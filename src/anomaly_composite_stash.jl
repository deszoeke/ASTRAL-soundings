# Stash of the anomaly-compositing experiment (removed from mean_sonde.ipynb).
# Each section below was one notebook cell, in this order. To reproduce step by step,
# paste the cells back into mean_sonde.ipynb (or `include` the sections one at a time)
# after running its cells 0-5 (needs `as`, `c`, `cmap`/`fontsize`, SondeAveraging, PythonPlot).
# Status: not convincing -- the θ′ std stayed large (3-6 K) after aligning layers; see the
# diagnostic in D. Kept so the procedure can be rebuilt, not as a recommended method.

# ======================================================================
# CELL A. first attempt: anomalies relative to each sounding's own shallow-cumulus-layer mean (superseded by B-D)
# ======================================================================
# composite anomalies: remove each sounding's own reference (its mean over the shallow cumulus layer)
# so the θ, q offsets between soundings don't swamp the inversion / mixed-layer structure
nper = SondeAveraging.DEFAULT_NPER
refl = 1:nper[1]                                 # composite levels used as the reference (shallow cumulus layer)

anomaly(X) = X .- [mean(filter(isfinite, X[r, refl])) for r in axes(X, 1)]

ca = map((:θ, :q)) do f
    A = anomaly(c.X[f])
    s = SondeAveraging.colstats(A)
    (A=A, mean=s.mean, std=s.std, n=s.n)
end
scale = (θ=1.0, q=1e3)                           # q in g/kg
units_ = (θ="θ′ [K]", q="q′ [g/kg]")

fig, axs = subplots(1, 2; figsize=(10, 6), sharey=true)
for (ax, f, a) in zip(axs, (:θ, :q), ca)
    μ, σ = scale[f] .* a.mean, scale[f] .* a.std
    ax.fill_betweenx(c.ζ, μ .- σ, μ .+ σ; alpha=0.3, label="± σ")
    ax.fill_betweenx(c.ζ, μ .- σ ./ sqrt.(a.n), μ .+ σ ./ sqrt.(a.n); alpha=0.6, label="± σ/√n")
    ax.plot(μ, c.ζ; label="mean")
    ax.set_xlabel(units_[f]; fontsize=fontsize)
    ax.set_yticks(0:length(nper))
    ax.tick_params(axis="both"; labelsize=fontsize)
end
axs[0].set_ylabel("normalized layer height"; fontsize=fontsize)
axs[0].legend(; fontsize=fontsize)
gcf()

# ======================================================================
# CELL B. functions: reference profile, anomalies, layers from θ′, anomaly composite
# ======================================================================
# anomaly compositing, prototyped here to move into sonde_averaging.jl later:
#   1. reference = mean profile of the valid soundings on a height-above-surface grid
#   2. anomaly = sounding - reference, in physical height
#   3. layers are found from the θ anomaly itself: high over low θ′ (dθ′/dz > 0) is a stable layer,
#      low over high θ′ (dθ′/dz < 0) is a mixed layer
#   4. anomalies are stretched onto the layer levels, then averaged (mean, std, n) over soundings
Inv = SondeAveraging.Inversions
HGRID = collect(0.0:50.0:12e3)                  # height above each sounding's surface [m]

"field f of a profile on HGRID (height above the surface = lowest level)"
function on_hgrid(prof, f; maxgap=1500.0)
    h = prof.z .- prof.z[argmax(prof.p)]
    Inv.interp1(h, SondeAveraging.field_values(prof, f), HGRID; maxgap=maxgap)
end

"mean (reference) profile of f on HGRID, NaN where fewer than nmin soundings; X[sounding, level] on HGRID"
function reference_profile(profs, f; nmin=10, maxgap=1500.0)
    X = permutedims(reduce(hcat, [on_hgrid(p, f; maxgap=maxgap) for p in profs]))
    s = SondeAveraging.colstats(X)
    (h=HGRID, mean=ifelse.(s.n .>= nmin, s.mean, NaN), n=s.n, X=X)
end

"""
anomaly_layers(h, a; ...) -> (z0, sci_b, sci_t, ml_b, ml_t, mi_t) or nothing
Layer boundaries [m above surface] from one θ anomaly profile a(h) alone:
stable layer: smoothed dθ′/dz ≥ γs [K/m] with θ′ rise ≥ Δmin; mixed layer: dθ′/dz ≤ -γm, ≥ mldmin deep.
SCI = lowest stable layer with base in sci_zlims; MML = lowest mixed layer above it (within tol) that
has a stable layer (MI) starting within `gap` above its top.
"""
function anomaly_layers(h, a; γs=2e-3, γm=2e-3, Δmin=1.0, mldmin=400.0, smooth=100.0,
                        mergegap=100.0, mlgap=200.0, sci_zlims=(300.0, 4000.0), tol=200.0, gap=500.0)
    dz = h[2] - h[1]
    g = Inv.ddz(a, dz; nk=max(1, round(Int, smooth/dz)))
    stab = [(max(i-1, 1), min(j+1, length(h))) for (i, j) in
            Inv.runs(g .>= γs; mergegap=round(Int, mergegap/dz))]
    filter!(r -> (d = a[r[2]] - a[r[1]]; isfinite(d) && d >= Δmin), stab)
    mix = filter(r -> h[r[2]] - h[r[1]] >= mldmin, Inv.runs(g .<= -γm; mergegap=round(Int, mlgap/dz)))
    sci = filter(r -> sci_zlims[1] <= h[r[1]] <= sci_zlims[2], stab)
    isempty(sci) && return nothing
    s = sci[argmin([h[r[1]] for r in sci])]
    ngap = round(Int, gap/dz)
    for m in sort(mix; by=first)
        h[m[1]] >= h[s[2]] - tol || continue
        k = findfirst(r -> r[2] > m[2] && r[1] <= m[2] + ngap + 1, stab)
        isnothing(k) && continue
        return (z0=0.0, sci_b=h[s[1]], sci_t=h[s[2]], ml_b=max(h[m[1]], h[s[2]]),
                ml_t=h[m[2]], mi_t=h[stab[k][2]])
    end
    nothing
end

"""
anomaly_composite(profs, fields=(:θ,); layerfield=:θ, nper, nfreemin=15, nmin=10, maxgap=1500, layer_kw=(;))
Reference profile and anomalies of each field on HGRID, layers from the `layerfield` anomaly,
anomalies stretched to the layer levels. Returns (ζ, zmean, mean, std, n, A, ref, anom, edges, ok, bs, profs, nrejected):
mean/std/n/A (A[f][sounding, level]) are composites of the anomalies keyed by field; ref[f] is the reference
profile; anom[f][sounding, h] the anomalies on HGRID; ok indexes `profs` (the soundings with nfree ≥ nfreemin).
"""
function anomaly_composite(profs, fields=(:θ,); layerfield=:θ, nper=SondeAveraging.DEFAULT_NPER,
                           nfreemin=15, nmin=10, maxgap=1500.0, layer_kw=(;))
    profs = filter(p -> p.nfree >= nfreemin, profs)
    fields = Tuple(unique((layerfield, fields...)))
    refs = map(f -> reference_profile(profs, f; nmin=nmin, maxgap=maxgap), fields)
    anom = map(r -> r.X .- r.mean', refs)                       # anomalies in physical height
    il = findfirst(==(layerfield), fields)
    bs = [anomaly_layers(HGRID, anom[il][r, :]; layer_kw...) for r in axes(anom[il], 1)]
    ok = findall(!isnothing, bs)
    nl = sum(nper)
    zs = fill(NaN, length(ok), nl)
    edges = fill(NaN, length(ok), 6)
    A = map(_ -> fill(NaN, length(ok), nl), anom)
    for (r, i) in enumerate(ok)
        zs[r, :] = SondeAveraging.layer_levels(bs[i], nper)
        edges[r, :] = collect(SondeAveraging.layer_edges(bs[i]))
        for k in eachindex(anom)
            A[k][r, :] = Inv.interp1(HGRID, anom[k][i, :], zs[r, :]; maxgap=maxgap)
        end
    end
    st = map(SondeAveraging.colstats, A)
    nt(x) = NamedTuple{fields}(Tuple(x))
    (ζ=SondeAveraging.layer_ζ(nper), zmean=SondeAveraging.colstats(zs).mean,
     mean=nt(map(s -> s.mean, st)), std=nt(map(s -> s.std, st)), n=nt(map(s -> s.n, st)),
     A=nt(A), ref=nt(refs), anom=nt(anom), edges=edges, ok=ok, bs=bs, profs=profs,
     nrejected=length(profs) - length(ok))
end

# ======================================================================
# CELL C. run anomaly_composite on the May Arabian Sea soundings and plot
# ======================================================================
# composite of θ and q anomalies from the mean (reference) profile, layers found from θ′
ac = anomaly_composite(as, (:θ, :q))
println(length(ac.ok), " soundings used, ", ac.nrejected, " rejected (no stable layer / mixed layer found in θ′)")
for (k, name) in enumerate(("surface", "SCI base", "SCI top", "MML base", "MML top", "MI top"))
    println(rpad(name, 9), round(mean(ac.edges[:, k]); digits=0), " m above surface")
end

fig, axs = subplots(1, 3; figsize=(18, 7))
# a few individual θ′ profiles (offset) with the detected layer boundaries
off = 4.0
for (k, r) in enumerate(round.(Int, range(1, length(ac.ok); length=6)))
    a = ac.anom.θ[ac.ok[r], :]
    axs[0].plot(a .+ k*off, HGRID)
    axs[0].plot(Inv.interp1(HGRID, a, ac.edges[r, :]) .+ k*off, ac.edges[r, :], "ko"; markersize=4)
end
axs[0].set_xlabel("θ′ [K] (offset)"; fontsize=fontsize); axs[0].set_ylabel("height above surface [m]"; fontsize=fontsize)
axs[0].set_ylim(0, 8000)

scale = (θ=1.0, q=1e3)                           # q in g/kg
for (ax, f, lab) in zip(axs[1:2], (:θ, :q), ("θ′ [K]", "q′ [g/kg]"))
    μ, σ, n = scale[f] .* ac.mean[f], scale[f] .* ac.std[f], ac.n[f]
    ax.fill_betweenx(ac.ζ, μ .- σ, μ .+ σ; alpha=0.3, label="± σ")
    ax.fill_betweenx(ac.ζ, μ .- σ ./ sqrt.(n), μ .+ σ ./ sqrt.(n); alpha=0.6, label="± σ/√n")
    ax.plot(μ, ac.ζ; label="mean")
    ax.axvline(0; color="k", linewidth=0.5)
    ax.set_xlabel(lab; fontsize=fontsize)
    ax.set_yticks(0:length(SondeAveraging.DEFAULT_NPER))
end
axs[1].set_ylabel("normalized layer height"; fontsize=fontsize)
axs[1].legend(; fontsize=fontsize)
for ax in axs; ax.tick_params(axis="both"; labelsize=fontsize); end
gcf()

# ======================================================================
# CELL D. diagnostic: stretched vs unstretched vs offset-removed std of θ′
# ======================================================================
# how much of the θ′ spread is inversion/layer position, and how much is a sounding-wide offset?
#   unstretched : θ′ at the composite mean heights, no layer alignment
#   stretched   : θ′ at the layer levels (what the composite uses)
#   offset removed : stretched, minus each sounding's own column-mean θ′
ok = ac.ok
Au = permutedims(reduce(hcat, [Inv.interp1(HGRID, ac.anom.θ[i, :], ac.zmean; maxgap=1500.0) for i in ok]))
As = ac.A.θ
Ao = As .- [mean(filter(isfinite, As[r, :])) for r in axes(As, 1)]
curves = (("unstretched, at mean heights", Au), ("stretched to layer levels", As), ("stretched, column-mean offset removed", Ao))

fig, ax = subplots(figsize=(7, 7))
for (lab, X) in curves
    ax.plot(SondeAveraging.colstats(X).std, ac.ζ; label=lab)
end
ax.set_xlabel("std of θ′ [K]"; fontsize=fontsize); ax.set_ylabel("normalized layer height"; fontsize=fontsize)
ax.set_yticks(0:length(SondeAveraging.DEFAULT_NPER)); ax.tick_params(axis="both"; labelsize=fontsize)
ax.legend(; fontsize=fontsize-4)
for (lab, X) in curves
    println(rpad(lab, 40), "rms std = ", round(sqrt(mean(filter(isfinite, SondeAveraging.colstats(X).std) .^ 2)); digits=2), " K")
end
gcf()
