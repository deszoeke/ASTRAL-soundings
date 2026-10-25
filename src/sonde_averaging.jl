# Three ways to average soundings from all_profiles() (sonde_profiles.jl):
#   1. histogram of T, Td realizations (one per sounding) in 0.1 K, 10 hPa bins
#   2. mean on a 10 hPa grid: average high-resolution, interpolate low-resolution soundings
#   3. composite in normalized-height layers bounded by the shallow cumulus inversion (SCI),
#      the midtropospheric mixed layer (MML), and the midtropospheric inversion (MI)
#
#   include("sonde_averaging.jl"); using .SondeAveraging
#   as   = filter(...)                   # choose soundings in the notebook
#   h    = hist_T_Td(as)                 # method 1
#   m    = mean_on_pgrid(as, :θ)         # method 2
#   c    = composite_layers(as, (:θ, :u, :v))   # method 3
#
# Profiles are the NamedTuples from sonde_profiles.jl:
#   (source, station, time, lat, lon, p, z, θ, rh, u, v, z850, nfree)
# T and Td [°C] and q [kg/kg] are derived from p, θ, rh (Buck 1981). Fields are given as Symbols (:θ, :T, :Td, :rh, :u, :v).

module SondeAveraging

import Statistics, Dates
include(joinpath(@__DIR__, "inversions.jl")); using .Inversions

export T_Td, field_values, colstats,
       to_pgrid, profiles_on_pgrid, mean_on_pgrid, hist_pgrid, hist_field, hist_T_Td,
       layer_bounds, layer_levels, layer_ζ, composite_layers, LAYERS, DEFAULT_NPER

const κ = 287.04/1005.7          # Rd/cp, as in sonde_profiles.jl

# ---------------- common helpers ----------------

const ε = 0.622                  # Rd/Rv

"Buck (1981) saturation vapor pressure over water [hPa], T in °C, at pressure p [hPa]"
saturation_vp(T_C, p) = (1.0007 + 3.46e-6p) * 6.1121exp(17.502T_C/(240.97 + T_C))

"temperature [°C] from pressure [hPa] and potential temperature [K]"
temperature(p, θ) = θ * (p/1000)^κ - 273.15

"vapor pressure [hPa] from specific humidity q [kg/kg] and pressure p [hPa]"
vapor_pressure(q, p) = q*p / (ε + (1 - ε)*q)

"specific humidity [kg/kg] from vapor pressure e and pressure p [hPa]"
specific_humidity(e, p) = ε*e / (p - (1 - ε)*e)

"""
dewpoint(T_C, p; rh=NaN, q=NaN) dewpoint [°C] from temperature [°C], pressure [hPa], and
either relative humidity rh [%] or specific humidity q [kg/kg] (q is used if both are given).
Buck (1981): e_w = f(p) 6.1121 exp(17.502 T/(240.97 + T)), f = 1.0007 + 3.46e-6 p, inverted
exactly for Td; q is converted by e = q p/(ε + (1 - ε) q). NaN if humidity is missing or ≤ 0.
"""
function dewpoint(T_C, p; rh=NaN, q=NaN)
    if isfinite(q)
        q > 0 || return NaN
        e = vapor_pressure(q, p)
    elseif isfinite(rh)
        rh > 0 || return NaN
        e = rh/100 * saturation_vp(T_C, p)
    else
        return NaN
    end
    l = log(e / ((1.0007 + 3.46e-6p) * 6.1121))
    240.97l / (17.502 - l)
end

"(T, Td) [°C] at the levels of a profile"
function T_Td(prof)
    T = temperature.(prof.p, prof.θ)
    T, map((t, p, rh) -> dewpoint(t, p; rh=rh), T, prof.p, prof.rh)
end

"values of field f (:T, :Td, :q [kg/kg], or any profile field such as :θ, :rh, :u, :v) at the profile's levels"
function field_values(prof, f::Symbol)
    f === :T  && return T_Td(prof)[1]
    f === :Td && return T_Td(prof)[2]
    f === :q  && return specific_humidity.(saturation_vp.(T_Td(prof)[1], prof.p) .* prof.rh ./ 100, prof.p)
    float.(getproperty(prof, f))
end

"""
colstats(X) -> (mean, std, n) over soundings (rows) for each level (column) of X,
ignoring non-finite values. NaN where there are no (or, for std, fewer than 2) values.
"""
function colstats(X)
    nlev = size(X, 2)
    μ = fill(NaN, nlev); σ = fill(NaN, nlev); n = zeros(Int, nlev)
    for j in 1:nlev
        v = filter(isfinite, X[:, j])
        n[j] = length(v)
        n[j] > 0 && (μ[j] = Statistics.mean(v))
        n[j] > 1 && (σ[j] = Statistics.std(v))
    end
    (mean=μ, std=σ, n=n)
end

# ---------------- methods 1 and 2: pressure-grid averaging ----------------

"""
to_pgrid(p, x, pcenters; halfwidth=5, maxgap=Inf) -> x on the pressure grid `pcenters` [hPa].
Bins with at least 2 samples in [c-halfwidth, c+halfwidth) get their mean (higher resolution
than the grid); other bins are linearly interpolated to the bin center (lower resolution),
NaN outside the sampled range and across gaps wider than maxgap [hPa].
"""
function to_pgrid(p, x, pcenters; halfwidth=5.0, maxgap=Inf)
    ok = isfinite.(p) .& isfinite.(x)
    ps, xs = float.(p[ok]), float.(x[ok])
    o = sortperm(ps); ps, xs = ps[o], xs[o]
    xi = interp1(ps, xs, collect(float.(pcenters)); maxgap=maxgap)
    for (k, c) in enumerate(pcenters)
        i1 = searchsortedfirst(ps, c - halfwidth)
        i2 = searchsortedfirst(ps, c + halfwidth) - 1
        i2 - i1 >= 1 && (xi[k] = Statistics.mean(@view xs[i1:i2]))
    end
    xi
end

"matrix X[sounding, pressure level] of field f on the pressure grid"
function profiles_on_pgrid(profs, f::Symbol; pcenters=10.0:10.0:1020.0, kw...)
    X = fill(NaN, length(profs), length(pcenters))
    for (i, prof) in enumerate(profs)
        X[i, :] = to_pgrid(prof.p, field_values(prof, f), pcenters; kw...)
    end
    X
end

"""
mean_on_pgrid(profs, f; pcenters=10:10:1020, halfwidth=5, maxgap=Inf)
method 2: mean, std, n at each pressure level over soundings, plus the matrix X and pcenters.
"""
function mean_on_pgrid(profs, f::Symbol; pcenters=10.0:10.0:1020.0, kw...)
    X = profiles_on_pgrid(profs, f; pcenters=pcenters, kw...)
    merge(colstats(X), (X=X, pcenters=pcenters))
end

"""
hist_pgrid(X; dx=0.1, xlims=(-90, 40)) -> (counts[xbin, level], xedges)
2-D histogram of the columns of X[sounding, level]: one count per sounding per level.
"""
function hist_pgrid(X; dx=0.1, xlims=(-90.0, 40.0))
    nx = round(Int, (xlims[2] - xlims[1])/dx)
    H = zeros(Int, nx, size(X, 2))
    for j in axes(X, 2), i in axes(X, 1)
        x = X[i, j]
        isfinite(x) || continue
        k = floor(Int, (x - xlims[1])/dx) + 1
        1 <= k <= nx && (H[k, j] += 1)
    end
    H, collect(range(xlims[1], xlims[2]; length=nx+1))
end

"""
hist_field(profs, f; dx=nothing, xlims=nothing, edges=nothing, pcenters=10:10:1020, halfwidth=5, maxgap=Inf)
histogram of field f (:θ, :T, :Td, :q, :rh, :u, :v, :z, or any profile field) at each pressure level,
one realization per sounding, in dx bins over xlims. By default xlims is the range of the data and dx
is 1/100 of it. Pass increasing `edges` instead for unequal bins (e.g. geometric edges for q); values
outside the edges are not counted. Returns (H, xedges, pcenters, pedges, n) with H[xbin, level].
"""
function hist_field(profs, f::Symbol; dx=nothing, xlims=nothing, edges=nothing, pcenters=10.0:10.0:1020.0,
                    halfwidth=5.0, kw...)
    X = profiles_on_pgrid(profs, f; pcenters=pcenters, halfwidth=halfwidth, kw...)
    pedges = collect(range(first(pcenters) - halfwidth, last(pcenters) + halfwidth; length=length(pcenters)+1))
    if edges !== nothing
        xedges = collect(float.(edges))
        H = zeros(Int, length(xedges) - 1, size(X, 2))
        for j in axes(X, 2), i in axes(X, 1)
            x = X[i, j]
            isfinite(x) || continue
            k = searchsortedlast(xedges, x)
            1 <= k < length(xedges) && (H[k, j] += 1)
        end
        return (H=H, xedges=xedges, pcenters=pcenters, pedges=pedges, n=length(profs))
    end
    if xlims === nothing
        v = filter(isfinite, X)
        lo, hi = isempty(v) ? (0.0, 1.0) : extrema(v)
        hi > lo || (lo -= 0.5; hi += 0.5)
        xlims = (lo, hi + 1e-6*(hi - lo))        # include the maximum in the last bin
    end
    dx === nothing && (dx = (xlims[2] - xlims[1])/100)
    H, xedges = hist_pgrid(X; dx=dx, xlims=xlims)
    (H=H, xedges=xedges, pcenters=pcenters, pedges=pedges, n=length(profs))
end

"""
hist_T_Td(profs; dT=0.1, Tlims=(-90, 40), pcenters=10:10:1020, halfwidth=5, maxgap=Inf)
method 1: histograms of T and Td [°C] (dT bins) and specific humidity q [kg/kg] (dq = 1e-4
bins) at each pressure level, one realization per sounding, returned as
(T, Td, q, Tedges, qedges, pcenters, pedges, n). q is conserved, so it averages linearly;
Td, which has a sensible range as T and saturation q fall with height, is the better axis
for display. See hist_field for any other variable.
"""
function hist_T_Td(profs; dT=0.1, Tlims=(-90.0, 40.0), dq=1e-4, qlims=(0.0, 0.03),
                   pcenters=10.0:10.0:1020.0, halfwidth=5.0, kw...)
    kws = (; pcenters=pcenters, halfwidth=halfwidth, kw...)
    hT  = hist_field(profs, :T;  dx=dT, xlims=Tlims, kws...)
    hTd = hist_field(profs, :Td; dx=dT, xlims=Tlims, kws...)
    hq  = hist_field(profs, :q;  dx=dq, xlims=qlims, kws...)
    (T=hT.H, Td=hTd.H, q=hq.H, Tedges=hT.xedges, qedges=hq.xedges,
     pcenters=pcenters, pedges=hT.pedges, n=length(profs))
end

# ---------------- method 3: layer composite ----------------

"layers from the surface up: shallow cumulus layer, SCI, transition, MML, MI"
const LAYERS = (:SCL, :SCI, :TRN, :MML, :MI)
"normalized levels per layer in LAYERS"
const DEFAULT_NPER = (20, 5, 8, 20, 5)

"""
layer_bounds(prof; nfreemin=15, sci_zlims=(300, 4000), tol=200, sci_kw=(;), ml_kw=(;))
Layer boundaries of one sounding from inversions.jl, or `nothing` if the sounding is too
coarse (nfree < nfreemin) or a layer is not found. Heights [m] are geopotential.

  SCI   lowest inversion with base within sci_zlims above the surface (find_inversions)
  MML   lowest mixed layer starting within tol of the SCI top or above it, and having a top
        inversion (find_mixed_layers)
  MI    that top inversion

Returns (z0, sci_b, sci_t, ml_b, ml_t, mi_t, sci_Δθ, mi_Δθ). The shallow cumulus layer is
z0 to sci_b; the transition layer lies between sci_t and ml_b. sci_kw and ml_kw override
the detection thresholds (but not zlims).
"""
function layer_bounds(prof; nfreemin=15, sci_zlims=(300.0, 4000.0), tol=200.0,
                      sci_kw=(;), ml_kw=(;))
    prof.nfree >= nfreemin || return nothing
    z, θ = prof.z, prof.θ
    z0 = z[argmax(prof.p)]
    invs = find_inversions(z, θ; zlims=(z0 + sci_zlims[1], z0 + sci_zlims[2]), sci_kw...)
    isempty(invs) && return nothing
    sci = invs[argmin([i.zb for i in invs])]
    mls = filter(m -> m.ml_b >= sci.zt - tol && isfinite(m.top_t) && m.top_t > m.ml_t,
                 find_mixed_layers(z, θ; zlims=(z0 + sci_zlims[1], 10e3), ml_kw...))
    isempty(mls) && return nothing
    ml = mls[argmin([m.ml_b for m in mls])]
    (z0=z0, sci_b=sci.zb, sci_t=sci.zt, ml_b=max(ml.ml_b, sci.zt), ml_t=ml.ml_t, mi_t=ml.top_t,
     sci_Δθ=sci.Δθ, mi_Δθ=ml.top_Δθ)
end

layer_edges(b) = (b.z0, b.sci_b, b.sci_t, b.ml_b, b.ml_t, b.mi_t)

"normalized coordinate of the layer levels: layer l (1-based) spans l-1 to l"
layer_ζ(nper=DEFAULT_NPER) = reduce(vcat, [(l - 1) .+ ((1:n) .- 0.5)/n for (l, n) in enumerate(nper)])

"heights z of the normalized levels for boundaries b: nper levels equally spaced inside each layer"
function layer_levels(b, nper=DEFAULT_NPER)
    e = layer_edges(b)
    reduce(vcat, [e[l] .+ ((1:n) .- 0.5)/n .* (e[l+1] - e[l]) for (l, n) in enumerate(nper)])
end

"""
composite_layers(profs, fields=(:θ,); nper=DEFAULT_NPER, maxgap=1500, kw...)
method 3: interpolate each sounding with layer_bounds to the normalized layer levels, then
average over soundings. kw are passed to layer_bounds. Returns
(ζ, mean, std, n, X, zmean, edges, ok, nrejected), where mean/std/n/X are NamedTuples keyed by
`fields` (X[f][sounding, level] holds the composite values, for histograms), zmean is the mean height of each level, edges[sounding, 6] the layer boundaries,
and ok the indices of the soundings used.
"""
function composite_layers(profs, fields=(:θ,); nper=DEFAULT_NPER, maxgap=1500.0, kw...)
    fields = Tuple(fields)
    bs = [layer_bounds(p; kw...) for p in profs]
    ok = findall(!isnothing, bs)
    nl = sum(nper)
    zs = fill(NaN, length(ok), nl)
    edges = fill(NaN, length(ok), 6)
    X = Dict(f => fill(NaN, length(ok), nl) for f in fields)
    for (r, i) in enumerate(ok)
        zs[r, :] = layer_levels(bs[i], nper)
        edges[r, :] = collect(layer_edges(bs[i]))
        for f in fields
            X[f][r, :] = interp1(profs[i].z, field_values(profs[i], f), zs[r, :]; maxgap=maxgap)
        end
    end
    st = map(f -> colstats(X[f]), fields)
    (ζ=layer_ζ(nper),
     mean=NamedTuple{fields}(map(s -> s.mean, st)),
     std=NamedTuple{fields}(map(s -> s.std, st)),
     n=NamedTuple{fields}(map(s -> s.n, st)),
     X=NamedTuple{fields}(map(f -> X[f], fields)),
     zmean=colstats(zs).mean, edges=edges, ok=ok, nrejected=length(profs) - length(ok))
end

end # module SondeAveraging
