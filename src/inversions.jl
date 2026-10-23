# Inversion (capping stable layer) diagnostics for single soundings.
# Works on irregular, possibly sparse levels (10-m RAOBE, IMD, UWyo significant levels).
#
#   include("inversions.jl"); using .Inversions
#   invs = find_inversions(z, θ)
#   d    = inversion_diagnostics(z, θ, u, v; lat=12.0, rh=rh)
#
# Layers are defined by potential temperature θ only.

module Inversions

import Statistics

export interp1, regrid, find_inversions, find_mixed_layers, jumps, inversion_diagnostics,
       margules_tilt, primary_inversion, relative_composite

const g = 9.8        # m/s^2
const Ω = 7.292e-5   # 1/s

"mean that returns NaN for an empty collection"
mean(x) = isempty(x) ? NaN : Statistics.mean(x)

"""
interp1(z, x, zi; maxgap=Inf) linear interpolation of x(z) to zi.
Drops non-finite (z, x) pairs, sorts by z, averages duplicate z.
NaN outside the data range and inside data gaps wider than maxgap [m].
"""
function interp1(z, x, zi; maxgap=Inf)
    ok = isfinite.(z) .& isfinite.(x)
    zs, xs = float.(z[ok]), float.(x[ok])
    p = sortperm(zs); zs, xs = zs[p], xs[p]
    # average duplicate levels
    zu = unique(zs)
    if length(zu) < length(zs)
        xs = [mean(xs[zs .== zz]) for zz in zu]; zs = zu
    end
    xi = fill(NaN, length(zi))
    length(zs) < 2 && return xi
    for (k, zz) in enumerate(zi)
        (zz < zs[1] || zz > zs[end]) && continue
        j = searchsortedlast(zs, zz)
        if j == length(zs)
            xi[k] = xs[end]
        elseif zs[j+1] - zs[j] <= maxgap
            w = (zz - zs[j]) / (zs[j+1] - zs[j])
            xi[k] = (1-w)*xs[j] + w*xs[j+1]
        end
    end
    xi
end

"regrid(z, x; dz=50, zlims=(0, 12e3), maxgap) -> (zg, xg) on a uniform grid."
function regrid(z, x; dz=50.0, zlims=(0.0, 12e3), maxgap=Inf)
    zg = zlims[1]:dz:zlims[2]
    zg, interp1(z, x, zg; maxgap=maxgap)
end

"centered derivative dx/dz over ±nk grid points of a uniform grid (NaN at ends)"
function ddz(x, dz; nk=1)
    d = fill(NaN, length(x))
    for k in 1+nk:length(x)-nk
        d[k] = (x[k+nk] - x[k-nk]) / (2nk*dz)
    end
    d
end

"contiguous true runs of a Bool vector as (first, last) index pairs, merging gaps ≤ mergegap points"
function runs(mask; mergegap=0)
    r = Tuple{Int,Int}[]
    k = 1
    while k <= length(mask)
        if mask[k]
            j = k
            while j < length(mask) && mask[j+1]; j += 1; end
            if !isempty(r) && k - r[end][2] - 1 <= mergegap
                r[end] = (r[end][1], j)
            else
                push!(r, (k, j))
            end
            k = j + 1
        else
            k += 1
        end
    end
    r
end

"""
find_inversions(z, θ; kw...) -> Vector of NamedTuples, one per stable layer.

A stable layer is a contiguous height range where dθ/dz ≥ γinv, smoothed over
±smooth meters, with total θ increase Δθ ≥ Δθmin. Each layer reports

    zb, zt       base and top height [m]
    Δθ           θ(zt) - θ(zb) [K];  γ = Δθ/(zt-zb) [K/m]
    θb           θ at base [K]
    mld          depth [m] of the weakly stratified layer directly below zb: the
                 deepest layer (≤ mldmax) whose bulk dθ/dz and that of its bottom
                 mlwin meters are both ≤ γml; capping = mld ≥ mldmin
    zml          base of that mixed layer [m]

Defaults are for the tropical free troposphere where the moist adiabat has
dθ/dz ≈ 4-5 K/km near 5 km: γinv = 8 K/km, γml = 2.5 K/km.
"""
function find_inversions(z, θ; dz=50.0, zlims=(500.0, 10e3), maxgap=1500.0,
                         smooth=100.0, γinv=8e-3, Δθmin=1.0, mergegap=100.0,
                         γml=2.5e-3, mldmin=300.0, mldmax=3000.0, mlwin=500.0)
    zg, θg = regrid(z, θ; dz=dz, zlims=zlims, maxgap=maxgap)
    nk = max(1, round(Int, smooth/dz))
    dθ = ddz(θg, dz; nk=nk)
    dθ1 = ddz(θg, dz; nk=1)   # unsmoothed, to place the layer edges

    out = NamedTuple[]
    for (i, j) in runs(coalesce.(dθ .>= γinv, false); mergegap=round(Int, mergegap/dz))
        # the smoothed gradient widens the layer by ±nk; trim to the unsmoothed stable part
        while i < j && !(dθ1[i] >= γinv); i += 1; end
        while j > i && !(dθ1[j] >= γinv); j -= 1; end
        i, j = max(i-1, 1), min(j+1, length(zg))
        zb, zt = zg[i], zg[j]
        Δθ = θg[j] - θg[i]
        (isfinite(Δθ) && Δθ >= Δθmin) || continue

        # weakly stratified layer below the base, by bulk gradients (robust to noise)
        mld = 0.0
        for k0 in i-1:-1:max(1, i - round(Int, mldmax/dz))
            kw = min(i, k0 + round(Int, mlwin/dz))
            bulk = (θg[i] - θg[k0]) / (zg[i] - zg[k0])
            bot  = (θg[kw] - θg[k0]) / (zg[kw] - zg[k0])
            (isfinite(bulk) && isfinite(bot)) || break
            bulk <= γml && bot <= γml && (mld = zg[i] - zg[k0])
        end

        push!(out, (zb=zb, zt=zt, Δθ=Δθ, γ=Δθ/max(zt-zb, dz), θb=θg[i],
                    mld=mld, zml=zb - mld, capping=mld >= mldmin))
    end
    out
end

"""
find_mixed_layers(z, θ; kw...) -> Vector of NamedTuples, one per mixed layer.

A mixed layer is a contiguous height range where dθ/dz, smoothed over ±smooth
meters, is ≤ γml (gaps ≤ mlgap merged), at least mldmin deep, with bulk
dθ/dz ≤ γml. Its bounding inversions are the stable layers (smoothed
dθ/dz ≥ γs, gaps ≤ mergegap merged) that begin within `gap` above the mixed-layer
top (top inversion) and end within `gap` below its base (base inversion).

    ml_b, ml_t, mld, θml   mixed-layer base, top, depth [m], mean θ [K]
    top_b, top_t, top_Δθ   top inversion base and top [m], and
                           θ(top_t) - θ(ml_t) [K]; NaN, NaN, 0 if there is none
    base_b, base_t, base_Δθ  base inversion base and top [m], and
                           θ(ml_b) - θ(base_b) [K]; NaN, NaN, 0 if there is none
    top_Δθ500, base_Δθ500  threshold-free strengths: θ(ml_t+500) - θ(ml_t) and
                           θ(ml_b) - θ(ml_b-500) [K]

The default γs = 8 K/km is well above the moist-adiabatic dθ/dz (3-6 K/km at
2-7 km), so free-tropospheric layers don't merge with the inversions.
"""
function find_mixed_layers(z, θ; dz=50.0, zlims=(500.0, 10e3), maxgap=1500.0, smooth=100.0,
                           γml=2.5e-3, mlgap=200.0, mldmin=400.0,
                           γs=8e-3, mergegap=100.0, gap=500.0)
    zg, θg = regrid(z, θ; dz=dz, zlims=zlims, maxgap=maxgap)
    dθ = ddz(θg, dz; nk=max(1, round(Int, smooth/dz)))
    dθ1 = ddz(θg, dz; nk=1)
    # stable cores from the smoothed gradient, edges extended while the unsmoothed one stays stable
    stab = map(runs(coalesce.(dθ .>= γs, false); mergegap=round(Int, mergegap/dz))) do (a, b)
        while a > 1 && coalesce(dθ1[a-1] >= γs, false); a -= 1; end
        while b < length(zg) && coalesce(dθ1[b+1] >= γs, false); b += 1; end
        (max(a-1, 1), min(b+1, length(zg)))
    end
    ngap = round(Int, gap/dz)

    out = NamedTuple[]
    for (i, j) in runs(coalesce.(dθ .<= γml, false); mergegap=round(Int, mlgap/dz))
        mld = zg[j] - zg[i]
        (mld >= mldmin && (θg[j] - θg[i])/mld <= γml) || continue
        # top inversion: stable run starting just above the mixed-layer top
        kt = findfirst(r -> r[2] > j && r[1] <= j + ngap + 1, stab)
        # base inversion: stable run ending just below the mixed-layer base
        kb = findlast(r -> r[1] < i && r[2] >= i - ngap - 1, stab)
        top  = isnothing(kt) ? (NaN, NaN, 0.0) :
               (zg[stab[kt][1]], zg[stab[kt][2]], θg[stab[kt][2]] - θg[j])
        base = isnothing(kb) ? (NaN, NaN, 0.0) :
               (zg[stab[kb][1]], zg[stab[kb][2]], θg[i] - θg[stab[kb][1]])
        θat(zz) = interp1(z, θ, [zz])[1]
        push!(out, (ml_b=zg[i], ml_t=zg[j], mld=mld, θml=mean(θg[i:j]),
                    top_b=top[1], top_t=top[2], top_Δθ=top[3],
                    base_b=base[1], base_t=base[2], base_Δθ=base[3],
                    top_Δθ500=θat(zg[j] + 500) - θg[j], base_Δθ500=θg[i] - θat(zg[i] - 500)))
    end
    out
end

"""
jumps(z, x, zb, zt; δ=300, zmin=-Inf) -> (below, above, Δ): means of x over
[max(zb-δ, zmin), zb] and [zt, zt+δ].
"""
function jumps(z, x, zb, zt; δ=300.0, dz=25.0, zmin=-Inf)
    below = zb > zmin ? mean(filter(isfinite, interp1(z, x, max(zb-δ, zmin):dz:zb))) : NaN
    above = mean(filter(isfinite, interp1(z, x, zt:dz:(zt+δ))))
    (below=below, above=above, Δ=above-below)
end

"""
margules_tilt(Δu, Δv, Δθ, θ0, lat) -> (dhdy, dhdx) slope of a balanced interface.

Two-layer thermal wind (Margules) across a θ jump with g′ = gΔθ/θ0:
    u_above - u_below =  (g′/f) ∂h/∂y
    v_above - v_below = -(g′/f) ∂h/∂x
"""
function margules_tilt(Δu, Δv, Δθ, θ0, lat)
    f = 2Ω*sind(lat)
    gp = g*Δθ/θ0
    (dhdy = f*Δu/gp, dhdx = -f*Δv/gp)
end

"least-squares slope dx/dz of x(z) over [z1, z2]"
function slope(z, x, z1, z2; dz=25.0)
    zi = z1:dz:z2
    xi = interp1(z, x, zi)
    ok = isfinite.(xi)
    sum(ok) < 3 && return NaN
    zz = zi[ok] .- mean(zi[ok]); xx = xi[ok]
    sum(zz .* (xx .- mean(xx))) / sum(zz.^2)
end

"""
inversion_diagnostics(z, θ, u, v; lat, rh=nothing, δ=300, zshear=-Inf, kw...)

find_inversions plus, for each layer, jumps across it (averaged over δ below
the base and above the top) in θ, u, v, rh; bulk Richardson number across the
layer; Margules tilt implied by the wind jump; and the mean shear in the mixed
layer below. Wind jumps and shears use only levels above zshear [m], to exclude
the frictional (ageostrophic) boundary layer, e.g. zshear = z(850 hPa).
kw are passed to find_inversions.
"""
function inversion_diagnostics(z, θ, u, v; lat, rh=nothing, δ=300.0, zshear=-Inf, kw...)
    map(find_inversions(z, θ; kw...)) do inv
        jθ = jumps(z, θ, inv.zb, inv.zt; δ=δ)
        ju = jumps(z, u, inv.zb, inv.zt; δ=δ, zmin=zshear)
        jv = jumps(z, v, inv.zb, inv.zt; δ=δ, zmin=zshear)
        jrh = isnothing(rh) ? (Δ=NaN,) : jumps(z, rh, inv.zb, inv.zt; δ=δ)
        θ0 = (jθ.below + jθ.above)/2
        H = inv.zt - inv.zb + δ   # distance between the centers of the averaging windows
        Ri = g*jθ.Δ/θ0 * H / (ju.Δ^2 + jv.Δ^2)
        tilt = margules_tilt(ju.Δ, jv.Δ, jθ.Δ, θ0, lat)
        mlz1 = max(inv.zml, inv.zb - 2000, zshear)   # mixed-layer shear: ≤ 2 km below base, above zshear
        dudz_ml = inv.zb - mlz1 >= 300 ? slope(z, u, mlz1, inv.zb) : NaN
        dvdz_ml = inv.zb - mlz1 >= 300 ? slope(z, v, mlz1, inv.zb) : NaN
        merge(inv, (Δθj=jθ.Δ, u_b=ju.below, u_a=ju.above, Δu=ju.Δ,
                    v_b=jv.below, v_a=jv.above, Δv=jv.Δ, Δrh=jrh.Δ,
                    Ri=Ri, dhdy=tilt.dhdy, dhdx=tilt.dhdx,
                    dudz_ml=dudz_ml, dvdz_ml=dvdz_ml))
    end
end

"""
primary_inversion(invs; zlims=(3000, 7000), capping=true)
strongest (largest Δθ) inversion with base in zlims, or nothing.
"""
function primary_inversion(invs; zlims=(3000.0, 7000.0), capping=true)
    c = filter(i -> zlims[1] <= i.zb <= zlims[2] && (!capping || i.capping), invs)
    isempty(c) ? nothing : c[argmax([i.Δθ for i in c])]
end

"""
relative_composite(zs, xs, zref; zrel=-2000:50:2000)
mean and count of profiles x(z - zref) over soundings. zs, xs, zref are vectors
(one entry per sounding); soundings with non-finite zref are skipped.
"""
function relative_composite(zs, xs, zref; zrel=-2000.0:50.0:2000.0)
    s = zeros(length(zrel)); n = zeros(Int, length(zrel))
    for (z, x, z0) in zip(zs, xs, zref)
        isfinite(z0) || continue
        xi = interp1(z, x, z0 .+ zrel)
        ok = isfinite.(xi)
        s[ok] .+= xi[ok]; n[ok] .+= 1
    end
    zrel, s ./ n, n
end

end # module Inversions
