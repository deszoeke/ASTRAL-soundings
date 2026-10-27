# Idealized May Arabian Sea sounding for CRM initial conditions (see case_plan.md).
# Standalone: parameters were fit in mean_sonde.ipynb and copied here; this file does not depend on the notebook.
#
#   using IdealizedSounding            # Pkg.develop(path="IdealizedSounding") once
#   s = idealized_sounding(p)          # (; p, θ, T, Td, rh, r, q, u, v) at arbitrary pressures p [hPa]
#   v = idealized_vertices()           # the same at the vertices (significant levels) of the profile
#   using PythonPlot                   # loads the plotting extension (ext/PlotExt.jl)
#   fig = plot_idealized_sounding()    # profiles with the vertices marked
#   keyword arguments override the parameters, e.g. idealized_sounding(p; θtop2=310.5)

module IdealizedSounding

using ForwardDiff

export sounding_params, idealized_sounding, idealized_vertices, plot_idealized_sounding

# ============================== parameters ==============================

"""
Default parameters of the idealized sounding (mean_sonde.ipynb, May Arabian Sea soundings). Pressures [hPa], θ [K].

Each mixture layer pbot ≥ p ≥ ptop is θ(p) = Σ dF max(θbot, θm(p)), θm the moist adiabat through θtop at ptop.
The branch at fit-grid level pv leaves θm as a dry adiabat θbot = θm(pv) downward, with weight dF (Σ dF = 1).
"""
const sounding_params = (
    # lower mixture, 810-1010 hPa (notebook: dry-adiabat mixture fit, "nonnegative least squares ...")
    pbot2=1010.0, ptop2=810.0, θtop2=309.98,
    pv2=[1010.0, 950.0, 930.0, 890.0, 870.0, 850.0, 830.0, 810.0],                  # fit-grid levels pc
    dF2=[0.501, 0.017, 0.200, 0.014, 0.117, 0.053, 0.077, 0.021],
    # upper mixture, 590-750 hPa (same notebook cell); θ linear in p in the gap 810-750 hPa
    pbot1=750.0, ptop1=590.0, θtop1=321.52,
    pv1=[750.0, 710.0, 690.0, 650.0, 630.0, 590.0],                                 # fit-grid levels pc
    dF1=[0.698, 0.012, 0.024, 0.043, 0.049, 0.174],
    # above ptop1 (notebook: "# idealize the mean profile above 590 hPa ..."): θ on the moist adiabat labeled
    #   θes = θes1 + dθes1dp (p - ptop1)                                        for pexp ≤ p < ptop1
    #   θes = θes(pexp) + dθes2dp (p - pexp) + Aexp (exp((pexp - p)/Hexp) - 1)  for ptop ≤ p < pexp
    θes1=347.8828150278081, dθes1dp=-0.008026126039451558,                           # [K], [K/hPa]
    pexp=400.0, dθes2dp=-0.024167949346775845, Aexp=0.00098085422382202, Hexp=30.0,  # [hPa], [K/hPa], [K], [hPa]
    ptop=10.0,
    # humidity and wind at the vertices ph, linear in log p between them
    # (notebook: "# plot dewpoint depression and RH profile", "# plot wind u,v profiles")
    ph=[1010.0, 1000, 950, 850, 810, 650,  500, 440, 400,  330, 210,  100,   70,   30,    10],
    dd=[   4.0,    4, 3.5,   6,   8,   8, 14.5,  16,  16, 13.8,  12, 12.6, 13.8,   24,    31],   # T - Td [K]
    uh=[   2.6,  3.1, 6.6, 7.1, 5.9, 3.2,  4.6, 3.8, 2.9,  1.5, -3.6, -16.2, -9.3, -13.7, -18.6],  # [m/s]
    vh=[  -1.1, -1.0, -1.5, -1.8, -1.9, -1.4, -1.3, -0.8, -0.6, -0.4, 1.3, -0.1, -0.1, -0.3, 1.0],  # [m/s]
)

# ============================== construction ==============================

const κ_id = 287.04 / 1005.7          # Rd/cp, as in sonde_averaging.jl

"Bolton (1980) saturation vapor pressure [hPa] at T [K]"
es_bolton(T) = 6.112exp(17.67(T - 273.15) / (T - 29.65))

"Bolton (1980) saturation θe [K] at T [K], p [hPa] (eq. 39, Bolton es)"
function θes_bolton(T, p)
    e = es_bolton(T)
    r = 0.622e / (p - e)
    T * (1000 / (p - e))^κ_id * exp((3036 / T - 1.78) * r * (1 + 0.448r))
end

"T [K] on the moist adiabat labeled θes [K] at p [hPa]: Newton with ForwardDiff"
function T_moist(θes, p; tol=1e-4)
    T = θes * (p / 1000)^κ_id - 30.0
    for _ in 1:20
        dT = (θes_bolton(T, p) - θes) / ForwardDiff.derivative(t -> θes_bolton(t, p), T)
        T -= dT
        abs(dT) < tol && break
    end
    T
end

"θ [K] on the moist adiabat labeled θes at p [hPa]"
θ_moist(θes, p) = T_moist(θes, p) * (1000 / p)^κ_id

"θes [K] labeling the moist adiabat through potential temperature θ [K] at p [hPa]"
θes_through(p, θ) = θes_bolton(θ * (p / 1000)^κ_id, p)

"mixture θ(p) [K]: dry adiabats leaving the moist adiabat θes at levels pv, weights dF"
function mix(pv, dF, θes, p)
    θm = θ_moist(θes, p)
    sum(w * max(θ_moist(θes, pb), θm) for (pb, w) in zip(pv, dF))
end

"linear interpolation of x(p) in log p; NaN outside the data"
function interp_logp(p, x, pp)
    o = sortperm(p); lp, x = log.(p[o]), x[o]
    l = log(pp)
    (l < lp[1] || l > lp[end]) && return NaN
    j = min(searchsortedlast(lp, l), length(lp) - 1)
    x[j] + (x[j+1] - x[j]) * (l - lp[j]) / (lp[j+1] - lp[j])
end

"θ [K] at pressure p [hPa] for parameters s; NaN below pbot2 or above ptop"
function θ_profile(p, s)
    θesm1, θesm2 = θes_through(s.ptop1, s.θtop1), θes_through(s.ptop2, s.θtop2)     # capping moist adiabats
    θb, θt = mix(s.pv1, s.dF1, θesm1, s.pbot1), mix(s.pv2, s.dF2, θesm2, s.ptop2)   # gap end points
    θesx = s.θes1 + s.dθes1dp * (s.pexp - s.ptop1)                                  # θes at pexp
    !(s.ptop <= p <= s.pbot2) ? NaN :
    p >= s.ptop2  ? mix(s.pv2, s.dF2, θesm2, p) :
    p > s.pbot1   ? θb + (θt - θb) * (p - s.pbot1) / (s.ptop2 - s.pbot1) :
    p >= s.ptop1  ? mix(s.pv1, s.dF1, θesm1, p) :
    θ_moist(p >= s.pexp ? s.θes1 + s.dθes1dp * (p - s.ptop1) :
            θesx + s.dθes2dp * (p - s.pexp) + s.Aexp * (exp((s.pexp - p) / s.Hexp) - 1), p)
end

"""
profile(p, s) -> (; p, θ, T, Td, rh, r, q, u, v) at pressures p [hPa] for parameters s.
Td = T - dd, dd interpolated linearly in log p from ph; rh = e/es(T), mixing ratio r and specific humidity q [kg/kg]
from e = es(Td) (Bolton 1980 es); u, v [m/s] linear in log p from ph.
"""
function profile(p, s)
    θ = [θ_profile(pp, s) for pp in p]
    T = θ .* (p ./ 1000) .^ κ_id
    Td = T .- [interp_logp(s.ph, s.dd, pp) for pp in p]
    e = es_bolton.(Td)                                                              # vapor pressure [hPa]
    (p=p, θ=θ, T=T, Td=Td, rh=e ./ es_bolton.(T),
     r=0.622 .* e ./ (p .- e), q=0.622 .* e ./ (p .- 0.378 .* e),
     u=[interp_logp(s.ph, s.uh, pp) for pp in p], v=[interp_logp(s.ph, s.vh, pp) for pp in p])
end

# ============================== user functions ==============================

"""
idealized_sounding(p; kw...) -> (; p, θ, T, Td, rh, r, q, u, v)
Idealized sounding at arbitrary pressures p [hPa]; kw override sounding_params.
θ, T, Td [K]; rh [0-1]; r, q [kg/kg]; u, v [m/s].
"""
idealized_sounding(p; kw...) = profile(collect(float.(p)), merge(sounding_params, kw))

"""
idealized_vertices(; kw...) -> (; p, θ, T, Td, rh, r, q, u, v)
Idealized sounding at its vertices (descending p): the mixture branch levels pv1, pv2, the layer ends, pexp, ptop,
and the humidity/wind levels ph. kw override sounding_params.
"""
function idealized_vertices(; kw...)
    s = merge(sounding_params, kw)
    pv = sort(unique(vcat(s.pv2, s.pbot2, s.ptop2, s.pv1, s.pbot1, s.ptop1, s.pexp, s.ptop, s.ph)); rev=true)
    profile(pv, s)
end

"plot_idealized_sounding(; kw...) -> fig: methods in ext/PlotExt.jl, loaded with `using PythonPlot`"
function plot_idealized_sounding end

end # module IdealizedSounding
