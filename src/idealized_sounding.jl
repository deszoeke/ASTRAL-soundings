# Idealized May Arabian Sea sounding for CRM cases (see case_plan.md).
#   include("idealized_sounding.jl")
#   θ, T = idealized_sounding(p)          # all parameters have defaults

using ForwardDiff

const κ_id = 287.04 / 1005.7          # Rd/cp, as in sonde_averaging.jl

"Bolton (1980) saturation θe [K] at T [K], p [hPa] (eq. 39, Bolton es)"
function θes_bolton(T, p)
    e = 6.112exp(17.67(T - 273.15) / (T - 29.65))
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

"mixture of dry adiabats θbot with weights dF capped by the moist adiabat θes: Σ dF max(θbot, θm(p))"
mix(θbot, dF, θes, p) = sum(dF .* max.(θbot, θ_moist(θes, p)))

"θes [K] labeling the moist adiabat through potential temperature θ [K] at p [hPa]"
θes_through(p, θ) = θes_bolton(θ * (p / 1000)^κ_id, p)

"""
idealized_sounding(p; kw...) -> (θ, T) [K] at pressures p [hPa]

Standalone idealized May Arabian Sea sounding (defaults from mean_sonde.ipynb, May AS soundings).
Layers from the surface up:
  lower mixture  ps ≥ p ≥ pm2:   θ = Σ dF2 max(θbot2, θm2(p)), θm2 the moist adiabat through θm2top at pm2
  gap            pm2 > p > pb1:  θ linear in p between the two mixtures
  upper mixture  pb1 ≥ p ≥ p1:   θ = Σ dF1 max(θbot1, θm1(p)), θm1 the moist adiabat through θm1top at p1
  above p1       p < p1:         θes = up.b0 + up.b1 (p - p1) for p ≥ up.pq,
                                 θes(pq) + up.B (pq - p) + up.A (exp((pq - p)/up.H) - 1) above
Moist adiabats: Bolton (1980) θes, inverted for T by Newton/ForwardDiff (the notebook fit used metpy moist_lapse).
"""
function idealized_sounding(p;
        ps=1010.0,
        pm2=810.0, θm2top=309.98,                                                   # lower mixture
        θbot2=[298.38, 301.60, 302.72, 305.03, 306.23, 307.45, 308.70, 309.98],
        dF2=[0.501, 0.017, 0.200, 0.014, 0.117, 0.053, 0.077, 0.021],
        pb1=750.0,                                                                  # upper mixture
        p1=590.0, θm1top=321.52,
        θbot1=[310.27, 312.91, 314.27, 317.08, 318.53, 321.52],
        dF1=[0.698, 0.012, 0.024, 0.043, 0.049, 0.174],
        up=(pq=400.0, b0=missing, b1=missing, B=missing, A=missing, H=missing))      # above p1: run cell 13
    θes1, θes2 = θes_through(p1, θm1top), θes_through(pm2, θm2top)
    θb, θt = mix(θbot1, dF1, θes1, pb1), mix(θbot2, dF2, θes2, pm2)                # gap end points
    θ = map(p) do pp
        pp > ps   ? NaN :
        pp >= pm2 ? mix(θbot2, dF2, θes2, pp) :
        pp > pb1  ? θb + (θt - θb) * (pp - pb1) / (pm2 - pb1) :
        pp >= p1  ? mix(θbot1, dF1, θes1, pp) :
        θ_moist(pp >= up.pq ? up.b0 + up.b1 * (pp - p1) :
                (up.b0 + up.b1 * (up.pq - p1)) + up.B * (up.pq - pp) + up.A * (exp((up.pq - pp) / up.H) - 1), pp)
    end
    θ, θ .* (p ./ 1000) .^ κ_id
end

"linear interpolation of x(p) in log p; NaN outside the data"
function interp_logp(p, x, pp)
    o = sortperm(p); lp, x = log.(p[o]), x[o]
    l = log(pp)
    (l < lp[1] || l > lp[end]) && return NaN
    j = min(searchsortedlast(lp, l), length(lp) - 1)
    x[j] + (x[j+1] - x[j]) * (l - lp[j]) / (lp[j+1] - lp[j])
end

"""
fit_θes_upper(p, θes; pa=590, pq=400, pz=10, Hs=10:5:500) -> (; pa, pq, b0, b1, B, A, H)
Fit θes(p) above pa: linear for pq ≤ p ≤ pa, θes = b0 + b1 (p - pa); linear + exponential for pz ≤ p < pq,
θes = θq + B (pq - p) + A (exp((pq - p)/H) - 1), continuous at pq (θq = b0 + b1 (pq - pa)).
b0, b1, then B, A by least squares; H from a scan over Hs [hPa]. Evaluate with θes_upper(fit, p).
"""
function fit_θes_upper(p, θes; pa=590.0, pq=400.0, pz=10.0, Hs=10.0:5.0:500.0)
    kl, kq = findall(pq .<= p .<= pa), findall(pz .<= p .< pq)
    b0, b1 = [ones(length(kl)) (p[kl] .- pa)] \ θes[kl]
    y = θes[kq] .- (b0 + b1 * (pq - pa))
    fits = map(Hs) do H
        X = [(pq .- p[kq]) (exp.((pq .- p[kq]) ./ H) .- 1)]
        B, A = X \ y
        (B=B, A=A, H=H, r=sum(abs2, X * [B, A] .- y))
    end
    f = argmin(f -> f.r, fits)
    (pa=pa, pq=pq, b0=b0, b1=b1, B=f.B, A=f.A, H=f.H)
end

"θes [K] at p [hPa] (p ≤ pa) from a fit_θes_upper result"
function θes_upper(f, p)
    θq = f.b0 + f.b1 * (f.pq - f.pa)
    p >= f.pq ? f.b0 + f.b1 * (p - f.pa) : θq + f.B * (f.pq - p) + f.A * (exp((f.pq - p) / f.H) - 1)
end
