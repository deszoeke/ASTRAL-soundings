# Plan for CRM cases

Analysis of Arabian Sea soundings in May are in mean_sonde.ipynb.

Here we specify a CRM case reference density (potential temperature) and humidity structure.
The reference is mainly defined by its stability relative to the moist adiabat (labeled by theta_e), d(theta-dtheta_e)/dp or dtheta_es/dp. Moist convection sets the stability of the tropical troposphere.

Above pRC \approx 300 hPa, the profile is an idealized radiative convective equilibrium mean observed profile. Below 300-400 hPa, the profile is slightly more stable than the moist adiabat.

**The parameters for p>=pRC are:**

theta_es and theta_e are have a minimum at p1 approx 590 hPa. Below this minimum (p>p1), theta and theta_es increase downward by mixing in a distribution of dry adiabats (increasing downward faster than the moist adiabat). 
    - The stratification approaches moist adiabatic at a tangent point at about 850 hPa.

Between p1 and p2,
theta_es increases with height slightly faster than the moist adiabat.

Between p2 and pRC, theta_es follows a moist adiabat at theta_es = theta_w2.

**To do:**
Plan a function that models idealized soundings according the above specifications.

First list the parameters of the function here.

## Function: `idealized_sounding` → (p, θ, T, q, z)

```julia
idealized_sounding(p; ps=1010.0, pt=850.0, p1=550.0, p2=430.0, pRC=380.0,   # hPa (p1, p2: θes extrema of the mean)
                      θes1=347.79, θes2=349.5,                               # moist adiabats (θes labels) at p1, p2 [K]
                      θbot=missing, dF=missing,                              # dry-adiabat mixture below p1 (from fit_dry_mixture)
                      upper=:observed)                                       # profile above pRC
```
(p1, p2: local θes min and max of the May Arabian Sea mean sounding, mean_sonde.ipynb.)

Pressure ordering: surface ps = 1010 > tangent pt ≈ 850 > p1 = 550 > p2 = 430 > pRC = 380 hPa.
θm(p; θes) denotes θ along the moist adiabat labeled by θes (θw = `t0_moist(θes)` at 1000 hPa), from metpy `moist_lapse`.

**Parameters**

| name | meaning | first guess / source |
|---|---|---|
| `ps` | surface pressure | 1010 hPa |
| `p1` | θes local minimum | 550 hPa |
| `θes1` | moist adiabat through the θes minimum at p1 | 347.79 K |
| `pt` | tangent point: dry-adiabat mixture becomes moist adiabatic | 850 hPa |
| `θs_dist` | distribution of surface θ (dry adiabats) mixed below p1; mean = E[max(θs, θm(p; θes1))] | Gaussian μ ≈ 302 K, σ ≈ 2.75 K, or the isotonic F from mean_sonde |
| `p2` | θes local maximum; top of the slightly-stable layer | 430 hPa |
| `θes2` | moist adiabat followed from p2 to pRC; θes rises from θes1 at p1 to θes2 at p2 | 349.5 K |
| `pRC` | base of the RCE profile | 380 hPa |
| `upper` | profile above pRC: observed mean θ (May AS) or an idealized RCE profile | observed mean |
| humidity | q (or RH) per layer: below p1 from θe of the mixture; θe minimum value at p1; RH above p1 | TBD |

**Derived, not free:** θes between p1 and p2 (linear in p between the two end values), T from θ, z hydrostatic.
**Checks:** pt follows from `θs_dist` and θes1 (warmest member meets θm near pt), so either fix `pt` or the distribution's upper tail, not both. Continuity of θ at p1, p2, pRC.

**Open questions:** humidity specification per layer; which CRM (input format, winds).
