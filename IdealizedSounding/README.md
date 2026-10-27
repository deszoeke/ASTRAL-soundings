# IdealizedSounding

An idealized thermodynamic and wind sounding for the eastern Arabian Sea in May, built to set initial conditions for cloud-resolving model (CRM) runs. It is fit to the mean of May soundings over and next to the Arabian Sea (inland stations excluded) in `../src/mean_sonde.ipynb`.

The package is standalone. The fit parameters are copied into `sounding_params`, so it does not depend on the notebook, on Python, or on the sounding data. Its only dependency is ForwardDiff. Plotting is an optional extension that loads with PythonPlot.

## Install

From the `soundings` environment:

```julia
using Pkg
Pkg.develop(path="IdealizedSounding")
```

## Use

```julia
using IdealizedSounding

s = idealized_sounding(1010:-10:10)   # at any pressures [hPa]
v = idealized_vertices()              # at the profile's vertices (significant levels)

s.p, s.θ, s.T, s.Td                   # [hPa], [K], [K], [K]
s.q, s.r, s.rh                        # specific humidity, mixing ratio [kg/kg], relative humidity [0-1]
s.u, s.v                              # [m/s]

using PythonPlot                      # loads the plotting extension
fig = plot_idealized_sounding()       # T/Td, θ, q/rh, u/v against log p, with dots at the vertices
```

Both functions return a NamedTuple of vectors `(; p, θ, T, Td, rh, r, q, u, v)`. Pressures outside 1010–10 hPa give `NaN`.

To change any parameter, pass it as a keyword. Keywords override `sounding_params` in all three functions:

```julia
s = idealized_sounding(p; θtop2=310.5, dd=dd_dry)
```

### Exports

| function | returns |
|---|---|
| `idealized_sounding(p; kw...)` | the sounding at pressures `p` [hPa] |
| `idealized_vertices(; kw...)` | the sounding at its vertices, in descending `p` |
| `plot_idealized_sounding(; p, fontsize=14, kw...)` | a PythonPlot figure; needs `using PythonPlot` |
| `sounding_params` | the default parameters (a NamedTuple) |

`idealized_vertices` gives the levels to hand to a model that interpolates linearly between significant levels. `idealized_sounding` gives the full-resolution θ structure, which is curved between vertices where it follows moist adiabats.

## How the profile is built

Units: pressure in hPa and θ in K.

### Potential temperature

Moist adiabats use Bolton's (1980) saturation equivalent potential temperature θes. They are inverted for T by Newton's method with ForwardDiff derivatives. From the surface up:

| layer | pressures | θ(p) |
|---|---|---|
| lower mixture | 1010–810 | Σ dF₂ max(θbot, θm₂(p)) |
| gap | 810–750 | linear in p between the two mixtures |
| upper mixture | 750–590 | Σ dF₁ max(θbot, θm₁(p)) |
| upper troposphere | 590–400 | on the moist adiabat with θes = θes1 + dθes1dp (p − 590) |
| top | 400–10 | on the moist adiabat with θes = θes(400) + dθes2dp (p − 400) + Aexp (exp((400 − p)/Hexp) − 1) |

**Mixture layers.** Each is an idealized mixture of profiles. Every member follows the moist adiabat θm through θtop at the top of the layer, `ptop`. At a branch level `pv`, a member leaves θm and goes straight down a dry adiabat with θbot = θm(pv). That gives

  θ(p) = Σ dF · max(θbot, θm(p))

with weights dF ≥ 0 and Σ dF = 1. The branch levels `pv1` and `pv2` are levels of the notebook's 20 hPa fit grid. The weights come from a nonnegative least-squares fit to the mean θ, with rms misfit 0.02 K (590–750 hPa) and 0.04 K (810–1010 hPa). The branch levels are the kinks of θ(p).

**Above 590 hPa.** θes is linear in p up to 400 hPa. Above that it is linear plus exponential, fit to the mean θes up to 10 hPa. The rms misfit is 2.9 K in θes and 2.8 K in θ, small compared with θ's range of about 300–800 K.

### Humidity and wind

These are specified at the vertices `ph` and interpolated linearly in log p:
- **Humidity:** the dewpoint depression `dd` = T − Td is read off the mean sounding on a skew-T. It is chosen to give three moisture layers: surface–850, 850–650 and 650–440 hPa.
- **Wind:** `uh` and `vh` are the mean u and v interpolated to `ph`, rounded to 0.1 m/s.

Td = T − dd, where T comes from the full-resolution θ profile. Moisture uses Bolton's saturation vapor pressure, e = es(Td):
- rh = e/es(T)
- r = 0.622 e/(p − e)
- q = 0.622 e/(p − 0.378 e)

### Parameters (`sounding_params`)

| name | default | meaning |
|---|---|---|
| `pbot2`, `ptop2`, `θtop2` | 1010, 810, 309.98 | lower mixture: bottom, top, θ at top |
| `pv2`, `dF2` | 8 levels, weights | lower mixture branch levels and weights |
| `pbot1`, `ptop1`, `θtop1` | 750, 590, 321.52 | upper mixture: bottom, top, θ at top |
| `pv1`, `dF1` | 6 levels, weights | upper mixture branch levels and weights |
| `θes1`, `dθes1dp` | 347.88 K, −0.00803 K/hPa | θes at 590 hPa and its slope up to `pexp` |
| `pexp` | 400 | base of the exponential θes layer |
| `dθes2dp`, `Aexp`, `Hexp` | −0.0242 K/hPa, 0.00098 K, 30 hPa | slope, amplitude and e-folding pressure above `pexp` |
| `ptop` | 10 | top of the profile |
| `ph` | 15 levels | humidity and wind vertices |
| `dd` | K | dewpoint depression T − Td at `ph` |
| `uh`, `vh` | m/s | wind at `ph` |

## Where the parameters come from

The fits are in `../src/mean_sonde.ipynb`. The notebook does not load this package, except to test it.

| parameters | notebook cell (first line) |
|---|---|
| `pv`, `dF`, `θtop`, layer ends | `"nonnegative least squares …"`, the dry-adiabat mixture fit; it prints the `pv`/`dF` lines |
| `θes1` … `Hexp` | `# idealize the mean profile above 590 hPa …`; it prints a parameters line |
| `ph`, `dd` | `# plot dewpoint depression and RH profile` |
| `uh`, `vh` | `# plot wind u,v profiles` |

To update a parameter, rerun its cell and paste the printed values into `sounding_params`.

## Reference

Bolton, D., 1980: The computation of equivalent potential temperature. *Mon. Wea. Rev.*, **108**, 1046–1053.
