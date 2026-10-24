# tests for sonde_averaging.jl on synthetic soundings:  julia --project=.. test_sonde_averaging.jl
using Test, Dates
include("sonde_averaging.jl"); using .SondeAveraging

const κ = 287.04/1005.7
es(T, p) = SondeAveraging.saturation_vp(T, p)

# θ with shallow cumulus layer to 1500 m, SCI 1500-1700 m (+3 K), moist-adiabatic-ish
# free troposphere to 3000 m, MML 3000-4500 m, MI 4500-4800 m (+3 K)
function θsyn(z)
    z < 1500 && return 299 + 1e-3*z
    z < 1700 && return 299 + 1.5 + 3*(z-1500)/200
    z < 3000 && return 303.5 + 4.5e-3*(z-1700)
    z < 4500 && return 303.5 + 4.5e-3*1300 + 1e-3*(z-3000)
    z < 4800 && return 303.5 + 4.5e-3*1300 + 1.5 + 3*(z-4500)/300
    return 303.5 + 4.5e-3*1300 + 1.5 + 3 + 4.5e-3*(z-4800)
end

"synthetic profile in the all_profiles() format; Td is T - dTd"
function syn(z; dTd=5.0, time=DateTime(2024, 5, 10), lon=72.0)
    p = 1013.25 .* exp.(-z ./ 7500)
    θ = θsyn.(z)
    T = θ .* (p ./ 1000) .^ κ .- 273.15
    rh = 100 .* es.(T .- dTd, p) ./ es.(T, p)
    (source="syn", station="syn", time=time, lat=10.0, lon=lon, p=p, z=collect(float.(z)), θ=θ, rh=rh,
     u=fill(5.0, length(z)), v=fill(0.0, length(z)), z850=1450.0,
     nfree=count(zz -> 1450 <= zz <= 7000, z))
end

zhi = 0.0:10.0:12000.0
hi = syn(zhi)

@testset "T_Td" begin
    T, Td = T_Td(syn(zhi; dTd=7.0))
    @test T .- Td ≈ fill(7.0, length(T)) atol=1e-6
    @test isnan(SondeAveraging.dewpoint(10.0, 1000.0; rh=0.0))
    @test isnan(SondeAveraging.dewpoint(10.0, 1000.0))
    # specific humidity: q at saturation gives Td = T; q from vapor pressure round-trips
    p0, T0 = 900.0, 20.0
    qs = SondeAveraging.specific_humidity(SondeAveraging.saturation_vp(T0, p0), p0)
    @test SondeAveraging.dewpoint(T0, p0; q=qs) ≈ T0 atol=1e-9
    @test SondeAveraging.dewpoint(T0, p0; rh=100.0) ≈ T0 atol=1e-9
    @test SondeAveraging.dewpoint(T0, p0; q=qs/2) < T0
    # Buck (1981): e_w(20 °C, 1000 hPa) = 23.4 hPa, about 14.8 g/kg
    @test SondeAveraging.saturation_vp(20.0, 1000.0) ≈ 23.4 atol=0.1
    @test SondeAveraging.dewpoint(20.0, 1000.0; rh=50.0) ≈ 9.3 atol=0.1
    # rh and the equivalent q agree, and q is a field
    pr = syn(zhi; dTd=7.0)
    @test field_values(pr, :q)[1] > 0
    T1, Td1 = T_Td(pr)
    @test map((t, p, q) -> SondeAveraging.dewpoint(t, p; q=q), T1, pr.p, field_values(pr, :q)) ≈ Td1 atol=1e-6
end

@testset "select_may_arabian_sea" begin
    ps = [syn(zhi; time=DateTime(2024, 5, 1)), syn(zhi; time=DateTime(2024, 6, 1)),
          syn(zhi; time=DateTime(2024, 5, 20), lon=80.0)]
    @test length(select_may_arabian_sea(ps)) == 1
end

@testset "to_pgrid: average high resolution, interpolate low resolution" begin
    pc = 100.0:10.0:1000.0
    # linear in p, 1 hPa samples centered in each bin: bin mean = value at the center
    phi = 95.5:1.0:1004.5
    @test to_pgrid(phi, 2 .* phi .+ 1, pc) ≈ 2 .* pc .+ 1
    # linear in p, mandatory levels every 50 hPa, listed from the surface up: interpolated exactly
    plo = collect(1000.0:-50.0:100.0)
    @test to_pgrid(plo, 2 .* plo .+ 1, pc) ≈ 2 .* pc .+ 1
    # NaN beyond the sampled range and in too-wide gaps
    @test isnan(to_pgrid([500.0, 600.0], [1.0, 2.0], [1000.0])[1])
    @test isnan(to_pgrid([500.0, 900.0], [1.0, 2.0], [700.0]; maxgap=100)[1])
end

@testset "mean_on_pgrid" begin
    pc = 220.0:10.0:1000.0      # the synthetic profiles reach 205 hPa
    m = mean_on_pgrid([hi, hi, syn(zhi; time=DateTime(2024, 5, 2))], :θ; pcenters=pc)
    @test all(m.n .== 3)
    @test all(m.std .< 1e-9)
    # below the SCI (p > 900 hPa), where θ is linear in z: bin mean ≈ θ at the center
    k = findall(pc .>= 910)
    zc = -7500 .* log.(pc[k] ./ 1013.25)
    @test m.mean[k] ≈ θsyn.(zc) atol=0.02
end

@testset "hist_T_Td" begin
    N = 5
    ps = [syn(zhi; dTd=5.0) for _ in 1:N]
    h = hist_T_Td(ps; pcenters=220.0:10.0:1000.0)
    @test size(h.T) == (length(h.Tedges) - 1, 79)
    @test all(sum(h.T; dims=1) .== N)            # one count per sounding per level
    @test maximum(h.T) == N                       # identical soundings fall in one bin
    @test maximum(h.Td) == N
    @test all(sum(h.q; dims=1) .== N) && maximum(h.q) == N
    @test h.qedges[2] - h.qedges[1] ≈ 1e-4
    # q in the histogram bin matches the profile's q at 900 hPa
    k = argmax(h.q[:, 69]); q900 = mean_on_pgrid(ps[1:1], :q; pcenters=220.0:10.0:1000.0).mean[69]
    @test h.qedges[k] <= q900 < h.qedges[k+1]
    # Td is 5 K colder than T
    kT = argmax(h.T[:, 5]); kTd = argmax(h.Td[:, 5])
    @test 4.8 <= (h.Tedges[kT] - h.Tedges[kTd]) <= 5.2
end

@testset "layer_bounds recovers planted layers" begin
    b = layer_bounds(hi)
    @test b !== nothing
    @test isapprox(b.sci_b, 1500; atol=150) && isapprox(b.sci_t, 1700; atol=150)
    @test isapprox(b.ml_b, 3000; atol=300) && isapprox(b.ml_t, 4500; atol=150)
    @test isapprox(b.mi_t, 4800; atol=200)
    @test b.sci_Δθ > 2 && b.mi_Δθ > 2
    # too coarse: mandatory levels only
    zlo = [0.0, 800.0, 1500.0, 3000.0, 5800.0, 9000.0]
    @test layer_bounds(syn(zlo)) === nothing
    # no inversions
    flat = merge(hi, (θ=300 .+ 4.5e-3 .* hi.z,))
    @test layer_bounds(flat) === nothing
end

@testset "composite_layers" begin
    ps = [hi, hi, syn([0.0, 800.0, 3000.0, 9000.0])]
    c = composite_layers(ps, (:θ, :T, :u))
    @test c.nrejected == 1 && c.ok == [1, 2]
    @test length(c.ζ) == sum(DEFAULT_NPER) && all(c.n.θ .== 2)
    # layer boundaries fall at integer ζ; ζ increases monotonically
    @test issorted(c.ζ) && c.ζ[1] > 0 && c.ζ[end] < length(LAYERS)
    # composite of identical soundings is the sounding at the layer heights
    b = layer_bounds(hi)
    @test c.mean.θ ≈ θsyn.(layer_levels(b)) atol=0.05
    @test all(c.mean.u .≈ 5.0)
    @test c.zmean ≈ layer_levels(b)
    @test issorted(c.zmean)
end
