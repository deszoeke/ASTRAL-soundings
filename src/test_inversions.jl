# tests for inversions.jl on synthetic profiles:  julia --project=.. test_inversions.jl
using Test, Random
include("inversions.jl"); using .Inversions

Random.seed!(1)
lat = 12.0
f = 2*7.292e-5*sind(lat)

# background θ with 4.5 K/km, mixed layer 3000-4800 m, 8 K stable cap 4800-5200 m
function θsyn(z; zml=3000.0, zb=4800.0, zt=5200.0, Δθ=8.0)
    θ0 = 300 + 4.5e-3*z
    z < zml && return θ0
    z < zb  && return 300 + 4.5e-3*zml
    z < zt  && return 300 + 4.5e-3*zml + Δθ*(z-zb)/(zt-zb)
    return 300 + 4.5e-3*zml + Δθ + 4.5e-3*(z-zt)
end
z = collect(20.0:10.0:12000.0)
θ = θsyn.(z)

@testset "interp1" begin
    @test interp1([0, 1, 2], [0, 10, 20], [0.5, 1.5]) ≈ [5, 15]
    @test isnan(interp1([0, 1], [0, 1], [2.0])[1])
    @test interp1([2, 0, 1, NaN], [20, 0, 10, 5], [0.5]) ≈ [5]       # unsorted, NaN dropped
    @test isnan(interp1([0, 10], [0, 1], [5.0]; maxgap=5)[1])         # gap too wide
end

@testset "clean capping inversion" begin
    invs = find_inversions(z, θ)
    @test length(invs) == 1
    i = invs[1]
    @test abs(i.zb - 4800) <= 100 && abs(i.zt - 5200) <= 100
    @test 7.5 <= i.Δθ <= 9.5
    @test i.capping && abs(i.mld - 1800) <= 200
end

@testset "noise robustness" begin
    θn = θ .+ 0.1 .* randn(length(z))    # 0.1 K sensor noise at 10 m
    invs = find_inversions(z, θn)
    @test count(i -> i.Δθ >= 2, invs) == 1
    i = primary_inversion(invs)
    @test !isnothing(i) && abs(i.zb - 4800) <= 150 && i.capping
end

@testset "no inversion" begin
    @test isempty(find_inversions(z, 300 .+ 4.5e-3 .* z))
    θn = 300 .+ 4.5e-3 .* z .+ 0.1 .* randn(length(z))
    @test all(i -> i.Δθ < 2, find_inversions(z, θn))
end

@testset "sparse significant levels" begin
    zs = [20, 1000, 3000, 4000, 4800, 5200, 7000, 9000, 12000.0]
    invs = find_inversions(zs, θsyn.(zs))
    @test length(invs) == 1
    @test abs(invs[1].zb - 4800) <= 100 && invs[1].capping
end

@testset "winds, Ri, Margules tilt signs" begin
    # interface rising northward 1 km/1000 km and eastward 0.5 km/1000 km
    gp = 9.8*8/((300 + 4.5e-3*3000) + 4)
    Δu =  gp/f * 1e-3
    Δv = -gp/f * 0.5e-3
    u = [zz < 5000 ? -8.0 : -8.0 + Δu for zz in z]
    v = [zz < 5000 ?  2.0 :  2.0 + Δv for zz in z]
    d = inversion_diagnostics(z, θ, u, v; lat=lat)
    @test length(d) == 1
    @test d[1].Δu ≈ Δu rtol=0.01
    @test d[1].dhdy ≈ 1e-3 rtol=0.1          # θ0 and Δθ from 300 m windows, so ~5-10% off
    @test d[1].dhdx ≈ 0.5e-3 rtol=0.1
    @test d[1].Ri > 1
    @test margules_tilt(1.0, 0.0, 1.0, 300.0, lat).dhdy > 0     # westerly jump ⇒ deeper to north
    @test margules_tilt(0.0, 1.0, 1.0, 300.0, lat).dhdx < 0     # southerly jump ⇒ deeper to west
    # shear in the mixed layer
    u2 = -1e-3 .* z
    d2 = inversion_diagnostics(z, θ, u2, 0 .* z; lat=lat)
    @test d2[1].dudz_ml ≈ -1e-3 rtol=1e-6
end

@testset "shear excluded below zshear" begin
    # mixed layer 1000-1800 m capped at 1800 m; strong frictional shear below 1000 m
    θb = θsyn.(z; zml=1000.0, zb=1800.0, zt=2200.0)
    u = [zz < 1000 ? 10 - 0.01*(1000 - zz) : 10 - 1e-3*(zz - 1000) for zz in z]
    d = inversion_diagnostics(z, θb, u, 0 .* z; lat=lat, zshear=1000.0, δ=1000.0)
    i = d[findfirst(x -> abs(x.zb - 1800) <= 100, d)]
    @test i.dudz_ml ≈ -1e-3 rtol=0.05          # free-tropospheric shear only
    @test i.u_b ≈ 10 - 1e-3*(1800 - 1000 - 400) atol=0.2   # window clipped to [1000, 1800]
    @test isnan(jumps(z, u, 900.0, 1100.0; zmin=1000.0).below)
end

# elevated mixed layer 3000-4600 m between a base inversion (2600-3000 m) and a cap (4600-5000 m)
function θml(z; Δbase=3.0, Δtop=8.0)
    z < 2600 && return 300 + 4.5e-3*z
    θ1 = 300 + 4.5e-3*2600
    z < 3000 && return θ1 + Δbase*(z-2600)/400
    z < 4600 && return θ1 + Δbase
    z < 5000 && return θ1 + Δbase + Δtop*(z-4600)/400
    θ1 + Δbase + Δtop + 4.5e-3*(z-5000)
end

@testset "mixed layer and bounding inversions" begin
    for (Δb, Δt) in ((4.0, 8.0), (6.0, 3.6))       # strong cap; strong base, weak cap (≥ 9 K/km)
        θm = θml.(z; Δbase=Δb, Δtop=Δt) .+ 0.05 .* randn(length(z))
        mls = find_mixed_layers(z, θm)
        m = mls[findfirst(x -> 2500 < x.ml_t < 5500, mls)]
        @test abs(m.ml_b - 3000) <= 150 && abs(m.ml_t - 4600) <= 150
        @test abs(m.top_t - 5000) <= 150 && abs(m.base_b - 2600) <= 150
        @test m.top_Δθ ≈ Δt atol=0.6
        @test m.base_Δθ ≈ Δb atol=0.6
        @test m.top_Δθ500 ≈ Δt + 4.5e-3*100 atol=0.6   # 400 m cap + 100 m above
        @test m.base_Δθ500 ≈ Δb + 4.5e-3*100 atol=0.6
    end
    # no base inversion when the layer below is just moist-adiabatic
    θn = [zz < 3000 ? 300 + 4.5e-3*zz : θml(zz) - θml(2999.0) + 300 + 4.5e-3*3000 for zz in z]
    m = find_mixed_layers(z, θn)[1]
    @test isnan(m.base_b) && m.base_Δθ == 0 && m.top_Δθ > 6
    @test isempty(find_mixed_layers(z, 300 .+ 4.5e-3 .* z))
end

@testset "relative composite" begin
    zr, m, n = relative_composite([z, z], [θ, θ], [4800.0, NaN]; zrel=-500.0:50.0:500.0)
    @test all(n .== 1)
    @test m[zr .== 0][1] ≈ θsyn(4800.0)
end
