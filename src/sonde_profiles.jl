# Readers that put soundings from every source into one profile format:
#   (source, station, time, lat, lon, p, z, θ, rh, u, v, z850, nfree)
# p [hPa], z geopotential height [m], θ [K], rh [%], u, v [m/s],
# z850 = z(850 hPa) [m], the lower limit for inversion and shear diagnostics;
# nfree = number of levels between 850 hPa and 7 km (vertical resolution; old
# Amini soundings may have only mandatory levels and can't resolve inversions).
#
#   include(abspath("sonde_profiles.jl"))
#   profs = all_profiles()            # Revelle 2023, TGT 2024, IMD, Amini
#
# IMD text files are parsed by the IMDSondes module (imd_sondes.jl).

using CSV, DataFrames, Dates, NCDatasets
include(abspath("inversions.jl")); using .Inversions: interp1
include(abspath("imd_sondes.jl")); using .IMDSondes

const κ = 287.04/1005.7          # Rd/cp
θfun(T_C, p) = (T_C + 273.15) * (1000/p)^κ
"""
dθdz_moist(θ, p) saturated (pseudo)adiabatic dθ/dz [K/m] at potential temperature
θ [K] and pressure p [hPa]: (θ/T)(g/cp - Γm), with T = θ(p/1000)^κ and
Γm = g(1 + L qs/(Rd T)) / (cp + L² qs ε/(Rd T²)). About 4.5-5.5 K/km at 3-5 km
in the tropics, so the moist-adiabatic Δθ over 500 m is about 2.5 K.
"""
function dθdz_moist(θ, p)
    g, cp, Rd, L, ε = 9.8, 1005.7, 287.04, 2.501e6, 0.622
    T = θ * (p/1000)^κ
    es = 6.112exp(17.67(T - 273.15)/(T - 29.65))
    qs = ε*es/(p - (1 - ε)*es)
    Γm = g*(1 + L*qs/(Rd*T)) / (cp + L^2*qs*ε/(Rd*T^2))
    (θ/T) * (g/cp - Γm)
end

"Buck (1981) saturation vapor pressure over water [hPa], T in °C, p in hPa"
esat_buck(T_C, p) = (1.0007 + 3.46e-6p) * 6.1121exp(17.502T_C/(240.97 + T_C))

uwind(spd, dir) = -spd * sind(dir)
vwind(spd, dir) = -spd * cosd(dir)
fl(x) = Float64.(coalesce.(x, NaN))

datadir = abspath("../data")

"assemble a profile NamedTuple; drops levels without p, z, or T"
function profile(source, station, time, lat, lon, p, z, T_C, rh, u, v)
    p, z, T_C, rh, u, v = fl.((p, z, T_C, rh, u, v))
    ok = isfinite.(p) .& isfinite.(z) .& isfinite.(T_C) .& (p .> 0)
    p, z, T_C, rh, u, v = (x[ok] for x in (p, z, T_C, rh, u, v))
    z850 = interp1(p, z, [850.0])[1]
    (source=source, station=station, time=time, lat=Float64(lat), lon=Float64(lon),
     p=p, z=z, θ=θfun.(T_C, p), rh=rh, u=u, v=v,
     z850=z850, nfree=count(zz -> z850 <= zz <= 7000, z))
end

# ---------------- R/V Revelle 2023 RAOBE csv ----------------
function read_revelle(file)
    lat = lon = NaN
    df = open(file) do io
        ln = readline(io)
        while !occursin("RAOB/DATA", ln)
            m = match(r"^LAT,([\d.]*),([NS])", ln); m !== nothing && (lat = parse(Float64, m[1]) * (m[2]=="N" ? 1 : -1))
            m = match(r"^LON,([\d.]*),([EW])", ln); m !== nothing && (lon = parse(Float64, m[1]) * (m[2]=="E" ? 1 : -1))
            ln = readline(io)
        end
        CSV.read(io, DataFrame; delim=",", stripwhitespace=true, ignorerepeated=true, missingstring=["-999"])
    end
    t = DateTime(match(r"RAOBE(\d{10})", file)[1], "yyyymmddHH")
    profile("Revelle", "Revelle", t, lat, lon, df.PRES, df.GPM, df.TEMP, df.RH,
            uwind.(fl(df.SPEED), fl(df.WIND)), vwind.(fl(df.SPEED), fl(df.WIND)))
end
revelle_files() = sort(filter(f -> occursin(r"RAOBE\d{10}\.csv$", f),
                              readdir(joinpath(datadir, "Revelle"); join=true)))

# ---------------- R/V Thompson (TGT) 2024 L1 netCDF ----------------
function read_tgt(file)
    NCDataset(file) do ds
        lat = fl(ds["lats"][:]); lon = fl(ds["lons"][:])
        i = findfirst(isfinite, lat)
        profile("TGT", "TGT", ds["time"][1], lat[i], lon[i], ds["pressure"][:], ds["gp_height"][:],
                ds["T"][:], ds["rh"][:], ds["u"][:], ds["v"][:])
    end
end
tgt_files() = sort(filter(f -> occursin(r"TGT_\d{8}_\d{4}_L1\.nc$", f),
                          readdir(joinpath(datadir, "tgt/nc/L1"); join=true)))

# ---------------- Aminidivi (UWyo netCDF, significant levels) ----------------
const amini_latlon = (11.12, 72.73)
function read_amini(file)
    t = DateTime(match(r"(\d{8}_\d{4})", file)[1], "yyyymmdd_HHMM")
    NCDataset(file) do ds
        T = fl(ds["temperature"][:]); Td = fl(ds["dewpoint"][:]); p = fl(ds["pressure"][:])
        profile("Amini", "Amini", t, amini_latlon..., p, ds["height"][:], T,
                100 .* esat_buck.(Td, p) ./ esat_buck.(T, p), ds["u_wind"][:], ds["v_wind"][:])
    end
end
amini_files() = sort(filter(f -> endswith(f, ".nc"), readdir(joinpath(datadir, "uwyo/aminidivi"); join=true)))

# ---------------- IMD text soundings (parser in imd_sondes.jl) ----------------
imd_files(station) = soundingfiles(station; reldatapath=datadir)

function read_imd(station, file)
    df = read_sonde(joinpath(datadir, station, file))
    nrow(df) == 0 && return nothing
    p, T, Td, rh, z, wspd, wdir, u, v = get_sounding_cols(df)
    profile("IMD", station, parse_datetime_from_filename(file), imd_latlon[station]..., p, z, T, rh, u, v)
end

"try a reader; return nothing (and count the failure) if the file can't be read"
function tryread(f, args...; failures)
    try
        r = f(args...)
        (isnothing(r) || count(isfinite, r.θ) < 5) ? (push!(failures, join(string.(args), " ")); nothing) : r
    catch err
        push!(failures, join(string.(args), " ") * " : " * sprint(showerror, err)[1:min(end, 80)])
        nothing
    end
end

"all soundings from every source; returns (profiles, failures)"
function all_profiles(; sources=(:revelle, :tgt, :imd, :amini))
    profs = NamedTuple[]; failures = String[]
    :revelle in sources && append!(profs, filter(!isnothing, [tryread(read_revelle, f; failures) for f in revelle_files()]))
    :tgt     in sources && append!(profs, filter(!isnothing, [tryread(read_tgt, f; failures) for f in tgt_files()]))
    :amini   in sources && append!(profs, filter(!isnothing, [tryread(read_amini, f; failures) for f in amini_files()]))
    if :imd in sources
        for st in sort(collect(keys(imd_latlon)))
            isdir(joinpath(datadir, st)) || continue
            files = imd_files(st)
            # one file per launch time (directories hold duplicates in several formats)
            t = parse_datetime_from_filename.(files)
            keep = unique(i -> t[i], findall(!isnothing, t))
            append!(profs, filter(!isnothing, [tryread(read_imd, st, files[i]; failures) for i in keep]))
        end
    end
    profs, failures
end
