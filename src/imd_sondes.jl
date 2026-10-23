# Parser for IMD radiosonde text files ("Standard Summary" exports), from
# read_imd_txt_ASTRAL2025.ipynb. Handles the variants in file names, header
# position, footer lines, column names and missing-value markers.
#
#   include(abspath("imd_sondes.jl")); using .IMDSondes
#   files = soundingfiles("kochi")                          # in ../data/kochi
#   t     = parse_datetime_from_filename.(files)
#   df    = read_sonde(joinpath("../data/kochi", files[1]))
#   p, T, Td, rh, z, wspd, wdir, u, v = get_sounding_cols(df)

module IMDSondes

using CSV, DataFrames, Dates

export imd_stations, imd_latlon, parse_datetime_from_filename, soundingfiles, uniqueidx,
       stripjunk, read_sonde, get_column, get_sounding_cols

"IMD station directory names under ../data"
const imd_stations = ["ahmedabad", "bhubaneswar", "chennai", "karaikal", "kochi", "kolkata", "machilipatnam",
                      "mangalore", "minicoy", "portblair", "pune", "santacruz", "trivandrum", "visakhapatnam"]

"station (lat °N, lon °E), from ../data/IMDstations.kml"
const imd_latlon = Dict(
    "ahmedabad" => (23.0225, 72.5714), "bhubaneswar" => (20.2961, 85.8245), "chennai" => (13.0827, 80.2707),
    "karaikal" => (10.9254, 79.8662), "kochi" => (9.9312, 76.2711), "kolkata" => (22.5726, 88.3639),
    "mangalore" => (12.9141, 74.8560), "minicoy" => (8.3066, 73.0442), "portblair" => (11.6234, 92.7457),
    "pune" => (18.5204, 73.8567), "santacruz" => (19.0902, 72.8508), "visakhapatnam" => (17.6868, 83.2185))

const month_lookup = Dict(lowercase(Dates.monthname(i)[1:3]) => i for i in 1:12)

"Parse dates from filenames using regexs. Returns `nothing` if no pattern matches."
function parse_datetime_from_filename(filename::AbstractString)::Union{DateTime, Nothing}
    patterns = [
        # 0. ddmmyy or ddmmyyyy with HHMM
        r"(?<d>\d{2})(?<m>\d{2})(?<y>\d{2,4})[ _](?<H>\d{2})(?<M>\d{2})",
        # 1. 20230528-12-... -> YYYYMMDD-HH
        r"(?<y>\d{4})(?<m>\d{2})(?<d>\d{2})[-_](?<H>\d{2})",
        # 2. 43353-20250528-1200 -> YYYYMMDD-HHMM
        r"\d{5}[-_](?<y>\d{4})(?<m>\d{2})(?<d>\d{2})[-_](?<H>\d{2})(?<M>\d{2})",
        # 3. 2023-06-27T23 -> ISO
        r"(?<y>\d{4})-(?<m>\d{2})-(?<d>\d{2})[T_](?<H>\d{2})",
        # 4. 01MAY23Chennai or 02May2023...
        r"(?<d>\d{2})(?<mon>[A-Za-z]{3,})?(?<y>\d{2,4})",
    ]
    for pat in patterns
        m = match(pat, filename)
        m === nothing && continue
        try
            year = parse(Int, m["y"])
            year += (year < 100) ? 2000 : 0
            if haskey(m, "mon")   # pattern 4: month name, no hour
                m["mon"] === nothing && continue
                month = get(month_lookup, lowercase(m["mon"][1:3]), 0)
                return DateTime(year, month, parse(Int, m["d"]))
            end
            hour   = haskey(m, "H") ? parse(Int, m["H"]) : 0
            minute = haskey(m, "M") ? parse(Int, m["M"]) : 0
            return DateTime(year, parse(Int, m["m"]), parse(Int, m["d"]), hour, minute)
        catch err
            @warn "Failed to parse matched date" filename err
        end
    end
    return nothing
end

"Standard Summary .txt files in a station directory"
soundingfiles(stationdir; reldatapath="../data/") =
    filter(x -> occursin(r"(?i)(Standard_Summary|Standard-Summary|Standard Summary).*(\.txt)$(?-i)", x),
           readdir(joinpath(reldatapath, stationdir)))

"indices of the first occurrence of each unique value of v"
uniqueidx(v) = unique(i -> v[i], eachindex(v))

"single-character replacer handles tabs and ° degree character °"
stripjunk(f) = IOBuffer(replace(read(f), UInt8('\t') => UInt8(' '), UInt16('°') => UInt8('d')))

"read the file and return a DataFrame (header at the first line starting with \"Time,\")"
function read_sonde(file)
    lines = readlines(file)
    hdr = findfirst(x -> length(strip(x)) > 4 && strip(x)[1:5] == "Time,", lines)
    isnothing(hdr) && (hdr = length(lines) + 1)
    footskip = 0
    if !isempty(lines) && !isempty(lines[end]) && lines[end][1] == '-'
        footskip = 1
        length(lines) > 1 && !isempty(lines[end-1]) && lines[end-1][1] == '-' && (footskip = 2)
    end
    CSV.read(stripjunk(file), DataFrame, delim=",",
             stripwhitespace=true, ignorerepeated=true, header=hdr,
             footerskip=footskip,
             missingstring=["NaN", "NA", "-", "/", ".", "//", "./", ".//", "/.//", "//./", "//-", "///./"])
end

"""
get_column(df, keys) the first column of df whose name is in keys, as numbers
(text entries tryparse'd, unparseable -> NaN); all missing if no key matches.
"""
function get_column(df::DataFrame, keys::Vector{String})
    for k in keys
        if k in names(df)
            c = df[!, k]
            eltype(c) <: Union{Missing, Number} && return c
            return [x isa Number ? x : ismissing(x) ? missing : something(tryparse(Float64, string(x)), NaN) for x in c]
        end
    end
    @warn "No column matching $(keys)"
    return fill(missing, size(df, 1))
end

"""
p, T, Td, rh, z, wspd, wdir, u, v = get_sounding_cols(df)
extract columns by name, insensitive to station and variations in column names.
"""
function get_sounding_cols(df::DataFrame)
    # Row	Time	p(hPa)	T(\xb0C)	U(%)	Dew(\xb0C)	Geo(gpm)	Geo(gpft)	Wdir(\xb0)	Wspd(m/s)	...
    pressure_keys = ["P(hPa)", "p(hPa)"]
    temp_keys     = ["T(C)", "T(dC)", "T(Deg.C)", "T(\xb0C)", "T(°C)"]
    dew_keys      = ["Dew(C)", "Dew(dC)", "Dew(Deg.C)", "Dew(\xb0C)", "Dew(°C)"]
    rh_keys       = ["U(%)"]
    z_keys        = ["Geo(gpm)"]
    wspd_keys     = ["Wspd(m/s)"]
    wdir_keys     = ["Wdir(d)", "Wdir(Deg)", "Wdir(\xb0)", "Wdir(°)"]

    isempty(df) && return ntuple(_ -> Float64[], 9)
    p    = get_column(df, pressure_keys)
    T    = get_column(df, temp_keys)
    Td   = get_column(df, dew_keys)
    rh   = get_column(df, rh_keys)
    z    = get_column(df, z_keys)
    wspd = get_column(df, wspd_keys)
    wdir = get_column(df, wdir_keys)
    u = -wspd .* sin.(pi/180.0*wdir)
    v = -wspd .* cos.(pi/180.0*wdir)
    return p, T, Td, rh, z, wspd, wdir, u, v
end

end # module IMDSondes
