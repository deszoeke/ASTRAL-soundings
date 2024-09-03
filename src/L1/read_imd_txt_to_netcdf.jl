# Code to read and process IMD soundings for ASTraL IOP1
# Step 1: Standardize files
# Input: IMD raw txt files
# Output: Txt and NetCDF files with standard names and formats

# Alex Kinsella, August 2024, modified from code by Simon de Szoeke for ASraL Pilot, June 2023

env_path = "../.."
data_path = "/Users/kinsella/Documents/WHOI/Projects/ASTraL/2024_IOP1/ASTRAL-soundings/data/"

# Activate the environment
using Pkg; Pkg.activate(env_path)

using Revise
using CSV, DataFrames, Dates, NCDatasets

CamelStations = ["Ahmedabad","Bhubaneswar","Chennai","Karaikal","Kochi","Kolkata","Mangalore","Minicoy","PortBlair","Pune","SantaCruz","Visakhapatnam"]

# Function to parse date formats in file names
function imd_file_date(fn)
    if fn[1:6]=="43150-" && all(isnumeric, fn[[7:14; 16:17]]) # 43150 Visakhap. WMO station ID
        dt = Dates.DateTime(fn[7:17], "yyyymmdd-HH")
    elseif fn[10:12]=="VSK" && all(isnumeric, fn[[1:2; 6:9]])  # Visakhap.
        dt = Dates.DateTime(fn[1:9], "ddUUUyyyy")
    elseif all(isnumeric, fn[[1:8; 10:11]])
        dt = Dates.DateTime(fn[1:11], "yyyymmdd-HH")
    elseif all(isnumeric, fn[[1:4; 6:7; 9:10; 12:13; 15:16; 18:19]])
        dt = Dates.DateTime(fn[1:13], "yyyy-mm-dd_HH")
    else
        # expect undelimited 3-letter month
        if all(isnumeric, fn[[1:2; 6:7]])
            if isletter(fn[8])
                dt = DateTime( 2000+year(DateTime(fn[6:7], "yy")),
                    month(DateTime(fn[3:5], "U")),
                    day(DateTime(fn[1:2], "d")) )
            elseif all(isnumeric, fn[8:9])
                dt = DateTime( year(DateTime(fn[6:9], "yyyy")),
                    month(DateTime(fn[3:5], "U")),
                    day(DateTime(fn[1:2], "d")) )
            end
        end
    end
    return dt
end

soundingfiles(stationdir; data_path=data_path) = filter(x -> occursin(r"(?i)(Standard_Summary|Standard-Summary|Standard Summary).*(\.txt)$(?-i)",x),
                   readdir(joinpath(data_path, stationdir)) ) 

stripjunk(f) = IOBuffer(replace(read(f), UInt8('\t') => UInt8(' '), UInt16('\u00b0') => UInt8('d')) )

function extract_flight_start_time(file_path::String)
    open(file_path, "r") do io
        for line in eachline(io)
            if occursin("Flight DT:", line) || occursin("Filght DT:", line)
                # Use a regular expression to capture the datetime string, with or without space
                match = Base.match(r"(?:Flight|Filght) DT:\s*(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})", line)
                if match !== nothing
                    dt_str = match.captures[1]  # Extract the captured datetime string
                    #println("Extracted datetime string: $dt_str")  # Debugging output
                    try
                        return DateTime(dt_str, "yyyy-mm-dd HH:MM:SS")  # Return the parsed datetime
                    catch e
                        error("Could not parse flight start time from: $dt_str")
                    end
                end
            end
        end
    end
end

# Function to read and process sounding files
function read_sonde(file)
    flight_start_time = extract_flight_start_time(file)
    open( file, "r" ) do io
        lines = readlines( io )
        hdr = findfirst(x-> ~isempty(x)&&length(x)>4&&strip(x)[1:5]=="Time,", lines)
        if isnothing(hdr)
            hdr = length(lines)+1
        end
        footskip=0
        if lines[end][1]=='-'
            footskip=1
            if lines[end-1][1]=='-'
                footskip=2
            end
        end
        df = CSV.read( stripjunk( file ),
            DataFrame, delim=",", stripwhitespace=true, ignorerepeated=true, header=hdr[1],footerskip=footskip)
        return df, flight_start_time
    end
end

# Function to convert "MMM.ss:cc" to total seconds
function parse_time_str(time_str::AbstractString)
    time_str = String(time_str)  # Ensure it's a standard String
    parts = split(time_str, r"[:.]")  # Use a regex to split on both ':' and '.'
    minutes = parse(Int, parts[1])
    seconds = parse(Int, parts[2])
    centiseconds = parse(Int, parts[3])
    total_seconds = minutes * 60 + seconds + centiseconds / 1e2
    return total_seconds
end

# Function to extract data columns from DataFrame
function get_sounding_cols(df, station)
    if isempty(df)
        time_col = p = T = Td = rh = wspd = wdir = u = v = NaN*zeros(0)
        ascent = gp_height = msl_height = lat = lon = density = mix_ratio = NaN*zeros(0)
    else
        time_col = parse_time_str.(df[!,"Time"][:])  # Convert time column to total seconds
        if in(station, ["mangalore", "bhubaneswar", "visakhapatnam", "portblair"])
            if in("P(hPa)", names(df))
                p = df[!,"P(hPa)"][:]
                T = df[!,"T(C)"][:]
                Td = df[!,"Dew(C)"][:]
                rh = df[!,"U(%)"][:]
                gp_height = df[!,"Geo(gpm)"][:]
                wspd = df[!,"Wspd(m/s)"][:]
                wdir = df[!,"Wdir(d)"][:]
            elseif in("p(hPa)", names(df)) # slightly different keys
                p = df[!,"p(hPa)"][:]
                rh = df[!,"U(%)"][:]
                gp_height = df[!,"Geo(gpm)"][:]
                wspd = df[!,"Wspd(m/s)"][:]
                if in("T(dC)", names(df))
                    T = df[!,"T(dC)"][:]
                    Td = df[!,"Dew(dC)"][:]
                    wdir = df[!,"Wdir(d)"][:]
                elseif in("T(Deg.C)", names(df)) # another variation
                    T = df[!,"T(Deg.C)"][:]
                    Td = df[!,"Dew(Deg.C)"][:]
                    wdir = df[!,"Wdir(Deg)"][:]
                end
            end
        else # slightly different keys
            p = df[!,"p(hPa)"][:]
            rh = df[!,"U(%)"][:]
            gp_height = df[!,"Geo(gpm)"][:]
            wspd = df[!,"Wspd(m/s)"][:]
            if in("T(dC)", names(df))
                T = df[!,"T(dC)"][:]
                Td = df[!,"Dew(dC)"][:]
                wdir = df[!,"Wdir(d)"][:]
            elseif in("T(Deg.C)", names(df)) # another variation
                T = df[!,"T(Deg.C)"][:]
                Td = df[!,"Dew(Deg.C)"][:]
                wdir = df[!,"Wdir(Deg)"][:]
            end
        end
        # Additional fields that might not be present in all files:
        ascent = in("Asc(m/m)", names(df)) ? df[!,"Asc(m/m)"][:] : NaN*ones(size(p))
        msl_height = NaN*ones(size(p)) # Placeholder for geometric height
        lat = in("Lat", names(df)) ? df[!,"Lat"][:] : NaN*ones(size(p))
        lon = in("Lon", names(df)) ? df[!,"Lon"][:] : NaN*ones(size(p))
        density = NaN*ones(size(p)) # Placeholder for density, not in data
        mix_ratio = NaN*ones(size(p)) # Placeholder for mixing ratio, not in data
    end
    u = -wspd .* sin.(pi/180.0*wdir)
    v = -wspd .* cos.(pi/180.0*wdir)
    
    return time_col, p, T, Td, rh, wspd, wdir, u, v, ascent, gp_height, msl_height, lat, lon, density, mix_ratio
end

# Function to write standardized data to txt and netcdf
function write_standardized_soundings(time_col, flight_start_time, p, T, Td, rh, u, v, ascent, gp_height, msl_height, wspd, wdir, lat, lon, density, mix_ratio, station, nominal_start_time)
    # Convert elapsed time in time_col to UTC time by adding it to start_time
    utc_times = flight_start_time .+ Dates.Second.(time_col)
    
    # Convert the UTC time to seconds since 1970-01-01 00:00:00 UTC
    epoch = DateTime(1970, 1, 1)
    seconds_since_epoch = (utc_times .- epoch) ./ Dates.Second(1)

    # Define output file paths
    fileout = joinpath(data_path, "standardized_sondes/", "$(station)_" * Dates.format(nominal_start_time, "yyyymmdd-HH") * ".csv")
    fileoutnc = joinpath(data_path, "standardized_sondes/", "$(station)_" * Dates.format(nominal_start_time, "yyyymmdd-HH") * ".nc")
    
    # Create DataFrame with units in column headers
    dfo = DataFrame(
        "time (s)" => seconds_since_epoch,
        "p (hPa)" => p,
        "T (°C)" => T,
        "Td (°C)" => Td,
        "rh (%)" => rh,
        "u (m/s)" => u,
        "v (m/s)" => v,
        "ascent (m/min)" => ascent,
        "gp_height (m)" => gp_height,
        "msl_height (m)" => msl_height,
        "lat (degrees_north)" => lat,
        "lon (degrees_east)" => lon,
        "density (g/m³)" => density,
        "mix_ratio (g/kg)" => mix_ratio
    )

    # Write to CSV
    CSV.write(fileout, dfo)
   
    
    # Write to NetCDF
    ds = NCDataset(fileoutnc, "c")
    ds.attrib["nominal_start_time"] = Dates.format(nominal_start_time, "yyyymmdd-HH")
    defDim(ds, "time", length(p))
    
    time_var = defVar(ds, "time", seconds_since_epoch, ("time",))
    gp_var = defVar(ds, "gp_height", gp_height, ("time",))
    msl_var = defVar(ds, "msl", msl_height, ("time",))
    lat_var = defVar(ds, "lats", lat, ("time",))
    lon_var = defVar(ds, "lons", lon, ("time",))
    Pvar = defVar(ds, "P", p, ("time",))
    Tvar = defVar(ds, "T", T, ("time",))
    Tdvar = defVar(ds, "Td", Td, ("time",))
    qvar = defVar(ds, "q", mix_ratio, ("time",))
    rhvar = defVar(ds, "rh", rh, ("time",))
    rhostat = defVar(ds, "rho", density, ("time",))
    ws_var = defVar(ds, "ws", wspd, ("time",))
    wd_var = defVar(ds, "wd", wdir, ("time",))
    uvar = defVar(ds, "u", u, ("time",))
    vvar = defVar(ds, "v", v, ("time",))
    asc_var = defVar(ds, "asc", ascent, ("time",))

    time_var.attrib["units"] = "seconds since 1970-01-01 00:00:00 UTC"
    gp_var.attrib["units"] = "meter"
    msl_var.attrib["units"] = "meter"
    lat_var.attrib["units"] = "degrees_north"
    lon_var.attrib["units"] = "degrees_east"
    Pvar.attrib["units"] = "hPa"
    Tvar.attrib["units"] = "C"
    Tdvar.attrib["units"] = "C"
    qvar.attrib["units"] = "g/kg"
    rhvar.attrib["units"] = "%"
    rhostat.attrib["units"] = "g/m^3"
    ws_var.attrib["units"] = "m/s"
    wd_var.attrib["units"] = "degrees from the North"
    uvar.attrib["units"] = "m/s"
    vvar.attrib["units"] = "m/s"
    asc_var.attrib["units"] = "m/min"
    
    close(ds)
end

# Main loop to process each station's sounding files
stationnames = lowercase.(CamelStations)

for station in stationnames
    println("Processing station: $station")
    files = soundingfiles(station)
    for file in files
        df, flight_start_time = read_sonde(joinpath(data_path, station, file))
        if !isempty(df)
            time_col, p, T, Td, rh, wspd, wdir, u, v, ascent, gp_height, msl_height, lat, lon, density, mix_ratio = get_sounding_cols(df, station)
            nominal_start_time = imd_file_date(file)
            write_standardized_soundings(time_col, flight_start_time, p, T, Td, rh, u, v, ascent, gp_height, msl_height, wspd, wdir, lat, lon, density, mix_ratio, station, nominal_start_time)
        end
    end
end