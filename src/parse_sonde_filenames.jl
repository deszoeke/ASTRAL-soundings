using Dates

function parse_datetime_from_filename(filename::String)::Union{DateTime, Nothing}
    patterns = [
        # 1. 20230528-12-... -> YYYYMMDD-HH
        r"(?<y>\d{4})(?<m>\d{2})(?<d>\d{2})[-_](?<H>\d{2})",
        # 2. 43353-20250528-1200 -> YYYYMMDD-HHMM
        r"\d{5}[-_](?<y>\d{4})(?<m>\d{2})(?<d>\d{2})[-_](?<H>\d{2})(?<M>\d{2})",
        # 3. 2023-06-27T23 -> ISO
        r"(?<y>\d{4})-(?<m>\d{2})-(?<d>\d{2})[T_](?<H>\d{2})",
        # 4. 01MAY23Chennai or 02May2023...
        r"(?<d>\d{2})(?<mon>[A-Za-z]{3,})?(?<y>\d{2,4})",
    ]

    # Map 3-letter month abbreviations to numbers
    month_lookup = Dict{String, Int}()
    for (i, m) in enumerate( map(x -> lowercase(x[1:3]), Dates.monthname.(1:12)) )
        month_lookup[lowercase(m[1:3])] = i
    end

    for pat in patterns
        m = match(pat, filename)
        if m !== nothing
            try
                if haskey(m, "mon")
                    # Handle '01MAY23' style
                    mon = lowercase(m["mon"])
                    month = get(month_lookup, mon, 0)
                    year = parse(Int, m["y"])
                    year += (year < 100) ? 2000 : 0
                    return DateTime(year, month, parse(Int, m["d"]))
                else
                    year  = parse(Int, m["y"])
                    month = parse(Int, m["m"])
                    day   = parse(Int, m["d"])
                    hour  = haskey(m, "H") ? parse(Int, m["H"]) : 0
                    minute = haskey(m, "M") ? parse(Int, m["M"]) : 0
                    return DateTime(year, month, day, hour, minute)
                end
            catch err
                @warn "Failed to parse matched date" filename err
            end
        end
    end

    return nothing  # No match found
end


# test
using Printf

for filename in readlines("file_variety.txt")
    dt = parse_datetime_from_filename(filename)
    @printf("%s, %s\n", filename, Dates.format(dt, dateformat"yyyy-mm-dd HH"))
end
