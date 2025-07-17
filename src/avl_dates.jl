using Pkg; Pkg.activate("..")
using Dates
using HTTP

# get available Aminidivi (station no 43311) soundings
#
# urls that return valid soundings for each year
site = "https://weather.uwyo.edu"
urls = [
"/wsgi/sounding?datetime=2025-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2024-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2022-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2021-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2020-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2019-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2018-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2017-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2016-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2015-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2014-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2013-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2012-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2011-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2010-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2009-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2008-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2007-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2006-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2005-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2004-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2003-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2002-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2001-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=2000-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1999-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1998-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1997-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1996-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1995-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1994-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1993-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1992-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1991-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1990-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1989-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1988-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1987-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1986-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1985-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1984-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1983-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1982-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1981-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1980-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1979-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1978-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1977-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1976-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1975-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1974-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
"/wsgi/sounding?datetime=1973-01-01 00:00:00&id=43311&type=INVENTORY&src=FM35"
]

# h = "https://weather.uwyo.edu/wsgi/sounding?datetime=2025-01-01%2000:00:00&id=43311&type=INVENTORY&src=FM35"
url_encode(u) = replace(u, " " => "%20")

# regex matches date parts of valid soundings
re = r"<TD><a href=\"\/wsgi\/sounding\?src=FM35&datetime=(\d{4})-(\d{2})-(\d{2}) (\d{2}):(\d{2}):(\d{2})"

dts = DateTime[] # initial empty
for u in urls
      html = String(HTTP.get(site*url_encode(u)).body)
      for m in eachmatch(re, html)
	   dt = DateTime(parse.(Int, m.captures)...)
	   println(dt)
	   push!(dts, dt)
      end
end

open("dates.txt", "w") do io
      [print(io, dt,"\n") for dt in sort(dts)]
end
