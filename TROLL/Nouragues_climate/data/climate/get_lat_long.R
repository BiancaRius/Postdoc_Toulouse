library(nominatimlite)
library(lutz)

local <- "Nouragues, French Guyana"
resultado <- geo_lite_sf(local)

#get latitude and longitude
print(resultado$geometry)

latitude <- 4.0
longitude <- -52.0

#obtaining time zone
time_zone <- tz_lookup_coords(latitude, longitude)
time_zone
