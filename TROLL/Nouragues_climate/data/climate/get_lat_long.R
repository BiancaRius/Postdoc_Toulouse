library(nominatimlite)
library(lutz)

local <- "Nouragues, French Guyana"
resultado <- geo_lite_sf(local)

#get latitude and longitude
print(resultado$geometry)

latitude <- 4.060414
longitude <- -52.75468

#obtaining time zone
time_zone <- tz_lookup_coords(latitude, longitude)
time_zone
