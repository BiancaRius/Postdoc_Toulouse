library(terra) # to read the netCDF files
library(lubridate) # to deal with dates and times
library(dplyr) # to wrangle and tidy the data
library(tidyr) # to wrangle and tidy the data
library(ggplot2) # to make statis maps
library(gganimate) # to make a temporal gif of climate variation
library(rcontroll) # to generate TROLL climate inputs

library(ecmwfr) # to request data from Copernicus
library(osmdata) # to get bounding box from the study area
library(lutz) # to get time zone
library(nominatimlite) # to get coordinates from the study area
library(leaflet) # to make interactive maps
library(sf) # to extract coordinates from spatial objects
library(readr) # to save the data in txt format

wf_set_key(
  user = "6aea1719-7d4d-4e0b-8232-eccd84c72369",
  key = "eb550e20-12ce-4f20-84f7-f34c6c40f6d2"
)

request <- list(
  "dataset_short_name" = "reanalysis-era5-land-monthly-means",
  "format" = "netcdf",
  "product_type" = "monthly_averaged_reanalysis_by_hour_of_day",
  "variable" = c(
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "2m_dewpoint_temperature",
    "2m_temperature",
    "surface_pressure",
    "total_precipitation",
    "surface_solar_radiation_downwards"
  ),
  "month" = sprintf("%02d", 1:12),
  "time" = sprintf("%02d:00", 0:23),
  "year" = as.character(2022),
  "target" = "ERA5land_hr_Nouragues_2022.nc",
  "area" = "3.960414/-52.85468/4.160414/-52.65468"
)

ncfile <- wf_request(
  user = "6aea1719-7d4d-4e0b-8232-eccd84c72369",
  request = request,
  transfer = TRUE,
  path = ".",
  verbose = FALSE
)

request <- list(
  "dataset_short_name" = "reanalysis-era5-land-monthly-means",
  "format" = "netcdf",
  "product_type" = "monthly_averaged_reanalysis", "time" = "00:00",
  "variable" = c(
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "2m_dewpoint_temperature",
    "2m_temperature",
    "surface_pressure",
    "total_precipitation",
    "surface_solar_radiation_downwards"
  ),
  "month" = sprintf("%02d", 1:12),
  "time" = sprintf("%02d:00", 0:23),
  "year" = as.character(2021:2022),
  "target" = "ERA5land_mth_Nouragues_2021_2022.nc",
  "area" = "3.960414/-52.85468/4.160414/-52.65468"
)

#ncfile <- wf_request(
#user = "6aea1719-7d4d-4e0b-8232-eccd84c72369",
#request = request,
#transfer = TRUE,
#path = ".",
#verbose = FALSE
#)

test_r <- suppressWarnings(rast(
  system.file("extdata",
              "ERA5land_hr_Nouragues_2022.nc",
              package = "rcontroll"
  )
))

test <- suppressWarnings(as.data.frame(test_r, xy = TRUE)) %>%
  gather("variable", "value", -x, -y) %>%
  group_by(x, y) %>%
  mutate(date = rep(as_date(terra::time(test_r)))) %>%
  separate(variable, c("variable", "t"), sep = "_(?=\\d)") %>%
  select(-t) %>%
  separate(variable, c("variable", "expver"), sep = "_expver=") %>%
  group_by(x, y, date, variable) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  spread(variable, value) %>%
  arrange(date)


ggplot(test, aes(date, tp)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  xlab("") +
  ylab("Total precipitation")

# get the nc files that you downloaded:

path_era5land_hour = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/Nouragues_climate/data/climate/ERA5land_hr_Nouragues_2022_2.nc"

path_era5land_month = "~/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/Nouragues_climate/data/climate/ERA5land_mth_Nouragues_2021_2022_2.nc"

#climate <- generate_climate(
# x = -52.75468,
# y = 4.060414,
# tz = "America/Cayenne",
# era5land_hour = system.file("extdata", "ERA5land_hr_Nouragues_2022.nc",
#                             package = "rcontroll"
# ),
# era5land_month = system.file("extdata", "ERA5land_mth_Nouragues_2021_2022.nc",
#                              package = "rcontroll"
# )
#)

climate <- generate_climate(
  x = -52.75468,
  y = 4.060414,
  tz = "America/Cayenne",
  era5land_hour = path_era5land_hour,
  era5land_month = path_era5land_month
  )

write_tsv(climate$daytimevar, "/Users/biancarius/Desktop/Postdoc_Toulouse/Postdoc_Toulouse/TROLL/Nouragues_climate/data/climate/ERA5land_daytimevar.txt")
write_tsv(climate$climatedaytime12, "ERA5land_climatedaytime12.txt")


#verifying the files:
library(ncdf4)

# Abra o arquivo e veja suas variáveis
nc_data_local <- nc_open(path_era5land_hour)
print(nc_data_local)

nc_data_r_path = system.file("extdata", "ERA5land_hr_Nouragues_2022.nc",
                                                  package = "rcontroll")
nc_data_r = nc_open(nc_data_r_path)
print(nc_data_r)

# Liste as variáveis disponíveis
variables_local <- names(nc_data_local$var)
print(variables_local)

variables_r <- names(nc_data_r$var)
print(variables_r)

# Feche o arquivo
nc_close(nc_data)