# Script to download nc files from CDS - model ERA5-LAND
# Script based on the API request section at (end of the webpage) https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-monthly-means?tab=download
# Additionally, some rationale and user configuration from "/rcontroll-main/vignettes/climate.Rmd"
# ATTENTION: the file ~/.cdsapirc needs to be well configured. See instructions at: https://cds.climate.copernicus.eu/how-to-api

import cdsapi

# You need two types of product: (1) monthly averages by hour of day and (2) monthly averages at 00:00 (Schmitt, climate.Rmd)


dataset = "reanalysis-era5-land"
request = {
    "variable": [
        "2m_dewpoint_temperature"
    ],
    "year": "2004",
    "month": "01",
    "day": [
        "01"
    ],
    "time": [
        "00:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [4, -52, 4, -52]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()

# product = input("which product? type 1 for monthly_averaged_reanalysis_by_hour_of_day or 2 for monthly averages at 00:00  ")


# if product == "1":

# # Request (1): Month means for each hour of day for 12 months 
# # Each month contains 24 means (one mean for each hour). For the 12 months in a year, you will have 
# # 24×12=288 mean values.
#     request = {
#         "product_type": ["monthly_averaged_reanalysis_by_hour_of_day"],
#         "variable": ["10m_u_component_of_wind",
#             "10m_v_component_of_wind",
#             "2m_dewpoint_temperature",
#             "2m_temperature",
#             "surface_pressure",
#             "total_precipitation",
#             "surface_solar_radiation_downwards"],
#         "year": ["2022"],
#         "month": ["01", "02", "03","04", "05", "06","07", "08", "09","10", "11", "12"],
#         "time": ["00:00", "01:00", "02:00",
#                 "03:00", "04:00", "05:00",
#                 "06:00", "07:00", "08:00",
#                 "09:00", "10:00", "11:00",
#                 "12:00", "13:00", "14:00",
#                 "15:00", "16:00", "17:00",
#                 "18:00", "19:00", "20:00",
#                 "21:00", "22:00", "23:00"],
#         "data_format": "netcdf",
#         "download_format": "unarchived",
#         "area": [3.960414,-52.85468,4.160414,-52.65468]
#     }
   
#     client = cdsapi.Client()
#     client.retrieve(dataset, request).download("ERA5land_hr_Nouragues_2022_2.nc")

# else:
# # Request (2): Month means for the variable value in an specific time 00:00 UTC
# # Each month will have only one value
#     request = {
#         "product_type": ["monthly_averaged_reanalysis"],
#         "variable": ["10m_u_component_of_wind",
#             "10m_v_component_of_wind",
#             "2m_dewpoint_temperature",
#             "2m_temperature",
#             "surface_pressure",
#             "total_precipitation",
#             "surface_solar_radiation_downwards"],
#         "year": ["2021","2022"],
#         "month": ["01", "02", "03","04", "05", "06","07", "08", "09","10", "11", "12"],
#         "time": ["00:00"],
#         "data_format": "netcdf",
#         "download_format": "unarchived",
#         "area": [3.960414,-52.85468,4.160414,-52.65468]
#     }
    
#     client = cdsapi.Client()
#     client.retrieve(dataset, request).download("ERA5land_mth_Nouragues_2021_2022_2.nc")

