import pandas as pd
import ee
import numpy as np

# ee.Authenticate()
ee.Initialize(project='ee-biancafaziorius')
# Initialize Earth Engine
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')

# Read the sites file
sites = pd.read_table("sites_tab.csv", delimiter=";")
sites['site_plot'] = sites['site'] + "_" + sites['plot']



# Select the desired site
site = "Nouragues_plot1"  # Modify as necessary
sites = sites[sites["site_plot"] == site]

# Define the region of interest
longitude, latitude = sites["longitude"].values[0], sites["latitude"].values[0]
point = ee.Geometry.Point([longitude, latitude])

# ERA5-LAND data collection
ic = (
    ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
    .filterBounds(point)
    .filterDate('2001-01-01', '2001-02-28')  # Adjust as necessary
)

# Reduce the collection to the point
def extract_at_point(image):
    values = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point,
        scale=1000,
        maxPixels=1e6
    )
    return ee.Feature(None, values).set({'time': image.date().format()})

# Apply the extraction and create a feature collection
features = ic.map(extract_at_point).getInfo()


# Read the sites file
sites = pd.read_table("sites_tab.csv", delimiter=";")
sites['site_plot'] = sites['site'] + "_" + sites['plot']

# Select the desired site
site = "Nouragues_plot1"  # Modify as necessary
sites = sites[sites["site_plot"] == site]

# Define the region of interest
longitude, latitude = sites["longitude"].values[0], sites["latitude"].values[0]
point = ee.Geometry.Point([longitude, latitude])

# Filter the ERA5-LAND data collection
ic = (
    ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
    .filterBounds(point)
    .filterDate('2001-01-01', '2001-02-28')  # Adjust as necessary
    .select(['dewpoint_temperature_2m','temperature_2m','surface_solar_radiation_downwards','total_precipitation_hourly', 'u_component_of_wind_10m', 'v_component_of_wind_10m','surface_pressure'])  # Select only the desired variables
)

# Extract data at the selected point
def extract_at_point(image):
    values = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point,
        scale=1000,
        maxPixels=1e6
    )
    return ee.Feature(None, values).set({'time': image.date().format()})

# Map and collect the data
features = ic.map(extract_at_point).getInfo()

# Create DataFrame with only the filtered variables
data = []
for feature in features['features']:
    properties = feature['properties']
    properties['time'] = pd.to_datetime(properties['time'])
    data.append(properties)

df = pd.DataFrame(data)

# Calculate wind speed
df['ws'] = np.sqrt(df['u_component_of_wind_10m']**2 + df['v_component_of_wind_10m']**2)

# Add site and plot information
df.insert(0, "site", sites["site"].values[0])
df.insert(0, "plot", sites["plot"].values[0])



# Save to CSV
output_file = "output_daily.csv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Filtered data saved to {output_file}")


#monthly:ic = (
#     ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")
#     .filterBounds(point)
#     .filterDate('2001-01-01', '2001-02-28')  # Adjust as necessary
# )