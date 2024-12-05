import pandas as pd
import ee
import numpy as np

# ee.Authenticate()
ee.Initialize(project='ee-biancafaziorius')
# Inicializar Earth Engine
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')

# Ler o arquivo de sites
sites = pd.read_table("sites_tab.csv", delimiter=";")
sites['site_plot'] = sites['site'] + "_" + sites['plot']

# Selecionar o site desejado
site = "Nouragues_plot1"  # Modifique conforme necessário
sites = sites[sites["site_plot"] == site]

# Definir a região de interesse
longitude, latitude = sites["longitude"].values[0], sites["latitude"].values[0]
point = ee.Geometry.Point([longitude, latitude])

# Coleção de dados ERA5-LAND
ic = (
    ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
    .filterBounds(point)
    .filterDate('2001-01-01', '2001-02-28')  # Ajuste conforme necessário
)

# Reduzir a coleção ao ponto
def extract_at_point(image):
    values = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point,
        scale=1000,
        maxPixels=1e6
    )
    return ee.Feature(None, values).set({'time': image.date().format()})

# Aplicar a extração e criar uma coleção de features
features = ic.map(extract_at_point).getInfo()


# Ler o arquivo de sites
sites = pd.read_table("sites_tab.csv", delimiter=";")
sites['site_plot'] = sites['site'] + "_" + sites['plot']

# Selecionar o site desejado
site = "Nouragues_plot1"  # Modifique conforme necessário
sites = sites[sites["site_plot"] == site]

# Definir a região de interesse
longitude, latitude = sites["longitude"].values[0], sites["latitude"].values[0]
point = ee.Geometry.Point([longitude, latitude])


# Filter the ERA5-LAND data collection
ic = (
    ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
    .filterBounds(point)
    .filterDate('2001-01-01', '2001-02-28')  # Modifique conforme necessário
    .select(['total_precipitation_hourly', 'u_component_of_wind_10m', 'v_component_of_wind_10m'])  # Seleciona apenas as variáveis desejadas
)

# Extrair dados no ponto selecionado
def extract_at_point(image):
    values = image.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=point,
        scale=1000,
        maxPixels=1e6
    )
    return ee.Feature(None, values).set({'time': image.date().format()})

# Mapear e coletar os dados
features = ic.map(extract_at_point).getInfo()

# Criar DataFrame com apenas as variáveis filtradas
data = []
for feature in features['features']:
    properties = feature['properties']
    properties['time'] = pd.to_datetime(properties['time'])
    data.append(properties)

df = pd.DataFrame(data)

# Calcular velocidade do vento
df['ws'] = np.sqrt(df['u_component_of_wind_10m']**2 + df['v_component_of_wind_10m']**2)

# Adicionar informações do site e plot
df.insert(0, "site", sites["site"].values[0])
df.insert(0, "plot", sites["plot"].values[0])

# Salvar em CSV
output_file = "output_filtered_data.csv"
df.to_csv(output_file, sep="\t", index=False)

print(f"Dados filtrados salvos em {output_file}")
