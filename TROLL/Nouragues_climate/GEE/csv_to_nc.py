import pandas as pd
import xarray as xr
import numpy as np

# Carregar o CSV
csv_file = "nouragues_wind_prec_gee.csv"
df = pd.read_csv(csv_file, delimiter="\t")

# Adicionar coordenadas fixas
df['latitude'] = 4.060414
df['longitude'] = -52.75468

# Certificar-se de que 'time' está no formato datetime
df['time'] = pd.to_datetime(df['time'])

# Criar um pivot table para estruturar os dados como 3D
df_pivot = df.pivot_table(
    index='time', columns=['latitude', 'longitude'], values='ws'
)

# Converter para xarray.Dataset
ds = xr.Dataset(
    {
        "ws": (["time", "latitude", "longitude"], df_pivot.values.reshape(-1, 1, 1))
    },
    coords={
        "time": df_pivot.index,
        "latitude": df_pivot.columns.get_level_values(0).unique(),
        "longitude": df_pivot.columns.get_level_values(1).unique(),
    }
)

# Salvar como arquivo NetCDF
ds.to_netcdf("output.nc")
print("Arquivo NetCDF criado com sucesso!")
