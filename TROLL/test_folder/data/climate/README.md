# Climate data
Sylvain Schmitt -
Mar 5, 2024

This folder contains climate data used to force the TROLL model at
hourly or half-hourly time steps depending on the product:

- Guyaflux at Paracou provided by Damien Bonal (D. Bonal pers. com.)
  used mainly for validation and soil water content (m3/m3)
- FLUXNET 2015 at Paracou and Tapajos was retrieved freely from FLUXNET
  with the manifest.txt file tracing this retrieve
- ERA5-Land at Paracou was retrieved from the Climate Data Store (CDS)
  using the R package rcontroll
- ERA5-Land at Tapajos was retrieved from
  [Open-Meteo](https://open-meteo.com/) using mamba environment defined
  in get_era.yml and python code defined in get_era.py ran with the code
  below

``` bash
mamba env create -f get_era.yml 
conda activate get-era
python get_era.py
```

``` r
fs::dir_tree()
```

    .
    ├── FLX_BR-Sa1_FLUXNET2015_SUBSET_HR_2002-2011_1-4.csv
    ├── FLX_GF-Guy_FLUXNET2015_SUBSET_HH_2004-2014_2-4.csv
    ├── GX-METEO^MEDDY-2004-2022.csv
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    ├── era5_Paracou.tsv
    ├── era5_tapajos.csv
    ├── get_era.py
    ├── get_era.yml
    └── manifest.txt
