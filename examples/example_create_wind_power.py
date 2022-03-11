import pandas as pd
from pathlib import Path
import rewemo.windpower as wp
from rewemo import era5
from yaml import safe_load
import logging

# Modify these:
file_wpp_locations = "wpp_locations.yaml"
file_power_curves = "ncep_reanalysis/wind_powercurves_tradewind.csv"
files_wind_data = "C:/Users/hsven/code/energy_timeseries/era5/data_europe/era5data_month=2021-*.grib"
output_path = Path("./output_wind_timeseies")

with open(file_wpp_locations, "r", encoding="utf8") as f:
    wpp_locations = pd.DataFrame(safe_load(f)).T
power_curves = pd.read_csv(file_power_curves, index_col=0)

logging.info("Extracting wind data from ERA5 GRIB files (may take some time)")
wind_data = era5.extract_wind_speed_at_locations(files_wind_data, wpp_locations)

# Create output folder if it doesn't already exist:
Path(output_path).mkdir(parents=True, exist_ok=True)

for i, wpp in wpp_locations.iterrows():
    logging.info(f"Location {i}")
    pcurve_ref = wpp["power_curve"]
    pcurve = power_curves[pcurve_ref]
    wind_scaling = wpp["wind_scaling"]

    windp = wp.compute_wind_power(df_wind=wind_data[i], df_power_curve=pcurve, wind_scaling=wind_scaling)
    wind_data[i]["power"] = windp
    # Save to file
    csv_filename = output_path / f"windpower_{i}.csv"
    wind_data[i].to_csv(csv_filename)
