import pandas as pd
from pathlib import Path
import rewemo.windpower as wp
from rewemo import era5
from yaml import safe_load
import logging
from helper_functions import (calculate_indicators_windp)

# Modify these:
file_wpp_locations = "wpp_locations.yaml"
file_power_curves = "../src/rewemo/ncep_reanalysis/wind_powercurves_tradewind.csv"
files_wind_data = "../data/data_europe_era5/era5data_month=2021-11.grib"
output_path = Path("./output_wind_timeseries")

with open(file_wpp_locations, "r", encoding="utf8") as f:
    wpp_locations = pd.DataFrame.from_dict(safe_load(f), orient="index")
power_curves = pd.read_csv(file_power_curves, index_col=0)

logging.info("Extracting wind data from ERA5 GRIB files (may take some time)")
wind_data = era5.extract_wind_speed_at_locations(files_wind_data, wpp_locations)

# Initiate dataframe to hold diagnostic indicator values
key_indicators = pd.DataFrame(index=wind_data.keys(), columns = [
    'capacity factor (%)', 'full load hours (h)',
    'average power', 'median power', 
    'time fraction: zero power (%)',
    'time fraction: zero power low wind (%)',
    'time fraction: zero power high wind (%)',
    'time fraction: max power (%)',
    'min power', 'min power excl. zero power', 'max power', 'stdev power'])

# Create output folder if it doesn't already exist:
Path(output_path).mkdir(parents=True, exist_ok=True)

for i, wpp in wpp_locations.iterrows():
    logging.info(f"Location {i}")
    pcurve_ref = wpp["power_curve"]
    pcurve = power_curves[pcurve_ref]
    wind_scaling = wpp["wind_scaling"]

    windp = wp.compute_wind_power(df_wind=wind_data[i], df_power_curve=pcurve, wind_scaling=wind_scaling)
    wind_data[i]["power"] = windp
    
    # Calculate diagnostic indicator values
    key_indicators = calculate_indicators_windp(
        wind_data, pcurve, key_indicators, i)
    
    # Save to file
    csv_filename = output_path / f"windpower_{i}.csv"
    wind_data[i].to_csv(csv_filename)
