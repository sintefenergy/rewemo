import pandas as pd
from pathlib import Path
import rewemo.solarpower as sol
from rewemo import era5
from yaml import safe_load
import logging
import numpy as np
import xarray as xr

# Modify these:
file_pv_locations = "pv_locations.yaml"
files_era5_data = "C:/Users/hsven/code/energy_timeseries/data/data_europe/era5data_month=2021-*.grib"
output_path = Path("./output_solar_timeseies")

with open(file_pv_locations, "r", encoding="utf8") as f:
    pv_locations = pd.DataFrame.from_dict(safe_load(f), orient="index")
logging.info("Extracting solar radiation data from ERA5 GRIB files (may take some time)")
# ds_ssrd = xr.open_mfdataset(files_era5_data, engine="cfgrib", backend_kwargs={"filter_by_keys": {"shortName": "ssrd"}})
pv_data = era5.extract_solar_radiation_at_locations(files_era5_data, pv_locations)

# Create output folder if it doesn't already exist:
Path(output_path).mkdir(parents=True, exist_ok=True)

for i, pv in pv_locations.iterrows():
    logging.info(f"Location {i}")

    solarp = sol.compute_solar_power(
        df_rad=pv_data[i],
        lat=pv["lat"],
        lon=pv["lon"],
        panel_slope=pv["slope"] * np.pi / 180,
        panel_azimuth=pv["azimuth"] * np.pi / 180,
        albedo=pv["albedo"],
        eta_el=pv["eta_el"],
        tracking=pv["tracking"],
    )
    pv_data[i]["power"] = solarp
    # Save to file
    csv_filename = output_path / f"solarpower_{i}.csv"
    pv_data[i].to_csv(csv_filename)
