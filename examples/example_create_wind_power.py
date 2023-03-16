import pandas as pd
from pathlib import Path
import rewemo.windpower as wp
from rewemo import era5
from yaml import safe_load
import logging
from helper_functions import (
    calculate_indicators_windp, compute_wind_power2)

# Modify these:
areas = ['VESTSYD', 'NORGEMIDT'] # Larger (aggregate) areas considered
file_power_curves = "../src/rewemo/ncep_reanalysis/wind_powercurves_tradewind_with_expansion.csv"
files_wind_data = "../data/data_europe_era5/era5data_month=20*.grib"
output_path = "output_wind_timeseries/"

power_curves = pd.read_csv(file_power_curves, index_col=0)

# Create output folder if it doesn't already exist:
Path(output_path).mkdir(parents=True, exist_ok=True)

file_wpp_locations = []
for area_index, area in enumerate(areas):
    file_wpp_locations.append("wpp_locations_"+area+".yaml")
    
    with open(file_wpp_locations[area_index], "r", encoding="utf8") as f:
        wpp_locations = pd.DataFrame.from_dict(safe_load(f), orient="index")
    
    logging.info("Extracting wind data from GRIB files (may take some time)")
    wind_data = era5.extract_wind_speed_at_locations(
        files_wind_data, wpp_locations)
    
    # Initiate dataframe to hold diagnostic indicator values
    KEY_INDICATORS_NAMES = [
        'capacity factor (%)', 'full load hours (h)',
        'average power', 'median power', 
        'time fraction: zero power (%)',
        'time fraction: zero power low wind (%)',
        'time fraction: zero power high wind (%)',
        'time fraction: max power (%)',
        'min power', 'min power excl. zero power', 'max power', 'stdev power',
        'epsilon_cf', 'alpha', 'beta']
    key_indicators = pd.DataFrame(index=wind_data.keys(),
        columns=KEY_INDICATORS_NAMES)
    
    for i, wpp in wpp_locations.iterrows():
        logging.info(f"Location {i}")
        pcurve_ref = wpp["power_curve"]
        pcurve = power_curves[pcurve_ref]
        
        # If capacity factor target not defined in input data, add assumption
        if "cf_target" not in wpp.index:
            if wpp["power_curve"] == "offshore":
                wpp["cf_target"] = 0.42
            else:
                wpp["cf_target"] = 0.34
        
        # Use either option 1 or option #2 below
        # Option 1: predefined scaling factor for wind speed
        # wind_scaling = wpp["wind_scaling"]
        # windp = wp.compute_wind_power(df_wind=wind_data[i],
        #     df_power_curve=pcurve, wind_scaling=wind_scaling)
        # Option 2: scale to predefined average capacity factor
        (windp, epsilon_cf, alpha, beta) = compute_wind_power2(
            df_wind=wind_data[i], df_power_curve=pcurve,
            cf_target=wpp["cf_target"])
        key_indicators['epsilon_cf'].loc[i] = epsilon_cf
        key_indicators['alpha'].loc[i] = alpha
        key_indicators['beta'].loc[i] = beta
        
        wind_data[i]["power"] = windp
        
        # Calculate diagnostic indicator values
        key_indicators = calculate_indicators_windp(
            wind_data, pcurve, key_indicators, i)
        
        # Save to file
        csv_filename = output_path +"windpower_"+str(i)+".csv"
        wind_data[i].to_csv(csv_filename, sep=';')
    
    # Save to file
    csv_filename = output_path +"key_indicators_"+str(area)+".csv"
    key_indicators.to_csv(csv_filename, sep=';')
