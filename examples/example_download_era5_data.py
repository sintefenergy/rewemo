from energy_timeseries import era5
import pandas as pd
from pathlib import Path
import logging

output_data_path = Path("./data_europe")
area_europe = [70, -10, 35, 30]  # [north, west, south, east]
years = range(2020, 2022)  # [2020,2021]
months = range(1,13)  # 1-12
era5.download_era5_data_area(area_europe, years=years, months=months, data_path=output_data_path)
logging.info("Done")

