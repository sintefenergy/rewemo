from rewemo import era5
import pandas as pd
from pathlib import Path
import logging

output_data_path = Path("../data/data_europe")
area_europe = [70, -10, 35, 30]  # [north, west, south, east]
# years = range(2020, 2022)  # [2020,2021]
# months = range(1,13)  # 1-12
months = list(range(12, 0, -1))
years = list(range(2021, 2016, -1))  # last 5 years
# years = list(range(2021,1978, -1)) # all available years
era5.download_era5_data_area(area_europe, years=years, months=months, data_path=output_data_path)
logging.info("Done")
