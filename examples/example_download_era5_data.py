from rewemo import era5
from pathlib import Path
import logging

# Data download takes some time, so try first for a smaller time range (e.g. a month or a single year)
# 1 month of data for Europe requires about 300 MB, so think about where to store the data
# Data download requires a CDS api token (register here:)

output_data_path = Path("../data/data_europe")
area_europe = [70, -10, 35, 30]  # [north, west, south, east]
# years = range(2020, 2022)  # [2020,2021]
# months = range(1,13)  # 1-12
months = list(range(12, 0, -1))
years = list(range(2021, 2016, -1))  # last 5 years
# years = list(range(2021,1978, -1)) # all available years
era5.download_era5_data_area(area_europe, years=years, months=months, data_path=output_data_path)
logging.info("Done")
