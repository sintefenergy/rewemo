import xarray as xr
import pandas as pd
import numpy as np
import cdsapi
import xarray as xr
from pathlib import Path
from typing import List
import logging


def find_nearest_datapoint(lat, lon, ds):
    """Find the point in the dataset closest to the given latitude and longitude"""
    datapoint_lats = ds.coords.indexes["latitude"]
    datapoint_lons = ds.coords.indexes["longitude"]
    lat_nearest = min(datapoint_lats, key=lambda x: abs(x - lat))
    lon_nearest = min(datapoint_lons, key=lambda x: abs(x - lon))
    return lat_nearest, lon_nearest


def _download_era5_data_single_location(cds, era5_variables, lat, lon, year, file_path):
    """Download data for a specific lat/lon location for a range of years"""
    # It does not seem to work to pick a single lat/lon. Need an area large enough
    # that we pick up two points (or more). Since resolution is 0.25 deg, we use that
    epsilon = 0.25
    params = {
        "product_type": "reanalysis",
        "format": "grib",
        "variable": era5_variables,
        "date": f"{year}-01-01/{year}-12-31",
        "time": [f"{x:02d}:00" for x in range(24)],  # hours 00-23
        "area": [lat + epsilon, lon, lat, lon + epsilon],  # [north, west, south, east]
    }
    cds.retrieve("reanalysis-era5-single-levels", params, file_path)


def _download_era5_atm_grid(cds, file_path):
    """Download data for a single hour to get lat/lon coordinates"""
    params = {
        "product_type": "reanalysis",
        "format": "grib",
        "variable": "100m_u_component_of_wind",  # arbitrary variable
        "date": "2022-01-01",  # arbitrary date
        "time": "00:00",  # arbitrary time
    }
    cds.retrieve("reanalysis-era5-single-levels", params, file_path)


def get_default_variable_list():
    variables = [
        "100m_u_component_of_wind",  # u100
        "100m_v_component_of_wind",  # v100
        "10m_u_component_of_wind",  # "nice to have" for wind scaling?
        "10m_v_component_of_wind",
        "surface_pressure",
        "surface_solar_radiation_downwards",  # ssrd
        "total_sky_direct_solar_radiation_at_surface",  # fdir (diffuse = ssrd - fdir)
        "2m_temperature",
        "total_precipitation",
    ]
    return variables

# TODO: Remove this - inefficient way of downloading data.
def INEFFICIENT_download_era5_data(
    locations: pd.DataFrame, years: List, data_path: Path, replace_existing=False, era5_variables=None
):
    """Download selected ERA5 data to grib files (one for each location)"""

    if era5_variables is None:
        era5_variables = get_default_variable_list()

    cds = cdsapi.Client()

    # make path if it does not exist:
    data_path.mkdir(parents=True, exist_ok=True)

    # Get grid - download entire area for a single hour
    era5_atm_grid_file_name = "era5_atm_grid.grib"
    file_grid = data_path / era5_atm_grid_file_name
    if (not Path(file_grid).is_file()) or replace_existing:
        logging.info(f"Grid file: {file_grid}")
        _download_era5_atm_grid(cds, file_grid)
    ds_grid = xr.open_dataset(file_grid, engine="cfgrib")

    # TODO - check what is more efficient
    # Request single location, or area containing all locations in single request?
    # Request all years, single year or single month per request?
    # 2 locations (x4 poins) and 2 years -> Copernicus running time ca 40 min (+queing) (but because server is busy?)
    for i in locations.index:
        lat_req = locations.loc[i, "lat"]
        lon_req = locations.loc[i, "lon"]
        lat_data, lon_data = find_nearest_datapoint(lat_req, lon_req, ds_grid)
        filename_part1 = f"era5data_lat={lat_data}_lon={lon_data}"
        locations.loc[i, "datafile"] = filename_part1  # excluding the last part (eg. "_year=2021.grib")
        for year in years:
            filename = f"{filename_part1}_year={year}.grib"
            file_data = data_path / filename
            logging.info(f"{i}: Data file: {file_data}")
            if (not Path(file_data).is_file()) or replace_existing:
                logging.info(f"Downloading data to: {file_data}")
                _download_era5_data_single_location(cds, era5_variables, lat_data, lon_data, year, file_data)
    return locations


def _download_era5_data_single_month(cds, era5_variables, area, year, month, file_path):
    params = {
        "product_type": "reanalysis",
        "format": "grib",
        "variable": era5_variables,
        "year": year,
        "month": month,
        "day": [x + 1 for x in range(31)],  # days 1-31
        "time": [f"{x:02d}:00" for x in range(24)],  # hours 00-23
        "area": area,  # [north, west, south, east]
    }
    cds.retrieve("reanalysis-era5-single-levels", params, file_path)


def download_era5_data_area(
    area: List, years: List, months, data_path: Path, replace_existing=False, era5_variables=None
):
    """Download selected ERA5 data to grib files (one for each location)"""

    if era5_variables is None:
        era5_variables = get_default_variable_list()

    cds = cdsapi.Client()

    # make path if it does not exist:
    data_path.mkdir(parents=True, exist_ok=True)

    for year in years:
        for month in months:
            file_data = data_path / f"era5data_month={year}-{month:02d}.grib"
            logging.info(file_data)
            if (not Path(file_data).is_file()) or replace_existing:
                _download_era5_data_single_month(cds, era5_variables, area, year, month, file_data)
    return


def dataarray_to_dataframe(da):
    """Change forecast time and steps to a single hourly time index"""
    if da.step.size == 1:
        df = da.to_pandas()
    if da.step.size > 1:
        df = pd.DataFrame(da.values, index=da.time.values, columns=da.step.values).stack()
        ind2 = df.index.get_level_values(0) + df.index.get_level_values(1)
        df.index = ind2
    return df


def extract_solar_radiation_at_locations(file_pattern, locations):
    """Extract radiation data from era5 grib files for selected locations

    Nearest datapoint is used

    file_pattern : pathlib.Path
        which files to read
    locations : pandas.DataFrame
        columns "lat" and "lon" give locations of panel

    Returns : dictionary of pandas.DataFrames. The keys are the indices in the locations input dataframe
    """
    ds_ssrd = xr.open_mfdataset(file_pattern, engine="cfgrib", backend_kwargs={"filter_by_keys": {"shortName": "ssrd"}})
    ds_fdir = xr.open_mfdataset(file_pattern, engine="cfgrib", backend_kwargs={"filter_by_keys": {"shortName": "fdir"}})
    data_dict = {}
    for i, row in locations.iterrows():
        logging.info(i)
        lat = locations["lat"]
        lon = locations["lon"]
        data_lat, data_lon = find_nearest_datapoint(lat, lon, ds_ssrd)
        da_ssrd = ds_ssrd.sel(latitude=data_lat, longitude=data_lon).ssrd
        da_fdir = ds_fdir.sel(latitude=data_lat, longitude=data_lon).fdir
        df_rad = pd.DataFrame()
        df_rad["fdir"] = dataarray_to_dataframe(da_fdir)
        df_rad["ssdr"] = dataarray_to_dataframe(da_ssrd)
        df_rad["diffuse"] = df_rad["ssdr"] - df_rad["fdir"]
        data_dict[i] = df_rad
    return data_dict


def extract_wind_speed_at_locations(file_pattern, locations):
    """Extract wind speed data from era5 grib files for selected locations

    Nearest datapoint is used

    file_pattern : pathlib.Path
        which files to read
    locations : pandas.DataFrame
        columns "lat" and "lon" give locations of wind power plant

    Returns : dictionary of pandas.DataFrames. The keys are the indices in the locations input dataframe
    """
    ds_u_wind = xr.open_mfdataset(
        file_pattern, engine="cfgrib", backend_kwargs={"filter_by_keys": {"shortName": "100u"}}
    )
    ds_v_wind = xr.open_mfdataset(
        file_pattern, engine="cfgrib", backend_kwargs={"filter_by_keys": {"shortName": "100v"}}
    )
    data_dict = {}
    for i, row in locations.iterrows():
        logging.info(i)
        lat = row["lat"]
        lon = row["lon"]
        data_lat, data_lon = find_nearest_datapoint(lat, lon, ds_u_wind)
        da_u100 = ds_u_wind.sel(latitude=data_lat, longitude=data_lon).u100
        da_v100 = ds_v_wind.sel(latitude=data_lat, longitude=data_lon).v100
        df = pd.DataFrame()
        df["u100"] = dataarray_to_dataframe(da_u100)
        df["v100"] = dataarray_to_dataframe(da_v100)
        data_dict[i] = df
    return data_dict
