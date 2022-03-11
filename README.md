# rewemo - Renewable energy time series from numerical weather model data

This Python package includes methods for downloading weater model data from open data repositories, and for creating wind and PV energy time series based on these data.

# 1. Download weather model data

## ERA5 data

NOTE! To download data using these scripts, you need to have a CDS API key. Read more here: https://cds.climate.copernicus.eu/api-how-to

Documentation about the
[ERA5 dataset](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation).


The data downloaded is:
* See the `examples/example_download_era5_data.py` script.

Solar radiation:
* ssrd - Surface solar radiation downwards (https://apps.ecmwf.int/codes/grib/param-db?id=169)
    - the amount of solar radiation (also known as shortwave radiation) that reaches a horizontal plane at the surface of the Earth. This parameter comprises both direct and diffuse solar radiation.
* fdir - Total sky direct solar radiation at surface (https://apps.ecmwf.int/codes/grib/param-db?id=228021)
    -  the amount of direct solar radiation (also known as shortwave radiation) reaching the surface of the Earth. It is the amount of radiation passing through a horizontal plane, not a plane perpendicular to the direction of the Sun.
* diffuse radiation = ssrd - fdir

Wind speed:
* u100 - 100 metre U wind component (https://apps.ecmwf.int/codes/grib/param-db?id=228246)
* v100 - 100 metre V wind component (https://apps.ecmwf.int/codes/grib/param-db?id=228247)

These are not used at presented but also downloaded:

* u10 - 10 metre U wind component (https://apps.ecmwf.int/codes/grib/param-db?id=165)
* v10 - 10 metre V wind component (https://apps.ecmwf.int/codes/grib/param-db?id=166)
* sp - Surface pressure (https://apps.ecmwf.int/codes/grib/param-db?id=134)
* 2t - 2 metre temperature (https://apps.ecmwf.int/codes/grib/param-db?id=167)
* tp - Total precipitation (https://apps.ecmwf.int/codes/grib/param-db?id=228)

Data are downloaded as GRIB files.
Note that data download can be quite slow, depending on traffic on the server. 
- See request queue: https://cds.climate.copernicus.eu/live/queue

Read more about selection of parameters for efficient downloading: https://cds.climate.copernicus.eu/toolbox/doc/how-to/4_how_to_use_output_widgets/4_how_to_use_output_widgets.html#output-ct-output-download

Read more about solar radiation quantities: https://www.ecmwf.int/sites/default/files/elibrary/2015/18490-radiation-quantities-ecmwf-model-and-mars.pdf



### Additional notes about ERA5 data

Related quantities (not used):
* dsrp - direct solar radiation (perpendicular to sun) https://apps.ecmwf.int/codes/grib/param-db/?id=47
* vddsf - Visible Diffuse Downward Solar Flux (https://apps.ecmwf.int/codes/grib/param-db/?id=260347)
* vbdsf - Visible Beam Downward Solar Flux (W/m2) (https://apps.ecmwf.int/codes/grib/param-db/?id=260346)



## NCEP/NCAR Reanalysis data
This is found in the ncep_reanalysis folder

## MERRA-2 data (Not implemented yet)
Data access: https://disc.gsfc.nasa.gov/datasets?project=MERRA-2
1980-present
Wind: 
 -Read FAQ: https://disc.gsfc.nasa.gov/information/faqs?title=I%20am%20interested%20in%20using%20MERRA-2%20for%20wind%20analysis
 -Data collection: M2T1NXSLV (e.g. u&v wind at 50 m above surface)

# 2. Create energy time series

## Using ERA5 data (GRIB files)
1. extract data at selected coordinates using nearest datapoint and store in Pandas data frames
2. convert to energy time series
    * WIND: Use *effective* wind power curve to convert wind speed at hub height to power output
    * PV: Use information about panel location, orientation and tracking system (if any) to convert radiation data to power output

See `examples/example_create_wind_power.py`and `examples/example_create_solar_power.py`for how to do this.

## Using NCEP/NCAR Reanalysis data (NetCDF files):
* See ncep_reanalysis folder for the scripts
* See Memo (http://hdl.handle.net/11250/2468143)



# Literature

Comparing ERA-Interim and MERRA: https://doi.org/10.1016/j.renene.2014.09.042
