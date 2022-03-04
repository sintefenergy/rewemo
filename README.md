# ReWeMo - Renewable energy time series from numerical weather model data


Documentation about the
[ERA5 dataset](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation).


## Solar data

Radiation quantities: https://www.ecmwf.int/sites/default/files/elibrary/2015/18490-radiation-quantities-ecmwf-model-and-mars.pdf


* ssrd - Surface solar radiation downwards (https://apps.ecmwf.int/codes/grib/param-db?id=169)
    - the amount of solar radiation (also known as shortwave radiation) that reaches a horizontal plane at the surface of the Earth. This parameter comprises both direct and diffuse solar radiation.
* fdir - Total sky direct solar radiation at surface (https://apps.ecmwf.int/codes/grib/param-db?id=228021)
    -  the amount of direct solar radiation (also known as shortwave radiation) reaching the surface of the Earth. It is the amount of radiation passing through a horizontal plane, not a plane perpendicular to the direction of the Sun.
* diffuse radiation = ssrd - fdir



Related quantities (not used):
* dsrp - direct solar radiation (perpendicular to sun) https://apps.ecmwf.int/codes/grib/param-db/?id=47
* vddsf - Visible Diffuse Downward Solar Flux (https://apps.ecmwf.int/codes/grib/param-db/?id=260347)
* vbdsf - Visible Beam Downward Solar Flux (W/m2) (https://apps.ecmwf.int/codes/grib/param-db/?id=260346)

## Wind data
* u100 - 100 metre U wind component (https://apps.ecmwf.int/codes/grib/param-db?id=228246)
* v100 - 100 metre V wind component (https://apps.ecmwf.int/codes/grib/param-db?id=228247)
* u10 - 10 metre U wind component (https://apps.ecmwf.int/codes/grib/param-db?id=165)
* v10 - 10 metre V wind component (https://apps.ecmwf.int/codes/grib/param-db?id=166)
* sp - Surface pressure (https://apps.ecmwf.int/codes/grib/param-db?id=134)

Other data
* 2t - 2 metre temperature (https://apps.ecmwf.int/codes/grib/param-db?id=167)
* tp - Total precipitation (https://apps.ecmwf.int/codes/grib/param-db?id=228)


# NCEP/NCAR Reanalysis data
This is found in the ncep_reanalysis package

# MERRA-2 data
Data access: https://disc.gsfc.nasa.gov/datasets?project=MERRA-2
1980-present
Wind: 
 -Read FAQ: https://disc.gsfc.nasa.gov/information/faqs?title=I%20am%20interested%20in%20using%20MERRA-2%20for%20wind%20analysis
 -Data collection: M2T1NXSLV (e.g. u&v wind at 50 m above surface)


# Literature

Comparing ERA-Interim and MERRA: https://doi.org/10.1016/j.renene.2014.09.042
