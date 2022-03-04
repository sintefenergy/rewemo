Geographical coverage:
area_europe = [70, -10, 35, 30] #[north, west, south, east]

Data from
https://cds.climate.copernicus.eu/api/v2/resources/reanalysis-era5-single-levels

Variables:
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
