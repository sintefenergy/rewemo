import numpy as np


def decl(n):
    """Declination of the sun at day number n

    Parameters
    ----------
    n : int
        day of year (0-365)
    Returns
    -------
    Declination of the sun (radians)
    """
    return np.pi / 180 * 23.45 * np.sin(2 * np.pi * (284 + n) / 365)


def hourangle(h, lon):
    """Compute the hourangle of the sun at a given time

    Parameters
    ----------
    h : float
        hour of day (universal time, 0-24)
    lon : float
        longitude (degrees), east is positive, west is negative
    Returns
    -------
    hour angle in radians in the range (-pi,pi)
    """
    w = np.pi / 180 * (15 * (h - 12) + lon)
    # shift to (-pi,pi) region
    w2 = np.array(w)
    mask1 = w2 > np.pi  # too large
    mask2 = w2 < -np.pi  # too small
    w2[mask1] = w2[mask1] - 2 * np.pi
    w2[mask2] = w2[mask2] + 2 * np.pi
    return w2


def hourangle_sunset(n, lat):
    """Compute the hourangle where sun sets (rad)

    Parameters
    ----------
    n : int
        day of year (0-356)
    lat : float
        latitude (degrees)
    Returns
    -------
    hour angle in radians
    """
    cosws = -np.tan(lat * np.pi / 180) * np.tan(decl(n))
    hourangle = np.arccos(cosws)
    # special cases:
    hourangle[cosws > 1] = -1  # no sunset, summer light
    hourangle[cosws < -1] = np.pi
    return hourangle


def zenithangle(h, n, lat, lon):
    """Compute the angle between zenith and the sun

    Parameters
    ----------
    n : int
        day of year (0-356)
    h : float
        hour of the day (0-24)
    lat,lon : float
        latitude,longitude (degrees)
    Returns
    -------
    zenith angle of the sun in radians (0-pi)
    """
    lat = lat * np.pi / 180
    z = np.arccos(np.sin(lat) * np.sin(decl(n)) + np.cos(lat) * np.cos(decl(n)) * np.cos(hourangle(h, lon)))
    return z


def cos_incidenceangle(slope, delta_phi, theta_z):
    """
    Cosine of the incidence angle between sun and PV panel

    Parameters
    ----------
    slope : float
        slope of PV panel (radians)
    delta_phi : float
        azimuth angle difference between panel and sun
    theta_z : float
        solar zenith angle

    Returns
    -------
    cos(theta) : float
        cos of incidence angle
    """
    costh = np.cos(theta_z) * np.cos(slope) + np.sin(theta_z) * np.sin(slope) * np.cos(delta_phi)
    costh[costh > 1] = 1  # may occur due to approximations?
    costh[costh < 0] = 0  # costh<0 may occur close to sunrise/sunset
    return costh


def panel_angles(tracker, slope, azimuth, lat, lon, hour_of_day, day_of_year):
    """Computes sun-panel azimuth angle difference, and panel slope"""
    if (tracker is None) or (tracker == "fixed"):
        # no tracking
        delta_phi = hourangle(hour_of_day, lon - azimuth)
    elif tracker == "azimuth":
        # panel azimuth angle follows the sun, slope is fixed
        delta_phi = 0
    elif tracker == "2-axis":
        # panel is always facing the sun, slope equals sun angle from zenith
        delta_phi = 0
        slope = zenithangle(hour_of_day, day_of_year, lat, lon)
    else:
        raise ValueError(f"Tracker must be '2-axis','azimuth','fixed' (None). Got {tracker}")
    return {"slope": slope, "delta_phi": delta_phi}


def compute_solar_power(df_rad, lat, lon, panel_slope, panel_azimuth, albedo, eta_el, tracking=None):
    """Compute solar power from radiation data for a single panel location

    df_rad : pandas.DataFrame with UTC timestamp index and these columns:
        ssrd : surface solar radiation downward (J/m2s) onto horizontal sufrace
        fdir : surface direct radiation (beam) onto horizontal surface

    lat, lon : latitude and longitude location of panel

    panel_slope : float (radians)
        panel angle from horizontal position (horizonal: slope=0, vertical: slope = pi/2)

    panel_azimuth : float (radians)
        panel orientation, north=0, east=pi/2, south=pi

    albedo : float
        ground albedo

    eta_el : float
        efficiency of electrical system from panel to grid connection point

    tracking : str
        panel tracking type (None, "azimuth", "2-axis")

    Returns
        pandas.Series - panel power output (relative to installed capacity)


    """

    # Maximum ratio of beam radion on tilted panel vs horizontal surface
    # (to deal with bad approximation when zenith angle is close to 90 degrees)
    # (used in Rb calculation)
    MAX_BEAM_PANEL_RATIO = 10
    MIN_BEAM_PANEL_RATIO = 0

    diffuse_horizontal = df_rad["ssdr"] - df_rad["fdir"]
    beam_horizontal = df_rad["fdir"]

    day_of_year = np.array(df_rad.index.dayofyear)
    hour_of_day = np.array(df_rad.index.hour)

    # 1. compute angle between panel and the sun
    angles = panel_angles(tracking, panel_slope, panel_azimuth, lat, lon, hour_of_day, day_of_year)
    slope = angles["slope"]
    delta_phi = angles["delta_phi"]
    # Compute R_b factor
    theta_z = zenithangle(hour_of_day, day_of_year, lat, lon)
    cos_theta = cos_incidenceangle(slope, delta_phi, theta_z)
    rb_factor = cos_theta / np.cos(theta_z)

    # Correct for poor approximation close to sunrise/sunset at high latitudes:
    # TODO: Double check that this is reasonable
    solar_hour_angle_sunset = hourangle_sunset(day_of_year, lat)
    mask1 = solar_hour_angle_sunset < 30 * np.pi / 180  # sun sets very early
    mask2 = theta_z > 85 / 180 * np.pi  # solar zenith is low on horizon
    rb_factor[mask1 & mask2] = 0
    # Make sure the rb factor does not blow up - this is also a sunrise/sunset issue
    # TODO: Double check that this is reasonable
    rb_factor = np.clip(rb_factor, a_min=MIN_BEAM_PANEL_RATIO, a_max=MAX_BEAM_PANEL_RATIO)

    cosb = np.cos(slope)
    rad_direct = rb_factor * beam_horizontal
    rad_diffuse = (1 + cosb) / 2 * diffuse_horizontal
    rad_reflect = albedo * (1 - cosb) / 2 * (beam_horizontal + diffuse_horizontal)
    rad_on_tilted_surface = rad_direct + rad_diffuse + rad_reflect

    # Radiation is accumulated J/m**2 at over one hour. Convert to W/m2:
    rad_on_tilted_surface_w_per_m2 = rad_on_tilted_surface / 3600

    # Def panel capacity: output under standard test conditions: 1000 W/m2 irradiation
    # Power output per installed capacity (so electrical system efficiency is relevant,
    # but not panel conversion efficiency which is instead reflected in panel power capacity)
    solar_power = rad_on_tilted_surface_w_per_m2 * eta_el / 1000

    return solar_power
