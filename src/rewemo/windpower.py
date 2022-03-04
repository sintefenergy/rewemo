import numpy as np


def compute_wind_power(df_wind, df_power_curve, wind_scaling=1.0):
    """Compute solar power from radiation data for a single panel location

    df_wind : pandas.DataFrame with UTC timestamp index and these columns:
        u100 : u component wind speed at 100 m height
        v100 : v component wind speed at 100 m height

    power_curve : numpy array with wind speed vs power output

    wind_scaling : float - wind scaling factor (e.g. to adjust for height)

    """

    wind_speed = np.sqrt(df_wind["u100"] ** 2 + df_wind["v100"] ** 2)
    wind_speed = wind_speed * wind_scaling
    wind_power = np.interp(wind_speed, df_power_curve.index, df_power_curve.values, left=0, right=0)
    return wind_power
