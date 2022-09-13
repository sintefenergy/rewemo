""" Functions

"""

# %% Import
import numpy as np


# %% Constants, including unit conversion factors
CONV_DECIMAL_TO_PERCENT = 1e2
NO_HOURS_IN_YEAR = 8760
SMALL_NUMBER = 1e-9


# %% Functions

def calculate_indicators_windp(wind_data, power_curve, key_indicators, area):
    
    """ 
    Args:
        wind_data: dataframe with u, v and p for given site
        power_curve: array with power curve for given site
        key_indicators: dataframe with diagnostic indicator values
        area: name of site
        
    Returns
    -------
    Dictionary:
        key_indicators: dataframe with diagnostic indicator values, updated
            
    Notes
    -------
    
    """
    
    # Capacity factor, full load hours, average power, median power
    key_indicators['capacity factor (%)'].loc[area] = np.average(
        wind_data[area]['power']*CONV_DECIMAL_TO_PERCENT)
    key_indicators['full load hours (h)'].loc[area] = np.average(
        wind_data[area]['power']*NO_HOURS_IN_YEAR)
    key_indicators['average power'].loc[area] = (np.average(
        (wind_data[area]['power'])))
    key_indicators['median power'].loc[area] = (np.median(
        (wind_data[area]['power'])))
    
    no_time_steps = len(wind_data[area])
    wind_speed = np.sqrt(
        wind_data[area]['u100']**2+wind_data[area]['v100']**2)
        
    # Get first wind speed value from power curve which has power > 0
    # "(i-1)" because of interpolation for winds between i-1 and i
    wind_speed_power_start = (SMALL_NUMBER + power_curve.index\
        [next((i-1) for i, ws in enumerate(power_curve) if ws>0)])
            
    # Get last wind speed value from power curve which has power > 0
    # (data.iloc[::-1] is to reverse order)
    wind_speed_power_cut = (power_curve.iloc[::-1].index\
        [next((i-1) for i, ws in enumerate(power_curve.iloc[::-1]) if ws>0)])
        
    # Time fractions
    key_indicators['time fraction: zero power (%)'].loc[area] = ( 
        np.count_nonzero(np.array(wind_data[area]['power']==0))/ 
        no_time_steps*CONV_DECIMAL_TO_PERCENT)
    key_indicators['time fraction: zero power high wind (%)'].loc[area] = (
        np.count_nonzero(np.array(wind_speed>=wind_speed_power_cut))
        /no_time_steps*CONV_DECIMAL_TO_PERCENT)
    key_indicators['time fraction: max power (%)'].loc[area] = (
        np.count_nonzero(np.array(wind_data[area]['power']==np.max(
        wind_data[area]['power'])))/no_time_steps*CONV_DECIMAL_TO_PERCENT)
    key_indicators['time fraction: zero power low wind (%)'].loc[area] = (
        np.count_nonzero(np.array(wind_speed<wind_speed_power_start))/
        no_time_steps*CONV_DECIMAL_TO_PERCENT)
        
    # Max, min, standard deviation power
    key_indicators['max power'].loc[area] = (np.max(
        np.array((wind_data[area]['power']))))
    key_indicators['min power'].loc[area] = (np.min(
        np.array((wind_data[area]['power']))))
    key_indicators['min power excl. zero power'].loc[area] = (np.min(
        wind_data[area]['power'].loc[wind_data[area]['power']>0]))
    key_indicators['stdev power'].loc[area] = (np.std(
        np.array((wind_data[area]['power']))))
        
    return key_indicators
        