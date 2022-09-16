""" Functions

"""

# %% Import
from yaml import safe_load
import numpy as np
import pandas as pd


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
        
    Returns:
        key_indicators: dataframe with diagnostic indicator values, updated
            
    Notes
 
    
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


# %% if 
if __name__ == '__main__':
    
    """             
    If code executed as program, do several things:
        Load previously calculated wind power outputs from data files
        Group data for individual locations by larger areas
        Calculate key diagnostic indicator values
        Aggregate data for locations to larger areas

    """
    
           
    # """ Preparations """
    from example_create_wind_power import (KEY_INDICATORS_NAMES,
        file_power_curves)
    areas = ['VESTSYD', 'NORGEMIDT'] # Larger (aggregate) areas considered
    path_root = 'output_wind_timeseries/windpower_' # Data path, root
    wind_data_grouped_by_area = {}
    key_indicators_grouped_by_area = {}
    power_curves = pd.read_csv(file_power_curves, index_col=0)
    
    # """ Do several tings for all areas """
    for area_index, area in enumerate(areas):
            
        # """ Import data from files and group by area """ 
        file_wpp_locations = "wpp_locations_"+area+".yaml"
        with open(file_wpp_locations, "r", encoding="utf8") as f:
            wpp_locations = pd.DataFrame.from_dict(safe_load(f),
                orient="index")     
        wind_data_grouped_by_area[area] = {}
        for i, wpp in wpp_locations.iterrows():
            wind_data_grouped_by_area[area][i] = pd.read_csv(
                path_root+i+'.csv', sep=',', index_col='time')
            
        # """ Calculate key diagnostic indicator values """
        key_indicators_grouped_by_area[area] = pd.DataFrame(
            index=wind_data_grouped_by_area[area].keys(),
            columns=KEY_INDICATORS_NAMES)
        if area_index==0:  
            # Assume same time steps defined for all locations
            time_steps = wind_data_grouped_by_area[area][
                list(wind_data_grouped_by_area[area].keys())[0]].index
        for i, wpp in wpp_locations.iterrows():
            pcurve = power_curves[wpp["power_curve"]]
            key_indicators_grouped_by_area[area] = (
                calculate_indicators_windp(wind_data_grouped_by_area[area],
                pcurve, key_indicators_grouped_by_area[area], i))
            
        # """ Aggregate data for larger areas """
        if area_index==0:
            wind_data_aggregated_by_area = pd.DataFrame(index=time_steps,
                columns=areas) # Initiate dataframe
        wind_data_power_in_array = np.zeros(
            [len(time_steps), len(wpp_locations)]) # Initiate temporary array
        for index, wpp_name in enumerate(wpp_locations.T):
            wind_data_power_in_array[:,index] = np.array(
                wind_data_grouped_by_area[area][wpp_name]['power'])
        wind_data_aggregated_by_area[area] = np.average(
            wind_data_power_in_array, axis=1)
        del wind_data_power_in_array
        
        

            
            
            
            
            
            
        
            
        
        
        
    
    