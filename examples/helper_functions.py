""" Functions

"""

# %% Import
from yaml import safe_load
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# %% Constants, including unit conversion factors
CONV_DECIMAL_TO_PERCENT = 1e2
NO_HOURS_IN_YEAR = 8760
SMALL_NUMBER = 1e-9


# %% Functions

def calculate_indicators_windp(wind_data, power_curve, key_indicators, area):
    
    """ Calculate diagnostic indicator values
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


def sort_data(data, variable_names=['u100', 'v100', 'power']):
            
    """ Sort data for the purpose of duration curves/analyses
    Args:
        data: dictionary with dataframes
        variable_names: names of variables to sort
    
    Returns:
        data_sorted: dictionary with sorted data
                
    Notes:
        Sorting done for all dataframes in dict 'data' and all variable names, 
        keeping track of time steps
    
    """
    
    # Initiate dictionary to hold sorted values
    data_sorted = {}
    
    # For all dataframes ('df') in data input
    for df in data:
        
        data_sorted[df] = {} # Initiate dict
        
        for var in variable_names:
            
            # Determine sorted values (minuses to reverse to descending order)
            values_temp = -np.sort(-data[df][var])
            
            # Determine indexes of sorted values
            index_temp = np.argsort(np.array(-data[df][var]))
            time_steps_sorted = data[df][var].index[tuple([index_temp])]
            
            # Determine sorted values 
            data_sorted[df][var] = pd.DataFrame(values_temp, 
                index=time_steps_sorted, columns=[var])
            del values_temp, index_temp
            
            # Set index name
            data_sorted[df][var].index.name = 'time'
            
    return data_sorted


def plot_duration_curve(data_sorted, area, variable='power'):
    
    """ Plot duration curves with sorted data
    Args:
        data_sorted: dictionary with dataframes with already sorted values
        area: string, name of area
        variable: string, name of variable
        
     """         
    
    fig, ax = plt.subplots()
    for wpp in data_sorted:
        ax.plot(np.array(data_sorted[wpp][variable]), label=wpp)
    ax.set_ylabel('normalized power output [nominal power=1]')
    ax.legend()
    ax.set_title('duration_curve_'+area)
    plt.gcf().set_dpi(300)
    plt.show()
    fig.savefig('duration_curve_'+area+'.png', dpi=300,  bbox_inches='tight')

      
# %%

def compute_wind_power2(df_wind, df_power_curve, cf_target=None):
    """Compute solar power from radiation data for a single panel location
       Wind speed scaling is performed iteratively to match a target value
       for average capacity factor (cf). The scaling follows the approach
       outlined in section 3.2 'Wind speed correction bias' in Staffell and
       Pfenninger (2016), DOI: 10.1016/j.energy.2016.08.068

    df_wind : pandas.DataFrame with UTC timestamp index and these columns:
        u100 : u component wind speed at 100 m height
        v100 : v component wind speed at 100 m height

    power_curve : numpy array with wind speed vs power output

    cf_target : float with target average capacity factor value)

    """
    
    # Wind speed from u and v wind speed components
    wind_speed = np.sqrt(df_wind["u100"] ** 2 + df_wind["v100"] ** 2)
    
    # Wind power
    wind_power = np.interp(wind_speed, df_power_curve.index,
        df_power_curve.values, left=0, right=0)
    
    if cf_target is None:
        print('cf_target is None')
        return (wind_power, None, None)
    
    # Average capacity factor and deviation from target value
    cf_calculate = np.average(wind_power)
    cf_deviation = cf_calculate-cf_target
    
    # Ratio of target capacity factor to calculated capacity factor
    # epsilon, alpha and beta defined following Staffell and Pfenninger (2016)
    epsilon_cf = cf_target/cf_calculate # eq (2) Staffell and Pfenninger (2016)
    alpha = 0.6*epsilon_cf + 0.2        # eq (6) Staffell and Pfenninger (2016)
    beta = 0
        
    # Scale wind speed to match target value for average capacity factor
    ITERATIONS_MAX = 1000
    TOLERANCE = 0.0005
    BETA_INCREMENT = 0.01
    iterations_counter = 0
    if abs(cf_deviation)<=TOLERANCE:
        print('cf_devation: ', cf_deviation)
        return (wind_power, None, None)
    while (abs(cf_deviation)>TOLERANCE and iterations_counter<=ITERATIONS_MAX):
        
        # Update beta
        beta = beta + BETA_INCREMENT*np.sign(-cf_deviation)
        
        # New adjusted wind speed, wind power and capacity factor deviation
        df_wind_new = alpha*df_wind + beta # eq (5) 
        wind_speed_new = np.sqrt(
            df_wind_new["u100"] ** 2 + df_wind_new["v100"] ** 2)
        wind_power_new = np.interp(wind_speed_new, df_power_curve.index,
            df_power_curve.values, left=0, right=0)
        cf_deviation = np.average(wind_power_new)-cf_target
        
        # Counter
        iterations_counter = iterations_counter + 1

    wind_power = wind_power_new
    
    print(iterations_counter, wind_power, alpha, beta)
    return wind_power, alpha, beta


# %%  
if __name__ == '__main__':
    
    """             
    If code executed as program, do several things:
        Load previously calculated wind power outputs from data files
        Group data for individual locations by larger areas
        Calculate key diagnostic indicator values
        Aggregate data for locations to larger areas
    Assume same time steps defined for all locations

    """
           
    # """ Preparations """
    from example_create_wind_power import (KEY_INDICATORS_NAMES,
        file_power_curves)
    areas = ['VESTSYD', 'NORGEMIDT'] # Larger (aggregate) areas considered
    path_root = 'output_wind_timeseries/windpower_' # Data path, root
    wind_data_grouped_by_area = {}
    wind_data_grouped_by_area_and_sorted = {}
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
        
        
        # """ Sort data and plot duration curves """
        wind_data_grouped_by_area_and_sorted[area] = sort_data(
            wind_data_grouped_by_area[area])
        plot_duration_curve(wind_data_grouped_by_area_and_sorted[area], area,
            variable='power')

        

            
            
            
            
            
            
        
            
        
        
        
    
    