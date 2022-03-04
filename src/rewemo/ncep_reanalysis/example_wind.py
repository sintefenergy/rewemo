# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 10:06:30 2017

@author: hsven
"""

import pandas as pd
import wind
import matplotlib.pyplot as plt

datarange = {'lat': [55, 60], 'lon': [2.5, 8],'year': [2010,2016]}

powercurves = pd.read_csv('wind_powercurves_tradewind.csv',
                          index_col=0)

datapath='J:/DOK/12/hsven/reanalysis_data/nc_files'


# scaling factor - 40m to 90m height, roughness
# power law:
scaling1 = (90/40)**0.143 # = 1.12

#https://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XLII-4-W1/217/2016/isprs-archives-XLII-4-W1-217-2016.pdf
z0=0.0002 #surface roughness in m
d=0 #zero plane displacement
scaling2 = pd.np.log((90-d)/z0)/pd.np.log((40-d)/z0) # = 1.06

myWind = wind.WindPower(datarange)

(wind,wind_onland) = myWind.makePowerTimeSeries(nc_path=datapath,
                                   outpath='output',
                                   powercurves=powercurves,
                                   windspeed_scaling = 1.3,
                                   remove_29feb=True)

myWind.makePlots(wind,wind_onland,outpath="output",k_plot=(57.5,5))
