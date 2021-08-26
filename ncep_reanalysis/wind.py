# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 10:06:30 2017

@author: hsven
"""

import pandas as pd
from matplotlib import pyplot as plt
#from mpl_toolkits.basemap import Basemap
#import mpl_toolkits.basemap
import netCDF4
from scipy.interpolate import interp1d
import datetime
import os
#import pyresample
import numpy as np
import timeseries_tools as tt


#datarange = {'lat': [40, 75], 'lon': [-10, 35],'year': [1948,2016]}
#powercurves = pd.read_csv('wind_powercurves_tradewind.csv',
#                          index_col=0)
#datapath='J:/DOK/12/hsven/reanalysis_data/nc_files'

class WindPower:

    def __init__(self,datarange):
        self.datarange = datarange
        self.nc_lats = None
        self.nc_lons = None
        self.lats_selected = None
        self.lons_r = None
        self.lons_selected = None
        self.coords = None
        self.powercurve_refs = None # if specified, list of which powercurve to use


    def computepower(self,wind,powercurve):
        f = interp1d(powercurve.index, powercurve.values, kind='linear',
                     fill_value=(0,0),bounds_error=False)
        power = f(wind)
        return power


    def getAndValidateNcData(self,nc_path,year):
        print('{}: Year={}'.format(datetime.datetime.now(),year))

        datadict = {}
        nc_vwnd = netCDF4.Dataset('{}/vwnd.sig995.{}.nc'.format(nc_path,year))
        hv = nc_vwnd.variables['vwnd']
        nc_uwnd = netCDF4.Dataset('{}/uwnd.sig995.{}.nc'.format(nc_path,year))
        hu = nc_uwnd.variables['uwnd']
        nc_temp = netCDF4.Dataset('{}/air.sig995.{}.nc'.format(nc_path,year))
        hT = nc_temp.variables['air']
        nc_land_sea = netCDF4.Dataset('{}/land.nc'.format(nc_path))

        datadict['hv'] = hv
        datadict['hu'] = hu
        datadict['hT'] = hT
        datadict['nc_land_sea'] = nc_land_sea

        times = nc_vwnd.variables['time']
        if self.nc_lats is None:
            self.nc_lats = nc_vwnd.variables['lat']
            self.nc_lons = nc_vwnd.variables['lon']
        lats = self.nc_lats
        lons = self.nc_lons

        #Check that coordinates are the same:
        if not (nc_uwnd.variables['time'].actual_range
                == nc_vwnd.variables['time'].actual_range).all():
            raise Exception("Different times.")
        if not (nc_temp.variables['time'].actual_range
                == nc_vwnd.variables['time'].actual_range).all():
            raise Exception("Different times.")
        if list(nc_vwnd.variables['lat']) != list(lats):
            raise Exception("Different latitudes.")
        if list(nc_vwnd.variables['lon']) != list(lons):
            raise Exception("Different longitudes.")
        if list(nc_uwnd.variables['lat']) != list(lats):
            raise Exception("Different latitudes.")
        if list(nc_uwnd.variables['lon']) != list(lons):
            raise Exception("Different longitudes.")
        if list(nc_land_sea.variables['lon']) != list(lons):
            raise Exception("Different longitudes.")
        if list(nc_land_sea.variables['lat']) != list(lats):
            raise Exception("Different laitudes.")
        if list(nc_temp.variables['lat']) != list(lats):
            raise Exception("Different latitudes.")
        if list(nc_temp.variables['lon']) != list(lons):
            raise Exception("Different longitudes.")
        times = nc_vwnd.variables['time']

        jd = netCDF4.num2date(times[:],times.units,
            only_use_cftime_datetimes=False,only_use_python_datetimes=True)
        datadict['jd'] = jd

        # if data coordinates have not yet been specified (first run)
        if self.coords is None:
            # if no lat,lon points have been selected (use entire grid)
            if self.lats_selected is None:
                # Use all lat,lon pairs in the range
                self.lats_selected = [lat for lat in lats
                                 if ((lat>=self.datarange['lat'][0]) and
                                     (lat<=self.datarange['lat'][-1]))]
                self.lons_r = ([lon if lon<=180 else lon-360 for lon in lons])
                self.lons_selected = sorted([lon for lon in self.lons_r
                                 if ((lon>=self.datarange['lon'][0]) and
                                     (lon<=self.datarange['lon'][-1]))])
                self.coords = [(lat,lon) for lat in self.lats_selected
                               for lon in self.lons_selected]
            else:
                self.coords = pd.DataFrame({'lat':self.lats_selected,
                                        'lon':self.lons_selected})
            #= [(self.lats_selected[i],self.lons_selected[i])
            #                for i in range(len(self.lats_selected))]

        return datadict

    def makePowerTimeSeries(self,nc_path,outpath,powercurves,
                            windspeed_scaling=1.0,
                            remove_29feb=True,gzip=False):
        '''
        make wind power time series

        Parameters
        ==========
        nc_path : str
            Path where Reanalysis NC files are located
        outpath : str
            Path where generated time series should be placed
        powercurves : pandas dataframe
            Table giving powercurve
        windspeed_scaling : float or dict
            Scaling factor for wind speed, can be single number, or dictionary
            with different values for different (lat,lon) keys
        remove_29feb : bool
            Whether leap year days should be removed
        gzip : bool
            Whether to zip the output csv files
        '''
        wind = {}
        wind_onland = {}

        if not os.path.exists(outpath):
            print("Making path for output data:\n{}".format(outpath))
            os.makedirs(os.path.abspath(outpath))

        datarange = self.datarange

        for year in range(datarange['year'][0],datarange['year'][1]+1):
            datadict = self.getAndValidateNcData(nc_path,year)
            hu = datadict['hu']
            hv = datadict['hv']
            hT = datadict['hT']
            jd = datadict['jd']
            nc_land_sea = datadict['nc_land_sea']
            # Loop thorugh geo range of interest
            for k in self.coords:
                (lat,lon) = k
                ind_lat = list(self.nc_lats).index(lat)
                ind_lon = list(self.lons_r).index(lon)

                #hs = pd.Series(h[:,ind_lat,ind_lon],index=jd)
                thisone = pd.DataFrame({'vwnd':hv[:,ind_lat,ind_lon],
                                        'uwnd':hu[:,ind_lat,ind_lon],
                                        'T':hT[:,ind_lat,ind_lon]},index=jd)
                thisone['v'] = np.sqrt(np.square(thisone['uwnd'])
                                          +np.square(thisone['vwnd']))

                if k in wind:
                    wind[k] = pd.concat([wind[k],thisone],axis=0)
                else:
                    wind[k] = thisone

                if k not in wind_onland:
                    wind_onland[k] = nc_land_sea.variables['land'][
                                                0,ind_lat,ind_lon]

        print('{}: Done'.format(datetime.datetime.now()))

        print("Interploating to hourly wind speed values...")
        for k in wind:
            wind[k] = wind[k].resample('1H').interpolate('linear')
            # Add missing hours at the end 19,20,21,22,23:
            missinghours = 23-wind[k].index[-1].hour
            if missinghours!=5:
                raise Exception("Something is wrong")
            missingindx = pd.date_range(start=wind[k].index[-1]+1,
                                      periods=missinghours,freq='H')
            missingwind = wind[k].iloc[[-1]*missinghours].set_index(missingindx)
            wind[k] = wind[k].append(missingwind)

        # TODO: If necessary, 2d-interpolate to selected coordinates

        print("Computing power and exporting results...")
        summary = pd.DataFrame(self.coords,columns=['lat','lon'])
        self.powercurve_refs = pd.Series(index=wind.keys())
        for k in wind:
            print(k,end=" ")
            if wind_onland[k]:
                pc='avg_lowland_upland'
            else:
                pc='offshore'
            self.powercurve_refs.loc[k] = pc


            if isinstance(windspeed_scaling,dict):
                scaling = windspeed_scaling[k]
            else:
                scaling = windspeed_scaling
            wind[k]['power'] = self.computepower(
                    scaling*wind[k]['v'],powercurves[pc])
            timeseries_data=wind[k]
            if remove_29feb:
                mask = ((wind[k].index.month==2) & (wind[k].index.day==29))
                timeseries_data = wind[k][~mask]

            if gzip:
                timeseries_data.to_csv(
                        '{}/wind_{}_pc{}.csv.gz'.format(outpath,k,pc),
                        compression="gzip")
            else:
                timeseries_data.to_csv(
                        '{}/wind_{}_pc{}.csv'.format(outpath,k,pc))

            for c in timeseries_data.keys():
                summary.loc[self.coords.index(k),c] = timeseries_data[c].mean()
        summary.to_csv('{}/wind_SUMMARY.csv'.format(outpath))

        return (wind,wind_onland)

    def computeInterpolationWeights(self,Ninterpolate=3,
                                    plotWeights=False,outpath=''):
        '''Compute geographical interpolation weights

        parameters
        ----------
        Ninterpolate : int
          number of points to interpolate (N=1: nearest neighbour)
        plotWeights : bool
          whether to create a plot illustrating interpolation
        outpath : string
          path where plot (if any) is saved

        returns
        -------
        df_weight : pandas.DataFrame
          NxM table of interpolation weights. The columns are the M selected
          points of interest, and rows are tne N grid points.
        '''

        if self.nc_lats is None:
            raise Exception("No data Reanalysis has yet been retrieved.")

#        df_weight = pd.DataFrame(columns=self.lats_selected.index)
        print("Determining geographical interpolation weights")
        g_lats = self.nc_lats[:]
        g_lons = [l if l<180 else l-360 for l in self.nc_lons[:]]
        repl = {-90:-89.999, 90:89.999}
        g_lats = [repl[x] if x in repl else x for x in g_lats]
        # this is a list of all (lat,lon) coordinates for which there is data:
        latlon_grid = np.array([[lat,lon] for lat in g_lats for lon in g_lons ])
        latlon_grid = pd.DataFrame(data=latlon_grid,columns=['lat','lon'])
        #latlon_grid['lat'] = g_lats[:]
        #latlon_grid['lon'] = g_lons[:]
        # this is a list of (lat,lon) coordinates of interest:
#        latlon_select = np.array([self.lats_selected,
#                                  self.lons_selected]).transpose()
        latlon_select = pd.DataFrame()
        latlon_select['lat'] = self.lats_selected
        latlon_select['lon'] = self.lons_selected

        latlon_grid = pd.DataFrame(latlon_grid,columns=['lat','lon'])
        latlon_select = pd.DataFrame(latlon_select,columns=['lat','lon'])
        df_weight = tt.computeInterpolationWeights(
                latlon_grid,latlon_select,Ninterpolate)

#        for i in range(self.lats_selected.size):
#            latlon = latlon_select[i,:]
#            dist=pd.Series(_spherical_dist(latlon_grid,latlon))
#            df_weight.loc[:,i] = [0]*latlon_grid.shape[0]
#            if (dist==0).any():
#                # weight one where dist=0, zero elsewhere
#                df_weight.loc[dist==0,i] = 1
#            else:
#                nsmallest = dist.nsmallest(n=Ninterpolate)
#                # weight proportioinal to inverse distance
#                nweight = 1/nsmallest
#                # normalise
#                nweight = nweight/nweight.sum()
#                df_weight.loc[nsmallest.index,i] = nweight
#
        if plotWeights:
            tt.plotWeights(latlon_grid,latlon_select,df_weight,plotlabels=True,
                           filename=outpath+"/interpolate_wind.png")
#            plt.figure()
#            lon_min = latlon_select[:,1].min()-3
#            lon_max = latlon_select[:,1].max()+3
#            lat_min = latlon_select[:,0].min()-3
#            lat_max = latlon_select[:,0].max()+3
#            m = mpl_toolkits.basemap.Basemap(projection='mill',resolution='l',
#                        llcrnrlon=lon_min,llcrnrlat=lat_min,
#                        urcrnrlon=lon_max,urcrnrlat=lat_max,
#                        lat_0=60,lon_0=10)
#            m.drawcoastlines(linewidth=0.25)
#            m.drawcountries(linewidth=0.25)
#            m.drawmeridians(np.arange(0,360,10),labels=[True,True,True,True])
#            m.drawparallels(np.arange(-90,90,10),labels=[True,True,True,True])
#            xg,yg = m(latlon_grid[:,1],latlon_grid[:,0])
#            xp,yp = m(latlon_select[:,1],latlon_select[:,0])
#            for wcol in df_weight.columns:
#                neighbours = df_weight.loc[df_weight[wcol]!=0,wcol]
#                for ni,nv in neighbours.iteritems():
#                    plt.plot([xp[wcol],xg[ni]],[yp[wcol],yg[ni]],
#                             linewidth=nv*10,color='k')
#            m.scatter(xg,yg,marker='o',color='r',zorder=10)
#            m.scatter(xp,yp,marker='o',color='b',zorder=11)
#
#
#            plt.gcf().set_size_inches(6,6)
#            plt.savefig(os.path.join(outpath,"fig_geo_interpol.png"),
#                        bbox_inches = 'tight')
        return df_weight

    def makePowerTimeSeriesSelection(self,nc_path,outpath,powercurves,
                            windspeed_scaling=1.0,
                            remove_29feb=True,gzip=False,plotWeights=False,
                            Ninterpolate=3):
        '''
        make wind power time series for selected lat,lons, using interpolation

        Parameters
        ==========
        nc_path : str
            Path where Reanalysis NC files are located
        outpath : str
            Path where generated time series should be placed
        powercurves : pandas dataframe
            Table giving powercurve
        windspeed_scaling : float or dict
            Scaling factor for wind speed, can be single number, or dictionary
            with different values for different (lat,lon) keys
        remove_29feb : bool
            Whether leap year days should be removed
        gzip : bool
            Whether to zip the output csv files
        Ninterpolate : int
            number of nearest neighbours used for interpolation
        plotWeights : bool
            whether to plot interpolation weights in a figure
        '''
        wind = {}

        if not os.path.exists(outpath):
            print("Making path for output data:\n{}".format(outpath))
            os.makedirs(os.path.abspath(outpath))

        years = self.datarange['year']

        # If not getting windspeed for grid points, compute 2D interpolation
        # factors
        df_weight = None

        for year in range(years[0],years[1]+1):

            datadict = self.getAndValidateNcData(nc_path,year)

            hu = datadict['hu']
            hv = datadict['hv']
            hT = datadict['hT']
            jd = datadict['jd']
            vdim = hu.shape[0]
            # 1464x73x144 (time x lat x lon) -> 1464x10512 (time x latlon_grid)
            hv_2d = hv[:].reshape((vdim,-1),order="C")
            hu_2d = hu[:].reshape((vdim,-1),order="C")
            hT_2d = hT[:].reshape((vdim,-1),order="C")

            # If not done before, calculate geographical interpolation weights:
            if df_weight is None:
                df_weight = self.computeInterpolationWeights(
                        plotWeights=plotWeights,outpath=outpath)
            for i in df_weight.columns:
                thisone = pd.DataFrame(
                        {'vwnd':hv_2d.dot(df_weight[i]),
                         'uwnd':hu_2d.dot(df_weight[i]),
                         'T':hT_2d.dot(df_weight[i])},
                         index=pd.to_datetime(jd))
                thisone['v_calc'] = np.sqrt(np.square(thisone['uwnd'])
                                          +np.square(thisone['vwnd']))
                thisone['v'] = np.sqrt(np.square(hu_2d)
                                          +np.square(hv_2d)).dot(df_weight[i])
                #print("thisone = ",thisone.index)
                if i in wind:
                    wind[i] = pd.concat([wind[i],thisone],axis=0)
                else:
                    wind[i] = thisone

        print('{}: Done'.format(datetime.datetime.now()))

        print("Interploating to hourly wind speed values...")
        for k in wind:
            wind[k] = wind[k].resample('1H').interpolate('linear')
            # Add missing hours at the end 19,20,21,22,23:
            missinghours = 23-wind[k].index[-1].hour
            if missinghours!=5:
                raise Exception("Something is wrong")
            missingindx = pd.date_range(start=wind[k].index[-1]
                                        +datetime.timedelta(hours=1),
                                      periods=missinghours,freq='H')
            missingwind = wind[k].iloc[[-1]*missinghours].set_index(missingindx)
            wind[k] = wind[k].append(missingwind)

        print("Computing power and exporting results...")
        summary = pd.DataFrame(self.coords,columns=['lat','lon'])
        for k in wind:
            print(k,end=" ")
            pc = self.powercurve_refs[k]

            if isinstance(windspeed_scaling,dict):
                scaling = windspeed_scaling[k]
            else:
                scaling = windspeed_scaling
            wind[k]['power'] = self.computepower(
                    scaling*wind[k]['v'],powercurves[pc])
            timeseries_data=wind[k]
            if remove_29feb:
                mask = ((wind[k].index.month==2) & (wind[k].index.day==29))
                timeseries_data = wind[k][~mask]

            if gzip:
                timeseries_data.to_csv(
                        '{}/wind_{}_pc{}.csv.gz'.format(outpath,k,pc),
                        compression="gzip")
            else:
                timeseries_data.to_csv(
                        '{}/wind_{}_pc{}.csv'.format(outpath,k,pc))

            tmp = pd.DataFrame(timeseries_data.mean(),columns=[k]).T
            summary = summary.combine_first(tmp)
        summary.to_csv('{}/wind_SUMMARY.csv'.format(outpath))

        return wind


    def plotTimeseries(self,outpath,windpower,wind_onland=None,k_plot=None):

        df_coords = pd.DataFrame(self.coords,columns=['lat','lon'])
        if wind_onland is not None:
            df_coords['onland'] = [wind_onland[k] for k in self.coords]
        if k_plot is None:
            k_plot = list(windpower.keys())[0]

        fig = plt.figure(figsize=(6,3))
        #ax = fig.add_subplot(111)
        ax1 = windpower[k_plot]['v'].plot(label="Speed",title='{} at {}'.format('wind',k_plot))
        ax2 = windpower[k_plot]['power'].plot(label="Power",ax=ax1,secondary_y=True)
        ax1.set_ylabel('m/s')
        ax2.set_ylabel('p.u.')
        plt.xlabel("Time")
        lines = ax1.get_lines() + ax2.get_lines()
        ax1.legend(lines, [line.get_label() for line in lines])
        #plt.ylabel("MW")
        #plt.gcf().set_size_inches(6,3)
        plt.savefig(os.path.join(outpath,"fig_wind_speedpower.png"),
                    bbox_inches = 'tight')

        fig = plt.figure(figsize=(12,4))
        ax = fig.add_subplot(111)
        windpower[k_plot]['T'].plot(ax=ax,title='{} at {}'.format('Temperature (sig995)',k_plot))
        plt.ylabel("K")
        plt.gcf().set_size_inches(6,3)
        plt.xlabel("Time")
        plt.savefig(os.path.join(outpath,"fig_temperature_sig995.png"),
                    bbox_inches = 'tight')

    def plotPoints(self,outpath,wind_onland=None,labels=True):
        df_coords = pd.DataFrame(self.coords,columns=['lat','lon'])
        if wind_onland is not None:
            df_coords['onland'] = [wind_onland[k] for k in self.coords]
        plt.figure()
        # projections: cyl, mill
        m = mpl_toolkits.basemap.Basemap(projection='mill',resolution='l',
                    llcrnrlon=-11,llcrnrlat=39,urcrnrlon=36,urcrnrlat=76,
                    lat_0=60,lon_0=10)
        m.drawcoastlines(linewidth=0.25)
        m.drawcountries(linewidth=0.25)
        #m.fillcontinents(color='coral',lake_color='aqua')
        m.drawmeridians(np.arange(0,360,10),labels=[True,True,True,True])
        m.drawparallels(np.arange(-90,90,10),labels=[True,True,True,True])
        x,y = m(df_coords['lon'].values,df_coords['lat'].values)
#        if wind_onland is None:
#            colours = 'r'
#        else:
#            colours = ['r' if k==1 else 'b' for k in df_coords['onland']]
        labs, levs = self.powercurve_refs.factorize()
        m.scatter(x,y,marker='o',c=labs,cmap=plt.cm.Set1,vmin=0,vmax=8)
        if labels:
            for i in df_coords.index:
                plt.text(x[i], y[i],i)
        plt.gcf().set_size_inches(6,6)
        plt.savefig(os.path.join(outpath,"fig_datapoints_sig995.png"),
                    bbox_inches = 'tight')


    def plotContour(self,outpath,windpower):
        df_coords = pd.DataFrame(self.coords,columns=['lat','lon'])
        plt.figure()
        # projections: cyl, mill
        m = mpl_toolkits.basemap.Basemap(projection='mill',resolution='l',
                    llcrnrlon=-11,llcrnrlat=39,urcrnrlon=36,urcrnrlat=76,
                    lat_0=60,lon_0=10)
        m.drawcoastlines(linewidth=0.25)
        m.drawcountries(linewidth=0.25)
        #labels = [left,right,top,bottom]
        m.drawmeridians(np.arange(0,360,10),labels=[True,False,True,True])
        m.drawparallels(np.arange(-90,90,10),labels=[True,False,True,True])
        x,y = m(df_coords['lon'].values,df_coords['lat'].values)
        nx=len(self.lons_selected)
        ny=len(self.lats_selected)
        data=[w['power'].mean() for i,w in windpower.items()]
        data2=np.array(data).reshape(ny,nx)
        cs = m.contourf(x.reshape(ny,nx),y.reshape(ny,nx),data2,cmap=plt.cm.viridis)
        m.scatter(x,y,marker='o',edgecolor='k',c=data,cmap=plt.cm.jet)
        cbar = m.colorbar(cs)
        cbar.set_label('Mean capacity factor')
        plt.gcf().set_size_inches(6,6)
        plt.savefig(os.path.join(outpath,"fig_wind_capacityfactor_mean.png"),
                    bbox_inches = 'tight')

def _spherical_dist(pos1, pos2, r=6378.137):
    """Calculates the distance between locations, based on each point's longitude and latitude.
    (Harvesine formula)
    Args:
        pos1: array of 1D arrays containing the longitude and latitude of a point.
            example: [9.1, 45.3] or [[9.1, 45.3], [19.3, 20.3]]
        pos2: array of 1D arrays containing the longitude and latitude of a point.
            example: [9.1, 45.3] or [[9.1, 45.3], [19.3, 20.3]]
        r: Multiplicator to trasform the distance from radians to miles or kilometers.

    Returns:
        The distances between the given points in radians(r=1), miles or kilometers.
    """
    ilat=0
    ilon=1
    pos1 = pos1 * np.pi / 180
    pos2 = pos2 * np.pi / 180
    cos_lat1 = np.cos(pos1[..., ilat])
    cos_lat2 = np.cos(pos2[..., ilat])
    cos_lat_d = np.cos(pos1[..., ilat] - pos2[..., ilat])
    cos_lon_d = np.cos(pos1[..., ilon] - pos2[..., ilon])
    return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))
