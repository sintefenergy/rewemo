# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 10:06:30 2017

@author: hsven
"""

import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import netCDF4
import datetime
import solar as sol
plt.close("all")

#coordinates= {'lat': 63.394, 'lon': 10.415, 'year':[2016]}
datarange = {'lat': [62, 65], 'lon': [9, 12],'year': [2016,2016]}

'''
For the whole geographical range, one year takes about 3 min
I.e. all years take about 3.5 hours
'''
remove_29feb = True
datapath='J:/DOK/12/hsven/reanalysis_data/nc_files'

albedo = 0.2
eta_el = 0.65 #=0.87*0.744
tracker = 'fixed'
slope = 50*pd.np.pi/180 #vertical
azimuth = -24 #degrees south-east


# SCRIPT ---------------------------------------------------

solar = {}

def checkSameVariables(datasets,variables):
    for d in datasets:
        for k,v in variables.items():
            if not list(d.variables[k])==list(v):
                print("ERROR: d={},k={},v={}".format(d,k,v))
                raise Exception("Inconsistent variables ({},{})"
                                .format(d.variables[k],k))

lats = None
temperatures = {}

for year in range(datarange['year'][0],datarange['year'][1]+1):
    print('{}: Year={}'.format(datetime.datetime.now(),year))
    nc_nbdsf = netCDF4.Dataset('{}/nbdsf.sfc.gauss.{}.nc'.format(datapath,year))
    nc_nddsf = netCDF4.Dataset('{}/nddsf.sfc.gauss.{}.nc'.format(datapath,year))
    nc_vbdsf = netCDF4.Dataset('{}/vbdsf.sfc.gauss.{}.nc'.format(datapath,year))
    nc_vddsf = netCDF4.Dataset('{}/vddsf.sfc.gauss.{}.nc'.format(datapath,year))
    nc_temp = netCDF4.Dataset('{}/air.2m.gauss.{}.nc'.format(datapath,year))
    h_nb = nc_nbdsf.variables['nbdsf']
    h_nd = nc_nddsf.variables['nddsf']
    h_vb = nc_vbdsf.variables['vbdsf']
    h_vd = nc_vddsf.variables['vddsf']
    hT = nc_temp.variables['air']
    
    # Do this only once (first year)
    if lats is None:
        lats = nc_nbdsf.variables['lat']
        lons = nc_nbdsf.variables['lon']
    times = nc_nbdsf.variables['time']
    timesT = nc_temp.variables['time']
    #Check that coordinates are the same:
    checkSameVariables([nc_nbdsf,nc_nddsf,nc_vbdsf,nc_vddsf],
                       variables={'lat':lats,'lon':lons,'time':times})
    checkSameVariables([nc_nbdsf,nc_temp],variables={'lat':lats,'lon':lons})

    units = h_nb.units
    jd = netCDF4.num2date(times[:],times.units)
    jdT = netCDF4.num2date(timesT[:],timesT.units)

    lats_selected = [lat for lat in lats 
                     if ((lat>=datarange['lat'][0]) and
                         (lat<=datarange['lat'][-1]))]
    lons_r = ([lon if lon<=180 else lon-360 for lon in lons])
    lons_selected = sorted([lon for lon in lons_r 
                     if ((lon>=datarange['lon'][0]) and 
                         (lon<=datarange['lon'][-1]))])
    coords = [(lat,lon) for lat in lats_selected for lon in lons_selected]
    # Loop thorugh geo range of interest
    for k in coords:
        (lat,lon) = k
        #print("({},{})".format(lat,lon),end=" ")
        ind_lat = list(lats).index(lat)
        ind_lon = list(lons_r).index(lon)

        #hs = pd.Series(h[:,ind_lat,ind_lon],index=jd)
        # Sum visible and near infraread, sum per day W -> Wh
        Hb = 24*(h_nb[:,ind_lat,ind_lon] + h_vb[:,ind_lat,ind_lon])
        Hd = 24*(h_nd[:,ind_lat,ind_lon] + h_vd[:,ind_lat,ind_lon])

        print("{}:({},{})".format(coords.index(k),lat,lon),end=" ")
        
        thisone = pd.DataFrame()
        for i in range(jd.shape[0]):
            t = jd[i]
            nday = t.timetuple().tm_yday
            hours = pd.np.linspace(0,23,24)
            hours_indx = pd.date_range(start=t,periods=24,freq='H')
            rad_tot = sol.H_tiltedpanel(slope=slope,dph=tracker,h=hours,n=nday,
                                        lat=lat,lon=lon,Hb=Hb[i],Hd=Hd[i],
                                        albedo=albedo,azimuth=azimuth)
            solP = sol.pvPower(rad_tot,eta_el)
            

            radiation =  pd.DataFrame({'rad_tot':rad_tot,'power':solP},
                                      index=hours_indx)
            thisone = pd.concat([thisone,radiation])
        temperature = pd.DataFrame(hT[:,ind_lat,ind_lon],index=jdT)        
        
        if k in solar:
            solar[k] = pd.concat([solar[k],thisone],axis=0)
        else:
            solar[k] = thisone

        if k in temperatures:
            temperatures[k] = pd.concat([temperatures[k],temperature],axis=0)
        else:
            temperatures[k] = temperature

#Interpolate temperature (cf how it is done for wind)
for k in temperatures:
    temperatures[k] = temperatures[k].resample('1H').interpolate('linear')
    # Add missing hours at the end 19,20,21,22,23:
    missinghours = 23-temperatures[k].index[-1].hour
    if missinghours!=5:
        raise Exception("Something is wrong")
    missingindx = pd.date_range(start=temperatures[k].index[-1]+1,
                              periods=missinghours,freq='H')
    missingwind = temperatures[k].iloc[[-1]*missinghours].set_index(missingindx)
    temperatures[k] = temperatures[k].append(missingwind)
    solar[k]['T']=temperatures[k]

print('{}: Done'.format(datetime.datetime.now()))
   
print("Exporting results...")
summary = pd.DataFrame(coords,columns=['lat','lon'])
for k in coords:
    print(k,end=" ")
    timeseries_data = solar[k]
    if remove_29feb:
        mask = ((solar[k].index.month==2) & (solar[k].index.day==29))
        timeseries_data = solar[k][~mask]     
    timeseries_data.to_csv('krokstien/solar_({:.3f}, {:.3f}).csv.gz'
                           .format(k[0],k[1]),compression="gzip")
    for c in timeseries_data.keys():
        summary.loc[coords.index(k),c] = solar[k][c].mean()
summary.to_csv('krokstien/solar_SUMMARY.csv')
    

if True:
    df_coords = pd.DataFrame(coords,columns=['lat','lon'])

    plt.figure()
    # projections: cyl, mill
    m = Basemap(projection='mill',resolution='l',
                llcrnrlon=4,llcrnrlat=55,urcrnrlon=20,urcrnrlat=70,
                lat_0=64,lon_0=10)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.25)
    #m.fillcontinents(color='coral',lake_color='aqua')
    m.drawmeridians(pd.np.arange(0,360,10),labels=[True,True,True,True])
    m.drawparallels(pd.np.arange(-90,90,10),labels=[True,True,True,True])
    x,y = m(df_coords['lon'].values,df_coords['lat'].values)
    m.scatter(x,y,marker='o',color='b')
    plt.gcf().set_size_inches(6,6)
    plt.savefig("krokstien/fig_datapoints_surface.png", bbox_inches = 'tight')
    
#    plt.figure()
#    # projections: cyl, mill
#    m = Basemap(projection='mill',resolution='l',
#                llcrnrlon=-11,llcrnrlat=39,urcrnrlon=36,urcrnrlat=76,
#                lat_0=60,lon_0=10)
#    m.drawcoastlines(linewidth=0.25)
#    m.drawcountries(linewidth=0.25)
#    #labels = [left,right,top,bottom]
#    m.drawmeridians(pd.np.arange(0,360,10),labels=[True,False,True,True])
#    m.drawparallels(pd.np.arange(-90,90,10),labels=[True,False,True,True])
#    x,y = m(df_coords['lon'].values,df_coords['lat'].values)
#    nx=len(lons_selected)
#    ny=len(lats_selected)
#    data=[w['power'].mean() for i,w in solar.items()]
#    data2=pd.np.array(data).reshape(ny,nx)
#    clevs = pd.np.arange(0.30,0.40,0.01)
#    cs = m.contourf(x.reshape(ny,nx),y.reshape(ny,nx),data2,cmap=plt.cm.jet,
#                    vmin=0.07, vmax=0.20,levels=pd.np.linspace(0.07,0.20,num=14))
#    m.scatter(x,y,marker='o',edgecolor='k',c=data,cmap=plt.cm.jet,
#              vmin=0.07, vmax=0.20)
#    cbar = m.colorbar(cs)
#    cbar.set_label('Mean capacity factor')
#    plt.gcf().set_size_inches(6,6)
#    plt.savefig("krokstien/fig_solar_capacityfactor_mean.png", bbox_inches = 'tight')

#    plt.figure()
#    # projections: cyl, mill
#    m = Basemap(projection='mill',resolution='l',
#                llcrnrlon=-11,llcrnrlat=39,urcrnrlon=36,urcrnrlat=76,
#                lat_0=60,lon_0=10)
#    m.drawcoastlines(linewidth=0.25)
#    m.drawcountries(linewidth=0.25)
#    #labels = [left,right,top,bottom]
#    m.drawmeridians(pd.np.arange(0,360,10),labels=[True,False,True,True])
#    m.drawparallels(pd.np.arange(-90,90,10),labels=[True,False,True,True])
#    x,y = m(df_coords['lon'].values,df_coords['lat'].values)
#    nx=len(lons_selected)
#    ny=len(lats_selected)
#    data=[(w['T']-273.15).mean() for i,w in solar.items()]
#    data2=pd.np.array(data).reshape(ny,nx)
#    clevs = pd.np.arange(0.30,0.40,0.01)
#    cs = m.contourf(x.reshape(ny,nx),y.reshape(ny,nx),data2,cmap=plt.cm.jet,
#                    vmin=-4, vmax=22,levels=pd.np.linspace(-4,22,num=14))
#    m.scatter(x,y,marker='o',edgecolor='k',c=data,cmap=plt.cm.jet,
#              vmin=-4, vmax=22)
#    cbar = m.colorbar(cs)
#    cbar.set_label('Mean temperature (degC)')
#    plt.gcf().set_size_inches(6,6)
#    plt.savefig("krokstien/fig_solar_temperature_mean.png", bbox_inches = 'tight')

if True:    
    k_plot = list(solar.keys())[0]
    #k_plot = (60,10)
    k_plot = coords[1] # = (59.9986, 11.25) = close to Oslo
    
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(111)
    solar[k_plot]['power'].plot(ax=ax,title='{} at {}'.format('solar',k_plot))
    ax.set_ylabel(units)
    #solar[k_plot]['power'].plot(ax=ax,secondary_y=True)
    plt.ylabel("p.u.")
    plt.gcf().set_size_inches(6,3)
    plt.xlabel("Time")
    plt.savefig("krokstien/fig_solar_power.png", bbox_inches = 'tight')
    
    fig = plt.figure(figsize=(12,4))
    ax = fig.add_subplot(111)
    (solar[k_plot]['T']-273.15).plot(ax=ax,title='{} at {}'.format('Temperature (2m)',k_plot))
    plt.ylabel("degC")
    plt.gcf().set_size_inches(6,3)
    plt.xlabel("Time")
    plt.savefig("krokstien/fig_temperature_2m.png", bbox_inches = 'tight')
    
    fig = plt.figure(figsize=(12,4))
    dd = pd.DataFrame(solar[k]['power'])
    dd['hour']=dd.index.hour
    dd_gr=dd.groupby(by='hour')
    plt.fill_between(dd_gr.mean().index,
                     (dd_gr.mean()-dd_gr.std())['power'].clip(lower=0),
                     (dd_gr.mean()+dd_gr.std())['power'],
                     color='lightgrey',label='$\pm$1 std')
    dd_gr.mean().plot(color='r',ax=plt.gca())
    plt.legend()
    plt.title("Daily profile at {}".format(k))
    plt.gcf().set_size_inches(6,3)
    plt.xlabel("Time")
    plt.ylabel("Mean power")
    plt.savefig("krokstien/fig_dailymean.png", bbox_inches = 'tight')
    
    plt.clf()
    davg=solar[list(solar.keys())[1]].resample("D").mean()
    davg['power'].plot()
    plt.gcf().set_size_inches(6,3)
    plt.xlabel("Time")
    plt.ylabel("Mean power")
    plt.savefig("krokstien/fig_seasonal.png", bbox_inches = 'tight')

    fig = plt.figure(figsize=(10,4))
    df = solar[coords[1]].copy()
    df['hour'] = df.index.hour
    df['date'] = df.index.date
    df2 = df.pivot(index='date',columns='hour',values='power')
    vmin=0.0
    vmax=1.0
    im=plt.imshow(df2.T,vmin=vmin,vmax=vmax,aspect=2,origin="lower")
    plt.gcf().set_size_inches(10,4)
    plt.xlabel("Day")
    plt.ylabel("Hour")
    plt.savefig("krokstien/fig_solar_power_heatmap.png", bbox_inches = 'tight')

if False:
    # Radiation components on PV panel
    Hb = 1#0.01 #example for illustration
    Hd = 0.5#2.5 #example for illustration
    albedo = 0.2
    lat=coords[1][0]
    lon=coords[1][1]
    slope=90*pd.np.pi/180
    day=int(1/2*365)
    #n=np.linspace(0,365,200)
    #d=180/np.pi*decl(n)
    h=pd.np.linspace(0,24,200)
    
    for track in ['south','fixed','azimuth','2-axis']:
        azimuth=0
        if track=="fixed":
            azimuth=-24
        plt.figure()
        plt.title('Radiation on PV panel ({})'.format(track))
        (rad_tot,rad_b,rad_d,rad_r) = sol.H_tiltedpanel_components(
                Hb=Hb,Hd=Hd,albedo=albedo,lat=lat,lon=lon,n_day=day,hour=h,
                tracker=track,slope=slope,azimuth=azimuth)
        plt.plot(h,rad_b,label='direct')
        plt.plot(h,rad_d,label='diffuse')
        plt.plot(h,rad_r,label='reflected')
        plt.plot(h,rad_tot,'--k',label='total',linewidth=2)
        #plt.plot(h,rt,label='rt')
        #plt.plot(h,rd,label='rd')
        plt.xlim([0,24])
        #plt.ylim([0,0.20])
        plt.legend()
        plt.gcf().set_size_inches(5,4)
        plt.savefig('krokstien/fig_pv_rad_{}.png'.format(track), bbox_inches = 'tight')


