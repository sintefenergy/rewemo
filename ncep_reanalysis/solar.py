"""
2017-12-19
@author: Harald G Svendsen, SINTEF Energy Research

This is a module to compute hourly solar radiation on tilted surfaces
from daily average beam and direct solar radiation (obtained e.g. from
Reanalysis datasets)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4
import datetime
import timeseries_tools as tt

# Maximum ratio of beam radion on tilted panel vs horizontal surface
# (to deal with bad approximation when zenith angle is close to 90 degrees)
# (used in Rb calculation)
MAX_BEAM_PANEL_RATIO = 10
MIN_BEAM_PANEL_RATIO = 0

def decl(n):
    '''Declination of the sun at day number n

    Parameters
    ----------
    n : int
        day of year (0-365)
    Returns
    -------
    Declination of the sun (radians)
    '''
    return np.pi/180*23.45*np.sin(2*np.pi*(284+n)/365)


def hourangle(h,lon):
    '''Compute the hourangle of the sun at a given time

    Parameters
    ----------
    h : float
        hour of day (universal time, 0-24)
    lon : float
        longitude (degrees), east is positive, west is negative
    Returns
    -------
    hour angle in radians in the range (-pi,pi)
    '''
    w = np.pi/180*(15*(h-12)+lon)
    #shift to (-pi,pi) region
    idx= w>np.pi
    w[idx]=w[idx]-np.pi*2
    idx = w<-np.pi
    w[idx]=w[idx]+np.pi*2
    return w


def hourangle_sunset(n,lat):
    '''Compute the hourangle where sun sets (rad)

    Parameters
    ----------
    n : int
        day of year (0-356)
    lat : float
        latitude (degrees)
    Returns
    -------
    hour angle in radians
    '''
    cosws = -np.tan(lat*np.pi/180)*np.tan(decl(n))
    if cosws > 1:
        return -1 # no sunset, summer light
    elif cosws < -1:
        return np.pi
    else:
        return np.arccos(cosws)


def zenithangle(h,n,lat,lon):
    '''Compute the angle between zenith and the sun

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
    '''
    lat = lat*np.pi/180
    z = np.arccos(np.sin(lat)*np.sin(decl(n))
            + np.cos(lat)*np.cos(decl(n))*np.cos(hourangle(h,lon)))
    return z


def cpr(h,n,lat,lon):
    '''Collares-Pereira-Rabl (CPR) factor for decomposition of solar radiation

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
    CPR factor
    '''

    ws = hourangle_sunset(n,lat)
    a = 0.4090 + 0.5016*np.sin(ws-np.pi/3)
    b = 0.6609 - 0.4767*np.sin(ws-np.pi/3)
    w = hourangle(h=h,lon=lon)
    f= (np.cos(w)-np.cos(ws))/(np.sin(ws)-ws*np.cos(ws))
    f[w>ws] = 0
    f[w<-ws] = 0
    return [a,b,f]


def r_td(h,n,lat,lon):
    '''
    Compute decomposition factors r_t(h) and r_d(h)
    i.e. ratio of hourly to daily radiation:

    Parameters
    ----------
    n : int
        day of year (0-356)
    h : float
        hour of the day (0-24)
    lat, lon : float
        latitude, longitude (degrees)
    Returns
    -------
    [rt,rd] : float
        decomposition factors r_t(h) and r_d(h)
    '''
    [a,b,f] = cpr(h=h,n=n,lat=lat,lon=lon)
    w = hourangle(h=h,lon=lon)
    rt = np.pi/24*f*(a+b*np.cos(w)) # total=global
    rd = np.pi/24*f                 # diffuse
    rt[rt<0]=0
    rd[rd<0]=0
    return [rt, rd]


def cos_incidenceangle(slope,dph,h,n,lat,lon,azimuth):
    '''
    Cosine of the incidence angle between sun and PV panel

    Parameters
    ----------
    slope : float
        slope of PV panel (radians)
    dph : str
        tracker type ("south","fixed","azimuth","2-axis")
    n : int
        day of year (0-356)
    h : float
        hour of the day (universal time, 0-24)
    lat, lon : float
        latitude, longitude (degrees)
    azimuth : float
        azimuth position of PV panel (degrees)

    Returns
    -------
    cos(theta) : float
        cos of incidence angle
    '''
    thz= zenithangle(h,n,lat,lon)
    if dph=='south':
        # south facing, fixed azimuth angle
        dph_a = hourangle(h,lon)
    elif dph=='fixed':
        # fixed position, may be different from south
        dph_a = hourangle(h,lon-azimuth)
    elif dph=='tracker':
        # azimuth tracking
        dph_a = 0
    elif dph=='2-axis':
        # fixed azimuth angle different from south (not very useful)
        dph_a = dph
    else:
        raise Exception("undefined tracker type")
    costh = np.cos(thz)*np.cos(slope) + np.sin(thz)*np.sin(slope)*np.cos(dph_a)
    costh[costh>1] = 1 # may occur due to approximations?
    costh[costh<0] = 0 # costh<0 may occur close to sunrise/sunset
    return costh


def Rb(slope,dph,h,n,lat,lon,azimuth):
    '''
    Ratio of beam radiation onto the PV panel relative to horizontal surface

    Parameters
    ----------
    slope : float
        slope of PV panel (radians)
    dph : str
        tracker type ("south","fixed","azimuth","2-axis")
    n : int
        day of year (0-356)
    h : float
        hour of the day (universal time, 0-24)
    lat, lon : float
        latitude, longitude (degrees)
    azimuth : float
        azimuth position of PV panel (degrees)

    Returns
    -------
    Rb : float
        Rb factor
    '''
    thz = zenithangle(h,n,lat,lon)
    costh = cos_incidenceangle(slope,dph,h,n,lat,lon,azimuth)
    y = costh/np.cos(thz)

    # if zenith angle at midday is large (close to polar night)...
    # (or equivalently sunset angle is small)
    #if (zenithangle(np.array([12]),n,lat,0)>85*np.pi/180):
    if (hourangle_sunset(n,lat)< 30*np.pi/180):
        # poor approximation close to sunrise/sunset, si set to zero
        y[thz>85/180*np.pi] = 0


    #TODO: Determine how to deal with bad approximation close
    # to polar night in high latitudes

    # previous (un-tested) option
    # clip ratio to deal with bad approximation when zenith angle is
    # close to 90 degrees
    #y = np.clip(y,a_min=MIN_BEAM_PANEL_RATIO,a_max=MAX_BEAM_PANEL_RATIO)

    return y


def H_tiltedpanel(slope,dph,h,n,lat,lon,Hb,Hd,albedo,azimuth=0,use_rt=True):
    '''
    Hourly total irradiane on tilted surface

    Parameters
    ----------
    slope : float
        slope of PV panel (radians)
    dph : str
        tracker type ("south", "fixed", "azimuth", "2-axis")
    n : int
        day of year (0-356)
    h : float
        hour(s) of the day (universal time, 0-24)
    lat, lon : float
        latitude, longitude (degrees)
    Hb, Hd : float
        daily beam and diffuce irradiation on horizontal surface (Wh/m2)
    albedo : float
        albedo of the ground (reflection factor)
    azimuth : float
        azimuth position of PV panel (degrees)
    use_rt : boolean
        see code.

    Returns
    -------
    Radiation on tilted surface (PV panel)
    '''
    cosb = np.cos(slope)
    [rt, rd] = r_td(h=h,n=n,lat=lat,lon=lon)
    Hht = rt*(Hb+Hd)

    '''Note:
    If use_rt is true:
    Using the factor rt for both total, beam and diffuse radiation in order
    to avoid negative values for beam radiation (cf subtraction below).
    This is different from RETscreen, but seems reasonable instead of ad-hoc
    setting negative beam radiaion values to zero.
    '''
    if use_rt:
        Hhd = rt*Hd # changed from rd to rt (difference from RETscreen)
    else:
        Hhd = rd*Hd

    Hhb = Hht-Hhd
    # replace negative values with zero (may happen if use_rt=False)
    Hhb[Hhb<0]=0
    # scale beam and diffuse radiation so daily total matches
    if Hhb.sum()>0:
        scale_b = Hb/Hhb.sum()
        Hhb = Hhb * scale_b
    if Hhd.sum()>0:
        scale_d = Hd/Hhd.sum()
        Hhd = Hhd * scale_d

    # Liu and Jordan model
    # Liu, B.; Jordan, R. Daily insolation on surfaces tilted towards
    # equator. ASHRAE J. 1961, 10, 526–541
    y = ( Rb(slope,dph,h,n,lat,lon,azimuth)*Hhb + Hhd*(1+cosb)/2
            + (Hhb+Hhd)*albedo*(1-cosb)/2)

    # The following was used in the matlab code to get right values
    # for horizontal panels (slope=0). However, the calculation of
    # horizontal panel calculation was wrong, so the correction
    # factors were also wrong. Actually, correctly calculated factors
    # are nearly 1.00, so such a correction is not necessary
#    finalCorrect = False
#    if finalCorrect:
#        # same, zero slope:
#        y_hor = ( Rb(0,dph,h,n,lat,lon,azimuth)*Hhb + Hhd*(1+1)/2
#                + (Hhb+Hhd)*albedo*(1-1)/2)
#        if y_hor.sum()>0:
#            correctionFactor = (Hb+Hd)/y_hor.sum()
#            y = y * correctionFactor
#            print("Hb+Hd={}, y_hor={}, Correction factor = {}".format(
#                    Hb+Hd,y_hor.sum(),correctionFactor))

    return y


def H_tiltedpanel_components(Hb,Hd,albedo,lat,lon,n_day,hour,tracker='south',
                   slope=0,azimuth=0,use_rt=True):
    '''
    Hourly irradiance on tilted surface (decomposed)

    Parameters
    ----------
    Hb, Hd : float
        beam and diffuce irradiation on horizontal surface (W/m2)
    albedo : float
        albedo of the ground (reflection factor)
    n : int
        day of year (0-356)
    hour : float
        hour(s) of the day (universal time, 0-24)
    lat, lon : float
        latitude, longitude (degrees)
    tracker : str
        tracker type ("south", "fixed", "azimuth", "2-axis")
    slope : float
        slope of PV panel (radians)
    azimuth : float
        azimuth position of PV panel (degrees)

    Returns
    -------
    (rad_tot,rad_b,rad_d,rad_r) : floats
        Radiation on tilted surface (PV panel), total, beam, diffuce, reflected
        (W/m2)
    '''
    if tracker=='south':
        cosb = np.cos(slope) #fixed slope azimuth=0
        track_ew='south'
        if azimuth !=0:
            raise ValueError("With south tracking, azimuth must be zero")
    elif tracker=='fixed':
        cosb = np.cos(slope) #fixed slope and azimuth
        track_ew='fixed'
    elif tracker=='azimuth':
        cosb = np.cos(slope) #fixed slope
        track_ew='tracker'
    elif tracker=='2-axis':
        thz = zenithangle(h=hour,n=n_day,lat=lat,lon=lon) # = cos thz // altitude tracker
        cosb = np.cos(thz)
        slope = thz
        track_ew='tracker'
    else:
        raise ValueError("Tracker must be 'south','azimuth' or '2-axis'")
    [rt, rd] = r_td(hour,n=n_day,lat=lat,lon=lon)
    Hht = rt*(Hb+Hd)
    if use_rt:
        Hhd = rt*Hd # changed from rd to rt (difference from RETscreen)
    else:
        Hhd = rd*Hd

    Hhb = Hht-Hhd
    # replace negative values with zero (may happen if use_rt=False)
    Hhb[Hhb<0]=0
    # scale beam and diffuse radiation so daily total matches
    if Hhb.sum()>0:
        scale_b = Hb/Hhb.sum()
        Hhb = Hhb * scale_b
    if Hhd.sum()>0:
        scale_d = Hd/Hhd.sum()
        Hhd = Hhd * scale_d

    Rb1 = Rb(slope=slope,dph=track_ew,h=hour,n=n_day,lat=lat,lon=lon,azimuth=azimuth)
    rad_b = Hhb*Rb1
    rad_d = Hhd*(1+cosb)/2
    rad_r = (Hhb+Hhd)*albedo*(1-cosb)/2
    rad_tot = H_tiltedpanel(slope=slope,dph=track_ew,h=hour,n=n_day,
               lat=lat,lon=lon,Hb=Hb,Hd=Hd,albedo=albedo,azimuth=azimuth)
    return (rad_tot,rad_b,rad_d,rad_r)


def pvPower(H,eta_el):
    '''
    Hourly electrical powre from tilted PV panel

    Parameters
    ----------
    H : array of floats
        irradiation onto PV surface (W/m2)
    eta_el : float
        efficiency of electrical systems (ratio of exported power to
        power generated in PV panels)

    Returns
    -------
    Daily profile for electrical power from tilted PV panel (p.u.)
    '''
    p = eta_el*H/1000
    return p


def _checkSameVariables(datasets,variables):
    for d in datasets:
        for k,v in variables.items():
            if not list(d.variables[k])==list(v):
                print("ERROR: d={},k={},v={}".format(d,k,v))
                raise Exception("Inconsistent variables ({},{})"
                                .format(d.variables[k],k))

class SolarPower:
    def __init__(self,datarange,locations):
        self.datarange = datarange
        self.nc_lats = None
        self.nc_lons = None
        self.lons_r = None
        self.coords = None
        self.df_weight = None
        self.locations = locations
        self.locations["sol_slope_rad"] = self.locations["sol_slope"]*np.pi/180
        self.lats_selected = locations["latitude"]
        self.lons_selected = locations["longitude"]


    def getDataCoords(self):
        """Return latitudes and longitudes of Reanalysis data points"""
        return (self.nc_lats, self.nc_lons)

    def getAndValidateNcData(self,nc_path,year):
        datadict = {}

        nc_nbdsf = netCDF4.Dataset('{}/nbdsf.sfc.gauss.{}.nc'.format(nc_path,year))
        nc_nddsf = netCDF4.Dataset('{}/nddsf.sfc.gauss.{}.nc'.format(nc_path,year))
        nc_vbdsf = netCDF4.Dataset('{}/vbdsf.sfc.gauss.{}.nc'.format(nc_path,year))
        nc_vddsf = netCDF4.Dataset('{}/vddsf.sfc.gauss.{}.nc'.format(nc_path,year))
        nc_temp = netCDF4.Dataset('{}/air.2m.gauss.{}.nc'.format(nc_path,year))
        h_nb = nc_nbdsf.variables['nbdsf']
        h_nd = nc_nddsf.variables['nddsf']
        h_vb = nc_vbdsf.variables['vbdsf']
        h_vd = nc_vddsf.variables['vddsf']
        hT = nc_temp.variables['air']

        datadict['h_nb'] = h_nb # near, beam
        datadict['h_nd'] = h_nd # near, diffuse
        datadict['h_vb'] = h_vb # visible, beam
        datadict['h_vd'] = h_vd # visible, diffuse
        datadict['hT'] = hT # temperature


        if self.nc_lats  is None:
            self.nc_lats  = nc_nbdsf.variables['lat']
            self.nc_lons  = nc_nbdsf.variables['lon']
        lats = self.nc_lats
        lons = self.nc_lons
        times = nc_nbdsf.variables['time']
        timesT = nc_temp.variables['time']

        #Check that coordinates are the same:
        _checkSameVariables([nc_nbdsf,nc_nddsf,nc_vbdsf,nc_vddsf],
                           variables={'lat':lats,'lon':lons,'time':times})
        _checkSameVariables([nc_nbdsf,nc_temp],variables={'lat':lats,'lon':lons})

        units = h_nb.units
        jd = netCDF4.num2date(times[:],times.units,
            only_use_cftime_datetimes=False,only_use_python_datetimes=True)
        jdT = netCDF4.num2date(timesT[:],timesT.units,
            only_use_cftime_datetimes=False,only_use_python_datetimes=True)

        datadict['jd'] = jd
        datadict['jdT'] = jdT

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

        return datadict

    def MakePowerTimeSeriesSelection(self,nc_path,datarange,
            slope=0,azimuth=0,tracker="south",albedo=0.2, eta_eff=0.65,
            rad_scaling=1.0,Ninterpolate=3):
        '''
        nc_path : string
            Path to where Reanalyis data  (NC) files are kept
        datarange: dict
            Dictionary containing range of years, latitudes, longitudes,
            and (optionally) points of interest. If points are given,
        '''

        solar = {}
        temperatures = {}

        for year in range(datarange['year'][0],datarange['year'][1]+1):
            print('{}: Year={}'.format(datetime.datetime.now(),year))

            datadict = self.getAndValidateNcData(nc_path,year)
            h_nb = datadict['h_nb']
            h_nd = datadict['h_nd']
            h_vb = datadict['h_vb']
            h_vd = datadict['h_vd']
            hT = datadict['hT']
            jd = datadict['jd']
            jdT = datadict['jdT']

            # 1464x73x144 (time x lat x lon) -> 1464x10512 (time x latlon_grid)
            vdim = h_vb.shape[0]
            h_vb_2d = h_vb[:].reshape((vdim,-1),order="C")
            h_vd_2d = h_vd[:].reshape((vdim,-1),order="C")
            h_nb_2d = h_nb[:].reshape((vdim,-1),order="C")
            h_nd_2d = h_nd[:].reshape((vdim,-1),order="C")
            hT_2d = hT[:].reshape((hT.shape[0],-1),order="C")

            # If not done before, calculate geographical interpolation weights:
            if self.df_weight is None:
                print("Computing geographical interpolation weights")
                latlon_grid = tt.getDataGrid(self.nc_lats,self.nc_lons)
                latlon_select = pd.DataFrame({
                    'lats':self.lats_selected,
                    'lon':self.lons_selected})
                self.df_weight = tt.computeInterpolationWeights(
                        latlon_grid,latlon_select,Ninterpolate=Ninterpolate)

            for k,row in self.locations.iterrows():
            #for k,row in locations.loc[[5],:].iterrows():
                # weighted average of neighbour grid points. - sum per day: W -> Wh
                Hb = 24*(h_nb_2d+h_vb_2d).dot(self.df_weight[k])
                Hd = 24*(h_nd_2d+h_vd_2d).dot(self.df_weight[k])
                HT = hT_2d.dot(self.df_weight[k])

                lat = row['latitude']
                lon = row['longitude']
                slope = row['sol_slope_rad']
                eta_eff = row['sol_adj']
                #print("{}:({},{})".format(k,lat,lon),end=" ")
                thisone = pd.DataFrame()

                # NOTE! For compatibility with previous matlab code:
                # add 0.5 hours to get datapints in the middle of the interval
                # use use_rt=True in decomposition of radiation.
                use_rt=False
                h_add=0 #previously used 0.5 (middle of hourly time intervals)

               # loop over all days:
                for i in range(jd.shape[0]):
                    t = jd[i]
                    nday = t.timetuple().tm_yday
                    hours = np.linspace(0,23,24) + h_add
                    hours_indx = pd.date_range(start=t,periods=24,freq='H')
                    rad_tot = H_tiltedpanel(slope=slope,dph=tracker,h=hours,n=nday,
                                                lat=lat,lon=lon,Hb=Hb[i],Hd=Hd[i],
                                                albedo=albedo,azimuth=azimuth,
                                                use_rt=use_rt)
                    solP = pvPower(rad_tot,eta_eff)
                    # radiation (beam/diffuse) on horizontal surface hour-by-hour
        #            #rad_comp = (rad_tot,rad_b,rad_d,rad_r)
                    (rad0_tot,rad0_b,rad0_d,rad0_r) = H_tiltedpanel_components(
                            slope=0,tracker=tracker,
                            hour=hours,n_day=nday,
                            lat=lat,lon=lon,Hb=Hb[i],Hd=Hd[i],
                            albedo=albedo,azimuth=azimuth)
                    correctionFactor = 1
                    if (rad0_tot.sum()>0):
                        correctionFactor = (Hb[i]+Hd[i])/rad0_tot.sum()
                        solP = solP * correctionFactor


                    radiation =  pd.DataFrame({'rad_panel':rad_tot,'power':solP,
                                               'rad_beam':rad0_b,
                                               'rad_diffuse':rad0_d},
                                              index=hours_indx)
                    thisone = pd.concat([thisone,radiation])
        #        temperature = pd.DataFrame(hT[:,ind_lat,ind_lon],index=jdT)
                temperature = pd.DataFrame(HT,index=jdT)

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
            missingindx = pd.date_range(start=temperatures[k].index[-1]
                                        +datetime.timedelta(hours=1),
                                      periods=missinghours,freq='H')
            missingwind = temperatures[k].iloc[[-1]*missinghours].set_index(missingindx)
            temperatures[k] = temperatures[k].append(missingwind)
            solar[k]['T']=temperatures[k]

        print('{}: Done'.format(datetime.datetime.now()))
        return solar



def makePlotsForIllustration():
    n=np.linspace(0,365,200)
    d=180/np.pi*decl(n)
    h=np.linspace(0,24,200)

    z1=90-180/np.pi*zenithangle(h,n=81,lat=60,lon=15)
    z2=90-180/np.pi*zenithangle(h,n=172,lat=60,lon=15)
    z3=90-180/np.pi*zenithangle(h,n=354,lat=60,lon=15)
    z1b=90-180/np.pi*zenithangle(h,n=81,lat=30,lon=15)
    z2b=90-180/np.pi*zenithangle(h,n=172,lat=30,lon=15)
    z3b=90-180/np.pi*zenithangle(h,n=354,lat=30,lon=15)


    [rt,rd] = r_td(h+0.5,n=172,lat=30,lon=15)
    [rt2,rd2] = r_td(h+0.5,n=354,lat=30,lon=15)
    [rt3,rd3] = r_td(h+0.5,n=172,lat=60,lon=15)
    [rt4,rd4] = r_td(h+0.5,n=354,lat=60,lon=15)
    #rb = rt-rd #this is incorrect

    # Beam and diffuse radiation values:
    Hb = 1#0.01 #example for illustration
    Hd = 0.5#2.5 #example for illustration

    Hpv = H_tiltedpanel(slope=0*np.pi/180,dph='south',h=h,n=172,
               lat=60,lon=15,Hb=Hb,Hd=Hd,albedo=0.2)
    Hpv2 = H_tiltedpanel(slope=40*np.pi/180,dph='south',h=h,n=172,
                lat=60,lon=15,Hb=Hb,Hd=Hd,albedo=0.2)
    Hpv3 = H_tiltedpanel(slope=70*np.pi/180,dph='south',h=h,n=172,
                lat=60,lon=15,Hb=Hb,Hd=Hd,albedo=0.2)

    Hpvt = H_tiltedpanel(slope=0*np.pi/180,dph='tracker',h=h,n=172,
               lat=60,lon=15,Hb=Hb,Hd=Hd,albedo=0.2)
    Hpvt2 = H_tiltedpanel(slope=40*np.pi/180,dph='tracker',h=h,n=172,
                lat=60,lon=15,Hb=Hb,Hd=Hd,albedo=0.2)
    Hpvt3 = H_tiltedpanel(slope=70*np.pi/180,dph='tracker',h=h,n=172,
                lat=60,lon=15,Hb=Hb,Hd=Hd,albedo=0.2)

    plt.figure()
    plt.plot(n,d)
    plt.plot([0,365],[0,0],'--k')
    plt.xlim(0,365)
    plt.xlabel('n (day of year)')
    plt.ylabel('declination')
    plt.gcf().set_size_inches(5,2)
    plt.savefig('fig_declination.png', bbox_inches = 'tight')

    plt.figure()
    plt.plot([0,24],[0,0],'--k')
    plt.plot(h,z1,label='v')
    plt.plot(h,z2,label='s')
    plt.plot(h,z3,label='w')
    plt.xlim(0,24)
    plt.ylim(-90,90)
    plt.legend()
    plt.xlabel('hour')
    plt.title('Solar altitude angle, lat=60, lon=15')
    plt.gcf().set_size_inches(5,4)
    plt.savefig('fig_zenithangles60.png', bbox_inches = 'tight')

    plt.figure()
    plt.plot([0,24],[0,0],'--k')
    plt.plot(h,z1b,label='v')
    plt.plot(h,z2b,label='s')
    plt.plot(h,z3b,label='w')
    plt.xlim(0,24)
    plt.ylim(-90,90)
    plt.legend()
    plt.xlabel('hour')
    plt.title('Solar altitude angle, lat=30, lon=15')
    plt.gcf().set_size_inches(5,4)
    plt.savefig('fig_zenithangles30.png', bbox_inches = 'tight')


    plt.figure()
    plt.plot(h,rt,'g',label='summer')
    #plt.plot(h,rd,'--g',label='s diffuse')
    plt.plot(h,rt2,'r',label='winter')
    #plt.plot(h,rd2,'--r,',label='w diffuse')
    plt.legend()
    plt.title('Daily radiation profile, lat=30,lon=15')
    plt.xlabel('hour')
    plt.xlim(0,24)
    plt.ylim(0,0.30)
    plt.gcf().set_size_inches(5,4)
    plt.savefig('fig_dailyradiation30.png', bbox_inches = 'tight')

    plt.figure()
    plt.plot(h,rt3,'g',label='summer')
    #plt.plot(h,rd3,'--g',label='s diffuse')
    plt.plot(h,rt4,'r',label='winter')
    #plt.plot(h,rd4,'--r,',label='w diffuse')
    plt.legend()
    plt.title('Daily radiation profile, lat=60,lon=15')
    plt.xlabel('hour')
    plt.xlim(0,24)
    plt.ylim(0,0.30)
    plt.gcf().set_size_inches(5,4)
    plt.savefig('fig_dailyradiation60.png', bbox_inches = 'tight')


    plt.figure()
    plt.plot(h,Hpv,label='0')
    plt.plot(h,Hpv2,label='40')
    plt.plot(h,Hpv3,label='70')
    plt.xlim(0,24)
    #plt.ylim(0,0.3)
    plt.legend()
    plt.title('Radiation on PV panel, south facing')
    plt.xlabel('hour')
    plt.gcf().set_size_inches(5,4)
    plt.savefig('fig_pv_radiation1.png', bbox_inches = 'tight')

    plt.figure()
    plt.plot(h,Hpvt,label='0')
    plt.plot(h,Hpvt2,label='40')
    plt.plot(h,Hpvt3,label='70')
    plt.xlim(0,24)
    #plt.ylim(0,0.3)
    plt.legend()
    plt.xlabel('hour')
    plt.title('Radiation on PV panel, tracking')
    plt.gcf().set_size_inches(5,4)
    plt.savefig('fig_pv_radiation2.png', bbox_inches = 'tight')


def makePlotsForIllustration2():
    # Radiation components on PV panel
    Hb = 1#0.01 #example for illustration
    Hd = 0.5#2.5 #example for illustration
    albedo = 0.2
    lat=40
    lon=15
    #slope=35*np.pi/180
    slope=np.pi/2
    day=172
    #n=np.linspace(0,365,200)
    #d=180/np.pi*decl(n)
    h=np.linspace(0,24,200)

    for track in ['south','fixed','azimuth','2-axis']:
        azimuth=0
        if track=="fixed":
            azimuth=-30
        plt.figure()
        plt.title('Radiation on PV panel ({})'.format(track))
        (rad_tot,rad_b,rad_d,rad_r) = H_tiltedpanel_components(
                Hb=Hb,Hd=Hd,albedo=albedo,lat=lat,lon=lon,n_day=day,hour=h,
                tracker=track,slope=slope,azimuth=azimuth)
        plt.plot(h,rad_b,label='direct')
        plt.plot(h,rad_d,label='diffuse')
        plt.plot(h,rad_r,label='reflected')
        plt.plot(h,rad_tot,'--k',label='total',linewidth=2)
        #plt.plot(h,rt,label='rt')
        #plt.plot(h,rd,label='rd')
        plt.xlim([0,24])
        plt.ylim([0,0.20])
        plt.legend()
        plt.gcf().set_size_inches(5,4)
        plt.savefig('fig_pv_rad_{}.png'.format(track), bbox_inches = 'tight')




if __name__=="__main__":
    plt.close('all')
    makePlotsForIllustration()
    makePlotsForIllustration2()
