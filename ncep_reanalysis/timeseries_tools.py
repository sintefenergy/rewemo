import numpy as np
import matplotlib.pyplot as plt
import geopandas
import shapely
import pandas as pd
import pyproj
import os
import sys
try:
    import geoplot
    import cartopy.crs as ccrs
except ImportError:
    print("Warning: geoplot and/or cartopy could not be imported.")
    pass

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


def computeInterpolationWeights(latlon_grid,latlon_select,Ninterpolate=3):
    '''Compute geographical interpolation weights

    parameters
    ----------
    latlon_grid : pandas.DataFrame
      latitude and longitude of data grid points
    latlon_select : pandas.DataFrame
      latitude and longitude of points of interest
    Ninterpolate : int
      number of points to interpolate (N=1: nearest neighbour)

    returns
    -------
    df_weight : pandas.DataFrame
      NxM table of interpolation weights. The columns are the M selected
      points of interest, and rows are tne N grid points.
    '''

    df_weight = pd.DataFrame(0,columns=latlon_select.index,index=latlon_grid.index)
    for i,row in latlon_select.iterrows():
        print(i,end=" ")
        latlon = row.values
        dist=pd.Series(_spherical_dist(latlon_grid.values,latlon))
        #df_weight.loc[:,i] = [0]*latlon_grid.shape[0]
        if (dist==0).any():
            # weight one where dist=0, zero elsewhere
            df_weight.loc[dist==0,i] = 1
        else:
            nsmallest = dist.nsmallest(n=Ninterpolate)
            # weight proportioinal to inverse distance
            nweight = 1/nsmallest
            # normalise
            nweight = nweight/nweight.sum()
            df_weight.loc[nsmallest.index,i] = nweight
    print(".")
    return df_weight

def getDataGrid(nc_lats,nc_lons):
    """Return dataframe with all lat,lon in reanalysis data"""

    if (nc_lats is None) or (nc_lons is None):
        raise Exception("No Reanalysis data has yet been retrieved.")

    g_lats = nc_lats[:]
    g_lons = [l if l<180 else l-360 for l in nc_lons[:]]
    repl = {-90:-89.999, 90:89.999}
    g_lats = [repl[x] if x in repl else x for x in g_lats]
    # this is a list of all (lat,lon) coordinates for which there is data:
    latlon_grid = np.array([[lat,lon] for lat in g_lats for lon in g_lons ])
    latlon_grid = pd.DataFrame(data=latlon_grid,columns=['lat','lon'])
#        latlon_grid = pd.DataFrame(latlon_grid,columns=['lat','lon'])
    return latlon_grid


def plotWeights(wind_obj,plotlabels=True,filename=None,
        proj=geoplot.crs.Mercator()):
    '''Draw grid points and data points on a map, with lines showing
    interpolation weights

    parameters
    ----------
    plotlabels : boolean
        wheter to show labels
    filename : str
        name of file for saving figure
    '''
    (nc_lats,nc_lons) = wind_obj.getDataCoords()
    latlon_grid = getDataGrid(nc_lats,nc_lons)
    latlon_select = pd.DataFrame({
        'lat':wind_obj.lats_selected,
        'lon':wind_obj.lons_selected})
    df_weight = wind_obj.df_weight

    gdf_grid = geopandas.GeoDataFrame(latlon_grid,geometry=[
            shapely.geometry.Point(x, y)
            for x, y in zip(latlon_grid.lon, latlon_grid.lat)]
        )
    gdf_select = geopandas.GeoDataFrame(latlon_select,geometry=[
            shapely.geometry.Point(x, y)
            for x, y in zip(latlon_select.lon, latlon_select.lat)]
        )
    extent = gdf_select.total_bounds

    world = geopandas.read_file(geoplot.datasets.get_path('world'))
    #proj = geoplot.crs.Mercator()
    #proj = geoplot.crs.Miller()
    #ax = geoplot.webmap(world, projection=proj)

    dfw = pd.melt(df_weight.rename_axis("point_grid").reset_index(),
                  id_vars="point_grid",var_name="point_select")
    dfw = dfw[dfw['value']!=0]
    dfw['from'] = dfw.merge(gdf_grid,left_on="point_grid",right_index=True)['geometry']
    dfw['to'] = dfw.merge(gdf_select,left_on="point_select",right_index=True)['geometry']
    lines = [shapely.geometry.LineString(xy) for xy in zip(dfw['from'], dfw['to'])]
    dfw = geopandas.GeoDataFrame(dfw,geometry=lines)
    dfw['w']=dfw['value']

    ax = geoplot.polyplot(world, projection=proj, figsize=(6,8), extent=extent,
                          facecolor="lightgray",edgecolor="white")
    geoplot.pointplot(gdf_grid,color="black",marker=".",ax=ax,extent=extent,
                      label="data grid",projection=proj)
    geoplot.pointplot(gdf_select,color="red",marker='o',ax=ax,extent=extent,
                      label="selected point",projection=proj)
    geoplot.sankey(dfw,ax=ax,scale="w",limits=[0,10],zorder=0,projection=proj)

    if plotlabels:
        for x, y, label in zip(gdf_select.geometry.x, gdf_select.geometry.y,
                gdf_select.index):
            plt.text(x, y, label, horizontalalignment='left',
                transform=ccrs.Geodetic())

    plt.legend(loc="upper left", facecolor='white', framealpha=1)

    if filename is not None:
        plt.savefig(filename,bbox_inches='tight',dpi=150)

def plotMap(wind_obj,windp,filename=None,proj=geoplot.crs.Mercator(),
    variable="capacityfactor"):
    '''Plot summary data on map

    parameters
    ----------
    windp : dict
        object holding power series (output from makePowerTimeSeriesSelection)
    plotlabels : boolean
        wheter to show labels
    filename : str
        name of file for saving figure
    '''
    latlon_select = pd.DataFrame({
        'lat':wind_obj.lats_selected,
        'lon':wind_obj.lons_selected})
    gdf_select = geopandas.GeoDataFrame(latlon_select,geometry=[
            shapely.geometry.Point(x, y)
            for x, y in zip(latlon_select.lon, latlon_select.lat)]
        )
    margin = 2
    extent = gdf_select.total_bounds+[-margin,-margin,margin,margin]

    world = geopandas.read_file(geoplot.datasets.get_path('world'))
    #proj = geoplot.crs.Mercator()
    #proj = geoplot.crs.Miller()
    #ax = geoplot.webmap(world, projection=proj)

    plotvalues=None
    if variable == "capacityfactor":
        plotvalues = pd.concat(windp).loc[:,"power"].unstack(0).mean()
    else:
        variable = None
    gdf_select[variable] = plotvalues

    ax = geoplot.polyplot(world, projection=proj, figsize=(6,8), extent=extent,
                          facecolor="lightgray",edgecolor="white")
    geoplot.pointplot(gdf_select,marker='o',ax=ax,extent=extent,
                      projection=proj, hue=variable,legend=True)

    #plt.legend(loc="upper left", facecolor='white', framealpha=1)

    if filename is not None:
        plt.savefig(filename,bbox_inches='tight',dpi=150)



def plotWeights_OLD(latlon_grid,latlon_select,df_weight,plotlabels=True,
                filename="fig_geo_interpol.png"):
    '''Draw grid points and data points on a map, with lines showing
    interpolation weights

    parameters
    ----------
    latlon_grid : pandas.DataFrame
      latitude and longitude of data grid points
    latlon_select : pandas.DataFrame
      latitude and longitude of points of interest
    df_weight : pandas.DataFrame
      NxM table of interpolation weights. The columns are the M selected
      points of interest, and rows are tne N grid points.
    '''
    gdf_grid = geopandas.GeoDataFrame(latlon_grid,geometry=[
            shapely.geometry.Point(x, y)
            for x, y in zip(latlon_grid.lon, latlon_grid.lat)]
        )
    gdf_select = geopandas.GeoDataFrame(latlon_select,geometry=[
            shapely.geometry.Point(x, y)
            for x, y in zip(latlon_select.lon, latlon_select.lat)]
        )
    # lat/lon:
    gdf_grid.crs = {'init' :'epsg:4326'}
    gdf_select.crs = {'init' :'epsg:4326'}

    xy_max = pyproj.transform(p1=pyproj.Proj(init='epsg:4326'),
                              p2=pyproj.Proj(init='epsg:3395'),
                              x=gdf_select['lon'].max()+2,
                              y=gdf_select['lat'].max()+2)
    xy_min = pyproj.transform(p1=pyproj.Proj(init='epsg:4326'),
                              p2=pyproj.Proj(init='epsg:3395'),
                              x=gdf_select['lon'].min()-2,
                              y=gdf_select['lat'].min()-2)

    fig, ax = plt.subplots(figsize=(6, 8))
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world = world[(world.name != "Antarctica")
                    & (world.name != "Fr. S. Antarctic Lands")]
    #world = geopandas.read_file("custom.geo.json")
    world = world.to_crs({'init': 'epsg:3395'})
    #world.plot(ax=ax)
    europe = world[world.continent == "Europe"]
    europe.plot(ax=ax)
    plt.xlim(xy_min[0],xy_max[0])
    plt.ylim(xy_min[1],xy_max[1])

#    num_lines = (df_weight!=0).sum().sum()
#    weightlines = pd.DataFrame(index=range(num_lines),columns=['geo','w'])
    weightlines = pd.DataFrame(columns=['geo','w'])
#    i=0
    for wcol in df_weight.columns:
        neighbours = df_weight.loc[df_weight[wcol]!=0,wcol]
        for ni,nv in neighbours.iteritems():
            wline = shapely.geometry.LineString(
                    [gdf_select.loc[wcol,'geometry'],
                     gdf_grid.loc[ni,'geometry']])
#            weightlines.loc[i,'geo'] = wline
#            weightlines.loc[i,'w'] = nv*10
            weightlines = weightlines.append({'geo':wline,'w':nv*10},
                                             ignore_index=True)
#            i=i+1
    weightlines = geopandas.GeoDataFrame(weightlines,geometry="geo")
    weightlines.crs = {'init' :'epsg:4326'}
    weightlines = weightlines.to_crs({'init': 'epsg:3395'})
    for i,wline in weightlines.iterrows():
        geopandas.plotting.plot_linestring_collection(
                ax,[wline['geo']],color='cyan',linewidth=wline['w'],zorder=1)

    # Change to mercator projection:
    gdf_grid=gdf_grid.to_crs({'init': 'epsg:3395'})
    gdf_select=gdf_select.to_crs({'init': 'epsg:3395'})

    gdf_grid.plot(ax=ax,color="k",marker=".",label="data grid",zorder=2)
    gdf_select.plot(ax=ax,color="r",marker="o",label="selected points",zorder=3)

    plt.legend(loc="upper left", facecolor='white', framealpha=1)
    plt.xticks(ticks=[])
    plt.yticks(ticks=[])

    if filename is not None:
        plt.savefig(filename,bbox_inches = 'tight',dpi=150)

def _maxmin(gdf):
    xlim = [min(g.x for g in gdf.geometry),max(g.x for g in gdf.geometry)]
    ylim = [min(g.y for g in gdf.geometry),max(g.y for g in gdf.geometry)]
    xmarg = (xlim[1]-xlim[0])/20
    ymarg = (ylim[1]-ylim[0])/20
    xlim = (xlim[0]-xmarg, xlim[1]+xmarg)
    ylim = (ylim[0]-ymarg, ylim[1]+ymarg)
    return {'x':xlim,'y':ylim}

def plotPoints(coords,outpath,wind_onland=None,labels=True,
               projection='epsg:3395',color=None):
    latlon_select = pd.DataFrame(coords,columns=['lat','lon'])
    if wind_onland is not None:
        latlon_select['onland'] = [wind_onland[k] for k in coords]

    gdf_select = geopandas.GeoDataFrame(latlon_select,geometry=[
            shapely.geometry.Point(x, y)
            for x, y in zip(latlon_select.lon, latlon_select.lat)]
        )
    gdf_select.crs = {'init' :'epsg:4326'} # lat/lon

    #projection = 'epsg:3395' # Mercator
    #projection = 'epsg:4326' # Lat/lon
    world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    world = world.to_crs({'init': projection})
    gdf_select = gdf_select.to_crs({'init': projection})

    fig, ax = plt.subplots(figsize=(6, 8))
    world.plot(ax=ax)
    if color is None:
        color="red"
    else:
        color= coords.groupby(color).ngroup()
    gdf_select.plot(ax=ax,marker="o",c=color)


    maxmin = _maxmin(gdf_select)
    plt.xlim(maxmin['x'])
    plt.ylim(maxmin['y'])
    if labels:
        for x,y,label in zip(gdf_select.geometry.x,gdf_select.geometry.y,
                             gdf_select.index):
            ax.annotate(label, xy=(x,y),xytext=(3,3),textcoords="offset points")




#    x,y = m(df_coords['lon'].values,df_coords['lat'].values)
#        if wind_onland is None:
#            colours = 'r'
#        else:
#            colours = ['r' if k==1 else 'b' for k in df_coords['onland']]
#    labs, levs = self.powercurve_refs.factorize()
#    m.scatter(x,y,marker='o',c=labs,cmap=plt.cm.Set1,vmin=0,vmax=8)
#    if labels:
#        for i in df_coords.index:
#            plt.text(x[i], y[i],i)

    plt.savefig(os.path.join(outpath,"fig_datapoints.png"),
                bbox_inches = 'tight')

def plotWithValues(gdf,column=None,labels=False,outpath=None):
    world = geopandas.read_file(geoplot.datasets.get_path('world'))
    proj = geoplot.crs.Mercator()
    #proj = geoplot.crs.Miller()
    #ax = geoplot.webmap(world, projection=proj)
    extent = gdf.total_bounds + pd.np.array([-1,-1,1,1])*2
    ax = geoplot.polyplot(world,projection=proj,extent=extent, figsize=(8,8),
                          facecolor="lightgray",edgecolor="white")
    geoplot.pointplot(gdf,hue=column,ax=ax,extent=extent)
    if labels:
        crs_latlon = ccrs.PlateCarree() #latlon proj
        for x,y,lab in zip(gdf.geometry.x,gdf.geometry.y,gdf.index):
            xy=ax.projection.transform_point(x,y,crs_latlon)
            ax.annotate(lab, xy=xy, xytext=(2,2),textcoords="offset points")
    if outpath is not None:
        plt.savefig(os.path.join(outpath,"fig_values_{}.png".format(column)),
                    bbox_inches = 'tight')

def saveToFile(tseries,outpath,remove_29feb=True,gzip=False,emps_format=False):
    """Save power time series to files

    Parameters
    ==========
    tseries : dict
        power time series data (output from makePowerTimeSeriesSelection)
    outpath : str
        Path where generated time series should be placed
    remove_29feb : bool
        Whether leap year days should be removed
    gzip : bool
        Whether to zip the output csv files
    emps_format : bool
        Create EMPS format files
    """
    print("Exporting to file")

    if emps_format:
        saveToFileEmps(tseries,outpath)
    else:
        if not os.path.exists(outpath):
            print("Making path for output data:\n{}".format(outpath))
            os.makedirs(os.path.abspath(outpath))

    for k,timeseries_data in tseries.items():
        print(k,end=" ")
        if remove_29feb:
            mask = ((timeseries_data.index.month==2) &
                    (timeseries_data.index.day==29))
            timeseries_data = timeseries_data[~mask]
        if gzip:
            timeseries_data.to_csv(
                    '{}/wind_{}_pc{}.csv.gz'.format(outpath,k,pc),
                    compression="gzip")
        else:
            timeseries_data.to_csv(
                    '{}/wind_{}_pc{}.csv'.format(outpath,k,pc))


def saveToFileEmps(tseries,outpath):
    """Save power time series to files in EMPS text file format

    Parameters
    ==========
    tseries : dict
        power time series data (output from makePowerTimeSeriesSelection)
    outpath : str
        Path where generated time series should be placed
    """
    print("Exporting to file")
    if not os.path.exists(outpath):
        print("Making path for output data:\n{}".format(outpath))
        os.makedirs(os.path.abspath(outpath))

    for k,timeseries_data in tseries.items():
        print(k,end=" ")

        # remove leap year days (29 Feb)
        mask = ((timeseries_data.index.month==2) &
                (timeseries_data.index.day==29))
        timeseries_data = timeseries_data[~mask]

        # remove last day (31 dec) to get 52*7=364 days
        mask = ((timeseries_data.index.month==12) &
                (timeseries_data.index.day==31))
        timeseries_data = timeseries_data[~mask]

        dfx=timeseries_data[["power"]].copy()
        y_num= len(dfx.index.year.unique())
        y_start=dfx.index.year[0]
        dfx["hour"] = dfx.index.hour
        dfx["date"] = dfx.index.date
        dfx=dfx.set_index(["date","hour"]).unstack("hour")
        headers= [
            "Antall aar;Startaar;Antall uker;Startuke;Sluttuke;Startdogn;"
            +"Type data(Vind=1, Tilsig=2);"
            +"Type oppløsning(Uke=1, Dogn=2, Time=3);",
            "; ".join([str(x) for x in
                        [y_num,y_start,52,1,52,0,1,3,"VIND"]]),
            "Vindserier på timenivå;"]
        fname = '{}/{}.csv'.format(outpath,k)
        with open(fname, 'w') as empsfile:
            for line in headers:
                empsfile.write(line+"\n")
            dfx.to_csv(empsfile,decimal=".",sep=";",
                float_format="{0:9.3f}".format,header=False,
                index=False,line_terminator="\n")
