import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys
from netCDF4 import Dataset
from numpy import *
import numpy as np
from pylab import *
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from matplotlib.colors import ListedColormap
import nclcmaps
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from wrf import *
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap, maskoceans



### Read in wrf outfiles

#Thomspon
thomp = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/thompson_mp/wrfout_d03_2013-12-11_00:00:00')
tp_thomp = getvar(thomp, 'RAINNC', timeidx = ALL_TIMES)
elevation = getvar(thomp, 'HGT')
xland = getvar(thomp, 'XLAND')



#Morriosn2m
morr = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/morrison2m_mp/wrfout_d03_2013-12-11_00:00:00')
tp_morr = getvar(morr, 'RAINNC', timeidx = ALL_TIMES)

#Goddard
godd = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/goddard_mp/wrfout_d03_2013-12-11_00:00:00')
tp_godd = getvar(godd, 'RAINNC', timeidx = ALL_TIMES)


#wrf_file = "/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/thompson_mp/wrfout_d03_2013-12-11_00:00:00"
#fh = Dataset(wrf_file)
#lon = fh.variables['XLONG'][:]
#lat = fh.variables['XLAT'][:]
#precip = fh.variables['RAINNC'][:]

#Radar Precip
radar = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/precip_data/Tom_KTYX_QPE_2013121100_2013121200.nc')
lon_radar = radar.variables['Lon'][:]
lat_radar = radar.variables['Lat'][:]
tp_radar = radar.variables['QPE'][:,:]


#Smooth radar data

# Applies a much stronger smoother (radius = 10) to the area of cluter than the whole area of precip
#tp_radar_clutter = smoother_leah.smoother(tp_radar[:200,200:], 'rect', 10)
#tp_radar[165:184,425:455] = tp_radar_clutter[165:184,225:255]
#tp_radar = smoother_leah.smoother(tp_radar, 'crect', 4)
#np.save('/uufs/chpc.utah.edu/common/home/u1013082/lake_effect/numpy_arrays/tp_radar', tp_radar)

tp_radar = np.load('/uufs/chpc.utah.edu/common/home/u1013082/lake_effect/numpy_arrays/tp_radar.npy')



#%%


### Define colorbar using functions and data in nclcmaps module

colors1 = np.array(nclcmaps.colors['WhiteBlueGreenYellowRed'])#perc2_9lev'])
colors_int = colors1.astype(int)
colors = list(colors_int)
cmap_precip = nclcmaps.make_cmap(colors, bit=True)

colors1_t = np.array(nclcmaps.colors['OceanLakeLandSnow'])
colors_int_t = colors1_t.astype(int)
colors_t = list(colors_int_t)
cmap_terrain = nclcmaps.make_cmap(colors_t, bit=True)


#%%


#############   Plots ########################################
    


fig = plt.figure(num=None, figsize=(12,16), dpi=400, facecolor='w', edgecolor='k')
for j in range(1,4):
    subplot = 310 + j
    
    #Levels
    lev_el = np.arange(135,1050,25)
    lev_ellab = np.arange(0,1000,5)
    lev_tp = np.arange(5,60.01,2.5)
    lev_tplab = np.arange(5,60.01,5)
    lev_water = [1.5,2.5]
    
    #Map
    latlon = [-77.9, 43.1, -74.3, 44.3]
    map = Basemap(projection='merc',llcrnrlon=latlon[0],llcrnrlat=latlon[1],urcrnrlon=latlon[2],urcrnrlat=latlon[3],resolution='h')
    
    
    #Plot
    ax = plt.subplot(subplot,aspect = 'equal')
    plt.subplots_adjust(left=0.04, bottom=0.01, right=0.88, top=0.95, wspace=0.1, hspace=0)
    #plt.axis('off')
    lats, lons = latlon_coords(tp_thomp)
    x, y = map(to_np(lons), to_np(lats))
    #Label to loop over runs
    run = ['tp_thomp', 'tp_morr', 'tp_godd']
    model_run = eval(run[j-1])
    tp = map.contourf(x,y,model_run[144,:,:], lev_tp, cmap = cmap_precip, zorder = 3, alpha = 1, vmin = -2.5)  
    el = map.contourf(x,y,elevation[:,:], lev_el, cmap = cm.Greys, zorder = 2, alpha = 1)
    #el2 = map.contourf(x,y,elevation[:,:], lev_water, colors = ('lightsteelblue'), zorder = 2, alpha = 1)
    #el2 = map.contour(x,y,elevation[:,:], lev_water, colors = ('k'), zorder = 2, alpha = 1)
    #water = map.contourf(x,y,xland[:,:], lev_water, colors = ('lightsteelblue'), zorder = 2, alpha = 1) 
    map.drawcoastlines(linewidth = 1)
    map.drawstates()  
    map.fillcontinents(lake_color='lightsteelblue', color = 'white')



    #Labels
    sub_title = ['Thompson', 'Morrison2m', 'Goddard']
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(10000, 165000, sub_title[j-1], fontsize = 25, bbox = props, zorder = 5)

    #Title
    if j == 1:
        ax.set_title("Total Precipitation", fontsize = 30, y = 1.01) 

#Colorbar
cbaxes = fig.add_axes([0.915, 0.2, 0.035, 0.55])             
cbar = plt.colorbar(tp, cax = cbaxes, ticks = lev_tplab)
cbar.ax.tick_params(labelsize=16)
plt.text(-0.27, -0.05, 'mm', fontsize = 25)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/wrf/mp_sensitivity_nwp.png")
plt.close(fig)



#%%

#######################  Radar Plot ###########################################

fig = plt.figure(num=None, figsize=(12,6), dpi=400, facecolor='w', edgecolor='k')

subplot = 111

#Levels
lev_el = np.arange(135,1050,25)
lev_ellab = np.arange(0,1000,5)
lev_tp = np.arange(5,60.01,2.5)
lev_tplab = np.arange(5,60.01,5)
lev_water = [1.5,2.5]

#Map
latlon = [-77.9, 43.1, -74.3, 44.3]
map = Basemap(projection='merc',llcrnrlon=latlon[0],llcrnrlat=latlon[1],urcrnrlon=latlon[2],urcrnrlat=latlon[3],resolution='h')


#Plot
ax = plt.subplot(subplot,aspect = 'equal')
plt.subplots_adjust(left=0.04, bottom=0.01, right=0.8, top=0.95, wspace=0.1, hspace=0)
#plt.axis('off')

#Lat lon grid
lon = np.zeros((301,601))
lat = np.zeros((301,601))
for i in range(601):
    lat[:,i] = lat_radar
for i in range(301):
    lon[i,:] = lon_radar
x, y = map(lon, lat)

#Label to loop over runs
tp = map.contourf(x,y,tp_radar[:,:], lev_tp, cmap = cmap_precip, zorder = 3, alpha = 1, vmin = -2.5) 

#Lat lon for elevation grid
lats, lons = latlon_coords(tp_thomp)
x, y = map(to_np(lons), to_np(lats))
    
el = map.contourf(x,y,elevation[:,:], lev_el, cmap = cm.Greys, zorder = 2, alpha = 1)
#el2 = map.contourf(x,y,elevation[:,:], lev_water, colors = ('lightsteelblue'), zorder = 2, alpha = 1)
#el2 = map.contour(x,y,elevation[:,:], lev_water, colors = ('k'), zorder = 2, alpha = 1)
#water = map.contourf(x,y,xland[:,:], lev_water, colors = ('lightsteelblue'), zorder = 2, alpha = 1) 
map.drawcoastlines(linewidth = 1)
map.drawstates()  
map.fillcontinents(lake_color='lightsteelblue', color = 'white')

#Title
ax.set_title("Radar Derived Total Precipitation", fontsize = 25, y = 1.01) 

#Colorbar
cbaxes = fig.add_axes([0.815, 0.2, 0.025, 0.55])             
cbar = plt.colorbar(tp, cax = cbaxes, ticks = lev_tplab)
cbar.ax.tick_params(labelsize=13)
plt.text(-0.27, -0.07, 'mm', fontsize = 18)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/wrf/mp_sensitivity_nwp_radar_precip.png")
plt.close(fig)

