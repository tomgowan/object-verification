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
import matplotlib.patches as mpatches
import nclcmaps
import smoother_leah
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
from wrf import *
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap, maskoceans
import scipy.signal



### Read in wrf outfiles

#Thompson
thomp = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/thompson_mp/wrfout_d03_2013-12-11_00:00:00_RAINNC_regrid.nc')
lon_thomp = thomp.variables['lon'][:]
lat_thomp = thomp.variables['lat'][:]
tp_thomp = thomp.variables['RAINNC'][144,:,:]

#Morrison
morr = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/morrison2m_mp/wrfout_d03_2013-12-11_00:00:00_RAINNC_regrid.nc')
lon_morr = morr.variables['lon'][:]
lat_morr = morr.variables['lat'][:]
tp_morr = morr.variables['RAINNC'][144,:,:]

#Goddard
godd = Dataset('/uufs/chpc.utah.edu/common/home/steenburgh-group7/tom/wrf/wrf3.9.1b/WRFV3/run/goddard_mp/wrfout_d03_2013-12-11_00:00:00_RAINNC_regrid.nc')
lon_godd = godd.variables['lon'][:]
lat_godd = godd.variables['lat'][:]
tp_godd = godd.variables['RAINNC'][144,:,:]

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

#Combine model data
tp_model = np.zeros((3,301,601))
tp_model[0,:,:] = tp_thomp
tp_model[1,:,:] = tp_morr
tp_model[2,:,:] = tp_godd
loc_save = zeros((3,2))

#%%
####################################################################################################################
####### Object-based verifficaiton (http://www.cawcr.gov.au/projects/verification/CRA/CRA_verification.html) #######
####################################################################################################################
obj_model = np.zeros((3,200,340))
cc = np.zeros((3,10000,3))
best_model_anal = np.zeros((3,200,340))
best_model_plot = np.zeros((3,200,340))
orig_model_anal = np.zeros((3,200,340))
orig_model_plot = np.zeros((3,200,340))
loc_save = zeros((3,2))

#Loop over model run
for mod in range(3):
    #Select region of interest (lake-effect band)
    obj_model[mod,:,:] = np.array(tp_model[mod,:200,200:540])
    obj_radar = np.array(tp_radar[:200,200:540])
    
    #Choose threshold (>15mm)
    cc_model = np.array(obj_model[mod,:,:])
    cc_radar = np.array(obj_radar)
    
    cc_model[cc_model < 15] = 0
    cc_radar[cc_radar < 15] = 0
    

    w = 0
    #To find best match, I shift object 100 gridpoints in every direction and calc CC
    for i in range(100):
        for j in range(100):
            cc_match_model = zeros((200,340))
            cc_match_model[i:,j:] = cc_model[:200-i,:340-j]
            #Unravel into 1D array to compute CC
            x = cc_match_model.ravel()
            y = cc_radar.ravel()
            #Compute Correlation Coeffcient
            cc[mod,w,0] = np.corrcoef(x,y)[0,1]
            cc[mod,w,1] = i
            cc[mod,w,2] = j
            w = w + 1
            
    #New Array for best fit
    m = np.nanmax(cc[mod,:,0])
    loc = np.where(cc[mod,:,0] == m)
    locint = int(loc[0])
    ii = int(cc[mod,locint,1])
    jj = int(cc[mod,locint,2])
    loc_save[mod,:] = [ii,jj]
    
    ### Best Match ###
    best_model_anal[mod,ii:,jj:] =  cc_model[:200-ii,:340-jj]  #For analysis
    best_model_plot[mod,ii:,jj:] =  obj_model[mod,:200-ii,:340-jj] #For plotting
    
    ### Original ###
    orig_model_anal[mod,:,:] =  cc_model[:,:]  #For analysis
    orig_model_plot[mod,:,:] =  obj_model[mod,:,:] #For plotting
    
    
    
    
    
    
#%%
###############################################################################
#######################  Calculations  ########################################
###############################################################################

## Create Pandas data frame to store data
rows = ["Number of Gridpoints >15mm", "Average (mm)", "Max (mm)", "Volume (km^3)", "Displacement NS (km)", "Displacement EW (km)"]
cols = ["Thompson", "Morrison2m", "Goddard", "Radar Derived"]

r = pd.Index(["Number of Gridpoints >15mm", "Average (mm)", "Max (mm)", "Volume (km^3)", "Displacement NS (km)", "Displacement EW (km)"], name="rows")
c = pd.Index(["Thompson", "Morrison2m", "Goddard", "Radar Derived"], name="columns")
sum_calcs = pd.DataFrame(data=np.zeros((6,4)), index=r, columns=c)
sum_calcs.set_value(rows[0],cols[0],3)

##################### Number of Gridpoints >15mm ##############################
calc = 0
for mod in range(3):
    x = np.where(orig_model_anal[mod,:,:] > 0)
    num = len(x[0])
    sum_calcs.set_value(rows[calc],cols[mod],num)
# Radar  
x = np.where(cc_radar[:,:] > 0)
num = len(x[0])
sum_calcs.set_value(rows[calc],cols[3],num)


########################### Average (mm) ######################################
calc = 1
for mod in range(3):
    x = np.average(orig_model_anal[mod,:,:], weights=(orig_model_anal[mod,:,:]> 0))
    sum_calcs.set_value(rows[calc],cols[mod],x)
# Radar  
x = np.average(cc_radar, weights =(cc_radar>0))
sum_calcs.set_value(rows[calc],cols[3],x)


########################### Max (mm) ##########################################
calc = 2
for mod in range(3):
    x = np.max(orig_model_anal[mod,:,:],)
    sum_calcs.set_value(rows[calc],cols[mod],x)
# Radar  
x = np.max(cc_radar)
sum_calcs.set_value(rows[calc],cols[3],x)


########################### Volume (km^3) ##########################################
calc = 3
for mod in range(3):
    avg = np.average(orig_model_anal[mod,:,:], weights=(orig_model_anal[mod,:,:]> 0))/10**6
    grpts = np.where(orig_model_anal[mod,:,:] > 0)
    area = len(grpts[0])*1.2321 #Approx area of gridpoint
    vol = avg*area
    sum_calcs.set_value(rows[calc],cols[mod],vol)
# Radar  
avg = np.average(cc_radar, weights =(cc_radar>0))/10**6
grpts = np.where(cc_radar[:,:] > 0)
area = len(grpts[0])*1.2321 #Approx area of gridpoint
vol = avg*area
sum_calcs.set_value(rows[calc],cols[3],vol)



########################### Displacement NS (km) ##############################
calc = 4
for mod in range(3):
    dist = loc_save[mod,0]*1.11#km in .001 degree
    sum_calcs.set_value(rows[calc],cols[mod],dist)
# Radar  
sum_calcs.set_value(rows[calc],cols[3],np.NaN)

    
########################### Displacement EW (km) ##############################
calc = 5
for mod in range(3):
    dist = loc_save[mod,1]*1.11#km in .001 degree
    sum_calcs.set_value(rows[calc],cols[mod],dist)
# Radar  
sum_calcs.set_value(rows[calc],cols[3],np.NaN)





###############################################################################
#######################  Error Calcualtions  ##################################
###############################################################################

    
rows = ["RMS Error (mm)", "Correlation Coefficient","MSE Displacement (mm^2)", "MSE Volume (mm^2)", "MSE Pattern (mm^2)", "Displacement Error (%)", "Volume Error (%)", "Pattern Error (%)"]
cols = ["Thompson Original", "Morrison2m Original", "Goddard Original","Thompson 'Best Match'", "Morrison2m 'Best Match'", "Goddard 'Best Match'"]
r = pd.Index(["RMS Error (mm)", "Correlation Coefficient", "MSE Displacement (mm^2)", "MSE Volume (mm^2)", "MSE Pattern (mm^2)", "Displacement Error (%)", "Volume Error (%)", "Pattern Error (%)"], name="rows")
c = pd.Index(["Thompson Original", "Morrison2m Original", "Goddard Original","Thompson 'Best Match'", "Morrison2m 'Best Match'", "Goddard 'Best Match'"], name="columns")
error_calcs = pd.DataFrame(data=np.zeros((8,6)), index=r, columns=c)
    
    
########################### RMS Error (mm) ##############################
calc = 0
for mod in range(3):
    diff_orig = orig_model_anal[mod,:,:]-cc_radar
    diff_best = best_model_anal[mod,:,:]-cc_radar
    rms_orig = np.sqrt(np.mean(diff_orig**2))
    rms_best = np.sqrt(np.mean(diff_best**2))
    
    error_calcs.set_value(rows[calc],cols[mod],rms_orig)
    error_calcs.set_value(rows[calc],cols[mod+3],rms_best)
 

########################### Correlation Coefficient ###########################
calc = 1
for mod in range(3):
    cc_orig = cc[mod,0,0]
    cc_best = np.nanmax(cc[mod,:,0])
    
    error_calcs.set_value(rows[calc],cols[mod],cc_orig)
    error_calcs.set_value(rows[calc],cols[mod+3],cc_best)
    
    
########################### MSE Displacement (mm^2) ############################
calc = 2
for mod in range(3):
    diff_orig = orig_model_anal[mod,:,:]-cc_radar
    diff_best = best_model_anal[mod,:,:]-cc_radar
    rms_orig = np.mean(diff_orig**2)
    rms_best = np.mean(diff_best**2)
    
    mse = rms_orig-rms_best
    
    error_calcs.set_value(rows[calc],cols[mod],mse)
    error_calcs.set_value(rows[calc],cols[mod+3],np.nan)
    

########################### MSE Volume (mm^2) ##################################
calc = 3
for mod in range(3):
    f = np.mean(best_model_anal[mod,:,:])
    x = np.mean(cc_radar)
    mse = (f-x)**2
    
    
    error_calcs.set_value(rows[calc],cols[mod],mse)
    error_calcs.set_value(rows[calc],cols[mod+3],np.nan)
    

########################### MSE Pattern (mm^2) #################################
calc = 4
for mod in range(3):
    diff_best = best_model_anal[mod,:,:]-cc_radar
    rms_best = np.mean(diff_best**2)
    mse = rms_best - error_calcs.get_value(rows[3],cols[mod])
    
    error_calcs.set_value(rows[calc],cols[mod],mse)
    error_calcs.set_value(rows[calc],cols[mod+3],np.nan)
    

########################### Displacement Error (%) ############################
calc = 5
for mod in range(3):
    diff_orig = orig_model_anal[mod,:,:]-cc_radar
    diff_best = best_model_anal[mod,:,:]-cc_radar
    rms_orig = np.mean(diff_orig**2)
    rms_best = np.mean(diff_best**2)
    
    mse = ((rms_orig-rms_best)/((error_calcs.get_value(rows[0],cols[mod]))**2)*100)
    
    error_calcs.set_value(rows[calc],cols[mod],mse)
    error_calcs.set_value(rows[calc],cols[mod+3],np.nan)
    

########################### Volume Error (%) ##################################
calc = 6
for mod in range(3):
    f = np.mean(best_model_anal[mod,:,:])
    x = np.mean(cc_radar)
    mse = (((f-x)**2)/((error_calcs.get_value(rows[0],cols[mod]))**2)*100)
    
    
    error_calcs.set_value(rows[calc],cols[mod],mse)
    error_calcs.set_value(rows[calc],cols[mod+3],np.nan)
    

########################### Pattern Error (%) #################################
calc = 7
for mod in range(3):
    mse = 100 - (error_calcs.get_value(rows[2],cols[mod]) + error_calcs.get_value(rows[3],cols[mod]))
    
    error_calcs.set_value(rows[calc],cols[mod],mse)
    error_calcs.set_value(rows[calc],cols[mod+3],np.nan)
    
## Save to csv file
sum_calcs.to_csv('nwp_obj_sum_calcs.csv')  
error_calcs.to_csv('nwp_obj_error_calcs.csv')  
   



     
#%%


############   Plots ########################################
#Correlation Coeffcients

fig = plt.figure(num=None, figsize=(16,4.6), dpi=400, facecolor='w', edgecolor='k')

size = 1.5
plot = 131
ax1 = fig.add_subplot(plot)
y = best_model_anal[0,:,:].ravel()
y2 = orig_model_anal[0,:,:].ravel()
plt.grid(True)
x = cc_radar.ravel()
plt.scatter(x,y2, s = size, c = 'green')
plt.scatter(x,y, s = size, c = 'blue')
plt.xlim(15, 65)
plt.ylim(15, 65)
#Labels
sub_title = 'Thompson'
props = dict(boxstyle='square', facecolor='white', alpha=1)
ax1.text(17.5, 60, sub_title, fontsize = 17, bbox = props, zorder = 5)
ax1.text(37, 19, 'Original CC = %0.3f' % cc[0,0,0], fontsize = 12, color = 'red', zorder = 5)
ax1.text(37, 16, 'Best-Match CC = %0.3f' % np.nanmax(cc[0,:,0]), fontsize = 12, color = 'blue', zorder = 5)
ax1.set_ylabel('WRF Forecast Precip. (mm)', fontsize  = 11, labelpad = 10)
ax1.set_xlabel('Radar Derived Precip. (mm)', fontsize  = 11, labelpad = 10)

plot = 132
ax1 = fig.add_subplot(plot)
y = best_model_anal[1,:,:].ravel()
y2 = orig_model_anal[1,:,:].ravel()
plt.grid(True)
x = cc_radar.ravel()
plt.scatter(x,y2, s = size, c = 'green')
plt.scatter(x,y, s = size, c = 'blue')
plt.xlim(15, 65)
plt.ylim(15, 65)
sub_title = 'Morrison2m'
props = dict(boxstyle='square', facecolor='white', alpha=1)
ax1.text(17.5, 60, sub_title, fontsize = 17, bbox = props, zorder = 5)
ax1.text(37, 19, 'Original CC = %0.3f' % cc[1,0,0], fontsize = 12, color = 'red', zorder = 5)
ax1.text(37, 16, 'Best-Match CC = %0.3f' % np.nanmax(cc[1,:,0]), fontsize = 12, color = 'blue', zorder = 5)
ax1.set_xlabel('Radar Derived Precip. (mm)', fontsize  = 11, labelpad = 10)

plot = 133
ax1 = fig.add_subplot(plot)
y = best_model_anal[2,:,:].ravel()
y2 = orig_model_anal[2,:,:].ravel()
plt.grid(True)
x = cc_radar.ravel()
plt.scatter(x,y2, s = size, c = 'green')
plt.scatter(x,y, s = size, c = 'blue')
plt.xlim(15, 65)
plt.ylim(15, 65)
sub_title = 'Goddard'
props = dict(boxstyle='square', facecolor='white', alpha=1)
ax1.text(17.5, 60, sub_title, fontsize = 17, bbox = props, zorder = 5)
ax1.text(37, 19, 'Original CC = %0.3f' % cc[2,0,0], fontsize = 12, color = 'red', zorder = 5)
ax1.text(37, 16, 'Best-Match CC = %0.3f' % np.nanmax(cc[2,:,0]), fontsize = 12, color = 'blue', zorder = 5)
ax1.set_xlabel('Radar Derived Precip. (mm)', fontsize  = 11, labelpad = 10)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/wrf/mp_sensitivity_nwp_correlation.png")
plt.close(fig)


#%%
############   Plots ########################################
#3 Panel with all three outcomes overlayed    


fig = plt.figure(num=None, figsize=(12,15), dpi=400, facecolor='w', edgecolor='k')
for mod in range(1,4):
    subplot = 310 + mod
    
    
    #Levels
    lev_el = np.arange(135,1050,25)
    lev_ellab = np.arange(0,1000,5)
    lev_tp = [15,100]
    lev_tp2 = np.arange(15,25,10)
    lev_tplab = np.arange(5,60.01,5)
    lev_water = [1.5,2.5] 
    
    #Map
    latlon = [-77.9, 43.2, -74.3, 44.2]
    map = Basemap(projection='merc',llcrnrlon=latlon[0],llcrnrlat=latlon[1],urcrnrlon=latlon[2],urcrnrlat=latlon[3],resolution='h')
    
    #Plot
    ax = plt.subplot(subplot,aspect = 'equal')
    plt.subplots_adjust(left=0.07, bottom=0.1, right=0.92, top=0.95, wspace=0.1, hspace=0)
    
    #Lat lon grid
    lon = np.zeros((301,601))
    lat = np.zeros((301,601))
    for i in range(601):
        lat[:,i] = lat_radar
    for i in range(301):
        lon[i,:] = lon_radar
    x, y = map(lon[:200,200:540], lat[:200,200:540])

    
    #Contours
    tp = map.contourf(x,y,orig_model_plot[mod-1,:,:], lev_tp, colors = ('blue'), zorder = 3, alpha = 0.4) 
    tp2 = map.contour(x,y,orig_model_plot[mod-1,:,:], lev_tp2, linewidths = 1.8, colors = ('blue'), zorder = 5, alpha = 1) 
    tp = map.contourf(x,y,best_model_plot[mod-1,:,:], lev_tp, colors = ('green'), zorder = 3, alpha = 0.4) 
    tp2 = map.contour(x,y,best_model_plot[mod-1,:,:], lev_tp2, linewidths = 1.8, colors = ('green'), zorder = 5, alpha = 1) 
    tp1 = map.contourf(x,y,obj_radar, lev_tp, colors = ('red'), zorder = 3, alpha = 0.4) 
    tp3 = map.contour(x,y,obj_radar, lev_tp2,  linewidths = 1.8, colors = ('red'), zorder = 5, alpha = 1)  
    #el = map.contourf(x,y,elevation[:,:], lev_el, cmap = cm.Greys, zorder = 2, alpha = 1)
    map.drawcoastlines(linewidth = 1)
    map.drawstates()  
    map.fillcontinents(lake_color='lightsteelblue', color = 'white')


    #Labels
    sub_title = ['Thompson', 'Morrison2m', 'Goddard']
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(10000, 133000, sub_title[mod-1], fontsize = 25, bbox = props, zorder = 5)

    #Title
    if mod == 1:
        ax.set_title("Precipitation Objects (>15mm)", fontsize = 30, y = 1.02) 

#Legend
blue_patch = mpatches.Patch(color='blue', label='WRF Original',alpha = 0.4, edgecolor="red")
green_patch = mpatches.Patch(color='green', label='WRF "Best-Match"', alpha = 0.4, edgecolor="red")
red_patch = mpatches.Patch(color='red', label='Radar Derived', alpha = 0.4, edgecolor="red")
plt.legend(bbox_to_anchor=(1.03, -0.09), handles=[blue_patch, green_patch, red_patch], fontsize = 21, ncol = 3)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/wrf/mp_sensitivity_nwp_objects_3panel.png")
plt.close(fig)








#%%

#6 Panel with all thre outcomes overlayed    

mod = 0
mod_best = 0

fig = plt.figure(num=None, figsize=(16,10), dpi=400, facecolor='w', edgecolor='k')


for p in range(1,7):
    subplot = 320 + p
    
    
    #Levels
    lev_el = np.arange(135,1050,25)
    lev_ellab = np.arange(0,1000,5)
    lev_tp = [15,100]
    lev_tp2 = np.arange(15,100,15)
    lev_tplab = np.arange(5,60.01,5)
    lev_water = [1.5,2.5] 
    
    #Map
    latlon = [-77.9, 43.2, -74.3, 44.2]
    map = Basemap(projection='merc',llcrnrlon=latlon[0],llcrnrlat=latlon[1],urcrnrlon=latlon[2],urcrnrlat=latlon[3],resolution='h')
    
    #Plot
    ax = plt.subplot(subplot,aspect = 'equal')
    plt.subplots_adjust(left=0.07, bottom=0.1, right=0.93, top=0.93, wspace=0.1, hspace=0)
    
    #Lat lon grid
    lon = np.zeros((301,601))
    lat = np.zeros((301,601))
    for i in range(601):
        lat[:,i] = lat_radar
    for i in range(301):
        lon[i,:] = lon_radar
    x, y = map(lon[:200,200:540], lat[:200,200:540])

    
    #Contours
    if p == 1 or p == 3 or p ==  5:
        print(mod)
        print(p)
        tp = map.contourf(x,y,obj_model[mod,:,:], lev_tp, colors = ('blue'), zorder = 3, alpha = 0.4) 
        tp2 = map.contour(x,y,obj_model[mod,:,:], lev_tp2, linewidths = 1.8, colors = ('blue'), zorder = 5, alpha = 1)
        mod = mod + 1

    else:
        tp = map.contourf(x,y,best_model_plot[mod_best,:,:], lev_tp, colors = ('green'), zorder = 3, alpha = 0.4) 
        tp2 = map.contour(x,y,best_model_plot[mod_best,:,:], lev_tp2, linewidths = 1.8, colors = ('green'), zorder = 5, alpha = 1) 
        mod_best = mod_best + 1
    tp1 = map.contourf(x,y,obj_radar, lev_tp, colors = ('red'), zorder = 3, alpha = 0.4) 
    tp3 = map.contour(x,y,obj_radar, lev_tp2,  linewidths = 1.8, colors = ('red'), zorder = 5, alpha = 1)  
    #el = map.contourf(x,y,elevation[:,:], lev_el, cmap = cm.Greys, zorder = 2, alpha = 1)
    map.drawcoastlines(linewidth = 1)
    map.drawstates()  
    map.fillcontinents(lake_color='lightsteelblue', color = 'white')


    #Labels
    sub_title = ['Thompson', 'Thompson','Morrison2m', 'Morrison2m','Goddard','Goddard']
    props = dict(boxstyle='square', facecolor='white', alpha=1)
    ax.text(10000, 132000, sub_title[p-1], fontsize = 20, bbox = props, zorder = 5)

#Title
plt.suptitle("Precipitation Objects (>15mm)", fontsize = 30, y = 0.97) 

#Legend
blue_patch = mpatches.Patch(color='blue', label='WRF Original',alpha = 0.4, edgecolor="red")
green_patch = mpatches.Patch(color='green', label='WRF "Best-Match"', alpha = 0.4, edgecolor="red")
red_patch = mpatches.Patch(color='red', label='Radar Derived', alpha = 0.4, edgecolor="red")
plt.legend(bbox_to_anchor=(0.8, -0.09), handles=[blue_patch, green_patch, red_patch], fontsize = 21, ncol = 3)


plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/wrf/mp_sensitivity_nwp_objects_6panel.png")
plt.close(fig)




#%%

### Quick plot
#mod = 2
#fig = plt.figure(num=None, figsize=(8,5), dpi=200, facecolor='w', edgecolor='k')
#subplot = 111
##Levels
#lev_el = np.arange(135,1050,25)
#lev_ellab = np.arange(0,1000,5)
#lev_tp = [15,100]
#lev_tp2 = np.arange(15,100,10)
#lev_tplab = np.arange(5,60.01,5)
#lev_water = [1.5,2.5]    
##Map
#latlon = [-77.9, 43.1, -74.3, 44.3]
#map = Basemap(projection='merc',llcrnrlon=latlon[0],llcrnrlat=latlon[1],urcrnrlon=latlon[2],urcrnrlat=latlon[3],resolution='h')
##Plot
#ax = plt.subplot(111)
##plt.axis('off')
#lon = np.zeros((301,601))
#lat = np.zeros((301,601))
#for i in range(601):
#    lat[:,i] = lat_radar
#for i in range(301):
#    lon[i,:] = lon_radar
#    
#x, y = map(lon[:200,200:540], lat[:200,200:540])
##Label to loop over runs
#tp = map.contourf(x,y,obj_model[mod,:,:], lev_tp, colors = ('blue'), zorder = 3, alpha = 0.5) 
#tp2 = map.contour(x,y,obj_model[mod,:,:], lev_tp2, linewidths = 0.4, colors = ('blue'), zorder = 5, alpha = 1) 
#tp = map.contourf(x,y,best_model_plot[mod,:,:], lev_tp, colors = ('green'), zorder = 3, alpha = 0.5) 
#tp2 = map.contour(x,y,best_model_plot[mod,:,:], lev_tp2, linewidths = 0.4, colors = ('green'), zorder = 5, alpha = 1) 
#tp1 = map.contourf(x,y,obj_radar, lev_tp, colors = ('red'), zorder = 3, alpha = 0.5) 
#tp3 = map.contour(x,y,obj_radar, lev_tp2,  linewidths = .4, colors = ('red'), zorder = 5, alpha = 1) 
#plt.savefig("/uufs/chpc.utah.edu/common/home/u1013082/public_html/phd_plots/wrf/mp_sensitivity_nwp_obj_test.png")
#plt.close(fig)

#
#
##%%
#
#
#### Define colorbar using functions and data in nclcmaps module
#
#colors1 = np.array(nclcmaps.colors['WhiteBlueGreenYellowRed'])#perc2_9lev'])
#colors_int = colors1.astype(int)
#colors = list(colors_int)
#cmap_precip = nclcmaps.make_cmap(colors, bit=True)
#
#colors1_t = np.array(nclcmaps.colors['OceanLakeLandSnow'])
#colors_int_t = colors1_t.astype(int)
#colors_t = list(colors_int_t)
#cmap_terrain = nclcmaps.make_cmap(colors_t, bit=True)
#
#
##%%




