import numpy as np
import pandas as pd  
from netCDF4 import Dataset, num2date # netCDF4 handles netCDF files
import array
import gsw

import cartopy.crs as ccrs
import cartopy

import matplotlib.pyplot as plt # matplotlib contains functions for graphics and plot manipulation

import glob

#########

files1 = sorted(glob.glob('19*/profiles/*.nc'))
files2 = sorted(glob.glob('29*/profiles/*.nc'))
files3 = sorted(glob.glob('39*/profiles/*.nc'))
files4 = sorted(glob.glob('49*/profiles/*.nc'))
files5 = sorted(glob.glob('59*/profiles/*.nc'))
files6 = sorted(glob.glob('69*/profiles/*.nc'))
files7 = sorted(glob.glob('79*/profiles/*.nc'))

files = files1+files2+files3+files4+files5+files6+files7;


#nc = Dataset(files,'r')

#file = "1900204/profiles/D1900204_120.nc"

#nc = Dataset(file,'r')

###########

lon_grid = range(-35,-40,1);
lat_grid = range(-35,-40,1);

longg = [];
latt = [];

for i in range(len(files)): 
    #a = print(files[i])
    nc = Dataset(files[i])
    #print(nc.variables['LATITUDE'][:])
    
    lon = nc.variables['LONGITUDE'][:]
    lat = nc.variables['LATITUDE'][:]
    
    longg= np.append(longg,lon);
    latt= np.append(latt,lat);
    
    nc.close()
    
    for x in lon_grid:
        print(lon_grid)
        if lon[:] < (lon_grid[x]+3):
            if lon[:] > (lon_grid[x]-1):
                print('hi')
                

############

extent = [-60, -30, -55, -30]
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent)
ax.stock_img()
ax.coastlines(resolution='50m')

gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=2, color='gray', alpha=0.5, linestyle='--', draw_labels=True);
gl.xlabels_top = False
gl.ylabels_left = False
gl.ylabels_right=True
gl.xlines = True
#gl.xlocator = mticker.FixedLocator([120, 140, 160, 180, -160, -140, -120])
#gl.ylocator = mticker.FixedLocator([0, 20, 40, 60])
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER

ax.plot(longg,latt,'b.',markersize=1);

###########

siggy = gsw.density.sigma2(Temp[0:1590267],Sal)

Temp1 = Temp[0:1590267];
#Temp = Temp1[siggy>5.1]
Temp = Temp1[siggy<5.23]

#Sal = Sal[siggy>5.1]
Sal = Sal[siggy<5.23]

############

plt.show()
plt.hist(Sal,bins=np.arange(34,37,0.09),edgecolor='black');
plt.title("Salinity Histogram");
plt.ylabel("Number of points");
plt.xlabel("g/kg");

plt.show()
plt.hist(Temp,bins=np.arange(1,5,0.1),facecolor='red',edgecolor='black');
plt.title("Temperature Histogram");
plt.ylabel("Number of points");
plt.xlabel("Deg C");

plt.show()
plt.hist(siggy,bins=np.arange(3,23,0.5),facecolor='pink',edgecolor='black');
plt.title("Sigma2 Histogram");
plt.ylabel("Number of points");
plt.xlabel("g/kg");


##############
