############## script for MIMOC data ##############
#
# Written by: Rosalind Echols, UW Oceanography (rechols@uw.edu) 
#
####################################################

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import gsw
import cmocean

def import_data(file_name):
    '''This function works regardless of the file called as long as a list of 
    variables can be identified easily in the file.'''
    
    fh = Dataset(file_name, mode = 'r')
    data = {}
    
    for vname in list(fh.variables):
        data[str(vname)] = fh.variables[vname][:]
        
    fh.close()
            
    return data


def lat_lon_indices(latitudes,longitudes,latlonrange):
    '''This function identifies the index values for latitude and longitude
    that correspond to the desired bounds for studying a particular region so 
    that the data can be sliced before processing. It returns data in the form
    [lon_min,lon_max,lat_min,lat_max] based on a latlonrange input of the
    same value.'''
    
    index=np.zeros(4,dtype='int') 
    
    index[0] = int(next(i for i,j in enumerate(longitudes) if j == latlonrange[0]))
    index[1] = int(next(i for i,j in enumerate(longitudes) if j == latlonrange[1]))
    index[2] = int(next(i for i,j in enumerate(latitudes) if j == latlonrange[2]))
    index[3] = int(next(i for i,j in enumerate(latitudes) if j == latlonrange[3]))
                            
    return index


#MIMOC VARIABLES: ['LATITUDE', 'PRESSURE', 'POTENTIAL_TEMPERATURE', 'LONGITUDE', 'SALINITY']
#OR:
    #MIMOC VARIABLES: ['LATITUDE', 'PRESSURE', 'CONSERVATIVE_TEMPERATURE', 'LONGITUDE', 'ABSOLUTE_SALINITY']
#Shape of salinity: (81, 341, 720) = depth, latitude, longitude

filename='Provide the path to your MIMOC file'

#import data:
data=import_data(filename)

#define the latitude-longitude range you want to look at:
llrange=['longitude1','longitude2','latitude1','latitude2']  #example: [120,135,30,40] captures the Korean Peninsula area

#find the relevant section of the data to look at based on your latitude/longitude range
ll_index=lat_lon_indices(data['LATITUDE'],data['LONGITUDE'],llrange)

fig=plt.figure()
m = Basemap(llcrnrlon=data['LONGITUDE'][ll_index[0]],llcrnrlat=data['LATITUDE'][ll_index[2]],
            urcrnrlon=data['LONGITUDE'][ll_index[1]],urcrnrlat=data['LATITUDE'][ll_index[3]],projection='mill')
x,y=np.meshgrid(data['LONGITUDE'][ll_index[0]:ll_index[1]],data['LATITUDE'][ll_index[2]:ll_index[3]])

#make a contour plot of the data you want to look at; you'll need to give it an index for the pressure, 
#which I think is spaced out by 5 meters near the surface. Therefore, if you pick 0 for plevel, you'll get 0m,
#1 for plevel is 5m, 2 for plevel is 10m, and so on. 
plevel='PICK A NUMBER CORRESPONDING TO THE PRESSURE YOU WANT TO LOOK AT'
print('You are looking at a depth of ', data['PRESSURE'][plevel])
m.contourf(x,y,data['NAME A VARIABLE'][plevel,ll_index[2]:ll_index[3],ll_index[0]:ll_index[1]],10,latlon='True')
plt.colorbar()
plt.show()
