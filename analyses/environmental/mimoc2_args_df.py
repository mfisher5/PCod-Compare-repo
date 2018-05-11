############## PYTHON CODE TO SUBSET AND PLOT MIMOC DATA ################
#
# Written by Rosalind Echols, UW Oceanography (rechols@uw.edu) 
#
#
# MF Edited 5/10/2018 to take command line arguments
#
###########################################################################

######## arguments #########
import argparse 

parser = argparse.ArgumentParser(description="subset a netCDF file from the MIMOC database, and print out data. THIS SHOULD ONLY BE RUN FOR LAT/LONG GRIDS WITH DATA COVERAGE. Use plotting script to ID grids to extract.")

parser.add_argument("-i", "--input", help="input file; should be a .nc.gz file from MIMOC")
parser.add_argument("-o", "--output", help="output text file for data frame; do not include `.txt`")
parser.add_argument("-lat", "--latitude", help="start, end latitude. 120,135 Korean peninsula [start,end]")
parser.add_argument("-long", "--longitude", help="start, end longitude. 30,40 Korean peninsula [start,end]")
parser.add_argument("-d", "--depth", help="depth to retrieve data from. Coded by seq(0,x,by=5) so d1 = 0m, d2 = 5m, d2 = 10m, etc. enter as the start and end of list ['0,3'] = '0,1,2,3")
parser.add_argument("-var", "--variable", help="Variable to plot [CONSERVATIVE_TEMPERATURE / ABSOLUTE_SALINITY]")

args = parser.parse_args()

########## packages ############

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import gsw
import cmocean


########### functions ##############
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

# define function that creates a range using floats
def frange(x, y, jump):
	while x < y:
	    yield x
	    x += jump


############ subset data ###############
print "subsetting data..."
#MIMOC VARIABLES: ['LATITUDE', 'PRESSURE', 'POTENTIAL_TEMPERATURE', 'LONGITUDE', 'SALINITY']
#OR:
    #MIMOC VARIABLES: ['LATITUDE', 'PRESSURE', 'CONSERVATIVE_TEMPERATURE', 'LONGITUDE', 'ABSOLUTE_SALINITY']
#Shape of salinity: (81, 341, 720) = depth, latitude, longitude

filename=args.input

#import data:
data=import_data(filename)


#define the latitude-longitude range you want to look at. example: [120,135,30,40] captures the Korean Peninsula area
lon_start = float(args.longitude.split(",")[0])
lon_end = float(args.longitude.split(",")[1])
lat_start = float(args.latitude.split(",")[0])
lat_end = float(args.latitude.split(",")[1])
llrange=[lat_start,lat_end,lon_start,lon_end]  
print llrange

#find the relevant section of the data to look at based on your latitude/longitude range
ll_index=lat_lon_indices(data['LATITUDE'],data['LONGITUDE'],llrange)


# subset data frame and save as new object "mydata"; you'll need to give it an index for the pressure, which I think is spaced out by 5 meters near the surface. 
#     Therefore, if you pick 0 for plevel, you'll get 0m,
#     1 for plevel is 5m, 2 for plevel is 10m, and so on. 
depth_start = int(args.depth.split(",")[0])
depth_stop = int(args.depth.split(",")[1])
depth_list = range(depth_start, depth_stop, 1)
outfile = open(args.output, "w")

# create and write file header (longitudes)
header_list = list(frange(lat_start, lat_end, 0.5))
header = "\t".join([str(i) for i in header_list])
outfile.write("lat\t" + header + "\tdepth\n")

for depth in depth_list:
	plevel=int(depth)
	print('You are looking at a depth of ', data['PRESSURE'][plevel])
	mydata = data[args.variable][plevel,ll_index[2]:ll_index[3],ll_index[0]:ll_index[1]]

	############ write data out to file ###########
	print "writing subset to file for depth ", depth, "..."

	# create list of latitudes
	row_list = list(frange(lon_start, lon_end, 0.5))

	# if-else to make sure row names match dimensions of data frame; 
	#     then write out data (one row per 0.5 degrees of latitude)
	if len(row_list) != len(mydata):
		print "Data and Lat/Long coordinates do not match!"
		print len(row_list)
		print row_list
		print len(mydata)
		print mydata
	else:
	    for i in range(0, len(row_list)):
	        outfile.write(str(row_list[i]) + "\t" + "\t".join([str(i) for i in mydata[i]]) + "\t" + str(data['PRESSURE'][plevel]) + "\n")


# close output file
outfile.close()


print "done."
