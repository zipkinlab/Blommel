#-------------------------------------------#
#----Extract by mask LULC to site buffer----#
#----Created by Matthew Farr----------------#
#-------------------------------------------#

import arcpy
from arcpy.sa import *

#Working environment
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Herbivore\\GIS\\SiteBuffer\\"
#Activate extensions
arcpy.CheckOutExtension("Spatial")
#List of shapefiles(i.e., sites)
shlist = arcpy.ListFeatureClasses()
#Raster file
raster = "C:\\Users\\farrm\\Documents\\GitHub\\Herbivore\\GIS\\LULC\\kenya_lsat_20090612_palsar_20090613_hh_hv__pca1_rf.tif"
#Calculate centroid for each site
for fc in shlist:
    name = fc.split('.')
    filename = "C:\\Users\\farrm\\Documents\\GitHub\\Herbivore\\GIS\\SiteLULC\\" + name[0] + "LULC.tif"
    outExtractByMask = ExtractByMask(raster, fc)
    outExtractByMask.save(filename)
