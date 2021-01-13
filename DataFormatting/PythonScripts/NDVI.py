#----------------------------------------------#
#----Management and extraction of NDVI data----#
#----from MODIS via USGS/NASA.-----------------#
#----Created by Matthew Farr-------------------#
#----------------------------------------------#
import arcpy,shutil,os,csv
from arcpy.sa import *

#Working environment within HDF folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\HDF\\"
#List of HDF files
rasterlist = arcpy.ListRasters("*" "HDF")
#UTM projection
projection = arcpy.SpatialReference(4326)
#Extract NDVI from HDF files
for raster in rasterlist:
	name = raster.split('.')
	filename = "..\\NDVI\\NDVI"+str(name[1][1:]) + ".tif"
	arcpy.ExtractSubDataset_management(raster, filename, "0")

#Working environment within site folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\Site\\"
#List of shapefiles(i.e., sites)
shlist = arcpy.ListFeatureClasses()
#Buffer sites 650 meters
for fc in shlist:
	name=fc.split('.')
	filename = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\SiteBuffer\\" + name[0] + ".shp"
	arcpy.Buffer_analysis(fc, filename, "1000 meters", "FULL", "FLAT")

#Working environment within sitebuffer folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\SiteBuffer\\"
#List of shapefiles(i.e., site buffers)
shlist = arcpy.ListFeatureClasses()
#Working environment within NDVI folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\NDVI\\"
#List of TIF files (i.e., NDVI)
rasterlist = arcpy.ListRasters("*" "TIF")
#Working environment within GIS folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\"
#Clipping each NDVI by the site buffer
for raster in rasterlist:
        raster1 = "NDVI\\" + raster
        name1 = raster.split(".")
	for fc in shlist:
		clip = "SiteBuffer\\" + fc
		name2 = fc.split(".")
		filename = "SiteNDVI\\" + name1[0] + name2[0] + ".tif"
		arcpy.Clip_management(raster1, "#", filename, clip, "", "ClippingGeometry", "")

#Working environment within SiteNDVI
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\SiteNDVI\\"
#List of TIF files
rasterlist = arcpy.ListRasters("*" "TIF")
#Calculate mean NDVI for site and convert to NDVI scale (i.e., -1 to 1)
NDVI = []
Year = []
Day = []
Site = []
for raster in rasterlist:
    means = arcpy.GetRasterProperties_management(raster, "Mean")
    NDVI.append(str(float(means[0])*0.0001))
    names = raster.split("Site")
    Year.append(names[0][4:8])
    Day.append(names[0][8:])
    Site.append(names[1].split(".")[0])

#Write to file
myfile = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\RawData\\NDVI.csv"
with open(myfile, 'w') as f:
    writer = csv.writer(f, delimiter=',',lineterminator = '\n')
    writer.writerows(zip(NDVI,Year,Day,Site))
quit














	
	
	
	
	
