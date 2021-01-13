#-------------------------------------------------#
#----Extraction of distance to nearest feature----#
#----for each site. Creates Dst2Water and---------#
#----Dst2Bound.-----------------------------------#
#----Created by Matthew Farr----------------------#
#-------------------------------------------------#
import arcpy,shutil,os,csv
from arcpy.sa import *

#Working environment within SiteBuffer folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\SiteBuffer\\"
#List of shapefiles(i.e., sites)
shlist = arcpy.ListFeatureClasses()
#Calculate centroid for each site
for fc in shlist:
    name = fc.split('.')
    filename = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\Centroid\\" + name[0] + "centroid.shp"
    arcpy.FeatureToPoint_management(fc, filename, "CENTROID")

#Working environment within SiteBuffer folder
arcpy.env.workspace = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\GIS\\Centroid\\"
#List of shapefiles(i.e., sites)
shlist = arcpy.ListFeatureClasses()
#Rivers and boundary shapefiles
rivers = "..\\DstFeatures\\MaraRivers.shp"
boundary = "..\\DstFeatures\\MMNR.shp"
#Field name to extract
field = "NEAR_DIST"
#Distance values
RiverDst = []
BoundDst = []
#Calculate nearest distnace to rivers and boundary
for fc in shlist:
    arcpy.Near_analysis(fc, rivers)
    cursor = arcpy.SearchCursor(fc)
    for row in cursor:
        RiverDst.append(row.getValue(field))
    arcpy.Near_analysis(fc, boundary)
    cursor = arcpy.SearchCursor(fc)
    for row in cursor:
        BoundDst.append(row.getValue(field))

#Write to file
myfile = "C:\\Users\\farrm\\Documents\\GitHub\\Ungulate\\RawData\\Dst.csv"
with open(myfile, 'w') as f:
    writer = csv.writer(f, delimiter=',',lineterminator = '\n')
    writer.writerows(zip(RiverDst, BoundDst))
quit
