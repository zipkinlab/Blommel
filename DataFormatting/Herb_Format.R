#-----------------------#
#-Set Working directory-#
#-----------------------#

setwd("C:/Users/cblom/Documents/ZQE_Lab/HerbData")


#----------------#
#-Load libraries-#
#----------------#

library(dplyr)
library(tidyr)
library(rgdal)
library(sp)
library(raster)
library(sf)

#------------#
#-Import CSV-#
#------------#

raw <- read.csv("~/ZQE_Lab/HerbData/Herbivore Utilization Complete.csv", header=TRUE)
raw <- tbl_df(raw)

#Sort by species
raw <- arrange(raw, Animal)

#List species
levels(raw$Animal)

#---------------------#
#-Filter for migrants-#
#---------------------#

data <- filter(raw,
               Animal == "Buffalo" |
                 Animal == "Eland" |
                 Animal == "Elephant" |
                 Animal == "Giraffe" |
                 Animal == "Grants" |
                 Animal == "Hartebeest" |
                 Animal == "Hippo" |
                 Animal == "Impala" |
                 Animal == "Thomsons" |
                 Animal == "Topi" |
                 Animal == "Warthog" |
                 Animal == "Waterbuck" |
                 Animal == "Wildebeest" | 
                 Animal == "Zebra" |
                 Animal == "Cattle" |
                 Animal == "Shoat" |
                 Animal == "Lion" |
                 Animal == "Hyena" |
                 Animal == "BlackBackedJackal")

#Remove incomplete obs
data <- data[-which(is.na(data$Count)|data$Count<1),]
data
N <- length(data$Animal)


u1 <- data$AdjEasting
u2 <- data$AdjNorthing

#----------------------------------#
#Assign each DS point to a transect#
#----------------------------------#
 
#Directory for transects by site shapefile
d.dir <- "~/ZQE_Lab/HerbData/Shapefiles/"
setwd(d.dir)
#setwd for distance sampling shapefiles  
setwd("./DS")
#read in DS data 
DS <- st_read(dsn = ".", layer = "DS_10kmpersite") 
DS

#read in data as an sf object
data <- st_as_sf(data, coords = c("AdjEasting", "AdjNorthing"), crs = st_crs(DS))

#create distance matrix
ds_matrix <- st_distance(data, DS)
ds_matrix #visualize 
dim(ds_matrix) #23443 by 18, are there 18 transects? - CB

#assign transect to each observation 
ds_obs_transects <- apply(ds_matrix, 1, which.min)
length(ds_obs_transects) #23443, should be correct, matches number of obs 


#how to add transect back to the data set? convert back to dataframe? - CB
#creat observation array or assign distance classes first? - CB

#add transect (Site) for each observation back to the dataset 
data <- as.data.frame(data)
data$Site = ds_obs_transects
data$Distance_to_transect = min_dst

#-----------------------#
#Assign distance classes#
#-----------------------#

#ID for distance class
di <- seq(0,1000,25)

#Distance class
dclass <- rep(NA, N) #should this be rep 0? - CB

#Number of distance classes
nG <- length(di) - 1

#Minimum distance to assigned transect
min_dst <- apply(ds_matrix, 1, min)
length(min_dst)
dst <- data$Distance_to_transect

for(i in 1:N){
  for(k in 1:nG){
    if(di[k] < dst[i] && dst[i] <= di[k+1]) #why the k+1 argument? 
      dclass[i] <- k
      
  }
}

head(dclass)
length(dclass)
data$dclass = dclass

#--------------------#
#DS Observation array#
#--------------------#

length(unique(data$Animal)) #19 species? - CB
length(unique(data$Site)) #18 sites? - CB
length(unique(data$Month)) #11 months, monthly replicates? - CB
#how many replicates? 16? Does replicate column exist in the data set?- CB

y <- array(rep(0), dim = c(11,18,19))

dim(data)[1]

#suggested loop from Matt
#for(i in dim(data)[1]){ 
#  y[rep[i], site[i], spec[i]] <- y[i,i,i] + data$count[i] 
#}



#----------------------------------#
#Assign each TC point to a transect#
#----------------------------------#

#setwd for transect count shapefiles 
#d.dir <- "~/ZQE_Lab/HerbData/Shapefiles/"
#setwd(d.dir)
#setwd("./Transects")
#read in transect count data
#TC <- st_read(dsn = ".", layer = "Transects")  #TC for transect count
#TC

#back a directory and read in TC data
#setwd("C:/Users/cblom/Documents/ZQE_Lab/HerbData")

#read in TC data as an sf object 
#TCdata <- read.csv("~/ZQE_Lab/HerbData/tblPreyCensus_2012to2014.csv", header=TRUE)
#length(unique(TCdata$transect))
#TCdata <- st_as_sf(TCdata, coords = c(), crs = st_crs(TC))

#---------------#
#-Simulate Data-#
#---------------#

#Remove obs over 1000 #cutoff in distance samplng, detection gets too low 
#data$dclass <- dclass
#data$site <- site
data <- data[-(which(dst>1000)),]
u1 <- u1[-(which(dst>1000))]
u2 <- u2[-(which(dst>1000))]

#Sample ID
sample_ID <- rep(NA, length(data[,1]))
data <- data.frame(data, sample_ID)

data

data$sample_ID[data$Territory == "North"] <- filter(data,Territory=="North")%>%group_by(Year, Month)%>%group_indices()
data$sample_ID[data$Territory == "South"] <- filter(data,Territory=="South")%>%group_by(Year, Month)%>%group_indices()
data$sample_ID[data$Territory == "West"] <- filter(data,Territory=="West")%>%group_by(Year, Month)%>%group_indices()


#Filter out unnecessary data
D <- data.frame(data$Animal, data$Site, data$sample_ID, data$Count, data$dclass)
colnames(D) <- c("Animal", "Site_ID", "Sample_ID", "Count", "dclass")

#Vector of species
name <- c("Buffalo", "Eland", "Elephant", "Giraffe", "Grants", "Hartebeest", "Hippo", "Impala", "Thomsons", "Topi", "Warthog", "Waterbuck", "Wildebeest", "Zebra")

H <- D %>%filter(Animal %in% name)


#Initialize observation array (rep x site x species)
#create arracy with rep, number of sites, number of species 
# could probably replace it with 
#not sure about the dimensions of this - CB
length(unique(H$Sample_ID))
length(unique(H$Site_ID)) #18 sites?- CB
length(unique(H$Animal))
y <- array(NA, dim = c(16,17,14))

#Generate observation array
for(s in 1:14){ #for each species
  A <- (filter(H, Animal == name[s])) #why is this necessary? 
  A <- group_by(A, Site_ID, Sample_ID, Animal)%>%summarize(n()) #summarized data within one species 
  W <- data.frame(rep(1:17, rep(16, 17)), rep(1:16, 17)) #done to get around the problem of inequal sampling in different regions 
  colnames(W) <- c("Site_ID", "Sample_ID")
  B <- full_join(W, A, by.x = c("Site_ID", "Sample_ID"), by.y = c("Site_ID", "Sample_ID"))
  B$`n()`[is.na(B$`n()`)] = 0
  C <- split(B$`n()`, f = B$Site_ID) #not sure what this does either 
  C <- do.call(cbind, C)
  for(j in 14:17){
    for(t in 14:16){
      C[t,j] = NA
    }
  }
  y[,,s] <- C
} #getting consistent indexing errors/length mismatches - CB

#-------------------------#
#-Create distance classes-#
#-------------------------#

#Width of distance classes
v <- 25 #meters

#Transect half-width
B <- 1000 #meters

#Site ID
site <- H$Site_ID

#Replicate ID
rep <- H$Sample_ID

#Species ID
spec <- as.integer(droplevels(H$Animal))

#Distance class midpoint ID
mdpt <- seq(12.5, 1000, 25)

#Number of sites
nsites <- 17

#Number of reps
nreps <- as.vector(c(rep(16, 13), rep(13, 4)))

#Number of speceis
nherb <- length(unique(spec))

#Number of distance classes
nD <- length(mdpt)

#Group size
gs <- H$Count

#Dclass
dclass <- H$dclass

#----------------------------#
#-Offset for transect length-# 
#----------------------------#

offset <- as.vector(c(1, 1, 1, 1, 1, 1, 1, 1.080,
                      0.878, 1, 1, 1, 1, 1, 1.100,
                      1.300, 1.237)) #not all transects are the same length 

#-------------------------#
#-Create Region Covariate-#
#-------------------------#

region <- c(rep(0, 13), rep(1, 4))

#-----------------------#
#-Create NDVI Covariate-#  
#-----------------------#

#stopped here, I think I need more data

NDVIdata <- read.csv("RawData/NDVI.csv", header = FALSE) #Normalized Difference Vegetation Index, essentially represents greenness - how much vegetation there is
colnames(NDVIdata) <- c("NDVI", "Year", "Day", "Site")

NDVIdata$Date <- as.Date(NDVIdata$Day, origin = paste(NDVIdata$Year,"-01-01", sep = ""))


Date <- as.data.frame(as.Date(with(filter(data, Animal %in% name), paste(Year, Month, Day, sep = "-")), "%Y-%m-%d"))
Date$site <- site
Date$rep <- rep
colnames(Date) <- c("date","site","rep")

NDVI <- array(NA, dim = c(16,17))

#this is a way to 'match' NDVI to the sampling data we have (spatially and temporally) 
for(j in 1:nsites){
  for(t in 1:nreps[j]){
    tmp1 <- filter(Date, site==j&rep==t)%>%summarize(min(date)+1)
    tmp2 <- as.vector(filter(NDVIdata, Site==j)%>%select(Date))
    tmp3 <- which.min(abs(as.numeric(tmp2[,]-tmp1[,]))) #what does this do?  
    NDVI[t,j] <- NDVIdata[NDVIdata$Site==j&NDVIdata$Date==tmp2[tmp3,],1] 
  }
}

NDVI <- (NDVI - mean(NDVI, na.rm = TRUE))/sd(NDVI, na.rm = TRUE)

#-------------------------------------------#
#-Create Livestock and Carnivore Covariates-#
#-------------------------------------------#

#Vector of species
name <- c("BlackBackedJackal", "Cattle", "Hyena", "Lion", "Shoat")

H <- D %>%filter(Animal %in% name)

#Initialize observation array (rep x site x species)
covar <- array(NA, dim = c(16,17,5))

#Generate observation array 
#CB: I need to review why/how this is done 
for(s in 1:5){
  A <- (filter(H, Animal == name[s]))
  A <- group_by(A, Site_ID, Sample_ID, Animal)%>%summarize(n())
  W <- data.frame(rep(1:17, rep(16, 17)), rep(1:16, 17))
  colnames(W) <- c("Site_ID", "Sample_ID")
  Q <- full_join(W, A, by.x = c("Site_ID", "Sample_ID"), by.y = c("Site_ID", "Sample_ID"))
  Q$`n()`[is.na(Q$`n()`)] = 0
  C <- split(Q$`n()`, f = Q$Site_ID)
  C <- do.call(cbind, C)
  for(j in 14:17){
    for(t in 14:16){
      C[t,j] = NA
    }
  }
  covar[,,s] <- C
}

#Site ID
site <- c(site, H$Site_ID)

#Replicate ID
rep <- c(rep, H$Sample_ID)

#Species ID
spec <- c(spec, 14 + as.integer(droplevels(H$Animal)))

#Number of observations
nobs <- length(spec)

#Group size
gs <- c(gs, H$Count)

#Dclass
dclass <- c(dclass, H$dclass)

#Number of total species
nspec <- 19
#standard error? 
BBJ <- (Carnivore$mean$GSrep[,,1] - 
          mean(Carnivore$mean$GSrep[,,1], na.rm = TRUE))/sd(Carnivore$mean$GSrep[,,1], na.rm = TRUE) #what is GSrep? 
Hyena <- (Carnivore$mean$GSrep[,,2] - 
            mean(Carnivore$mean$GSrep[,,2], na.rm = TRUE))/sd(Carnivore$mean$GSrep[,,2], na.rm = TRUE)
Lion <- (Carnivore$mean$GSrep[,,3] - 
           mean(Carnivore$mean$GSrep[,,3], na.rm = TRUE))/sd(Carnivore$mean$GSrep[,,3], na.rm = TRUE)
Cattle <- cbind(matrix(0, nrow = 16, ncol = 13), 
                rbind(Livestock$mean$GSrep[,,1], 
                      matrix(NA, nrow = 3, ncol = 4)))
Cattle <- (Cattle - mean(Cattle, na.rm = TRUE))/sd(Cattle, na.rm = TRUE)
Shoat <- cbind(matrix(0, nrow = 16, ncol = 13), 
                rbind(Livestock$mean$GSrep[,,2], 
                      matrix(NA, nrow = 3, ncol = 4)))
Shoat <- (Shoat - mean(Shoat, na.rm = TRUE))/sd(Shoat, na.rm = TRUE)

#---------------------#
#-Distance Covariates-#
#---------------------#

Dstdata <- read.csv("RawData/Dst.csv", header = FALSE)
colnames(Dstdata) <- c("RiverDst", "BoundDst")

riverDst <- Dstdata$RiverDst #distance from river 
riverDst <- (riverDst - mean(riverDst, na.rm = TRUE))/sd(riverDst, na.rm = TRUE)

boundDst <- Dstdata$BoundDst #distance from the border of the park or region? 
boundDst <- (boundDst - mean(boundDst, na.rm = TRUE))/sd(boundDst, na.rm = TRUE)

#----------------#
#-LULC Covariate-# #Land Use Land Cover 
#----------------#

files <- list.files(path = "~/GitHub/Herbivore/GIS/SiteLULC", pattern = "tif$", full.names = TRUE)
LULCstack <- list()
for(i in 1:length(files)){
  LULCstack[[i]] <- raster(files[i])
}

LULC <- matrix(0, ncol = 7, nrow = 17)
for(i in 1:17){
  for(j in 1:7){
    LULC[i,j] <- sum(values(LULCstack[[i]])==j, na.rm = TRUE)/sum(!is.na(values(LULCstack[[i]])))
  }
}
LULC <- scale(LULC)
LULC <- LULC[,-6]

#---------------------#
#-Migration Covariate-#
#---------------------#

Migrant <- read.csv("./RawData/Migration.csv", header = FALSE)
Migrant <- as.matrix(Migrant)
Migrant1 <- array(0, dim = c(16,17,14))
Migrant1[,,13:14] <- Migrant
Migrant <- Migrant1
Migrant[14:16,14:17,] <- NA

#--------------#
#-Combine data-#
#--------------#

DSdata <- list(y, dclass, v, B, site, rep, spec, mdpt, nsites, nreps, nherb, nspec, nobs, 
               nD, gs, region, offset, NDVI, covar, riverDst, boundDst, LULC, Migrant, 
               BBJ, Hyena, Lion, Cattle, Shoat)
heads <- c("y", "dclass", "v", "B", "site", "rep", "spec", "mdpt", "nsites", "nreps", 
           "nherb", "nspec", "nobs", "nD", "gs", "region", "offset", "NDVI", "covar",
           "riverDst", "boundDst", "LULC", "Migrant", "BBJ", "Hyena", "Lion", "Cattle",
           "Shoat")
DSdata <- setNames(DSdata, nm = heads)

#-------------#
#-Export data-#
#-------------#

save(DSdata, file = "C:/Users/farrm/Documents/GitHub/Herbivore/herbdata.R")

#----------------#
#-Visualize data-#
#----------------#

#Easting
xlim <- c(715304, 752393)

#Northing
ylim <- c(9831970, 9857296)

#Visualize
#could we do this in a loop? Especially with simple features 
#maybe don't try to adapt this, just start fresh 
plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, 
     yaxt = "n", xaxt = "n", ylab = "", xlab = "")
plot(Site1, add=T, col="red")
plot(Site2, add=T, col="darkorange")
plot(Site3, add=T, col="gold")
plot(Site4, add=T, col="darkolivegreen2")
plot(Site5, add=T, col="forestgreen")
plot(Site6, add=T, col="aquamarine")
plot(Site7, add=T, col="cadetblue3")
plot(Site8, add=T, col="cornflowerblue")
plot(Site9, add=T, col="blue")
plot(Site10, add=T, col="blueviolet")
plot(Site11, add=T, col="darkmagenta")
plot(Site12, add=T, col="deeppink4")
plot(Site13, add=T, col="darkred")
plot(Site14, add=T, col="black")
plot(Site15, add=T, col="coral2")
plot(Site16, add=T, col="brown")
plot(Site17, add=T, col="grey")

points(cbind(u1[site==1], u2[site==1]), col="red", lwd = 1)
points(cbind(u1[site==2], u2[site==2]),  col="darkorange", lwd = 1)
points(cbind(u1[site==3], u2[site==3]), col="gold", lwd = 1)
points(cbind(u1[site==4], u2[site==4]), col="darkolivegreen2", lwd = 1)
points(cbind(u1[site==5], u2[site==5]), col="forestgreen", lwd = 1)
points(cbind(u1[site==6], u2[site==6]), col="aquamarine", lwd = 1)
points(cbind(u1[site==7], u2[site==7]), col="cadetblue3", lwd = 1)
points(cbind(u1[site==8], u2[site==8]), col="cornflowerblue", lwd = 1)
points(cbind(u1[site==9], u2[site==9]), col="blue", lwd = 1)
points(cbind(u1[site==10], u2[site==10]), col="blueviolet", lwd = 1)
points(cbind(u1[site==11], u2[site==11]), col="darkmagenta", lwd = 1)
points(cbind(u1[site==12], u2[site==12]), col="deeppink4", lwd = 1)
points(cbind(u1[site==13], u2[site==13]), col="darkred", lwd = 1)
points(cbind(u1[site==14], u2[site==14]), col="black", lwd = 1)
points(cbind(u1[site==15], u2[site==15]), col="coral2", lwd = 1)
points(cbind(u1[site==16], u2[site==16]), col="brown", lwd = 1)
points(cbind(u1[site==17], u2[site==17]), col="grey", lwd = 1)

plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, 
     yaxt = "n", xaxt = "n", ylab = "", xlab = "")
plot(Site1, add=T, col="red")
plot(Site2, add=T, col="darkorange")
plot(Site3, add=T, col="gold")
plot(Site4, add=T, col="darkolivegreen2")
plot(Site5, add=T, col="forestgreen")
plot(Site6, add=T, col="aquamarine")
plot(Site7, add=T, col="cadetblue3")
plot(Site8, add=T, col="cornflowerblue")
plot(Site9, add=T, col="blue")
plot(Site10, add=T, col="blueviolet")
plot(Site11, add=T, col="darkmagenta")
plot(Site12, add=T, col="deeppink4")
plot(Site13, add=T, col="darkred")
plot(Site14, add=T, col="black")
plot(Site15, add=T, col="coral2")
plot(Site16, add=T, col="brown")
plot(Site17, add=T, col="grey")


points(cbind(u1[data$Animal=="Buffalo"], u2[data$Animal=="Buffalo"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Eland"], u2[data$Animal=="Eland"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Elephant"], u2[data$Animal=="Elephant"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Giraffe"], u2[data$Animal=="Giraffe"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Grants"], u2[data$Animal=="Grants"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Hartebeest"], u2[data$Animal=="Hartebeest"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Hippo"], u2[data$Animal=="Hippo"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Impala"], u2[data$Animal=="Impala"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Topi"], u2[data$Animal=="Topi"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Waterbuck"], u2[data$Animal=="Waterbuck"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Zebra"], u2[data$Animal=="Zebra"]), col="red", lwd = 1)
points(cbind(u1[data$Animal=="Wildebeest"], u2[data$Animal=="Wildebeest"]), col="black", lwd = 1)
points(cbind(u1[data$Animal=="Thomsons"], u2[data$Animal=="Thomsons"]), col="grey", lwd = 1)

#--------------------------#
#-Visualize covariate data-#
#--------------------------#

#1 Buffalo
plot(y[,,1]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,1])~as.vector(NDVI)))
#2 Elend
plot(y[,,2]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,2])~as.vector(NDVI)))
#3 Elephant
plot(y[,,3]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,3])~as.vector(NDVI)))
#4 Giraffe
plot(y[,,4]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,4])~as.vector(NDVI)))
#5 Grant's
plot(y[,,5]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,5])~as.vector(NDVI)))
#6 Hartebeest
plot(y[,,6]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,6])~as.vector(NDVI)))
#7 Hippo
plot(y[,,7]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,7])~as.vector(NDVI)))
#8 Impala
plot(y[,,8]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,8])~as.vector(NDVI)))
#9 Thomsons
plot(y[,,9]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,9])~as.vector(NDVI)))
#10 Topi
plot(y[,,10]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,10])~as.vector(NDVI)))
#11 Warthog
plot(y[,,11]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,11])~as.vector(NDVI)))
#12 Waterbuck
plot(y[,,12]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,12])~as.vector(NDVI)))
#13 Wildebeest
plot(y[,,13]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,13])~as.vector(NDVI)))
#14 Zebra
plot(y[,,14]~NDVI, xlim = c(-1,1))
abline(lm(as.vector(y[,,14])~as.vector(NDVI)))

#1 Buffalo
plot(y[,,1]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,1])~as.vector(cov[,,1])))
#2 Elend
plot(y[,,2]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,2])~as.vector(cov[,,1])))
#3 Elephant
plot(y[,,3]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,3])~as.vector(cov[,,1])))
#4 Giraffe
plot(y[,,4]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,4])~as.vector(cov[,,1])))
#5 Grant's
plot(y[,,5]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,5])~as.vector(cov[,,1])))
#6 Hartebeest
plot(y[,,6]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,6])~as.vector(cov[,,1])))
#7 Hippo
plot(y[,,7]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,7])~as.vector(cov[,,1])))
#8 Impala
plot(y[,,8]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,8])~as.vector(cov[,,1])))
#9 Thomsons
plot(y[,,9]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,9])~as.vector(cov[,,1])))
#10 Topi
plot(y[,,10]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,10])~as.vector(cov[,,1])))
#11 Warthog
plot(y[,,11]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,11])~as.vector(cov[,,1])))
#12 Waterbuck
plot(y[,,12]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,12])~as.vector(cov[,,1])))
#13 Wildebeest
plot(y[,,13]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,13])~as.vector(cov[,,1])))
#14 Zebra
plot(y[,,14]~cov[,,1], xlim = c(0,10))
abline(lm(as.vector(y[,,14])~as.vector(cov[,,1])))

#1 Buffalo
plot(y[,,1]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,1])~as.vector(cov[,,2])))
#2 Elend
plot(y[,,2]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,2])~as.vector(cov[,,2])))
#3 Elephant
plot(y[,,3]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,3])~as.vector(cov[,,2])))
#4 Giraffe
plot(y[,,4]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,4])~as.vector(cov[,,2])))
#5 Grant's
plot(y[,,5]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,5])~as.vector(cov[,,2])))
#6 Hartebeest
plot(y[,,6]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,6])~as.vector(cov[,,2])))
#7 Hippo
plot(y[,,7]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,7])~as.vector(cov[,,2])))
#8 Impala
plot(y[,,8]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,8])~as.vector(cov[,,2])))
#9 Thomsons
plot(y[,,9]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,9])~as.vector(cov[,,2])))
#10 Topi
plot(y[,,10]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,10])~as.vector(cov[,,2])))
#11 Warthog
plot(y[,,11]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,11])~as.vector(cov[,,2])))
#12 Waterbuck
plot(y[,,12]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,12])~as.vector(cov[,,2])))
#13 Wildebeest
plot(y[,,13]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,13])~as.vector(cov[,,2])))
#14 Zebra
plot(y[,,14]~cov[,,2], xlim = c(0,10))
abline(lm(as.vector(y[,,14])~as.vector(cov[,,2])))

#1 Buffalo
plot(y[,,1]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,1])~as.vector(cov[,,3])))
#2 Elend
plot(y[,,2]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,2])~as.vector(cov[,,3])))
#3 Elephant
plot(y[,,3]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,3])~as.vector(cov[,,3])))
#4 Giraffe
plot(y[,,4]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,4])~as.vector(cov[,,3])))
#5 Grant's
plot(y[,,5]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,5])~as.vector(cov[,,3])))
#6 Hartebeest
plot(y[,,6]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,6])~as.vector(cov[,,3])))
#7 Hippo
plot(y[,,7]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,7])~as.vector(cov[,,3])))
#8 Impala
plot(y[,,8]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,8])~as.vector(cov[,,3])))
#9 Thomsons
plot(y[,,9]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,9])~as.vector(cov[,,3])))
#10 Topi
plot(y[,,10]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,10])~as.vector(cov[,,3])))
#11 Warthog
plot(y[,,11]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,11])~as.vector(cov[,,3])))
#12 Waterbuck
plot(y[,,12]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,12])~as.vector(cov[,,3])))
#13 Wildebeest
plot(y[,,13]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,13])~as.vector(cov[,,3])))
#14 Zebra
plot(y[,,14]~cov[,,3], xlim = c(0,10))
abline(lm(as.vector(y[,,14])~as.vector(cov[,,3])))