#-----------------------#
#-Set Working directory-#
#-----------------------#

#MTF: for version controlling you should make sure the local directory
#folders match the remote repository on github. I would also recommend
#using R projects. It seets the path root to your project directory.
#prevents you from having to change directories

#setwd("C:/Users/cblom/Documents/ZQE_Lab/HerbData") #CB
#delete - test

#----------------#
#-Load libraries-#
#----------------#

library(tidyverse)
# library(raster) I don't (anymore) add the raster library because it interfers with tidyverse functions. Instead I call raster functions using raster::
library(sf)

#------------------------#
#-Distance sampling data-#
#------------------------#

DS <- read.csv("~/ZQE_Lab/HerbData/Herbivore Utilization Complete.csv", header=TRUE) #CB
DS <- read.csv("./RawData/Herbivore Utilization Complete.csv") #MTF

#-------------#
#-Filter data-#
#-------------#

DS$Animal <- as.factor(DS$Animal)

DS <- filter(DS,
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
             Animal == "BlackBackedJackal") %>% droplevels()

DS$Animal <- factor(DS$Animal,
                    levels = c("Buffalo", "Eland", "Elephant",
                               "Giraffe", "Grants", "Hartebeest",
                               "Hippo", "Impala", "Thomsons",
                               "Topi", "Warthog", "Waterbuck",
                               "Wildebeest","Zebra", "Cattle", 
                               "Shoat", "Lion", "Hyena", "BlackBackedJackal"))

#Remove incomplete obs
DS <- DS[-which(is.na(DS$Count)|DS$Count<1),]

#-------------------------------#
#-Format distance sampling data-#
#-------------------------------#

#read in DS data 
DSshape <- st_read(dsn = "./Shapefiles", layer = "DS_10kmpersite") #CB
DSshape <- st_read(dsn = "./RawData/Shapefiles", layer = "DS_10kmpersite") #MTF

#read in data as an sf object
DS <- st_as_sf(DS, coords = c("AdjEasting", "AdjNorthing"), crs = st_crs(DSshape))

#create distance matrix
ds_matrix <- st_distance(DS, DSshape)

#assign transect to each observation 
DS$site <- apply(ds_matrix, 1, which.min)

#add corrected distances
DS$dst <- apply(ds_matrix, 1, min)
DS <- DS %>% filter(dst <= 1000)

#-----------------------#
#Assign distance classes#
#-----------------------#

#Number of observations for distance sampling
nobs <- NULL
nobs[1] <- dim(DS)[1]

#ID for distance class
di <- seq(0,1000,25)

#Distance class
dclass <- rep(NA, nobs[1])

#Number of distance classes
nG <- length(di) - 1

#Minimum distance to assigned transect
dst <- DS$dst

for(i in 1:nobs[1]){
  for(k in 1:nG){
    if(di[k] < dst[i] && dst[i] <= di[k+1]) #why the k+1 argument? 
      dclass[i] <- k
    
  }
}

DS$dclass <- dclass

#Replicate
DS$reps[DS$Territory == "North"] <- filter(DS,Territory=="North")%>%group_by(Year, Month)%>%group_indices()
DS$reps[DS$Territory == "South"] <- filter(DS,Territory=="South")%>%group_by(Year, Month)%>%group_indices()
DS$reps[DS$Territory == "West"] <- filter(DS,Territory=="West")%>%group_by(Year, Month)%>%group_indices()

#Filter out unnecessary data
DS <- DS %>% select(Animal, site, reps, Count, dclass, Day, Month, Year, Territory)

DS$spec <- as.numeric(DS$Animal)

#---------------#
#-Transect data-#
#---------------#

#read in TC data
TC <- read.csv("~/ZQE_Lab/HerbData/tblPreyCensus_2012to2014.csv", header=TRUE) #CB
TC <- read.csv("./RawData/tblPreyCensus_2012to2014.csv") #MTF

#read in transect count shapefiles
TCshape <- st_read(dsn = "./RawData/Shapefiles", layer = "Transects")

#reshape TC to long format
TC <- reshape2::melt(TC, id=colnames(TC)[1:8])

#rename columns
TC <- TC %>% rename("Animal" = "variable", "Count" = "value")

#filter species
TC <- filter(TC,
               Animal == "buffalo" |
               Animal == "eland" |
               Animal == "elephant" |
               Animal == "giraffe" |
               Animal == "grants" |
               Animal == "hartebeest" |
               Animal == "hippo" |
               Animal == "impala" |
               Animal == "thomsons" |
               Animal == "topi" |
               Animal == "warthog" |
               Animal == "waterbuck" |
               Animal == "wildebeest" | 
               Animal == "zebra" |
               Animal == "cow" |
               Animal == "shoat" |
               Animal == "lion" |
               Animal == "spotted_hyena" |
               Animal == "bb_jackal") %>% droplevels()

#set transects to factor
TC$transect <- as.factor(TC$transect)

#filter transects
TC <- filter(TC,
               transect == "WLOW"|
               transect == "WHIGH"|
               transect == "W3"|
               transect == "North"|
               transect == "S1"|
               transect == "S2"|
               transect == "RSP"|
               transect == "SST"|
               transect == "HZT") %>% droplevels()

#select relevant columns
TC <- TC %>% select(transect, date, Animal, Count)

#split date to month, day, year
TC <- separate(data = TC, col = date, into = c("month", "day", "year"), sep = "/")

#set columns to numeric values
TC$year <- as.numeric(TC$year)
TC$month <- as.numeric(TC$month)
TC$day <- as.numeric(TC$day)
TC$Count <- as.numeric(TC$Count)

#determine 1st or 2nd replicate of the month
TC <- TC %>% group_by(Animal, transect, month, year) %>%
  mutate(min_day = min(day),
         max_day = max(day),
         n_day = n_distinct(day)) %>%
  filter(!(n_day == 3 & max_day == day)) %>%
  mutate(k = ifelse(n_day == 1 & day < 16, 1,
                    ifelse(n_day == 2 & day == min_day, 1,
                           ifelse(n_day == 3 & day == min_day, 1,0)))) %>%
  select(transect, month, day, year, Animal, Count, k)

#function to generate replicates
reps <- function(month, day, year, k){
  yr <- year - 2011
  reps <- ((month * 2) - k) + yr * 24 - 24
  return(reps)
}

#generate replicate, site, and species IDs
TC <- TC %>% mutate(reps = reps(month, day, year, k),
                    site = as.numeric(transect),
                    spec = as.numeric(Animal))

#--------------------#
#-Observation arrays-#
#--------------------#

#number of sites
nsites <- NULL
nsites[1] <- max(DS$site) #distance sampling
nsites[2] <- max(TC$site) #transect counts

#adjust transect counts site ID
TC <- TC %>% mutate(site = site + nsites[1])

#number of replicates
nreps <- NULL
nreps[1] <- max(DS$reps) #distance sampling
nreps[2] <- max(TC$reps) #transect counts

#start and end replicate for each site
nstart <- c(rep(1,nsites[1]), rep(nreps[1]+1,nsites[2]))
nend <- c(DS %>% group_by(site) %>% summarise(nend = max(reps)) %>% select(nend) %>% .$nend, rep(nreps[1]+nreps[2],nsites[2]))

#adjust transect counts replicate ID
TC <- TC %>% mutate(reps = reps + nreps[1])

#number of species
nspec <- 14

#observation array
y <- array(NA, dim = c(sum(nreps), sum(nsites), nspec))

#set distance sampling counts to 0 as sampling always occuried 
y[1:nreps[1],1:nsites[1],] <- 0

for(i in 1:dim(DS)[1]){
  if(DS$spec[i] > 14){
    next
  }
  y[DS$reps[i],DS$site[i],DS$spec[i]] <- as.numeric(DS$Count[i]) + y[DS$reps[i],DS$site[i],DS$spec[i]]
}

for(i in 1:dim(TC)[1]){
  if(TC$spec[i] > 14){
    next
  }
  y[TC$reps[i],TC$site[i],TC$spec[i]] <- as.numeric(TC$Count[i])
}

#-------------------#
#-Remaining indices-#
#-------------------#

#Site ID
site <- DS %>% filter(spec < 15) %>% select(site) %>% .$site

#Spec ID
spec <- DS %>% filter(spec < 15) %>% select(spec) %>% .$spec

#Distance class
dclass <- DS %>% filter(spec < 15) %>% select(dclass) %>% .$dclass

#Number of observations
nobs <- length(site)

#Width of distance classes
v <- 25 #meters

#Transect half-width
B <- 1000 #meters

#Distance class midpoint ID
mdpt <- seq(v/2, B, v)

#---------------------------------#
#-Area offset for transect length-# 
#---------------------------------#

#area of transects (m^2)
area <- as.numeric(c(st_length(DSshape)*1000*2, st_length(TCshape)*200*2))

#set baseline unit as 1 km^2
offset <- area/1E6

#------------#
#-Covariates-#
#------------#

#-Region-#

region <- c(DS %>% group_by(site) %>% summarize(region = max(Territory)) %>% mutate(region = ifelse(region == "West", 1, 0)) %>% select(region) %>% .$region,
            TC %>% group_by(site) %>% summarize(region = unique(transect)) %>% mutate(region = ifelse(region == "W3"|region == "WHIGH"|region == "WLOW", 1, 0)) %>% select(region) %>% .$region)
  
#-Migration-#

migration <- c(DS %>% mutate(migration = ifelse(Month %in% 7:11, 1, 0)) %>% group_by(reps) %>% summarise(migration = max(migration)) %>% select(migration) %>% .$migration,
            TC %>% mutate(migration = ifelse(month %in% 7:11, 1, 0)) %>% group_by(reps) %>% summarise(migration = max(migration)) %>% select(migration) %>% .$migration)

tmp <- array(0, dim = dim(y)[-2])
tmp[,9] <- migration
tmp[,13] <- migration
tmp[,14] <- migration
migration <- tmp

#-NDVI-#

#-Carnivores-#

#-Livestock-#

#-LULC-#

#--------------#
#-Compile data-#
#--------------#

Data <- list(y, dclass, v, B, mdpt, nG, nobs,
             nreps, nstart, nend, nsites, nspec, site, spec,
             offset, region, migration)

heads <- c("y", "dclass", "v", "B", "mdpt", "nG", "nobs",
           "nreps", "nstart", "nend", "nsites", "nspec", "site", "spec",
           "offset", "region", "migration")
           
Data <- setNames(Data, nm = heads)

#-------------#
#-Export data-#
#-------------#

save(Data, file = "./DataFormatting/FormattedData.Rdata")

