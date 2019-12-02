####################################################################################################
#
# Elif Kardas - Dec 2019 - ENMeval for Apidae of Puerto Rico Island
# from https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html
# and bioclim var : https://cran.r-project.org/web/packages/dismo/dismo.pdf
#
####################################################################################################
# SPECIES 1: APIS MELLIFERA #
####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval 
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Apis mellifera', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for A. mellifera

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occsamgbif <- as.data.frame(bv$gbif$data$Apis_mellifera[,2:3])
occsambison <- as.data.frame(bv$bison$data$Apis_mellifera[,2:3])
occsam <- rbind(occsambison,occsamgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsam <- occsam[!duplicated(occsam),]

####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsam.sp <- SpatialPoints(occsam)

# Get the bounding box of the points
bb <- bbox(occsam.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.
bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsam)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsam, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsam, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
evalam <- ENMevaluate(occ=occsam, env=envs, bg.coords=bg, method='checkerboard2', 
                               RMvalues=c(1,2), fc=c('L','LQ','LQP'), parallel=FALSE, 
                               algorithm='maxent.jar')
evalam@results
evalam@results[which(evalam@results$delta.AICc==0),]
evalam@predictions

plot(evalam@predictions[[which(evalam@results$delta.AICc==0)]], main="Relative occurrence rate of Apis mellifera", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))

# Let’s see how model complexity changes the predictions in our example. We’ll compare the model predictions of the model with only linear feature classes and with the highest regularization multiplier value we used (i.e., fc=‘L’, RM=2) versus the model with all feature class combination and the lowest regularization multiplier value we used (i.e., fc=‘LQP’, RM=1).
# bisect the plotting area to make two columns
par(mfrow=c(1,2), mar=c(2,2,1,0))

plot(evalam@predictions[['L_1']], xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), legend=T, main='L_1 prediction')

plot(evalam@predictions[['LQP_1']], xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), legend=T, main='LQP_1 prediction')


####################################################################################################
# SPECIES 2: CENTRIS HAEMORRHOIDALIS
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Centris haemorrhoidalis', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris haemorrhoidalis

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occschgbif <- as.data.frame(bv$gbif$data$Centris_haemorrhoidalis[,2:3])
occschbison <- as.data.frame(bv$bison$data$Centris_haemorrhoidalis[,2:3])
occsch <- rbind(occschbison,occschgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsch <- occsch[!duplicated(occsch),]

####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsch.sp <- SpatialPoints(occsch)

# Get the bounding box of the points
bb <- bbox(occsch.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsch)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsch, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalch <- ENMevaluate(occ=occsch, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plot(evalch@predictions[[which(evalch@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Centris haemorrhoidalis", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))

####################################################################################################
# SPECIES 3: ANTHOPHORA KRUGII
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Anthophora krugii', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris decolorata

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occsakgbif <- as.data.frame(bv$gbif$data$Anthophora_krugii[,2:3])
occsakbison <- as.data.frame(bv$bison$data$Anthophora_krugii[,2:3])
occsak <- rbind(occsakbison,occsakgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsak <- occsak[!duplicated(occsak),]
#
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))

####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsak.sp <- SpatialPoints(occsak)

# Get the bounding box of the points
bb <- bbox(occsak.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsak)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsak, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalak <- ENMevaluate(occ=occsak, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plot(eval2@predictions[[which(eval2@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Anthophora krugii", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))

####################################################################################################
# SPECIES 4: CENTRIS DECOLORATA
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Centris decolorata', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris decolorata

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occscdgbif <- as.data.frame(bv$gbif$data$Centris_decolorata[,2:3])
occscdbison <- as.data.frame(bv$bison$data$Centris_decolorata[,2:3])
occscd <- rbind(occscdbison,occscdgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
# occscd <- occscd[!duplicated(occscd),]
# not for this one
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occscd.sp <- SpatialPoints(occscd)

# Get the bounding box of the points
bb <- bbox(occscd.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occscd)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occscd, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalcd <- ENMevaluate(occ=occscd, env=envs, bg.coords=bg, method='checkerboard2', 
                               RMvalues=c(1,2), fc=c('L','LQ','LQP'), parallel=FALSE, 
                               algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plot(evalcd@predictions[[which(evalcd@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Centris decolorata", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))

####################################################################################################
# SPECIES 5: XYLOCOPA MORDAX
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Xylocopa mordax', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris decolorata

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occsxmgbif <- as.data.frame(bv$gbif$data$Xylocopa_mordax[,2:3])
occsxmbison <- as.data.frame(bv$bison$data$Xylocopa_mordax[,2:3])
occsxm <- rbind(occsxmbison,occsxmgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsxm <- occsxm[!duplicated(occsxm),]
#
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsxm.sp <- SpatialPoints(occsxm)

# Get the bounding box of the points
bb <- bbox(occsxm.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsxm)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsxm, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalxm <- ENMevaluate(occ=occsxm, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plotxm <- plot(evalxm@predictions[[which(evalxm@results$delta.AICc==0)]], 
               main="Relative occurrence rate of Xylocopa mordax", 
               xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))


####################################################################################################
# SPECIES 6: EXOMALOPSIS SIMILIS 
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Exomalopsis similis', from=c('gbif', 'bison'), 
          limit=500, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris decolorata

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occses <- as.data.frame(bv$gbif$data$Exomalopsis_similis[,2:3])
#occsesbison <- as.data.frame(bv$bison$data$Exomalopsis_similis[,2:3]) #nothing for bison
#occses <- rbind(occsesbison,occsesgbif) 
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
#occses <- occses[!duplicated(occses),]
#
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))

####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occses.sp <- SpatialPoints(occses)

# Get the bounding box of the points
bb <- bbox(occses.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occses)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')


####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evales <- ENMevaluate(occ=occses, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plot(evales@predictions[[which(evales@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Exomalopsis similis", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.005))

####################################################################################################
# SPECIES 7: MELISSODES TRIFASCIATA 
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Melissodes trifasciata', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris decolorata

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occsmtgbif <- as.data.frame(bv$gbif$data$Melissodes_trifasciata[,2:3])
occsmtbison <- as.data.frame(bv$bison$data$Melissodes_trifasciata[,2:3])
occsmt <- rbind(occsmtbison,occsmtgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsmt <- occsmt[!duplicated(occsmt),]
#
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsmt.sp <- SpatialPoints(occsmt)

# Get the bounding box of the points
bb <- bbox(occsmt.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsmt)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsmt, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalmt <- ENMevaluate(occ=occsmt, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plotmt <- plot(evalmt@predictions[[which(evalmt@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Melissodes trifasciata", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))

####################################################################################################
# SPECIES 8: XEROMELECTA TIBIALIS 
####################################################################################################

####################################################################################################
#
# Elif Kardas - Nov 2019 - Test ENMeval
# from https://cran.r-project.org/web/packages/ENMeval/vignettes/ENMeval-vignette.html
# and bioclim var : https://cran.r-project.org/web/packages/dismo/dismo.pdf
#
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Xeromelecta tibialis', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates
## only for Centris decolorata

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occsxtgbif <- as.data.frame(bv$gbif$data$Xeromelecta_tibialis[,2:3])
occsxtbison <- as.data.frame(bv$bison$data$Xeromelecta_tibialis[,2:3])
occsxt <- rbind(occsxtbison,occsxtgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsxt <- occsxt[!duplicated(occsxt),]
#
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsxt.sp <- SpatialPoints(occsxt)

# Get the bounding box of the points
bb <- bbox(occsxt.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsxt)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')

####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsxt, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalxt <- ENMevaluate(occ=occsxt, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plot(eval2@predictions[[which(eval2@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Xeromelecta tibialis", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00009, 0.005))

####################################################################################################
# SPECIES 9: NOMADA KRUGII 
####################################################################################################

####################################################################################################
# 
#
# A. Acquisition and pre-processing of input data for ENMeval
#
#
####################################################################################################

####################################################################################################
# Occurence dataset for the bees downloading 
####################################################################################################

if (!require('spocc')) install.packages('spocc', repos="http://cran.us.r-project.org")
if (!require('ENMeval')) install.packages('ENMeval', repos="http://cran.us.r-project.org")

library(spocc) # We are using spocc to download occurrence records for species
library(ENMeval)

# Search GBIF for occurrence data+ rest
pr <- c(-67.3, 17.5, -65.3, 18.6) # with a bigger extent, rounded:
# pr <- c(-68, 17, -66, 19)
# apsp <- c('Apis mellifera', 'Centris versicolor', 'Centris haemorrhoidalis','Centris decolorata', 
#           'Anthophora krugii','Xylocopa mordax','Exomalopsis similis','Melissodes trifasciata',
#           'Xeromelecta tibialis', 'Nomada krugii')
bv <- occ(query='Nomada krugii', from=c('gbif', 'bison'), 
          limit=5000, has_coords=TRUE, geometry = pr) 
## has_coords=TRUE to only consider data with coordinates

# Get the latitude/coordinates for each locality. Also convert the tibble that occ() outputs
# to a data frame for compatibility with ENMeval functions.
occsnkgbif <- as.data.frame(bv$gbif$data$Nomada_krugii[,2:3])
occsnkbison <- as.data.frame(bv$bison$data$Nomada_krugii[,2:3])
occsnk <- rbind(occsnkbison,occsnkgbif)
# the [,2:3] are selecting only the coordinates

# Remove duplicate rows (Note that you may or may not want to do this). 
occsnk <- occsnk[!duplicated(occsnk),]
#
####################################################################################################
# Raster predictors for PR from worldclim
####################################################################################################

library(raster)

files <- getData("worldclim",var="bio",res=0.5, lon=-66, lat=18) # higher resolution 

# Put the rasters into a RasterStack:
envs <- stack(files)
envs
#plot(envs) # to see all plots
## it means we we will have all these predators in a big file ("stack")
plot(envs[[1]], xlim=c(-67.3,-65.3), ylim=c(17.5,18.6))
####################################################################################################
# Buffer the species limits
####################################################################################################

library(sp)
# Make a SpatialPoints object
occsnk.sp <- SpatialPoints(occsnk)

# Get the bounding box of the points
bb <- bbox(occsnk.sp)
bb

# Add 5 degrees to each bound by stretching each bound by 10, as the resolution is 0.5 degree.

bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)

envs.backg <- crop(envs, bb.buf)

####################################################################################################
# Crop the environment to these species limits
####################################################################################################

if (!require('maptools')) install.packages('maptools', repos="http://cran.us.r-project.org")
if (!require('rgeos')) install.packages('rgeos', repos="http://cran.us.r-project.org")
library(maptools)
library(rgeos)
# Get a SIMPLE world countries polygon
data(wrld_simpl)

# Get polygons for PR
prr <- wrld_simpl[wrld_simpl@data$NAME=='Puerto Rico',]
plot(prr)
# Both spatial objects have the same geographic coordinate system with slightly 
# different specifications, so just name the coordinate reference system (crs) 
# for ca.sa with that of envs.backg to ensure smooth geoprocessing.
crs(envs.backg) <- crs(prr)
crs(envs.backg)

# Mask envs by this polygon after buffering a bit to make sure not to lose coastline.
prr <- gBuffer(prr)
envs.backg <- mask(envs.backg, prr)
# Let's check our work. We should see only PR.
plot(envs.backg[[1]], main=names(envs.backg)[1], xlab = "Longitude", ylab="Latitude", 
     xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(occsxt)


####################################################################################################
# sample 10000 random points from the background
####################################################################################################

library(dismo)

# Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
bg <- randomPoints(envs.backg[[1]], n=5000)
bg <- as.data.frame(bg)

# Notice how we have pretty good coverage (every cell).
plot(envs.backg[[1]], legend=FALSE, xlim=c(-67.3,-65.2), ylim=c(17.8,18.6))
points(bg, col='red')
####################################################################################################
# Partitioning Occurrences for Evaluation
####################################################################################################
# checkerboard2 method: 
check2 <- get.checkerboard2(occsnk, envs, bg, aggregation.factor=c(5,5))

plot(envs.backg[[1]], col='gray', legend=FALSE, xlim=c(-67.5,-65.5),ylim=c(17.5,19))
points(bg, pch=21, bg=check2$bg.grp)
points(occsch, pch=21, bg=check2$occ.grp, col='white', cex=1.5)

####################################################################################################
# 
#
# B. Running ENMeval
#
#
####################################################################################################
library(rJava)
library(maxnet)
eval2 <- evalnk <- ENMevaluate(occ=occsnk, env=envs, bg.coords=bg, method='checkerboard2', RMvalues=c(1,2),
                               fc=c('L','LQ','LQP'), parallel=FALSE, algorithm='maxent.jar')
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
eval2@predictions
plot(eval2@predictions[[which(eval2@results$delta.AICc==0)]], main="Relative occurrence rate of 
     Xeromelecta tibialis", xlim=c(-67.3,-65.2), ylim=c(17.8,18.6), zlim=c(0.00001, 0.0015))








