

############################### save data ###################################
dir()
# save.image(file="./Data/all_data_velutina.RData")
# load(file="./Data/all_data_velutina.RData")

##########################################################################


# Exploring bias

tiff(filename="./Output/Plots/Geographic Space.tif", width = 45, height = 15, units = 'cm', res = 300)

par(mfrow=c(1, 2), pty='s')

# plot in geographic space
plot(current$wc2.1_30s_bio_1, main='Geographic Space')
plot(wrld_simpl, add=TRUE, border='gray45')
points(records$longitude,
  records$latitude,
  pch=21,
  bg=alpha('red', 0.5)
)

# plot in environmental space
plot(randomBg[, c("wc2.1_30s_bio_1")],
  randomBg[, c("wc2.1_30s_bio_12")],
  pch=16,
  col=alpha('red', 0.5),
  xlab='MAT (deg)',
  ylab='MAP (mm)',
  main='Environmental Space'
)

points(records$wc2.1_30s_bio_1,
  records$wc2.1_30s_bio_12,
  col=alpha('black'),
  bg='black',
  pch=20
)

legend('topleft',
  legend=c('Background', 'Presences'),
  col=c('red', 'black'),
  pch=16,
  cex=0.6
)

dev.off()

#################### Backround points ###########################################

coordinates(velutina) <- ~longitude+latitude
projection(velutina) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

library(dismo)
# create circles using an arbitrary radius of 50 km =  d=50000
# 50 km to small
x     <- circles(velutina, d=1000000, lonlat=TRUE) 
pol   <- rgeos::gUnaryUnion(x@polygons)
# sample randomly from all circles
samp1 <- spsample(pol, 10000, type="random", iter=50)
# get unique cells
cells <- cellFromXY(climate$wc2.1_30s_bio_1, samp1)
length(cells)
cells <- unique(cells)
length(cells)
xy    <- xyFromCell(climate$wc2.1_30s_bio_1, cells)

# Plot to inspect the results

plot(climate$wc2.1_30s_bio_1, axes=TRUE)
plot(pol, add=T)
points(xy, cex=0.1, pch=20, col="blue")
points(velutina, cex=1, pch=20, col="red")

targetSites <-xy

###########################################################################
dev.off()
# plot background sites on map

tiff(filename="./Output/Plots/Target_BG.tif", width = 25, height = 25, units = 'cm', res = 300)

par(pty='s')

plot(env$wc2.1_30s_bio_1, axes=TRUE)
points(randomBg[,1], randomBg[,2],  pch=20, col="gray50")
plot(pol, add=T)
points(targetSites, cex=0.1, pch=20, col="blue")
points(velutina, cex=1, pch=20, col="red")

legend('topright',
       legend=c('Random', 'Target', 'presence'),
       col=c('gray', 'blue', 'red'),
       pch=c(16, 16, 16),
       pt.cex=c(1, 1, 1),
       bg='white')


dev.off()

###############################################################################

# extract environmental data at target background sites
targetEnv <- as.data.frame(raster::extract(climate, targetSites))

# remove target background sites with no environmental data (fall in ocean)
outside <- which(is.na(rowSums(targetEnv))) # if a record has missing data (NA) then the sum of its rows is also NA... find these (if any)
if (length(outside) > 0) {
  targetSites <- targetSites[-outside, ]
  targetEnv <- targetEnv[-outside, ]
}

# combine coordinates and environmental data
targetBg <- cbind(targetSites, targetEnv)
names(targetBg)[1:2] <- c('longitude', 'latitude')

# how many target records do we have?
# if we had >10000 we could use a random sample of 10000 for expediency
nrow(targetBg)

getwd()
# setwd( "C:/Users/floed/Dropbox/vkm/Vespa velutina/Vespa velutina/")
# save target background sites
save(targetBg, file='./Data/Background Sites/Target Background Sites Drawn from All data Downloaded from GBIF.Rdata', compress=TRUE)

#########################################

dev.off()

tiff(filename="./Output/Plots/Target_BG_scatter.tif", width = 25, height = 25, units = 'cm', res = 300)

# plot in environmental space
plot(randomBg[, c("wc2.1_30s_bio_1")],
     randomBg[, c("wc2.1_30s_bio_12")],
     pch=16,
     col='green')

points(targetBg[, c("wc2.1_30s_bio_1")],
     targetBg[, c("wc2.1_30s_bio_12")],
     pch=16,
     col='black')

points(records$wc2.1_30s_bio_1,
       records$wc2.1_30s_bio_12,
       col='red',
       pch=16)

legend('topright',
       legend=c('Random', 'Target', 'presence'),
       col=c('green', 'black', 'red'),
       pch=c(16, 16, 16),
       pt.cex=c(1, 1, 1),
       bg='white')

dev.off()


# subsample  --------------------------------------------------------------
test <- records[,3:2]
coordinates(test) <- ~latitude +longitude
# create a RasterLayer with the extent of record
r <- raster(test)
# set the resolution of the cells to (for example) 1 degree
res(r) <- 0.2
# expand (extend) the extent of the RasterLayer a little
r <- extend(r, extent(r)+1)
# sample:
sample <- gridSample(test, r, n=1)
# to illustrate the method and show the result
p <- rasterToPolygons(r)
plot(p, border='gray')
points(test)
# selected points in red
points(sample, cex=1, col='red', pch='x')


dim(sample)
dim(records)


# Model with target background  -------------------------------------------


# make training data frame with predictors and vector of 1/0 for presence/background
trainData <- rbind(records[ , predictors], targetBg[ , predictors]) 
presBg <- c(rep(1, nrow(records)), rep(0, nrow(na.omit(targetBg))))

# create output directory for model object and rasters
dir.create('./Models/Model 03 Bias Correction - Target Background',
           recursive=TRUE, showWarnings=FALSE)

# model species
targetBgModel <- maxnet::maxnet(p=presBg, data=na.omit(trainData))

# save model
save(targetBgModel,
     file='./Models/Model 03 Bias Correction - Target Background/Model.Rdata',
     compress=TRUE)

# write raster
targetBgMap <- predict(climate[[predictors]], targetBgModel,
  filename='./Models/Model 03 Bias Correction - Target Background/maxentPrediction1970to2000', 
  format='GTiff',  overwrite=TRUE,  type='cloglog', progress='text')

# plot
dev.off()

tiff(filename="./Output/Plots/Target_VS_random_BG_EC.tif", width = 25, height = 25, units = 'cm', res = 300)

par(mfrow=c(1, 2), pty='s')
plot(manualSelectMap, main='Random Background & Ecologist-Chosen Predictors')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude, pch=21, bg='mediumseagreen')

plot(targetBgMap, main='Targer Background & Ecologist-Chosen Predictors')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude, pch=21, bg='mediumseagreen')

dev.off()
