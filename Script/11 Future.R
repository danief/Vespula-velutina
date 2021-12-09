# lib ---------------------------------------------------------------------

library(dismo)
library(tidymodels)

# read data ---------------------------------------------------------------

future <- raster::stack(
  list.files(
    path="./Data/Study region future/",
    full.names=TRUE,
    pattern='.tif'
  )
)

load(file="./Data/velutina_presence.RData")
velutina <- Europe

# clean elevation ---------------------------------------------------------

elevationraster <- raster("./Data/worldclim/vegetation and other/elevation_1KMmd_GMTEDmd.tif")

velutina$elevation <- raster::extract(elevationraster, cbind(velutina$longitude, velutina$latitude))

velutina <- filter(velutina, elevation < 1500)

rm(list=setdiff(ls(), c("velutina", "future")))

# data --------------------------------------------------------------------

presence <- data.frame(velutina[, c("longitude", "latitude")])
# rasss <- current[[candidates]]
rasss <- future

# generate 10,000 background sites
library(dismo)
# create circles using an arbitrary radius of 50 km =  d=50000
# 50 km to small
x     <- circles(velutina[,3:2], d=1000000, lonlat=TRUE) 
pol   <- rgeos::gUnaryUnion(x@polygons)
# sample randomly from all circles
samp1 <- spsample(pol, 10000, type="random", iter=50)
# get unique cells
cells <- cellFromXY(rasss$mp85bi701, samp1)
length(cells)
cells <- unique(cells)
length(cells)
xy    <- xyFromCell(rasss$mp85bi701, cells)

# Plot to inspect the results

plot(rasss$mp85bi701, axes=TRUE)
plot(pol, add=T)
points(xy, cex=0.1, pch=20, col="blue")
points(velutina[,3:2], cex=1, pch=20, col="red")

targetSites <-xy

# Background --------------------------------------------------------------

back <- raster::extract(rasss, cbind(targetSites[,1], targetSites[,2]))
back <- as.data.frame(back)
back <- cbind(targetSites, back)
head(back)
back <- na.omit(back)
back_df <- back
back <- back[,1:2]
names(back) <- c("longitude", "latitude")

# thinn -------------------------------------------------------------------

test <- presence
#names(test) <- c("x", "y")
coordinates(test) <- ~longitude+latitude
# create a RasterLayer with the extent of acgeo
r <- raster(test)
# set the resolution of the cells to (for example) 1 degree
res(r) <- .4
# expand (extend) the extent of the RasterLayer a little
r <- extend(r, extent(r)+ 1)
# sample:
pres <- gridSample(test, r, n=1)
# to illustrate the method and show the result
p <- rasterToPolygons(r)
plot(p, border='gray')
points(test)
# selected points in red
points(pres, cex=1, col='red', pch='x')
pres

# presence 
env <- raster::extract(rasss, cbind(pres[,1], pres[,2]))
env <- as.data.frame(env)
pres <- cbind(pres, env)
head(pres)
pres <- na.omit(pres)
pres_df <- pres
pres <- pres[,1:2]
names(pres) <- c("longitude", "latitude")

# plot
plot(rasss$mp85bi701)
points(back, cex=.5, pch=20, col="red")
points(pres, cex=1, pch=20, col="blue")
#points(subset(pres, longitude > 60), cex=1, pch=20, col="black")
#pres <- subset(pres, longitude > 60)

# dataframe ---------------------------------------------------------------

candidates <- c(
  "mp85bi701",  "mp85bi7010", "mp85bi7011", "mp85bi7012",
  "mp85bi7013", "mp85bi7014", "mp85bi7015", "mp85bi7016",
  "mp85bi7017", "mp85bi7018", "mp85bi7019", "mp85bi702", 
  "mp85bi703",  "mp85bi704",  "mp85bi705",  "mp85bi706", 
  "mp85bi707",  "mp85bi708",  "mp85bi709" )


candidates <- c("mp85bi701", "mp85bi703", "mp85bi704", "mp85bi709", "mp85bi7012")

trainData <- rbind(
  pres_df[ ,candidates],
  back_df[ ,candidates]
)

dim(trainData)

presBg <- c(rep(1, nrow(pres_df)), rep(0, nrow(back_df)))
length(presBg)
presBg <- data.frame(presBg)
names(presBg) <- "pres"

df_all <- bind_cols(presBg, trainData)
df_all <- na.omit(df_all)

# Test combinations -------------------------------------------------------

enmSdm::trainMaxNet(data=df_all, regMult=seq(.1, 1, by=.1), classes='lqht', verbose=T, testClasses = TRUE)

regMult numPres classes numClasses numCoeff    logLik      AICc   deltaAICc       relLike     aicWeight
1      0.9     450      lh          2       44 -3258.188  6614.153    0.000000  1.000000e+00  6.701435e-01
2      1.0     450      lh          2       42 -3261.478  6615.830    1.676624  4.324400e-01  2.897968e-01
3      1.0     450     lqh          3       44 -3261.629  6621.035    6.881580  3.203937e-02  2.147097e-02

regMult numPres classes numClasses numCoeff    logLik     AICc  deltaAICc      relLike    aicWeight
1      0.6     450     lqh          3       26 -3303.474 6662.268   0.000000 1.000000e+00 3.929066e-01
2      0.7     450      lh          2       26 -3304.104 6663.526   1.258512 5.329883e-01 2.094146e-01
3      0.6     450      lh          2       27 -3303.509 6664.600   2.332322 3.115607e-01 1.224143e-01

# folds -------------------------------------------------------------------

fold <- kfold(pres, k=5)
prestest <-  pres[fold == 1, ]
prestrain <- pres[fold != 1, ]

# Model -------------------------------------------------------------------

maxent_future <- dismo::maxent(x=rasss[[candidates]], p=prestrain, a=back, removeDuplicates=TRUE,  
                          
                          args=c("-J", 
                                 "-P", 
                                 "cloglog", 
                                 "betamultiplier=0.8", 
                                 "replicatetype=crossvalidate",
                                 "replicates=5",
                                 "randomseed=TRUE",
                                 "threads=3",
                                 
                                 "writeclampgrid=TRUE",
                                 "pictures=TRUE",
                                 "writemess=TRUE",
                                 "extrapolate=TRUE",
                                 
                                 "autofeature=FALSE",
                                 "product=FALSE",    # describe pairwise interactions between predictors 
                                 "hinge=TRUE",       # combine linear and step functions
                                 "threshold=FALSE",   # simple step functions
                                 "quadratic=FALSE",  # features use squared predictor values
                                 "linear=TRUE"       # simple linear coefficients for each predictor
                          ))

show(maxent_future)

# Predict -----------------------------------------------------------------

preda_future<- predict(maxent_future, rasss[[candidates]], progress="text")

# Plot --------------------------------------------------------------------

future_map <- mean(preda_future$layer.1, preda_future$layer.2, preda_future$layer.3, preda_future$layer.4, preda_future$layer.5) 


writeRaster(future_map,
            "./Output/Raster/future_map", 
            format='GTiff',
            datatype='INT2S',
            overwrite = TRUE)

# eval --------------------------------------------------------------------

mod_eval_train <- dismo::evaluate(p = prestest, a = back, x= rasss[[candidates]], model = maxent_future)

plot(mod_eval_train, 'ROC')
print(mod_eval_train)
# calculate thresholds of models
thd1 <- dismo::threshold(mod_eval_train, "no_omission")  # 0% omission rate 
thd2 <- dismo::threshold(mod_eval_train, "spec_sens")  # highest TSS

# plotting points that are above the previously calculated
# thresholded value

tiff(filename="./Maxent/threshold_future_europe_subset_lh_b.7_auc089_10fold.tif", width = 30, height = 25, units = 'cm', res = 300)
plot(future_map >= thd1)
dev.off()


# Base plot ---------------------------------------------------------------


tiff(filename="./Maxent/future_europe_subset_lh_b.7_auc089_10fold.tif", width = 30, height = 25, units = 'cm', res = 300)
plot(future_map, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
dev.off()


# End ---------------------------------------------------------------------


