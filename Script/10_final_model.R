# Library -----------------------------------------------------------------

library(dismo)
library(tidymodels)

# read data ---------------------------------------------------------------


current <- raster::stack(
  list.files(
    path="./Data/Study region/",
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

rm(list=setdiff(ls(), c("velutina", "current")))

# data --------------------------------------------------------------------

presence <- data.frame(velutina[, c("longitude", "latitude")])
# rasss <- current[[candidates]]
rasss <- current

# generate 10,000 background sites
library(dismo)
# create circles using an arbitrary radius of 50 km =  d=50000
# 50 km to small
x     <- circles(velutina[,3:2], d=1000000, lonlat=TRUE) 
pol   <- rgeos::gUnaryUnion(x@polygons)
# sample randomly from all circles
samp1 <- spsample(pol, 10000, type="random", iter=50)
# get unique cells
cells <- cellFromXY(current$wc2.1_30s_bio_1, samp1)
length(cells)
cells <- unique(cells)
length(cells)
xy    <- xyFromCell(current$wc2.1_30s_bio_1, cells)

# Plot to inspect the results

plot(current$wc2.1_30s_bio_1, axes=TRUE)
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
plot(current$wc2.1_30s_bio_9)
points(back, cex=.5, pch=20, col="red")
points(pres, cex=1, pch=20, col="blue")

# dataframe ---------------------------------------------------------------

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2",  "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5",
                 "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_15",
                 "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19")

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
1      1.0     452      lh          2       51 -3300.323  6715.907    0.000000  1.000000e+00  7.912187e-01
2      0.8     452      lh          2       59 -3291.807  6719.675    3.767548  1.520153e-01  1.202773e-01
3      0.7     452      lh          2       63 -3287.046  6720.876    4.969293  8.335501e-02  6.595204e-02

# folds -------------------------------------------------------------------

fold <- kfold(pres, k=5)
prestest <-  pres[fold == 1, ]
prestrain <- pres[fold != 1, ]

# models ------------------------------------------------------------------

maxent_1 <- dismo::maxent(x=rasss[[candidates]], p=prestrain,  a=back, removeDuplicates=TRUE,
                          
                                        args=c("-J", 
                                               "-P", 
                                               "cloglog", 
                                               "betamultiplier=1", 
                                               "replicatetype=crossvalidate",
                                               "replicates=10",
                                               "randomseed=TRUE",
                                               "threads=3",
                                               
                                               "writeclampgrid=TRUE",
                                               "pictures=TRUE",
                                               "writemess=TRUE",
                                               "extrapolate=TRUE",
                                               
                                               "autofeature=FALSE",
                                               "product=FALSE",   # describe pairwise interactions between predictors 
                                               "hinge=TRUE",      # combine linear and step functions
                                               "threshold=FALSE", # simple step functions
                                               "quadratic=FALSE",  # features use squared predictor values
                                               "linear=TRUE"      # simple linear coefficients for each predictor
                                        ))


show(maxent_1)

preda_ia <- predict(maxent_1, rasss[[candidates]], progress="text")


final_map <- mean(preda_ia$layer.1, preda_ia$layer.2)

# plot the averaged map
plot(final_map)
dev.off()
tiff(filename="./Output/Plots/all_asia_standardMaxent_hql_auc923.tif", width = 30, height = 25, units = 'cm', res = 300)
plot(final_map, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
dev.off()


# Final model  ------------------------------------------------------------

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_4",
                  "wc2.1_30s_bio_5", "wc2.1_30s_bio_11", "wc2.1_30s_bio_17")

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

fold <- kfold(pres, k=5)
prestest <-  pres[fold == 1, ]
prestrain <- pres[fold != 1, ]


maxent_2 <- dismo::maxent(x=rasss[[candidates]], p=prestrain, a=back, removeDuplicates=TRUE,  
                          
                          args=c("-J", 
                                 "-P", 
                                 "cloglog", 
                                 "betamultiplier=1", 
                                 "replicatetype=crossvalidate",
                                 "replicates=5",
                                 "randomseed=TRUE",
                                 "threads=3",
                                 
                                 "writeclampgrid=TRUE",
                                 "pictures=TRUE",
                                 "writemess=TRUE",
                                 "extrapolate=TRUE",
                                 
                                 "autofeature=FALSE",
                                 "product=FALSE",   # describe pairwise interactions between predictors 
                                 "hinge=TRUE",      # combine linear and step functions
                                 "threshold=FALSE", # simple step functions
                                 "quadratic=FALSE",  # features use squared predictor values
                                 "linear=TRUE"      # simple linear coefficients for each predictor
                          ))

show(maxent_2)


# Predict -----------------------------------------------------------------

preda_ia_2 <- predict(maxent_2, rasss[[candidates]], progress="text")

# Evaluate  ---------------------------------------------------------------

mod_eval_train <- dismo::evaluate(p = prestest, a = back, x= rasss[[candidates]], model = maxent_2)

plot(mod_eval_train, 'ROC')
print(mod_eval_train)
# calculate thresholds of models
thd1 <- dismo::threshold(mod_eval_train, "no_omission")  # 0% omission rate 
thd2 <- dismo::threshold(mod_eval_train, "spec_sens")  # highest TSS

# plotting points that are above the previously calculated
# thresholded value
plot(preda_ia_2 >= thd1)

final_map_2 <- mean(preda_ia_2$layer.1, preda_ia_2$layer.2, preda_ia_2$layer.3, preda_ia_2$layer.4, preda_ia_2$layer.5) # Take a mean of all the prediction


writeRaster(final_map_2,
            "./Output/Raster/final_map_2", 
            format='GTiff',
            datatype='INT2S',
            overwrite = TRUE)

# Base plot  --------------------------------------------------------------

# plot the averaged map
tiff(filename="./Maxent/Final_europe_subset_lq_b1_auc085.tif", width = 30, height = 25, units = 'cm', res = 300)
plot(final_map_2, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
dev.off()
