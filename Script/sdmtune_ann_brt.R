
library(SDMtune)
library(zeallot)

# https://consbiol-unibern.github.io/SDMtune/articles/articles/evaluation_strategies.html


current <- raster::stack(
  list.files(
    path="./Data/Study region/",
    full.names=TRUE,
    pattern='.tif'
  )
)

load(file="./Data/velutina_presence.RData")
velutina <- Europe
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
#points(subset(pres, longitude > 60), cex=1, pch=20, col="black")
#pres <- subset(pres, longitude > 60)


# dataframe ---------------------------------------------------------------

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

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_4",
                 "wc2.1_30s_bio_5", "wc2.1_30s_bio_11", "wc2.1_30s_bio_17")
# prep data  --------------------------------------------------------------

data <- prepareSWD(species = "velutina", p = pres, a = back, env = rasss)

#c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE, seed = 25)
c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE, seed = 12345)

# train models  -----------------------------------------------------------

output <- SDMtune::train(method = c("ANN", "BRT", "RF"), data = data, size = 10)

# ANN ---------------------------------------------------------------------

selected_variables_model <- varSel(output$ANN, metric = "auc", test = test,
                                   bg4cor = test, method = "spearman",
                                   cor_th = 0.7, permut = 10)

getTunableArgs(output$ANN)

# size: integer. Number of the units in the hidden layer.
# decay numeric. Weight decay, default is 0.
# rang numeric. Initial random weights, default is 0.7.
# maxit integer. Maximum number of iterations, default is 100.

h <- list(size = 40:100, decay = seq(0.01, 1, by=0.1), rang=seq(0.01, 1, by=0.1) ,  maxit = seq(5000, 10000, by=1000))

om <- optimizeModel(selected_variables_model, hypers = h, metric = "auc", seed = 25, test=test)

best_model <- om@models[[1]]
om@results[1, ]

#folds <- randomFolds(train, k = 5)

final_model <- train("ANN", data = train, size = 47, decay = om@results[1, 2], rang= om@results[1, 3], maxit = om@results[1, 4])
plotROC(final_model, test = test)

map <- predict(final_model, data = rasss,  file = "ann_map_aurope_selctvar_auc976", format = "GTiff", overwrite=TRUE)

tiff(filename="./Output/Plots/ann_map_aurope_selctvar_auc976.tif", width = 25, height = 25, units = 'cm', res = 300)
plot(map,  col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
title(main = "Artificial Neural Networks", 
      sub = "Presence data for Europe and selected variables")

dev.off()


# BRT ---------------------------------------------------------------------

getTunableArgs(output$BRT)
# distribution: character. Name of the distribution to use, default is "bernoulli".
# n.trees: integer. Maximum number of tree to grow, default is 100.
# interaction.depth: integer. Maximum depth of each tree, default is 1.
# shrinkage: numeric. The shrinkage parameter, default is 0.1.
# bag.fraction: numeric. Random fraction of data used in the tree expansion, default is 0.5.

h <- list(distribution="bernoulli", n.trees=c(1500, 2000, 3000, 5000), interaction.depth=c(5, 7, 9), shrinkage=c(0.01, 0.05 , 0.1, 0.2), bag.fraction=0.5)

om <- optimizeModel(output$BRT, metric = "auc",  test=test, hypers = h)

om@results[1, ]

#folds <- randomFolds(train, k = 5)

final_model <- train("BRT", data = train, size = om@results[1, 1], decay = om@results[1, 2], maxit = om@results[1, 4])
plotROC(final_model, test = test)

map <- predict(final_model, data = rasss,  file = "brt_map_aurope_19_auc895", format = "GTiff", overwrite=TRUE)

tiff(filename="./Output/Plots/BRT_19bioclim_europe_auc895.tif", width = 25, height = 25, units = 'cm', res = 300)
plot(map,  col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
title(main = "Boosted regression trees - 5 fold cross validation", 
      sub = "Presence data for Europe and all 19 bioclim variables")

dev.off():dev.off()

writeRaster(map,
            "brt_pred", 
            format='GTiff',
            datatype='INT2S',
            overwrite = TRUE)


# RF ----------------------------------------------------------------------

getTunableArgs(output$RF)
# mtry: integer. Number of variable randomly sampled at each split, default is floor(sqrt(number of variables)).
# ntree: integer. Number of tree to grow, default is 500.
# nodesize: integer. Minimum size of terminal nodes, default is 1.

h <- list(ntree=seq(from = 500, to = 1000, by=500), mtry=c(3:10), nodesize=seq(from = 1, to = 100, by=1))

om <- optimizeModel(output$RF, metric = "auc", test=test, hypers = h)

om@results[1, ]

#folds <- randomFolds(train, k = 5)

final_model <- train("BRT", data = train, size = om@results[1, 1], decay = om@results[1, 2], maxit = om@results[1, 4])
plotROC(final_model, test = test)

map <- predict(final_model, data = rasss,  file = "brt_map_aurope_19_auc895", format = "GTiff", overwrite=TRUE)

tiff(filename="./Output/Plots/BRT_19bioclim_europe_auc895.tif", width = 25, height = 25, units = 'cm', res = 300)
plot(map,  col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
title(main = "Boosted regression trees - 5 fold cross validation", 
      sub = "Presence data for Europe and all 19 bioclim variables")

dev.off():dev.off()
