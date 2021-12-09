
library(tidymodels)
library(ggplot2)
library(purrr)
library(dplyr)
library(spatialsample)
library(reshape2)
library(dismo)

current <- raster::stack(
  list.files(
    path="./Data/Study region/",
    full.names=TRUE,
    pattern='.tif'
  )
)

other <- raster::stack(
  list.files(
    path="./Data/Study Region vegetation and other/",
    full.names=TRUE,
    pattern='.tif'
  )
)

current<-stack(other, current)
rm(other)

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
res(r) <- .3
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

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2",  "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5",
                 "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_15",
                 "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19")


pres_df <- cbind(pres, pres_df[ ,candidates])
back_df <- cbind(back, back_df[ ,candidates])

trainData <- rbind(
  pres_df,
  back_df
)

dim(trainData)

presBg <- c(rep(1, nrow(pres_df)), rep(0, nrow(back_df)))
length(presBg)
presBg <- data.frame(presBg)
names(presBg) <- "pres"

pres_all <- bind_cols( presBg, trainData)
pres_all <- na.omit(pres_all)

# spatial folds -----------------------------------------------------------

# https://spatialsample.tidymodels.org/
folds <- spatial_clustering_cv(df_all, coords = c("Latitude", "Longitude"), v = 10)

folds

plot_splits <- function(split) {
  p <- analysis(split) %>%
    mutate(analysis = "Analysis") %>%
    bind_rows(assessment(split) %>%
                mutate(analysis = "Assessment")) %>%
    ggplot(aes(Longitude, Latitude, color = analysis)) + 
    geom_point(alpha = 0.5) +
    labs(color = NULL)
  print(p)
}

walk(folds$splits, plot_splits)


# tidymodels  -------------------------------------------------------------

rec <-  recipe(pres ~ . , data = df_all) 

mod <- logistic_reg(penalty = tune(), mixture = tune()) %>% 
  set_mode("regression") %>%
  set_engine("glmnet")

wf <- workflow(rec, mod)


grid <- grid_regular(penalty(range = c(-5, 0)), levels = 20)

doParallel::registerDoParallel()
rs <-
  tune_grid(
    wf,
    folds,
    grid = grid
  )

rs

autoplot(rs)


# predict -----------------------------------------------------------------

mean_pred <- predict(rs, new_data = new_points)

#