
# http://www.earthenv.org/landcover 

####################### Load library ######################################

library(raster)
library(maptools)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)


# Set exstent  ------------------------------------------------------------

data(wrld_simpl)

eurasia <- ne_countries(continent = c('europe', "asia"), scale = "medium", returnclass = "sp" )
extent(bbox(eurasia))

extent <- c(-25, 180, -21.36904, 81.8542)

max_map_df <- fortify(eurasia) %>% mutate(long = ifelse(long < -100, long + 180*2, long))

ggplot() +
  geom_map(data=max_map_df, map=max_map_df, aes(map_id=id, x=long, y=lat), fill="black") +
 # geom_point(data=data, aes(x=longitude, y=latitude , col=name)) +
  scale_x_continuous(limits = c(-25, 180), expand = c(0, 0))

############################  current worldclim  variables ######################################
getwd()

dir.create("./Data/Study Region/", recursive=TRUE, showWarnings=FALSE)

referance <- raster("./Data/worldclim/bioclim/wc2.1_30s_bio_1.tif")
#referance <- raster("C:/Users/floed/Dropbox/Min PC (MSI)/Desktop/worldclim/present/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")

library(raster)
library(doParallel)  # Foreach Parallel Adaptor 
library(foreach)     # Provides foreach looping construct

#Define how many cores you want to use
UseCores <- detectCores() -1
#Register CoreCluster
cl       <- makeCluster(UseCores)
registerDoParallel(cl)

# list of variavle names for loop
stack_list <- list.files(path = "./Data/worldclim/bioclim", full.names = TRUE, pattern = "\\.tif$", recursive=F)

# for each BIOCLIM raster
#Use foreach loop and %dopar% command
foreach(i=1:length(stack_list), .packages = c("raster")) %dopar% {
  
  files <- list.files(path = "./Data/worldclim/bioclim", full.names = TRUE, pattern = "\\.tif$", recursive=F)
  
  rast <- raster(files[i])
  
  filename <- (paste(basename(files[i]), sep = ""))
  
  # crop the raster and update tt, 
  rast <- setMinMax(rast)
  projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  rast <-raster::resample(rast, referance, method="bilinear")
  # Crop the raster
  rast <- crop(rast, extent, snap = "out")
  
  rast<- mask(rast, eurasia)
  
  # save raster
  writeRaster(rast,
              paste0("./Data/Study Region/", filename, i), 
              format='GTiff',
              datatype='INT2S',
              overwrite = TRUE)
  
  
}

#end cluster
stopCluster(cl)

# visually check
current <- raster::stack(
  list.files(
    path="./Data/Study region/",
    full.names=TRUE,
    pattern='.tif'
  )
)


dir("./Data/")
other <- raster::stack(
  list.files(
    path="./Data/Study Region vegetation and other/",
    full.names=TRUE,
    pattern='.tif'
  )
)

current<-stack(other, current)
rm(other)

################################### CMIP5 8.5 50 years #####################################

getwd()

dir.create("./Data/Study Region future/", recursive=TRUE, showWarnings=T)

referance <- raster("./Data/worldclim/bioclim/wc2.1_30s_bio_1.tif")
#referance <- raster("C:/Users/floed/Dropbox/Min PC (MSI)/Desktop/worldclim/present/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")

library(raster)
library(doParallel)  # Foreach Parallel Adaptor 
library(foreach)     # Provides foreach looping construct

#Define how many cores you want to use
UseCores <- detectCores() -1
#Register CoreCluster
cl       <- makeCluster(UseCores)
registerDoParallel(cl)

# list of variavle names for loop
stack_list <- list.files(path = "./Data/worldclim/bioclim_MPI-ESM-LR_rcp85_2070_30s", full.names = TRUE, pattern = "\\.tif$", recursive=F)

# for each BIOCLIM raster
#Use foreach loop and %dopar% command
foreach(i=1:length(stack_list), .packages = c("raster")) %dopar% {
  
  files <- list.files(path = "./Data/worldclim/bioclim_MPI-ESM-LR_rcp85_2070_30s", full.names = TRUE, pattern = "\\.tif$", recursive=F)
  
  rast <- raster(files[i])
  
  filename <- (paste(basename(files[i]), sep = ""))
  
  # crop the raster and update tt, 
  rast <- setMinMax(rast)
  projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  rast <-raster::resample(rast, referance, method="bilinear")
  # Crop the raster
  rast <- crop(rast, extent, snap = "out")
  
  rast<- mask(rast, eurasia)
  
  # save raster
  writeRaster(rast,
              paste0("./Data/Study Region future/", filename, i), 
              format='GTiff',
              datatype='INT2S',
              overwrite = TRUE)
  
}

#end cluster
stopCluster(cl)

# visually check
future <- raster::stack(
  list.files(
    path="./Data/Study Region future/",
    full.names=TRUE,
    pattern='.tif'
  )
)

plot(future)

######################### difference #########################


par(mfrow=c(1, 3))

currentTemp <- raster('C:/Users/floed/Dropbox/Min PC (MSI)/Desktop/worldclim/present/wc2.1_2.5m_bio/Study Region/01.tif')

futureTemp <- raster('C:/Users/floed/Dropbox/Min PC (MSI)/Desktop/worldclim/Future/cmip5/2_5m/Study Region/01.tif')

deltaTemp <- futureTemp - currentTemp

plot(deltaTemp, main='Change in Temperature')
plot(currentTemp, main='currentTempe')
plot(futureTemp, main='futureTemp')

#################################################################

setwd("C:/Users/floed/Dropbox/Karl")
dir()
# save.image("bioclim.RData")
# load("bioclim.RData")
