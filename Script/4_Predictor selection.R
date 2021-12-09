
library(dismo)

# Predictor selection

memory.limit()
memory.limit(size=90000)

current
climate <- env

########################### RANDOM BG #####################################

# generate 10,000 background sites
randomBgSites <- randomPoints(climate, 10000)

# extract environment at sites
randomBgEnv <- raster::extract(climate, randomBgSites)
randomBgEnv <- as.data.frame(randomBgEnv)

# remove any sites with NA for at least one variable
isNa <- is.na(rowSums(randomBgEnv))
if (any(isNa)) {
  randomBgSites <- randomBgSites[-which(isNa), ]
  randomBgEnv <- randomBgEnv[-which(isNa), ]
}

# combine with coordinates and rename coordinate fields
randomBg <- cbind(randomBgSites, randomBgEnv)
names(randomBg)[1:2] <- c('longitude', 'latitude')
head(randomBg)

dev.off()
plot(climate[['wc2.1_30s_bio_1']], main='Mean temp') # plot first raster in climate stack
points(randomBgSites, pch='.')

dir.create('./Data/Background Sites/', recursive=TRUE, showWarnings=FALSE)
save(randomBg, file='./Data/Background Sites/Random Background.Rdata',  compress=TRUE)

################## model random bg

dim(records)
dim(randomBgEnv)

trainDataAuto <- rbind(records[,4:32], randomBgEnv)

presBg <- c(rep(1, nrow(records)), rep(0, nrow(randomBgEnv)))

library(maxnet)
autoSelectModel <- maxnet(p=presBg, data=trainDataAuto)

getwd()
dir.create('./Models/Model 01 Predictors - Automated Selection', recursive=TRUE,
           showWarnings=FALSE)

save(autoSelectModel,
     file='./Models/Model 01 Predictors - Automated Selection/Model.Rdata',
     compress=TRUE)

autoSelectMap <- predict(climate, autoSelectModel,
                         filename='./Models/Model 01 Predictors - Automated Selection/maxentPrediction1970to2000',
                         format='GTiff', overwrite=TRUE, type='cloglog',  progress='text')


# plot output
#x11()
dev.off()
tiff(filename="./Output/Plots/Random_backround.tif", width = 25, height = 24, units = 'cm', res = 300)
plot(autoSelectMap, main="Random backround")
plot(wrld_simpl, add=TRUE)
points(records$longitude, records$latitude, pch=18, col="red")
dev.off()
###########################################################

dev.off()

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5" ,
 "wc2.1_30s_bio_6" , "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" ,  "wc2.1_30s_bio_10",
 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", 
 "Barren", "Cultivated_and_Managed_Vegetation", "Deciduous_Broadleaf_Trees", "Deciduous_Needleleaf_Trees", 
  "Evergreen_Broadleaf_Trees", "Herbaceous_Vegetation", "Regularly_Flooded_Vegetation")

#candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2",  "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5" ,
#                 "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
#                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_15",
#                 "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19") 

correl <- cor(records[ , candidates],  method='spearman')
print(correl, digits=2)

pos <- correl > 0.7
neg <- correl < -0.7

tiff(filename="./Output/Plots/cor_all.tif", width = 25, height = 24, units = 'cm', res = 300)

legendary::spoke(
  pos=pos,
  neg=neg,
  lwdPos=2,
  lwdNeg=2,
  colPos='black',
  colNeg='red',
  pty='s'
)

dev.off()
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))

# BIO3 = Isothermality (BIO2/BIO7) (?100)
# BIO4 = Temperature Seasonality (standard deviation ?100

# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)

# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter

# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)

# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

dev.off()



# no response = 2, 3, 4, 5, 6, 7, 10, 11, 13, 14, 16, 17, 18, 19

#candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5", "wc2.1_30s_bio_8",
#                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14" )




candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2",  "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5" ,
                 "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_15",
                 "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19")


correl <- cor(records[ , candidates],  method='spearman')
print(correl, digits=2)

pos <- correl > 0.7
neg <- correl < -0.7

tiff(filename="./Output/Plots/cor_bio.tif", width = 25, height = 24, units = 'cm', res = 300)

legendary::spoke(
  pos=pos,
  neg=neg,
  lwdPos=2,
  lwdNeg=2,
  colPos='black',
  colNeg='red',
  pty='s'
)

dev.off()



# subset  -----------------------------------------------------------------

#Throw out det one that are hard to enterperate 
# BIO3 = Isothermality (BIO2/BIO7) (?100)
# BIO4 = Temperature Seasonality (standard deviation ?100
# BIO15 = Precipitation Seasonality (Coefficient of Variation)

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_5" ,
                 "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", 
                 "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19")

correl <- cor(records[ , candidates],  method='spearman')
print(correl, digits=2)

pos <- correl > 0.7
neg <- correl < -0.7

tiff(filename="./Output/Plots/cor_bio_sub.tif", width = 25, height = 24, units = 'cm', res = 300)

legendary::spoke(
  pos=pos,
  neg=neg,
  lwdPos=2,
  lwdNeg=2,
  colPos='black',
  colNeg='red',
  pty='s'
)

dev.off()


# cor ---------------------------------------------------------------------



# ecospat package 
library(ecospat)

records %>% dplyr::select("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_5" ,
                     "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                     "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", 
                     "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19") %>% ecospat.cor.plot()


library(usdm)
records %>% dplyr::select("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_5" ,
                          "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                          "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", 
                          "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19") %>% vif()

# We can see that almost all the variables are above a value of 10.0. 
# Usually values from 5 to 10 are considered as critical for multi-variable correla- tion


records %>% dplyr::select("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_5" ,
                          "wc2.1_30s_bio_6" , "wc2.1_30s_bio_7",  "wc2.1_30s_bio_8",  "wc2.1_30s_bio_9" , "wc2.1_30s_bio_10",
                          "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", 
                          "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19") %>% vifcor( th=.7) 


records %>% dplyr::select("wc2.1_30s_bio_1", "wc2.1_30s_bio_2", "wc2.1_30s_bio_7",  "wc2.1_30s_bio_13",
                          "wc2.1_30s_bio_17") %>% ecospat.cor.plot()



candidates <- c("wc2.1_30s_bio_2", "wc2.1_30s_bio_5" , "wc2.1_30s_bio_8",  
                "wc2.1_30s_bio_11", "wc2.1_30s_bio_19",
                 
                 "Barren", "Cultivated_and_Managed_Vegetation", 
                 "Deciduous_Broadleaf_Trees", "Deciduous_Needleleaf_Trees", 
                 "Evergreen_Broadleaf_Trees", "Herbaceous_Vegetation", 
                 "Regularly_Flooded_Vegetation")

correl <- cor(records[ , candidates],  method='spearman')
print(correl, digits=2)

pos <- correl > 0.7
neg <- correl < -0.7

tiff(filename="./Output/Plots/cor_bio_sub.tif", width = 25, height = 24, units = 'cm', res = 300)

legendary::spoke(
  pos=pos,
  neg=neg,
  lwdPos=2,
  lwdNeg=2,
  colPos='black',
  colNeg='red',
  pty='s'
)

dev.off()

############################


predictors <-  c("wc2.1_30s_bio_2", "wc2.1_30s_bio_5" , "wc2.1_30s_bio_8",  
                 "wc2.1_30s_bio_11", "wc2.1_30s_bio_19",
                 
                 "Barren", "Cultivated_and_Managed_Vegetation", 
                 "Deciduous_Broadleaf_Trees", "Deciduous_Needleleaf_Trees", 
                 "Evergreen_Broadleaf_Trees", "Herbaceous_Vegetation", 
                 "Regularly_Flooded_Vegetation")




# select chosen predictor variables

trainDataManual <- rbind(records[ , predictors], randomBg[ , predictors])

dim(trainDataManual)
  
presBg <- c(rep(1, nrow(records)), rep(0, nrow(randomBg)))

length(presBg)

# model
manualSelectModel <- maxnet(p=presBg, data=trainDataManual)

# save model
save(manualSelectModel,
     file="./Models/Manual_Selection_Model.Rdata",
     compress=TRUE)

# subset new raster containing only variables that are uncorrelated and gives response 
climateSelect <- subset(climate,  c("wc2.1_30s_bio_2", "wc2.1_30s_bio_5" , "wc2.1_30s_bio_8",  
                                    "wc2.1_30s_bio_11", "wc2.1_30s_bio_19",
                                    
                                    "Barren", "Cultivated_and_Managed_Vegetation", 
                                    "Deciduous_Broadleaf_Trees", "Deciduous_Needleleaf_Trees", 
                                    "Evergreen_Broadleaf_Trees", "Herbaceous_Vegetation", 
                                    "Regularly_Flooded_Vegetation"))
getwd()
# speeds things up to use just the necessary rasters
manualSelectMap <- predict(climateSelect, manualSelectModel,
                           filename='./maxentPrediction1970to2000',
                           format='GTiff', overwrite=TRUE, type='cloglog', progress='text')

dev.off():dev.off()
# x11()


tiff(filename="./Output/Plots/Model__and_Ecologist__random_BG.tif", width = 45, height = 15, units = 'cm', res = 300)

par(mfrow=c(1, 2))
plot(autoSelectMap, main='Model-Selected Predictors with random BG') 
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude, pch=21, bg=alpha('mediumseagreen', 0.5), cex=1.4)

plot(manualSelectMap, main='Ecologist-Chosen Predictors with random BG')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude, pch=21, bg=alpha('mediumseagreen', 0.5), cex=1.4)
dev.off()

