
#########################  Model Tuning #########################  

# get min/max value of each predictor across study region
minPred <- minValue(climate[[predictors]])
maxPred <- maxValue(climate[[predictors]])
names(minPred) <- names(maxPred) <- names(climate[[predictors]])

# get median value of each predictor across species' thinned presences
medianPred <- apply(records[ , predictors], 2, median)

# make data frame with median value of each predictor
env <- as.data.frame(medianPred)
env <- t(env)
env <- env[rep(1, 100), ]
row.names(env) <- 1:nrow(env)
head(env)

par(mfrow=c(5, 5), mai=c(0.3,0.3,0.3,0.3))

# for each predictor...
for (pred in predictors) {
  
  # make copy of data frame
  thisEnv <- env
  
  # now vary focal predictor from min to max value... 
  # all other predictors keep median value
  thisEnv[ , pred] <- seq(minPred[pred], maxPred[pred], length.out=100)
  
  # make prediction using this data frame
  prediction <- predict(targetBgModel, thisEnv, type='cloglog')
  
  # plot
  plot(x=thisEnv[ , pred], y=prediction, ylim=c(0, 1), xlab=pred,
       ylab='Suitability', main=pred, type='l', col='blue',
       lty='dotted', lwd=2)
  
  # add species' presences (top rug)
  rug(records[ , pred], side=3, col='mediumseagreen')
  
  # add background sites (bottom rug)
  rug(randomBg[ , pred], side=1, col='red')
  
}

#########################  manual Feature selectio n#########################  

# select only linear and quadratic features ! 

trainData <- rbind(
  records[ , predictors],
  targetBg[ , predictors]
)

presBg <- c(rep(1, nrow(records)), rep(0, nrow(na.omit(targetBg))))

# create output directory for model object and rasters
dir.create('./Models/Model 07 Model Tuning - Feature Selection', recursive=TRUE, showWarnings=FALSE)

#  model using only lqp features
f <- maxnet.formula(p=as.vector(presBg), data=trainData, classes="lqh")
lpqModel <- maxnet(p=presBg, data=na.omit(trainData), f=f)

# save model
save(lpqModel,
     file='./Models/Model 07 Model Tuning - Feature Selection/Model - lpq Features.Rdata',
     compress=TRUE)



par(mfrow=c(5, 5), mai=c(0.3,0.3,0.3,0.3))

# for each predictor
for (pred in predictors) {

    # make copy of data frame
    thisEnv <- env
    
    # now vary focal predictor from min to max value... all other predictors keep median value
    thisEnv[ , pred] <- seq(minPred[pred], maxPred[pred], length.out=100)
    
    # make prediction using this data frame
    predictionTargetBg <- predict(targetBgModel, thisEnv, type='cloglog')
    predictionlpq <- predict(lpqModel, thisEnv, type='cloglog')
    
    # plot
    plot(x=thisEnv[ , pred],
         y=predictionTargetBg,
         ylim=c(0, 1),
         xlab=pred,
         ylab='Suitability',
         main=pred,
         type='l',
         col='blue',
         lty='dotted',
         lwd=2
    )
    
    lines(x=thisEnv[ , pred], y=predictionlpq, col='red', lwd=1,
          lty='solid')
    
    legend('topright',
    legend=c('all features', 'lpq'),
    lty=c('dotted', 'solid'),
    col=c('blue', 'red'),
      lwd=2,
    cex=0.6,
    bty='n'
  )
    
    # add species' presences (top rug)
rug(records[ , pred], side=3, col='mediumseagreen')

# add background sites (bottom rug)
rug(targetBg[ , pred], side=1, col='red')

}

#########################  lqh prdict #########################  

lpqMap <- predict(climate[[predictors]], lpqModel,
  filename='./Models/Model 07 Model Tuning - Feature Selection/maxentPrediction1970to2000', 
  format='GTiff',  overwrite=TRUE,  type='cloglog', progress='text')

# plot
dev.off()

tiff(filename="./Output/Plots/All_vs_lqh.tif", width = 45, height = 15, units = 'cm', res = 300)

par(mfrow=c(1, 2), pty='s')

plot(targetBgMap, main='All features')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude)

plot(lpqMap, main='linear and quadratic features')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude)

dev.off()

#########################  autotune #########################  

# create output directory for model object and rasters
dir.create('./Models/Model 08 Model Tuning - Beta Parameter', recursive=TRUE, showWarnings=FALSE)

trainData <- cbind(presBg, na.omit(trainData))

# remotes::install_github('adamlilith/omnibus', dependencies=TRUE)

# to find the model that has the lowest AICc
# note: glmnet causes weird error if run more than once... restart R to run again

library(enmSdm)
library(glmnet)
# av og til rar feilmelding. bygg traindata p? nytt..
tunedModel <- enmSdm::trainMaxNet(data=trainData,  regMult=seq(.1, 1.5, by=.1),  verbose=TRUE)
#  0.6     170     lqh          3       23 -1346.428 2746.418  0.000000 1.000000e+00 8.661566e-01

# get just model object from output
save(tunedModel,
     file='./Models/Model 08 Model Tuning - Beta Parameter/Model.Rdata',
     compress=TRUE)


dev.off()
par(mfrow=c(5, 5), mai=c(0.3,0.3,0.3,0.3))

# for each predictor
for (pred in predictors) {
  
  # make copy of data frame
  thisEnv <- env
  
  # now vary focal predictor from min to max value
  # all other predictors keep median value
  thisEnv[ , pred] <- seq(minPred[pred], maxPred[pred], length.out=100)
  
  # make prediction using this data frame
  predictionUntuned <- predict(lpqModel, thisEnv, type='cloglog')
  predictionTuned <- predict(tunedModel, thisEnv, type='cloglog')
  
  # plot
  plot(thisEnv[ , pred],
       predictionUntuned,
       ylim=c(0, 1),
       xlab=pred,
       ylab='Suitability',
       main=pred,
       type='l',
       col='blue',
       lty='dashed',
       lwd=2
  )
  
  lines(x=thisEnv[ , pred], y=predictionTuned, col='black', lwd=2)
  
  legend('topright',
         legend=c(' lpq tuned', 'auto tuned'),
         lty=c('dashed', 'solid'),
         col=c('blue', 'black'),
         lwd=2,
         cex=1,
         bty='n'
  )
  
  # add species' presences (top rug)
  rug(records[ , pred], side=3, col='mediumseagreen')
  
  # add background sites (bottom rug)
  rug(targetBg[ , pred], side=1, col='red')
  
}

#########################  predict auttune  #########################  

# predict to raster
tunedMap <- predict(climate[[predictors]], tunedModel,
  filename='./Models/Model 08 Model Tuning - Beta Parameter/maxentPrediction1970to2000',
  format='GTiff', overwrite=TRUE, type='cloglog', progress='text')

# plot
dev.off()
#x11()

tiff(filename="./Output/Plots/lqh_vs_Auto_AIC.tif", width = 45, height = 15, units = 'cm', res = 300)

par(mfrow=c(1, 2), pty='s')

plot(lpqMap, main='lqp only Model')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude)

plot(tunedMap, main='Auto AIC tuned Model')
plot(wrld_simpl, add=TRUE, border='gray45')
#points(records$longitude, records$latitude)
dev.off()

#########################  vinner #########################  

# save.image(file="./Data/all_data_velutina.RData")

# convert to a df for plotting in two steps,
# First, to a SpatialPointsDataFrame
lpqMap_pts <- rasterToPoints(tunedMap, spatial = TRUE)

# Then to a 'conventional' dataframe
lpqMap_df  <- data.frame(lpqMap_pts)
head(lpqMap_df)

library(viridis)
memory.limit()
memory.limit(size=100000)

library(ggplot2)

#tiff("lpqMap_df.tiff", units="cm", width=30, height=30, res=600) 

ggplot() +
  geom_raster(data = lpqMap_df , aes(x = x, y = y, fill = maxentPrediction1970to2000 )) + 
  geom_point(data=all_presence, aes(x=longitude, y=latitude), col="black", size=.5, alpha=.5) +
  scale_fill_gradientn(colours=plasma(7),  breaks=c(0, 0.5, 1), labels=c(0 ,0.5, 1), limits=c(0, 1)) +
  labs(title = "Geithams",
       subtitle = "Current climate - Ecologically-chosen Bioclim predictors")

dev.off()

#########################  save #########################  

#save.image(file="./Data/all_data_velutina.RData")

  
  
