
# cross-validation --------------------------------------------------------

backup <- records

# calculate k-folds for presences and background sites
kPres  <- kfold(x=records, k=10)
kBg    <- kfold(x=targetBg, k=10)

head(kPres)

head(kBg)

# map
plot(climate$wc2.1_30s_bio_1, main='k-fold #1')
points(records$longitude, records$latitude)
points(records$longitude[kPres==1],
       records$latitude[kPres==1],
       bg='red',
       pch=21
)

legend('topright',
       legend=c('Training presence', 'Test presence'),
       pch=c(1, 16),
       col=c('black', 'red'),
       bg='white',
       cex=0.8
)


# dopar -------------------------------------------------------------------

# create output directory for model object and rasters
dir.create('./Models/Model 09 Model Evaluation - Random K-Folds', recursive=TRUE, showWarnings=FALSE)

# for storing AUC and CBI
aucRandom <- cbiRandom <- rep(NA, 5)

library(omnibus)
library(enmSdm)

library(doParallel)  #Foreach Parallel Adaptor 
library(foreach)     #Provides foreach looping construct

#Define how many cores you want to use
UseCores <- detectCores() -1
#Register CoreCluster
cl       <- makeCluster(UseCores)
registerDoParallel(cl)


#Use foreach loop and %dopar% command
foreach(i=1:10, .packages = c("omnibus", "enmSdm")) %dopar% {
  
  say('K-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0 for
  # presence/background... using only points not in this k-fold
  envData <- rbind(
    records[kPres!=i, predictors],
    targetBg[kBg!=i, predictors]
  )
  
  presBg <- c(rep(1, sum(kPres!=i)), rep(0, sum(kBg!=i)))
  
  trainData <- cbind(presBg, envData)
  
  # tuned model: just using one regularization value for expdiancy!
  model <- enmSdm::trainMaxNet(
    data=na.omit(trainData),
    regMult=seq(.1, 1.5, by=.1),
    classes='lpqht',
    verbose=T
  )
  
  # model <- model$model
  
  save(model,
       file=paste0('./Models/Model 09 Model Evaluation - Random K-Folds/Model ', i, '.Rdata'),
       compress=TRUE)
  
  # predict to presences and background sites
  predPres <- raster::predict(model, newdata=records[kPres==i, ], type='cloglog')
  predBg   <- raster::predict(model, newdata=randomBg[kBg==i, ],  type='cloglog')
  
  # evaluate and remember result
  thisEval <- evaluate(p=as.vector(predPres), a=as.vector(predBg))
  
  thisAuc <- thisEval@auc
  thisCbi <- contBoyce(pres=predPres, bg=predBg)
  
  say(': AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2))
  
  aucRandom[i] <- thisAuc
  cbiRandom[i] <- thisCbi
  
  
}

#end cluster
stopCluster(cl)

# AUC model --------------------------------------------------------------

say('Mean AUC:', round(mean(aucRandom), 2))
##  Mean AUC:0.84
say('Mean CBI:', round(mean(cbiRandom), 2))

#Lobo et al  2008
# 0.5-0.6, fail;
# 0.6-0.7, poor;
# 0.7-0.8, fair; 
# 0.8-0.9, good;
# 0.9-1, excellent.


# Predict -----------------------------------------------------------------

## predict to presences and background sites of non-k-folded sites
predPres <- predict(tunedModel, records, type='cloglog')
predBg <- predict(tunedModel, randomBg, type='cloglog')

# evaluate and remember result
thisEval <- evaluate(p=as.vector(predPres), a=as.vector(predBg))

aucResub <- thisEval@auc
cbiResub <- contBoyce(pres=predPres, bg=predBg, numBins=100)

# PLot --------------------------------------------------------------------

par(mfrow=c(1, 2), pty='m')
barplot(c(aucResub, mean(aucRandom)),
        names.arg=c('resubstituted', 'k-fold'),
        cex.lab=0.7,
        col=c('white', 'cornflowerblue'), ylab='AUC', main='AUC'
)


barplot(c(cbiResub, mean(cbiRandom)),
        names.arg=c('re-substituted', 'k-fold'),
        cex.lab=0.7,
        col=c('white', 'cornflowerblue'), ylab='CBI', main='CBI'
)






# calculate g-folds -------------------------------------------------------

gPres <- geoFold(
  x=records,
  k=3,
  minIn=5,
  minOut=10,
  longLat=c('longitude', 'latitude')
)

# maps
par(mfrow=c(1, 3), pty='s')
for (i in 1:3) {
  
  plot(climate$wc2.1_30s_bio_1, main=paste0('g-fold #', i))
  points(records$longitude, records$latitude)
  points(records$longitude[gPres==i],
         records$latitude[gPres==i],
         bg='red',
         pch=21
  )
  
  legend('topright',
         legend=c('Training presence', 'Test presence'),
         pch=c(1, 16),
         col=c('black', 'red'),
         bg='white',
         cex=0.8
  )
  
}





# initialize vectors to store g-fold assignments
gTrainBg <- rep(NA, nrow(na.omit(targetBg)))
gTestBg <- rep(NA, nrow(na.omit(randomBg)))

# for each TRAINING background site, find closest TRAINING presence and use it's g-fold designation
for (i in 1:nrow(na.omit(targetBg))) {
  dists <- geosphere::distCosine(cbind(na.omit(targetBg)$longitude[i], na.omit(targetBg)$latitude[i]), cbind(records$longitude, records$latitude))
  closest <- which.min(dists)
  gTrainBg[i] <- gPres[closest]
}

# for each TEST background site, find closest TEST presence and use it's g-fold designation
for (i in 1:nrow(na.omit(randomBg))) {
  dists <- geosphere::distCosine(cbind(na.omit(randomBg)$longitude[i], na.omit(randomBg)$latitude[i]), cbind(records$longitude, records$latitude))
  closest <- which.min(dists)
  gTestBg[i] <- gPres[closest]
}



# create output directory
dirCreate('./Models/Model 10 Model Evaluation - Geographic K-Folds')

# for storing model evaluation metrics for each k-fold
aucGeog <- cbiGeog <- numeric()

# cycle through each k-fold
for (i in 1:3) {
  
  say('G-fold ', i, ':', post=0)
  
  # make training data frame with predictors and vector of 1/0
  # for presence/background
  envData <- rbind(
    records[gPres!=i, predictors],
    targetBg[gTrainBg!=i, predictors]
  )
  
  presBg <- c(rep(1, sum(gPres!=i)), rep(0, sum(gTrainBg!=i)))
  
  trainData <- cbind(presBg, envData)
  
  # tuned model
  model <- trainMaxNet(
    data=trainData,
    regMult=1,
    classes='lpq',
    out='model'
  )
  
  # save
  save(model,
       file=paste0('./Models/Model 10 Model Evaluation - Geographic K-Folds/Model ', i, '.Rdata'),
       compress=TRUE)
  
  # predict to presences and background sites
  predPres <- predict(model, records[gPres==i, ], type='cloglog')
  predBg <- predict(model, randomBg[gTestBg==i, ], type='cloglog')
  
  # evaluate and remember result
  thisEval <- evaluate(p=as.vector(predPres), a=as.vector(predBg))
  
  thisAuc <- thisEval@auc
  thisCbi <- contBoyce(pres=predPres, bg=predBg, numBins=100)
  
  say(': AUC (geog) = ', round(thisAuc, 2), ' | CBI (geog) = ', round(thisCbi, 2))
  
  aucGeog <- c(aucGeog, thisAuc)
  cbiGeog <- c(cbiGeog, thisCbi)
  
}
##  G-fold 1:  : AUC (geog) = 0.8 | CBI (geog) = 0.94 
##  G-fold 2:  : AUC (geog) = 0.96 | CBI (geog) = 0.6 
##  G-fold 3:  : AUC (geog) = 0.71 | CBI (geog) = 0.75
say('Mean AUC (geog): ', round(mean(aucGeog), 2))
##  Mean AUC (geog): 0.82
say('Mean CBI (geog): ', round(mean(cbiGeog), 2))
##  Mean CBI (geog): 0.76



  # plot
par(mfrow=c(1, 2), pty='m')
barplot(c(aucResub, mean(aucRandom), mean(aucGeog)),
        names.arg=c('resubstituted', 'k-fold', 'geo-fold'),
        cex.lab=0.7,
        col=c('white', 'cornflowerblue', 'forestgreen'),
        ylab='AUC',
        main='AUC'
)

barplot(c(cbiResub, mean(cbiRandom), mean(cbiGeog)),
        names.arg=c('resubstituted', 'k-fold', 'geo-fold'),
        cex.lab=0.7,
        col=c('white', 'cornflowerblue', 'forestgreen'),
        ylab='CBI',
        main='CBI'
)

# save 
save.image(file="./Data/all_data_velutina.RData")

