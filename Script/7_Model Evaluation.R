
# Model Evaluation 

randomBg

library(raster)
predPres <- raster::extract(lpqMap, cbind(records$longitude, records$latitude))
predBg   <- raster::extract(lpqMap, cbind(randomBg$longitude, randomBg$latitude))

minTrainPresThreshold <- min(predPres)
minTrainPresThreshold
## [1] 0.07852979
# sensitivity (true positive rate)
sum(predPres >= minTrainPresThreshold) / length(predPres)
## [1] 1
# specificity (true negative rate)
sum(predBg < minTrainPresThreshold) / length(predBg)

# calculate an "evaluation" object
eval <- evaluate(p=as.vector(predPres), a=as.vector(predBg), tr=seq(0, 1, by=0.01))

# see some evaluation statistics
eval 

equalSeSpThreshold <- eval@t[which.min(abs(eval@TPR - eval@TNR))]
equalSeSpThreshold

# sensitivity (true positive rate)
sum(predPres >= equalSeSpThreshold) / length(predPres)

# specificity (true negative rate)
sum(predBg < equalSeSpThreshold) / length(predBg)   

# sensitivity (true positive rate)
sum(predPres >= maxSeSpThreshold) / length(predPres)

# specificity (true negative rate)
sum(predBg < maxSeSpThreshold) / length(predBg) 



sensitivity <- FPR <- numeric()

# use 21 threshold values from 0 to 1
for (threshold in seq(0, 1, length.out=21)) {
  thisSens <- sum(predPres >= threshold) / length(predPres)
  thisFPR <- sum(predBg >= threshold) / length(predBg)
  sensitivity <- c(sensitivity, thisSens)
  FPR <- c(FPR, thisFPR)
}

# plot ROC curve
plot(FPR,
     sensitivity,
     type='l',
     xlab='False positive rate (=1 - specificity)',
     ylab='True positive rate (= sensitivity)',
     main='ROC Curve',
     lwd=2,
     col='blue'
)

abline(a=0, b=1, lty='dashed')

legend('bottomright',
       legend=c('ROC Curve', 'Random'),
       lwd=c(2, 1),
       lty=c('solid', 'dashed'),
       col=c('blue', 'black'),
       bty='n'
)

#########################################################################


