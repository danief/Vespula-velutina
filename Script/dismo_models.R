# Library -----------------------------------------------------------------

library(dismo)
library(maptools)
data(wrld_simpl)

# Data --------------------------------------------------------------------

predictors <- raster::stack(
  list.files(
    path="./Data/Study region/",
    full.names=TRUE,
    pattern='.tif'
  )
)

load(file="./Data/velutina_presence.RData")
bradypus <- Europe
rm(list=setdiff(ls(), c("bradypus", "predictors")))


# data prep ---------------------------------------------------------------

presvals <- extract(predictors, bradypus[,3:2])
set.seed(0)
backgr <- randomPoints(predictors, 10000)
absvals <- extract(predictors, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))

ext <- extent(predictors)

set.seed(0)
group <- kfold(bradypus[,3:2], 5)
pres_train <- bradypus[,3:2][group != 1, ]
pres_test <- bradypus[,3:2][group == 1, ]

set.seed(10)
backg <- randomPoints(predictors, n=10000, ext=ext, extf = 1.25)
colnames(backg) = c( "longitude", "latitude")
group <- kfold(backg, 5)
backg_train <- backg[group != 1, ]
backg_test <- backg[group == 1, ]

r <- raster(predictors, 1)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(ext, add=TRUE, col='red', lwd=2)
points(backg_train, pch='-', cex=0.5, col='yellow')
points(backg_test, pch='-',  cex=0.5, col='black')
points(pres_train, pch= '+', col='green')
points(pres_test, pch='+', col='blue')

train <- rbind(pres_train, backg_train)
pb_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(predictors, train)
envtrain <- data.frame( cbind(pa=pb_train, envtrain) )
head(envtrain)
testpres <- data.frame( extract(predictors, pres_test) )
testbackg <- data.frame( extract(predictors, backg_test) )

# glm ---------------------------------------------------------------------

gm1 <- glm(pa ~ wc2.1_30s_bio_4 + wc2.1_30s_bio_5 + wc2.1_30s_bio_11 + wc2.1_30s_bio_17 , family = gaussian(link = "identity"), data=envtrain)
summary(gm1)

ge1 <- evaluate(testpres, testbackg, gm1)
ge1

pg <- predict(predictors, gm1, ext=ext)
tr <- threshold(ge2, 'spec_sens')


tiff(filename="./Output/Plots/glm.tif", width = 45, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(pg, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1), main='GLM/gaussian, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(pg > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
#points(pres_train, pch='+')
#points(backg_train, pch='-', cex=0.25)
dev.off()

# randomForest ------------------------------------------------------------

names(sdmdata)
library(randomForest)
## randomForest 4.6-14
## Type rfNews() to see new features/changes/bug fixes.
model <- factor(pa) ~ wc2.1_30s_bio_4 + wc2.1_30s_bio_5 + wc2.1_30s_bio_11 + wc2.1_30s_bio_17 
rf1 <- randomForest(model, data=na.omit(envtrain))
erf <- evaluate(testpres, testbackg, rf1)
erf

pr <- predict(predictors, rf1, ext=ext)
tr <- threshold(erf, 'spec_sens')


tiff(filename="./Output/Plots/rf.tif", width = 45, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(pr, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1), main='Random Forest, regression')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(pr > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
#points(pres_train, pch='+')
#points(backg_train, pch='-', cex=0.25)
dev.off()

# support vector ----------------------------------------------------------

library(kernlab)

svm <- ksvm(pa ~  wc2.1_30s_bio_4 + wc2.1_30s_bio_5 + wc2.1_30s_bio_11 + wc2.1_30s_bio_17 , data=envtrain)
esv <- evaluate(testpres, testbackg, svm)
esv

ps <- predict(predictors, svm, ext=ext)
tr <- threshold(esv, 'spec_sens')




tiff(filename="./Output/Plots/svm.tif", width = 45, height = 20, units = 'cm', res = 300)
par(mfrow=c(1,2))
plot(ps, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1), main='Support Vector Machine')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(ps > tr, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
#points(pres_train, pch='+')
#points(backg_train, pch='-', cex=0.25)
dev.off()

# Combining model predictions ---------------------------------------------

models <- stack(pg, pr, ps)
names(models) <- c("glm", "rf", "svm")

tiff(filename="./Output/Plots/three_models.tif", width = 25, height = 25, units = 'cm', res = 300)
plot(models,col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
dev.off()

auc <- sapply(list(ge1, erf, esv), function(x) x@auc)
w <- (auc-0.5)^2
m2 <- weighted.mean( models[[c("glm", "rf", "svm")]], w)

m2.min = cellStats(m2, "min")
m2.max = cellStats(m2, "max")

m2.scale <- (m2 - m2.min) / (m2.max - m2.min)

tiff(filename="./Output/Plots/weighted_mean_three_models.tif", width = 25, height = 25, units = 'cm', res = 300)
plot(m2.scale, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1), main='weighted mean of three models')
plot(wrld_simpl, add=TRUE, border='dark grey')
dev.off()

# end ---------------------------------------------------------------------

