
library(ecospat)
library(tidyverse)

plot(current$wc2.1_30s_bio_1)
points(subset(clean_data[3:2], longitude < 60))


# invasive range 
inv <- subset(clean_data, longitude < 60)
plot(current$wc2.1_30s_bio_1)
points(inv[3:2])

inv <- cbind(inv[3:2], raster::extract(current, inv[3:2]))
# inv <- raster::extract(current, inv[3:2], df=TRUE)
# setting random seed to always create the same
# random set of points for this example
set.seed(0)
backgr <- randomPoints(current, 10000)

absvals <- cbind(backgr, raster::extract(current, backgr))
colnames(absvals)[1:2] <- c("longitude", "latitude")
#absvals <- extract(current, backgr)
pb <- c(rep(1, nrow(inv)), rep(0, nrow(absvals)))
inv <- data.frame(cbind(pb, rbind(inv, absvals)))
head(inv)


# nattural range 
plot(current$wc2.1_30s_bio_1)
points(subset(clean_data[3:2], longitude > 60))

nat <- subset(clean_data, longitude > 60)

nat <- cbind(nat[3:2], raster::extract(current, nat[3:2]))

#nat <- raster::extract(current, nat[3:2])
# setting random seed to always create the same
# random set of points for this example
set.seed(0)
backgr <- randomPoints(current, 10000)
absvals <- cbind(backgr, raster::extract(current, backgr))
colnames(absvals)[1:2] <- c("longitude", "latitude")
pb <- c(rep(1, nrow(nat)), rep(0, nrow(absvals)))
nat <- data.frame(cbind(pb, rbind(nat, absvals)))
head(nat)

dim(inv)
dim(nat)

nat <-na.omit(nat)
inv <- na.omit(inv)
dim(inv)
dim(nat)

################################################################################
pca.env <- dudi.pca(rbind(nat, inv)[,4:22], scannf=F, nf=2)


ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)


# PCA scores for the whole study area
scores.globclim <- pca.env$li
# PCA scores for the species native distribution
scores.sp.nat <- suprow(pca.env, nat[which(nat[,1]==1),4:22])$li
# PCA scores for the species invasive distribution
scores.sp.inv <- suprow(pca.env, inv[which(inv[,1]==1),4:22])$li
# PCA scores for the whole native study area
scores.clim.nat <- suprow(pca.env, nat[,4:22])$li
# PCA scores for the whole invaded study area
scores.clim.inv <- suprow(pca.env, inv[,4:22])$li


grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.nat,
                                       sp=scores.sp.nat,
                                       R=100,
                                       th.sp=0)


grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv,
                                       sp=scores.sp.inv,
                                       R=100,
                                       th.sp=0)

D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D
D.overlap


eq.test <- ecospat.niche.equivalency.test(grid.clim.nat, grid.clim.inv, rep=1000, alternative = "greater")

ecospat.plot.overlap.test(eq.test, "D", "Equivalency")


sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv, rep=1000, alternative = "greater", rand.type=2)

ecospat.plot.overlap.test(sim.test, "D", "Similarity")

niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = 0.1)


ecospat.plot.niche.dyn(grid.clim.nat, 
                       grid.clim.inv, 
                       quant=0.25, interest=2,
                       title= "Niche Overlap",
                       name.axis1="PC1",
                       name.axis2="PC2")

ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)


# gridding the native niche
grid.clim.t.nat <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat, inv)[,4]),
                                         glob1=as.data.frame(nat[,4]),
                                         sp=as.data.frame(nat[which(nat[,1]==1),4]),
                                         R=1000, th.sp=0)

# gridding the invaded niche
grid.clim.t.inv <- ecospat.grid.clim.dyn(glob=as.data.frame(rbind(nat,inv)[,4]),
                                         glob1=as.data.frame(inv[,4]),
                                         sp=as.data.frame(inv[which(inv[,1]==1),4]),
                                         R=1000, th.sp=0)

t.dyn<-ecospat.niche.dyn.index (grid.clim.t.nat, grid.clim.t.inv, intersection=0.1)


ecospat.plot.niche.dyn(grid.clim.t.nat, grid.clim.t.inv, quant=0,
                       interest=2, 
                       title= "Niche Overlap",
                       name.axis1="Average temperature")

###################################################################################
