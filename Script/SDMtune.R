
#https://cran.r-project.org/web/packages/SDMtune/vignettes/var-selection.html 

# SDMtune -----------------------------------------------------------------

library(SDMtune)
library(zeallot)
library(dismo)

# Load data  -------------------------------------------------------------

load("./Data/velutina_presence.RData")

velutina <- data.frame(Europe[,3:2])

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_jitter(data=velutina, aes(x=longitude, y=latitude), col="red") 

# remove everything except
rm(list=setdiff(ls(), c("velutina")))

current <- stack(
  list.files(
    path="./Data/Study region/",
    full.names=TRUE,
    pattern='.tif'
  )
)

# Target BP ---------------------------------------------------------------

# create circles using an arbitrary radius of 50 km =  d=50000
# 50 km to small
x     <- circles(velutina, d=1000000, lonlat=TRUE) 
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
points(velutina, cex=1, pch=20, col="red")

targetSites <-xy
randomBgSites <- targetSites

names(randomBgSites) <- names(velutina)

# Data  -------------------------------------------------------------------

data <- prepareSWD(species = "velutina", p = velutina, a = randomBgSites, env = current)

#c(train, test) %<-% trainValTest(data, test = 0.2, only_presence = TRUE, seed = 25)
c(train, val, test) %<-% trainValTest(data, val = 0.2, test = 0.2, only_presence = TRUE, seed = 12345)

# Set parameters ----------------------------------------------------------

# Train cross validation model
h <- list(reg = seq(0.1, 3, 0.1), fc = c("l", "lq", "lp", "lh", "lt", "qp", "qh", "qt", "ph", "pt", "ht", "lqp", "lqph", "lqpt", "lqpht"))
folds <- randomFolds(data, k = 10, only_presence = TRUE, seed = 25)

# Train model -------------------------------------------------------------

cv_model <- train("Maxent", data = data, folds = folds)

# Data-driven variable selection ------------------------------------------

plotCor(data, method = "spearman", cor_th = 0.7)

corVar(data, method = "spearman", cor_th = 0.7)

selected_variables_model <- varSel(cv_model, metric = "auc", test = test,
                                   bg4cor = test, method = "spearman",
                                   cor_th = 0.7, permut = 10)

cat("Testing AUC before: ", auc(cv_model, test = test))
cat("Testing AUC after: ", auc(selected_variables_model, test = test))

varImp(selected_variables_model, permut = 10)

# Test all the possible combinations with gridSearch ----------------------

gs <- gridSearch(selected_variables_model, hypers = h, metric = "auc", test = val)
head(gs@results[order(-gs@results$test_AUC), ])  # Best combinations
plot(gs, title = "10-fold gridSearch")

index <- which.max(gs@results$test_AUC)  # Index of the best model in the experiment
new_train <- gs@models[[index]]@data  # New train dataset containing only the selected variables
merged_data <- mergeSWD(new_train, test, only_presence = TRUE) # Merge only presence data

# Final model  ------------------------------------------------------------

final_model <- train("Maxnet", data = merged_data, fc = gs@results[index, 1], reg = gs@results[index, 2])
auc(final_model, test = test)

# Subset variables  -------------------------------------------------------

drop<-c( "wc2.1_30s_bio_10", "wc2.1_30s_bio_11", 
         "wc2.1_30s_bio_12", "wc2.1_30s_bio_14", 
         "wc2.1_30s_bio_15", "wc2.1_30s_bio_16", 
         "wc2.1_30s_bio_18", "wc2.1_30s_bio_19", 
         "wc2.1_30s_bio_3", "wc2.1_30s_bio_5", 
         "wc2.1_30s_bio_6", "wc2.1_30s_bio_7")

sub<- dropLayer(current, drop)

# Predict and plot final model  -------------------------------------------

map <- predict(final_model, data = sub, type = "cloglog", file = "my_map", format = "GTiff", hr=T, progress="text")
plotPred(map)

tiff(filename="./Output/Plots/sdmtune_asia_and_europe_data.tif", width = 25, height = 25, units = 'cm', res = 300)

plot(map, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))

dev.off()

# Plot tgreshold map ------------------------------------------------------

ths <- thresholds(final_model, type = "cloglog")
thresholds(final_model, type = "cloglog")
memory.limit()
memory.limit(size=300000)

tiff(filename="./Output/Plots/sdmtune_asia_and_europe_data_thresholds.tif", width = 25, height = 25, units = 'cm', res = 300)

plotPA(map, th = ths[3, 2], hr=T)

dev.off()

# Model report ------------------------------------------------------------

save.image(file="./Data/all_data_velutina.RData")

modelReport(final_model, type = "cloglog", folder = "velutina SDMtune", test = test,
            response_curves = TRUE, only_presence = TRUE, jk = F,
            env = sub)

# Variable importance -----------------------------------------------------

model@model@results

vi <- maxentVarImp(model)
vi

# variable importance
plotVarImp(vi[, 1:2])

# or the permutation importance with:
plotVarImp(vi[, c(1,3)])


# Permutation importance --------------------------------------------------

maxnet_model <- train("Maxnet", data = train)

vi_maxnet <- varImp(maxnet_model, permut = 5)
vi_maxnet

plotVarImp(vi_maxnet)

# Compute the permutation importance
vi_maxent <- varImp(model, permut = 10)
# Print it
vi_maxent
# Compare with Maxent output
maxentVarImp(model)


# Jackknife test for variable importance ----------------------------------

jk <- doJk(maxnet_model, metric = "auc", test = test)
jk

plotJk(jk, type = "train", ref = auc(maxnet_model))

# and the Jackknife test for the testing AUC:
plotJk(jk, type = "test", ref = auc(maxnet_model, test = test))

# Response curves ---------------------------------------------------------

plotResponse(model, var = "wc2.1_30s_bio_1", type = "cloglog", only_presence = TRUE,
             marginal = FALSE, rug = TRUE)

plotResponse(cv_model, var = "wc2.1_30s_bio_1", type = "cloglog", only_presence = TRUE,
             marginal = TRUE, fun = mean, rug = TRUE)


