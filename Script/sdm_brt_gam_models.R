

library(sdm)

names(pres_df)
names(back_df)

names(back_df)[1] <- "longitude"
names(back_df)[2] <- "latitude"

presBg <- c(rep(1, nrow(pres_df)), rep(0, nrow(back_df)))
length(presBg)
presBg <- data.frame(presBg)
names(presBg) <- "pres"

trainData <- rbind(
  pres_df,
  back_df
)

df_all <- bind_cols(presBg, trainData)
df_all <- na.omit(df_all)

df_all <- df_all %>% select( pres,  longitude, latitude,
                             wc2.1_30s_bio_1, wc2.1_30s_bio_4,
                             wc2.1_30s_bio_5, wc2.1_30s_bio_11, wc2.1_30s_bio_17)

fold <- kfold(df_all, k=5)
prestest <-  df_all[fold == 1, ]
prestrain <- df_all[fold != 1, ]

getmethodNames('sdm')

data <- sdmData(formula= pres  ~. , train=prestrain[,c(1, 4:8)], test=prestest[,c(1, 4:8)])

m1 <- sdm(pres ~ . , data=data, methods=c("brt", "gam"))

tiff(filename="./Output/Plots/BRT_GAM_roc.tif", width = 30, height = 25, units = 'cm', res = 300)
roc(m1)
dev.off()

candidates <- c( "wc2.1_30s_bio_1", "wc2.1_30s_bio_4",
                 "wc2.1_30s_bio_5", "wc2.1_30s_bio_11", "wc2.1_30s_bio_17")

p1 <- predict(m1, newdata=rasss[[candidates]]) 

tiff(filename="./Output/Plots/BRT_GAM.tif", width = 45, height = 20, units = 'cm', res = 300)
plot(p1)
dev.off()

# ensemble ----------------------------------------------------------------

e2 <- ensemble(m1, newdata=rasss[[candidates]],  ncore=3,
               setting=list(method='weighted', stat='TSS'))

e2
roc(e2)

tiff(filename="./Output/Plots/BRT_GAM_weighted_ensemble.tif", width = 45, height = 20, units = 'cm', res = 300)
plot(e2, col=rev(terrain.colors(10)), breaks=seq(0, 1, by=0.1))
dev.off()
