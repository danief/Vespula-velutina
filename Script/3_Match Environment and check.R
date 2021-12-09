
# Data --------------------------------------------------------------------
memory.limit()
memory.limit(325750)
load(file="./Data/velutina_presence.RData")
velutina <-china_presence
rm(list=setdiff(ls(), c("velutina", "current")))

env <- current

# Matche presence to environment ------------------------------------------

# match species' records with environment at each location
envSpecies <- raster::extract(env, cbind(velutina$longitude, velutina$latitude))
envSpecies <- as.data.frame(envSpecies)
records <- cbind(velutina, envSpecies)
head(records)

#  remove any records that might have NA's as environmental data 
if (any(is.na(rowSums(envSpecies)))) records <-  records[-which(is.na(rowSums(envSpecies))), ]

par(mfrow=c(1, 2))
hist(records$wc2.1_30s_bio_1, breaks=20, xlab='Celcius', main='BIO1 = Annual Mean Temperature')
hist(records$wc2.1_30s_bio_12, breaks=20, xlab='mm', main='BIO12 = Annual Precipitation')

records %>% 
  ggplot() +
  geom_point(aes(wc2.1_30s_bio_1, wc2.1_30s_bio_12))

dev.off()
# Low temp records may be outliers
records %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=wc2.1_30s_bio_1, col=wc2.1_30s_bio_12))  

# records in the asia? at a much lower temp
records %>% 
filter(longitude > 60) %>% 
  ggplot() +
  geom_point(aes(x=longitude, y=as.factor(wc2.1_30s_bio_1), col=as.factor(wc2.1_30s_bio_1)))+
  theme(legend.position = "none")


# Consider removing US presence points
ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_point(data=records, aes(x=longitude, y=latitude, col=as.factor(wc2.1_30s_bio_1)), size=2) 



ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_point(data=filter(records, longitude>0), aes(x=longitude, y=latitude, col=as.factor(wc2.1_30s_bio_1)), size=2) 

library(enmSdm)
recordsNoDups <- elimCellDups(x=records, r=env$wc2.1_30s_bio_1, longLat=c('longitude', 'latitude'))

nrow(records)
## [1] 333
nrow(recordsNoDups)

par(mfrow=c(2, 2))
hist(records$wc2.1_30s_bio_1, breaks=20, xlab='Celcius', main='BIO1 = Annual Mean Temperature')
hist(records$wc2.1_30s_bio_12, breaks=20, xlab='mm', main='BIO12 = Annual Precipitation')

hist(recordsNoDups$wc2.1_30s_bio_1, breaks=20, xlab='Celcius', main='BIO1 = Annual Mean Temperature - no dups')
hist(recordsNoDups$wc2.1_30s_bio_12, breaks=20, xlab='mm', main='BIO12 = Annual Precipitation - no dups')


records <- recordsNoDups
getwd()
saveRDS(records, "./Data/Species_presence.rds")

