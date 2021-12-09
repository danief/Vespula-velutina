# Load Packages -----------------------------------------------------------

library(tidyverse)
library(kgc)
library(dplyr)
library(taxize)
library(rgbif)
library(spocc)
library(scrubr)
library(raster)
library(spocc)

# Species names  ----------------------------------------------------------

# Subspecies 
sub<-name_suggest(q='Vespa velutina', rank='subspecies')
sub
# Species 
splist<-c("Vespa velutina")

# bind 
splist<- c(splist, sub$data$canonicalName)

# Resolve names using Global Names Resolv
Species<-taxize::gnr_resolve(splist, best_match_only=T)

Species

# make new vektor of species to check gbif backbone occurance 
splist<-Species$matched_name
# split list for plotting 

# fill in your gbif.org credentials 
user <- "danief"  
pwd <- "ZXasqw12" 
email <- "floedaniel@gmail.com" 

# use only eccpeted names 
test<- splist %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  purrr::imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows()  %>% # combine all data.frames into one
  filter(matchtype == "EXACT" & status == "ACCEPTED")

test

# bind and download 
dat2 <- occ_data(test$specieskey, limit=100000)
dd <- data.table::rbindlist(lapply(dat2, function(z) z$data), fill = TRUE, use.names = TRUE)
dd

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_point(data=dd, aes(x=decimalLongitude, y=decimalLatitude, col=species)) 

# georef missing lat long with geocode ------------------------------------

georef <- subset(dd, (is.na(decimalLongitude) | is.na(decimalLatitude)) & ! is.na(locality) )
georef

library(ggmap)
ggmap::register_google(key = "AIzaSyA46QZE1JSsQy_r26YnXbJS5py2rFAT9v0")

geolatlong<-geocode(georef$locality)

# 
t <- georef %>% as_tibble()

# change values in column with new lat long 
t[,4:3] <- geolatlong

# check 
t

# bind new values to original df. 
# dont care about doubles because NA`s will disrepair anyway
dim(dd)
dim(t)
complete <- bind_rows(t, dd)
dim(complete)


# Select and plot ---------------------------------------------------------

# selct
dat <- complete %>% dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                            gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, basisOfRecord, 
                            institutionCode, datasetName)

# plot 
library(maptools)
library(ggplot2)
data(wrld_simpl)

wrld_simpl@data$id <- wrld_simpl@data$NAME
wrld <- subset(wrld_simpl, wrld_simpl@data$NAME != "Antarctica") # ta vekk Antarctica
wrld <- fortify(wrld)

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_point(data=dat, aes(x=decimalLongitude, y=decimalLatitude, col=species)) 

################################## Clean GBIF data #########################################################################################

# Improving data quality using GBIF meta-data -----------------------------

hist(dat$coordinateUncertaintyInMeters / 1000, breaks = 20)

dat <- dat %>%  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# Delete NA coords 
dat <- dat[which(!is.na(dat$decimalLatitude) & !is.na(dat$decimalLongitude)),]


# Use CoordinateCleaner to automatically flag problematic records ---------
library(CoordinateCleaner)

dat<-cc_dupl(x = dat,
                lon = "decimalLongitude",
                lat = "decimalLatitude", 
                species = "species",
                value = "clean")

dat$countryCode <-  countrycode::countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')

#flag problems
flags <- clean_coordinates(x = dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions", "zeros", "countries", "seas"))


summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

#Exclude problematic records
dat_cl <- dat[flags$.summary,]

#The flagged records
dat_fl <- dat[!flags$.summary,]

# Removes or flags records that are temporal outliers based on interquantile ranges.
# Temporal outliers
flags <- cf_age(x = dat_cl,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                taxon = "species", 
                min_age = "year", 
                max_age = "year", 
                value = "flagged")


dat_cl[!flags, "year"]

dat_cl <- dat_cl[flags, ]

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_point(data=dat_cl, aes(x=decimalLongitude, y=decimalLatitude, col=species)) 


# basis of records --------------------------------------------------------

table(dat_cl$basisOfRecord)

# no neede 
dat_cl <- filter(dat_cl, basisOfRecord != "FOSSIL_SPECIMEN")

#Individual count
table(dat_cl$individualCount)

table(dat_cl$year)

hist(dat_cl$year)

dat_cl <- dat_cl %>%  filter(year > 1945) 

table(dat_cl$family)

clean_gbif <- dat_cl

clean_gbif

# occ ---------------------------------------------------------------------

# check additional  multipple data sources 
out <- occ(query= splist, from=c("bison","idigbio", "inat", "obis", "ecoengine", "ala"),  has_coords = T, limit=5000)

# inspect individual elements
out$inat
out$ala

## Coerce to combined data.frame, selects minimal set of
## columns (name, lat, long, provider, date, occurrence key)
out_df<-occ2df(out)
# only 10
 
 # select, fix names and bind 
 out_df <- dplyr::select(out_df, name, latitude, longitude)
 
 out_df$latitude<-as.numeric( out_df$latitude)
 out_df$longitude<-as.numeric( out_df$longitude)
 
 out_df
 
 # delete all NA
 data<- out_df %>% filter(!is.na(latitude)) %>% filter(!is.na(longitude))
 
 data
 
 data_cl<-cc_dupl(x = data,
                 lon = "longitude",
                 lat = "latitude", 
                 species = "name",
                 value = "clean")
 
 flags <- clean_coordinates(x = data_cl,
                            lon = "longitude",
                            lat = "latitude",
                            species = "name",
                            tests = c("capitals", "centroids", "equal", "gbif", "institutions", "zeros", "seas"))
 
 
 
 # 
 summary(flags)
 plot(flags, lon = "longitude", lat = "latitude", size=2)
 
 # Exclude problematic records
 clean_data <- data_cl[flags$.summary,]
 
 # The flagged records
 flaged_data <- data_cl[!flags$.summary,]
 
 
 ggplot() +
   geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
   geom_point(data=clean_data, aes(x=longitude, y=latitude, col=name), size=2) 
 
 
 clean_occ_data <- clean_data
 
 ####################################################################################################################
 
 clean_gbif<- clean_gbif %>%
   filter(!countryCode=="SWE") %>% 
   filter(!countryCode=="GBR") %>% 
   filter(!countryCode=="FIN") %>% 
   filter(!countryCode=="DNK") 
   
 
 clean_occ_data
 
 data <- bind_rows(clean_occ_data, dplyr::select(clean_gbif, name=species, longitude=decimalLongitude, latitude =decimalLatitude))
 
 unique(data$name)
 
 ########################################################################################################################
 
# get rid of DNA codes
data <- dplyr::filter(data, !grepl("BOLD", name))
 
# set all to lover case
data$name <- tolower(data$name)
 
# grep only the first 3 words
#data$name  <-stringr::word(data$name, 1, 2)
 
unique(data$name)
 
match <- tolower(splist)
 
splist
 
# data <- data %>% filter(name %in% match)
 
# plot final data 
ggplot() +
   geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
   geom_jitter(data=data, aes(x=longitude, y=latitude, col=name)) 
 
all_presence <- data

# Europe ------------------------------------------------------------------

Europe <- data %>% filter(longitude < 60)


Europe <- Europe %>% 
  cc_outl(lon = "longitude", 
          lat = "latitude",
          species = "name",
          method="mad",
          mltpl=.1,
          value = "clean",
          thinning=T)

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_jitter(data=Europe, aes(x=longitude, y=latitude, col=name)) 

# nat range ---------------------------------------------------------------
 
ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_point(data=filter(data, longitude  > 60) , aes(x=longitude, y=latitude, col=name)) +
  coord_cartesian()

 require(raster)
 china <- getData('GADM', country='china', level=0)
 plot(china)
 
 china@data$id <- rownames(china@data)
 mapa.df     <- fortify(china)
 mapa.df     <- plyr::join(mapa.df, china@data, by="id")
 
 head(mapa.df)
 
 summary(mapa.df)
 
 mapa.df <- mapa.df %>%  filter(piece ==1)
 
 
 china <- mapa.df
 
 # fra nedre venstre mot h?yre 
 range <- Polygon(cbind(c(china$long), # x
                        c(china$lat)))  # y
 
 

 range <- Polygons(list(range), ID = c("A"))
 wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
 range <- SpatialPolygons(list(range), proj4string = CRS(wgs84))
 ## Warning in showSRID(uprojargs, format = "PROJ", multiline = "NO"): Discarded datum WGS_1984 in CRS definition,
 ##  but +towgs84= values preserved
 df <- data.frame(species = c("Vespa velutina"), row.names = c("A"))
 nat_range <- SpatialPolygonsDataFrame(range, data = as.data.frame(df))
 
 
 # Visualize range
 plo <- fortify(nat_range)
 ## Regions defined for each Polygons
 
 ggplot() +
   borders("world", colour="gray50")+
   geom_polygon(data = plo, aes(x = long, y = lat, group = id), fill="green", alpha=.5)+
   geom_point(data=data, aes(x=longitude, y=latitude)) +
   theme_bw()+
   theme(legend.position = "none",
         axis.title = element_blank())
 
 
 # run cc_iucn()
 range_flags <- cc_iucn(x = dat_cl,
                        range = nat_range,
                        lon = "decimalLongitude",
                        lat = "decimalLatitude",
                        value = "flagged")
 
 
 dat_fin <- dat_cl[range_flags, ]
 
 ggplot() +
   borders("world", colour="gray50")+
   geom_polygon(data = plo, aes(x = long, y = lat, group = id), fill="green", alpha=.5)+
   geom_point(data=data, aes(x=longitude , y=latitude ), col="red") +
   geom_point(data=dat_fin, aes(x=decimalLongitude, y=decimalLatitude)) +
   theme_bw()+
   theme(legend.position = "none",
         axis.title = element_blank())
 
 dim(dat_fin)
 
 china_presence <- dat_fin %>% dplyr::select(name=species, latitude=decimalLatitude, longitude=decimalLongitude)
 

# Final data  -------------------------------------------------------------

all_presence
 
 ggplot() +
   geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
   geom_jitter(data=all_presence, aes(x=longitude, y=latitude, col=name)) 
 
china_presence

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_jitter(data=china_presence, aes(x=longitude, y=latitude, col=name)) 

Europe

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_jitter(data=Europe, aes(x=longitude, y=latitude, col=name)) 

velutina <- bind_rows(china_presence, Europe)

ggplot() +
  geom_map(data=wrld, map=wrld, aes(map_id=id, x=long, y=lat), fill="black") +
  geom_jitter(data=velutina, aes(x=longitude, y=latitude, col=name)) 

# Save --------------------------------------------------------------------

getwd()
# save.image("./Data/velutina_presence.RData")
# load("./Data/velutina_presence.RData")

# End ---------------------------------------------------------------------


