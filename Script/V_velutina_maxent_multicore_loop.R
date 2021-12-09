
# http://rstudio-pubs-static.s3.amazonaws.com/3911_c549d81ac2814caf9919dd0ec1b5829a.html

require(tidyverse)
require(dismo)
require(reshape)
library(rgbif)

memory.limit(size= 100000)
memory.limit()


splist  <- c("Vespa velutina", "Vespa mandarina")


# Resolve names using Global Names Resolv
species<-taxize::gnr_resolve(splist, best_match_only=T)

# make new vektor of species to check gbif backbone occurance 
splist<-species$matched_name

# fill in your gbif.org credentials 
user <- "danief"  
pwd <- "ZXasqw12" 
email <- "floedaniel@gmail.com" 

splist<- splist %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  purrr::imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() 


splist


eu <-occ_search(scientificName="Vespa velutina", continent="europe", limit = 99999, fields="minimal", hasGeospatialIssue=F, hasCoordinate=T)
eu

asia <-occ_search(scientificName="Vespa mandarina", continent="asia", limit = 99999, fields="minimal", hasGeospatialIssue=F, hasCoordinate=T)
asia

############################



library(spocc)

## Loading required package: ggplot2
out <- occ(query = splist$scientificname, from =c("gbif", "bison", "obis", "ala", "inat", "idigbio"), limit = 5000)

df_out <- occ2df(out)

as_tibble(df_out)

unique(df_out$prov)

df_out <- dplyr::select(df_out, name, longitude, latitude)

names(df_out)<-c("Species","Longitude", "Latitude")

df_out$Longitude <- as.numeric(df_out$Longitude)
df_out$Latitude <- as.numeric(df_out$Latitude)

library(scrubr)
df <-df_out %>%
  coord_impossible() %>%
  coord_incomplete() %>% 
  coord_unlikely()

df <- df %>% dedup() 

df$Species <- tolower(df$Species)
df <- dplyr::filter(df, !grepl("BOLD", Species))
df$Species <- word(df$Species, start=1, end=2)

df <- na.omit(df)

# ONLY ASIA 
df <-df %>%  filter(Longitude > 60) 

df <- df %>%  filter(tolower(Species) %in% tolower(splist))  

##################################################################################

# Bring in the env layers from the dismo package
predictors <- getData("worldclim", var="bio", res=5)
# view one predictor
predictors[[1]]
plot(predictors[[1]])
points(df)

############################# dopar ###########################################

# must have tis structure 
df <- df %>% dplyr::select(x=Longitude, y=Latitude, Species=Species)

df <- as.data.frame(df)



library(foreach)
library(doParallel)

# create list of unique counties to loop over 
sp_list <- unique(df$Species)

numCores <- detectCores()
numCores
registerDoParallel(numCores)  # use multicore, set to the number of our cores



max <- foreach::foreach (i=seq_along(sp_list), .packages = c("dismo", "rJava")) %dopar% {

  
  data=subset(df, df$Species==sp_list[i])
  
  occ <- data[, -3]
  
  # witholding a 20% sample for testing
  fold <- kfold(occ, k = 3)
  occtest <- occ[fold == 1, ]
  occtrain <- occ[fold != 1, ]
  
  me <- maxent(predictors, occtrain, removeDuplicates=T,
               args=c(
                 "cumulative", 
                 "outputformat=cloglog", 
                 "threads=3",
                 "betamultiplier=1.5", # !
                 
                 "autofeature=FALSE",
                 "product=FALSE",  # describe pairwise interactions between predictors 
                 "hinge=TRUE",    # combine linear and step functions
                 "threshold=TRUE", # simple step functions
                 "quadratic=TRUE", # features use squared predictor values
                 "linear=TRUE"     # simple linear coefficients for each predictor
                 
               ))
  
  # project into geographic space
  r <- predict(me, predictors)
  
}

stopImplicitCluster() 
##################################################################################

# stack all the outputs and name them
stack.out <- stack(max)
names(stack.out) <- sp_list
stack.out
#plot(stack.out)

mean <- calc(stack.out, fun = mean, na.rm = T)

plot(mean)

##################################################################################

for (i in 1:nlayers(stack.out)) {
  
  p=subset(df, df$Species==sp_list[i])
  
  plot(stack.out[[i]], main=sp_list[i])
  
  points(p[,1:2], type = "p", pch=16, size=1)
  
}
