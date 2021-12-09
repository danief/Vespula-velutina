

####################### Load library ######################################

library(raster)
library(maptools)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)

memory.limit()
memory.limit(size=100000)

# Set exstent  ------------------------------------------------------------

data(wrld_simpl)
eurasia <- ne_countries(continent = c('europe', "asia"), scale = "medium", returnclass = "sp" )
extent(bbox(eurasia))

extent <- c(-25, 180, -21.36904, 81.8542)

max_map_df <- fortify(eurasia) %>% mutate(long = ifelse(long < -100, long + 180*2, long))
ggplot() +
  geom_map(data=max_map_df, map=max_map_df, aes(map_id=id, x=long, y=lat), fill="black") +
  scale_x_continuous(limits = c(-25, 180), expand = c(0, 0))

############################  current worldclim  variables ######################################

referance <- raster("./Data/worldclim/bioclim/wc2.1_30s_bio_1.tif")

file <- "C:/Users/DAFL/Dropbox/vkm/Climat layers files/Beck_KG_V1/Beck_KG_V1_present_0p0083.tif"
rast <- raster(file)
projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
rast <-raster::resample(rast, referance, method="bilinear")
rast <- crop(rast, extent, snap = "out")
rast<- mask(rast, eurasia)

plot(rast)

# save raster
writeRaster(rast,
            "./Output/Raster/KG", 
            format='GTiff',
            datatype='INT2S',
            overwrite = TRUE)


kg_sub <- raster::extract(rast, cbind(velutina$longitude, velutina$latitude))
unique(kg_sub)
vec <-c(7, 8, 9, 14, 15, 26, 27)

kg_pts <- rasterToPoints(rast, spatial = TRUE)
kg_df  <- data.frame(kg_pts)
head(kg_df)

mapWorld <- borders("world", colour="gray50", fill=NA) # create a layer of borders


tiff("./Output/Plots/Köppen-Geiger.tif", units="cm", width=30, height=20, res=600) 

ggplot() +
  geom_raster(data = filter(kg_df, Beck_KG_V1_present_0p0083 %in% vec), aes(x = x, y = y, fill = factor(Beck_KG_V1_present_0p0083))) + 
  scale_fill_manual(na.value = "white", 
                    values  = c("#FEDB63", "#FEFE00", "#C7C700", "#C7FE4F", "#63FE31","#36C7FE" ,"#007C7C"),
                    breaks = c(7, 8, 9, 14, 15, 26, 27),
                    labels = c("BSk","Csa","Csb","Cfa","Cfb","Dfb","Dfc"),
                    name="Climate groups") +
  mapWorld +
  coord_quickmap(xlim = c(-15, 170), ylim = c(-10, 71)) +
  labs(title = "Köppen-Geiger climate classification",
       subtitle = "1‑km resolution current climate (1980–2016)")  +
  theme(panel.background = element_rect("white"),
        panel.grid = element_blank())

dev.off()


"#FFDC64", "#ffff00", "#C80000", "#50C8FF", "#FF5064","#37C8C8" ,"#7D7D00"
"BSk","Csa","Csb","Cfa","Cfb","Dfb","Dfc"

# disse haer feil hex 
#7:  BSk  Arid, steppe, cold                    [255 220 100] #FFDC64
#8:  Csa  Temperate, dry summer, hot summer     [255 255 0] #ffff00
#9:  Csb  Temperate, dry summer, warm summer    [200 200 0] #C80000
#14: Cfa  Temperate, no dry season, hot summer  [200 255 80] #50C8FF
#15: Cfb  Temperate, no dry season, warm summer [100 255 80] #FF5064
#26: Dfb  Cold, no dry season, warm summer      [55 200 255] #37C8C8
#27: Dfc  Cold, no dry season, cold summer      [0 125 125] #7D7D00





# future ------------------------------------------------------------------

referance <- raster("./Data/worldclim/bioclim/wc2.1_30s_bio_1.tif")

file <- "C:/Users/DAFL/Dropbox/vkm/Climat layers files/Beck_KG_V1/Beck_KG_V1_future_0p0083.tif"
rast <- raster(file)
projection(rast) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
rast <-raster::resample(rast, referance, method="bilinear")
rast <- crop(rast, extent, snap = "out")
rast<- mask(rast, eurasia)

plot(rast)

kg_sub <- raster::extract(rast, cbind(velutina$longitude, velutina$latitude))
sort(unique(kg_sub))
vec <-c(6, 7, 8, 9, 14, 15, 18, 26)

kg_pts_fut <- rasterToPoints(rast, spatial = TRUE)
kg_df_fut  <- data.frame(kg_pts_fut)
head(kg_df_fut)

# 6:  BSh  Arid, steppe, hot                     [245 165 0]
#7:  BSk  Arid, steppe, cold                    [255 220 100] #FFDC64
#8:  Csa  Temperate, dry summer, hot summer     [255 255 0] #ffff00
#9:  Csb  Temperate, dry summer, warm summer    [200 200 0] #C80000
#14: Cfa  Temperate, no dry season, hot summer  [200 255 80] #50C8FF
#15: Cfb  Temperate, no dry season, warm summer [100 255 80] #FF5064
# 18: Dsb  Cold, dry summer, warm summer         [200 0 200]
#26: Dfb  Cold, no dry season, warm summer      [55 200 255] #37C8C8

tiff("./Output/Plots/Future_Köppen-Geiger.tif", units="cm", width=30, height=20, res=600) 

ggplot() +
  geom_raster(data = filter(kg_df_fut, Beck_KG_V1_future_0p0083 %in% vec), aes(x = x, y = y, fill = factor(Beck_KG_V1_future_0p0083))) + 
  scale_fill_manual(na.value = "white", 
                    values  = c("#F5A301", "#FEDB63", "#FEFE00", "#C7C700", "#C7FE4F", "#63FE31", "#C600C7", "#36C7FE"),
                    breaks = c(6, 7, 8, 9, 14, 15, 18, 26),
                    labels = c("BSh", "BSk","Csa","Csb","Cfa","Cfb","Dsb","Dfb"),
                    name="Climate groups") +
  mapWorld +
  coord_quickmap(xlim = c(-15, 170), ylim = c(-10, 71)) +
  labs(title = "Köppen-Geiger climate classification",
       subtitle = "1‑km resolution projected future conditions (2071–2100) under climate change)")  +
  theme(panel.background = element_rect("white"),
        panel.grid = element_blank())

dev.off()









