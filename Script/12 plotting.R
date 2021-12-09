
library(viridis)
library(ggplot2)
library(raster)
memory.limit()
memory.limit(size=100000)

mapWorld <- borders("world", colour="gray50", fill=NA) # create a layer of borders


future_map <- raster("./Output/Raster/future_map.tif")
plot(future_map)

present_map <- raster("./Output/Raster/final_map_2.tif")
plot(present_map)

diff <- future_map-present_map

plot(diff)

# Future  ------------------------------------------------------------------

future_map <- raster("./Output/Raster/future_map.tif")

future_pts <- rasterToPoints(future_map, spatial = TRUE)

# Then to a  dataframe
future_df  <- data.frame(future_pts)

head(future_df)

tiff("./Output/Plots/future_europe_subset_lh_b.7_auc089_10fold.tif", units="cm", width=30, height=20, res=600) 

ggplot() +
  geom_raster(data = future_df , aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradientn(colors=rev(terrain.colors(10)),  breaks=seq(0, 1, by=0.1), labels=seq(0, 1, by=0.1), limits=c(0, 1)) +
  mapWorld +
  coord_quickmap(xlim = c(-15, 170), ylim = c(-10, 71)) +
  labs(title = "MaxEnt model for Vespa velutina",
       subtitle = "5-fold cross-validation with future climate (MPI-ESM-LR rcp85 2070)")  +
  theme(panel.background = element_rect("white"),
        panel.grid = element_blank())

dev.off()

# present ------------------------------------------------------------------

present_pts <- rasterToPoints(final_map_2, spatial = TRUE)

# Then to a dataframe
present_df  <- data.frame(present_pts)
head(present_df)

tiff("./Output/Plots/Final_europe_subset_lq_b1_auc085_df.tiff", units="cm", width=30, height=20, res=600) 

ggplot() +
  geom_raster(data = present_df , aes(x = x, y = y, fill = layer)) + 
  scale_fill_gradientn(colors=rev(terrain.colors(10)),  breaks=seq(0, 1, by=0.1), labels=seq(0, 1, by=0.1), limits=c(0, 1)) +
  mapWorld +
  coord_quickmap(xlim = c(-15, 170), ylim = c(-10, 71)) +
  labs(title = "MaxEnt model for Vespa velutina",
       subtitle = "5-fold cross-validation with current climate")  +
  theme(panel.background = element_rect("white"),
        panel.grid = element_blank())

dev.off()
