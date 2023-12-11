

library(raster)
library(sf)
shp <- st_read(dsn = "C:/Old__CCAFSToolkit/NGhana/schema_shapefile.shp")
S1 = stack("C:/Users/youedraogo/Downloads/spam2017V2r1_SSA_A_MAIZ_A.tif")
plot(shp)
LZ_names=paste0("LZ",sprintf("%03d", 1:length(shp$CellID)))
Mean <- as.data.frame(t(raster::extract(S1, shp,sum,na.rm = T )))

colnames(Mean)= as.vector(LZ_names)
length(colnames(Mean))


nb=length(shp$CellID)

colnames(Mean)[1:nb]=as.vector(shp$CellID)

M = t(Mean)
M = as.data.frame(M)
colnames(M) = "vA"
M$CellID = shp$CellID
#M$vA = M$vA/5000

dd <- merge(M, shp, by = "CellID", all = TRUE)

# Convert the 'File' dataframe to a spatial dataframe using 'st_as_sf'
dd_sf <- st_as_sf(dd)

# Set tmap mode to "plot"
library(tmap)
tmap_mode("plot")




p2 <- tm_shape(dd_sf) +
  tm_polygons(col = "vA", title = "Mask ",style = "quantile",palette= "Greens") +
  tm_layout(legend.outside = F)
p2

p2 <- tm_shape(dd_sf) +
  tm_polygons(col = "vA", title = "Mask",  palette = "Greens") +
  tm_layout(legend.outside = FALSE)
p2

