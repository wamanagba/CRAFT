
library(raster)
library(sf)
library(sp)
library(raster)
library(rgeos)
library(tmap)
source_folder <- "C:/CCAFSToolkit/CRAFT_Schema/"
source_folder <- "C:/Old__CCAFSToolkit/BF/CRAFT_Schema_ByR/"






Grid_shape = function(source_folder){
  fileloc <- paste0(source_folder)
  projs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  ## define cellid for the whole world in 5 arc-min
  rwrd <- raster(res = 1/12) 
  rwrd[] <- 1:ncell(rwrd)
  
  start_time <- Sys.time()
  dt = data.frame()
  i=1
  j=1
  setwd(paste0(fileloc, "/Level", 1))
  ## read the shp file
  tmp <- shapefile(paste0("Shape/",list.files(path = "Shape/", pattern = "*.shp$")))
  #plot(tmp)
  ## spatial projection
  tmp <- spTransform(tmp, CRSobj = projs)
  ## attribute of Level names, here by the define of sname, we can select shp in just one district
  sname <- paste0("Level",1:i,"Name")
  ## return name of the attribute
  rname <- tmp@data[j,sname[i]]
  ## select
  tmp2 <- tmp[tmp@data[,sname[i]] == rname,]
  ## build a raster for this polygon, and give cell ids 
  rlay <- raster(xmn = floor(xmin(tmp2)), xmx = ceiling(xmax(tmp2)), ymn = floor(ymin(tmp2)), ymx = ceiling(ymax(tmp2)), res = 1/12)
  rlay <- crop(rwrd,rlay)
  ## build fishnet
  fishnet <- rasterToPolygons(rlay)
  #plot(fishnet)
  fishnet <- spTransform(fishnet, CRSobj = proj4string(tmp2))
  fish <- intersect(tmp2,fishnet)
  fishfine <- fishnet[fishnet@data$layer %in% fish@data$layer,]
  #plot(fish)
  order_fish <- order(fish@data$layer)
  fish <- fish[order_fish, ]
  # Order fishfine according to layer
  order_fishfine <- order(fishfine@data$layer)
  fishfine <- fishfine[order_fishfine, ]
  plot(fishfine)
  fish_sf <- st_as_sf(fish, coords = c("lon", "lat"), crs = st_crs(projs))
  columns_to_remove = c("Level1Name", "Id", "ObjectID")
  fish_sf <- fish_sf[, !(names(fish_sf) %in% columns_to_remove)]
  names(fish_sf)[names(fish_sf) == "layer"] <- "CellID"
  # Define the path and file name for your shapefile
  output_shapefile <-paste0("C:/CCAFSToolkit/Grid_Shape/",rname,"/Shape.shp")
  output_shapefil <-paste0("C:/CCAFSToolkit/Grid_Shape/",rname)
  dir.create(output_shapefil, recursive = TRUE,showWarnings = FALSE)
  st_write(fish_sf, output_shapefile,append=FALSE)
  return(rname)
}
#rname=Grid_shape(source_folder)

SpamData = function(rname){
  shp <- st_read(dsn = paste0("C:/CCAFSToolkit/Grid_Shape/",rname,"/Shape.shp"))
  S1 = stack("C:/Users/youedraogo/Downloads/spam2017V2r1_SSA_A_SORG_A.tif")
  
  LZ_names=paste0("LZ",sprintf("%03d", 1:length(shp$CellID)))
  Sum <- as.data.frame(t(raster::extract(S1, shp,sum,na.rm = T )))
  nb =length(Sum)
  colnames(Sum)= as.vector(LZ_names)
  colnames(Sum)[1:nb]=as.vector(shp$CellID)
  M = t(Sum)
  M = as.data.frame(M)
  colnames(M) = "PhysicalA"
  M$CellID = shp$CellID
  M$Percentage = M$PhysicalA/10000
  
  dd <- merge(M, shp, by = "CellID", all = TRUE)
  return(dd)
}


#rname=Grid_shape(source_folder)
# Convert the 'File' dataframe to a spatial dataframe using 'st_as_sf'
dd_sf <- st_as_sf(SpamData(Grid_shape(source_folder)))


# Set tmap mode to "plot"

tmap_mode("plot")
dd_sf <- st_make_valid(dd_sf)

# Création du tracé avec l'objet spatial corrigé
p2 <- tm_shape(dd_sf) +
  tm_polygons(col = "PhysicalA", title = "Mask ", style = "quantile", palette = "Greens") +
  tm_layout(legend.outside = FALSE)

# Affichage du tracé
p2

n_intervals <- 10  # Définir le nombre d'intervalles souhaités
col_name <- "PhysicalA"  # Remplacez par le nom de votre colonne à colorier

# Créer la carte en utilisant le style "pretty"
p2 <- tm_shape(dd_sf) +
  tm_fill(col = col_name, style = "pretty", n = n_intervals, palette = "Greens", title = "Mask") +
  tm_layout(legend.outside = FALSE)

# Afficher la carte
p2


p2 <- tm_shape(dd_sf) +
  tm_polygons(col = "PhysicalA", title = "Mask ",style = "quantile",palette= "Greens") +
  tm_layout(legend.outside = F)
p2

p2 <- tm_shape(dd_sf) +
  tm_polygons(col = "Percentage", title = "Mask",n=5,  palette = "Greens") +
  tm_layout(legend.outside = FALSE)
p2


















DownloadShape = function(ISO){
  # Define the URL of the ZIP file you want to download
  url <- paste0("https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_",ISO,"_shp.zip")
  
  # Define the path where you want to save the ZIP file
  destination <- paste0("C:/CCAFSToolkit/CRAFT_Schema/gadm41_",ISO,"_shp.zip")
  
  # Define the destination folder where you want to extract the files
  source_folder <- "C:/CCAFSToolkit/CRAFT_Schema/"
  dir.create(source_folder, recursive = TRUE,showWarnings = FALSE)
  
  # Download the ZIP file
  download.file(url, destination, mode = "wb")
  
  # Unzip the ZIP file into the destination folder
  unzip(destination, exdir = source_folder)
  
  # Remove the downloaded ZIP file to clean up
  file.remove(destination)
  
  # Print a message to indicate that the download and extraction are complete
  cat("Download and extraction completed.\n")
  
}
