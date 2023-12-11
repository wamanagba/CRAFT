



source_folder <- "C:/CCAFSToolkit/CRAFT_Schema/"

# Record the start time
start_time <- Sys.time()
library(sp)
library(raster)
library(rgeos)
InputFile = function(source_folder){
  
  ## define the location of the file (absolute　path)
  
  fileloc <- paste0(source_folder)
  ##make projection of Shape file as defaulted projection in r
  projs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  ## define cellid for the whole world in 5 arc-min
  rwrd <- raster(res = 1/12) 
  ## equals to 
  ## rwrd <- raster(res = 5/60)
  rwrd[] <- 1:ncell(rwrd)
  
  dt = data.frame()
  i=1 
  j=1
  
  setwd(paste0(fileloc, "/Level", i))
  
  ## read the shp file
  tmp <- shapefile(paste0("Shape/",list.files(path = "Shape/", pattern = "*.shp$")))
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
  
  fishnet <- spTransform(fishnet, CRSobj = proj4string(tmp2))
  plot(fishnet)
  fish <- intersect(tmp2,fishnet)
  fishfine <- fishnet[fishnet@data$layer %in% fish@data$layer,]
  
  # Order fish according to layer
  order_fish <- order(fish@data$layer)
  fish <- fish[order_fish, ]
  # Order fishfine according to layer
  order_fishfine <- order(fishfine@data$layer)
  fishfine <- fishfine[order_fishfine, ]
  #plot(fish)
  #plot(fishfine)
  
  schema <- data.frame(
    ID = fishfine@data$layer,
    Latitude = coordinates(fishfine)[, 2],  # Latitude is the second column of the coordinates
    Longitude = coordinates(fishfine)[, 1]  # Longitude is the first column of the coordinates
  )
  #schema <- schema[order(schema[,1], decreasing = TRUE),]
  schema <- schema[order(schema[,1]),]
  
  ## save schema to a relevant file
  fileout <- paste0("C:/Old__CCAFSToolkit/CCAFSToolkit/InputFile/Input.csv") 
  dir.create("C:/CCAFSToolkit/InputFile/", recursive = TRUE,showWarnings = FALSE)
  write.table(schema, file = fileout, row.names = FALSE, col.names = TRUE,quote = FALSE,sep = ",")
  
}


InputFile(source_folder)
