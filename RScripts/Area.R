


rm(list = ls())
source("C:/Users/youedraogo/Desktop/CRAFT/RScripts/All.R")
source_folder <- "C:/CCAFSToolkit/CRAFT_Schema/"

Area = function(source_folder){
  
  ## define the location of the file (absolute　path)
  fileloc <- paste0(source_folder)
  ##make projection of Shape file as defaulted projection in r
  projs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  ## define cellid for the whole world in 5 arc-min
  rwrd <- raster(res = 1/12) 
  ## equals to 
  ## rwrd <- raster(res = 5/60)
  rwrd[] <- 1:ncell(rwrd)
  #plot(rwrd)
  dt = data.frame()
  i=1
  setwd(paste0(fileloc, "/Level", i))
    
  # Check if the folder have a shape file or if the folder shape exist
  if(CheckFolder(source_folder,i)==1){
      ## read the shp file
      tmp <- shapefile(paste0("Shape/",list.files(path = "Shape/", pattern = "*.shp$")))
      if(i==3){ tmp@data$Level3Name = paste(tmp@data$Level2Name, tmp@data$Level3Name, sep = "_")}
      
      ## spatial projection
      tmp <- spTransform(tmp, CRSobj = projs)
      ## attribute of Level names, here by the define of sname, we can select shp in just one district
      sname <- paste0("Level",1:i,"Name")
      ## make loops for all the polygon in for Level 'i'
      #j=1
      cpt=0
      j=1
      rname <- tmp@data[j,sname[i]]
      ## return of name ahead of the selected Level, this is used to name the result and write it out
      tname <- tmp@data[j, sname[1:(i - 1)]]
      tname <- paste(tname,collapse = "_")
      if (i == 1) { tname <- NULL}
        ## select
      tmp2 <- tmp[tmp@data[,sname[i]] == rname,]
      ## build a raster for this polygon, and give cell ids 
      rlay <- raster(xmn = floor(xmin(tmp2)), xmx = ceiling(xmax(tmp2)), ymn = floor(ymin(tmp2)), ymx = ceiling(ymax(tmp2)), res = 1/12)
      rlay <- crop(rwrd,rlay)
        
      ## build fishnet
      fishnet <- rasterToPolygons(rlay)
      fishnet <- spTransform(fishnet, CRSobj = proj4string(tmp2))
      #plot(fishnet)
      fish <- intersect(tmp2,fishnet)
      fishfine <- fishnet[fishnet@data$layer %in% fish@data$layer,]
      plot(fish)
      # Order fish according to layer
      order_fish <- order(fish@data$layer)
      fish <- fish[order_fish, ]
      # Order fishfine according to layer
      order_fishfine <- order(fishfine@data$layer)
      fishfine <- fishfine[order_fishfine, ]
        
      centroids <- st_centroid(st_as_sf(fishfine))
        
      # Retrieving centroid coordinates
      centroids_coords <- st_coordinates(centroids)

      Areea_fishfine = (area(fishfine) / 1e6)
      # Adding centroid coordinates to the schema dataframe
      LAT <- centroids_coords[, 2]  # Latitude
      LON <- centroids_coords[, 1]  # Longitude
        
      schema <- data.frame(CELLID = fishfine@data$layer,LAT = LAT, LONG= LON, elevation = -99.00,Area=Areea_fishfine,LandeUse=1)
      schema <- schema[order(schema[,1], decreasing = TRUE),]
      ## save schema to a relevant file
      if(is.null(tname)){fileout <- paste0(getwd(), "/Schema/", "5m", tname,"_", rname, ".txt") }else{fileout <- paste0(getwd(), "/Schema/", "5m_", tname,"_", rname, ".txt")}
      dir.create(paste0(getwd(), "/Schema/"), recursive = TRUE,showWarnings = FALSE)
      write.table(schema, file = fileout, row.names = FALSE, col.names = TRUE,quote = FALSE,sep = "\t")

      output_file <- "C:/CCAFSToolkit/Schema_By_R/SchemaGenerationDetails.csv"
      # Write the 'data' dataframe to the CSV file
      write.table(dt, file = output_file, sep = ",", row.names = FALSE, col.names = FALSE)
      
    }else{
      cat("No shape file on the Level",i ,".   Check in the folder Shape \n" )
    }
    cat("Schema generation in-progress: Schema generation Successful for Level ",i, "\n")
  
}

start_time <- Sys.time()
Area(source_folder)
# Record the end time
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
