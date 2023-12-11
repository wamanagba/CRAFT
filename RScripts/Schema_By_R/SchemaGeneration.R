


   ### Schema Generation & World Grid & Input File Generation    ###



##  1-    Schema Generation

##################################################################################################

#This script downloads a shapefile from the GADM online site in the form of a ZIP file, 
#extracts its contents into the "CRAFT_Schema" folder, and deletes the downloaded ZIP file upon completion. 
#Users must ensure they have an active internet connection. 
#If the destination folder "C:/CCAFSToolkit/CRAFT_Schema/" does not exist, 
#it will be created automatically before running the script. 
#This script takes the "ISO" as a parameter, which represents the ISO code of 
#the country from which you want to extract the SCHEMA.

rm(list = ls())
#ISO = "BFA"

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






##########################################################################################
"This script organizes Shapefile files into separate folders based on their 'Level.' 
It scans a source folder containing Shapefiles, and for each Shapefile, 
it determines the 'Level' based on the last character of the filename. 
It then moves each Shapefile to a corresponding folder named 'Level1,' 'Level2,' or 'Level3' 
depending on its 'Level.' The script considers Shapefiles with different extensions
such as '.shp,' '.shx,' and '.dbf' and ensures they are placed in the appropriate folders."

##################################################################################


CreateLevel= function(source_folder){
  # Path to the source folder where your Shapefiles are located
  Country= "Shape"
  
  # Create the three folders Level1, Level2, and Level3
  for (i in 1:6) {
    dir.create(file.path(source_folder, paste0("Level", i)), showWarnings = FALSE)
  }
  
  # List of file extensions to consider (shp, shx, dbf,cpg,prj,)
  extensions <- c("shp", "shx", "dbf","cpg","prj","csv")
  
  # Loop through the extensions
  for (ext in extensions) {
    # List of files with the current extension in the source folder
    files <- list.files(path = source_folder, pattern = paste0("\\.", ext, "$"), full.names = TRUE)
    
    # Loop through the files and move them to the appropriate folders
    for (file in files) {
      # Get the file name without extension
      file_name <- tools::file_path_sans_ext(basename(file))
      
      # Get the last character of the file name (0, 1, or 2)
      last_character <- substr(file_name, nchar(file_name), nchar(file_name))
      
      # Calculate the destination folder based on the last character
      last_character <- as.numeric(last_character) + 1
      destination_folder <- file.path(source_folder, paste0("Level", last_character))
      
      # Move the file to the destination folder with the current extension
      file.rename(file, file.path(destination_folder, paste0(file_name, ".", ext)))
    }
  }
  
  
}

  



##########################################

ShapeProcessing = function(source_folder){
  library(utils)
  library(sf)
  Country= "Shape"
  setwd(source_folder)
  # Create the three folders Level1, Level2, and Level3
  CRAFT_folder <- paste0(source_folder)
  
  # Check if the CRAFT directory exists, if not, create it
  if (!file.exists(CRAFT_folder)) {
    dir.create(CRAFT_folder, recursive = TRUE)
  }
  
  for (i in 1:3) {
    dir.create(paste0(CRAFT_folder, paste0("Level", i)), showWarnings = FALSE)
  }
  
# All the variable that we don't need
columns_to_remove = c("GID_1", "GID_0", "VARNAME_1", "NL_NAME_1", "TYPE_1", "ENGTYPE_1", 
                        "CC_1", "HASC_1", "ISO_1","GID_0","ID_0", "ISO",  "OBJECTID_1", "ISO3",     
                        "NAME_ENGLI", "NAME_ISO",   "NAME_FAO",   "NAME_LOCAL", "NAME_OBSOL",
                        "NAME_VARIA", "NAME_NONLA", "NAME_FRENC", "NAME_SPANI", "NAME_RUSSI",
                        "NAME_ARABI", "NAME_CHINE", "WASPARTOF" , "CONTAINS"   ,"SOVEREIGN" ,
                        "ISO2"       ,"WWW"    ,    "FIPS"   ,    "ISON"    ,   "VALIDFR"   ,
                        "VALIDTO" ,   "POP2000"  ,  "SQKM"  ,     "POPSQKM"  ,  "UNREGION1" ,
                        "UNREGION2" , "DEVELOPING", "CIS"   ,     "Transition", "OECD"      ,
                        "WBREGION",   "WBINCOME" ,  "WBDEBT" ,    "WBOTHER"    ,"CEEAC"     ,
                        "CEMAC" ,     "CEPLG" ,     "COMESA" ,    "EAC"      ,  "ECOWAS"    ,
                        "IGAD"  ,     "IOC"   ,     "MRU"    ,    "SACU"       ,"UEMOA"     ,
                        "UMA"  ,      "PALOP"  ,    "PARTA" ,     "CACM"   ,    "EurAsEC"   ,
                        "Agadir" ,    "SAARC" ,     "ASEAN" ,     "NAFTA"  ,    "GCC"       ,
                        "CSN"  ,      "CARICOM",    "EU"     ,    "CAN"     ,   "ACP"       ,
                        "Landlocked" ,"AOSIS" ,     "SIDS"    ,   "Islands"  ,  "LDC",
                        "GID_2", "GID_0",   "GID_1", "NL_NAME_1","CC_2",   "HASC_2", "VARNAME_2",
                        "NL_NAME_2",   "TYPE_2", "ENGTYPE_2","ID_1","ID_2")
  
  
  # Extract the Level 1
  directory_path=file.path(source_folder, paste0("Level", 1))
  # List all Shapefiles in the directory
  Shapefile_list <- list.files(directory_path, pattern = "\\.shp$", full.names = TRUE)
  
  # Check if any Shapefiles were found
  if (length(Shapefile_list) > 0) {
    # Read the first Shapefile in the list
    Shape <- st_read(Shapefile_list[1])
    Level1 =paste0(Country)
    names(Shape)[names(Shape) == "COUNTRY"] <- "Level1Name"
    names(Shape)[names(Shape) == "NAME_0"] <- "Level1Name"
    
    # delete column
    Shape <- Shape[, !(names(Shape) %in% columns_to_remove)]
    Shape$ObjectID = 1
    
    Shape$Level1Name <- iconv(Shape$Level1Name, to = "ASCII//TRANSLIT")
    
    dir.create(paste0(CRAFT_folder, paste0("Level", 1),"/Shape"), showWarnings = FALSE)
    st_write(Shape,paste0(CRAFT_folder,"Level1/Shape/",Level1,"_01.shp"), delete_layer = T)
    
  } else {cat("No Shapefiles found in the Level 1 folder \n")}
  

  
  
  # Extract the Level 2
  directory_path=file.path(source_folder, paste0("Level", 2))
  # List all Shapefiles in the directory
  Shapefile_list <- list.files(directory_path, pattern = "\\.shp$", full.names = TRUE)
  # Check if any Shapefiles were found
  if (length(Shapefile_list) > 0) {
    # Read the first Shapefile in the list
    Shape <- st_read(Shapefile_list[1])
    # delete some column
    Shape <- Shape[, !(names(Shape) %in% columns_to_remove)]
    # Renommer la colonne 'Level3Name' en 'Level2Name'
    names(Shape)[names(Shape) == "NAME_1"] <- "Level2Name"
    names(Shape)[names(Shape) == "COUNTRY"] <- "Level1Name"
    names(Shape)[names(Shape) == "NAME_0"] <- "Level1Name"
    
    Shape$ObjectID = 0
    
    Shape$Level1Name <- iconv(Shape$Level1Name, to = "ASCII//TRANSLIT")
    Shape$Level2Name <- iconv(Shape$Level2Name, to = "ASCII//TRANSLIT")
    
    dir.create(paste0(CRAFT_folder, paste0("Level", 2),"/Shape"), showWarnings = FALSE)
    st_write(Shape,paste0(CRAFT_folder,"Level2/Shape/",Level1,"_02.shp"), delete_layer = T)
    
    } else {cat("No Shapefiles found in the Level 2 folder.\n")}
    
  
  
  
  # Extract the Level 3
  directory_path=file.path(source_folder, paste0("Level", 3))
  # List all Shapefiles in the directory
  Shapefile_list <- list.files(directory_path, pattern = "\\.shp$", full.names = TRUE)
  
  # Check if any Shapefiles were found
  if (length(Shapefile_list) > 0) {
    # Read the first Shapefile in the list
    Shape <- st_read(Shapefile_list[1])
    # delete some column
    Shape <- Shape[, !(names(Shape) %in% columns_to_remove)]
    # Rename the columns
    names(Shape)[names(Shape) == "NAME_2"] <- "Level3Name"
    names(Shape)[names(Shape) == "NAME_1"] <- "Level2Name"
    names(Shape)[names(Shape) == "COUNTRY"] <- "Level1Name"
    names(Shape)[names(Shape) == "NAME_0"] <- "Level1Name"
    Shape$ObjectID = 0
    
    Shape$Level1Name <- iconv(Shape$Level1Name, to = "ASCII//TRANSLIT")
    Shape$Level2Name <- iconv(Shape$Level2Name, to = "ASCII//TRANSLIT")
    Shape$Level3Name <- iconv(Shape$Level3Name, to = "ASCII//TRANSLIT")
    
    dir.create(paste0(CRAFT_folder, paste0("Level", 3),"/Shape"), showWarnings = FALSE)
    st_write(Shape,paste0(CRAFT_folder,"Level3/Shape/",Level1,"_03.shp"), delete_layer = T)
    
    } else {cat("No Shapefiles found in the Level 3 folder. \n")}
    
}



####################################  Delete the firts Shape file in the folder
DeleteOldFolder = function(source_folder){
  # List of directory names to delete
  directories_to_delete <- c("Level4","Level5","Level6")
  
  # Loop through the directories and delete them
  for (dir_name in directories_to_delete) {
    directory_path <- file.path(source_folder, dir_name)
    
    if (file.exists(directory_path)) {
      # Delete the directory and its contents recursively
      unlink(directory_path, recursive = TRUE)
      cat("Deleted directory:", directory_path, "\n")
    } else {
      cat("Directory does not exist:", directory_path, "\n")
    }
  }
}









DeleteOldFiles <- function(){
  

# Specify the folders Level1, Level2, and Level3
folders <- c("Level1", "Level2", "Level3")

# File extensions to delete
extensions <- c("shp", "shx", "dbf", "cpg", "prj", "csv")

# Loop through the folders
for (folder in folders) {
  # Create the full path of the folder
  folder_path <- file.path(source_folder, folder)
  
  # List all files in the folder
  files <- list.files(path = folder_path, full.names = TRUE)
  
  # Filter files with the specified extensions
  files_to_delete <- files[grep(paste(extensions, collapse = "|"), files, ignore.case = TRUE)]
  
  # Delete the files
  if (length(files_to_delete) > 0) {
    file.remove(files_to_delete)
    cat("Files deleted in folder", folder, ":", length(files_to_delete), "file(s)\n")
  } else {
    cat("No files to delete in folder", folder, "\n")
  }
}
}

################################################################################

CheckFolder <- function(source_folder,i){
  fileloc <- paste0(source_folder)
  # Extract the Level 2
  directory_path=file.path(fileloc, paste0("Level", i,"/Shape/"))
  # List all Shapefiles in the directory
  Shapefile_list <- list.files(directory_path, pattern = "\\.shp$", full.names = TRUE)
  # Check if any Shapefiles were found
  if (length(Shapefile_list) > 0) {check <- 1} else {check <- 0}
  return(check)
}
###################################################################


####################################################

Message_start <- function(i){
  # Specify the folder name for messages
  messages_folder <- paste0(source_folder,"messages/")
  
  # Check if the messages folder exists, if not, create it
  if (!dir.exists(messages_folder)) {
    dir.create(messages_folder)
  }
  
  # Message at the beginning of the script
  start_message <- paste0("Schema generation in-progress: Attempting to generate schema for Level ", i, "\n")
  start_file <- paste0(messages_folder, "start_Level",i,".txt")
  writeLines(start_message, start_file)
}


Message_End <- function(i){
  # Specify the folder name for messages
  messages_folder <- paste0(source_folder,"messages/")
  
  # Check if the messages folder exists, if not, create it
  if (!dir.exists(messages_folder)) {
    dir.create(messages_folder)
  }
  
  
  start_message <- paste0("Schema generation in-progress: Schema generation Successful for Level ",i)
  start_file <- paste0(messages_folder, "End_Level",i,".txt")
  writeLines(start_message, start_file)
}


################################################
dt = data.frame()
File = function(dt,i,cpt){
  data <- data.frame(
    x = c(paste0("Level",i)),
    y = c(cpt)
  )
  dt =rbind(dt,data)
}
#############################################################################
##################  2 - World Grid function #################################
#############################################################################



WorldGrid = function(fileloc){
  j=1
  i=1
  ##make projection of Shape file as defaulted projection in r
  projs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  ## define cellid for the whole world in 5 arc-min
  rwrd <- raster(res = 1/12) 
  ## equals to 
  ## rwrd <- raster(res = 5/60)
  rwrd[] <- 1:ncell(rwrd)
  
  ## set working directory here for the level 1
  setwd(paste0(fileloc, "/Level", i))
  ## read the shp file
  tmp <- shapefile(paste0("shape/",list.files(path = "Shape/", pattern = "*.shp$")))
  ## spatial projection
  tmp <- spTransform(tmp, CRSobj = projs)
  ## attribute of level names, here by the define of sname, we can select shp in just one district
  sname <- paste0("Level",1:i,"Name")
  
  ## return name of the attribute
  rname <- tmp@data[j,sname[i]]
  ## return of name ahead of the selected level, this is used to name the result and write it out
  tname <- tmp@data[j, sname[1:(i - 1)]]
  tname <- paste(tname,collapse = "_")
  if (i == 1) { tname <- NULL}
  ## select
  tmp2 <- tmp[tmp@data[,sname[i]] == rname,]
  ## build a raster for this polygon, and give cell ids 
  rlay <- raster(xmn = (xmin(tmp2))-1/12, xmx = (xmax(tmp2))+1/12, ymn = (ymin(tmp2))-1/12, ymx = (ymax(tmp2))+1/12, res = 1/12)
  rlay <- crop(rwrd,rlay)
  
  # Convert the raster layer to polygons
  fishnet_polygons <- rasterToPolygons(rlay)
  # Assuming 'fishnet' is your polygon layer
  centroid_coords <- coordinates(fishnet_polygons)
  # Add the cell ID information to the polygons
  names(fishnet_polygons)[names(fishnet_polygons) == "layer"] <- "CELLID"
  fishnet_polygons$id <- 0
  
  # Add separate columns for LAT and LON
  fishnet_polygons$LAT <- centroid_coords[, 2]  # Assuming the latitude is in the second column of centroid_coords
  fishnet_polygons$LON <- centroid_coords[, 1]  # Assuming the longitude is in the first column of centroid_coords
  
  # Define the file path where you want to save the shapefile
  shapefile_path <- "C:/CCAFSToolkit/WorldGridData/"
  # Save the polygons as a shapefile
  dir.create(paste0(shapefile_path), recursive = TRUE,showWarnings = FALSE)
  writeOGR(fishnet_polygons, dsn = shapefile_path, layer = "WorldGrid", driver = "ESRI Shapefile",overwrite = TRUE)
  
  # close all files
  closeAllConnections()
  
}


##############################################################################
###########                                    ###############################
########### Input files for weather generation ###############################
###########                                    ###############################
##############################################################################

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
    CellID = fishfine@data$layer,
    Latitude = coordinates(fishfine)[, 2],  # Latitude is the second column of the coordinates
    Longitude = coordinates(fishfine)[, 1]  # Longitude is the first column of the coordinates
  )
  #schema <- schema[order(schema[,1], decreasing = TRUE),]
  schema <- schema[order(schema[,1]),]
  colnames(schema)[1] = "ID"
  ## save schema to a relevant file
  fileout <- paste0("C:/CCAFSToolkit/Schema_By_Python/Input.csv") 
  dir.create("C:/CCAFSToolkit/Schema_By_Python/", recursive = TRUE,showWarnings = FALSE)
  write.table(schema, file = fileout, row.names = FALSE, col.names = TRUE,quote = FALSE,sep = ",")
  
}

################################  Schema creation ####################################
## the script here is used to generate schema for Shape file among different
## Levels
######################################################################################

Schema = function(source_folder){
  
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
  
  
  # Record the start time
  start_time <- Sys.time()
  dt = data.frame()
  #i=1 
  ## from now on, do schema generation Level by Level
  for (i in 1:3) { ## note 3 is the Levels of Shape file you have
    ## set working directory here for the Level 'i'
    
    cat("Schema generation in-progress: Attempting to generate schema for Level ", i,"\n")
    #Message_start(i)
    setwd(paste0(fileloc, "/Level", i))
    
    # Check if the folder have a shape file or if the folder shape exist
    CheckFolder(source_folder,i)
    if(CheckFolder(source_folder,i)==1){
      
      ## read the shp file
      tmp <- shapefile(paste0("Shape/",list.files(path = "Shape/", pattern = "*.shp$")))
      ## spatial projection
      tmp <- spTransform(tmp, CRSobj = projs)
      ## attribute of Level names, here by the define of sname, we can select shp in just one district
      sname <- paste0("Level",1:i,"Name")
      ## make loops for all the polygon in for Level 'i'
      #j=1
      cpt=0
      for (j in 1:nrow(tmp@data)) {
        ## return name of the attribute
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
        #plot(fishnet)
        #plot(tmp2,add=T)
        #Sys.sleep(2)
        fishnet <- spTransform(fishnet, CRSobj = proj4string(tmp2))
        #plot(fishnet)
        fish <- intersect(tmp2,fishnet)
        fishfine <- fishnet[fishnet@data$layer %in% fish@data$layer,]
        #plot(fish)
        #plot(fishfine)
        
        plot(fish)
        #Sys.sleep(2)
        ##### New code
        # Order fish according to layer
        order_fish <- order(fish@data$layer)
        fish <- fish[order_fish, ]
        # Order fishfine according to layer
        order_fishfine <- order(fishfine@data$layer)
        fishfine <- fishfine[order_fishfine, ]
        
        
        ## the share percent was calculated by the area of fish/fishfine,
        ## through the comparison with result from CRAFT, the difference might just comes from the round() function
        share <- round(area(fish)/area(fishfine)*100,2)
        schema <- data.frame(CELLID = fishfine@data$layer,SHAREPERCENT = share)
        schema <- schema[order(schema[,1], decreasing = TRUE),]
        ## save schema to a relevant file
        if(is.null(tname)){fileout <- paste0(getwd(), "/Schema/", "5m", tname,"_", rname, ".txt") }else{fileout <- paste0(getwd(), "/Schema/", "5m_", tname,"_", rname, ".txt")}
        dir.create(paste0(getwd(), "/Schema/"), recursive = TRUE,showWarnings = FALSE)
        write.table(schema, file = fileout, row.names = FALSE, col.names = TRUE,quote = FALSE,sep = "\t")
        cpt=cpt+1
      }
      dt = File(dt,i,cpt)
      output_file <- "C:/CCAFSToolkit/Schema_By_R/SchemaGenerationDetails.csv"
      # Write the 'data' dataframe to the CSV file
      write.table(dt, file = output_file, sep = ",", row.names = FALSE, col.names = FALSE)
      
    }else{
      cat("No shape file on the Level",i ,".   Check in the folder Shape \n" )
    }
    cat("Schema generation in-progress: Schema generation Successful for Level ",i, "\n")
    #Message_End(i)
  }
  
  
}

# Record the start time
start_time <- Sys.time()
# prepare the inputs
source_folder <- "C:/CCAFSToolkit/CRAFT_Schema/"
#ISO=rio::import("C:/CCAFSToolkit/Schema_By_R/ISO.csv")
#ISO=ISO[1,1]
#DownloadShape(ISO)

# packages requirements
library(sp)
library(raster)
library(rgeos)
library(rgdal)


# Run the functions
CreateLevel(source_folder)

ShapeProcessing(source_folder)

DeleteOldFolder(source_folder)

DeleteOldFiles()

Schema(source_folder)

worldGrid(source_folder)

InputFile(source_folder)


# Record the end time
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)


