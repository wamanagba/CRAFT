



rm(list = ls())
source_folder <- "C:/CCAFSToolkit/CRAFT_Schema/"

# Record the start time
start_time <- Sys.time()


library(sp)
library(raster)
library(rgeos)
fileloc <- paste0(source_folder)
##make projection of Shape file as defaulted projection in r
projs <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
## define cellid for the whole world in 5 arc-min
rwrd <- raster(res = 1/12) 

rwrd[] <- 1:ncell(rwrd)

#plot(rwrd)


# Record the start time
start_time <- Sys.time()
dt = data.frame()
i=1
j=1

setwd(paste0(fileloc, "/Level", 2))

## read the shp file
tmp <- shapefile(paste0("Shape/",list.files(path = "Shape/", pattern = "*.shp$")))
plot(tmp)
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
plot(fishnet)
fishnet <- spTransform(fishnet, CRSobj = proj4string(tmp2))
plot(fish)
fish <- intersect(tmp2,fishnet)
fishfine <- fishnet[fishnet@data$layer %in% fish@data$layer,]



#Sys.sleep(2)

order_fish <- order(fish@data$layer)
fish <- fish[order_fish, ]
# Order fishfine according to layer
order_fishfine <- order(fishfine@data$layer)
fishfine <- fishfine[order_fishfine, ]

plot(fishfine)

library(sf)

fish_sf <- st_as_sf(fish, coords = c("lon", "lat"), crs = st_crs(projs))
columns_to_remove = c("Level1Name", "Id", "ObjectID")
fish_sf <- fish_sf[, !(names(fish_sf) %in% columns_to_remove)]
names(fish_sf)[names(fish_sf) == "layer"] <- "CellID"
# Define the path and file name for your shapefile
output_shapefile <- "C:/Old__CCAFSToolkit/NGhana2/schema_shapefile.shp"

# Save the data frame as a shapefile
st_write(fish_sf, output_shapefile,append=FALSE)
# Define the path to your shapefile
#shapefile_path <-"C:/CCAFSToolkit/InputFile/schema_shapefile.shp"

# Import the shapefile
fish_shapefile <- st_read(output_shapefile)
plot(fish_shapefile)

#Data = rio::import("C:/Users/youedraogo/Documents/Schema_Northern Ghana.txt")

library(data.table)

# Regular expression pattern to extract numbers from the filename
pattern <- "\\d+\\.WHT$"
dossier = "C:/Old__CCAFSToolkit/CCAFSToolkit/Schema_By_Python/Data"
# Create a list to store the data frames
dta <- list()

# Iterate through files in the directory
files <- list.files(dossier)
for (nom_fichier in files) {
  # Define the full file path
  chemin_du_fichier <- file.path(dossier, nom_fichier)
  
  # Read the WHT file as a data frame
  df <- fread(chemin_du_fichier, sep = ' ', skip = 7, header = FALSE, col.names = c("DATE","T2M","TMIN","TMAX","TDEW","RHUM","RAIN2","WIND","SRAD","RAIN"))
  # Calculate the sum of RAIN for each year
  df$RAIN <- as.numeric(df$RAIN)
  average_precipitation_per_year <- sum(df$RAIN, na.rm = TRUE)
  CELLID <- as.integer(substr(nom_fichier, 1, nchar(nom_fichier) - 4))
  
  # Create a data frame for each CELLID and average precipitation
  df <- data.frame(CellID = CELLID, average_prec = average_precipitation_per_year)
  
  # Add the data frame to the list
  dta <- c(dta, list(df))
}

# Concatenate the data frames in the list into one data frame
final_df <- do.call(rbind, dta)


shapefile <- fish_shapefile
DataShapefile <- as.data.frame(shapefile)

dd <- merge(DataShapefile, final_df, by = "CellID", all.x = TRUE)

# Convert the 'File' dataframe to a spatial dataframe using 'st_as_sf'
dd_sf <- st_as_sf(dd)

# Set tmap mode to "plot"
library(tmap)
tmap_mode("plot")




p2 <- tm_shape(dd_sf) +
  tm_polygons(col = "average_prec", title = "Precipitation ",
              style = "quantile") +
  tm_layout(legend.outside = F)
p2

















merged_data <- merge.data.frame(fish_shapefile, Data, by = "CellID")


# Convert the 'File' dataframe to a spatial dataframe using 'st_as_sf'
dd_sf <- st_as_sf(merged_data)

# Set tmap mode to "plot"
library(tmap)
tmap_mode("view")
CV <- tm_shape(dd_sf) +
  tm_polygons(col = "Mean", title = "Mean", palette= "YlOrBr") 

CV
tmap_save(CV, filename = file.path("C:/CCAFSToolkit/InputFile/UpperEast/", "Climatologyhtl.html"))

  
  