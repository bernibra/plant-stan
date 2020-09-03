
#Download WSL data if not there. These are only bioclim variables at approx 1km of resolution
download_climatic <- function(res=0.5, meanlat, meanlong){
  
  bioclim.data <- getData(name = "worldclim",
                          var = "bio",
                          res = res,
                          lat=meanlat,
                          lon=meanlong,
                          path = "../data/raw/climate-data/")

}

#Convert Lat and Lon to N and E
transform_coordinates <- function(lon, lat, data_system = "+init=epsg:4326", map_system = "+proj=somerc +init=world:CH1903"){
  #Alternatively, you can use the following projection for Antoine's data: map_system <- "+init=epsg:21781"
  
  # Defining current coordinates
  cord.dec = SpatialPoints(cbind(lon, lat), proj4string=CRS(data_system))

  # Chaning coordinates
  cord.UTM <- spTransform(cord.dec, map_system)

  return(list(x=cord.UTM@coords[,1], y=cord.UTM@coords[,2]))

}

#Check if the site is in the Rechalp area
is.rechalp <- function(folder, lon, lat){

  # Change to Northing Easting
  NorthEeast <- transform_coordinates(lon=lon,lat=lat)
  
  # Import raster example
  file <- paste(folder, "biovars/yearly/bio1_tmean1981_Rechalp_ngb5_mwconic.tif", sep="")
  bioclim.data <- stack(file)
  crs(bioclim.data)<-"+proj=somerc +init=world:CH1903"

  # Check if the coordinates are there
  return(!is.na(extract(bioclim.data, data.frame(NorthEeast))))
}

little.check.rechalp <- function(folder){

  # Generate grid
  lat_ <- seq(46, 46.8, length.out = 100)
  lon_<- seq(6.5, 8, length.out = 100)

  # Import raster example
  file <- paste(folder, "biovars/yearly/bio1_tmean1981_Rechalp_ngb5_mwconic.tif", sep="")
  bioclim.data <- stack(file)
  crs(bioclim.data)<-"+proj=somerc +init=world:CH1903"
  
  # Check that things would work as expected
  dat <- expand.grid(lat=lat_, lon=lon_)
  NorthEeast <- transform_coordinates(lon=dat$lon, lat=dat$lat)
  dat$value <- as.vector(!is.na(extract(bioclim.data, data.frame(NorthEeast))))
  plot(dat$lon, dat$lat, col=as.integer(dat$value)+2, pch=20)
  
}

start.rechalp.file <- function(filenames, rows, columns){
  
  # Start first all R objects to optimize the code
  data <- data.frame(matrix(NA, nrow=length(rows), ncol=length(columns)))
  colnames(data) <- columns
  rownames(data) <- rows
  for(filename in filenames){
    saveRDS(data, file = filename)
  }

}

write.value.data <- function(filenames, row, column, z){
  
  #Write a value in the R object
  for(x in 1:length(filenames)){
    data <- readRDS(file = filenames[x])
    data[row, column] <- z[x]
    saveRDS(data, file = filenames[x])
  }
  
}

extract.wsl <- function(folder, data, start_files=TRUE){
  
  #Two folders
  folder_A <- paste(folder, "wsl/", sep="")
  folder_B <- paste(folder, "wc0.5/", sep="")
  
  #Translate all coordinates
  Coordinates_chelsa_and_worldclim <- transform_coordinates(lat=data$Lat_WGS84, lon=data$Long_WGS84, map_system = "+proj=longlat +datum=WGS84 +no_defs")
  
  #Download additional data
  download_climatic(res=0.5, meanlat = mean(data$Lat_WGS84), meanlong = mean(data$Long_WGS84))
  
  #Find files for all climatic variables
  files_A <- list.files(folder_A, pattern = "tif$", full.names = FALSE)
  files_B <- paste("worldclim", list.files(folder_B, pattern = "bil$", full.names = FALSE), sep="_")
  if (length(files_A)==0) {
    stop("You didn't download the WSL data. You can do this using the 'wget' command in the terminal to download the data from 'https://envidatrepo.wsl.ch/uploads/chelsa/'")
  }
  
  #Find out all variables that we need
  files_ <- c(files_A, files_B)
  info_files <- as.data.frame(do.call(rbind, strsplit(files_, split = "_")))
  info_files$files <- files_
  info_files$bioclim <- NA
  chelsa <- as.character(info_files$V1)=="CHELSA"
  info_files$bioclim[chelsa] <- paste("bio", as.numeric(gsub("\\.tif","", info_files$V3[chelsa])), sep="")
  info_files$bioclim[!chelsa] <- as.character(info_files$V2[!chelsa])
  info_files$files[chelsa] <- paste(folder_A, info_files$files[chelsa], sep="")
  info_files$files[!chelsa] <- paste(folder_B, gsub("worldclim_", "", info_files$files[!chelsa]), sep = "")

  #Deine output dimensions  
  rows <- unique(info_files$V1)
  columns <- unique(info_files$bioclim)
  
  #Start all files to store the data
  if(start_files){
    start.rechalp.file(data$filename, rows, columns) 
  }
  
  #Loop over all variables and years to extract the information
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  for(i in 1:length(info_files$files)){
    bioclim.data <- stack(info_files$files[i])
    crs(bioclim.data)<-"+proj=longlat +datum=WGS84 +no_defs"
    z <- extract(bioclim.data, data.frame(Coordinates_chelsa_and_worldclim))
    setTxtProgressBar(pb, i)
    write.value.data(filenames = data$filename, row = as.character(info_files$V1[i]), column = as.character(info_files$bioclim[i]), z = as.numeric(z))
  }
  close(pb) 
  
}

extract.rechalp <- function(folder, data, additional_folder, start_files=TRUE){
  
  #Translate all coordinates
  NorthEeast <- transform_coordinates(lat=data$Lat_WGS84, lon=data$Long_WGS84)
  
  #Find out all variables that we need
  folder <- paste(folder, "biovars/yearly/", sep="")
  files_ <- list.files(folder, pattern = "tif$", full.names = FALSE)
  files <- paste(folder, files_, sep="")
  info_files <- as.data.frame(do.call(rbind, strsplit(files_, split = "_")))
  info_files$files <- files
  info_files$bioclim <- info_files$V1
  info_files$variable <- gsub("[[:digit:]]", "", info_files$V2)
  info_files$year <- gsub("[^[:digit:]]", "", info_files$V2)
  
  #Defining output dimensions
  rows <- unique(info_files$year)
  columns <- unique(info_files$bioclim)
  
  #Addin wordlcim and Chelsa to the output data
  rows <- c(as.character(rows), "CHELSA", "worldclim")
  
  #Start all files to store the data if start_files=TRUE
  if(start_files){
    start.rechalp.file(data$filename, rows, columns) 
  }
  
  #Loop over all variables and years to extract the information
  pb <- txtProgressBar(min = 0, max = length(files), style = 3)
  for(i in 1:length(info_files$files)){
    bioclim.data <- stack(info_files$files[i])
    crs(bioclim.data)<-"+proj=somerc +init=world:CH1903"
    z <- extract(bioclim.data, data.frame(NorthEeast))
    setTxtProgressBar(pb, i)
    write.value.data(filenames = data$filename, row = as.character(info_files$year[i]), column = as.character(info_files$bioclim[i]), z = as.numeric(z))
  }
  close(pb)
  
  #Loop over Chelsa and Worldclim
  extract.wsl(folder = additional_folder, data = data, start_files = FALSE)
  
}

write.readme.rechalp <- function(folder, data){
  
  for (i in 1:nrow(data)){
    
    # Writing things down: Readme, csv, zip...
    filename <- paste(processed_folder, "sites/", as.character(data$Codes[i]), "/", as.character(data$Codes[i]), sep = "")
    file.copy("../data/raw/readme-templates/rechalp.md", paste(filename, ".md", sep = ""), overwrite = T)
    mat <- readRDS(file = paste(filename, ".Rds", sep = ""))
    write.table(mat, file = paste(filename, ".csv", sep = ""), col.names=NA, quote = F, sep = ",")
    
    # For some reason this doesn't work very well in mac... use "zip -r outputfile.zip folder" in the command line instead
    # foldername <- paste(processed_folder, "sites/", as.character(data$Codes[i]), sep = "")
    # zip(paste(foldername, ".zip", sep=""), foldername)
    
  }
}

write.readme.wsl <- function(folder, data){
  
  for (i in 1:nrow(data)){

    # Writing things down: Readme, csv, zip...
    filename <- paste(processed_folder, "sites/", as.character(data$Codes[i]), "/", as.character(data$Codes[i]), sep = "")
    file.copy("../data/raw/readme-templates/other.md", paste(filename, ".md", sep = ""), overwrite = T)
    mat <- readRDS(file = paste(filename, ".Rds", sep = ""))
    write.table(mat, file = paste(filename, ".csv", sep = ""), col.names=NA, quote = F, sep = ",")

    # For some reason this doesn't work very well in mac... use "zip -r outputfile.zip folder" in the command line instead
    # foldername <- paste(processed_folder, "sites/", as.character(data$Codes[i]), sep = "")
    # zip(paste(foldername, ".zip", sep=""), foldername)

    }
}
