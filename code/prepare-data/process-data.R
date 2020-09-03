## Load library
# library(readtext)
# library(stringr)
# library(dplyr)
# library(rgdal)
# library(dismo)
library(raster)
# library(FD)
# library(FactoMineR);library(factoextra)

## Local library
source("useful-tools.R")

plot_cordinates <- function(x){
  
  # Load the data to use for our base map
  data(wrld_simpl)
  # Plot base map
  extra <- 0
  plot(wrld_simpl, 
       xlim = c(min(x$long)-extra, max(x$long)+extra),
       ylim = c(min(x$lat)-extra, max(x$lat)+extra),
       axes = TRUE, 
       col = "grey95")
  
  # Add original observations
  points(x$long, x$lat, col = "olivedrab", pch = 20, cex = 0.75)
  box()
}

prepare_data_SDM <- function(){
  ## Load data
  dat.weighted <- read.csv("../../data/raw/distribution/grassland.abundance.csv", sep=";")

  ## Extract coordinates
  x.weighted <- dat.weighted[, "X"]
  y.weighted <- dat.weighted[, "Y"]
  elevation <- dat.weighted[, "Elevation"]

  ## Changing the coordinate system
  dat <- data.frame(Easting=x.weighted, Northing=y.weighted)
  coords <- cbind(Easting=x.weighted, Northing=y.weighted)
  swissgrid = "+proj=somerc +init=world:CH1903"
  
  dat_SP <- SpatialPointsDataFrame(coords,
                                   data = dat,
                                   proj4string = CRS(swissgrid))
  
  ## Extract coocurrence data
  dat.weighted <- dat.weighted[,7:ncol(dat.weighted)]
  
  ## Extract all places and species names to loop over, and rename stuff
  sp.weighted <- colnames(dat.weighted)
  sp.codes <- cbind(sp.weighted, 1:length(sp.weighted),colSums(1*(dat.weighted>0)))
  places.codes <- cbind(rownames(dat.weighted), 1:length(rownames(dat.weighted)), rowSums(1*(dat.weighted>0)), dat_SP$Easting, dat_SP$Northing)
  
  ## Rename rows and columns
  colnames(dat.weighted) <- 1:length(sp.weighted)
  rownames(dat.weighted) <- 1:length(rownames(dat.weighted))

  ## Save species and places codes
  write.table(sp.codes, file = "../../data/properties/codes/sp_codes.csv", quote = F, row.names = F, col.names = c("sp", "id", "range"), sep=",")
  write.table(places.codes, file = "../../data/properties/codes/places_codes.csv", quote = F, row.names = F, col.names = c("place", "id", "richness", "easting", "northing"), sep=",")
  write.table(file="../../data/processed/distribution/distribution.csv", x=1*(dat.weighted>0), sep = ",",row.names = F, col.names = F, quote = F )

  ## loop over speces to generate presence and absence data
  for(i in 1:length(sp.weighted)){
    idx <- dat.weighted[,i]>0
    presence <- cbind(dat$long[idx], dat$lat[idx])
    if(length(presence)>0){
      write.table(presence, file = paste("../../data/processed/presence_lola/",as.character(i), ".csv", sep=""), quote = F, row.names = F, col.names = F, sep = ",")      
    }
  }
  for(i in 1:length(sp.weighted)){
    idx <- dat.weighted[,i]>0
    presence <- cbind(dat$Easting[idx], dat$Northing[idx])
    if(length(presence)>0){
      write.table(presence, file = paste("../../data/processed/presence_eano/",as.character(i), ".csv", sep=""), quote = F, row.names = F, col.names = F, sep = ",")      
    }
  }
}

download_climatic <- function(res=0.5, meanlat, meanlong){
  bioclim.data <- getData(name = "worldclim",
                          var = "bio",
                          res = res,
                          lat=meanlat,
                          lon=meanlong,
                          path = "../../data/raw/climatic-data/")
}

