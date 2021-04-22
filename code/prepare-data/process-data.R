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

## Prepare general data
prepare.data <- function(){
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

  # Define coordinates
  dat_SP <- SpatialPointsDataFrame(coords,
                                   data = dat,
                                   proj4string = CRS(swissgrid))
  
  ## Extract coocurrence data
  dat.weighted <- dat.weighted[,7:ncol(dat.weighted)]
  
  ## Extract all places and species names to loop over, and rename stuff
  sp.weighted <- colnames(dat.weighted)
  sp.codes <- cbind(sp.weighted, 1:length(sp.weighted),colSums(1*(dat.weighted>0)))
  places.codes <- cbind(rownames(dat.weighted), 1:length(rownames(dat.weighted)), rowSums(1*(dat.weighted>0)), dat_SP$Easting, dat_SP$Northing, as.vector(as.numeric(elevation)))
  
  
  ## Rename rows and columns
  colnames(dat.weighted) <- 1:length(sp.weighted)
  rownames(dat.weighted) <- 1:length(rownames(dat.weighted))

  ## Save species and places codes
  write.table(sp.codes, file = "../../data/properties/codes/sp_codes.csv", quote = F, row.names = F, col.names = c("sp", "id", "range"), sep=",")
  write.table(places.codes, file = "../../data/properties/codes/places_codes.csv", quote = F, row.names = F, col.names = c("place", "id", "richness", "easting", "northing", "elevation"), sep=",")
  write.table(file="../../data/processed/distribution/distribution.csv", x=1*(dat.weighted>0), sep = ",",row.names = F, col.names = F, quote = F )

  return(list(dat = dat.weighted, coord = dat_SP, places = places.codes, sp = sp.codes))
}

prepare.data.sdm <- function(absences=T){
  dat <- prepare.data()
  
  ## loop over speces to generate presence and absence dat
  for(i in 1:nrow(dat$sp)){
    idx <- dat$dat[,i]>0 | rep(absences, nrow(dat$dat))
    presence <- cbind(dat$coord$Easting[idx], dat$coord$Northing[idx], dat$dat[idx, i])

    if(length(presence)>0){
      write.table(presence, file = paste("../../data/processed/sdm/",as.character(i), ".csv", sep=""), quote = F, row.names = F, col.names = c("easting", "northing", "abundance"), sep = ",")
    }
  }
}

prepare.sliding <- function(){
  dat <- prepare.data()
  elevation <- as.numeric(dat$places[,6])

  pi <- rethinking::PI(elevation, prob = c((0:10)*0.1))
  noms <- unique(names(pi))
  pi <- unique(pi)
  names(pi) <- noms
  
  for(i in 1:16){
    condition <- (elevation>pi[i] & elevation<pi[i+5])
    ocu <- as.matrix(dat$dat[condition,])
    places <- dat$places[condition,]
    sp <- dat$sp
    Easting <- dat$coord$Easting[condition]
    Northing <- dat$coord$Northing[condition]
    
    sp.codes <- cbind(sp[,1], sp[,2],colSums(1*(ocu>0)))
    places.codes <- cbind(rownames(ocu), 1:length(rownames(ocu)), rowSums(1*(ocu>0)), places[,4], places[,5], places[,6])
    print(c(pi[i], pi[i+5],sum(colSums(1*(ocu>0))>=20)))
    
    ## Save species and places codes
    write.table(sp.codes, file = paste("../../data/properties/codes/sdm-sliding/sp_codes_",i,".csv",sep=""), quote = F, row.names = F, col.names = c("sp", "id", "range"), sep=",")
    write.table(places.codes, file = paste("../../data/properties/codes/sdm-sliding/places_codes_",i,".csv",sep=""), quote = F, row.names = F, col.names = c("place", "id", "richness", "easting", "northing", "elevation"), sep=",")

    ## loop over speces to generate presence and absence dat
    for(j in 1:nrow(sp.codes)){
      idx <- ocu[,j]>0 | rep(T, nrow(ocu))
      presence <- cbind(Easting[idx], Northing[idx], ocu[idx, j])
      
      if(length(presence)>0){
        write.table(presence, file = paste("../../data/processed/sdm-sliding/partition-",i,"/",as.character(j), ".csv", sep=""), quote = F, row.names = F, col.names = c("easting", "northing", "abundance"), sep = ",")
      }
    }
  }
}

