library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

prepare.data <- function(variables = c("bio5_", "bio6_","bio12_")){
        
        # Determine geographic extent of our data
        places <- read.csv(file = "../../data/properties/codes/places_codes.csv", header = T)
        max.lat <- max(places$northing)
        min.lat <- min(places$northing)
        max.lon <- max(places$easting)
        min.lon <- min(places$easting)
        geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
        
        # Load data for all species
        files <- as.numeric(gsub(".csv", "", list.files("../../data/processed/sdm/")))
        
        # Check that there aren't unnexpected files
        if(!all(sort(files)==1:length(files))){
            stop("Odd files in the folder")    
        }
        
        # Prepare main file
        obs.data <- data.frame()
        
        # Read observations
        for(idx in 1:length(files)){
             obs.data_ <- read.csv(file = paste("../../data/processed/sdm/", as.character(idx), ".csv", sep = ""), header = T)                
             obs.data_$id <- idx
             obs.data <- rbind(obs.data, obs.data_)
        }
        
        obs.data$obs <- 1*(obs.data$abundance>0)

        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        files <- files[grepl(paste(variables,collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-projection
        
        return(list(obs.data = obs.data, bioclim.data = bioclim.data, xlim = c(min.lon, max.lon), ylim = c(min.lat, max.lat)))
}

####
# Run species distribution model for a given species
####
species_distribution.data <- function(variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_"), pca=F, ndim=2){
        
        if(pca){
                variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")  
        }
        
        # Load all data
        dat <- prepare.data(variables = variables)
        
        # extract environmental data
        clim <- as.data.frame(raster::extract(dat$bioclim.data,data.frame(easting = dat$obs.data$easting, northing = dat$obs.data$northing)))
        
        # Find main axes
        if(pca){
                pca.clim <- prcomp(clim, center = TRUE, scale = TRUE) 
                # Prepare full dataset
                dataset = cbind(dat$obs.data, pca.clim$x[,c(1:ndim)])
        }else{
                # Prepare full dataset
                dataset = cbind(dat$obs.data, clim)
        }
        
        # Standarize environmental variables
        for(i in c(1:ncol(dataset))[-c(1,2,3,4,5)]){dataset[,i] <- scale(dataset[,i])}
        
        return(dataset)
}

