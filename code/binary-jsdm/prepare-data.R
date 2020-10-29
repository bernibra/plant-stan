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

# Generate fake data to test the extent to which the model works
simulated.data <- function(simulated.type="linear"){
        
        if(simulated.type=="linear"){
                # Define system dimensions
                N <- 20
                sites <- 200
                
                # Environmental predictors for each site
                e1 <- rnorm(sites)
                e2 <- rnorm(sites)
                
                # uncorrelated coefficients for each species
                alpha <- rnorm(N)
                z1 <- rnorm(N)
                z2 <- rnorm(N)
                sigma1 <- 0.3
                mean1 <- -1
                sigma2 <- 1.3
                mean2 <- 1.5
                
                # Generate correlations
                rho <- 0.4
                beta1 <- z1 * sigma1 + mean1
                beta2 <- (rho*z1 + sqrt(1-rho^2)*z2)*sigma2 + mean2
                
                # Simulate data
                dataset <- expand.grid(site=1:sites, id=1:N)
                dataset$S1 <- e1[dataset$site]
                dataset$S2 <- e2[dataset$site]
                dataset$beta1 <- beta1[dataset$id] 
                dataset$beta2 <- beta2[dataset$id] 
                dataset$alpha <- alpha[dataset$id]
                dataset$p <- inv_logit(alpha[dataset$id] + beta1[dataset$id] * dataset$S1  + beta2[dataset$id] * dataset$S2 )
                dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
                dataset <- data.frame(id=dataset$id, obs=dataset$obs, alpha=dataset$alpha, beta1=dataset$beta1, beta2=dataset$beta2, S1=dataset$S1, S2=dataset$S2)                
        }else{
                # Define system dimensions
                N <- 20
                sites <- 200
                
                # Environmental predictors for each site
                e1 <- rnorm(sites)
                e2 <- rnorm(sites)
                
                # uncorrelated coefficients for each species
                alpha <- rnorm(N)
                z1 <- rnorm(N)
                z2 <- rnorm(N)
                sigma1 <- 0.3
                mean1 <- -1
                sigma2 <- 1.3
                mean2 <- 1.5
                
                # Generate correlations
                rho <- 0.4
                beta1 <- z1 * sigma1 + mean1
                beta2 <- (rho*z1 + sqrt(1-rho^2)*z2)*sigma2 + mean2
                
                # Simulate data
                dataset <- expand.grid(site=1:sites, id=1:N)
                dataset$S1 <- e1[dataset$site]
                dataset$S2 <- e2[dataset$site]
                dataset$beta1 <- beta1[dataset$id] 
                dataset$beta2 <- beta2[dataset$id] 
                dataset$alpha <- alpha[dataset$id]
                dataset$p <- inv_logit(alpha[dataset$id] + beta1[dataset$id] * dataset$S1  + beta2[dataset$id] * dataset$S2 )
                dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
                dataset <- data.frame(id=dataset$id, obs=dataset$obs, alpha=dataset$alpha, beta1=dataset$beta1, beta2=dataset$beta2, S1=dataset$S1, S2=dataset$S2)
        }

        return(dataset)
}

####
# Run species distribution model for a given species
####
species_distribution.data <- function(variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_"),
                                      pca=F, ndim=2,
                                      simulated=F, simulated.type="linear"){
        
        if(simulated){
                if(!(simulated.type %in% c("linear", "gaussian"))){
                        stop(paste("'", simulated.type, "' is not a valid 'simulated.type'", sep=""))
                }
                dataset <- simulated.data(simulated.type=simulated.type)
                return(dataset)
        }else{
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
}

