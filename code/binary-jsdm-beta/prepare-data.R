library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")
library(MASS)
library(rethinking)

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

prepare.data <- function(variables = c("bio5_", "bio6_","bio12_"), min.occurrence=0){
        
        # Determine geographic extent of our data
        places <- read.csv(file = "../../data/properties/codes/places_codes.csv", header = T)
        max.lat <- max(places$northing)
        min.lat <- min(places$northing)
        max.lon <- max(places$easting)
        min.lon <- min(places$easting)
        geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
        
        # Load data for all species
        files <- as.numeric(gsub(".csv", "", list.files("../../data/processed/sdm/")))
        sp_codes <- read.table("../../data/properties/codes/sp_codes.csv", sep=",", header = T)
        dictionary <- read.table("../../data/properties/codes/dictionary.csv", sep=",", header = T)
        correlation_matrix_ids <- read.table("../../data/properties/codes/correlation_matrix_ids.csv", sep="\t", header = F)
        denvironment <- as.matrix(read.table("../../data/properties/distance-matrices/environment.csv", sep=","))
        dvariation <- as.matrix(read.table("../../data/properties/distance-matrices/variation.csv", sep=","))
        
        # Check that there aren't unnexpected files
        if(!all(sort(files)==1:length(files))){
            stop("Odd files in the folder")    
        }
        
        # Prepare main file
        obs.data <- data.frame()
        name.idx <- c()
        kdx <- 1
        
        # Read observations
        for(idx in 1:length(files)){
             if(sp_codes$range[sp_codes$id==idx]<min.occurrence){
                     next
             }
             obs.data_ <- read.csv(file = paste("../../data/processed/sdm/", as.character(idx), ".csv", sep = ""), header = T)
             ### Ok so I think the order for the correlation matrix is the same, but I do need to double-check
             ### This is a very dumb way of doing just that. I didn't want to think.
             new.name <- as.character(dictionary$new.names[as.character(dictionary$old.names)==as.character(sp_codes$sp[idx])])
             name.idx <- c(name.idx,correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name])
             ###
             obs.data_$id <- kdx
             obs.data_$real.id <- idx
             obs.data_$mmsbm.id <- correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name]
             obs.data <- rbind(obs.data, obs.data_)
             kdx <- kdx+1
        }
        
        # reshape correlation matrices
        denvironment <- denvironment[,name.idx]
        denvironment <- denvironment[name.idx,]
        dvariation <- dvariation[,name.idx]
        dvariation <- dvariation[name.idx,]
        
        # rename cols and rows
        colnames(denvironment) <- 1:ncol(denvironment)
        rownames(denvironment) <- 1:nrow(denvironment)
        colnames(dvariation) <- 1:ncol(dvariation)
        rownames(dvariation) <- 1:nrow(dvariation)
        
        obs.data$obs <- 1*(obs.data$abundance>0)

        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        files <- files[grepl(paste(variables,collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-projection
        
        return(list(obs.data = obs.data, bioclim.data = bioclim.data, xlim = c(min.lon, max.lon), ylim = c(min.lat, max.lat), 
                    denv = denvironment, dvar = dvariation))
}

# Generate fake data to test the extent to which the model works
simulated.data <- function(){
        # Define system dimensions
        N <- 20
        sites <- 100
        
        # Environmental predictors for each site
        e1 <- rnorm(sites)
        e2 <- rnorm(sites)
        
        # uncorrelated coefficients for each species and parameters
        alpha <- rnorm(N)
        sigma1 <- 0.3 # sd beta1
        mean1 <- -1 # mean beta1
        sigma2 <- 0.4 # sd beta2
        mean2 <- 1.5 # mean beta2
        rho <- 0.4 # correlation between betas
        
        # Coefficients for generating the variance-covariance matrix
        nu <- 3
        s <- 0.5
        Dis <- (as.matrix(dist(1:N))/N)
        
        # coefficients for each species
        Sigma <- nu*exp(-1/(s*s)*(Dis^2)) + diag(N)*sigma1*sigma1
        z1 <- mvrnorm(mu = rep(mean1, times = N), Sigma = Sigma)
        Sigma <- nu*exp(-1/(s*s)*(Dis^2)) + diag(N)*sigma2*sigma2
        z2 <- mvrnorm(mu = rep(mean2, times = N), Sigma = Sigma)
        
        # Generate correlations
        beta1 <- z1
        beta2 <- z2
        sigma_beta1 <- abs(rnorm(N, 0,0.2))
        sigma_beta2 <- abs(rnorm(N, 0,0.2))
        
        vec <- c(1:round(N*0.5), 1:round(N*0.5))
        Dis_sigma <- as.matrix(dist(vec))+1-diag(N)
        Dis_sigma <- (Dis_sigma/max(Dis_sigma))
        
        Sigma <- 1*exp(-1/(0.3*0.3)*(Dis_sigma^2)) + diag(N)*0.2
        sigma_beta1 <- exp(mvrnorm(mu = rep(-1, times = N), Sigma = Sigma))
        Sigma <- 1*exp(-1/(0.2*0.2)*(Dis_sigma^2)) + diag(N)*0.1
        sigma_beta2 <- exp(mvrnorm(mu = rep(-2, times = N), Sigma = Sigma))

        # Simulate data
        dataset <- expand.grid(site=1:sites, id=1:N)
        dataset$S1 <- e1[dataset$site]
        dataset$S2 <- e2[dataset$site]
        dataset$beta1 <- beta1[dataset$id] 
        dataset$beta2 <- beta2[dataset$id] 
        dataset$alpha <- alpha[dataset$id]
        dataset$sigma_beta1 <- sigma_beta1[dataset$id] 
        dataset$sigma_beta2 <- sigma_beta2[dataset$id]
        dataset$p <- inv_logit(alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2 - sigma_beta2[dataset$id]*(beta2[dataset$id] - dataset$S2)**2)              
        
        dataset$obs <- rbinom(n = length(dataset$S1), size = 1, prob = dataset$p)
        dataset <- data.frame(id=dataset$id, obs=dataset$obs, alpha=dataset$alpha, beta1=dataset$beta1, beta2=dataset$beta2, sigma_beta1=dataset$sigma_beta1, sigma_beta2=dataset$sigma_beta2,  S1=dataset$S1, S2=dataset$S2)
        Y <- matrix(dataset$obs, sites, N)
        X1 <- matrix(dataset$S1, sites, N)[,1]
        X1 <- matrix(dataset$S2, sites, N)[,1]

        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, Y=Y, X1=X1, X2=X2))
}

####
# Run species distribution model for a given species
####
species_distribution.data <- function(variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_"),
                                      pca=F, ndim=2,
                                      simulated=F, min.occurrence=0){
        
        if(simulated){
                dataset <- simulated.data()
                return(dataset)
        }else{
                if(pca){
                        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")  
                }
                
                # Load all data
                dat <- prepare.data(variables = variables, min.occurrence=min.occurrence)
                
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
                for(i in c(1:ncol(dataset))[-c(1:(ncol(dataset)-ndim))]){dataset[,i] <- scale(dataset[,i])}
                
                return(list(dataset=dataset, corr=dat$denv, corr2=dat$dvar))
        }
}

playing.with.multivariate <- function(n, s, c, mu, N){
        vec <- c(1:round(N*0.5), 1:round(N*0.5))
        # vec <- 1:N
        Dis <- dist(vec)
        Sigma <- n*exp(-(1/(s*s))*((as.matrix(Dis)/max(Dis))^2)) + diag(N)*c
        z1 <- exp(mvrnorm(mu = rep(mu, times = N), Sigma = Sigma))
        plot(1:N, z1)
        return(z1)
}

