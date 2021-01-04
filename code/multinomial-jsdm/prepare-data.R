library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")
library(MASS)
library(MCMCpack)

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
        dtraits <- as.matrix(read.table("../../data/properties/distance-matrices/trait.csv", sep=","))
        
        indicator <- read.table("../../data/properties/codes/temperature_indicator.csv", sep=",")
        
        # Check that there aren't unnexpected files
        if(!all(sort(files)==1:length(files))){
            stop("Odd files in the folder")    
        }
        
        # Prepare main file
        obs.data <- data.frame()
        name.idx <- c()
        kdx <- 1
        Tind <- c()
        
        # Read observations
        for(idx in 1:length(files)){
             if(sp_codes$range[sp_codes$id==idx]<min.occurrence){
                     next
             }
             obs.data_ <- read.csv(file = paste("../../data/processed/sdm/", as.character(idx), ".csv", sep = ""), header = T)
             ### Ok so I think the order for the correlation matrix is the same, but I do need to double-check
             ### This is a very dumb way of doing just that. I didn't want to think.
             new.name <- as.character(dictionary$new.names[as.character(dictionary$old.names)==as.character(sp_codes$sp[idx])])
             Tind <- rbind(Tind, c(kdx, new.name, as.character(indicator$nflor.T[new.name==indicator$nflor.spnames])))
             name.idx <- c(name.idx,correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name])
             ###
             obs.data_$id <- kdx
             obs.data_$real.id <- idx
             obs.data_$mmsbm.id <- correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name]
             obs.data <- rbind(obs.data, obs.data_)
             kdx <- kdx+1
        }
        
        write.table(Tind, "../../data/properties/codes/temperature_indicator_reindexed.csv", sep=",")
        
        # reshape correlation matrices
        denvironment <- denvironment[,name.idx]
        denvironment <- denvironment[name.idx,]
        dvariation <- dvariation[,name.idx]
        dvariation <- dvariation[name.idx,]
        dtraits <- dtraits[,name.idx]
        dtraits <- dtraits[name.idx,]
        
        # rename cols and rows
        colnames(denvironment) <- 1:ncol(denvironment)
        rownames(denvironment) <- 1:nrow(denvironment)
        colnames(dvariation) <- 1:ncol(dvariation)
        rownames(dvariation) <- 1:nrow(dvariation)
        colnames(dtraits) <- 1:ncol(dtraits)
        rownames(dtraits) <- 1:nrow(dtraits)
        
        obs.data$obs <- factor(obs.data$abundance,
               levels = sort(unique(obs.data$abundance)),
               labels = 1:length(unique(obs.data$abundance)))

        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        files <- files[grepl(paste(variables,collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-projection
        
        return(list(obs.data = obs.data, bioclim.data = bioclim.data, xlim = c(min.lon, max.lon), ylim = c(min.lat, max.lat), 
                    denv = denvironment, dvar = dvariation, dtrait=dtraits))
}

# Generate fake data to test the extent to which the model works
simulated.data <- function(){
        
        # Define system dimensions
        N <- 50
        sites <- 100
        
        # Environmental predictors for each site
        e1 <- rnorm(sites)
        e2 <- rnorm(sites)
        
        # uncorrelated coefficients for each species and parameters
        levels <- 1:5
        alpha <- rnorm(N)
        b01 <- 1
        b02 <- 0.05
        b03 <- -0.05
        b04 <- -1
        alpha1 <- 1
        alpha2 <- 0.05
        alpha3 <- -0.05
        alpha4 <- -1

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
        dataset$alpha1 <- rep(alpha1, length(dataset$id))
        dataset$alpha2 <- rep(alpha2, length(dataset$id))
        dataset$alpha3 <- rep(alpha3, length(dataset$id))
        dataset$alpha4 <- rep(alpha4, length(dataset$id))
        
        dataset$sigma_beta1 <- sigma_beta1[dataset$id] 
        dataset$sigma_beta2 <- sigma_beta2[dataset$id]
        dataset$prob_2to5 <- inv_logit(alpha1 + alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2 - sigma_beta2[dataset$id]*(beta2[dataset$id] - dataset$S2)**2)              
        dataset$prob_3to5 <- inv_logit(alpha2 + alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2 - sigma_beta2[dataset$id]*(beta2[dataset$id] - dataset$S2)**2)              
        dataset$prob_4to5 <- inv_logit(alpha3 + alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2 - sigma_beta2[dataset$id]*(beta2[dataset$id] - dataset$S2)**2)              
        dataset$prob_5 <- inv_logit(alpha4 + alpha[dataset$id] - sigma_beta1[dataset$id]*(beta1[dataset$id] - dataset$S1)**2 - sigma_beta2[dataset$id]*(beta2[dataset$id] - dataset$S2)**2)              
        dataset$prob_1 <- 1 - dataset$prob_2to5
        dataset$prob_2 <- dataset$prob_2to5 - dataset$prob_3to5
        dataset$prob_3 <- dataset$prob_3to5 - dataset$prob_4to5
        dataset$prob_4 <- dataset$prob_4to5 - dataset$prob_5
        
        obs <- c()
        for (i in 1:length(dataset$prob_4)) {
                obs[i] <- sample(
                        x = c(1:5), 
                        size = 1, 
                        prob = c(dataset$prob_1[i], dataset$prob_2[i], dataset$prob_3[i], dataset$prob_4[i], dataset$prob_5[i])
                )
        }
        
        dataset$obs <- obs
        dataset <- data.frame(id=dataset$id, obs=dataset$obs, alpha=dataset$alpha, alpha1=dataset$alpha1, alpha2=dataset$alpha2, alpha3=dataset$alpha3, alpha4=dataset$alpha4,  beta1=dataset$beta1, beta2=dataset$beta2, sigma_beta1=dataset$sigma_beta1, sigma_beta2=dataset$sigma_beta2,  S1=dataset$S1, S2=dataset$S2)                
        return(list(dataset=dataset, corr=Dis, corr2=Dis_sigma, corr3=Dis_sigma))

}

####
# Run species distribution model for a given species
####
species_distribution.data <- function(variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_"),
                                      pca=F, ndim=2,
                                      simulated=F, simulated.type="linear.corr", min.occurrence=0){
        
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
                
                return(list(dataset=dataset, corr=dat$denv, corr2=dat$dvar, corr3=dat$dtrait))
        }
}


