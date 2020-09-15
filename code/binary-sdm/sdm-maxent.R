library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

prepare.data <- function(idx=1, variables = c("bio5_", "bio6_","bio12_")){
        
        # Determine geographic extent of our data
        places <- read.csv(file = "../../data/properties/codes/places_codes.csv", header = T)
        max.lat <- max(places$northing)
        min.lat <- min(places$northing)
        max.lon <- max(places$easting)
        min.lon <- min(places$easting)
        geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
        
        # Read observations
        obs.data <- read.csv(file = paste("../../data/processed/sdm/", as.character(idx), ".csv", sep = ""), header = T)

        # Define presence and absence
        presence <- as.numeric(obs.data$abundance)>0
        absence.data <- data.frame(easting = obs.data$easting[!(presence)], northing = obs.data$northing[!(presence)])        
        obs.data <- data.frame(easting = obs.data$easting[presence], northing = obs.data$northing[presence])
        
        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        files <- files[grepl(paste(variables,collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-projection
        
        # # Crop bioclim data to geographic extent of saguaro
        # bioclim.data <- crop(x = bioclim.data, y = geographic.extent)
        
        return(list(obs.data = obs.data, absence.data=absence.data, bioclim.data = bioclim.data, xlim = c(min.lon, max.lon), ylim = c(min.lat, max.lat)))
}

####
# Run species distribution model for a given species
####
species_distribution.maxent <- function(idx=1, view_plots=T, variables=c("bio5_", "bio6_","bio12_"), pseudoA = F){
        
        # # Regularizing factor
        # lr_ <- c(0,0.1,0.5,1,2,5)
        
        # Load all data
        dat <- prepare.data(idx=idx, variables = variables)

        if(view_plots){
                par(mfrow=c(2,1))
                
                # Basic map data (only for visualization purposes)
                data(wrld_simpl)
                
                # Change of projection
                wrld_simpl <- spTransform(wrld_simpl,crs(projection))
                
                # Plot the base map
                plot(wrld_simpl, 
                     xlim = dat$xlim, 
                     ylim = dat$ylim,
                     axes = TRUE, 
                     col = "grey95")
                
                # Add the points for individual observation
                points(x = dat$obs.data$easting,
                       y = dat$obs.data$northing,
                       col = "olivedrab", 
                       pch = 20, 
                       cex = 0.75)
                # And draw a little box around the graph
                box()
        }
        
        ### separate test and train data
        # Arbitrarily assign group 1 as the testing data group
        testing.group <- 1

        # Create vector of group memberships
        group.presence <- kfold(x = dat$obs.data, k = 4) # kfold is in dismo package
        group.absence <- kfold(x = dat$absence.data, k = 4) # kfold is in dismo package
        
        # Divide data between training and testing data
        presence.train <- dat$obs.data[group.presence != testing.group, ]
        presence.test <- dat$obs.data[group.presence == testing.group, ]
        
        # Generate absence test dataset
        if(pseudoA){
                absence.train <- randomPoints(mask = dat$bioclim.data[[1]],     # Provides resolution of sampling points
                                              n = 10000)
                absence.test <- randomPoints(mask = dat$bioclim.data[[1]],     # Provides resolution of sampling points
                                             n = nrow(presence.test))
                colnames(absence.test) <- colnames(presence.test)
                
        }else{
                absence.train <- dat$absence.data[group.absence != testing.group, ]
                absence.test <- dat$absence.data[group.absence == testing.group, ]
        }

        # Run maxent
        bc.model <- dismo::maxent(dat$bioclim.data, p = presence.train, a = absence.train, removeDuplicates=TRUE) # Maxent uses all the environment and presense data, regularisation from the line above. removeDuplicates is just a precaution, it is not neccessary.

        # Also run maxent with pseudo-absences for comparison purposes
        if(!(pseudoA)){
           bc.model_ <- dismo::maxent(dat$bioclim.data, p = presence.train, removeDuplicates=TRUE) # Maxent uses all the environment and presense data, regularisation from the line above. removeDuplicates is just a precaution, it is not neccessary.                
        }

        # Predict presence from model
        predict.presence <- dismo::predict(object = bc.model, x = dat$bioclim.data)
        
        if (view_plots){
                # # Plot base map
                plot(wrld_simpl,
                     xlim = dat$xlim,
                     ylim = dat$ylim,
                     axes = TRUE,
                     col = "grey95")
                
                # Add model probabilities
                plot(predict.presence, add = TRUE)
                
                # Redraw those country borders
                plot(wrld_simpl, add = TRUE, border = "grey5")
                
                # Add original observations
                points(dat$obs.data$easting, dat$obs.data$northing, col = "olivedrab", pch = 20, cex = 0.75)
                box()
        }
        
        # Prepare full dataset
        if(pseudoA){
                train.presence <- data.frame(obs=1, train=1, bc.model@presence)
                train.absence <- data.frame(obs=0, train=1, bc.model@absence)
                test.presence <- data.frame(obs=1, train=0, raster::extract(dat$bioclim.data, presence.test))
                test.absence <- data.frame(obs=0, train=0,raster::extract(dat$bioclim.data, absence.test))
        }else{
                train.presence <- data.frame(obs=1, train=1, raster::extract(dat$bioclim.data, presence.train))
                train.absence <- data.frame(obs=0, train=1, raster::extract(dat$bioclim.data, absence.train))
                test.presence <- data.frame(obs=1, train=0, raster::extract(dat$bioclim.data, presence.test))
                test.absence <- data.frame(obs=0, train=0,raster::extract(dat$bioclim.data, absence.test))
        }
        
        # Prepare data for stan model
        dataset = rbind(train.presence, train.absence, test.presence, test.absence)
        # Standarize environmental variables
        for(i in c(1:ncol(dataset))[-c(1,2)]){dataset[,i] <- scale(dataset[,i])}
        
        # Calculate AUC and ROC for maxent model
        auc_maxent <- dismo::evaluate(p = as.matrix(presence.test), a = as.matrix(absence.test), model = bc.model, x = dat$bioclim.data)
        roc_maxent <- roc(c(rep(1,nrow(presence.test)), rep(0,nrow(absence.test))), c(auc_maxent@presence, auc_maxent@absence))
        
        # In addition, I will calculate AUC and ROC for the maxent model using pseudo-absences for comparision purposes
        if(!(pseudoA)){
            auc_maxent_ <- dismo::evaluate(p = as.matrix(presence.test), a = as.matrix(absence.test), model = bc.model_, x = dat$bioclim.data)
            roc_maxent_ <- roc(c(rep(1,nrow(presence.test)), rep(0,nrow(absence.test))), c(auc_maxent_@presence, auc_maxent_@absence))
        }else{
            auc_maxent_ <- NULL
            roc_maxent_ <- NULL
        }

        return(list(dataset = dataset, model.maxent = bc.model, roc_maxent = roc_maxent,
                    auc_pseudoAmaxent=auc_maxent_, roc_pseudoAmaxent=roc_maxent_))
}

