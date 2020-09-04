library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")

####
# Run species distribution model for a given species
####
species_distribution.maxent <- function(idx=1, view_plots=T){
        
        # # Regularizing factor
        # lr_ <- c(0,0.1,0.5,1,2,5)
        
        # Define projection
        projection <- "+proj=somerc +init=world:CH1903"
        
        # Determine geographic extent of our data
        places <- read.csv(file = "../../data/properties/codes/places_codes.csv", header = T)
        max.lat <- max(places$northing)
        min.lat <- min(places$northing)
        max.lon <- max(places$easting)
        min.lon <- min(places$easting)
        geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
        
        # Read observations
        obs.data <- read.csv(file = paste("../../data/processed/sdm/", as.character(idx), ".csv", sep = ""), header = T)
        obs.data <- data.frame(easting = obs.data$easting, northing = obs.data$northing)
        
        if(view_plots){
                par(mfrow=c(2,1))
                
                # Basic map data (only for visualization purposes)
                data(wrld_simpl)
                
                # Change of projection
                wrld_simpl <- spTransform(wrld_simpl,crs(projection))
                
                # Plot the base map
                plot(wrld_simpl, 
                     xlim = c(min.lon, max.lon), 
                     ylim = c(min.lat, max.lat),
                     axes = TRUE, 
                     col = "grey95")
                
                # Add the points for individual observation
                points(x = obs.data$easting,
                       y = obs.data$northing,
                       col = "olivedrab", 
                       pch = 20, 
                       cex = 0.75)
                # And draw a little box around the graph
                box()
        }
        
        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        # files <- files[grepl(paste(c("bio4_","bio5_","bio6_","bio8_","bio12_","bio15_","bio17_"),collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-"+proj=somerc +init=world:CH1903"
        
        # Crop bioclim data to geographic extent of saguaro
        bioclim.data <- crop(x = bioclim.data, y = geographic.extent)

        
        ### separate test and train data
        # Arbitrarily assign group 1 as the testing data group
        testing.group <- 1

        # Create vector of group memberships
        group.presence <- kfold(x = obs.data, k = 4) # kfold is in dismo package

        # Divide data between training and testing data
        presence.train <- obs.data[group.presence != testing.group, ]
        presence.test <- obs.data[group.presence == testing.group, ]
        csv <- rbind(presence.train,presence.test)
        
        # Run maxent
        bc.model <- dismo::maxent(bioclim.data, presence.train, removeDuplicates=TRUE) # Maxent uses all the environment and presense data, regularisation from the line above. removeDuplicates is just a precaution, it is not neccessary.

        # Predict presence from model
        predict.presence <- dismo::predict(object = bc.model, x = bioclim.data)
        
        if (view_plots){
                # # Plot base map
                plot(wrld_simpl,
                     xlim = c(min.lon, max.lon),
                     ylim = c(min.lat, max.lat),
                     axes = TRUE,
                     col = "grey95")
                
                # Add model probabilities
                plot(predict.presence, add = TRUE)
                
                # Redraw those country borders
                plot(wrld_simpl, add = TRUE, border = "grey5")
                
                # Add original observations
                points(obs.data$easting, obs.data$northing, col = "olivedrab", pch = 20, cex = 0.75)
                box()
        }

        # 
        # if(test){
        #         #Background points
        #         test_loc <- csv[csv[,3]=="test",1:2]
        #         bg = randomPoints(predict.presence,n=1000)  # should be n=10000
        #         bg = as.data.frame(bg)
        #         test_prior = sapply(seq(dim(test_loc)[1]),function(x) extract(predict.presence,test_loc[x,])) #extract suitability index values for presense data on the west coast
        #         bg_prior = sapply(seq(dim(bg)[1]),function(x) extract(predict.presence,bg[x,])) # extract suitability index for pseudo-absenses
        #         prior_stat = SDM_Stat(test_prior,bg_prior) # calculate AUC & TSS.
        #         return(prior_stat)
        # }else{
        #         return(NA)
        # }

}

