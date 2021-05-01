library("sp")
library("raster")
library("maptools")
library("dismo")
library("rJava")
library(SDMTools)
source("../prepare-data/useful-tools.R")
library(MASS)
library(gridExtra)
library(grid)

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

prepare.data <- function(partition=1, variables = c("bio5_", "bio6_","bio12_"), min.occurrence=0, elevation=F){
        
        # Determine geographic extent of our data
        places <- read.csv(file = paste("../../data/properties/codes/sdm-sliding/places_codes_",as.character(partition),".csv", sep=""), header = T)
        # places$elevation <- as.vector(scale(places$elevation))
        max.lat <- max(places$northing)
        min.lat <- min(places$northing)
        max.lon <- max(places$easting)
        min.lon <- min(places$easting)
        geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))
        
        # Load data for all species
        files <- as.numeric(gsub(".csv", "", list.files(paste("../../data/processed/sdm-sliding/partition-",as.character(partition),"/", sep=""))))
        sp_codes <- read.table(paste("../../data/properties/codes/sdm-sliding/sp_codes_",as.character(partition),".csv", sep=""), sep=",", header = T)
        dictionary <- read.table("../../data/properties/codes/dictionary.csv", sep=",", header = T)
        correlation_matrix_ids <- read.table("../../data/properties/codes/correlation_matrix_ids.csv", sep="\t", header = F)
        denvironment <- as.matrix(read.table("../../data/properties/distance-matrices/environment.csv", sep=","))
        dvariation <- as.matrix(read.table("../../data/properties/distance-matrices/variation.csv", sep=","))
        dtraits <- as.matrix(read.table("../../data/properties/distance-matrices/trait.csv", sep=","))
        
        indicator <- read.table("../../data/properties/codes/temperature_indicator.csv", sep=",")
        neophytes <- read.table("../../data/properties/codes/neophytes-list.csv", sep=",")
        tendency <- read.table("../../data/properties/codes/change-tendency.csv", sep=",")
        competitive <- read.table("../../data/properties/codes/competitive_indicator.csv", sep=",")
        
        # Check that there aren't unnexpected files
        if(!all(sort(files)==1:length(files))){
            stop("Odd files in the folder")    
        }
        
        # Prepare main file
        obs.data <- data.frame()
        name.idx <- c()
        kdx <- 1
        Tind <- c()
        NEO <- c()
        Tend <- c()
        compet <- c()
        
        # Read observations
        for(idx in 1:length(files)){
             if(sp_codes$range[sp_codes$id==idx]<min.occurrence){
                     next
             }
             obs.data_ <- read.csv(file = paste("../../data/processed/sdm-sliding/partition-",as.character(partition),"/", as.character(idx), ".csv", sep = ""), header = T)
             ### Ok so I think the order for the correlation matrix is the same, but I do need to double-check
             ### This is a very dumb way of doing just that. I didn't want to think.
             new.name <- as.character(dictionary$new.names[as.character(dictionary$old.names)==as.character(sp_codes$sp[idx])])
             Tind <- rbind(Tind, c(kdx, new.name, as.character(indicator$nflor.T[new.name==indicator$nflor.spnames])))
             NEO <- rbind(NEO, c(kdx, new.name, neophytes$neo[new.name==neophytes$names]))
             Tend <- rbind(Tend, c(kdx, new.name, tendency$decrease[new.name==tendency$names], tendency$decrease.low[new.name==tendency$names], tendency$increase[new.name==tendency$names], tendency$other[new.name==tendency$names]))
             compet <- rbind(compet, c(kdx, new.name, as.character(competitive$nflor.KS[new.name==competitive$nflor.spnames])))
             
             name.idx <- c(name.idx,correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name])
             ###
             obs.data_$id <- kdx
             obs.data_$real.id <- idx
             obs.data_$mmsbm.id <- correlation_matrix_ids$V1[as.character(correlation_matrix_ids$V2)==new.name]
             if(elevation){
                     obs.data_$elevation <- places$elevation
             }
             obs.data <- rbind(obs.data, obs.data_)
             kdx <- kdx+1
        }
        
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
        
        obs.data$obs <- 1*(obs.data$abundance>0)
        obs.data$abundance <- factor(obs.data$abundance,
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

####
# Run species distribution model for a given species
####
species_distribution.data <- function(partition=1, variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_"),
                                      pca=F, ndim=1,
                                      simulated=F, simulated.type="linear.corr", min.occurrence=0, elevation=F){
        
        if(pca){
                        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")  
                }
                
                if(elevation){
                        variables=c("bio1_")  
                }
                
                # Load all data
                dat <- prepare.data(partition=partition, variables = variables, min.occurrence=min.occurrence, elevation=elevation)
                
                # extract environmental data
                clim <- as.data.frame(raster::extract(dat$bioclim.data,data.frame(easting = dat$obs.data$easting, northing = dat$obs.data$northing)))
                
                # Find main axes
                if(pca){
                        # colnames(clim) <- unlist(lapply(strsplit(colnames(clim), split = "_"), function(xx) xx[1]))
                        pca.clim <- prcomp(clim, center = TRUE, scale = TRUE) 
                        # Prepare full dataset
                        dataset = cbind(dat$obs.data, pca.clim$x[,c(1:ndim)])
                        
                        # p1 <- fviz_eig(pca.clim)+
                        #         scale_y_continuous(expand = expansion(add = c(0, 0)))+
                        #         theme_bw()+
                        #         theme(plot.title = element_blank(), text = element_text(size=10))
                        # p2 <- fviz_pca_var(pca.clim,
                        #              col.var = "contrib", # Color by contributions to the PC
                        #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        #              repel = T,     # Avoid text overlapping
                        # ) + theme(plot.title = element_blank(), text = element_text(size=10))
                        # 
                        # p <- grid.arrange(grobs=list(p1, p2), ncol=2, nrow=1, widths=c(0.65,1))
                        
                }else{
                        if(elevation){
                                # Prepare full dataset
                                dataset = cbind(dat$obs.data, data.frame(elevationstd=dat$obs.data$elevation))
                        }else{
                                # Prepare full dataset
                                dataset = cbind(dat$obs.data, clim)
                        }
                }
                
                # Standarize environmental variables
                for(i in c(1:ncol(dataset))[-c(1:(ncol(dataset)-ndim))]){dataset[,i] <- as.vector(scale(dataset[,i]))}
                
                return(list(dataset=dataset, corr=dat$denv, corr2=dat$dvar, corr3=dat$dtrait))
}


