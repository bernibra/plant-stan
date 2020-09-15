source("./models.R")
library(rethinking)
library(rstan)
library(gtools)
library(pROC)
library(raster)
library(shinystan)
library(dismo)

####
## Run poisson models for number of species per site
####

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

# Prepare all the data for the poission regressions
prepare.data <- function(variables=c("bio5_", "bio6_","bio12_"), test=F){
        # Load site data
        d <- read.csv("../../results/poisson/data/places_traits.csv", header = T)
        d <- d[d$richness>2,]
        
        # Load environmental data
        files <- list.files(path="../../data/raw/climatic-data/", pattern = "bil$", full.names = TRUE)
        files <- files[grepl(paste(variables,collapse="|"), files)]
        bioclim.data <- stack(files)
        crs(bioclim.data)<-projection
        
        # Extract climatic data
        bio <- raster::extract(bioclim.data, data.frame(d$easting, d$northing))
        colnames(bio) <- gsub("_", "", variables)
        
        # standarize variables
        for(i in 1:ncol(bio)){bio[,i] <- scale(bio[,i])}
        
        if(test){
                testing.group <- 1
                g <- kfold(1:length(d$richness), k=4)
                
                d.test <- d[g==1,]
                d.train <- d[g!=1,]
                bio.test <- bio[g==1,]
                bio.train <- bio[g!=1,]
        }else{
                d.test <- d
                d.train <- d
                bio.test <- bio
                bio.train <- bio
        }
        
        return(list(dat=d.train, bio=bio.train, dat.test=d.test, bio.test=bio.test))
}

# This first one is just a poisson regression, where all "variables" are used as predictors in a linear form.
poisson.stan.neutral <- function(d=NULL, variables=c("bio5_", "bio6_","bio12_"),  test=T){

        # Load the data
        if(is.null(d)){
            d <- prepare.data(variables = variables, test=test)
        }
        
        # create dataset for the model
        dat_1.0 <- list(N=length(d$dat$richness),
                        K=length(variables),
                        sp=as.numeric(d$dat$richness))
        
        # Set starting values for the parameters
        start_1.0 <- list(
                alpha = 3
                )
        
        # Initialize data structure
        n_chains_1.0 <- 3
        init_1.0 <- list()
        for ( i in 1:n_chains_1.0 ) init_1.0[[i]] <- start_1.0
        
        # Run stan model
        mfit_1.0 <- stan ( model_code=model1.0 ,
                           data=dat_1.0 ,
                           chains=n_chains_1.0 ,
                           cores= n_chains_1.0 ,
                           warmup=1000, iter=4000,
                           init=init_1.0 , control = list(adapt_delta = 0.95))
        
        # Make sure that priors produce sensible expectations. For example, higher sigmas will produce overrepresenation of ones and zeros.
        dens(exp(rnorm(n = 1e5, mean = 3,sd = 0.5)))
        
        # Extract posterior samples
        post <- extract.samples(mfit_1.0)
        bio <- d$bio.test
        
        # Build link function to make predictions        
        my_link <- function(bio, post){
                p <- as.matrix(bio) %*% matrix(0,ncol(bio),length(post$alpha))
                p <- exp(t(p + as.matrix(post$alpha)[col(p)]))
                return(p)
        }
        
        # Predict new points
        p_post <- my_link(post = post, bio=bio)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Visualize predictions
        par(mfrow=c(2,1))
        plot(p_mu, d$dat.test$richness, type="p", ylab="real", xlab="prediction",
             ylim=c(0, max(d$dat.test$richness)), xlim=c(0, max(d$dat.test$richness)))
        for(i in 1:length(p_mu)){
                lines(c(p_ci[1,i], p_ci[2,i]), c(d$dat.test$richness[i], d$dat.test$richness[i]),col=rangi2)
        }
        lines(c(0, max(d$dat.test$richness)), c(0, max(d$dat.test$richness)), lty=2)
        
        plot(precis(mfit_1.0, pars = "beta", depth = 2))

        return(list(mfit = mfit_1.0,d=d))
}



# This first one is just a poisson regression, where all "variables" are used as predictors in a linear form.
poisson.stan.environment <- function(d=NULL, variables=c("bio5_", "bio6_","bio12_"), test=T){
        
        # Load the data
        if(is.null(d)){
                d <- prepare.data(variables = variables, test=test)
        }
        
        # create dataset for the model
        dat_1.1 <- list(N=length(d$dat$richness),
                        K=length(variables),
                        sp=as.numeric(d$dat$richness),
                        bio=d$bio)
        
        # Set starting values for the parameters
        start_1.1 <- list(
                alpha = 3,
                beta = rep(0,length(variables))
        )
        
        # Initialize data structure
        n_chains_1.1 <- 3
        init_1.1 <- list()
        for ( i in 1:n_chains_1.1 ) init_1.1[[i]] <- start_1.1
        
        # Run stan model
        mfit_1.1 <- stan ( model_code=model1.1 ,
                           data=dat_1.1 ,
                           chains=n_chains_1.1 ,
                           cores= n_chains_1.1 ,
                           warmup=1000, iter=4000,
                           init=init_1.1 , control = list(adapt_delta = 0.95))

        # Extract posterior samples
        post <- extract.samples(mfit_1.1)
        bio <- d$bio.test
        
        # Build link function to make predictions        
        my_link <- function(bio, post){
                beta <- t(post$beta)
                p <- as.matrix(bio) %*% beta
                p <- exp(t(p + as.matrix(post$alpha)[col(p)]))
                # outputMatrix <- rpois(length(p), p)
                # dim(outputMatrix) <- dim(p)
                return(p)
        }
        
        # Predict new points
        p_post <- my_link(post = post, bio=bio)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Visualize predictions
        par(mfrow=c(2,1))
        plot(p_mu, d$dat.test$richness, type="p", ylab="real", xlab="prediction",
             ylim=c(0, max(d$dat.test$richness)), xlim=c(0, max(d$dat.test$richness)))
        for(i in 1:length(p_mu)){
             lines(c(p_ci[1,i], p_ci[2,i]), c(d$dat.test$richness[i], d$dat.test$richness[i]),col=rangi2)
        }
        lines(c(0, max(d$dat.test$richness)), c(0, max(d$dat.test$richness)), lty=2)
        
        plot(precis(mfit_1.1, pars = "beta", depth = 2))
        
        return(list(mfit = mfit_1.1, d))
}

# This first one is just a poisson regression, where all "variables" are used as predictors in a linear form.
poisson.stan.EandD <- function(d=NULL, variables=c("bio5_", "bio6_","bio12_"), test=T){
        
        # Load the data
        if(is.null(d)){
                d <- prepare.data(variables = variables, test=test)
        }
        
        # create dataset for the model
        dat_1.1 <- list(N=length(d$dat$richness),
                        K=length(variables),
                        sp=as.numeric(d$dat$richness),
                        bio=d$bio)
        
        # Set starting values for the parameters
        start_1.1 <- list(
                alpha = 3,
                beta = rep(0,length(variables))
        )
        
        # Initialize data structure
        n_chains_1.1 <- 3
        init_1.1 <- list()
        for ( i in 1:n_chains_1.1 ) init_1.1[[i]] <- start_1.1
        
        # Run stan model
        mfit_1.1 <- stan ( model_code=model1.1 ,
                           data=dat_1.1 ,
                           chains=n_chains_1.1 ,
                           cores= n_chains_1.1 ,
                           warmup=1000, iter=4000,
                           init=init_1.1 , control = list(adapt_delta = 0.95))
        
        # Extract posterior samples
        post <- extract.samples(mfit_1.1)
        bio <- d$bio.test
        
        # Build link function to make predictions        
        my_link <- function(bio, post){
                beta <- t(post$beta)
                p <- as.matrix(bio) %*% beta
                p <- exp(t(p + as.matrix(post$alpha)[col(p)]))
                # outputMatrix <- rpois(length(p), p)
                # dim(outputMatrix) <- dim(p)
                return(p)
        }
        
        # Predict new points
        p_post <- my_link(post = post, bio=bio)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Visualize predictions
        par(mfrow=c(2,1))
        plot(p_mu, d$dat.test$richness, type="p", ylab="real", xlab="prediction",
             ylim=c(0, max(d$dat.test$richness)), xlim=c(0, max(d$dat.test$richness)))
        for(i in 1:length(p_mu)){
                lines(c(p_ci[1,i], p_ci[2,i]), c(d$dat.test$richness[i], d$dat.test$richness[i]),col=rangi2)
        }
        lines(c(0, max(d$dat.test$richness)), c(0, max(d$dat.test$richness)), lty=2)
        
        plot(precis(mfit_1.1, pars = "beta", depth = 2))
        
        return(list(mfit = mfit_1.1, d))
}
