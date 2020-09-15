source("./sdm-maxent.R")
source("./models.R")
library(rethinking)
library(rstan)
library(gtools)
library(pROC)

####
## Run poisson models for number of species per site
####

# Define projection as global variable
projection <- "+proj=somerc +init=world:CH1903"

# Prepare all the data for the poission regressions
prepare.data <- function(variables=c("bio5_", "bio6_","bio12_")){
        # Load site data
        d <- read.csv("../../data/properties/codes/places_codes.csv", header = T)
        
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
        return(list(dat=d, bio=bio))
}

# This first one is just a poisson regression, where all "variables" are used as predictors in a linear form.
poisson.stan.neutral <- function(variables=c("bio5_", "bio6_","bio12_")){
        
        d <- prepare.data(variables = variables)
        
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
        priors <- function(bio, mean=0, sigma=1){
                N <- dim(bio)[1]
                K <- 100
                p <- matrix(0,K,N)
                for(j in 1:K){
                        for(i in 1:N){
                                p[j, i] <- inv.logit(rnorm(n = 1, mean = mean, sd = sigma) + sum(bio[i,] * rnorm(n = 3, mean = mean, sd = sigma)))
                        }
                }
                
                # # The same can be done much more efficiently
                # K <- 1000
                # L <- dim(bio)[2]
                # p <- as.matrix(bio) %*% matrix(rnorm(n = L*K, mean = mean, sd = sigma), nrow = L, ncol = K)
                # p <- inv.logit(t(p + rnorm(n=K, mean = mean, sd = sigma)[col(p)]))
                
                dens(as.vector(p))
        }
        
        # Build link function to make predictions        
        link_1.1 <- function(dat, post){
                N <- dim(post$alpha)[1]
                M <- nrow(dat)
                K <- ncol(dat)
                beta <- t(post$beta)
                
                p <- as.matrix(dat) %*% beta
                p <- inv.logit(t(p + as.matrix(post$alpha)[col(p)]))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[!(training.data)]
        dat <- d$dataset[!(training.data),3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_1.1)
        
        # Make predictions with test data
        p_post <- link_1.1(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Have a look at the estimates
        plot(precis(mfit_1.1,depth=2))
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc_obj <- auc(roc_obj)
        
        # Compare maxent and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("regular logistic regression")
        plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
        title("maxent model")
        
        return(list(mfit = mfit_1.1, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d))
}



# This first one is just a poisson regression, where all "variables" are used as predictors in a linear form.
poisson.stan.environment <- function(variables=c("bio5_", "bio6_","bio12_"), loglik=F){
        
        d <- prepare.data(variables = variables)
        
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

        # Make sure that priors produce sensible expectations. For example, higher sigmas will produce overrepresenation of ones and zeros.
        priors <- function(bio, mean=0, sigma=1){
                N <- dim(bio)[1]
                K <- 100
                p <- matrix(0,K,N)
                for(j in 1:K){
                        for(i in 1:N){
                                p[j, i] <- inv.logit(rnorm(n = 1, mean = mean, sd = sigma) + sum(bio[i,] * rnorm(n = 3, mean = mean, sd = sigma)))
                        }
                }
                
                # # The same can be done much more efficiently
                # K <- 1000
                # L <- dim(bio)[2]
                # p <- as.matrix(bio) %*% matrix(rnorm(n = L*K, mean = mean, sd = sigma), nrow = L, ncol = K)
                # p <- inv.logit(t(p + rnorm(n=K, mean = mean, sd = sigma)[col(p)]))
                
                dens(as.vector(p))
        }

        # Build link function to make predictions        
        link_1.1 <- function(dat, post){
                N <- dim(post$alpha)[1]
                M <- nrow(dat)
                K <- ncol(dat)
                beta <- t(post$beta)
                
                p <- as.matrix(dat) %*% beta
                p <- inv.logit(t(p + as.matrix(post$alpha)[col(p)]))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[!(training.data)]
        dat <- d$dataset[!(training.data),3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_1.1)
        
        # Make predictions with test data
        p_post <- link_1.1(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Have a look at the estimates
        plot(precis(mfit_1.1,depth=2))
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc_obj <- auc(roc_obj)
        
        # Compare maxent and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("regular logistic regression")
        plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
        title("maxent model")
        
        return(list(mfit = mfit_1.1, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d))
}
