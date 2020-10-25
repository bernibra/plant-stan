source("./sdm-maxent.R")
source("./models.R")
library(rethinking)
library(rstan)
library(gtools)
library(pROC)

####
## Run species distribution model for a given species considering ordered categorical observations.
## I wrote the following models:
# - Simple ordered-multinomial regression.
# - Zero-inflated ordered-multinomial regression.
####

# This first one is just a ordered categorical regression, all "variables" are used as predictors in a linear form.
ordered.categorical.stan <- function(d = NULL,idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F, splitdata=F, show.plot=T, ofolder="../../results/models/"){
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }

        # Split data
        if(splitdata){
            training.data <- d$dataset$train==1
            test.data <- !(training.data)
        }else{
            training.data <- d$dataset$train>-1
            test.data <- training.data
        }
        
        # Prepare training data for stan model
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_1.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(obs))-1,
                        obs=obs,
                        bio=bio)
        
        # Set starting values for the parameters
        start_1.1 <- list(
                alpha = 1:dat_1.1$L/dat_1.1$L-length(dat_1.1$L)*0.5,
                beta = rep(0,length(variables))
        )
        
        # Initialize data structure
        n_chains_1.1 <- 3
        init_1.1 <- list()
        for ( i in 1:n_chains_1.1 ) init_1.1[[i]] <- start_1.1
        
        # Should I calculate the likelihood?
        if(loglik){
           model_code=model1.1
        }else{
           model_code=model1.0     
        }
        
        # Run stan model
        mfit_1.1 <- stan ( model_code=model_code ,
                           data=dat_1.1 ,
                           chains=n_chains_1.1 ,
                           cores= n_chains_1.1 ,
                           warmup=1000, iter=4000,
                           init=init_1.1 , control = list(adapt_delta = 0.95))
        
        # Build link function to make predictions        
        link_1.1 <- function(dat, post){
                N <- dim(post$alpha)[1]
                L <- dim(post$alpha)[2]
                M <- nrow(dat)
                K <- ncol(dat)
                beta <- t(post$beta)
                
                p <- as.matrix(dat) %*% beta
                p <- lapply(1:L, function(x) inv.logit(t(as.matrix(post$alpha[,x])[col(p)]-p)))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[test.data]
        dat <- d$dataset[test.data,3:ncol(d$dataset)]
        ord <- order(d$dataset$bio6_tminc_8110)
        
        obs <- obs[ord]
        dat <- dat[ord,]
        
        # Sample from the posterior
        post <- extract.samples(mfit_1.1)

        if(show.plot){
                # Have a look at the estimates
                plot(precis(mfit_1.1,depth=2))
                
                #Visualize plot
                par(mfrow=c(3,3))
                condi1 <- c(-1,-1,-1,0,0,0,1,1,1)
                condi2 <- c(-1,0,1,-1,0,1,-1,0,1)
                text <- c("bio12=-1 bio5=-1", "bio12=-1 bio5=0", "bio12=-1 bio5=1",
                          "bio12=0 bio5=-1", "bio12=0 bio5=0", "bio12=0 bio5=1",
                          "bio12=1 bio5=-1", "bio12=1 bio5=0", "bio12=1 bio5=1")
                
                for(k in 1:9){
                        dat_ <- dat
                        dat_$bio12_p_8110 <- dat_$bio12_p_8110*0+condi1[k]
                        dat_$bio5_tmaxw_8110 <- dat_$bio5_tmaxw_8110*0+condi2[k]
                        
                        # Make predictions with test data
                        p_post <- link_1.1(post = post, dat=dat_)
                        p_mu <- lapply(1:length(p_post), function(x) apply( p_post[[x]] , 2 , mean ))
                        p_ci <- lapply(1:length(p_post), function(x) apply( p_post[[x]] , 2 , PI ))
                        
                        # Calculate AUC and stuff
                        plot( NULL , type="n" , xlab="tmin" , ylab="probability" ,
                              xlim=c(min(dat$bio6_tminc_8110),max(dat$bio6_tminc_8110)) , ylim=c(0,1) , yaxp=c(0,1,2) )
                        
                        # for ( s in 1:50 ) {
                        #         for ( i in 1:3 ) lines( dat$bio6_tminc_8110 , p_post[[i]][s,] , col=col.alpha("black",0.1) )
                        # }
                        for ( i in 1:3 ) lines( dat$bio6_tminc_8110 , p_mu[[i]] , col=rangi2)
                        for ( i in 1:3 ) shade(   p_ci[[i]],dat$bio6_tminc_8110  , col=col.alpha(rangi2,0.2))
                        
                        title(text[k])
                }
        }        
        
        out.object <- list(mfit = mfit_1.1, d=d)
        saveRDS(out.object, file = paste(ofolder, "sdm-ordered-categorical.rds", sep=""))
        return(out.object)
}

# This model is a bit more interesting. This is a zero-inflated ordered categorical distribution.
# In this model, for every site, there is a probability p that a plant was there but was not observed.
# To do so, we use a binomial distribution to decide whether or not the plant was there, and a second ordered categorical to decide whether or not we observed it and the abundance.
# The aim of a zero-inflated model here is to avoid having to use pseudo-absences. 
ordered.categorical.stan.zeroinflated <- function(d = NULL,idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F, splitdata=F, show.plot=T, ofolder="../../results/models/"){
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }
        
        # Split data
        if(splitdata){
                training.data <- d$dataset$train==1
                test.data <- !(training.data)
        }else{
                training.data <- d$dataset$train>-1
                test.data <- training.data
        }
        
        # Prepare training data for stan model
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_2.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(obs))-1,
                        obs=obs,
                        bio=bio)
        
        # Set starting values for the parameters
        start_2.1 <- list(
                gamma = 0,
                alpha = 1:dat_2.1$L/dat_2.1$L-length(dat_2.1$L)*0.5,
                beta = rep(0,length(variables))
        )
        
        # Initialize data structure
        n_chains_2.1 <- 3
        init_2.1 <- list()
        for ( i in 1:n_chains_2.1 ) init_2.1[[i]] <- start_2.1
        
        # Should I calculate the likelihood?
        if(loglik){
                model_code=model2.1
        }else{
                model_code=model2.0     
        }
        
        # Run stan model
        mfit_2.1 <- stan ( model_code=model_code ,
                           data=dat_2.1 ,
                           chains=n_chains_2.1 ,
                           cores= n_chains_2.1 ,
                           warmup=1000, iter=4000,
                           init=init_2.1 , control = list(adapt_delta = 0.95))
        
        # Build link function to make predictions        
        link_2.1 <- function(dat, post){
                N <- dim(post$alpha)[1]
                L <- dim(post$alpha)[2]
                M <- nrow(dat)
                K <- ncol(dat)
                beta <- t(post$beta)
                
                p <- as.matrix(dat) %*% beta
                p <- lapply(1:L, function(x) inv.logit(t(as.matrix(post$alpha[,x])[col(p)]-p)))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[test.data]
        dat <- d$dataset[test.data,3:ncol(d$dataset)]
        ord <- order(d$dataset$bio5_tmaxw_8110)
        
        obs <- obs[ord]
        dat <- dat[ord,]
        
        # Sample from the posterior
        post <- extract.samples(mfit_2.1)
        
        if(show.plot){
                # Have a look at the estimates
                par(mfrow=c(1,1))
                plot(precis(mfit_2.1,depth=2))
                
                #Visualize plot
                par(mfrow=c(3,3))
                condi1 <- c(-1,-1,-1,0,0,0,1,1,1)
                condi2 <- c(-1,0,1,-1,0,1,-1,0,1)
                text <- c("bio12=-1 bio6=-1", "bio12=-1 bio6=0", "bio12=-1 bio6=1",
                          "bio12=0 bio6=-1", "bio12=0 bio6=0", "bio12=0 bio6=1",
                          "bio12=1 bio6=-1", "bio12=1 bio6=0", "bio12=1 bio6=1")
                
                for(k in 1:9){
                        dat_ <- dat
                        dat_$bio12_p_8110 <- dat_$bio12_p_8110*0+condi1[k]
                        dat_$bio6_tminc_8110 <- dat_$bio6_tminc_8110*0+condi2[k]
                        
                        # Make predictions with test data
                        p_post <- link_2.1(post = post, dat=dat_)
                        p_mu <- lapply(1:length(p_post), function(x) apply( p_post[[x]] , 2 , mean ))
                        p_ci <- lapply(1:length(p_post), function(x) apply( p_post[[x]] , 2 , PI ))
                        
                        # Calculate AUC and stuff
                        plot( NULL , type="n" , xlab="tmax" , ylab="probability" ,
                              xlim=c(min(dat$bio5_tmaxw_8110),max(dat$bio5_tmaxw_8110)) , ylim=c(0,1) , yaxp=c(0,1,2) )
                        
                        # for ( s in 1:50 ) {
                        #         for ( i in 1:3 ) lines( dat$bio6_tminc_8110 , p_post[[i]][s,] , col=col.alpha("black",0.1) )
                        # }
                        for ( i in 1:3 ) lines( dat$bio5_tmaxw_8110 , p_mu[[i]] , col=rangi2)
                        for ( i in 1:3 ) shade(   p_ci[[i]],dat$bio5_tmaxw_8110  , col=col.alpha(rangi2,0.2))
                        
                        title(text[k])
                }
        }
        
        out.object <- list(mfit = mfit_2.1, d=d)
        saveRDS(out.object, file = paste(ofolder, "sdm-ordered-categorical-zeroinflated.rds", sep=""))
        return(out.object)
}

# This model might also be interesting. It is an ordered categorical regression, where all "variables" are used as predictors using gaussian RBF.
ordered.categorical.stan.gaussian <- function(d = NULL,idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F, splitdata=F, show.plot=T, ofolder="../../results/models/"){
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }
        
        # Split data
        if(splitdata){
                training.data <- d$dataset$train==1
                test.data <- !(training.data)
        }else{
                training.data <- d$dataset$train>-1
                test.data <- training.data
        }
        
        # Prepare training data for stan model
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_3.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(obs))-1,
                        obs=obs,
                        bio=bio)
        
        # Set starting values for the parameters
        start_3.1 <- list(
                gamma = 0,
                alpha = 1:dat_2.1$L/dat_2.1$L-length(dat_2.1$L)*0.5,
                mu_bar = rep(0,length(variables)),
                beta = rep(0,length(variables)),
                epsilon = rep(1,length(variables))
        )
        
        # Initialize data structure
        n_chains_3.1 <- 3
        init_3.1 <- list()
        for ( i in 1:n_chains_3.1 ) init_3.1[[i]] <- start_3.1
        
        # Should I calculate the likelihood?
        if(loglik){
                model_code=model3.1
        }else{
                model_code=model3.0     
        }
        
        # Run stan model
        mfit_3.1 <- stan ( model_code=model_code ,
                           data=dat_3.1 ,
                           chains=n_chains_3.1 ,
                           cores= n_chains_3.1 ,
                           warmup=1000, iter=4000,
                           init=init_3.1 , control = list(adapt_delta = 0.95))
        
        # Build link function to make predictions 
        link_3.1 <- function(dat, post){
                L <- dim(post$alpha)[2]
                p <- 0
                for(k in 1:dim(dat)[2]){
                        a <- as.matrix(dat[,k])
                        b <- t(post$mu_bar[,k])
                        c <- t(post$epsilon[,k])
                        e <- t(post$beta[,k])
                        p <- p + e[col(a)] * exp(-0.5*((a[,row(b)] - b[col(a), ])/c[col(a), ])**2)
                }
                p <- lapply(1:L, function(x) inv.logit(t(as.matrix(post$alpha[,x])[col(p)]-p)))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[test.data]
        dat <- d$dataset[test.data,3:ncol(d$dataset)]
        ord <- order(d$dataset$bio5_tmaxw_8110)
        
        obs <- obs[ord]
        dat <- dat[ord,]
        
        # Sample from the posterior
        post <- extract.samples(mfit_3.1)
        
        if(show.plot){
                # Have a look at the estimates
                par(mfrow=c(1,1))
                plot(precis(mfit_3.1,depth=2))
                
                #Visualize plot
                par(mfrow=c(3,3))
                condi1 <- c(-1,-1,-1,0,0,0,1,1,1)
                condi2 <- c(-1,0,1,-1,0,1,-1,0,1)
                text <- c("bio12=-1 bio6=-1", "bio12=-1 bio6=0", "bio12=-1 bio6=1",
                          "bio12=0 bio6=-1", "bio12=0 bio6=0", "bio12=0 bio6=1",
                          "bio12=1 bio6=-1", "bio12=1 bio6=0", "bio12=1 bio6=1")
                
                for(k in 1:9){
                        dat_ <- dat
                        dat_$bio12_p_8110 <- dat_$bio12_p_8110*0+condi1[k]
                        dat_$bio6_tminc_8110 <- dat_$bio6_tminc_8110*0+condi2[k]
                        
                        # Make predictions with test data
                        p_post <- link_3.1(post = post, dat=dat_)
                        p_mu <- lapply(1:length(p_post), function(x) apply( p_post[[x]] , 2 , mean ))
                        p_ci <- lapply(1:length(p_post), function(x) apply( p_post[[x]] , 2 , PI ))
                        
                        # Calculate AUC and stuff
                        plot( NULL , type="n" , xlab="tmax" , ylab="probability" ,
                              xlim=c(min(dat$bio5_tmaxw_8110),max(dat$bio5_tmaxw_8110)) , ylim=c(0,1) , yaxp=c(0,1,2) )
                        
                        # for ( s in 1:50 ) {
                        #         for ( i in 1:3 ) lines( dat$bio6_tminc_8110 , p_post[[i]][s,] , col=col.alpha("black",0.1) )
                        # }
                        for ( i in 1:3 ) lines( dat$bio5_tmaxw_8110 , p_mu[[i]] , col=rangi2)
                        for ( i in 1:3 ) shade(   p_ci[[i]],dat$bio5_tmaxw_8110  , col=col.alpha(rangi2,0.2))
                        
                        title(text[k])
                }
        }
        
        out.object <- list(mfit = mfit_3.1, d=d)
        saveRDS(out.object, file = paste(ofolder, "sdm-ordered-categorical-gaussian.rds", sep=""))
        return(out.object)
}



