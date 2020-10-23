source("./sdm-maxent.R")
source("./models.R")
library(rethinking)
library(rstan)
library(gtools)
library(pROC)

####
## Run species distribution model for a given species. I wrote the following models:
# - Simple binomial regression.
# - Binomial regression with quadratic terms.
# - Zero-inflated binomial regression.
# - Zero-inflated binomial regression with quadratic terms.
####

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan <- function(d = NULL,idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F){
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }

        # Prepare training data for stan model
        training.data <- d$dataset$train==1
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_1.1 <- list(N=length(obs),
                        K=length(variables),
                        obs=obs,
                        bio=bio)
        
        # Set starting values for the parameters
        start_1.1 <- list(
                alpha = 0,
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

# This second one is also just a logistic regression, but in this case we add fucking quadratic terms to the regression.
binomial.stan.quadratic <- function(d=NULL, idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F){
        if(pseudoA){
                warning("This will take quite some time to run.\n Unless you are running this on the cluster, I would think twice before wasting time on a model with stupid quadratic terms.")
        }
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }
        
        # Prepare training data for stan model
        training.data <- d$dataset$train==1
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        bio2 <- bio**2
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_2.1 <- list(N=length(obs),
                        K=length(variables),
                        obs=obs,
                        bio=bio,
                        bio2 = bio2)
        
        # Set starting values for the parameters
        start_2.1 <- list(
                alpha = 0,
                beta = rep(0,length(variables)),
                beta2 = rep(0,length(variables))
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
        
        priors <- function(bio, sigma1=1, sigma2=0.5){
                N <- dim(bio)[1]
                K <- 1000
                L <- dim(bio)[2]

                p <- as.matrix(bio) %*% matrix(rnorm(n = L*K, mean = mean, sd = sigma1), nrow = L, ncol = K)
                p <- p + as.matrix(bio**2) %*% matrix(rnorm(n = L*K, mean = mean, sd = sigma2), nrow = L, ncol = K)
                p <- inv.logit(t(p + rnorm(n=K, mean = mean, sd = sigma1)[col(p)]))

                dens(as.vector(p))
        }
        
        # Build link function to make predictions        
        link_2.1 <- function(dat, post){
                beta <- t(post$beta)
                beta2 <- t(post$beta2)
                
                p <- as.matrix(dat) %*% beta
                p <- p + as.matrix(dat**2) %*% beta2
                p <- inv.logit(t(p + as.matrix(post$alpha)[col(p)]))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[!(training.data)]
        dat <- d$dataset[!(training.data),3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_2.1)
        
        # Make predictions with test data
        p_post <- link_2.1(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Have a look at the different posteriors
        plot(precis(mfit_2.1,depth=2))
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc(roc_obj)
        
        # Compare maxent and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("logistic regression with quadratic terms")
        plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
        title("maxent model")
        return(list(mfit = mfit_2.1, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d))
}

# This model is a bit more interesting. This is a zero-inflated binommial distribution.
# In this model, for every site, there is a probability p that a plant was there but was not observed.
# To do so, we use a binomial distribution to decide whether or not the plant was there, and a second binomial to decide whether or not we observed it.
# The aim of a zero-inflated model here is to avoid having to use pseudo-absences. When adding quadratic terms, it seems
# to do the trick (see below). That said, some posterior distributions swing around quite a lot, which could mean that there is a fair bit of overfitting as well as
# correlations across predictor variables.
binomial.stan.zeroinflated <- function(d = NULL,idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F){
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }
        
        # Prepare training data for stan model
        training.data <- d$dataset$train==1
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_3.0 <- list(N=length(obs),
                        K=length(variables),
                        obs=obs,
                        bio=bio)
        
        # Set starting values for the parameters
        start_3.0 <- list(
                alpha1 = 0,
                alpha2 = 0,
                beta = rep(0,length(variables))
        )
        
        # Initialize data structure
        n_chains_3.0 <- 3
        init_3.0 <- list()
        for ( i in 1:n_chains_3.0 ) init_3.0[[i]] <- start_3.0
        
        # Should I calculate the likelihood?
        if(loglik){
                model_code=model3.1
        }else{
                model_code=model3.0     
        }
        
        # Run stan model
        mfit_3.0 <- stan ( model_code=model_code ,
                           data=dat_3.0 ,
                           chains=n_chains_3.0 ,
                           cores= n_chains_3.0 ,
                           warmup=1000, iter=4000,
                           init=init_3.0 , control = list(adapt_delta = 0.95))
        
        
        # Build link function to make predictions        
        link_3.0 <- function(dat, post){
                N <- dim(post$alpha1)[1]
                M <- nrow(dat)
                K <- ncol(dat)
                beta <- t(post$beta)
                
                p <- as.matrix(dat) %*% beta
                p <- inv.logit(t(p + as.matrix(post$alpha1)[col(p)])) * inv.logit(t(p*0 + as.matrix(post$alpha2)[col(p)]))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[!(training.data)]
        dat <- d$dataset[!(training.data),3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_3.0)
        
        # Make predictions with test data
        p_post <- link_3.0(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Have a look at the estimates
        plot(precis(mfit_3.0,depth=2))
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc_obj <- auc(roc_obj)
        
        # Compare maxent and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("zero-inflated logistic regression")
        plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
        title("maxent model")
        
        # Compare maxent with pseudo-absences and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("zero-inflated logistic regression")
        plot(d$roc_pseudoAmaxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_pseudoAmaxent))), cex = .8)
        title("maxent model with pseudo-absences")
        
        return(list(mfit = mfit_3.0, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d))
}

# This model is a bit more interesting. This is a zero-inflated binommial distribution with quadratic terms.
# In this model, for every site, there is a probability p that a plant was there but was not observed.
# To do so, we use a binomial distribution to decide whether or not the plant was there, and a second binomial to decide whether or not we observed it.
# The aim of a zero-inflated model here is to avoid having to use pseudo-absences. When adding quadratic terms, it seems
# to do the trick. That said, some posterior distributions swing around quite a lot, which could mean that there is a fair bit of overfitting as well as
# correlations across predictor variables.
binomial.stan.zeroinflated.quadratic <- function(d = NULL,idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F){
        if(pseudoA){
                warning("This will take quite some time to run.\n Unless you are running this on the cluster, I would think twice before wasting time on a model with stupid quadratic terms.")
        }
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }
        
        # Prepare training data for stan model
        training.data <- d$dataset$train==1
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]
        bio2 <- bio**2
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_4.0 <- list(N=length(obs),
                        K=length(variables),
                        obs=obs,
                        bio=bio,
                        bio2 = bio2)
        
        # Set starting values for the parameters
        start_4.0 <- list(
                alpha1 = 0,
                alpha2 = 0,
                beta = rep(0,length(variables)),
                beta2 = rep(0,length(variables))
        )
        
        # Initialize data structure
        n_chains_4.0 <- 3
        init_4.0 <- list()
        for ( i in 1:n_chains_4.0 ) init_4.0[[i]] <- start_4.0
        
        # Should I calculate the likelihood?
        if(loglik){
                model_code=model4.1
        }else{
                model_code=model4.0     
        }
        
        # Run stan model
        mfit_4.0 <- stan ( model_code=model_code ,
                           data=dat_4.0 ,
                           chains=n_chains_4.0 ,
                           cores= n_chains_4.0 ,
                           warmup=1000, iter=4000,
                           init=init_4.0 , control = list(adapt_delta = 0.95))
        
        
        # Build link function to make predictions        
        link_4.0 <- function(dat, post){
                beta <- t(post$beta)
                beta2 <- t(post$beta2)
                
                p <- as.matrix(dat) %*% beta
                p <- p + as.matrix(dat**2) %*% beta2
                p <- inv.logit(t(p + as.matrix(post$alpha1)[col(p)])) * inv.logit(t(p*0 + as.matrix(post$alpha2)[col(p)]))
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[!(training.data)]
        dat <- d$dataset[!(training.data),3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_4.0)
        
        # Make predictions with test data
        p_post <- link_4.0(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Have a look at the estimates
        plot(precis(mfit_4.0,depth=2))
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc_obj <- auc(roc_obj)
        
        # Compare maxent and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("zero-inflated logistic regression")
        plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
        title("maxent model")
        
        # Compare maxent with pseudo-absences and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("zero-inflated logistic regression")
        plot(d$roc_pseudoAmaxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_pseudoAmaxent))), cex = .8)
        title("maxent model with pseudo-absences")
        
        return(list(mfit = mfit_4.0, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d))
}


# This third one is also just a logistic regression, but in this case we add something that is much better than a
# quadratic term. Instead, we add a spline so that we find an optimal environment value and standard deviation from that mean value
# NOTE: while I like this option, it does not do very well. I am not sure why... but I like it for incorporating information about traits.
binomial.stan.gaussian <- function(d=NULL, idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F){
        if(pseudoA){
                warning("This will take quite some time to run.\n Unless you are running this on the cluster, I would think twice before wasting time on a model with stupid quadratic terms.")
        }
        
        # Load the data. This function will run maxent, but also process all the data for our stan model
        # The idea is that we can then compare the estimates from maxent and the estimates from our linear model
        # pseudoA defines whether or not we use pseudo absences or actual abasences.
        if(is.null(d)){
                d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = pseudoA)                
        }
        
        # Prepare training data for stan model
        training.data <- d$dataset$train==1
        obs <- d$dataset$obs[training.data]
        bio <- d$dataset[training.data,3:ncol(d$dataset)]

        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_5.0 <- list(N=length(obs),
                        K=length(variables),
                        obs=obs,
                        bio=bio)
        
        # Set starting values for the parameters
        start_5.0 <- list(
                alpha = 0,
                mu_bar = rep(0,length(variables)),
                beta = rep(0,length(variables)),
                epsilon = rep(1,length(variables))
        )
        
        # Initialize data structure
        n_chains_5.0 <- 3
        init_5.0 <- list()
        for ( i in 1:n_chains_5.0 ) init_5.0[[i]] <- start_5.0
        
        # Should I calculate the likelihood?
        if(loglik){
                model_code=model5.0
        }else{
                model_code=model5.0     
        }
        
        # Run stan model
        mfit_5.0 <- stan ( model_code=model_code ,
                           data=dat_5.0 ,
                           chains=n_chains_5.0 ,
                           cores= n_chains_5.0 ,
                           warmup=1000, iter=4000,
                           init=init_5.0 , control = list(adapt_delta = 0.95))
        
        priors <- function(bio, sigma1=1, sigma2=0.5){
                N <- dim(bio)[1]
                K <- 1000
                L <- dim(bio)[2]
                
                p <- as.matrix(bio) %*% matrix(rnorm(n = L*K, mean = mean, sd = sigma1), nrow = L, ncol = K)
                p <- p + as.matrix(bio**2) %*% matrix(rnorm(n = L*K, mean = mean, sd = sigma2), nrow = L, ncol = K)
                p <- inv.logit(t(p + rnorm(n=K, mean = mean, sd = sigma1)[col(p)]))
                
                dens(as.vector(p))
        }
        
        # Build link function to make predictions        
        link_5.0 <- function(dat, post){
                p <- 0
                for(k in 1:dim(dat)[2]){
                   a <- as.matrix(dat[,k])
                   b <- t(post$mu_bar[,k])
                   c <- t(post$epsilon[,k])
                   d <- t(post$beta[,k])
                   p <- p + d[col(a)] * exp(-0.5*((a[,row(b)] - b[col(a), ])/c[col(a), ])**2)
                }
                p <- inv.logit(t(p + as.matrix(post$alpha)[col(p)]))
                
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[!(training.data)]
        dat <- d$dataset[!(training.data),3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_5.0)
        
        # Make predictions with test data
        p_post <- link_5.0(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        # Have a look at the different posteriors
        plot(precis(mfit_5.0,depth=2))
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc(roc_obj)
        
        # Compare maxent and stan model
        par(mfrow=c(2,1))
        plot(roc_obj, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
        title("logistic regression with gaussian RBF")
        plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
        text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
        title("maxent model")
        return(list(mfit = mfit_5.0, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d))
}
