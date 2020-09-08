source("./sdm-maxent.R")
source("./models.R")
library(rethinking)
library(rstan)
library(gtools)

####
# Run species distribution model for a given species
####

# This first one is just a logistic regression, where we use the presence of species and pseudo-absences to estimate the probability of encountering a speces
# Those probabilities are not real probabilities, as the absences are not real absences.
binomial.stan <- function(idx=128, variables=c("bio5_", "bio6_","bio12_")){
        
        ## Run maxent model to make sure we are processing the same data as maxent
        # Alternatively, you could just load the presence data for species idx and the climatic data, and generate pseudo-absences. EASY
        d <- species_distribution.maxent(idx = idx, view_plots = F, variables=variables, pseudoA = F)

        # Prepare training data for stan model
        c <- d$dataset$train==1
        obs <- d$dataset$obs[c]
        bio <- d$dataset[c,3:ncol(d$dataset)]
        
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
        
        # Run stan model
        mfit_1.1 <- stan ( model_code=model1.0 ,
                           data=dat_1.1 ,
                           chains=n_chains_1.1 ,
                           cores= n_chains_1.1 ,
                           warmup=1000, iter=4000,
                           init=init_1.1 , control = list(adapt_delta = 0.95))

        priors <- function(bio, mean=0, sigma=1){
                N <- dim(bio)[1]
                K <- 100
                p <- matrix(0,K,N)
                for(j in 1:K){
                        for(i in 1:N){
                                p[j, i] <- inv.logit(rnorm(n = 1, mean = mean, sd = sigma) + sum(bio[i,] * rnorm(n = 3, mean = mean, sd = sigma)))
                        }
                }
                dens(as.vector(p))
        }

        # Build link function to make predictions        
        link_1.1 <- function(dat, post){
                N <- dim(post$alpha)[1]
                M <- nrow(dat)
                K <- ncol(dat)
                p <- matrix(0,N,M)
                for( i in 1:N){
                        for( j in 1:M){
                                p[i,j] <- post$alpha[i]
                                for( k in 1:K){
                                        p[i, j] <- p[i, j] + post$beta[i,k]*dat[j,k]
                                }
                                p[i,j] <- inv.logit(p[i,j])
                        }
                }
                return(p)
        }
        
        # Extract test data from the dataset
        obs <- d$dataset$obs[c]
        dat <- d$dataset[c,3:ncol(d$dataset)]
        
        # Sample from the posterior
        post <- extract.samples(mfit_1.1)
        
        # Make predictions with test data
        p_post <- link_1.1(post = post, dat=dat)
        p_mu <- apply( p_post , 2 , mean )
        p_ci <- apply( p_post , 2 , PI )
        
        plot(p_mu, obs)
        plot(precis(mfit_1.1,depth=2))
        
}

