source("./sdm-maxent.R")
source("./models.R")
library(rethinking)
library(rstan)
library(gtools)
library(pROC)

####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression.
####

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan <- function(d = NULL, idx=128, variables=c("bio5_", "bio6_","bio12_"), pseudoA=F, loglik=F, show.plot=T, ofolder="../../results/models/"){
        
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
        
        # Calculate AUC and stuff
        roc_obj <- roc(obs, p_mu)
        auc_obj <- auc(roc_obj)
        
        if (show.plot){
                # Have a look at the estimates
                plot(precis(mfit_1.1,depth=2))
                
                # Compare maxent and stan model
                par(mfrow=c(2,1))
                plot(roc_obj, xlim = c(1, 0), asp=NA)
                text(0.2, 0.2, paste("AUC =", as.character(auc(roc_obj))), cex = .8)
                title("regular logistic regression")
                plot(d$roc_maxent, xlim = c(1, 0), asp=NA)
                text(0.2, 0.2, paste("AUC =", as.character(auc(d$roc_maxent))), cex = .8)
                title("maxent model")
        }
        
        out.object <- list(mfit = mfit_1.1, roc_obj = roc_obj, roc_maxent = d$roc_maxent, d=d)
        saveRDS(out.object, file = paste(ofolder, "sdm-binomial.rds", sep=""))
        return(out.object)
}
