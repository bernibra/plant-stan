source("./prepare-data.R")
source("./models.R")
library(rethinking)
library(rstan)

####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression.
####

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=T, show.plot=T, pca=F, ndim=2, simulated=F, ofolder="../../results/models/"){
                
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "-simulated"
        }
        
        # Load the data
        if(is.null(d)){
                if(recompile){
                       d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="linear")
                       # rename variables
                       if(simulated){
                               variables <- c("S1", "S2")
                       }else{
                               if(pca){
                                       variables <- paste("PC", 1:ndim, sep="")
                               }
                       }
                       saveRDS(d, file = paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = ""))
                }else{
                       # rename variables
                        if(simulated){
                                variables <- c("S1", "S2")
                        }else{
                                if(pca){
                                        variables <- paste("PC", 1:ndim, sep="")
                                }
                        }
                        d <- readRDS(file = paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = ""))
                }
        }


        # Prepare training data for stan model
        obs <- d$obs
        id <- d$id
        bio <- d[,(ncol(d)-length(variables)+1):ncol(d)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_1.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(id)),
                        obs=obs,
                        bio=bio,
                        id=id)
        

        
        if (non.centered){
                # Set starting values for the parameters
                start_1.1 <- list(
                        zalpha = rep(0, dat_1.1$L),
                        alpha_bar = 0,
                        sigma_a = 0.1,
                        zbeta = matrix(0,dat_1.1$K,dat_1.1$L),
                        beta_bar = rep(0, dat_1.1$K),
                        sigma_b = rep(0.1, dat_1.1$K)
                )
                
                model_code=model1.2
                
        }else{
                # Set starting values for the parameters
                start_1.1 <- list(
                        zalpha = rep(0, dat_1.1$L),
                        alpha = 0,
                        sigma_a = 0.1,
                        beta = matrix(0,dat_1.1$L,dat_1.1$K),
                        zbeta = rep(0, dat_1.1$K),
                        sigma_b = rep(0.1, dat_1.1$K)
                )
                
                # Should I calculate the likelihood?
                if(loglik){
                        model_code=model1.1
                }else{
                        model_code=model1.0     
                }
        }
        
        # Initialize data structure
        n_chains_1.1 <- 3
        init_1.1 <- list()
        for ( i in 1:n_chains_1.1 ) init_1.1[[i]] <- start_1.1
        
        # Run stan model
        mfit_1.1 <- stan ( model_code=model_code ,
                           data=dat_1.1 ,
                           chains=n_chains_1.1 ,
                           cores= n_chains_1.1 ,
                           warmup=1000, iter=4000,
                           init=init_1.1 , control = list(adapt_delta = 0.95))
        
        saveRDS(mfit_1.1, file = paste(ofolder, "binomial-jsdm",extension,".rds", sep=""))
        return(mfit_1.1)
}

binomial.stan(recompile = F, pca = T, ndim = 3, ofolder="/cluster/scratch/bemora/plant-stan/")

