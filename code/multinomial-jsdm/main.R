# source("./prepare-data.R")
source("./models.R")
# library("gridExtra")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - Ordered categorical regression with gaussian RBFs to fit categorical data for many species.
####

ordered.categorical.stan.gauss.RBFs <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=F, show.plot=T, pca=T, ndim=2, simulated=F, gp_type=2, min.occurrence=10, ofolder="../../results/models/"){
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "-simulated"
        }
        extension <- paste(extension, as.character(gp_type), sep="")
        
        # Load the data
        if(is.null(d)){
                if(recompile){
                        d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="gauss.gauss", min.occurrence=min.occurrence)
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
        Dis_b <- d$corr
        Dis_g <- d$corr2
        dat <- d$dataset
        obs <- dat$obs
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        if (gp_type==1){
                dat_3.1 <- list(N=length(obs),
                                K=length(variables),
                                L=length(unique(id)),
                                obs=obs,
                                bio=bio,
                                id=id,
                                Dmat=Dis_b)
                
                # Set starting values for the parameters
                start_3.1 <- list(
                        zalpha = rep(0, dat_3.1$L),
                        zbeta = matrix(0, dat_3.1$K, dat_3.1$L),
                        zgamma = matrix(0, dat_3.1$K, dat_3.1$L),
                        alpha_bar = 0,
                        beta_bar = rep(0, dat_3.1$K),
                        gamma_bar = rep(0, dat_3.1$K),
                        sigma_a = 0.1,
                        sigma_b = rep(0.1, dat_3.1$K),
                        sigma_g = rep(0.1, dat_3.1$K),
                        etasq_b = rep(0.1, dat_3.1$K),
                        rhosq_b = rep(0.1, dat_3.1$K)
                )
                
                model_code=model3.0
        }else{
                dat_3.1 <- list(N=length(obs),
                                K=length(variables),
                                L=length(unique(id)),
                                obs=obs,
                                bio=bio,
                                id=id,
                                Dmat_b=Dis_b,
                                Dmat_g=Dis_g)
                
                # Set starting values for the parameters
                start_3.1 <- list(
                        zalpha = rep(0, dat_3.1$L),
                        zbeta = matrix(0, dat_3.1$K, dat_3.1$L),
                        zgamma = matrix(0, dat_3.1$K, dat_3.1$L),
                        alpha_bar = 0,
                        beta_bar = rep(0, dat_3.1$K),
                        gamma_bar = rep(0, dat_3.1$K),
                        sigma_a = 0.1,
                        sigma_b = rep(0.1, dat_3.1$K),
                        sigma_g = rep(0.1, dat_3.1$K),
                        etasq_b = rep(0.1, dat_3.1$K),
                        rhosq_b = rep(0.1, dat_3.1$K),
                        etasq_g = rep(0.1, dat_3.1$K),
                        rhosq_g = rep(0.1, dat_3.1$K)
                )
                
                model_code=model3.1
                
        }
        
        # Initialize data structure
        n_chains_3.1 <- 3
        init_3.1 <- list()
        for ( i in 1:n_chains_3.1 ) init_3.1[[i]] <- start_3.1
        
        # Run stan model
        mfit_3.1 <- stan ( model_code=model_code ,
                           data=dat_3.1 ,
                           chains=n_chains_3.1 ,
                           cores= n_chains_3.1 ,
                           warmup=1000, iter=2000,
                           init=init_3.1 , control = list(adapt_delta = 0.95))
        
        
        saveRDS(mfit_3.1, file = paste(ofolder, "binomial-stan-gauss-RBFs",extension,".rds", sep=""))
        return(mfit_3.1)
        
        
}

binomial.stan.gauss.RBFs(recompile = F, gp_type = 2, ofolder="/cluster/scratch/bemora/plant-stan/")

