source("./prepare-data.R")
source("./models.R")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression - with correlation between parameters.
# - Simple binomial regression - with covariance matrix for sp.
# - Simple binomial regression - with  correlation between parameters AND covariance matrix for sp.
# - Binomial regression with gaussian RBFs.
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
                       d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="linear.corr")
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
        d <- dataset$dataset
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

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan.corr <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=T, show.plot=T, pca=F, ndim=2, simulated=F, ofolder="../../results/models/"){
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "-simulated"
        }
        
        # Load the data
        if(is.null(d)){
                if(recompile){
                        d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="linear.gauss")
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
        Dis <- d$corr
        d <- d$dataset
        obs <- d$obs
        id <- d$id
        bio <- d[,(ncol(d)-length(variables)+1):ncol(d)]
        
        # dat_list <- list(
        #         obs = obs,
        #         L=length(unique(id)),
        #         id = id,
        #         bio1=bio$S1,
        #         bio2=bio$S2,
        #         Dmat=Dis)
        # 
        # m14.7 <- ulam(
        #         alist(
        #                 obs ~ dbinom(1, p),
        #                 logit(p) <- alpha[id]*sigma_a + alpha_bar + beta1[id] * bio1 + beta2[id] * bio2,
        #                 transpars> vector[L]: beta1 <<- L_SIGMA1 * zbeta1,
        #                 transpars> vector[L]: beta2 <<- L_SIGMA2 * zbeta2,
        #                 vector[L]:zbeta1 ~ normal( 0 , 1 ),
        #                 vector[L]:zbeta2 ~ normal( 0 , 1 ),
        #                 transpars> matrix[L,L]: L_SIGMA1 <<- cholesky_decompose( SIGMA1 ),
        #                 transpars> matrix[L,L]:SIGMA1 <- cov_GPL2( Dmat , etasq1 , rhosq1 , si1 ),
        #                 transpars> matrix[L,L]: L_SIGMA2 <<- cholesky_decompose( SIGMA2 ),
        #                 transpars> matrix[L,L]:SIGMA2 <- cov_GPL2( Dmat , etasq2 , rhosq2 , si2 ),
        #                 alpha[id]~dnorm(0,1),
        #                 alpha_bar ~ dnorm(0,1.3),
        #                 c(etasq1,etasq2) ~ dexp( 1 ),
        #                 c(rhosq1,rhosq2) ~ dexp( 0.5 ),
        #                 c(si1,si2) ~ dexp(1),
        #                 sigma_a ~ dexp( 1 )
        # ), data=dat_list , chains=3 , cores=3 , iter=3000 )
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_2.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(id)),
                        obs=obs,
                        bio=bio,
                        id=id,
                        Dmat=Dis)

        # Set starting values for the parameters
        start_2.1 <- list(
                zalpha = rep(0, dat_2.1$L),
                alpha_bar = 0,
                sigma_a = 0.1,
                beta1 = rep(0, dat_2.1$L),
                beta2 = rep(0, dat_2.1$L),
                beta_bar = rep(0, dat_2.1$K),
                sigma_b = rep(0.1, dat_2.1$K),
                etasq = rep(0.1, dat_2.1$K),
                rhosq = rep(0.1, dat_2.1$K)
        )

        model_code=model2.0

        # Initialize data structure
        n_chains_2.1 <- 3
        init_2.1 <- list()
        for ( i in 1:n_chains_2.1 ) init_2.1[[i]] <- start_2.1

        # Run stan model
        mfit_2.1 <- stan ( model_code=model_code ,
                           data=dat_2.1 ,
                           chains=n_chains_2.1 ,
                           cores= n_chains_2.1 ,
                           warmup=1000, iter=4000,
                           init=init_2.1 , control = list(adapt_delta = 0.95))
        

        saveRDS(mfit_1.1, file = paste(ofolder, "binomial-jsdm",extension,".rds", sep=""))
        return(mfit_1.1)
}


binomial.stan(recompile = F, pca = T, ndim = 3, ofolder="/cluster/scratch/bemora/plant-stan/")

