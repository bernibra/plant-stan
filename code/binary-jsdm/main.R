# source("./prepare-data.R")
source("./models.R")
# library("gridExtra")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression - with correlation between parameters.
# - Simple binomial regression - with covariance matrix for sp.
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
binomial.stan.gauss <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=T, show.plot=T, pca=F, ndim=2, simulated=F, ofolder="../../results/models/"){
        
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
                zbeta = matrix(0, dat_2.1$K, dat_2.1$L),
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
        

        saveRDS(mfit_2.1, file = paste(ofolder, "binomial-jsdm",extension,".rds", sep=""))
        return(mfit_2.1)
}

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan.gauss.RBFs <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=F, show.plot=T, pca=T, ndim=2, simulated=T, ofolder="../../results/models/"){
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "-simulated"
        }
        
        # Load the data
        if(is.null(d)){
                if(recompile){
                        d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="gauss.gauss")
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
        dat <- d$dataset
        obs <- dat$obs
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_3.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(id)),
                        obs=obs,
                        bio=bio,
                        id=id,
                        Dmat=Dis)
        
        # Set starting values for the parameters
        start_3.1 <- list(
                zalpha = rep(0, dat_3.1$L),
                alpha_bar = 0,
                sigma_a = 0.1,
                sigma_beta = matrix(0.1, dat_3.1$K, dat_3.1$L),
                zbeta = matrix(0, dat_3.1$K, dat_3.1$L),
                beta_bar = rep(0, dat_3.1$K),
                sigma_b = rep(0.1, dat_3.1$K),
                etasq = rep(0.1, dat_3.1$K),
                rhosq = rep(0.1, dat_3.1$K)
        )
        
        model_code=model3.0
        
        # Initialize data structure
        n_chains_3.1 <- 3
        init_3.1 <- list()
        for ( i in 1:n_chains_3.1 ) init_3.1[[i]] <- start_3.1
        
        # Run stan model
        mfit_3.1 <- stan ( model_code=model_code ,
                           data=dat_3.1 ,
                           chains=n_chains_3.1 ,
                           cores= n_chains_3.1 ,
                           warmup=1000, iter=4000,
                           init=init_3.1 , control = list(adapt_delta = 0.95))
        
        
        saveRDS(mfit_3.1, file = paste(ofolder, "binomial-stan-gauss-RBFs",extension,".rds", sep=""))
        return(mfit_3.1)
        
        
}

check_results_latest <- function(d, model){
        alphas <- precis(model, pars = "alpha", depth=2)
        betas <- precis(model, pars = "beta", depth=3)
        # sigmas <- precis(model, pars = "sigma_beta", depth=3)
        N <- length(unique(d$dataset$id))
        alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
        beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
        beta2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta2[1])
        d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , sd=c(rep(0,length(alpha_r)),alphas$sd))
        d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , sd=c(rep(0,length(beta1_r)),betas[1:N,]$sd))
        d_beta2 <- data.frame(N=1:N, id= c(rep("real", length(beta2_r)), rep("estimated", length(betas[(N+1):(N+N),]$mean))),value=c(beta2_r,betas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(beta2_r)),betas[(N+1):(N+N),]$sd))
        

        
        p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        p3 <- ggplot(d_beta2, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        figure <- grid.arrange(p1, p2, p3,
                            ncol = 1, nrow = 3)
}

binomial.stan.gauss.RBFs(recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")

