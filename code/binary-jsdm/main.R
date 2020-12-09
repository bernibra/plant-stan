# source("./prepare-data.R")
source("./models.R")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression - with correlation between parameters.
# - Simple binomial regression - with covariance matrix for sp.
# - Binomial regression with gaussian RBFs.
# - Binomial regression with gaussian RBFs - faster version
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
                       filename <- paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = "")
                       if (file.exists(filename)){
                               question <- askYesNo("Do you want to overwrite the file?", default = F, 
                                                    prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
                               if(is.na(question)){question <- F}
                               if(question){
                                       saveRDS(d, file=filename)
                               }else{
                                       stop("you should add 'recompile=T'")
                               }                                
                       }else{
                               saveRDS(d, file=filename)
                       }
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
                        filename <- paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = "")
                        if (file.exists(filename)){
                                question <- askYesNo("Do you want to overwrite the file?", default = F, 
                                                     prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
                                if(is.na(question)){question <- F}
                                if(question){
                                        saveRDS(d, file=filename)
                                }else{
                                        stop("you should add 'recompile=T'")
                                }                                
                        }else{
                                saveRDS(d, file=filename)
                        }
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
binomial.stan.gauss.RBFs <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=F, show.plot=T, pca=T, ndim=2, simulated=F, gp_type=2, min.occurrence=10, ofolder="../../results/models/"){
        
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
                        filename <- paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = "")
                        if (file.exists(filename)){
                                question <- askYesNo("Do you want to overwrite the file?", default = F, 
                                                     prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
                                if(is.na(question)){question <- F}
                                if(question){
                                        saveRDS(d, file=filename)
                                }else{
                                        stop("you should add 'recompile=T'")
                                }                                
                        }else{
                                saveRDS(d, file=filename)
                        }
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
                                bio=t(bio),
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
                
                model_code=model3.2
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

binomial.stan.gauss.RBFs.beta <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
        # Fixing some of the options
        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
        gp_type <- 2
        ndim <- 2
        pca <- T
        
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
                        filename <- paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = "")

                        if (file.exists(filename)){
                                question <- askYesNo("Do you want to overwrite the file?", default = F, 
                                                     prompts = getOption("askYesNo", gettext(c("Yes", "No", "Cancel"))))
                                if(is.na(question)){question <- F}
                                if(question){
                                        saveRDS(d, file=filename)
                                }else{
                                        stop("you should add 'recompile=T'")
                                }                                
                        }else{
                                saveRDS(d, file=filename)
                        }
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
        N <- sum(d$dataset$id==1)
        L <- length(unique(d$dataset$id))
        obs <- matrix(d$dataset$obs, N, L)
        dat <- d$dataset
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        X1 <- matrix(bio[,1], N, L)[,1]
        X2 <- matrix(bio[,2], N, L)[,1]
        
        dat_4.1 <- list(N=N,
                        L=L,
                        K=2,
                        Y=t(obs),
                        X1=X1,
                        X2=X2,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g
        )
        
        # Set starting values for the parameters
        start_4.1 <- list(
                zalpha = rep(0, dat_4.1$L),
                zbeta = matrix(0, dat_4.1$K, dat_4.1$L),
                zgamma = matrix(0, dat_4.1$K, dat_4.1$L),
                alpha_bar = 0,
                beta_bar = rep(0, dat_4.1$K),
                gamma_bar = rep(0, dat_4.1$K),
                sigma_a = 0.1,
                sigma_b = rep(0.1, dat_4.1$K),
                sigma_g = rep(0.1, dat_4.1$K),
                etasq_b = rep(0.1, dat_4.1$K),
                rhosq_b = rep(0.1, dat_4.1$K),
                etasq_g = rep(0.1, dat_4.1$K),
                rhosq_g = rep(0.1, dat_4.1$K)
        )
        
        model_code=model4.1
        
        # Initialize data structure
        n_chains_4.1 <- 3
        init_4.1 <- list()
        for ( i in 1:n_chains_4.1 ) init_4.1[[i]] <- start_4.1
        
        # Run stan model
        mfit_4.1 <- stan ( model_code=model_code ,
                           data=dat_4.1 ,
                           chains=n_chains_4.1 ,
                           cores= n_chains_4.1 ,
                           warmup=1000, iter=2000,
                           init=init_4.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
        
        
        saveRDS(mfit_4.1, file = paste(ofolder, "binomial-stan-gauss-RBFs-beta",extension,".rds", sep=""))
        return(mfit_4.1)
}

# model1 <- binomial.stan.gauss.RBFs(simulated=T, recompile = T, gp_type = 2, ofolder="~/Desktop/")
binomial.stan.gauss.RBFs.beta(simulated=F, recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")

