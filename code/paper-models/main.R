# source("./prepare-data.R")
source("./models.R")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - baseline model
####

baseline <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
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
        
        model_code=base.model
        
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
        
        
        saveRDS(mfit_4.1, file = paste(ofolder, "baseline-model", extension,".rds", sep=""))
        return(mfit_4.1)
}

baseline.traits <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
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
        Dis_t <- d$corr3
        N <- sum(d$dataset$id==1)
        L <- length(unique(d$dataset$id))
        obs <- matrix(d$dataset$obs, N, L)
        dat <- d$dataset
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        X1 <- matrix(bio[,1], N, L)[,1]
        X2 <- matrix(bio[,2], N, L)[,1]
        
        dat_5.1 <- list(N=N,
                        L=L,
                        K=2,
                        Y=t(obs),
                        X1=X1,
                        X2=X2,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g,
                        Dmat_t=Dis_t
        )
        
        # Set starting values for the parameters
        start_5.1 <- list(
                zalpha = rep(0, dat_5.1$L),
                zbeta = matrix(0, dat_5.1$K, dat_5.1$L),
                zgamma = matrix(0, dat_5.1$K, dat_5.1$L),
                alpha_bar = 0,
                beta_bar = rep(0, dat_5.1$K),
                gamma_bar = rep(0, dat_5.1$K),
                sigma_a = 0.1,
                etasq_a = 0.1,
                rhosq_a = 0.1,
                sigma_b = rep(0.1, dat_5.1$K),
                etasq_b = rep(0.1, dat_5.1$K),
                rhosq_b = matrix(0.1, dat_5.1$K, 2),
                sigma_g = rep(0.1, dat_5.1$K),
                etasq_g = rep(0.1, dat_5.1$K),
                rhosq_g = matrix(0.1, dat_5.1$K, 2)
        )
        
        model_code=base.model.traits
        
        # Initialize data structure
        n_chains_5.1 <- 3
        init_5.1 <- list()
        for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
        
        # Run stan model
        mfit_5.1 <- stan ( model_code=model_code ,
                           data=dat_5.1 ,
                           chains=n_chains_5.1 ,
                           cores= n_chains_5.1 ,
                           warmup=1000, iter=2000,
                           init=init_5.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
        
        
        saveRDS(mfit_5.1, file = paste(ofolder, "baseline-model-traits", extension,".rds", sep=""))
        return(mfit_5.1)
}

baseline.traits.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
        # Fixing some of the options
        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
        gp_type <- 2
        ndim <- 2
        pca <- T
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "-1d-simulated"
        }
        extension <- paste(extension, as.character(gp_type), sep="")
        
        if(min.occurrence==10){
                extension2 <- ""
        }else{
                extension2 <- "min30-"
        }
        
        
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
                        filename <- paste("../../data/processed/jsdm/", extension, paste(variables, collapse = ""), extension2, "data.rds", sep = "")
                        
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
        Dis_t <- d$corr3
        N <- sum(d$dataset$id==1)
        L <- length(unique(d$dataset$id))
        obs <- d$dataset$obs
        dat <- d$dataset
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        X1 <- matrix(bio[,1], N, L)[,1]

        dat_5.1 <- list(N=N,
                        L=L,
                        K=length(obs),
                        Y=obs,
                        X1=X1,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g,
                        Dmat_t=Dis_t
        )
        
        # Set starting values for the parameters
        start_5.1 <- list(
                zalpha = rep(0, dat_5.1$L),
                zbeta = rep(0, dat_5.1$L),
                zgamma = rep(0, dat_5.1$L),
                alpha_bar = 0,
                beta_bar = 0,
                gamma_bar = 0,
                sigma_a = 0.1,
                etasq_a = 0.1,
                rhosq_a = 0.1,
                sigma_b = 0.1,
                etasq_b = 0.1,
                rhosq_b = rep(0.1, 2),
                sigma_g = 0.1,
                etasq_g = 0.1,
                rhosq_g = rep(0.1, 2)
        )
        
        model_code=base.model.traits.1d.lent
        
        # Initialize data structure
        n_chains_5.1 <- 3
        init_5.1 <- list()
        for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
        
        # Run stan model
        mfit_5.1 <- stan ( model_code=model_code ,
                           data=dat_5.1 ,
                           chains=n_chains_5.1 ,
                           cores= n_chains_5.1 ,
                           warmup=1000, iter=2000,
                           init=init_5.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
        
        
        saveRDS(mfit_5.1, file = paste(ofolder, extension2, "baseline-model-traits-1d", extension,".rds", sep=""))
        return(mfit_5.1)
}

# model1 <- binomial.stan.gauss.RBFs(simulated=T, recompile = T, gp_type = 2, ofolder="~/Desktop/")
baseline.traits.1d(simulated=T, recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")

