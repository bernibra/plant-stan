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

baseline.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
        # Fixing some of the options
        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
        gp_type <- 2
        ndim <- 2
        pca <- T
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "1d-simulated"
        }
        extension <- paste(extension, as.character(gp_type), sep="")
        
        if(min.occurrence==10){
                extension2 <- ""
        }else{
                extension2 <- "min30-"
        }
        
        
        # Load the data
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
                if(is.null(d)){
                        d <- readRDS(file = paste("../../data/processed/jsdm/", extension, paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
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

        dat_5.1 <- list(N=N,
                        L=L,
                        Y=t(obs),
                        X1=X1,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g
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
                rhosq_b = 0.1,
                sigma_g = 0.1,
                etasq_g = 0.1,
                rhosq_g = 0.1
        )
        
        model_code=base.model.1d
        
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


        saveRDS(mfit_5.1, file = paste(ofolder, extension2, "baseline-model-1d", extension,".rds", sep=""))
        return(mfit_5.1)
}

categorical.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=30, ofolder="../../results/models/"){
        # Fixing some of the options
        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
        gp_type <- 2
        ndim <- 2
        pca <- T
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "1d-simulated"
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
                        d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="categorical", min.occurrence=min.occurrence)
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
                        d <- readRDS(file = paste("../../data/processed/jsdm/", extension, paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
                }
        }
        
        # Prepare training data for stan model
        Dis_b <- d$corr
        Dis_g <- d$corr2
        N <- sum(d$dataset$id==1)
        L <- length(unique(d$dataset$id))
        if(simulated){
                obs <- d$dataset$obs
        }else{
                obs <- d$dataset$abundance
        }
        dat <- d$dataset
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        X1 <- matrix(bio[,1], N, L)[,1]
        
        dat_5.1 <- list(N=N,
                        L=L,
                        K=length(obs),
                        M=length(unique(obs))-1,
                        Y=as.numeric(as.character(obs)),
                        X1=X1,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g
        )
        
        # Set starting values for the parameters
        start_5.1 <- list(
                phi = (1:dat_5.1$M)/dat_5.1$M-0.5,
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
                rhosq_b = 0.1,
                sigma_g = 0.1,
                etasq_g = 0.1,
                rhosq_g = 0.1
        )
        
        model_code=categorical.model.1d
        
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
        
        
        saveRDS(mfit_5.1, file = paste(ofolder, extension2, "categorical-model-traits-1d", extension,".rds", sep=""))
        return(mfit_5.1)
}

skew.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
        # Fixing some of the options
        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
        gp_type <- 2
        ndim <- 2
        pca <- T
        
        # File name extension
        extension <- ""
        
        # If we are dealing with simulated data
        if(simulated){
                extension <- "skew-simulated"
        }
        extension <- paste(extension, as.character(gp_type), sep="")
        
        if(min.occurrence==10){
                extension2 <- ""
        }else{
                extension2 <- paste("min", as.character(min.occurrence),"-",sep="")
        }
        
        
        # Load the data
        if(recompile){
                d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="skew", min.occurrence=min.occurrence)
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
                if(is.null(d)){
                        d <- readRDS(file = paste("../../data/processed/jsdm/", extension, paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
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
        
        dat_5.1 <- list(N=N,
                        L=L,
                        Y=t(obs),
                        X1=X1,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g
        )
        
        # Set starting values for the parameters
        start_5.1 <- list(
                zlambda = rep(0, dat_5.1$L),
                zalpha = rep(0, dat_5.1$L),
                zbeta = rep(0, dat_5.1$L),
                zgamma = rep(0, dat_5.1$L),
                lambda_bar = 0,
                alpha_bar = -1,
                beta_bar = 0,
                gamma_bar = 0,
                sigma_l = 0.1,
                sigma_a = 0.1,
                sigma_b = 0.1,
                etasq_b = 0.1,
                rhosq_b = 0.1,
                sigma_g = 0.1,
                etasq_g = 0.1,
                rhosq_g = 0.1
        )
        
        model_code=skew.model.traits.1d
        
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
                           init=init_5.1 , control = list(adapt_delta=0.95, max_treedepth = 15))
        
        
        saveRDS(mfit_5.1, file = paste(ofolder, extension2, "skew-model-traits-1d", extension,".rds", sep=""))
        return(mfit_5.1)
}

# model1 <- binomial.stan.gauss.RBFs(simulated=T, recompile = T, gp_type = 2, ofolder="~/Desktop/")
# baseline.1d(d=NULL, simulated=F, recompile = F, min.occurrence = 30, ofolder="/cluster/scratch/bemora/plant-stan/")

# d <- readRDS(file = "../../data/processed/jsdm/2PC1PC2min50-data.rds")
d <- readRDS(file = "../../data/processed/jsdm/skew-simulated2S1S2data.rds")
skew.1d(d=d, simulated=T, recompile = F, min.occurrence = 10, ofolder="/cluster/scratch/bemora/plant-stan/")
# # baseline.1d(d=d, simulated=T, recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")
# # 
# expose_stan_functions(model_r)
# 
# skew.1d(d=NULL, simulated=T, recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")
#
#
# skn <- function(x, alpha, sigma, beta, lambda){
# 
# 
#         delta <- lambda/sqrt(1+lambda**2)
#         # skewness <- 0.5*(4-pi)*(delta*sqrt(2/pi))**3/(1-2*delta**2/pi)**(3/2)
#         mu_z <- sqrt(2/pi)*delta
#         # sigma_z <- sqrt(1-mu_z**2)
#         # mode_x <- beta + 1/sqrt(2*sigma)*(mu_z- skewness*sigma_z*0.5-0.5*sign(lambda)*exp(-2*pi/abs(lambda)))
#         # maxy_ <- dsn(mode_x, xi=beta, omega=sqrt(1/(2*sigma)), alpha=lambda)
# 
#         maxy = 0.5 * ( 4 - pi ) * (delta * sqrt(2/pi))**3 / (1 - 2 * delta**2 / pi )**(3 / 2.0);
# 
#         maxy = beta + 1 / sqrt( 2 * sigma) * (mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5 - 0.5 * sign(lambda) * exp(- 2 * pi / abs(lambda) ))
#         print(c(maxy, 1 / sqrt( 2 * sigma), mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5,
#                                              0.5 * sign(lambda) * exp(- 2 * pi / abs(lambda) ),
#                 (mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5 - 0.5 * sign(lambda) * exp(- 2 * pi / abs(lambda) ))))
#         
#         maxy = exp(- sigma * (maxy - beta)**2) * (1 + pracma::erf((lambda * (maxy - beta)) * sqrt(sigma) ))
#         print(c(maxy,findmax(delta, beta, sigma, lambda)))
#         y <- exp(-log(findmax(delta, beta, sigma, lambda)+0.0001) - alpha - sigma * (beta - x)**2) * (1 + pracma::erf(lambda * (x-beta) * sqrt(sigma)))
#         # y <- exp(-alpha - sigma * (beta - x)**2) / (1 + lambda * (x-beta) * sqrt(2) * sqrt(sigma))
#         y
# }
# 
# 
# x <- seq(-10, 10, length.out = 2000)
# alpha=0
# lambda=15
# sigma_beta1=0.001
# beta1=-0
# lambda_hat <- lambda/sqrt(1+lambda**2)
# sigma_hat <- sigma_beta1 * (1 - (2*(lambda_hat**2))/pi)
# alpha_hat <- 1
# beta_hat <- beta1 - sqrt(1/(2*sigma_hat)) * lambda_hat * sqrt(2/pi)
# y <- skn(x, alpha=alpha, beta=beta_hat, sigma=sigma_hat, lambda=lambda)
# y_r <- rsn(n=6000, xi=beta_hat, omega=sqrt(1/(2*sigma_hat)), alpha=lambda)
# plot(x, y, type="l", ylim=c(0,1))
# abline(v=beta_hat, col="red")
# abline(v=mean(y_r), col="blue")
# lines(c(beta1-sqrt(1/(2*sigma_beta1)), beta1+sqrt(1/(2*sigma_beta1))), c(max(y)*0.5,max(y)*0.5), col="blue")
# lines(c(beta_hat-sqrt(1/(2*sigma_hat)), beta_hat+sqrt(1/(2*sigma_hat))), c(max(y)*0.4,max(y)*0.4), col="red")
# lines(c(beta1-sd(y_r), beta1+sd(y_r)), c(max(y)*0.45,max(y)*0.45), col="black")
# print(max(y))
