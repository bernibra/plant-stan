source("./prepare-data.R")
source("./models.R")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


categorical.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  gp_type <- 2
  ndim <- 2
  pca <- T
  
  # File name extension
  extension <- ""
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "simulated-categorical"
  }
  extension <- paste(extension, as.character(gp_type), sep="")
  
  if(min.occurrence==10){
    extension2 <- ""
  }else{
    extension2 <- "min30-"
  }
  
  
  # Load the data
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
    if(is.null(d)){
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
  obs <- as.numeric(as.character(obs))
  obs <- matrix(obs, N, L)
  dat <- d$dataset
  id <- dat$id
  bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
  X1 <- matrix(bio[,1], N, L)[,1]
  
  dat_5.1 <- list(N=N,
                  L=L,
                  M=length(unique(as.vector(obs)))-1,
                  Y=obs,
                  X1=X1,
                  minp=0.0000000001,
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
    sigma_b = 0.1,
    etasq_b = 0.1,
    rhosq_b = 0.1,
    sigma_g = 0.1,
    etasq_g = 0.1,
    rhosq_g = 0.1
  )
  
  model_code=categorical.model.1d
  
  # Initialize data structure
  n_chains_5.1 <- 15
  init_5.1 <- list()
  for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
  
  # Run stan model
  mfit_5.1 <- stan ( model_code=model_code ,
                     data=dat_5.1 ,
                     chains=n_chains_5.1 ,
                     cores= n_chains_5.1 ,
                     warmup=1000, iter=1200,
                     init=init_5.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
  
  
  saveRDS(mfit_5.1, file = paste(ofolder, extension2, "categorical-model-1d", extension,".rds", sep=""))
  return(mfit_5.1)
}

categorical.skew.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  gp_type <- 2
  ndim <- 2
  pca <- T
  
  # File name extension
  extension <- ""
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "simulated-categorical"
  }
  extension <- paste(extension, as.character(gp_type), sep="")
  
  if(min.occurrence==10){
    extension2 <- ""
  }else{
    extension2 <- "min30-"
  }
  
  
  # Load the data
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
    if(is.null(d)){
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
  obs <- as.numeric(as.character(obs))
  obs <- matrix(obs, N, L)
  dat <- d$dataset
  id <- dat$id
  bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
  X1 <- matrix(bio[,1], N, L)[,1]
  
  dat_5.1 <- list(N=N,
                  L=L,
                  M=length(unique(as.vector(obs)))-1,
                  Y=obs,
                  X1=X1,
                  minp=0.0000000001,
                  Dmat_b=Dis_b,
                  Dmat_g=Dis_g
  )
  
  # Set starting values for the parameters
  start_5.1 <- list(
    phi = (1:dat_5.1$M)/dat_5.1$M-0.5,
    zlambda = rep(0, dat_5.1$L),
    zalpha = rep(0, dat_5.1$L),
    zbeta = rep(0, dat_5.1$L),
    zgamma = rep(0, dat_5.1$L),
    lambda_bar = 0,
    alpha_bar = 0,
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
  
  model_code=categorical.model.skew.1d
  
  # Initialize data structure
  n_chains_5.1 <- 15
  init_5.1 <- list()
  for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
  
  # Run stan model
  mfit_5.1 <- stan ( model_code=model_code ,
                     data=dat_5.1 ,
                     chains=n_chains_5.1 ,
                     cores= n_chains_5.1 ,
                     warmup=1000, iter=1200,
                     init=init_5.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
  
  
  saveRDS(mfit_5.1, file = paste(ofolder, extension2, "categorical-model-skew-1d", extension,".rds", sep=""))
  return(mfit_5.1)
}

d <- readRDS(file = "../../data/processed/jsdm/1d-simulated2S1S2categorical-data.rds")
categorical.1d(d=d, simulated=T, recompile=F, min.occurrence = 10, ofolder="/cluster/scratch/bemora/plant-stan/")
# categorical.skew.1d(d=d, simulated=T, recompile=F, min.occurrence = 10, ofolder="/cluster/scratch/bemora/plant-stan/")
