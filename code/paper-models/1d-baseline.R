# source("./prepare-data.R")
source("./models-binomial.R")
library(rethinking)
library(rstan)
library(cmdstanr)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - baseline model
####

baseline.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, elevation=F, temperature=F, ofolder="../../results/models/"){

  if (temperature & elevation){
    stop("only one of 'temperature' and 'elevation' can be true")
  }else if( temperature ){
    # Fixing some of the options
    variables=c("bio1_")
    ndim <- 1
    pca <- F
  }else if( elevation ){
    # Fixing some of the options
    variables=c("elevation_")
    ndim <- 1
    pca <- F
  }else{
    # Fixing some of the options
    variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
    ndim <- 2
    pca <- T
  }
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "1d-simulated"
  }else{
    extension <- "1d"
  }
  
  if(min.occurrence==10){
    extension2 <- "-"
  }else{
    extension2 <- paste("min",min.occurrence,"-",sep="")
  }
  
  # Load the data
  if(recompile){
    d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, elevation = elevation, simulated=simulated, simulated.type="gauss.gauss", min.occurrence=min.occurrence)
    # rename variables
    if(simulated){
      variables <- c("S1", "S2")
    }else{
      if(pca){
        variables <- paste("PC", 1:ndim, sep="")
      }
    }
    filename <- paste("../../data/processed/jsdm/", extension, "-",paste(variables, collapse = ""), extension2, "data.rds", sep = "")
    
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
      d <- readRDS(file = paste("../../data/processed/jsdm/", extension,"-", paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
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
  if(elevation | temperature){
    X1 <- matrix(bio, N, L)[,1]
    if(elevation){
      extension <- gsub("1d","elevation",extension)
    }else{
      extension <- gsub("1d","temperature",extension)
    }
  }else{
    X1 <- matrix(bio[,1], N, L)[,1]
  }
  
  dat_5.1 <- list(N=N,
                  L=L,
                  minp=1e-100,
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


  saveRDS(mfit_5.1, file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  return(mfit_5.1)
}

generror.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  ndim <- 2
  pca <- T
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "1d-generror-simulated"
  }else{
    extension <- "1d-generror"
  }
  
  if(min.occurrence==10){
    extension2 <- "-"
  }else{
    extension2 <- paste("min",min.occurrence,"-",sep="")
  }
  
  # Load the data
  if(recompile){
    d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="generror", min.occurrence=min.occurrence)
    # rename variables
    if(simulated){
      variables <- c("S1", "S2")
    }else{
      if(pca){
        variables <- paste("PC", 1:ndim, sep="")
      }
    }
    filename <- paste("../../data/processed/jsdm/", extension, "-",paste(variables, collapse = ""), extension2, "data.rds", sep = "")
    
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
      d <- readRDS(file = paste("../../data/processed/jsdm/", extension,"-", paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
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
                  minp=1e-100,
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
    znu = rep(0, dat_5.1$L),
    alpha_bar = 0,
    nu_bar = 0,
    beta_bar = 0,
    gamma_bar = 0,
    sigma_a = 0.1,
    sigma_n = 0.1,
    sigma_b = 0.1,
    etasq_b = 0.1,
    rhosq_b = 0.1,
    sigma_g = 0.1,
    etasq_g = 0.1,
    rhosq_g = 0.1
  )
  
  # Initialize data structure
  n_chains_5.1 <- 3
  init_5.1 <- list()
  for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
  
  # model_code=base.model.generror.1d
  # 
  # # Run stan model
  # mfit_5.1 <- stan ( model_code=model_code ,
  #                    data=dat_5.1 ,
  #                    chains=n_chains_5.1 ,
  #                    cores= n_chains_5.1 ,
  #                    warmup=1000, iter=2000,
  #                    init=init_5.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
  # 
  # saveRDS(mfit_5.1, file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  # 
  
  dat_5.1$indices <- 1:L
  
  model_code = base.model.generror.1d.multithread
  generror1d <- cmdstan_model(write_stan_file(model_code), cpp_options = list(stan_threads = TRUE))
  mfit_5.1 <- generror1d$sample(data = dat_5.1,
                                init = init_5.1,
                                chains = 3,
                                threads_per_chain = 15,
                                parallel_chains = 3,
                                max_treedepth = 15,
                                max_depth = 15,
                                iter_sampling = 1000,
                                #adapt_delta = 0.95,
                                refresh = 500)
  mfit_5.1$save_object(file = paste(ofolder, extension2, "", extension,"-cdmstan.rds", sep=""))
  saveRDS(rstan::read_stan_csv(mfit_5.1$output_files()), file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  print(mfit_5.1$cmdstan_diagnose())
  
  return(mfit_5.1)
}


skew.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  ndim <- 2
  pca <- T
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "1d-skew-simulated"
  }else{
    extension <- "1d-skew"
  }
  
  if(min.occurrence==10){
    extension2 <- "-"
  }else{
    extension2 <- paste("min",min.occurrence,"-",sep="")
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
    filename <- paste("../../data/processed/jsdm/", extension,"-", paste(variables, collapse = ""), extension2, "data.rds", sep = "")
    
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
      d <- readRDS(file = paste("../../data/processed/jsdm/","-", extension, "-", paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
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
                  minp=1e-100,
                  Y=t(obs),
                  X1=X1,
                  Dmat_b=Dis_b,
                  Dmat_g=Dis_g
  )
  
  # Set starting values for the parameters
  start_5.1 <- list(
    lambda = rep(0, dat_5.1$L),
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
  
  # Initialize data structure
  n_chains_5.1 <- 3
  init_5.1 <- list()
  for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
  
  # model_code=skew.model.traits.1d
  # 
  # # Run stan model
  # mfit_5.1 <- stan ( model_code=model_code ,
  #                    data=dat_5.1 ,
  #                    chains=n_chains_5.1 ,
  #                    cores= n_chains_5.1 ,
  #                    warmup=1000, iter=2000,
  #                    init=init_5.1 , control = list(adapt_delta=0.95, max_treedepth = 15))
  # 
  # saveRDS(mfit_5.1, file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  
  dat_5.1$indices <- 1:L
  
  model_code = base.model.skew.1d.multithread
  generror1d <- cmdstan_model(write_stan_file(model_code), cpp_options = list(stan_threads = TRUE))
  mfit_5.1 <- generror1d$sample(data = dat_5.1,
                                init = init_5.1,
                                chains = 3,
                                threads_per_chain = 15,
                                parallel_chains = 3,
                                max_treedepth = 15,
                                max_depth = 15,
                                iter_sampling = 1000,
                                #adapt_delta = 0.95,
                                refresh = 500)
  mfit_5.1$save_object(file = paste(ofolder, extension2, "", extension,"-cdmstan.rds", sep=""))
  saveRDS(rstan::read_stan_csv(mfit_5.1$output_files()), file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  print(mfit_5.1$cmdstan_diagnose())
  
  return(mfit_5.1)
}

skew.traits.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  ndim <- 2
  pca <- T
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "1d-skew-traits-simulated"
  }else{
    extension <- "1d-skew-traits"
  }
  
  if(min.occurrence==10){
    extension2 <- "-"
  }else{
    extension2 <- paste("min",min.occurrence,"-",sep="")
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
    filename <- paste("../../data/processed/jsdm/", extension,"-", paste(variables, collapse = ""), extension2, "data.rds", sep = "")
    
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
      d <- readRDS(file = paste("../../data/processed/jsdm/","-", extension, "-", paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
    }
  }
  
  # Prepare training data for stan model
  Dis_b <- d$corr
  Dis_g <- d$corr2
  Dis_l <- d$corr
  N <- sum(d$dataset$id==1)
  L <- length(unique(d$dataset$id))
  obs <- matrix(d$dataset$obs, N, L)
  dat <- d$dataset
  id <- dat$id
  bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
  X1 <- matrix(bio[,1], N, L)[,1]
  
  dat_5.1 <- list(N=N,
                  L=L,
                  minp=1e-100,
                  Y=t(obs),
                  X1=X1,
                  Dmat_b=Dis_b,
                  Dmat_g=Dis_g,
                  Dmat_l=Dis_l
  )
  
  # Set starting values for the parameters
  start_5.1 <- list(
    lambda = rep(0, dat_5.1$L),
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
    rhosq_g = 0.1,
    etasq_l = 0.1,
    rhosq_l = 0.1
  )
  
  # Initialize data structure
  n_chains_5.1 <- 3
  init_5.1 <- list()
  for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1
  
  # model_code=skew.model.traits.1d
  # 
  # # Run stan model
  # mfit_5.1 <- stan ( model_code=model_code ,
  #                    data=dat_5.1 ,
  #                    chains=n_chains_5.1 ,
  #                    cores= n_chains_5.1 ,
  #                    warmup=1000, iter=2000,
  #                    init=init_5.1 , control = list(adapt_delta=0.95, max_treedepth = 15))
  # 
  # saveRDS(mfit_5.1, file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  
  dat_5.1$indices <- 1:L
  
  model_code = base.model.skew.traits.1d.multithread
  generror1d <- cmdstan_model(write_stan_file(model_code), cpp_options = list(stan_threads = TRUE))
  mfit_5.1 <- generror1d$sample(data = dat_5.1,
                                init = init_5.1,
                                chains = 3,
                                threads_per_chain = 15,
                                parallel_chains = 3,
                                max_treedepth = 15,
                                max_depth = 15,
                                iter_sampling = 1000,
                                #adapt_delta = 0.95,
                                refresh = 500)
  mfit_5.1$save_object(file = paste(ofolder, extension2, "", extension,"-cdmstan.rds", sep=""))
  saveRDS(rstan::read_stan_csv(mfit_5.1$output_files()), file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  print(mfit_5.1$cmdstan_diagnose())
  
  return(mfit_5.1)
}

skew.generror.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, elevation=F, temperature=F, ofolder="../../results/models/"){

  if (temperature & elevation){
    stop("only one of 'temperature' and 'elevation' can be true")
  }else if( temperature ){
    # Fixing some of the options
    variables=c("bio1_")
    ndim <- 1
    pca <- F
  }else if( elevation ){
    # Fixing some of the options
    variables=c("elevation_")
    ndim <- 1
    pca <- F
  }else{
    # Fixing some of the options
    variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
    ndim <- 2
    pca <- T
  }
  
  # If we are dealing with simulated data
  if(simulated){
    extension <- "1d-skew-generror-simulated"
  }else{
    extension <- "1d-skew-generror"
  }    
  
  if(min.occurrence==10){
    extension2 <- "-"
  }else{
    extension2 <- paste("min",min.occurrence,"-",sep="")
  }
  
  # Load the data
  if(recompile){
    d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, elevation = elevation, simulated.type="skew.generror", min.occurrence=min.occurrence)
    # rename variables
    if(simulated){
      variables <- c("S1", "S2")
    }else{
      if(pca){
        variables <- paste("PC", 1:ndim, sep="")
      }
    }
    filename <- paste("../../data/processed/jsdm/", extension, "-",paste(variables, collapse = ""), extension2, "data.rds", sep = "")
    
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
      d <- readRDS(file = paste("../../data/processed/jsdm/", extension,"-", paste(variables, collapse = ""), extension2, "data.rds", sep = ""))
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
  if(elevation | temperature){
    X1 <- matrix(bio, N, L)[,1]
    if(elevation){
      extension <- gsub("1d","elevation",extension)
    }else{
      extension <- gsub("1d","temperature",extension)
    }
  }else{
    X1 <- matrix(bio[,1], N, L)[,1]
  }
  
  dat_5.1 <- list(N=N,
                  L=L,
                  minp=1e-100,
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
    zlambda = rep(0, dat_5.1$L),
    znu = rep(0, dat_5.1$L),
    alpha_bar = 0,
    nu_bar = 0,
    lambda_bar = 0,
    beta_bar = 0,
    gamma_bar = 0,
    sigma_a = 0.1,
    sigma_n = 0.1,
    sigma_l = 0.1,
    sigma_b = 0.1,
    etasq_b = 0.1,
    rhosq_b = 0.1,
    sigma_g = 0.1,
    etasq_g = 0.1,
    rhosq_g = 0.1
  )

  # Initialize data structure
  n_chains_5.1 <- 3
  init_5.1 <- list()
  for ( i in 1:n_chains_5.1 ) init_5.1[[i]] <- start_5.1

  # model_code=base.model.skew.generror.1d
  # 
  # # Run stan model
  # mfit_5.1 <- stan ( model_code=model_code ,
  #                    data=dat_5.1 ,
  #                    chains=n_chains_5.1 ,
  #                    cores= n_chains_5.1 ,
  #                    warmup=1000, iter=2000,
  #                    init=init_5.1 , control = list(adapt_delta = 0.95, max_treedepth = 15))
  # 
  # saveRDS(mfit_5.1, file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  
  dat_5.1$indices <- 1:L
  
  model_code = base.model.skew.generror.1d.multithread
  generror1d <- cmdstan_model(write_stan_file(model_code), cpp_options = list(stan_threads = TRUE))
  mfit_5.1 <- generror1d$sample(data = dat_5.1,
                                init = init_5.1,
                                chains = 3,
                                threads_per_chain = 15,
                                parallel_chains = 3,
                                max_treedepth = 15,
                                max_depth = 15,
                                iter_sampling = 1000,
                                #adapt_delta = 0.95,
                                refresh = 500)
  mfit_5.1$save_object(file = paste(ofolder, extension2, "", extension,"-cdmstan.rds", sep=""))
  saveRDS(rstan::read_stan_csv(mfit_5.1$output_files()), file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  print(mfit_5.1$cmdstan_diagnose())
  
  return(mfit_5.1)
}


min.occurrence <- 20

# Basic analyses
d <- readRDS("../../data/processed/jsdm/1d-PC1PC2min20-data.rds")
# d <- readRDS("../../data/processed/jsdm/skew-generror-simulated-data.rds")
# d <- readRDS("../../data/processed/jsdm/skew-simulated-data.rds")
# d <- readRDS("../../data/processed/jsdm/generror-simulated-data.rds")

# skew.generror.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")
# skew.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")
# generror.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")
# baseline.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")

# Traits to explain lambda
# skew.traits.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")

# Elevation
# d <- readRDS("../../data/processed/jsdm/1d-skew-generror-elevation_min20-data.rds")
# skew.generror.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, elevation=T, temperature=F, ofolder="/cluster/scratch/bemora/plant-stan/")
# baseline.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, elevation=T, temperature=F, ofolder="/cluster/scratch/bemora/plant-stan/")

# Temperature
# d <- readRDS("../../data/processed/jsdm/1d-skew-generror-bio1_min20-data.rds")
# skew.generror.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, elevation=F, temperature=T, ofolder="/cluster/scratch/bemora/plant-stan/")
# baseline.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, elevation=F, temperature=T, ofolder="/cluster/scratch/bemora/plant-stan/")

