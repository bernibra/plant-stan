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

baseline.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  ndim <- 2
  pca <- T
  
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
    d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="gauss.gauss", min.occurrence=min.occurrence)
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
  
  model_code=base.model.generror.1d
  
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

  saveRDS(mfit_5.1, file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  return(mfit_5.1)
}

skew.generror.1d <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
  # Fixing some of the options
  variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
  ndim <- 2
  pca <- T
  
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
    d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, simulated.type="skew.generror", min.occurrence=min.occurrence)
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
  
  dat_5.1 <- list(N=N,
                  L=L,
                  minp=1e-100,
                  Y=t(obs),
                  indices=1:L,
                  X1=X1,
                  Dmat_b=Dis_b,
                  Dmat_g=Dis_g
  )
  model_code = base.model.skew.generror.1d.multithread
  generror1d <- cmdstan_model(write_stan_file(model_code), cpp_options = list(stan_threads = TRUE))
  mfit_5.1 <- generror1d$sample(data = dat_5.1,
                                init = init_5.1,
                                chains = 3,
                                threads_per_chain = 15,
                                parallel_chains = 3,
                                max_treedepth = 15,
                                adapt_delta = 0.95,
                                refresh = 500)
  mfit_5.1$save_object(file = paste(ofolder, extension2, "", extension,"-cdmstan.rds", sep=""))
  saveRDS(rstan::read_stan_csv(mfit_5.1$output_files()), file = paste(ofolder, extension2, "", extension,".rds", sep=""))
  
  return(mfit_5.1)
}


min.occurrence <- 20
# if(min.occurrence==10){
#   d <- readRDS(file = paste("../../data/processed/jsdm/1d-PC1PC2-data.rds", sep=""))
# }else{
#   d <- readRDS(file = paste("../../data/processed/jsdm/1d-PC1PC2min",min.occurrence,"-data.rds", sep=""))
# }
d <- readRDS("../../data/processed/jsdm/1d-PC1PC2min20-data.rds")
skew.generror.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")
# baseline.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")
# skew.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")
# generror.1d(d=d, simulated=F, recompile = F, min.occurrence = min.occurrence, ofolder="/cluster/scratch/bemora/plant-stan/")





####################################### OTHER STUFF ###############################################
####################################### OTHER STUFF ###############################################
####################################### OTHER STUFF ###############################################

# # 
# expose_stan_functions(mod)
# # 
# # skew.1d(d=NULL, simulated=T, recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")
# #
# #
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
#         # print(c(maxy, 1 / sqrt( 2 * sigma), mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5,
#         #                                      0.5 * sign(lambda) * exp(- 2 * pi / abs(lambda) ),
#         #         (mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5 - 0.5 * sign(lambda) * exp(- 2 * pi / abs(lambda) ))))
# 
#         maxy = exp(- sigma * (maxy - beta)**2) * (1 + pracma::erf((lambda * (maxy - beta)) * sqrt(sigma) ))
#         # print(c(maxy,findmax(delta, beta, sigma, lambda)))
#         y <- exp(-log(findmax(delta, beta, sigma, lambda)) - alpha - sigma * (beta - x)**2) * (1 + pracma::erf(lambda * (x-beta) * sqrt(sigma)))
#         # y <- exp(-alpha - sigma * (beta - x)**2) / (1 + lambda * (x-beta) * sqrt(2) * sqrt(sigma))
#         y
# }
# 
# 
# x <- seq(-5, 5, length.out = 2000)
# alpha=0
# lambda=10
# sigma_beta1=6
# beta1=1
# lambda_hat <- lambda/sqrt(1+lambda**2)
# sigma_hat <- sigma_beta1 * (1 - (2*(lambda_hat**2))/pi)
# beta_hat <- beta1 - sqrt(1/(2*sigma_hat)) * lambda_hat * sqrt(2/pi)
# alpha_hat <- log(findmax(lambda_hat, beta_hat, sigma_hat, lambda))+exp(alpha)
# 
# y <- skn(x, alpha=alpha, beta=beta_hat, sigma=sigma_hat, lambda=lambda)
# y_r <- rsn(n=6000, xi=beta_hat, omega=sqrt(1/(2*sigma_hat)), alpha=lambda)
# plot(x, y, type="l", ylim=c(0,1))
# abline(v=beta_hat, col="red")
# abline(v=mean(y_r), col="blue")
# lines(c(beta1-sqrt(1/(2*sigma_beta1)), beta1+sqrt(1/(2*sigma_beta1))), c(max(y)*0.5,max(y)*0.5), col="blue")
# lines(c(beta_hat-sqrt(1/(2*sigma_hat)), beta_hat+sqrt(1/(2*sigma_hat))), c(max(y)*0.4,max(y)*0.4), col="red")
# lines(c(beta1-sd(y_r), beta1+sd(y_r)), c(max(y)*0.45,max(y)*0.45), col="black")
# print(max(y))
# 

# generrskew <- function(x, a, mu, sigma, lambda, p){
#   v <- sqrt((pi*gamma(1/p))/(pi*(1+3*lambda**2)*gamma(3/p)-(16**(1/p))*lambda*lambda*(gamma(1/2+1/p)**2)*gamma(1/p)))
#   m <- (2**(2/p))*v*sigma*lambda*gamma(1/2+1/p)/sqrt(pi)
#   a*exp(-(abs(x-mu+m)/(v*sigma*(1+lambda*sign(x-mu+m))))**p)
# }
