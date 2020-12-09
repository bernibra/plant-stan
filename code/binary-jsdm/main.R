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
                        # saveRDS(d, file = paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = ""))
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
                           init=init_4.1 , control = list(adapt_delta = 0.95))
        
        
        saveRDS(mfit_4.1, file = paste(ofolder, "binomial-stan-gauss-RBFs-beta",extension,".rds", sep=""))
        return(mfit_4.1)
}


check_results_latest <- function(d, model){
        alphas <- precis(model, pars = "alpha", depth=2)
        betas <- precis(model, pars = "beta", depth=3)
        sigmas <- precis(model, pars = "gamma", depth=3)
        
        # sigmas <- precis(model, pars = "sigma_beta", depth=3)
        N <- length(unique(d$dataset$id))
        alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
        beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
        beta2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta2[1])
        sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
        sigma2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta2[1])
        d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , sd=c(rep(0,length(alpha_r)),alphas$sd))
        d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , sd=c(rep(0,length(beta1_r)),betas[1:N,]$sd))
        d_beta2 <- data.frame(N=1:N, id= c(rep("real", length(beta2_r)), rep("estimated", length(betas[(N+1):(N+N),]$mean))),value=c(beta2_r,betas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(beta2_r)),betas[(N+1):(N+N),]$sd))
        d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , sd=c(rep(0,length(sigma1_r)),sigmas[1:N,]$sd))
        d_sigma2 <- data.frame(N=1:N, id= c(rep("real", length(sigma2_r)), rep("estimated", length(sigmas[(N+1):(N+N),]$mean))),value=c(sigma2_r,sigmas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(sigma2_r)),sigmas[(N+1):(N+N),]$sd))
        
        
        p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        p3 <- ggplot(d_beta2, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        p4 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        p5 <- ggplot(d_sigma2, aes(x=N, y=value, group=id, color=id)) + 
                geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
        figure <- grid.arrange(p1, p2, p3,
                            ncol = 1, nrow = 3)
        print(figure)
        figure <- grid.arrange(p4, p5,
                               ncol = 1, nrow = 2)
        print(figure)
}

# model1 <- binomial.stan.gauss.RBFs(simulated=T, recompile = T, gp_type = 2, ofolder="~/Desktop/")
binomial.stan.gauss.RBFs.beta(simulated=T, recompile = F, ofolder="/cluster/scratch/bemora/plant-stan/")

