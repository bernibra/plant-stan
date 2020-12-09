source("./prepare-data.R")
source("./models.R")
library("gridExtra")
library(rethinking)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

####
## Run species distribution model for all species. I wrote the following models:
# - Binomial regression with gaussian RBFs.
####

binomial.stan.gauss.RBFs <- function(d = NULL, recompile = T, simulated=T, min.occurrence=10, ofolder="../../results/models/"){
        pca <- T
        ndim <- 2
        variables=c("bio5_", "bio6_","bio12_", "gdd5_", "bio1_","bio15_","bio17_", "bio8_", "TabsY_")
        gp_type <- 2
        
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
                        d <- species_distribution.data(variables=variables, pca=pca, ndim = ndim, simulated=simulated, min.occurrence=min.occurrence)
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
        Y <- d$Y
        X1 <- d$X1
        X2 <- d$X2
        dat <- d$dataset
        obs <- dat$obs
        id <- dat$id
        bio <- dat[,(ncol(dat)-length(variables)+1):ncol(dat)]
        
        dat_3.1 <- list(N=nrow(Y),
                        L=ncol(Y),
                        K=2,
                        Y=t(Y),
                        X1=X1,
                        X2=X2,
                        Dmat_b=Dis_b,
                        Dmat_g=Dis_g
                        )
        
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
        
        model_code=model3.beta

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

ptm <- proc.time()
model1 <- binomial.stan.gauss.RBFs(simulated=T, recompile = T, ofolder="~/Desktop/")
time.first <- proc.time() - ptm

