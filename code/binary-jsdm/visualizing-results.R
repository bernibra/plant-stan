# Import libraries
library(rethinking)
library(rstan)
library("gridExtra")
library(ggpubr)

####
# Visualizing the results of some of the models
####

plot.simulated.data <- function(beta=T, gp_type = 2){
  # load data
  d <- readRDS(file = paste("../../data/processed/jsdm/S1S2", "data.rds", sep = ""))
  if(beta){
    model <- readRDS(paste("../../results/models/binomial-stan-gauss-RBFs-beta-simulated",as.character(gp_type),".rds", sep=""))
  }else{
    model <- readRDS(paste("../../results/models/binomial-stan-gauss-RBFs-simulated",as.character(gp_type),".rds", sep=""))
  }

  # Extract variables from the model
  alphas <- precis(model, pars = "alpha", depth=2)
  betas <- precis(model, pars = "beta", depth=3)
  sigmas <- precis(model, pars = "gamma", depth=3)
  
  # Extract true values
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  beta2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta2[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
  sigma2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta2[1])
  
  # Build data.frames for the plots
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , sd=c(rep(0,length(alpha_r)),alphas$sd))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , sd=c(rep(0,length(beta1_r)),betas[1:N,]$sd))
  d_beta2 <- data.frame(N=1:N, id= c(rep("real", length(beta2_r)), rep("estimated", length(betas[(N+1):(N+N),]$mean))),value=c(beta2_r,betas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(beta2_r)),betas[(N+1):(N+N),]$sd))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , sd=c(rep(0,length(sigma1_r)),sigmas[1:N,]$sd))
  d_sigma2 <- data.frame(N=1:N, id= c(rep("real", length(sigma2_r)), rep("estimated", length(sigmas[(N+1):(N+N),]$mean))),value=c(sigma2_r,sigmas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(sigma2_r)),sigmas[(N+1):(N+N),]$sd))
  
  # Generate plot
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.title = element_blank())
  leg <- get_legend(p1)
  p1 <-  p1 + theme(legend.position = "none")
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p3 <- ggplot(d_beta2, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 2") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p4 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p5 <- ggplot(d_sigma2, aes(x=N, y=value, group=id, color=id)) + ggtitle("gamma 2") +
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none") 
  figure <- grid.arrange(p2, p3, p4, p5, p1,leg,
                         ncol = 2, nrow = 3)
  print(figure)
  return(figure)
}


