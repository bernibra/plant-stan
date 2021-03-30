# Import libraries
library(rethinking)
library(rstan)
library(gridExtra)
library(grid)
library(ggpubr)
library(hrbrthemes)
library(plotly)
library(cowplot)

####
# Visualizing the results of some of the models
####

transform.skew <- function(post){
  delta <- post$lambda/sqrt(1+post$lambda**2)
  beta_hat <- post$beta + sqrt(1/(2*post$gamma)) * delta * sqrt(2/pi)
  sigma_hat <- post$gamma * (1 - (2*(delta**2))/pi)^(-1)
  
  mu_z <- sqrt(2/pi)*delta
  maxy = 0.5 * ( 4 - pi ) * (delta * sqrt(2/pi))**3 / (1 - 2 * delta**2 / pi )**(3 / 2.0);
  maxy = post$beta + 1 / sqrt( 2 * sigma_hat) * (mu_z - maxy * sqrt(1 - mu_z**2 ) * 0.5 - 0.5 * sign(post$lambda) * exp(- 2 * pi / abs(post$lambda) ))
  maxy = exp(- post$gamma * (maxy - post$beta)**2) * (1 + pracma::erf((post$lambda * (maxy - post$beta)) * sqrt(post$gamma) ))
  
  alpha_hat <- -log((maxy+0.0001)) + post$alpha
  return(list(beta=beta_hat, sigma=sigma_hat, alpha=alpha_hat, lambda=post$lambda))
}

plot.simulated.compare <- function(){
  m1 <- readRDS(paste("../../../results/models/skew-model-traits-1dskew-simulated2.rds", sep=""))
  m2 <- readRDS(paste("../../../results/models/baseline-model-1d1d-simulated2.rds", sep=""))
  rethinking::compare(m1, m2)
  
  d <- readRDS(file = paste("../../../data/processed/jsdm/skew-simulated2S1S2", "data.rds", sep = ""))
  
  for (i in 1:2){
    if (i==1){
      model_r <- m1
      post <- extract.samples(model_r, pars=c("alpha", "beta", "gamma", "lambda"))
      
      post <- transform.skew(post)
      
      # Extract variables from the model
      alphas <- precis(model_r, pars = "alpha", depth=2, )
      betas <- precis(model_r, pars = "beta", depth=3)
      sigmas <- precis(model_r, pars = "gamma", depth=3)
      lambdas <- precis(model_r, pars = "lambda", depth=3)
      
      betas$mean <- sapply(1:dim(post$beta)[2], function(x) mean(post$beta[,x]), USE.NAMES = F)
      betas$`5.5%` <- sapply(1:dim(post$beta)[2], function(x) PI(post$beta[,x], prob = c(0.890))[1], USE.NAMES = F)
      betas$`94.5%` <- sapply(1:dim(post$beta)[2], function(x) PI(post$beta[,x], prob = c(0.890))[2], USE.NAMES = F)
      sigmas$mean <- sapply(1:dim(post$sigma)[2], function(x) mean(post$sigma[,x]), USE.NAMES = F)
      sigmas$`5.5%` <- sapply(1:dim(post$sigma)[2], function(x) PI(post$sigma[,x], prob = c(0.890))[1], USE.NAMES = F)
      sigmas$`94.5%` <- sapply(1:dim(post$sigma)[2], function(x) PI(post$sigma[,x], prob = c(0.890))[2], USE.NAMES = F)
      alphas$mean <- sapply(1:dim(post$alpha)[2], function(x) mean(post$alpha[,x]), USE.NAMES = F)
      alphas$`5.5%` <- sapply(1:dim(post$alpha)[2], function(x) PI(post$alpha[,x], prob = c(0.890))[1], USE.NAMES = F)
      alphas$`94.5%` <- sapply(1:dim(post$alpha)[2], function(x) PI(post$alpha[,x], prob = c(0.890))[2], USE.NAMES = F)
      
    }else{
      model_r <- m2
      
      # Extract variables from the model
      alphas <- precis(model_r, pars = "alpha", depth=2, prob = 0.99)
      betas <- precis(model_r, pars = "beta", depth=3, prob = 0.99)
      sigmas <- precis(model_r, pars = "gamma", depth=3, prob = 0.99)
    }
    
    
    # Extract true values
    N <- length(unique(d$dataset$id))
    alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
    beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
    sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
    sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
    
    # Build data.frames for the plots
    d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , ymin=c(alpha_r,alphas[,3]) , ymax=c(alpha_r,alphas[,4]))
    d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , ymin=c(beta1_r,betas[,3]), ymax=c(beta1_r,betas[,4]))
    d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , ymin=c(sigma1_r,sigmas[,3]), ymax=c(sigma1_r,sigmas[,4]))
    
    # Generate plot
    p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
      geom_pointrange(aes(ymin=ymin, ymax=ymax)) + theme_linedraw() + theme(legend.title = element_blank())
    if(i==1){
      lambda_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$lambda[1])
      d_lambda <- data.frame(N=1:N, id= c(rep("real", length(lambda_r)), rep("estimated", length(lambdas$mean))),value=c(lambda_r,lambdas$mean) , ymin=c(lambda_r,lambdas[,3]), ymax=c(lambda_r,lambdas[,4]))
      p4 <- ggplot(d_lambda, aes(x=N, y=value, group=id, color=id))  + ggtitle("lambda") + 
        geom_pointrange(aes(ymin=ymin, ymax=ymax)) + theme_linedraw() + theme(legend.position = "none")
    }else{
      p4 <- get_legend(p1)
      
    }
    p1 <-  p1 + theme(legend.position = "none")
    p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
      geom_pointrange(aes(ymin=ymin, ymax=ymax)) + theme_linedraw() + theme(legend.position = "none")
    p3 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
      geom_pointrange(aes(ymin=ymin, ymax=ymax)) + theme_linedraw() + theme(legend.position = "none")
    
    figure <- grid.arrange(p1, p2, p3, p4,
                           ncol = 2, nrow = 2)
    print(figure)
  }
  
}

plot.simulated.data <- function(beta=T, gp_type = 2){
  # load data
  d <- readRDS("../../../data/processed/jsdm/skew-simulated-data.rds")
  model_r <- readRDS(paste("../../../results/models/skew-simulated.rds", sep=""))
  
  # Extract variables from the model
  alphas <- precis(model_r, pars = "alpha_", depth=2, prob = 0.89)
  betas <- precis(model_r, pars = "beta_", depth=3, prob = 0.89)
  sigmas <- precis(model_r, pars = "gamma_", depth=3, prob = 0.89)
  lambdas <- precis(model_r, pars = "lambda", depth=3, prob = 0.89)
  precis(model_r, pars = "lambda_bar", depth=3)
  precis(model_r, pars = "sigma_l", depth=3)
  
  
  # Extract true values
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  lambda_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$lambda[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
  
  # Build data.frames for the plots
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , lower=c(alpha_r,alphas[,3]), upper=c(alpha_r,alphas[,4]))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , lower=c(beta1_r,betas[,3]), upper=c(beta1_r,betas[,4]))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , lower=c(sigma1_r,sigmas[,3]), upper=c(sigma1_r,sigmas[,4]))
  d_lambda <- data.frame(N=1:N, id= c(rep("real", length(lambda_r)), rep("estimated", length(lambdas$mean))),value=c(lambda_r,lambdas$mean) ,  lower=c(lambda_r,lambdas[,3]), upper=c(lambda_r,lambdas[,4]))
  
  
  # Generate plot
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.title = element_blank())
  leg <- get_legend(p1)
  p1 <-  p1 + theme(legend.position = "none")
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p3 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p4 <- ggplot(d_lambda, aes(x=N, y=value, group=id, color=id))  + ggtitle("lambda") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  
  
  figure <- grid.arrange(p1, p2, p3, p4,
                         ncol = 2, nrow = 2)
  print(figure)
  return(figure)
}

plot.simulated.data.tails <- function(beta=T, gp_type = 2){
  # load data
  d <- readRDS(file = paste("../../../data/processed/jsdm/generror-simulated-", "data.rds", sep = ""))
  model_r <- readRDS(paste("../../../results/models/-1d-generror-simulated.rds", sep=""))
  
  # Extract variables from the model
  alphas <- precis(model_r, pars = "alpha", depth=2)
  betas <- precis(model_r, pars = "beta", depth=3)
  sigmas <- precis(model_r, pars = "gamma", depth=3)
  nu <- precis(model_r, pars = "nu", depth=3)
  
  post <- extract.samples(model_r, n = 1000, pars=c("nu", "gamma", "beta", "alpha", "beta_bar"))
  post$gamma <- (post$gamma**2) * gamma(1/post$nu) / gamma(3/post$nu)
  
  # PI(exp(post$nu_bar + 0.5*post$sigma_n**2)+1)
  
  p <- 0.89
  alphas$mean <- sapply(1:dim(post$alpha)[2], function(x) median(post$alpha[,x]), USE.NAMES = F)
  alphas$`5.5%` <- sapply(1:dim(post$alpha)[2], function(x) PI(post$alpha[,x], prob = c(p))[1], USE.NAMES = F)
  alphas$`94.5%` <- sapply(1:dim(post$alpha)[2], function(x) PI(post$alpha[,x], prob = c(p))[2], USE.NAMES = F)
  betas$mean <- sapply(1:dim(post$beta)[2], function(x) median(post$beta[,x]), USE.NAMES = F)
  betas$`5.5%` <- sapply(1:dim(post$beta)[2], function(x) PI(post$beta[,x], prob = c(p))[1], USE.NAMES = F)
  betas$`94.5%` <- sapply(1:dim(post$beta)[2], function(x) PI(post$beta[,x], prob = c(p))[2], USE.NAMES = F)
  sigmas$mean <- sapply(1:dim(post$gamma)[2], function(x) median(post$gamma[,x]), USE.NAMES = F)
  sigmas$`5.5%` <- sapply(1:dim(post$gamma)[2], function(x) PI(post$gamma[,x], prob = c(p))[1], USE.NAMES = F)
  sigmas$`94.5%` <- sapply(1:dim(post$gamma)[2], function(x) PI(post$gamma[,x], prob = c(p))[2], USE.NAMES = F)
  nu$mean <- sapply(1:dim(post$nu)[2], function(x) median(post$nu[,x]), USE.NAMES = F)
  nu$`5.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$nu[,x], prob = c(p))[1], USE.NAMES = F)
  nu$`94.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$nu[,x], prob = c(p))[2], USE.NAMES = F)
  
  # Extract true values
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  nu_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$nu[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
  
  # Build data.frames for the plots
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , lower=c(alpha_r,alphas[,3]), upper=c(alpha_r,alphas[,4]))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , lower=c(beta1_r,betas[,3]), upper=c(beta1_r,betas[,4]))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , lower=c(sigma1_r,sigmas[,3]), upper=c(sigma1_r,sigmas[,4]))
  d_nu <- data.frame(N=1:N, id= c(rep("real", length(nu_r)), rep("estimated", length(nu$mean))),value=c(nu_r,nu$mean) ,  lower=c(nu_r,nu[,3]), upper=c(nu_r,nu[,4]))
  
  # Generate plot
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.title = element_blank())
  leg <- get_legend(p1)
  p1 <-  p1 + theme(legend.position = "none")
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p3 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p4 <- ggplot(d_nu, aes(x=N, y=value, group=id, color=id))  + ggtitle("nu") +
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  
  
  figure <- grid.arrange(p1, p2, p3, p4,
                         ncol = 2, nrow = 2)
  print(figure)
  return(figure)
}

transform.skew.generror <- function(post){
  post$gamma <- post$gamma * sqrt((pi*(1+3*post$lambda*post$lambda)*gamma(3/post$nu)-(16**(1/post$nu))*post$lambda*post$lambda*gamma(0.5+1/post$nu)*gamma(0.5+1/post$nu)*gamma(1/post$nu))/(pi*gamma(1/post$nu)))
  post$beta <- post$beta - 2**(2.0/post$nu)*post$lambda*gamma(0.5+1.0/post$nu)/(sqrt(pi)*post$gamma)
  return(post)
}

plot.simulated.data.tails.skew <- function(beta=T, gp_type = 2){
  # load data
  d <- readRDS(file = paste("../../../data/processed/jsdm/skew-generror-simulated-", "data.rds", sep = ""))
  model_r <- readRDS(paste("../../../results/models/skew-generror-simulated.rds", sep=""))

  # Extract variables from the model
  alphas <- precis(model_r, pars = "alpha", depth=2)
  betas <- precis(model_r, pars = "beta", depth=3)
  sigmas <- precis(model_r, pars = "gamma", depth=3)
  nu <- precis(model_r, pars = "nu", depth=3)
  lambda <- precis(model_r, pars = "lambda", depth=3)
  
  post <- extract.samples(model_r, n = 1000, pars=c("nu", "gamma", "beta", "alpha", "lambda"))
  post <- transform.skew.generror(post)
  # PI(exp(post$nu_bar + 0.5*post$sigma_n**2)+1)
  
  p <- 0.89
  alphas$mean <- sapply(1:dim(post$alpha)[2], function(x) median(post$alpha[,x]), USE.NAMES = F)
  alphas$`5.5%` <- sapply(1:dim(post$alpha)[2], function(x) PI(post$alpha[,x], prob = c(p))[1], USE.NAMES = F)
  alphas$`94.5%` <- sapply(1:dim(post$alpha)[2], function(x) PI(post$alpha[,x], prob = c(p))[2], USE.NAMES = F)
  betas$mean <- sapply(1:dim(post$beta)[2], function(x) median(post$beta[,x]), USE.NAMES = F)
  betas$`5.5%` <- sapply(1:dim(post$beta)[2], function(x) PI(post$beta[,x], prob = c(p))[1], USE.NAMES = F)
  betas$`94.5%` <- sapply(1:dim(post$beta)[2], function(x) PI(post$beta[,x], prob = c(p))[2], USE.NAMES = F)
  sigmas$mean <- sapply(1:dim(post$gamma)[2], function(x) median(post$gamma[,x]), USE.NAMES = F)
  sigmas$`5.5%` <- sapply(1:dim(post$gamma)[2], function(x) PI(post$gamma[,x], prob = c(p))[1], USE.NAMES = F)
  sigmas$`94.5%` <- sapply(1:dim(post$gamma)[2], function(x) PI(post$gamma[,x], prob = c(p))[2], USE.NAMES = F)
  nu$mean <- sapply(1:dim(post$nu)[2], function(x) median(post$nu[,x]), USE.NAMES = F)
  nu$`5.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$nu[,x], prob = c(p))[1], USE.NAMES = F)
  nu$`94.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$nu[,x], prob = c(p))[2], USE.NAMES = F)
  lambda$mean <- sapply(1:dim(post$lambda)[2], function(x) median(post$lambda[,x]), USE.NAMES = F)
  lambda$`5.5%` <- sapply(1:dim(post$lambda)[2], function(x) PI(post$lambda[,x], prob = c(p))[1], USE.NAMES = F)
  lambda$`94.5%` <- sapply(1:dim(post$lambda)[2], function(x) PI(post$lambda[,x], prob = c(p))[2], USE.NAMES = F)
  
  # Extract true values
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  nu_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$nu[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
  lambda_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$lambda[1])
  
  # Build data.frames for the plots
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , lower=c(alpha_r,alphas[,3]), upper=c(alpha_r,alphas[,4]))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , lower=c(beta1_r,betas[,3]), upper=c(beta1_r,betas[,4]))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , lower=c(sigma1_r,sigmas[,3]), upper=c(sigma1_r,sigmas[,4]))
  d_nu <- data.frame(N=1:N, id= c(rep("real", length(nu_r)), rep("estimated", length(nu$mean))),value=c(nu_r,nu$mean) ,  lower=c(nu_r,nu[,3]), upper=c(nu_r,nu[,4]))
  d_lambda <- data.frame(N=1:N, id= c(rep("real", length(lambda_r)), rep("estimated", length(lambda$mean))),value=c(lambda_r,lambda$mean) ,  lower=c(lambda_r,lambda[,3]), upper=c(lambda_r,lambda[,4]))
  
  # Generate plot
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.title = element_blank())
  leg <- get_legend(p1)
  p1 <-  p1 + theme(legend.position = "none")
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p3 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p4 <- ggplot(d_nu, aes(x=N, y=value, group=id, color=id))  + ggtitle("nu") +
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  p5 <- ggplot(d_lambda, aes(x=N, y=value, group=id, color=id))  + ggtitle("lambda") +
    geom_pointrange(aes(ymin=lower, ymax=upper)) + theme_linedraw() + theme(legend.position = "none")
  
  figure <- grid.arrange(p1, p2, p3, p4, p5, leg,
                         ncol = 2, nrow = 3)
  print(figure)
  return(figure)
}
