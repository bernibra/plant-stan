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
  m1 <- readRDS(paste("../../results/models/skew-model-traits-1dskew-simulated2.rds", sep=""))
  m2 <- readRDS(paste("../../results/models/baseline-model-1d1d-simulated2.rds", sep=""))
  rethinking::compare(m1, m2)
  
  d <- readRDS(file = paste("../../data/processed/jsdm/skew-simulated2S1S2", "data.rds", sep = ""))

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
  d <- readRDS(file = paste("../../data/processed/jsdm/skew-simulated-workingexample-S1S2", "data.rds", sep = ""))
  model_r <- readRDS(paste("../../results/models/skew-workingexample-simulated.rds", sep=""))

  # Extract variables from the model
  alphas <- precis(model_r, pars = "alpha", depth=2)
  betas <- precis(model_r, pars = "beta", depth=3)
  sigmas <- precis(model_r, pars = "gamma", depth=3)
  lambdas <- precis(model_r, pars = "lambda", depth=3)
  precis(model_r, pars = "lambda_bar", depth=3)
  precis(model_r, pars = "sigma_l", depth=3)
  
    
  # Extract true values
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  lambda_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$lambda[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])

  # Build data.frames for the plots
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , sd=c(rep(0,length(alpha_r)),alphas$sd))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , sd=c(rep(0,length(beta1_r)),betas[1:N,]$sd))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , sd=c(rep(0,length(sigma1_r)),sigmas[1:N,]$sd))
  d_lambda <- data.frame(N=1:N, id= c(rep("real", length(lambda_r)), rep("estimated", length(lambdas$mean))),value=c(lambda_r,lambdas$mean) , sd=c(rep(0,length(lambda_r)),lambdas$sd))
  
  
  # Generate plot
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.title = element_blank())
  leg <- get_legend(p1)
  p1 <-  p1 + theme(legend.position = "none")
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p3 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p4 <- ggplot(d_lambda, aes(x=N, y=value, group=id, color=id))  + ggtitle("lambda") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  
  
  figure <- grid.arrange(p1, p2, p3, p4,
                         ncol = 2, nrow = 2)
  print(figure)
  return(figure)
}

plot.simulated.data.tails <- function(beta=T, gp_type = 2){
  # load data
  d <- readRDS(file = paste("../../data/processed/jsdm/1d-generror-simulated-S1S2-", "data.rds", sep = ""))
  model_r <- readRDS(paste("../../results/models/-1d-generror-simulated.rds", sep=""))
  
  # Extract variables from the model
  alphas <- precis(model_r, pars = "alpha", depth=2)
  betas <- precis(model_r, pars = "beta", depth=3)
  sigmas <- precis(model_r, pars = "gamma", depth=3)
  nu <- precis(model_r, pars = "nu", depth=3)
  
  post <- extract.samples(model_r, n = 1000, pars=c("nu", "beta", "gamma"))
  post$gamma <- post$gamma
  
  nu$mean <- sapply(1:dim(post$nu)[2], function(x) mean(post$nu[,x]), USE.NAMES = F)
  nu$`5.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$nu[,x], prob = c(0.890))[1], USE.NAMES = F)
  nu$`94.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$nu[,x], prob = c(0.890))[2], USE.NAMES = F)
  sigmas$mean <- sapply(1:dim(post$nu)[2], function(x) mean(post$sigmas[,x]), USE.NAMES = F)
  sigmas$`5.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$sigmas[,x], prob = c(0.890))[1], USE.NAMES = F)
  sigmas$`94.5%` <- sapply(1:dim(post$nu)[2], function(x) PI(post$sigmas[,x], prob = c(0.890))[2], USE.NAMES = F)
  
  # Extract true values
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  nu_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$nu[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
  
  # Build data.frames for the plots
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , sd=c(rep(0,length(alpha_r)),alphas$sd))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , sd=c(rep(0,length(beta1_r)),betas[1:N,]$sd))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , sd=c(rep(0,length(sigma1_r)),sigmas[1:N,]$sd))
  d_nu <- data.frame(N=1:N, id= c(rep("real", length(nu_r)), rep("estimated", length(nu$mean))),value=c(nu_r,nu$mean) , sd=c(rep(0,length(nu_r)),nu$sd))
  
  
  # Generate plot
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + ggtitle("alpha") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.title = element_blank())
  leg <- get_legend(p1)
  p1 <-  p1 + theme(legend.position = "none")
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + ggtitle("beta 1") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p3 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id))  + ggtitle("gamma 1") + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + theme(legend.position = "none")
  p4 <- ggplot(d_nu, aes(x=N, y=value, group=id, color=id))  + ggtitle("nu") +
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw() + ylim(-1,4) + theme(legend.position = "none")
  
  
  figure <- grid.arrange(p1, p2, p3, p4,
                         ncol = 2, nrow = 2)
  print(figure)
  return(figure)
}

compare.models <- function(){
  m1 <- readRDS(paste("../../results/models/min30-skew-model-traits-1d2.rds", sep=""))
  m2 <- readRDS(paste("../../results/models/min30-baseline-model-1d2.rds", sep=""))
  m3 <- readRDS(paste("../../results/models/min30-1d-generror.rds", sep=""))
  comp <- rethinking::compare(m1, m2, m3, refresh = 1)
  
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed-30.csv", sep = " ", header=T)
  Tinvasive <- read.table("../../data/properties/codes/neophytes-list_reindexed-30.csv", sep = " ", header=T)
  competitive <- read.table("../../data/properties/codes/competitive_reindexed-30.csv", sep = " ", header=T)
  
  # extract samples
  post <- extract.samples(m1, pars=c("alpha", "beta", "gamma", "lambda"))
  post <- transform.skew(post)
  
  # alphas
  mu_lambda <- apply( post$lambda , 2 , mean )
  ci_lambda <- apply( post$lambda , 2 , PI )
  
  # betas
  mu_beta <- apply(post$beta,2,mean)
  ci_beta <- apply(post$beta,2,PI)
  
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,], post$lambda[x,]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(sort(c(as.vector(PI(corPI)),meanPI)), 2)
  print(c(meanPI, corPIci))
  plist <- list()
  
  # Generate empty list of plots
  idx <- sort(mu_beta,index.return=T)$ix
  Tind_ <- as.numeric(as.character(Tind[idx,3]))
  Tind_ <- Tind_[!(is.na(Tind_))]
  Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
  Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
  
  ylim = c(min(ci_beta), max(ci_beta))
  xlim = c(-1, length(mu_beta)+1)
  posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
  posx = (xlim[2]-xlim[1])*0.8 + xlim[1]
    
  mu_order <- sort(mu_beta,index.return=T)$ix
  
  # Beta
  plist[[1]] <- plot.ranking.x(mu_beta, ci_beta, color=colo[1], mu_order = mu_order, xlabel="species", ylabel=expression(beta), ylims=ylim, xlims=xlim, mar=margin(5.5,5.5,5.5,5.5), additional_label = Tind_, posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])))+
      coord_trans(x="identity")
  
  ylim = c(min(ci_lambda), max(ci_lambda))
  xlim = c(-1, length(mu_lambda)+1)
  posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
  posx = (xlim[2]-xlim[1])*0.8 + xlim[1]
  
  invasive <- Tinvasive$V3
  plist[[2]] <- plot.ranking.x(mu_lambda, ci_lambda, color=colo[2], mu_order = mu_order, xlabel="species", ylabel=expression(lambda), ylims=ylim, xlims=xlim, mar=margin(5.5,5.5,5.5,5.5), posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])), colorpoints=invasive)+
    coord_trans(x="identity")+ geom_hline(yintercept=0, linetype="dashed", color = "gray")
  
  grobs <- list()
  for (k in 1:length(plist)){
    grobs[[k]] <- ggplotGrob(plist[[k]])
  }   
  p <- grid.arrange(grobs=grobs, ncol=1, nrow=2#, vp=viewport(width=1, height=1, clip = TRUE),
                      # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)

  post <- extract.samples(m3, n = 1000, pars=c("nu_bar")) 
  
  plist <- list()
  nu_hat <- exp(post$nu_bar)+1
  nu_hat <- gamma(5/nu_hat)*gamma(1/nu_hat)/(gamma(3/nu_hat))**2-3
  plist[[1]] <- ggplot(data.frame(x=nu_hat), aes(x=x))+
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    theme_bw() +
    geom_vline(xintercept = 0, colour = "red", linetype="dashed") +
    ggtitle("(a)  average kurtosis of distributions")+
    # scale_fill_grey() +
    ylab("density")+
    xlab("kurtosis")+
    scale_x_continuous(limits = c(-0.5, 0.5))+
    theme(text = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          # legend.position="none",
          legend.position=c(0.70,.80),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(size=10))
  
  post <- extract.samples(m1, n = 1000, pars=c("lambda_bar")) 
  lambda_bar <- -post$lambda_bar
  lambda_bar <- lambda_bar/sqrt(1+lambda_bar**2)
  lambda_bar <- 0.5*(4-pi)*((lambda_bar*sqrt(2/pi))**3/(1-lambda_bar*lambda_bar*(2/pi))**(3/2))
  plist[[2]] <- ggplot(data.frame(x=lambda_bar), aes(x=x))+
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    theme_bw() +
    geom_vline(xintercept = 0, colour = "red", linetype="dashed") +
    ggtitle("(b)  average skewness of distributions")+
    # scale_fill_grey() +
    ylab(" ")+
    xlab("skewness")+
    scale_x_continuous(limits = c(-0.1, 0.75))+
    theme(text = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          # legend.position="none",
          legend.position=c(0.70,.80),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(size=10))
  p <- grid.arrange(grobs=plist, ncol=2, nrow=1#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )

}

plot.common.to.all <- function(p, fontsize=10, mar=margin(5.5,5.5,5.5,5.5)){
  p <- p +
    coord_cartesian(clip = 'off') +
    theme_bw() + 
    theme(legend.position = "none",
          text = element_text(size=fontsize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.4),
          plot.margin = mar
    )

  return(p)
}

plot.scatter <- function(mu, variance, color="black", xlabel="-", ylabel="-", tit=NULL, mar=margin(5.5,5.5,5.5,5.5), xlims=NA, ylims=NA){

  df <- data.frame(sp=1:length(mu), mean=mu, variance=variance)
  
  p <- ggplot(df, aes(x=variance, y=mean)) + 
    geom_point(color=color, alpha=0.7)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_x_continuous(expand = expansion(add = c(0, 0)), limits = xlims)+
    scale_y_continuous(expand = expansion(add = c(0, 0)), limits = ylims)
  
  p <- plot.common.to.all(p, mar=mar)
  if(!(is.null(tit))){
    p <- p + ggtitle(tit)
  }
  
  return(p)
}

plot.scatter2 <- function(mu, variance, label, color="black", alpha=NULL, xlabel="-", ylabel="-", tit=NULL, mar=margin(5.5,5.5,5.5,5.5)){
  
  if(is.null(alpha)){
    alpha <- rep(0.7, length(mu))
  }
  df <- data.frame(sp=1:length(mu), mean=mu, variance=variance, label=label, alpha=alpha)
  
  p <- ggplot(df, aes(x=variance, y=mean, color=label, alpha=alpha)) + 
    geom_point()+
    scale_color_manual(values=color)+
    guides(alpha = FALSE)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_x_continuous(expand = expansion(add = c(0, 0)))+
    scale_y_continuous(expand = expansion(add = c(0, 0)))
  
  p <- p +
    coord_cartesian(clip = 'off') +
    theme_bw() + 
    theme(legend.title = element_blank(),
          text = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.2),
          plot.margin = mar
    )
  if(!(is.null(tit))){
    p <- p + ggtitle(tit)
  }
  
  return(p)
}

plot.ranking.x <- function(mu, ci, color="black", xlabel="-", ylabel="-", mu_order=NULL, additional_label=c("", ""), xlims=NA, ylims=NA,mar=margin(5.5,5.5,5.5,5.5), posx=NULL, posy=NULL, colorpoints=NULL){
  # Sort data
  if(is.null(mu_order)){
    mu_order <- sort(mu,index.return=T)$ix
  }

  if(is.null(colorpoints)){
    colorpoints <- rep(0,length(mu_order))
  }
  
  # Generate data.frame  
  df <- data.frame(sp=1:length(mu), mean=mu[mu_order], low=ci[1,mu_order], high=ci[2,mu_order], color=colorpoints[mu_order])
  
  p <- ggplot(df, aes(x=sp, y=mean)) + 
    geom_segment(aes(x=sp, xend=sp, y=low, yend=high), data=df, color=color, alpha=0.5) +
    geom_point(color=color, alpha=0.7)+
    annotate("text", x = posx, y = posy, label = additional_label, colour = "#525252", size=3) +
    ylab(ylabel)+
    xlab(xlabel)+
    scale_x_continuous(expand = expansion(add = c(0, 0)),
                       limits = xlims,
                       breaks = seq(from=1, to=length(mu), length.out = 6))+
    scale_y_continuous(expand = expansion(add = c(0, 0)), limits = ylims)
  
  p <- plot.common.to.all(p, mar=mar)
  
  return(p)
}

plot.ranking.y <- function(mu, ci, color="black", xlabel="-", ylabel="-", mu_order=NULL, xlims=NA, ylims=NA, mar=margin(5.5,5.5,5.5,5.5)){
  # Sort data
  if(is.null(mu_order)){
    mu_order <- sort(mu,index.return=T)$ix 
  }
  
  # Generate data.frame  
  df <- data.frame(sp=1:length(mu),
                   mean=mu[mu_order],
                   low=ci[1,mu_order],
                   high=ci[2,mu_order])
  
  p <- ggplot(df, aes(y=sp, x=mean)) + 
    geom_segment(aes(y=sp, yend=sp, x=low, xend=high), data=df, color=color, alpha=0.5) +
    geom_point(color=color, alpha=0.7)+
    xlab(xlabel)+
    ylab(ylabel)+
    scale_x_continuous(expand = expansion(add = c(0, 0)), limits=xlims)+
    scale_y_continuous(expand = expansion(add = c(0, 0)),
                       limits = ylims,
                       breaks = seq(from=1, to=length(mu), length.out = 6))

  p <- plot.common.to.all(p, mar=mar)
  
  return(p)
}

link.model <- function(data, post){
  N <- ncol(post$alpha)
  M <- nrow(data)
  p <- list()
  a <- matrix(1,1,M)
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 1:N){
    setTxtProgressBar(pb, i)
    beta1 <- as.matrix(post$beta[,1,i]) %*% a
    beta2 <- as.matrix(post$beta[,2,i]) %*% a
    gamma1 <- as.matrix(post$gamma[,1,i]) %*% a
    gamma2 <- as.matrix(post$gamma[,2,i]) %*% a
    alpha <- t(as.matrix(post$alpha[,i]) %*% a)
    
    p[[i]] <- inv_logit(post$alpha[,i] - post$gamma[,1,i] * (t(data$x1[col(beta1)] - beta1))**2 - post$gamma[,2,i] * (t(data$x2[col(beta2)] - beta2))**2)
  }
  close(pb)
  return(p)
}

plot.actual.data <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")

  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/baseline-model.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  post$beta <- post$beta*(-1)
  
  # alphas
  mu_alpha <- apply( post$alpha , 2 , mean )
  ci_alpha <- apply( post$alpha , 2 , PI )
  
  # betas
  mu_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,mean))
  ci_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,PI))
  
  # gammas
  mu_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,mean))
  ci_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,PI))

  # Beta plots
  for(i in 1:1){
    # correlation
    corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,], post$gamma[x,i,]))
    meanPI <- mean(corPI)
    sdPI <- sd(corPI)
    corPIci <- round(sort(c(as.vector(PI(corPI)),meanPI)), 2)
    print(corPIci)
    
    # Generate empty list of plots
    plist = list()
    idx <- sort(mu_beta[[i]],index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
    ylim = c(min(ci_beta[[i]]), max(ci_beta[[i]]))
    xlim = c(-1, length(mu_beta[[i]])+1)
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posx = (xlim[2]-xlim[1])*0.2 + xlim[1]
    
    # Beta
    plist[[1]] <- plot.ranking.x(mu_beta[[i]], ci_beta[[i]], color=colo[1], xlabel="species", ylabel="", ylims=ylim, xlims=xlim, mar=margin(5.5,5.5,5.5,5.5), additional_label = Tind_, posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])))+
      coord_trans(x="identity")

    xlim = c(min(ci_gamma[[i]]), max(ci_gamma[[i]]))
    posx = exp((log(xlim[2])-log(xlim[1]))*0.796 + log(xlim[1]))
    
    # Scatter plot
    plist[[2]] <- plot.scatter(mu=mu_beta[[i]], variance=mu_gamma[[i]], color=colo[3], xlabel="", ylabel=expression(beta), mar=margin(5.5,5.5,5.5,5.5), ylims=ylim, xlims=xlim)
    breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
    plist[[2]] <- plist[[2]] + coord_trans(x="log")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    plist[[2]] <- plist[[2]] + annotate("text", x=posx, y = ylim[2]-(posy-ylim[1]), label=paste0(round(meanPI,2), " ± ", round(sdPI,2)), size=3)
    
    ylim = c(-1, length(mu_gamma[[i]])+1)
    
    # Gamma
    plist[[3]] <- plot.ranking.y(mu_gamma[[i]], ci_gamma[[i]], color=colo[2], xlabel=expression(gamma), ylabel="species", xlims=xlim, ylims=ylim, mar=margin(5.5,5.5,5.5,5.5))
    plist[[3]] <- plist[[3]] + coord_trans(x="log") +
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    
    hlay <- rbind(c(2,1),
                  c(3,NA))
    grobs <- list()
    widths <- list()
    for (k in 1:length(plist)){
      grobs[[k]] <- ggplotGrob(plist[[k]])
      widths[[k]] <- grobs[[k]]$widths[2:5]
    }    
    maxwidth <- do.call(grid::unit.pmax, widths)
    for (k in 1:length(grobs)){
      grobs[[k]]$widths[2:5] <- as.list(maxwidth)
    }
    p <- grid.arrange(grobs=grobs, ncol=2, nrow=2, heights=c(0.8,1), layout_matrix=hlay, widths=c(0.8,1)#, vp=viewport(width=1, height=1, clip = TRUE),
                      # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
                      )
    print(p)
   
  }

}

plot.actual.data.alpha <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/baseline-model.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  
  # alphas
  mu_alpha <- apply( post$alpha , 2 , mean )
  ci_alpha <- apply( post$alpha , 2 , PI )
  
  # betas
  mu_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,mean))
  ci_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,PI))
  
  # gammas
  mu_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,mean))
  ci_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,PI))
  
  # Beta plots
  for(i in 1:2){
    if(i==1){
      target <- post$beta[,i,]      
      mu_target <- mu_beta[[1]]
      ci_target <- ci_beta[[1]]
    }else{
      target <- post$gamma[,i,]
      mu_target <- mu_gamma[[1]]
      ci_target <- ci_gamma[[1]]
    }
    
    # correlation
    corPI <- sapply(1:dim(post$beta)[1], function(x) cor(target[x,], post$alpha[x,]))
    meanPI <- mean(corPI)
    sdPI <- sd(corPI)
    corPIci <- round(sort(c(as.vector(PI(corPI)),meanPI)), 2)
    
    # Generate empty list of plots
    plist = list()
    idx <- sort(mu_target,index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
    ylim = c(min(ci_target), max(ci_target))
    xlim = c(-1, length(mu_target)+1)
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posx = (xlim[2]-xlim[1])*0.8 + xlim[1]
    if(i==2){
      posy = exp((log(ylim[2])-log(ylim[1]))*0.1 + log(ylim[1]))
    }
    
    # Beta
    plist[[1]] <- plot.ranking.x(mu_target,ci_target, color=colo[1], xlabel="species", ylabel="", ylims=ylim, xlims=xlim, mar=margin(5.5,5.5,5.5,5.5), additional_label = Tind_, posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])))
    xlim = c(min(ci_alpha), max(ci_alpha))
    posx = (xlim[2]-xlim[1])*0.796 + xlim[1]
    
    if(i==2){
      plist[[2]] <- plot.scatter(mu=mu_target, variance=mu_alpha, color=colo[3], xlabel="", ylabel=expression(gamma), mar=margin(5.5,5.5,5.5,5.5), ylims=ylim, xlims=xlim)
      breaks <- round(seq(from=ylim[1], to=ylim[2], length.out = 4))+1
      plist[[1]] <- plist[[1]] + coord_trans(y="log")+
        scale_y_continuous(breaks = breaks, labels = breaks, limits = ylim)
      plist[[2]] <- plist[[2]] + coord_trans(y="log")+
      scale_y_continuous(breaks = breaks, labels = breaks, limits = ylim)
    }else{
      # Scatter plot
      plist[[2]] <- plot.scatter(mu=mu_target, variance=mu_alpha, color=colo[3], xlabel="", ylabel=expression(beta), mar=margin(5.5,5.5,5.5,5.5), ylims=ylim, xlims=xlim)
    }
    
    plist[[2]] <- plist[[2]] + annotate("text", x=posx, y = posy, label=paste0(round(meanPI,2), " ± ", round(sdPI,2)), size=3)
    
    ylim = c(-1, length(mu_alpha)+1)
    
    # Gamma
    plist[[3]] <- plot.ranking.y(mu_alpha, ci_alpha, color=colo[2], xlabel=expression(alpha), ylabel="species", xlims=xlim, ylims=ylim, mar=margin(5.5,5.5,5.5,5.5))

    hlay <- rbind(c(2,1),
                  c(3,NA))
    grobs <- list()
    widths <- list()
    for (k in 1:length(plist)){
      grobs[[k]] <- ggplotGrob(plist[[k]])
      widths[[k]] <- grobs[[k]]$widths[2:5]
    }    
    maxwidth <- do.call(grid::unit.pmax, widths)
    for (k in 1:length(grobs)){
      grobs[[k]]$widths[2:5] <- as.list(maxwidth)
    }
    p <- grid.arrange(grobs=grobs, ncol=2, nrow=2, heights=c(0.8,1), layout_matrix=hlay, widths=c(0.8,1)#, vp=viewport(width=1, height=1, clip = TRUE),
                      # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
    )
    print(p)
    
  }
}

plot.actual.data.means <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#e6ab02", "#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/baseline-model.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  Tend <- read.table("../../data/properties/codes/change-tendency_reindexed.csv", sep = ",")
  NEO <- read.table("../../data/properties/codes/neophytes-list_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  
  post$beta <- post$beta*(-1) 
  
  # alphas
  mu_alpha <- apply( post$alpha , 2 , mean )
  ci_alpha <- apply( post$alpha , 2 , PI )
  
  # betas
  mu_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,mean))
  ci_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,PI))
  
  # gammas
  mu_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,mean))
  ci_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,PI))
  
  mu_beta_dec <- as.vector(post$beta[,1,(Tend$V3+Tend$V4)==1])
  mu_beta_inc <- as.vector(post$beta[,1,Tend$V5==1])
  mu_gamma_dec <- as.vector(post$gamma[,1,(Tend$V3+Tend$V4)==1])
  mu_gamma_inc <- as.vector(post$gamma[,1,Tend$V5==1])
  
  plist <- list()
  df <- data.frame(group=c("other", "increasing")[c(rep(1, length(mu_beta_dec)), rep(2, length(mu_beta_inc)))],
                   x=c(mu_beta_dec, mu_beta_inc))
  plist[[1]] <- ggplot(df, aes(x=x, color=group, fill=group)) +
    geom_boxplot(alpha=0.3)
  df <- data.frame(group=c("other", "increasing")[c(rep(1, length(mu_gamma_dec)), rep(2, length(mu_gamma_inc)))],
                   x=c(mu_gamma_dec, mu_gamma_inc))
  plist[[2]] <- ggplot(df, aes(x=x, color=group, fill=group)) +
    geom_boxplot(alpha=0.3)
  p <- grid.arrange(grobs=plist, ncol=2, nrow=1)
  
  plist = list()
  ylab <- expression(beta)
  for(i in 1:2){
    idx <- sort(mu_beta[[i]],index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]

    ylim = c(min(mu_beta[[i]]), max(mu_beta[[i]]))
    xlim = c(min(mu_gamma[[i]]), max(mu_gamma[[i]]))
    posx = exp((log(xlim[2])-log(xlim[1]))*0.83 + log(xlim[1]))
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posy = c(posy, ylim[2]-(posy-ylim[1]))
    
    labels = ifelse(NEO$V3==1, "neophyte", "native")
    alphas = ifelse(NEO$V3==1, 0.7, 0.6)
    color <- c("#e7298a",colo[4])
    xlab=""
    breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
    
    if(i==1){
      correlations <- c()
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,NEO$V3==1], post$gamma[x,1,NEO$V3==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("neophytes", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,NEO$V3!=1], post$gamma[x,1,NEO$V3!=1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("native", corPIci))
    }
    
    p <- plot.scatter2(mu=mu_beta[[i]], variance=mu_gamma[[i]], label = labels, color=c(color[2],color[1]), alpha = alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
    # Scatter plot
    p <- p + coord_trans(x="log", clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    
    if(i==1){
      p <- p+annotate("text", x = posx, y = posy, label = Tind_, colour = "#525252", size=3)+
        ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
      title1 <- get_title(p)
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
      title2 <- get_title(p)
    }
    legend1 <- get_legend(p)
    plist[[2*(i-1)+1]] <- p+theme(legend.position = "none", plot.title = element_blank())
    # 
    labs <- c("decreasing", "decreasing low", "increasing", "stable")
    labels <- labs[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
    color <- colo
    alphas <- c(0.7,0.7,0.7,0.6)
    alphas <- alphas[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
    xlab=expression(gamma)
    
    if(i==1){
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V3==1], post$gamma[x,1,Tend$V3==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("decreasing", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V4==1], post$gamma[x,1,Tend$V4==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("decreasing low", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V5==1], post$gamma[x,1,Tend$V5==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("increasing", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V6==1], post$gamma[x,1,Tend$V6==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("stable", corPIci))
    }
    
    p <- plot.scatter2(mu=mu_beta[[i]], variance=mu_gamma[[i]], label = labels, color=color, alpha=alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
    legend2 <- get_legend(p)
    p <- p + coord_trans(x="log", clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    if(i==1){
      p <- p+ ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }
    plist[[2*(i-1)+2]] <- p+theme(legend.position = "none", plot.title = element_blank())
    ylab <- ""
  }
  plist[[5]] <- legend1
  plist[[6]] <- legend2
  plist[[7]] <- title1
  plist[[8]] <- title2
  
  hlay <- rbind(c(7,8,NA),
                c(1,3,5),
                c(2,4,6))

  p <- grid.arrange(grobs=plist, ncol=3, nrow=3, heights=c(0.1,1,1), layout_matrix=hlay, widths=c(1,1,0.4)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
  plist2 <- list()
  tab <- ggtexttable(correlations, rows = NULL,cols = c("", "mean", "sd"),
                            theme = ttheme("light",  base_size = 9))
  tab <- tab %>%
    tab_add_hline(at.row = 4, row.side = "top", linetype = 2)
  
  plist2[[1]] <- tab
  plist2[[2]] <- plist[[1]]
  plist2[[3]] <- plist[[2]] 
  plist2[[4]] <- plist[[5]]
  plist2[[5]] <- plist[[6]]
  
  hlay <- rbind(c(1,2,4),
                c(1,3,5))
  
  p <- grid.arrange(grobs=plist2, ncol=3, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,1,0.4)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
}

plot.actual.data.means2 <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#e6ab02", "#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/baseline-model.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  Tend <- read.table("../../data/properties/codes/change-tendency_reindexed.csv", sep = ",")
  NEO <- read.table("../../data/properties/codes/neophytes-list_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  
  post$beta <- post$beta*(-1) 
  post$alpha <- exp(-post$alpha)
 
  plist = list()
  ylab <- expression(beta)
  for(i in 1:2){
    if(i==1){
      # betas
      mu_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,mean))[[1]]
      ci_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,PI))[[1]]
      
      # gammas
      mu_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,mean))[[1]]
      ci_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,PI))[[1]]
      transformation <- "log"
    }else{
      # alphas
      mu_beta <- apply( post$alpha , 2 , mean )
      ci_beta <- apply( post$alpha , 2 , PI )
      
      mu_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,mean))[[1]]
      ci_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,PI))[[1]]
      transformation <- "log"
      
    }
    
    idx <- sort(mu_beta,index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
    ylim = c(min(mu_beta), max(mu_beta))
    xlim = c(min(mu_gamma), max(mu_gamma))
    posx = exp((log(xlim[2])-log(xlim[1]))*0.83 + log(xlim[1]))
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posy = c(posy, ylim[2]-(posy-ylim[1]))
    
    labels = ifelse(NEO$V3==1, "neophyte", "native")
    alphas = ifelse(NEO$V3==1, 0.7, 0.6)
    color <- c("#e7298a",colo[4])
    xlab=""
    
    if(i==1){
      breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
      
      correlations <- c()
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,NEO$V3==1], post$gamma[x,1,NEO$V3==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("neophytes", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,NEO$V3!=1], post$gamma[x,1,NEO$V3!=1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("native", corPIci))
    }else{
      breaks <- seq(from=0,to=1,length.out = 6)
    }
    
    p <- plot.scatter2(mu=mu_beta, variance=mu_gamma, label = labels, color=c(color[2],color[1]), alpha = alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
    # Scatter plot
    p <- p + coord_trans(x=transformation, clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    
    if(i==1){
      p <- p+annotate("text", x = posx, y = posy, label = Tind_, colour = "#525252", size=3)+
        ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
      title1 <- get_title(p)
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
      title2 <- get_title(p)
    }
    legend1 <- get_legend(p)
    plist[[2*(i-1)+1]] <- p+theme(legend.position = "none", plot.title = element_blank())
    # 
    labs <- c("decreasing", "decreasing low", "increasing", "stable")
    labels <- labs[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
    color <- colo
    alphas <- c(0.7,0.7,0.7,0.6)
    alphas <- alphas[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
  
    if(i==1){
      xlab=expression(gamma)
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V3==1], post$gamma[x,1,Tend$V3==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("decreasing", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V4==1], post$gamma[x,1,Tend$V4==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("decreasing low", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V5==1], post$gamma[x,1,Tend$V5==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("increasing", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,Tend$V6==1], post$gamma[x,1,Tend$V6==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("stable", corPIci))
    }else{
      xlab=expression(exp(-alpha))
    }
    
    p <- plot.scatter2(mu=mu_beta, variance=mu_gamma, label = labels, color=color, alpha=alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
    legend2 <- get_legend(p)
    p <- p + coord_trans(x=transformation, clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    if(i==1){
      p <- p+ ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }
    plist[[2*(i-1)+2]] <- p+theme(legend.position = "none", plot.title = element_blank())
    ylab <- ""
  }
  plist[[5]] <- legend1
  plist[[6]] <- legend2
  plist[[7]] <- title1
  plist[[8]] <- title2
  
  hlay <- rbind(c(7,8,NA),
                c(1,3,5),
                c(2,4,6))
  
  p <- grid.arrange(grobs=plist, ncol=3, nrow=3, heights=c(0.1,1,1), layout_matrix=hlay, widths=c(1,1,0.4)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
  plist2 <- list()
  tab <- ggtexttable(correlations, rows = NULL,cols = c("", "mean", "sd"),
                     theme = ttheme("light",  base_size = 9))
  tab <- tab %>%
    tab_add_hline(at.row = 4, row.side = "top", linetype = 2)
  
  plist2[[1]] <- tab
  plist2[[2]] <- plist[[1]]
  plist2[[3]] <- plist[[2]] 
  plist2[[4]] <- plist[[5]]
  plist2[[5]] <- plist[[6]]
  plist2[[6]] <- plist[[3]]
  plist2[[7]] <- plist[[4]]
  
  hlay <- rbind(c(1,2,6,4),
                c(1,3,7,5))
  
  p <- grid.arrange(grobs=plist2, ncol=4, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,1,1,0.4)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
}

plot.actual.data.means.pairwise <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#e6ab02", "#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/min30-skew-model-traits-1d2.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  Tend <- read.table("../../data/properties/codes/change-tendency_reindexed.csv", sep = ",")
  NEO <- read.table("../../data/properties/codes/neophytes-list_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta", "lambda")) 
  post <- transform.skew(post)
  
  # lambda
  mu_lambda <- apply( post$lambda , 2 , mean )
  ci_lambda <- apply( post$lambda , 2 , PI )
  
  # betas
  mu_beta <- apply(post$beta,2,mean)
  ci_beta <- apply(post$beta,2,PI)
  
  # gamma
  mu_gamma <- apply( post$sigma , 2 , mean )
  ci_gamma <- apply( post$sigma , 2 , PI )
  
  # alpha
  mu_alpha <- apply(post$alpha,2,mean)
  ci_alpha <- apply(post$alpha,2,PI)
  
  parameters_post <- list(beta=post$beta, gamma=post$sigma, alpha=post$alpha, lambda=post$lambda)
  parameters_mu <- list(beta=mu_beta, gamma=mu_gamma, alpha=mu_alpha, lambda=mu_lambda)
  parameters_ci <- list(beta=ci_beta, gamma=ci_gamma, alpha=ci_alpha, lambda=ci_lambda)

  transformation <- c("identity", "log", "log", "identity")
  label <- c(expression(beta), expression(gamma), expression(alpha), expression(lambda))
    
  plist = list()
  k <- 1
  for(i in 1:4){
    for(j in 1:4){
      # correlation
      corPI <- sapply(1:dim(parameters_post[[i]])[1], function(x) cor(parameters_post[[i]][x,], parameters_post[[j]][x,]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(sort(c(as.vector(PI(corPI)),meanPI)), 2)
      ylim = c(min(parameters_mu[[i]]), max(parameters_mu[[i]]))
      xlim = c(min(parameters_mu[[j]]), max(parameters_mu[[j]]))

      p <- plot.scatter(mu=parameters_mu[[i]], variance=parameters_mu[[j]], color=colo[i], xlabel=label[i], ylabel=label[j], xlims = xlim, ylims = ylim)
      # p <- p + coord_trans(x=transformation[i], y=transformation[j], clip = "off")
      plist[[k]] <- p
      k <- k+1
    }
  }
  p <- grid.arrange(grobs=plist, ncol=4, nrow=4 #, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
  for(i in 1:2){
    idx <- sort(mu_beta[[i]],index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
    ylim = c(min(mu_beta[[i]]), max(mu_beta[[i]]))
    xlim = c(min(mu_gamma[[i]]), max(mu_gamma[[i]]))
    posx = exp((log(xlim[2])-log(xlim[1]))*0.83 + log(xlim[1]))
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posy = c(posy, ylim[2]-(posy-ylim[1]))
    
    labels = ifelse(NEO$V3==1, "neophyte", "other")
    alphas = ifelse(NEO$V3==1, 0.7, 0.6)
    color <- c("#e7298a",colo[4])
    xlab=""
    breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
    
    p <- plot.scatter2(mu=mu_beta[[i]], variance=mu_gamma[[i]], label = labels, color=color, alpha = alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
    # Scatter plot
    p <- p + coord_trans(x="log", clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    
    if(i==1){
      p <- p+annotate("text", x = posx, y = posy, label = Tind_, colour = "#525252", size=3)+
        ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
      title1 <- get_title(p)
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
      title2 <- get_title(p)
    }
    legend1 <- get_legend(p)
    plist[[2*(i-1)+1]] <- p+theme(legend.position = "none", plot.title = element_blank())
    # 
    labs <- c("decreasing", "decreasing low", "increasing", "other")
    labels <- labs[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
    color <- colo
    alphas <- c(0.7,0.7,0.7,0.6)
    alphas <- alphas[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
    xlab=expression(gamma)
    
    p <- plot.scatter2(mu=mu_beta[[i]], variance=mu_gamma[[i]], label = labels, color=color, alpha=alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
    legend2 <- get_legend(p)
    p <- p + coord_trans(x="log", clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    if(i==1){
      p <- p+ ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }
    plist[[2*(i-1)+2]] <- p+theme(legend.position = "none", plot.title = element_blank())
    ylab <- ""
  }
  plist[[5]] <- legend1
  plist[[6]] <- legend2
  plist[[7]] <- title1
  plist[[8]] <- title2
  
  hlay <- rbind(c(7,8,NA),
                c(1,3,5),
                c(2,4,6))
  
  p <- grid.arrange(grobs=plist, ncol=3, nrow=3, heights=c(0.1,1,1), layout_matrix=hlay, widths=c(1,1,0.4)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
}

plot.actual.data.distribution <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/binomial-stan-gauss-RBFs-traits2.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  
  # alphas
  mu_alpha <- apply( post$alpha , 2 , mean )
  ci_alpha <- apply( post$alpha , 2 , PI )
  
  # betas
  mu_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,mean))
  ci_beta <- lapply(1:2, function(x) apply(post$beta[,x,],2,PI))
  
  # gammas
  mu_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,mean))
  ci_gamma <- lapply(1:2, function(x) apply(post$gamma[,x,],2,PI))
  
  plist = list()
  
  # Beta plots
  for(i in 1:2){
    xx <- seq(from = min(post$beta[,i,]), to = max(post$beta[,i,]), length.out = 100)
    beta <- lapply(1:1000, function(j) sapply(xx, function(x) sum((post$beta[j,i,]<=x)*1)/length(post$beta[j,i,])))
    beta <- matrix(unlist(beta), ncol = length(xx), byrow = TRUE)
    mu_beta <- apply( beta , 2 , mean )
    ci_beta <- apply( beta , 2 , PI )
    
    beta_normal <- matrix(rnorm(ncol(post$beta[,i,])*1000, mean=mean(post$beta[,i,]), sd=sd(post$beta[,i,])), 1000, ncol = ncol(post$beta[,i,]))
    beta_normal <- lapply(1:nrow(beta_normal), function(j) sapply(xx, function(x) sum((beta_normal[j,]<=x)*1)/ncol(beta_normal)))
    beta_normal <- matrix(unlist(beta_normal), ncol = length(xx), byrow = TRUE)
    mu_beta_normal <- apply( beta_normal , 2 , mean )
    ci_beta_normal <- apply( beta_normal , 2 , PI )
    
    beta_unif <- matrix(runif(ncol(post$beta[,i,])*1000, min =min(post$beta[,i,]), max=max(post$beta[,i,])), 1000, ncol = ncol(post$beta[,i,]))
    beta_unif <- lapply(1:nrow(beta_unif), function(j) sapply(xx, function(x) sum((beta_unif[j,]<=x)*1)/ncol(beta_unif)))
    beta_unif <- matrix(unlist(beta_unif), ncol = length(xx), byrow = TRUE)
    mu_beta_unif <- apply( beta_unif , 2 , mean )
    ci_beta_unif <- apply( beta_unif , 2 , PI )
    
    df <- data.frame(x=c(xx,xx,xx),
                     mean=c(mu_beta, mu_beta_normal, mu_beta_unif),
                     group=c(rep("real", length(mu_beta)), rep("normal", length(mu_beta_normal)), rep("unif", length(mu_beta_normal))),
                     lower=c(ci_beta[1,], ci_beta_normal[1,], ci_beta_unif[1,]),
                     upper=c(ci_beta[2,], ci_beta_normal[2,], ci_beta_unif[2,]))
    
    p <- ggplot(data=df, aes(x=x, y=mean, colour=group, fill=group)) +
         geom_line() +
        geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.3, colour = NA)+
        theme_bw()+
        theme(legend.title = element_blank(),
              legend.position = c(0.8, 0.2),
              legend.background = element_blank())
    print(p)
    # 
    # beta <- as.vector(post$beta[,i,])
    # beta_normal <- rnorm(length(beta),mean(beta),sd(beta))
    # beta_unif <- runif(length(beta),min=min(beta),max = max(beta))
    # gamma <- as.vector(post$gamma[,i,])
    # gamma_random <- exp(rnorm(length(gamma),mean(log(gamma)),sd(log(gamma))))
    # alpha <- as.vector(post$alpha)
    # alpha_random <- exp(rnorm(length(alpha),mean(log(alpha)),sd(log(alpha))))
    # xx <- seq(from = min(beta), to = max(beta), length.out = 100)
    # beta_c <- sapply(xx, function(x) sum((beta<=x)*1))
    # beta_c <- beta_c/max(beta_c)
    # beta_normal_c <- sapply(xx, function(x) sum((beta_normal<=x)*1))
    # beta_normal_c <- beta_normal_c/max(beta_normal_c)
    # beta_unif_c <- sapply(xx, function(x) sum((beta_unif<=x)*1))
    # beta_unif_c <- beta_unif_c/max(beta_unif_c)
    # probability_random <- sapply(xx, function(x) sum(rbinom(length(alpha_random), 1, prob = exp(- alpha_random - gamma_random * (x-beta_normal)^2))))
    # probability_unif <- sapply(xx, function(x) sum(rbinom(length(alpha_random), 1, prob = exp(- alpha_random - gamma_random * (x-beta_unif)^2))))
    # probability <- sapply(xx, function(x) sum(rbinom(length(alpha), 1, prob = exp(-alpha - gamma * (x-beta)^2))))
    
    # df=data.frame(group=c("real", "normal", "uniform")[c(rep(1,length(xx)),rep(2,length(xx)),rep(3,length(xx)))],
    #               x=c(xx,xx,xx),
    #               y=c(probability/max(probability),probability_random/max(probability_random), probability_unif/max(probability_unif)))
    # ggplot(df, aes(x=x, y=y, color=group)) +
    #   stat_ecdf(geom = "step")
    
    df=data.frame(group=c("real", "normal", "uniform")[c(rep(1,length(beta_c)),rep(2,length(beta_c)),rep(3,length(beta_c)))],
                  y=c(beta_c,beta_normal_c,beta_unif_c),
                  x=c(xx,xx,xx))
    ggplot(df, aes(x=x, y=y, color=group)) +
      geom_line()
  }
  
}

plot.distribution <- function(df, colors,legend, tit="", ylim=NULL){
  
  p <- ggplot(df, aes(x=x, color=type)) +
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    theme_bw() +
    ggtitle(tit)+
    # scale_fill_grey() +
    ylab("probability density")+
    xlab("value")+
    # guides(color = guide_legend(override.aes = list(shape=1,linetype = 1))) +
    scale_color_manual(values=colors)+
    theme(text = element_text(size=10),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          # legend.position="none",
          legend.position=c(0.70,.80),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(size=10))
  
  if(!legend){
    p <- p + theme(legend.position = "none")
  }
  if(!is.null(ylim)){
    p <- p + coord_cartesian(ylim=ylim)
  }
  return(p)
}

plot.distributions.gp <- function(model=NULL){

  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/binomial-stan-gauss-RBFs2.rds")
  }
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("etasq_a","etasq_b", "etasq_g", "rhosq_a", "rhosq_b", "rhosq_g", "sigma_b", "sigma_g", "sigma_a")) 
  
  sigmainfo <- c("#e6ab02", "#525252")
  
  plist <- list()
  df <- data.frame(x=c(post$sigma_b[,1], post$etasq_b[,1]),
                   type=c(rep("sigma", length(post$sigma_b[,1])),
                          rep("prior information", length(post$etasq_b[,1]))))
  plist[[1]] <- plot.distribution(df,colors=sigmainfo, legend=T,
                                  tit = paste("a)   ",expression(beta),sep=""), ylim=c(0,2.5))
  
  df <- data.frame(x=c(post$rhosq_b[,1,1], post$rhosq_b[,1,2]),
                   type=c(rep("indicator values", length(post$sigma_b[,1])),
                          rep("traits", length(post$etasq_b[,1]))))
  plist[[2]] <- plot.distribution(df,colors=c("#1b9e77", "#d95f02"), legend=T,
                                  tit = "")
  
  df <- data.frame(x=c(post$sigma_g[,1], post$etasq_g[,1]),
                   type=c(rep("sigma", length(post$sigma_g[,1])),
                          rep("information", length(post$etasq_g[,1]))))
  plist[[3]] <- plot.distribution(df,colors=sigmainfo, legend=F,
                                  tit = paste("b)   ",expression(gamma),sep=""))
  
  df <- data.frame(x=c(post$rhosq_g[,1,1], post$rhosq_g[,1,2]),
                   type=c(rep("indicator values", length(post$sigma_b[,1])),
                          rep("traits", length(post$etasq_b[,1]))))
  plist[[4]] <- plot.distribution(df,colors=c("#1b9e77", "#d95f02"), legend=F,
                                  tit = "")
  
  df <- data.frame(x=c(post$sigma_a, post$etasq_a),
                   type=c(rep("sigma", length(post$sigma_a)),
                          rep("information", length(post$etasq_a))))
  plist[[5]] <- plot.distribution(df,colors=sigmainfo, legend=F, 
                                  tit = paste("c)   ",expression(alpha),sep=""))
  
  p <- grid.arrange(grobs=plist, ncol=2, nrow=3, #, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
  plist = list()
  
  # Scatter plot
  plist[[1]] <- plot.scatter(mu=post$etasq_b[,1], variance=post$etasq_tb[,1], color="#fbb4ae", xlabel="traits", ylabel="indicator values", tit="beta 1")
  plist[[2]] <- plot.scatter(mu=post$etasq_b[,2], variance=post$etasq_tb[,2], color="#b3cde3", xlabel="traits", ylabel="indicator values", tit="beta 2")
  plist[[3]] <- plot.scatter(mu=post$etasq_g[,1], variance=post$etasq_tg[,1], color="#ccebc5", xlabel="traits", ylabel="indicator values", tit="gamma 1")
  plist[[4]] <- plot.scatter(mu=post$etasq_g[,2], variance=post$etasq_tg[,2], color="#decbe4", xlabel="traits", ylabel="indicator values", tit="gamma 2")
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, #, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
}

generate.association.matrix <- function(){
  comm <- readRDS("../../data/properties/communities/communities.rds")
  dsp <- readRDS("../../data/properties/communities/dictionary.rds")
  dictionary <- read.table("../../data/properties/codes/dictionary.csv", sep=",", header = T)
  
  # load data
  # d <- readRDS(file = paste("../../data/processed/jsdm/PC1PC2data_backup", ".rds", sep = ""))
  
  # Find dimensions
  L <- length(unique(d$dataset$id))
  N <- sum(d$dataset$id==1)
  
  # Generate base matrix
  results <- matrix(0,L,L)
  
  # indices
  real <- matrix(d$dataset$real.id, N, L)
  
  # Extract indices  
  indices = real[1,]
  dis_names <- dictionary[indices,]$new.names
  
  for (i in comm){
    if(sum(i %in% dis_names)>1){
      results[dis_names %in% i, dis_names %in% i] <- 1
    }
  }
  results <- results*(1-diag(L))
  return(results)
}


shade_curve <- function(MyDF, zstart, zend, fill = "red", alpha = .5, mean.1=0, sd.1=1){
  geom_area(data = subset(MyDF, x >= mean.1 + zstart*sd.1
                          & x < mean.1 + zend*sd.1),
            aes(y=y), fill = fill, color = NA, alpha = alpha)
}

prior_p <- function(N=2, sd_=1){
  set.seed(2)
  xmin=-2
  xmax=2
  a <- rnorm(1e6, 0,sd_)
  c <- rnorm(1e6, 0,sd_)
  b <- rnorm(N, 0,sd_)
  x <- seq(from=xmin, to=xmax, length.out = 1000)
  
  lims_a <- PI(a, prob=c(0.3333333333,1))
  lims_c <- PI(c, prob=c(0.33333333333,1))
  plist <- list()
  l <- 1
  
  for(i in 1:length(lims_a)){
    for(j in 2:(length(lims_c)+1)){
      if(i==1){
        if (j==(length(lims_c) + 1)){
          plist[[l]] <- ggplot() + theme_void()
        }else{
          xx <- seq(from=-5, to=5, length.out = 1000)
          MyDF <- data.frame(x = xx, y = dnorm(xx, mean = 0, sd = sd_))
          MyDF$y <- MyDF$y/max(MyDF$y)
          
          plist[[l]] <- ggplot(MyDF, aes(x = x, y = y)) + geom_line() +
            shade_curve(MyDF = MyDF, zstart = lims_c[(j-1)], zend = lims_c[(j)], fill = my_col, alpha = .3, sd.1 = sd_) +
            theme_bw()+
            annotate("text", x=0,y=0.4, label = expression(log(lambda)), size = 3)+
            theme(text = element_text(size=10),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title = element_blank(),
                  axis.text = element_blank(),
                  plot.margin = margin(0,0,0,30),
                  legend.position="none",
                  panel.border = element_blank())
        }
      }else if(j==(length(lims_c) + 1)){
        xx <- seq(from=-5, to=5, length.out = 1000)
        MyDF <- data.frame(x = xx, y = dnorm(xx, mean = 0, sd = sd_))
        MyDF$y <- MyDF$y/max(MyDF$y)
        
        plist[[l]] <- ggplot(MyDF, aes(x = x, y = y)) + geom_line() +
          shade_curve(MyDF = MyDF, zstart = lims_a[(i-1)], zend = lims_a[(i)], fill = my_col, alpha = .3, sd.1 = sd_) +
          theme_bw()+
          annotate("text", x=0,y=0.5, label = expression(log(alpha)), size = 3)+
          coord_flip()+
          theme(text = element_text(size=10),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                plot.margin = margin(0,0,10,0),
                legend.position="none",
                panel.border = element_blank())
      
      }else{
        a_ <- exp(sample(a[a>lims_a[(i-1)] & a<lims_a[(i)]], size = N))
        c_ <- exp(sample(c[c>lims_c[(j-1)] & c<lims_c[(j)]], size = N))
        
        dat <- data.frame(x=NULL, y=NULL, variable=NULL)
        for (k in 1:N){
          prob <- exp(-a_[k] - c_[k] * (x-b[k])**2)
          dat <- rbind(dat, data.frame(x=x, y=prob, variable=paste0("line",k)))
        }
        plist[[l]] <- ggplot(data=dat, aes(x=x, y=y)) + 
          geom_line(aes(colour=variable))+
          scale_color_brewer() +
          annotate("text", x = 1, y = 1.1, label="text") +
          scale_y_continuous(expand = expansion(add = c(0, 0)), limits = c(0,1))+
          scale_x_continuous(expand = expansion(add = c(0, 0)), limits = c(xmin, xmax))+
          theme_bw()+
          theme(text = element_text(size=10),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(colour = "black"),
                plot.margin = margin(5,5,0,0),
                legend.position="none",
                panel.border = element_rect(colour = "black", fill=NA, size=0.3))
        
        
        if(j!=2){
          plist[[l]] <- plist[[l]] + labs(y = "")
        }else{
          plist[[l]] <- plist[[l]] + labs(y = "probability")
        }
        if(i!=length(lims_a)){
          plist[[l]] <- plist[[l]] + labs(x = "")
        }else{
          plist[[l]] <- plist[[l]] + labs(x = expression(beta))
        }
      }
      l <- l+1
    }
  }

  p <- grid.arrange(grobs=plist, ncol=4, nrow=4, widths=c(1,1,1,0.3), heights=c(0.3,1,1,1))
  print(p)
}


