# Import libraries
library(rethinking)
library(rstan)
library(gridExtra)
library(grid)
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

plot.actual.data <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")

  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/binomial-stan-gauss-RBFs2.rds")
  }
  
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

  for(i in 1:2){
    
    # Sort beta and gamma
    beta_order <- sort(mu_beta[[i]],index.return=T)$ix
    gamma_order <- sort(mu_gamma[[i]],index.return=T)$ix
    
    # Generate empty list of plots
    plist = list()
      
    # Beta
    df <- data.frame(sp=1:length(mu_beta[[i]]),
                     mean=mu_beta[[i]][beta_order],
                     low=ci_beta[[i]][1,beta_order],
                     high=ci_beta[[i]][2,beta_order])
    
    plist[[1]] <- ggplot(df, aes(x=sp, y=mean)) + 
      geom_segment(aes(x=sp, xend=sp, y=low, yend=high), data=df, color=colo[1], alpha=0.5) +
      geom_point(color=colo[1], alpha=0.7)+
      ylab("mean")+
      xlab("species")+
      scale_x_continuous(expand = expansion(add = c(0, 0)),
                         limits = c(-extra, length(mu_beta[[i]])+extra),
                         breaks = seq(from=1, to=length(mu_beta[[i]]), length.out = 6))+
      scale_y_continuous(expand = expansion(add = c(0, 0)))+
      coord_cartesian(clip = 'off') +
      theme_bw() + 
      theme(legend.position = "none",
            text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.2)
            )
    
    beta_order <- sort(mu_beta[[i]],index.return=T)$ix
    gamma_order <- sort(mu_gamma[[i]],index.return=T)$ix
    
    # Scatter plot 
    df <- data.frame(sp=1:length(mu_gamma[[i]]),
                     mean=mu_beta[[i]],
                     variance=sqrt(1/2*mu_gamma[[i]]))

    plist[[2]] <- ggplot(df, aes(x=variance, y=mean)) + 
      geom_point(color=colo[3], alpha=0.7)+
      ylab("mean")+
      xlab("variance")+
      scale_x_continuous(expand = expansion(add = c(0, 0)))+
      scale_y_continuous(expand = expansion(add = c(0, 0)))+
      coord_cartesian(clip = 'off') +
      theme_bw() + 
      theme(legend.position = "none",
            text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.2)
      )

    # Gamma 
    df <- data.frame(sp=1:length(mu_gamma[[i]]),
                     mean=mu_gamma[[i]][gamma_order],
                     low=ci_gamma[[i]][1,gamma_order],
                     high=ci_gamma[[i]][2,gamma_order])
    df <- data.frame(sp=df$sp, mean=sqrt(1/(2*df$mean)), low=sqrt(1/(2*df$low)), high=sqrt(1/(2*df$high)))
    
    plist[[3]] <- ggplot(df, aes(y=sp, x=mean)) + 
      geom_segment(aes(y=sp, yend=sp, x=low, xend=high), data=df, color=colo[2], alpha=0.5) +
      geom_point(color=colo[2], alpha=0.7)+
      ylab("variance")+
      xlab("species")+
      scale_x_continuous(expand = expansion(add = c(0, 0)))+
      scale_y_continuous(expand = expansion(add = c(0, 0)),
                         limits = c(-extra, length(mu_beta[[i]])+extra),
                         breaks = seq(from=1, to=length(mu_beta[[i]]), length.out = 6))+
      coord_cartesian(clip = 'off') +
      theme_bw() + 
      theme(legend.position = "none",
            text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.2)
      )
    
    hlay <- rbind(c(1,2),
                  c(NA,3))
    p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(0.5, 1), layout_matrix=hlay, widths=c(1,0.5)#, vp=viewport(width=1, height=1, clip = TRUE),
                      # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
                      )
    print(p)
   
  }

}





