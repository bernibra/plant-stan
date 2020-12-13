# Import libraries
library(rethinking)
library(rstan)
library(gridExtra)
library(grid)
library(ggpubr)
library(hrbrthemes)
library(plotly)

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

plot.common.to.all <- function(p, fontsize=10){
  p <- p +
    coord_cartesian(clip = 'off') +
    theme_bw() + 
    theme(legend.position = "none",
          text = element_text(size=fontsize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.2)
    )
  return(p)
}

plot.scatter <- function(mu, variance, color="black", xlabel="-", ylabel="-"){

  df <- data.frame(sp=1:length(mu), mean=mu, variance=variance)
  
  p <- ggplot(df, aes(x=variance, y=mean)) + 
    geom_point(color=color, alpha=0.7)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_x_continuous(expand = expansion(add = c(0, 0)))+
    scale_y_continuous(expand = expansion(add = c(0, 0)))
  
  p <- plot.common.to.all(p)
  
  return(p)
}

plot.ranking.x <- function(mu, ci, color="black", xlabel="-", ylabel="-", extra=1, mu_order=NULL){
  # Sort data
  if(is.null(mu_order)){
    mu_order <- sort(mu,index.return=T)$ix
  }

  # Generate data.frame  
  df <- data.frame(sp=1:length(mu), mean=mu[mu_order], low=ci[1,mu_order], high=ci[2,mu_order])
  
  p <- ggplot(df, aes(x=sp, y=mean)) + 
    geom_segment(aes(x=sp, xend=sp, y=low, yend=high), data=df, color=color, alpha=0.5) +
    geom_point(color=color, alpha=0.7)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_x_continuous(expand = expansion(add = c(0, 0)),
                       limits = c(-extra, length(mu)+extra),
                       breaks = seq(from=1, to=length(mu), length.out = 6))+
    scale_y_continuous(expand = expansion(add = c(0, 0)))
  
  p <- plot.common.to.all(p)
  
  return(p)
}

plot.ranking.y <- function(mu, ci, color="black", xlabel="-", ylabel="-", extra=1, mu_order=NULL){
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
    scale_x_continuous(expand = expansion(add = c(0, 0)))+
    scale_y_continuous(expand = expansion(add = c(0, 0)),
                       limits = c(-extra, length(mu)+extra),
                       breaks = seq(from=1, to=length(mu), length.out = 6))

  p <- plot.common.to.all(p)
  
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
  
  ## Generate probabilities
  # data <- expand.grid(x1=seq(from=-3, to=3, length.out = 50), x2=seq(from=-3, to=3, length.out = 50))
  # p <- link.model(data, post)
  # p_mu <- lapply(1:length(p),function(x) apply(p[[x]], 1, mean))
  
  # Alpha plots
  plist = list()
  idx <- sort(mu_beta[[1]],index.return=T)$ix
  plist[[1]] <- plot.ranking.x(mu_beta[[1]], ci_beta[[1]], color=colo[1], xlabel="species", ylabel="beta 1", extra=1, mu_order = idx)
  plist[[2]] <- plot.ranking.x(sqrt(1/(2*mu_gamma[[1]])), sqrt(1/(2*ci_gamma[[1]])), color=colo[2], xlabel="species", ylabel="beta 2", extra=1, mu_order = idx)
  plist[[3]] <- plot.ranking.x(mu_alpha, ci_alpha, color=colo[3], xlabel="species", ylabel="alpha", extra=1, mu_order = idx)
  
  p <- grid.arrange(grobs=plist, ncol=1, nrow=3)
  print(p)
  
  p <- plot.scatter(mu=mu_beta[[1]], variance=mu_beta[[2]], color=colo[3], xlabel="variance", ylabel="mean")
  print(p)
  
  # Beta plots
  for(i in 1:2){
    # Generate empty list of plots
    plist = list()
      
    # Beta
    plist[[1]] <- plot.ranking.x(mu_beta[[i]], ci_beta[[i]], color=colo[1], xlabel="species", ylabel="mean", extra=1)

    # Scatter plot
    plist[[2]] <- plot.scatter(mu=mu_beta[[i]], variance=sqrt(1/2*mu_gamma[[i]]), color=colo[3], xlabel="variance", ylabel="mean")
    
    # Gamma
    plist[[3]] <- plot.ranking.y(sqrt(1/2*mu_gamma[[i]]), sqrt(1/2*ci_gamma[[i]]), color=colo[2], xlabel="variance", ylabel="species", extra=1)
    
    hlay <- rbind(c(2,1),
                  c(3,NA))
    p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(0.7,1), layout_matrix=hlay, widths=c(0.7,1)#, vp=viewport(width=1, height=1, clip = TRUE),
                      # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
                      )
    print(p)
   
  }

}

plot.distribution <- function(eta, sigma, tit=""){
  df <- data.frame(x=c(eta, sigma), type=c(rep("eta", length(eta)), rep("sigma", length(sigma))))
  
  p <- ggplot(df, aes(x=x, color=type)) +
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    theme_bw() +
    ggtitle(tit)+
    # scale_fill_grey() +
    ylab("probability density")+
    xlab("value")+
    # guides(color = guide_legend(override.aes = list(shape=1,linetype = 1))) +
    scale_color_manual(values=c("black", "#ce5c00"))+
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
          panel.border = element_rect(colour = "black", fill=NA, size=0.5))
}

plot.distributions.gp <- function(model=NULL){

  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/binomial-stan-gauss-RBFs2.rds")
  }
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("etasq_b", "etasq_g", "sigma_b", "sigma_g")) 
  
  plist <- list()
  plist[[1]] <- plot.distribution(post$etasq_b[,1], post$sigma_b[,1], tit = "beta 1")
  plist[[2]] <- plot.distribution(post$etasq_b[,2], post$sigma_b[,2], tit = "beta 2")
  plist[[3]] <- plot.distribution(post$etasq_g[,1], post$sigma_g[,1], tit = "gamma 1")
  plist[[4]] <- plot.distribution(post$etasq_g[,2], post$sigma_g[,2], tit = "gamma 2")
  
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, #, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
}

plot.direction <- function(){
  
}

plot.heatmaps <- function(type="trait", color="black", orderid=NULL, structure=NULL){
  # load data
  d <- readRDS(file = paste("../../data/processed/jsdm/PC1PC2", "data.rds", sep = ""))

  # Find dimensions
  L <- length(unique(d$dataset$id))
  N <- sum(d$dataset$id==1)
  
  # indices
  mmsbm <- matrix(d$dataset$mmsbm.id, N, L)
  
  # Extract indices  
  indices = data.frame(mmsbm = mmsbm[1,])

  # Distance matrices  
  dmat <- as.matrix(read.table(paste("../../data/properties/distance-matrices/", type, ".csv", sep=""), sep=","))
  
  # reshape correlation matrices
  dmat <- dmat[,indices$mmsbm]
  dmat <- dmat[indices$mmsbm,]
  
  # Normalize
  dmat <- dmat/max(dmat)
  
  if(!(is.null(structure))){
    dmat <- structure
  }
  
  # reorder
  if(!is.null(orderid)){
    dmat <- dmat[, orderid]
    dmat <- dmat[orderid,]
  }
  
  # Build dataset
  ij <- which(dmat>-1, arr.ind = T)
  data <- data.frame(x=ij[,1], y=ij[,2], z=log(1+dmat[ij]))
  p <- ggplot(data, aes(x, y, fill= z)) + 
    geom_tile()+
    ylab("species")+
    xlab("species")+
    scale_y_continuous(expand = expansion(add = c(0, 0)), trans = "reverse")+
    scale_x_continuous(expand = expansion(add = c(0, 0)))+
    scale_fill_gradient(low="white", high=color)+
    theme_bw()+
    theme(text = element_text(size=10),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          # legend.position="none",
          # legend.position=c(0.70,.80),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3))
  return(p)
}

generate.association.matrix <- function(){
  comm <- readRDS("../../data/properties/communities/communities.rds")
  dsp <- readRDS("../../data/properties/communities/dictionary.rds")
  dictionary <- read.table("../../data/properties/codes/dictionary.csv", sep=",", header = T)
  
  # load data
  d <- readRDS(file = paste("../../data/processed/jsdm/PC1PC2", "data.rds", sep = ""))
  
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
  
  for (i in comm.sp){
    if(sum(i %in% dis_names)>1){
      results[dis_names %in% i, dis_names %in% i] <- 1
    }
  }
  results <- results*(1-diag(L))
  return(results)
}

plot.some.heatmaps <- function(model=NULL){
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

  # Alpha  plots
  plist = list()
  idx <- sort(mu_alpha,index.return=T)$ix
  plist[[1]] <- plot.ranking.x(mu_alpha, ci_alpha, color=colo[3], xlabel="species", ylabel="alpha", extra=1, mu_order = idx)
  p <- plot.heatmaps(type="trait", color=colo[3], orderid=idx)
  plist[[3]] <- get_legend(p)
  plist[[2]] <- p + theme(legend.position="none")

  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,0.1))  
  
  # Beta 1  plots
  plist = list()
  idx <- sort(mu_beta[[1]],index.return=T)$ix
  plist[[1]] <- plot.ranking.x(mu_beta[[1]], ci_beta[[1]], color=colo[1], xlabel="species", ylabel="beta 1", extra=1, mu_order = idx)
  p <- plot.heatmaps(type="environment", color=colo[1], orderid=idx)
  plist[[3]] <- get_legend(p)
  plist[[2]] <- p + theme(legend.position="none")
  
  hlay <- rbind(c(1,NA),
                c(2,3))
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,0.1)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
  # Beta structure
  mat <- generate.association.matrix()
  plist = list()
  idx <- sort(mu_beta[[1]],index.return=T)$ix
  plist[[1]] <- plot.ranking.x(mu_beta[[1]], ci_beta[[1]], color=colo[1], xlabel="species", ylabel="beta 1", extra=1, mu_order = idx)
  p <- plot.heatmaps(type="environment", color=colo[1], orderid=idx, structure = mat)
  plist[[3]] <- get_legend(p)
  plist[[2]] <- p + theme(legend.position="none")
  

  hlay <- rbind(c(1,NA),
                c(2,3))
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,0.1)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
}



