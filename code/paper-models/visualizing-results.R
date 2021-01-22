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

plot.common.to.all <- function(p, fontsize=10, mar=margin(5.5,5.5,5.5,5.5)){
  p <- p +
    coord_cartesian(clip = 'off') +
    theme_bw() + 
    theme(legend.position = "none",
          text = element_text(size=fontsize),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.2),
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

plot.scatter2 <- function(mu, variance, label, color="black", xlabel="-", ylabel="-", tit=NULL, mar=margin(5.5,5.5,5.5,5.5)){
  
  df <- data.frame(sp=1:length(mu), mean=mu, variance=variance, label=label)
  
  p <- ggplot(df, aes(x=variance, y=mean, color=label)) + 
    geom_point(alpha=0.7)+
    scale_color_manual(values=color)+
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

plot.ranking.x <- function(mu, ci, color="black", xlabel="-", ylabel="-", mu_order=NULL, additional_label=c("", ""), xlims=NA, ylims=NA,mar=margin(5.5,5.5,5.5,5.5), posx=NULL, posy=NULL){
  # Sort data
  if(is.null(mu_order)){
    mu_order <- sort(mu,index.return=T)$ix
  }

  # Generate data.frame  
  df <- data.frame(sp=1:length(mu), mean=mu[mu_order], low=ci[1,mu_order], high=ci[2,mu_order])
  
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

  # Beta plots
  for(i in 1:2){
    # correlation
    corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,], post$gamma[x,i,]))
    meanPI <- mean(corPI)
    sdPI <- sd(corPI)
    corPIci <- round(sort(c(as.vector(PI(corPI)),meanPI)), 2)
    
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
    posx = (xlim[2]-xlim[1])*0.796 + xlim[1]
    
    # Beta
    plist[[1]] <- plot.ranking.x(mu_beta[[i]], ci_beta[[i]], color=colo[1], xlabel="species", ylabel="", ylim=ylim, xlim=xlim, mar=margin(5.5,5.5,5.5,5.5), additional_label = Tind_, posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])))

    xlim = c(min(sqrt(1/2*ci_gamma[[i]])), max(sqrt(1/2*ci_gamma[[i]])))
    posx = (xlim[2]-xlim[1])*0.796 + xlim[1]
    
    # Scatter plot
    plist[[2]] <- plot.scatter(mu=mu_beta[[i]], variance=sqrt(1/2*mu_gamma[[i]]), color=colo[3], xlabel="", ylabel="mean", mar=margin(5.5,5.5,5.5,5.5), ylim=ylim, xlim=xlim)
    plist[[2]] <- plist[[2]] + annotate("text", x=posx, y = posy, label=paste0(round(meanPI,2), " ± ", round(sdPI,2)), size=3)
    
    ylim = c(-1, length(sqrt(1/2*mu_gamma[[i]]))+1)
    
    # Gamma
    plist[[3]] <- plot.ranking.y(sqrt(1/2*mu_gamma[[i]]), sqrt(1/2*ci_gamma[[i]]), color=colo[2], xlabel="variance", ylabel="species", xlim=xlim, ylim=ylim, mar=margin(5.5,5.5,5.5,5.5))
    
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

plot.actual.data <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/binomial-stan-gauss-RBFs-traits2.rds")
  }
  
  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  Tend <- read.table("../../data/properties/codes/change-tendency_reindexed.csv", sep = ",")
  NEO <- read.table("../../data/properties/codes/neophytes-list_reindexed.csv", sep = ",")
  
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
  for(i in 1:2){
    idx <- sort(mu_beta[[i]],index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]

    ylim = c(min(mu_beta[[i]]), max(mu_beta[[i]]))
    xlim = c(min(sqrt(1/2*mu_gamma[[i]])), max(sqrt(1/2*mu_gamma[[i]])))
    posx = (xlim[2]-xlim[1])*0.83 + xlim[1]
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posy = c(posy, ylim[2]-(posy-ylim[1]))
    
    labels = ifelse(NEO$V3==1, "neophyte", "other")
    color <- c("red",colo[3])
    xlab=""

    p <- plot.scatter2(mu=mu_beta[[i]], variance=sqrt(1/2*mu_gamma[[i]]), label = labels, color=color, xlabel=xlab, ylabel="mean", mar=margin(5.5,5.5,5.5,5.5))
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
    labs <- c("decreasing", "decreasing", "increasing", "other")
    labels <- labs[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
    color <- colo
    xlab="variance"
    
    p <- plot.scatter2(mu=mu_beta[[i]], variance=sqrt(1/2*mu_gamma[[i]]), label = labels, color=color, xlabel=xlab, ylabel="mean", mar=margin(5.5,5.5,5.5,5.5))
    legend2 <- get_legend(p)
    if(i==1){
      p <- p+ ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }else{
      p <- p + ggtitle(label = "second axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
    }
    plist[[2*(i-1)+2]] <- p+theme(legend.position = "none", plot.title = element_blank())
    xlab="variance"
  }
  plist[[5]] <- legend1
  plist[[6]] <- legend2
  plist[[7]] <- title1
  plist[[8]] <- title2
  
  hlay <- rbind(c(7,8,NA),
                c(1,3,5),
                c(2,4,6))
  # grobs <- list()
  # widths <- list()
  # for (k in 1:length(plist)){
  #   grobs[[k]] <- ggplotGrob(plist[[k]])
  #   widths[[k]] <- grobs[[k]]$widths[2:5]
  # }    
  # maxwidth <- do.call(grid::unit.pmax, widths)
  # for (k in 1:length(grobs)){
  #   grobs[[k]]$widths[2:5] <- as.list(maxwidth)
  # }
  p <- grid.arrange(grobs=plist, ncol=3, nrow=3, heights=c(0.1,1,1), layout_matrix=hlay, widths=c(1,1,0.27)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
}

plot.actual.data.means <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../results/models/binomial-stan-gauss-RBFs-traits2.rds")
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
  


  # Generate empty list of plots
  plist = list()
    
  # Beta
  plist[[1]] <- plot.ranking.x(mu_beta[[1]], ci_beta[[1]], color=colo[1], xlabel="species", ylabel="mean", extra=1)
    
  # Scatter plot
  plist[[2]] <- plot.scatter(mu=mu_beta[[1]], variance=mu_beta[[2]], color=colo[3], xlabel="beta 2", ylabel="beta 1")
    
  # Gamma
  plist[[3]] <- plot.ranking.y(mu_beta[[2]], ci_beta[[2]], color=colo[2], xlabel="mean", ylabel="species", extra=1)
    
  hlay <- rbind(c(2,1),
                  c(3,NA))
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(0.7,1), layout_matrix=hlay, widths=c(0.7,1)#, vp=viewport(width=1, height=1, clip = TRUE),
                      # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
    )
  print(p)
    
}


plot.distribution <- function(eta, sigma, eta_prima, tit=""){
  df <- data.frame(x=c(sigma, eta, eta_prima),
                   type=c(rep("sigma", length(sigma)),
                          rep("indicator values", length(eta)),
                          rep("traits", length(eta_prima))))
  
  p <- ggplot(df, aes(x=x, color=type)) +
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    theme_bw() +
    ggtitle(tit)+
    # scale_fill_grey() +
    ylab("probability density")+
    xlab("value")+
    # guides(color = guide_legend(override.aes = list(shape=1,linetype = 1))) +
    scale_color_manual(values=c("#4daf4a", "black", "#e41a1c"))+
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
  post <- extract.samples(model, n = 1000, pars=c("etasq_b", "etasq_g", "sigma_b", "sigma_g", "etasq_tb", "etasq_tg")) 
  
  plist <- list()
  plist[[1]] <- plot.distribution(post$etasq_b[,1],
                                  post$sigma_b[,1],
                                  post$etasq_tb[,1],
                                  tit = "beta 1")
  plist[[2]] <- plot.distribution(post$etasq_b[,2],
                                  post$sigma_b[,2],
                                  post$etasq_tb[,2],
                                  tit = "beta 2")
  plist[[3]] <- plot.distribution(post$etasq_g[,1],
                                  post$sigma_g[,1],
                                  post$etasq_tg[,1],
                                  tit = "gamma 1")
  plist[[4]] <- plot.distribution(post$etasq_g[,2],
                                  post$sigma_g[,2],
                                  post$etasq_tg[,2],
                                  tit = "gamma 2")
  
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, #, vp=viewport(width=1, height=1, clip = TRUE),
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

plot.direction <- function(){
  
}

plot.heatmaps <- function(type="trait", color="black", orderid=NULL, structure=NULL){
  # load data
  # d <- readRDS(file = paste("../../data/processed/jsdm/PC1PC2data_backup", ".rds", sep = ""))

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

  #indicator values
  Tind <- read.table("../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  
  hlay <- rbind(c(1,NA),
                c(2,3))
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
  Tind_ <- as.numeric(as.character(Tind[idx,3]))
  Tind_ <- Tind_[!(is.na(Tind_))]
  Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
  Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
  
  plist[[1]] <- plot.ranking.x(mu_beta[[1]], ci_beta[[1]], color=colo[1], xlabel="species", ylabel="beta 1", extra=1, mu_order = idx, additional_label = Tind_)
  p <- plot.heatmaps(type="environment", color=colo[1], orderid=idx)
  plist[[3]] <- get_legend(p)
  plist[[2]] <- p + theme(legend.position="none")
  
  p <- grid.arrange(grobs=plist, ncol=2, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,0.1)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
  # Beta structure
  mat <- generate.association.matrix()
  plist = list()
  idx <- sort(mu_beta[[1]],index.return=T)$ix
  Tind_ <- as.numeric(as.character(Tind[idx,3]))
  Tind_ <- Tind_[!(is.na(Tind_))]
  Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
  Tind_ <- c("low elevation", "high elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
  
  plist[[1]] <- plot.ranking.x(mu_beta[[1]], ci_beta[[1]], color=colo[1], xlabel="species", ylabel="beta 1", extra=1, mu_order = idx, additional_label = Tind_)
  p <- plot.heatmaps(type="environment", color=colo[1], orderid=idx, structure = mat)
  # plist[[3]] <- get_legend(p)
  plist[[2]] <- p + theme(legend.position="none")

  p <- grid.arrange(grobs=plist, ncol=1, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1)#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
}


plot.actual.data <- function(model=NULL){
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
    beta <- as.vector(post$beta[,i,])
    beta_random <- rnorm(length(beta),0,2)
    gamma <- as.vector(post$gamma[,i,])
    gamma_random <- exp(rnorm(length(gamma),0,2))
    alpha <- as.vector(post$alpha)
    alpha_random <- rnorm(length(alpha),0,2)
    xx <- seq(from = min(beta), to = max(beta), length.out = 100)
    probability_random <- sapply(xx, function(x) sum(rbinom(length(alpha_random), 1, prob = inv_logit(alpha_random - gamma_random * (x-beta_random)^2))))
    probability <- sapply(xx, function(x) sum(rbinom(length(alpha), 1, prob = inv_logit(alpha - gamma * (x-beta)^2))))

    ggplot(data.frame(group=c(rep(1,length(xx)),rep(2,length(xx))), x=c(xx,xx), y=c(probability/max(probability),probability_random/max(probability_random))), aes(x=x, y=y, group=group)) +
          geom_line(aes(color=group))
    ggplot(data.frame(x=as.vector(post$beta[,i,])), aes(x=x)) +
      geom_density()
  }
  
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


