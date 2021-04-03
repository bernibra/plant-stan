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

compare.models <- function(){
  m1 <- readRDS(paste("../../../results/models/skew-min20-1d.rds", sep=""))
  m2 <- readRDS(paste("../../../results/models/baseline-min20-1d.rds", sep=""))
  m3 <- readRDS(paste("../../../results/models/generror-min20-1d.rds", sep=""))
  m4 <- readRDS(paste("../../../results/models/skew-generror-min20-1d.rds", sep=""))
  # comp <- rethinking::compare(m1, m2, m3, refresh = 1)
  
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed-30.csv", sep = " ", header=T)
  Tinvasive <- read.table("../../../data/properties/codes/neophytes-list_reindexed-30.csv", sep = " ", header=T)
  competitive <- read.table("../../../data/properties/codes/competitive_reindexed-30.csv", sep = " ", header=T)
  
  # extract samples
  post <- extract.samples(m1, pars=c("alpha", "beta", "gamma", "lambda"))
  # post <- transform.skew(post)
  post$beta <- post$betam
  
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
  
  post <- extract.samples(m3, n = 1000, pars=c("nu", "nu_bar"))
  post_ <- extract.samples(m4, n = 1000, pars=c("nu", "lambda_bar"))
  
  kurtosis.generror <- gamma(5/post$nu)*gamma(1/post$nu)/(gamma(3/post$nu))**2-3
  kurtosis.generror <- apply(kurtosis.generror, 1, median)
  
  nu_hat <- abs(post$nu_bar)+1
  nu_hat <- gamma(5/nu_hat)*gamma(1/nu_hat)/(gamma(3/nu_hat))**2-3
  
  
  plist <- list()
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
  
  post <- extract.samples(m4, n = 1000, pars=c("lambda_bar")) 
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

plot.actual.data <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  n <- 1
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS("../../../results/models/baseline-min20-1d.rds")
  }
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  # post$beta <- post$beta*(-1)
  # post$gamma <- 1/(2*post$gamma)
  
  # alphas
  mu_alpha <- apply( post$alpha , 2 , mean )
  ci_alpha <- apply( post$alpha , 2 , PI )
  
  # Beta plots
  for(i in 1:n){
    
    if(n==1){
      mu_beta <- apply( post$beta , 2 , mean )
      ci_beta <- apply( post$beta , 2 , PI )
      mu_gamma <- apply( post$gamma , 2 , mean )
      ci_gamma <- apply( post$gamma , 2 , PI )
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,], post$gamma[x,]))
    }else{
      mu_beta <- apply( post$beta[,i,] , 2 , mean )
      ci_beta <- apply( post$beta[,i,] , 2 , PI )
      mu_gamma <- apply( post$gamma[,i,] , 2 , mean )
      ci_gamma <- apply( post$gamma[,i,] , 2 , PI )
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,i,], post$gamma[x,i,]))
    }
    
    # correlation
    meanPI <- mean(corPI)
    sdPI <- sd(corPI)
    corPIci <- round(sort(c(as.vector(PI(corPI)),meanPI)), 2)
    print(corPIci)
    
    # Generate empty list of plots
    plist = list()
    idx <- sort(mu_beta,index.return=T)$ix
    Tind_ <- as.numeric(as.character(Tind[idx,3]))
    Tind_ <- Tind_[!(is.na(Tind_))]
    Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
    Tind_ <- c("high elevation", "low elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
    ylim = c(min(ci_beta), max(ci_beta))
    xlim = c(-1, length(mu_beta)+1)
    posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
    posx = (xlim[2]-xlim[1])*0.2 + xlim[1]
    
    # Beta
    plist[[1]] <- plot.ranking.x(mu_beta, ci_beta, color=colo[1], xlabel="species", ylabel="", ylims=ylim, xlims=xlim, mar=margin(5.5,5.5,5.5,5.5), additional_label = Tind_, posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])))+
      coord_trans(x="identity")
    
    xlim = c(min(ci_gamma), max(ci_gamma))
    posx = exp((log(xlim[2])-log(xlim[1]))*0.796 + log(xlim[1]))
    
    # Scatter plot
    plist[[2]] <- plot.scatter(mu=mu_beta, variance=mu_gamma, color=colo[3], xlabel="", ylabel=expression(beta), mar=margin(5.5,5.5,5.5,5.5), ylims=ylim, xlims=xlim)
    breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
    plist[[2]] <- plist[[2]] + coord_trans(x="log")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    plist[[2]] <- plist[[2]] + annotate("text", x=posx, y = ylim[2]-(posy-ylim[1]), label=paste0(round(meanPI,2), " ± ", round(sdPI,2)), size=3)
    
    ylim = c(-1, length(mu_gamma)+1)
    
    # Gamma
    plist[[3]] <- plot.ranking.y(mu_gamma, ci_gamma, color=colo[2], xlabel=expression(gamma), ylabel="species", xlims=xlim, ylims=ylim, mar=margin(5.5,5.5,5.5,5.5))
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
    model <- readRDS("../../../results/models/baseline-model.rds")
  }
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed.csv", sep = ",")
  
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
    model <- readRDS("../../../results/models/baseline-min20-1d.rds")
  }
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed-20.csv", sep = " ")
  Tend <- read.table("../../../data/properties/codes/change-tendency_reindexed-20.csv", sep = " ")
  NEO <- read.table("../../../data/properties/codes/neophytes-list_reindexed-20.csv", sep = " ")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  
  # post$beta <- post$beta*(-1) 
  # post$gamma <- 1/(2*post$gamma)
  
  # alphas
  mu_alpha <- apply( post$alpha , 2 , mean )
  ci_alpha <- apply( post$alpha , 2 , PI )
  
  # betas
  mu_beta <- apply(post$beta,2,mean)
  ci_beta <- apply(post$beta,2,PI)
  
  # gammas
  mu_gamma <- apply(post$gamma,2,mean)
  ci_gamma <- apply(post$gamma,2,PI)

  plist = list()
  ylab <- expression(beta)
  
  idx <- sort(mu_beta,index.return=T)$ix
  Tind_ <- as.numeric(as.character(Tind[idx,3]))
  Tind_ <- Tind_[!(is.na(Tind_))]
  Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
  Tind_ <- c("high elevation", "low elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
  ylim = c(min(mu_beta), max(mu_beta))
  xlim = c(min(mu_gamma), max(mu_gamma))
  posx = exp((log(xlim[2])-log(xlim[1]))*0.83 + log(xlim[1]))
  posy = (ylim[2]-ylim[1])*0.1 + ylim[1]
  posy = c(posy, ylim[2]-(posy-ylim[1]))
    
  labels = ifelse(NEO$V3==1, "neophyte", "native")
  alphas = ifelse(NEO$V3==1, 0.7, 0.6)
  color <- c("#e7298a",colo[4])
  xlab=""
  breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
    
  correlations <- c()
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,NEO$V3==1], post$gamma[x,NEO$V3==1]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(c(meanPI, sdPI), 2)
  correlations <- rbind(correlations, c("neophytes", corPIci))
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,NEO$V3!=1], post$gamma[x,NEO$V3!=1]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(c(meanPI, sdPI), 2)
  correlations <- rbind(correlations, c("native", corPIci))
    
  p <- plot.scatter2(mu=mu_beta, variance=mu_gamma, label = labels, color=c(color[2],color[1]), alpha = alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
  # Scatter plot
  p <- p + coord_trans(x="log", clip = "off")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    
  p <- p+annotate("text", x = posx, y = posy, label = Tind_, colour = "#525252", size=3)+
        ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
  title1 <- get_title(p)

  legend1 <- get_legend(p)
  plist[[1]] <- p+theme(legend.position = "none", plot.title = element_blank())
  labs <- c("decreasing", "decreasing low", "increasing", "stable")
  labels <- labs[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
  color <- colo
  alphas <- c(0.7,0.7,0.7,0.6)
  alphas <- alphas[Tend$V3+Tend$V4*2+Tend$V5*3+Tend$V6*4]
  xlab=expression(gamma)
    
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V3==1], post$gamma[x,Tend$V3==1]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(c(meanPI, sdPI), 2)
  correlations <- rbind(correlations, c("decreasing", corPIci))
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V4==1], post$gamma[x,Tend$V4==1]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(c(meanPI, sdPI), 2)
  correlations <- rbind(correlations, c("decreasing low", corPIci))
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V5==1], post$gamma[x,Tend$V5==1]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(c(meanPI, sdPI), 2)
  correlations <- rbind(correlations, c("increasing", corPIci))
  corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V6==1], post$gamma[x,Tend$V6==1]))
  meanPI <- mean(corPI)
  sdPI <- sd(corPI)
  corPIci <- round(c(meanPI, sdPI), 2)
  correlations <- rbind(correlations, c("stable", corPIci))

  p <- plot.scatter2(mu=mu_beta, variance=mu_gamma, label = labels, color=color, alpha=alphas, xlabel=xlab, ylabel=ylab, mar=margin(5.5,5.5,5.5,5.5))
  legend2 <- get_legend(p)
  p <- p + coord_trans(x="log", clip = "off")+
  scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
  p <- p+ ggtitle(label = "first axis") + theme(plot.title = element_text(hjust = 0.1, size=10))
  plist[[2]] <- p+theme(legend.position = "none", plot.title = element_blank())

  plist[[3]] <- legend1
  plist[[4]] <- legend2


  tab <- ggtexttable(correlations, rows = NULL,cols = c("", "mean", "sd"),
                     theme = ttheme("light",  base_size = 9))
  tab <- tab %>%
    tab_add_hline(at.row = 4, row.side = "top", linetype = 2)
  
  plist[[5]] <- tab

  hlay <- rbind(c(5,1,3),
                c(5,2,4))
  
  p <- grid.arrange(grobs=plist, ncol=3, nrow=2, heights=c(1,1), layout_matrix=hlay, widths=c(1,1,0.4)#, vp=viewport(width=1, height=1, clip = TRUE),
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
    model <- readRDS("../../../results/models/baseline-min20-1d.rds")
  }
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed-20.csv", sep = " ")
  Tend <- read.table("../../../data/properties/codes/change-tendency_reindexed-20.csv", sep = " ")
  NEO <- read.table("../../../data/properties/codes/neophytes-list_reindexed-20.csv", sep = " ")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=c("alpha", "gamma", "beta")) 
  
  # post$gamma <- 1/(2*post$gamma)
  post$alpha <- exp(-post$alpha)
  
  plist = list()
  ylab <- expression(beta)
  for(i in 1:2){
    if(i==1){
      # betas
      mu_beta <- apply(post$beta,2,mean)
      ci_beta <- apply(post$beta,2,PI)
      # gammas
      mu_gamma <- apply(post$gamma,2,mean)
      ci_gamma <- apply(post$gamma,2,PI)
      transformation <- "log"
    }else{
      # alphas
      mu_beta <- apply( post$alpha , 2 , mean )
      ci_beta <- apply( post$alpha , 2 , PI )
      
      mu_gamma <- apply(post$gamma,2,mean)
      ci_gamma <- apply(post$gamma,2,PI)
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
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,NEO$V3==1], post$gamma[x,NEO$V3==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("neophytes", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,NEO$V3!=1], post$gamma[x,NEO$V3!=1]))
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
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V3==1], post$gamma[x,Tend$V3==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("decreasing", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V4==1], post$gamma[x,Tend$V4==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("decreasing low", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V5==1], post$gamma[x,Tend$V5==1]))
      meanPI <- mean(corPI)
      sdPI <- sd(corPI)
      corPIci <- round(c(meanPI, sdPI), 2)
      correlations <- rbind(correlations, c("increasing", corPIci))
      corPI <- sapply(1:dim(post$beta)[1], function(x) cor(post$beta[x,Tend$V6==1], post$gamma[x,Tend$V6==1]))
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

