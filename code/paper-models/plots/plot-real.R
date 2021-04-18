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

kurtosis.skewness <- function(lambda){
  v <- ((pi-8*lambda*lambda+3*pi*lambda*lambda)/(2*pi))**(-0.5)
  kurtosis <- (v**4)*(-192*(lambda**4)+16*pi*(lambda**2)*(-5+lambda**2)+3*pi*pi*(1+10*lambda*lambda+5*(lambda**4)))/(4*(pi**2))-3
  return(kurtosis)
}
kurtosis.generror <- function(nu){
  kurtosis <- gamma(5/nu)*gamma(1/nu)/((gamma(3/nu))**2)-3
  return(kurtosis)
}
kurtosis.skew.generror <- function(nu, lambda){
  v <- ((pi*(1+3*lambda*lambda)*gamma(3/nu)-((16)**(1/nu))*lambda*lambda*(gamma(0.5+1/nu)**2)*gamma(1/nu))/(pi*gamma(1/nu)))**(-1/2)
  kurtosis <- ((v**4)/(pi*pi*gamma(1/nu)))*(-3*(256**(1/nu))*(lambda**4)*(gamma(0.5+1/nu)**4)*gamma(1/nu) +
                                                  3*(2**((4+nu)/nu))*pi*(lambda**2)*(1+3*(lambda**2))*(gamma(0.5+1/nu)**2)*gamma(3/nu)-
                                                  (2**(4+2/nu))*(pi**(3/2))*(lambda**2)*(1+lambda**2)*gamma(0.5+1/nu)*gamma(4/nu)+
                                                  (pi**2)*(1+10*(lambda**2)+5*(lambda**4))*gamma(5/nu))-3
  return(kurtosis)
}

skewness.skew <- function(lambda){
  v <- ((pi-8*lambda*lambda+3*pi*lambda*lambda)/(2*pi))**(-0.5)
  skewness <- (v**3)*lambda*(16*lambda*lambda-pi*(-1+5*lambda*lambda))/(pi**(3/2))
  return(skewness)
}
skewness.skew.generror <- function(nu, lambda){
  v <- ((pi*(1+3*lambda*lambda)*gamma(3/nu)-((16)**(1/nu))*lambda*lambda*(gamma(0.5+1/nu)**2)*gamma(1/nu))/(pi*gamma(1/nu)))**(-1/2)
  skewness <- (lambda*(v**3)/((pi**(3/2))*gamma(1/nu)))*((2**((6+nu)/nu))*(lambda**2)*(gamma(0.5+1/nu)**3)*gamma(1/nu)-
                                                           3*(4**(1/nu))*pi*(1+3*lambda*lambda)*gamma(0.5+1/nu)*gamma(3/nu)+
                                                           4*(pi**(3/2))*(1+lambda*lambda)*gamma(4/nu))
  return(skewness)
}


compare.models <- function(){
  # m1 <- readRDS(paste("../../../results/models/skew-min20-1d.rds", sep=""))
  # # m2 <- readRDS(paste("../../../results/models/baseline-min20-1d.rds", sep=""))
  # m3 <- readRDS(paste("../../../results/models/generror-min20-1d.rds", sep=""))
  # m4 <- readRDS(paste("../../../results/models/skew-generror-min20-1d.rds", sep=""))
  # # comp <- rethinking::compare(m1, m2, m3, refresh = 1)
  # 

  # Parameters for plots
  extra <- 1
  colo <- c("#1b9e77", "#d95f02", "#7570b3")
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed-30.csv", sep = " ", header=T)
  Tinvasive <- read.table("../../../data/properties/codes/neophytes-list_reindexed-30.csv", sep = " ", header=T)
  competitive <- read.table("../../../data/properties/codes/competitive_reindexed-30.csv", sep = " ", header=T)
  
  # extract samples
  post <- extract.samples(m4, pars=c("lambda", "beta"))
  # post <- transform.skew(post)
  # post$beta <- post$betam
  
  # alphas
  mu_lambda <- apply( post$lambda , 2 , mean )
  ci_lambda <- apply( post$lambda , 2 , PI )
  
  # betas
  mu_beta <- apply(post$beta,2,mean)
  ci_beta <- apply(post$beta,2,PI)
  
  corPI <- sapply(1:dim(post$nu)[1], function(x) cor(post$nu[x,], post$nu[x,]))
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
  
  post.m3 <- extract.samples(m3, n = 1000, pars=c("nu", "nu_bar"))
  post.m4 <- extract.samples(m4, n = 1000, pars=c("nu", "nu_bar", "lambda_bar", "lambda"))
  
  kurtosis.median.m3 <- kurtosis.generror(post.m3$nu)
  nu_bar.m3 <- apply(kurtosis.median.m3, 1, mean)

  kurtosis.median.m4 <- kurtosis.skew.generror(post.m4$nu, post.m4$lambda)
  nu_bar.m4 <- apply(kurtosis.median.m4, 1, mean)

  # nu_bar.m3 <- apply(post.m3$nu, 1, median)
  # nu_bar.m4 <- apply(post.m4$nu, 1, median)
  
  n <-  length(nu_bar.m3)
  dat <- data.frame(
    x = c(
      nu_bar.m3,
      nu_bar.m4
    ),
    label = c(
      rep("heavy-tails", n),
      rep("heavy-tails and skewed", n)
    )
  )
  dat$label <- factor(dat$label, levels = c("heavy-tails", "heavy-tails and skewed"))
  
  plist <- list()
  plist[[1]] <- ggplot(dat, aes(x=x, color = label, group = label))+
    geom_vline(xintercept = 0, colour = "#999999", linetype="dashed") +
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    scale_color_manual(values=c("black", "#E69F00")) +
    theme_bw() +
    ggtitle("(a)  average kurtosis of distributions")+
    # scale_fill_grey() +
    ylab("density")+
    xlab("kurtosis")+
    scale_x_continuous(limits = c(-0.75, 0.75))+
    theme(text = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          # legend.position="none",
          legend.position=c(0.25,.85),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(size=10))
  
  post.m1 <- extract.samples(m1, n = 1000, pars=c("lambda", "lambda_bar"))

  skewness.mean.m1 <- skewness.skew(post.m1$lambda)
  lambda_bar.m1 <- apply(skewness.mean.m1, 1, mean)
  
  skewness.mean.m4 <- skewness.skew.generror(post.m4$nu, post.m4$lambda)
  lambda_bar.m4 <- apply(skewness.mean.m4, 1, mean)
  
  # nu_bar.m3 <- apply(post.m3$nu, 1, median)
  # nu_bar.m4 <- apply(post.m4$nu, 1, median)
  
  n <-  length(lambda_bar.m4)
  dat <- data.frame(
    x = c(
      lambda_bar.m1,
      lambda_bar.m4
    ),
    label = c(
      rep("skewed", n),
      rep("heavy-tails and skewed", n)
    )
  )
  
  dat$label <- factor(dat$label, levels = c("skewed", "heavy-tails and skewed"))
  
  plist[[2]] <- ggplot(dat, aes(x=x, color = label, group = label))+
    geom_vline(xintercept = 0, colour = "#999999", linetype="dashed") +
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    scale_color_manual(values=c("black", "#E69F00")) +
    theme_bw() +
    ggtitle("(b)  average skewness of distributions")+
    # scale_fill_grey() +
    ylab(" ")+
    xlab("skewness")+
    scale_x_continuous(limits = c(-0.4, 0.4))+
    theme(text = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          # legend.position="none",
          legend.position=c(0.75,.85),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(size=10))
  p <- grid.arrange(grobs=plist, ncol=2, nrow=1#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
}

distributions.hyperparameters <- function(model=NULL){
  if(is.null(model)){
    model <- readRDS("../../../results/models/skew-traits-min20-1d.rds")
  }
  
  post <- extract.samples(model, n = 1000, pars=c("etasq_l", "sigma_l"))

  n <-  length(post$etasq_l)
  dat <- data.frame(
    x = c(
      post$etasq_l,
      post$sigma_l,
      rexp(n,rate = 1)
      ),
    label = c(
      rep("etasq", n),
      rep("sigma", n),
      rep("prior", n)
    )
  )
  dat$label <- factor(dat$label, levels = c("etasq", "sigma", "prior"))
  
  p <- ggplot(dat, aes(x=x, color = label, group = label))+
    stat_density(geom = "line", position = "identity") +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
    scale_color_manual(values=c("black", "#E69F00", "gray")) +
    theme_bw() +
    ggtitle("(a)  sigma vs etasq")+
    # scale_fill_grey() +
    ylab("density")+
    xlab("kurtosis")+
    scale_x_continuous(limits = c(0, 3))+
    theme(text = element_text(size=10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black"),
          legend.title = element_blank(),
          legend.spacing.x = unit(3, "pt"),
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
          # legend.position="none",
          legend.position=c(0.7,.85),
          legend.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.title = element_text(size=10))
 print(p)  
 # rethinking::pairs(model, pars=c("etasq_l", "sigma_l", "rhosq_l"))
 
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
  post$gamma <- 1/(2*post$gamma)
  
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
    posx = (xlim[2]-xlim[1])*0.8 + xlim[1]
    
    # Beta
    plist[[1]] <- plot.ranking.x(mu_beta, ci_beta, color=colo[1], xlabel="species", ylabel="", ylims=ylim, xlims=xlim, mar=margin(5.5,5.5,5.5,5.5), additional_label = Tind_, posx=posx, posy=c(posy, ylim[2]-(posy-ylim[1])))+
      coord_trans(x="identity")
    
    xlim = c(min(ci_gamma), max(ci_gamma))
    posx = exp((log(xlim[2])-log(xlim[1]))*0.2 + log(xlim[1]))
    
    # Scatter plot
    plist[[2]] <- plot.scatter(mu=mu_beta, variance=mu_gamma, color=colo[3], xlabel="", ylabel=expression(paste("mean ", beta)), mar=margin(5.5,5.5,5.5,5.5), ylims=ylim, xlims=xlim)
    breaks <- round(seq(from=xlim[1], to=xlim[2], length.out = 4))+1
    plist[[2]] <- plist[[2]] + coord_trans(x="log")+
      scale_x_continuous(breaks = breaks, labels = breaks, limits = xlim)
    plist[[2]] <- plist[[2]] + annotate("text", x=posx, y = posy, label=paste0(round(meanPI,2), " ± ", round(sdPI,2)), size=3)
    
    ylim = c(-1, length(mu_gamma)+1)
    
    # Gamma
    plist[[3]] <- plot.ranking.y(mu_gamma, ci_gamma, color=colo[2], xlabel=expression(variance %prop% gamma^{-1}), ylabel="species", xlims=xlim, ylims=ylim, mar=margin(5.5,5.5,5.5,5.5))
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
  post$gamma <- 1/(2*post$gamma)
  
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
  ylab <- expression(paste("mean ", beta))
  
  idx <- sort(mu_beta,index.return=T)$ix
  Tind_ <- as.numeric(as.character(Tind[idx,3]))
  Tind_ <- Tind_[!(is.na(Tind_))]
  Tind_ <- mean(Tind_[1:round(length(Tind_)*0.5)])>mean(Tind_[round(length(Tind_)*0.5):length(Tind_)])
  Tind_ <- c("high elevation", "low elevation")[c(Tind_*1+1, (!(Tind_))*1+1)]
    
  ylim = c(min(mu_beta), max(mu_beta))
  xlim = c(min(mu_gamma), max(mu_gamma))
  posx = exp((log(xlim[2])-log(xlim[1]))*0.2 + log(xlim[1]))
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
  xlab=expression(variance %prop% gamma^{-1})
    
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

plot.actual.data.pairs <- function(model=NULL){
  # Parameters for plots
  extra <- 1
  colo <- c("#e6ab02", "#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    readRDS(paste("../../../results/models/categorical-skew-generror-min20-1d.rds", sep=""))
  }
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed-20.csv", sep = " ")
  Tend <- read.table("../../../data/properties/codes/change-tendency_reindexed-20.csv", sep = " ")
  NEO <- read.table("../../../data/properties/codes/neophytes-list_reindexed-20.csv", sep = " ")
  
  params <- c("alpha", "gamma", "beta", "lambda", "nu")
  label <- c("amplitude", "variance", "mean", "skewness", "kurtosis")
  trans <- c("log", "log", "identity", "identity", "identity")
  
  com <- combn(c(1:length(params)), 2)
  # extract samples
  post <- extract.samples(model, n = 1000, pars=params) 
  post$gamma <- 1/(2*post$gamma)
  post$alpha <- exp(-post$alpha)
  
  plist <- list()
  hlay <- matrix(NA, nrow=length(params), ncol=length(params))
  correlations <- matrix(NA, nrow=length(params), ncol=length(params))
  
  for(i in 1:length(params)){
    correlations[i,i] <- label[i]
  }
  
  for( i in 1:ncol(com)){
    hlay[com[1,i], com[2,i]] <- i
    # betas
    mu_beta <- apply(post[[params[com[1,i]]]],2,mean)
    ci_beta <- apply(post[[params[com[1,i]]]],2,PI)
    transformation_y <- trans[com[1, i]]
    label_y <- ifelse(is.na(hlay[com[1,i], com[2,i]-1]), label[com[1,i]], "")
    
    # gammas
    mu_gamma <- apply(post[[params[com[2,i]]]],2,mean)
    ci_gamma <- apply(post[[params[com[2,i]]]],2,PI)
    transformation_x <- trans[com[2, i]]
    label_x <- ifelse(is.na(hlay[com[1,i], com[2,i]-1]), label[com[2,i]], "")
    
    corPI <- sapply(1:dim(post[[params[com[1,i]]]])[1], function(x) cor(post[[params[com[1,i]]]][x,NEO$V3==1], post[[params[com[2,i]]]][x,NEO$V3==1]))
    meanPI <- mean(corPI)
    sdPI <- sd(corPI)
    corPIci_NEO <- paste(round(c(meanPI, sdPI), 2), collapse = " ± ")
    correlations[com[1,i], com[2,i]] <- corPIci_NEO

    corPI <- sapply(1:dim(post[[params[com[1,i]]]])[1], function(x) cor(post[[params[com[1,i]]]][x,NEO$V3!=1], post[[params[com[2,i]]]][x,NEO$V3!=1]))
    meanPI <- mean(corPI)
    sdPI <- sd(corPI)
    corPIci_NATIVE <- paste(round(c(meanPI, sdPI), 2), collapse = " ± ")
    correlations[com[2,i],com[1,i]] <- corPIci_NATIVE
    
    posx <- ifelse(meanPI>0, 0.85, 0.15)
    
    if(transformation_x=="log"){
      posx <- exp((log(max(mu_gamma)) - log(min(mu_gamma)))*posx + log(min(mu_gamma)))
    }else{
      posx <- (max(mu_gamma) - min(mu_gamma))*posx + min(mu_gamma)
    }
    if(transformation_y=="log"){
      posy <- exp((log(max(mu_beta)) - log(min(mu_beta)))*0.15 + log(min(mu_beta)))
    }else{
      posy <- (max(mu_beta) - min(mu_beta))*0.15 + min(mu_beta)
    }

    labels = ifelse(NEO$V3==1, "neophyte", "native")
    alphas = ifelse(NEO$V3==1, 0.7, 0.6)
    color <- c("#e7298a",colo[4])
    
    p <- plot.scatter2(mu=mu_beta, variance=mu_gamma, label = labels, color=c(color[2],color[1]), alpha = alphas, xlabel=label_x, ylabel=label_y, mar=margin(5.5,5.5,5.5,5.5))
    # Scatter plot
    p <- p + coord_trans(x=transformation_x, y=transformation_y, clip = "off")
    p <- p + annotate("text", x=posx, y = posy, label=paste0("\n", corPIci_NEO), size=2, color="#e7298a")
    p <- p + annotate("text", x=posx, y = posy, label=paste0(corPIci_NATIVE, "\n"), size=2)
      # scale_x_continuous(breaks = breaks_x, labels = breaks_x, limits = xlim)+
      # scale_y_continuous(breaks = breaks_y, labels = breaks_y, limits = ylim)
    p <- p+theme(legend.position = "none", plot.title = element_blank())
    plist[[i]] <- p
  }
  
  # tab <- ggtexttable(correlations, rows = NULL,cols = NULL,
  #                    theme = ttheme("blank",  base_size = 9))

  
  hlay <- hlay[1:(nrow(hlay)-1), 2:ncol(hlay)]

  grobs <- list()
  widths <- list()
  heights <- list()
  
  for (k in 1:length(plist)){
    grobs[[k]] <- ggplotGrob(plist[[k]])
    widths[[k]] <- grobs[[k]]$widths[2:5]
    heights[[k]] <- grobs[[k]]$heights[2:5]
  }    
  maxwidth <- do.call(grid::unit.pmax, widths)
  maxheight <- do.call(grid::unit.pmax, heights)
  for (k in 1:length(grobs)){
    grobs[[k]]$widths[2:5] <- as.list(maxwidth)
    grobs[[k]]$heights[2:5] <- as.list(maxheight)
  }
  
  # grobs[[length(grobs)+1]] <- ggplotGrob(tab)
  # hlay[(ncol(hlay)/2+1):ncol(hlay), 1:(ncol(hlay)/2)] <- length(grobs)
  # 
  
  p <- grid.arrange(grobs=grobs, ncol=length(params)-1, nrow=length(params)-1, layout_matrix=hlay#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
}

plot.actual.data.distributions <- function(model=NULL){
  # Parameters for plots
  meansample=T
  extra <- 1
  colo <- c("#e6ab02", "#1b9e77", "#d95f02", "#7570b3")
  
  # Load the data if not added  
  if(is.null(model)){
    model <- readRDS(paste("../../../results/models/categorical-skew-generror-min20-1d.rds", sep=""))
  }
  
  #indicator values
  Tind <- read.table("../../../data/properties/codes/temperature_indicator_reindexed-20.csv", sep = " ")
  Tend <- read.table("../../../data/properties/codes/change-tendency_reindexed-20.csv", sep = " ")
  NEO <- read.table("../../../data/properties/codes/neophytes-list_reindexed-20.csv", sep = " ")
  
  params <- c("alpha", "gamma", "beta", "lambda", "nu")
  label <- c("amplitude", "variance", "mean", "skewness", "kurtosis")
  trans <- c("log", "log", "identity", "identity", "identity")
  
  # extract samples
  post <- extract.samples(model, n = 1000, pars=params) 
  post$gamma <- 1/(2*post$gamma)
  post$alpha <- exp(-post$alpha)
  
  plist <- list()
  for(i in 1:length(params)){

    label_y <- ifelse(i==1, "density", "")
    transformation_x <- trans[i]
    label_x <- label[i]
    
    # z <- matrix(NA, nrow(post[[params[i]]]), sum(NEO$V3==1))
    # for(j in 1:nrow(post[[params[i]]])){z[j,] <- post[[params[i]]][j,(1:length(NEO$V3)) %in% sample((1:length(NEO$V3))[NEO$V3!=1], sum(NEO$V3==1))]}

    x <- post[[params[i]]][,NEO$V3!=1]
    y <- post[[params[i]]][,NEO$V3==1]    
    if(meansample){
      x <- apply(x, 1, median)
      y <- apply(y, 1, median)
      # z <- apply(z, 1, mean)
    }else{
      x <- as.vector(x)
      y <- as.vector(y)
      # z <- as.vector(z)
    }
    
    dat <- data.frame(
      x = c(
        # z,
        x,
        y
      ),
      label = c(
        # rep("random", length(z)),
        rep("native", length(x)),
        rep("neophytes", length(y))
      )
    )
    
    # dat$label <- factor(dat$label, levels = c("random","native", "neophytes"))
    # colorins <- c("gray",colo[4], "#e7298a")
    colorins <- c(colo[4], "#e7298a")
    dat$label <- factor(dat$label, levels = c("native", "neophytes"))
    
    p <- ggplot(dat, aes(x=x, color = label, group = label, fill=label))+
      # geom_vline(xintercept = 0, colour = "#999999", linetype="dashed") +
      geom_density(position = "identity", alpha=0.2) +    # geom_vline(aes(xintercept=1), color="black", linetype="dashed", size=0.5) +
      scale_color_manual(values=colorins) +
      scale_fill_manual(values=colorins) +
      theme_bw() +
      # scale_fill_grey() +
      ylab(label_y)+
      xlab(label_x)+
      theme(text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(colour = "black"),
            legend.title = element_blank(),
            legend.spacing.x = unit(3, "pt"),
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            # legend.position="none",
            legend.position=c(0.75,.85),
            legend.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            plot.title = element_text(size=10))

    # p <- p + coord_trans(x=transformation_x, clip = "off")
    
    if(i!=length(params)){p <- p+theme(legend.position = "none", plot.title = element_blank())}
    plist[[i]] <- p
  }
  

  grobs <- list()
  widths <- list()
  heights <- list()
  
  for (k in 1:length(plist)){
    grobs[[k]] <- ggplotGrob(plist[[k]])
    widths[[k]] <- grobs[[k]]$widths[2:5]
    heights[[k]] <- grobs[[k]]$heights[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  maxheight <- do.call(grid::unit.pmax, heights)
  for (k in 1:length(grobs)){
    grobs[[k]]$widths[2:5] <- as.list(maxwidth)
    grobs[[k]]$heights[2:5] <- as.list(maxheight)
  }

  p <- grid.arrange(grobs=grobs, ncol=3, nrow=2#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  print(p)
  
}

