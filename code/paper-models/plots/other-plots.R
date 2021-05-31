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

distribution.type <- function(){
  x <- seq(-5, 5, length.out = 2000)
  generrskew <- function(x, a, mu, sigma, lambda, p){
    v <- sqrt((pi*gamma(1/p))/(pi*(1+3*lambda**2)*gamma(3/p)-(16**(1/p))*lambda*lambda*(gamma(1/2+1/p)**2)*gamma(1/p)))
    m <- (2**(2/p))*v*sigma*lambda*gamma(1/2+1/p)/sqrt(pi)
    return((p*sigma/(2*v*gamma(1/p)))*exp(-(sigma*abs(x-mu+m)/(v*(1+lambda*sign(x-mu+m))))**p))
  }
  
  plist <- list()
  
  for(i in 1:4){
    if(i==1){
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,0.6,1.5)
      lambda <- rep(0,3)
      nu <- rep(2,3)
      title <- "(a)  normal"
      labs <- c(expression(gamma==1), expression(gamma<1), expression(gamma>1))
    }else if (i==2){
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1,1)
      lambda <- rep(0,3)
      nu <- c(2,1.2,9)
      labs <- c(expression(nu==2), expression(nu<2), expression(nu>2))
      title <- "(b)  heavy-tails"
    }else if (i==3){
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1,1)
      lambda <- c(0,-0.5, 0.5)
      nu <- rep(2,3)
      labs <- c(expression(lambda==0), expression(lambda<0), expression(lambda>0))
      title <- "(c)  skewed"
    }else{
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1, 1)
      lambda <- c(0, -0.75, 0.75)
      nu <- c(2,1.2,9)
      title <- "(d)  heavy-tail and skewed"
      labs <- c(expression(list(lambda==0, nu==2)), expression(list(lambda<0, nu<2)), expression(list(lambda>0, nu>2)))
      
    }
    dat <- data.frame(x=x, y=generrskew(x, a[1], mu[1], sigma[1], lambda[1], nu[1]), color=rep(1, length(x)))
    for (j in 2:length(a)){
      dat <- rbind(dat, data.frame(x=x, y=generrskew(x, a[j], mu[j], sigma[j], lambda[j], nu[j]), color=rep(j, length(x))))
    }
    plist[[i]] <- ggplot() +
      geom_line(data = dat[dat$color==2,], aes(x=x, y=y,colour="2"), size=0.7, alpha=0.8)+
      geom_line(data = dat[dat$color==3,], aes(x=x, y=y,colour="3"), size=0.7, alpha=0.8)+
      geom_line(data = dat[dat$color==1,], aes(x=x, y=y,colour="1"), size=0.7, alpha=0.8)+
      theme_bw() +
      scale_color_manual(labels = labs, values=c("#7570b3", "#d95f02", "#1b9e77"))+
      ggtitle(title)+
      # scale_fill_grey() +
      ylab(" ")+
      xlab("elevation")+
      scale_x_continuous(limits = c(-4, 4))+
      theme(text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(colour = "black"),
            legend.title = element_blank(),
            legend.spacing.x = unit(3, "pt"),
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            # legend.position="none",
            legend.position=c(0.3,.85),
            legend.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            plot.title = element_text(size=10))
  }

  p <- grid.arrange(grobs=plist, ncol=4, nrow=1#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
}

distribution.type_2 <- function(){
  x <- seq(-5, 5, length.out = 2000)
  generrskew <- function(x, a, mu, sigma, lambda, p, normalized=T){
    v <- sqrt((pi*gamma(1/p))/(pi*(1+3*lambda**2)*gamma(3/p)-(16**(1/p))*lambda*lambda*(gamma(1/2+1/p)**2)*gamma(1/p)))
    m <- (2**(2/p))*v*sigma*lambda*gamma(1/2+1/p)/sqrt(pi)
    if(normalized){
      return((p*sigma/(2*v*gamma(1/p)))*exp(-(sigma*abs(x-mu+m)/(v*(1+lambda*sign(x-mu+m))))**p))
    }else{
      return(exp(-(sigma*abs(x-mu+m)/(v*(1+lambda*sign(x-mu+m))))**p))
    }
  }
  
  plist <- list()
  
  for(i in 1:6){
    normalization <- rep(1,3)
    if(i==2){
      a <- rep(0,3)
      mu <- c(-1.3,0,1.3)
      sigma <- c(1,1,1)
      lambda <- rep(0,3)
      nu <- rep(2,3)
      title <- "(b)  mean"
      labs <- c(expression(beta<0), expression(beta==0), expression(beta>0))
      xlab_ <- ""
      ylab_ <- ""
      normalized <- T
    }else if(i==1){
      normalization <- c(exp(-1),exp(-0.5),exp(0))
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1,1)
      lambda <- rep(0,3)
      nu <- rep(2,3)
      xlab_ <- ""
      ylab_ <- ""
      title <- "(a)  amplitude"
      labs <- c(expression(alpha==1), expression(alpha==0.5), expression(alpha==0))
      normalized <- F
    }else if(i==3){
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,0.6,1.5)
      lambda <- rep(0,3)
      nu <- rep(2,3)
      xlab_ <- ""
      ylab_ <- ""
      title <- "(c)  variance"
      labs <- c(expression(gamma==1), expression(gamma<1), expression(gamma>1))
      normalized <- T
    }else if (i==4){
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1,1)
      lambda <- rep(0,3)
      nu <- c(2,1.2,9)
      labs <- c(expression(nu==2), expression(nu<2), expression(nu>2))
      xlab_ <- ""
      ylab_ <- ""
      title <- "(d)  kurtosis"
      normalized <- T
    }else if (i==5){
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1,1)
      lambda <- c(0,-0.5, 0.5)
      nu <- rep(2,3)
      labs <- c(expression(lambda==0), expression(lambda<0), expression(lambda>0))
      xlab_ <- "elevation"
      ylab_ <- ""
      title <- "(e)  skewness"
      normalized <- T
    }else{
      a <- rep(0,3)
      mu <- rep(0,3)
      sigma <- c(1,1, 1)
      lambda <- c(0, -0.75, 0.75)
      nu <- c(2,1.2,9)
      xlab_ <- "elevation"
      ylab_ <- ""
      title <- "(f)  kurtosis and skewness"
      labs <- c(expression(list(lambda==0, nu==2)), expression(list(lambda<0, nu<2)), expression(list(lambda>0, nu>2)))
      normalized <- T
    }
    dat <- data.frame(x=x, y=normalization[1]*generrskew(x, a[1], mu[1], sigma[1], lambda[1], nu[1], normalized = normalized), color=rep(1, length(x)))
    for (j in 2:length(a)){
      dat <- rbind(dat, data.frame(x=x, y=normalization[j]*generrskew(x, a[j], mu[j], sigma[j], lambda[j], nu[j], normalized = normalized), color=rep(j, length(x))))
    }
    plist[[i]] <- ggplot() +
      geom_line(data = dat[dat$color==2,], aes(x=x, y=y,colour="2"), size=0.7, alpha=0.8)+
      geom_line(data = dat[dat$color==3,], aes(x=x, y=y,colour="3"), size=0.7, alpha=0.8)+
      geom_line(data = dat[dat$color==1,], aes(x=x, y=y,colour="1"), size=0.7, alpha=0.8)+
      theme_bw() +
      scale_color_manual(labels = labs, values=c("#7570b3", "#d95f02", "#1b9e77"))+
      ggtitle(title)+
      # scale_fill_grey() +
      ylab(ylab_)+
      xlab(xlab_)+
      scale_x_continuous(limits = c(-4, 4), breaks=c(-3,0,3))+
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+
      theme(text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(colour = "black"),
            legend.title = element_blank(),
            legend.spacing.x = unit(3, "pt"),
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            # legend.position="none",
            legend.position=c(0.3,.85),
            legend.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            legend.key.width = unit(12,"pt"),
            plot.title = element_text(size=10))
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
  
  p <- grid.arrange(grobs=plist, ncol=2, nrow=3#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
}

distribution.type_3 <- function(expressions=F){
  x <- seq(-4, 4, length.out = 2000)
  generrskew <- function(x, a, mu, sigma, lambda, p, normalized=T){
    v <- sqrt((pi*gamma(1/p))/(pi*(1+3*lambda**2)*gamma(3/p)-(16**(1/p))*lambda*lambda*(gamma(1/2+1/p)**2)*gamma(1/p)))
    m <- (2**(2/p))*v*sigma*lambda*gamma(1/2+1/p)/sqrt(pi)
    if(normalized){
      return((p*sigma/(2*v*gamma(1/p)))*exp(-(sigma*abs(x-mu+m)/(v*(1+lambda*sign(x-mu+m))))**p))
    }else{
      return(exp(-(sigma*abs(x-mu+m)/(v*(1+lambda*sign(x-mu+m))))**p))
    }
  }
  
  plist <- list()
  idx <- c(3,1,2)
  
  for(i in 1:4){
    normalization <- rep(1,3)
    if(i==1){
      normalization <- c(exp(-0.5),exp(0),exp(-1))[idx]
      a <- rep(0,3)[idx]
      mu <- c(0,2,-2)[idx]
      sigma <- c(1,2.5, 0.5)[idx]
      lambda <- rep(0,3)[idx]
      nu <- rep(2,3)[idx]
      xlab_ <- "elevation"
      ylab_ <- "probability"
      title <- "(a)  normal"
      labs <- c(expression(gamma==1), expression(gamma<1), expression(gamma>1))[idx]
      normalized <- F
    }else if (i==2){
      normalization <- c(exp(-0.5),exp(0),exp(-1))[idx]
      a <- rep(0,3)[idx]
      mu <- c(0,2,-2)[idx]
      sigma <- c(1,1,1)[idx]
      lambda <- rep(0,3)[idx]
      nu <- c(2,1.2,9)[idx]
      labs <- c(expression(nu==2), expression(nu<2), expression(nu>2))[idx]
      xlab_ <- "elevation"
      ylab_ <- ""
      title <- "(b)  fat-tailed"
      normalized <- F
    }else if (i==3){
      normalization <- c(exp(-0.5),exp(0),exp(-1))[idx]
      a <- rep(0,3)[idx]
      mu <- c(0,2,-2)[idx]
      sigma <- c(1,1,1)[idx]
      lambda <- c(0,-0.5, 0.5)[idx]
      nu <- rep(2,3)[idx]
      labs <- c(expression(lambda==0), expression(lambda<0), expression(lambda>0))[idx]
      xlab_ <- "elevation"
      ylab_ <- ""
      title <- "(c)  skewed"
      normalized <- F
    }else{
      normalization <- c(exp(-0.5),exp(0),exp(-1))[idx]
      a <- rep(0,3)[idx]
      mu <- c(0,2,-2)[idx]
      sigma <- c(1,1, 1)[idx]
      lambda <- c(0, -0.75, 0.75)[idx]
      nu <- c(2,1.2,9)[idx]
      xlab_ <- "elevation"
      ylab_ <- ""
      title <- "(d)  fat-tailed and skewed"
      labs <- c(expression(list(lambda==0, nu==2)), expression(list(lambda<0, nu<2)), expression(list(lambda>0, nu>2)))[idx]
      normalized <- F
    }
    dat <- data.frame(x=x, y=normalization[1]*generrskew(x, a[1], mu[1], sigma[1], lambda[1], nu[1], normalized = normalized), color=rep(1, length(x)))
    for (j in 2:length(a)){
      dat <- rbind(dat, data.frame(x=x, y=normalization[j]*generrskew(x, a[j], mu[j], sigma[j], lambda[j], nu[j], normalized = normalized), color=rep(j, length(x))))
    }
    if(expressions){
      labelsx <- c(expression(beta[1]), expression(beta[2]) , expression(beta[3]))
      labelsy <- c(expression(p[1]), expression(p[2]) , expression(p[3]))
      breaksy <- c(exp(-1),exp(-0.5),exp(0))
      breaksx <- c(-2,0,2)
    }else{
      labelsy <- c(exp(-1),exp(-0.5),exp(0))
      labelsx <- c(-2,0,2)
      breaksy <- c(exp(-1),exp(-0.5),exp(0))
      breaksx <- c(-2,0,2)
    }
    d=data.frame(x=c(-5,-5,-5,breaksx[1],breaksx[2], breaksx[3]), y=c(breaksy[1],breaksy[2], breaksy[3], -0.5, -0.5, -0.5),
                 vx=c(breaksx[1],breaksx[2], breaksx[3], breaksx[1],breaksx[2], breaksx[3]), vy=c(breaksy[1],breaksy[2], breaksy[3],breaksy[1],breaksy[2], breaksy[3]))

    plist[[i]] <- ggplot() +
      # geom_segment(data=d, mapping=aes(x=x, y=y, xend=vx, yend=vy), colour="gray80", linetype=3) +
      geom_line(data = dat[dat$color==2,], aes(x=x, y=y,colour="2"), size=0.7, alpha=0.8)+
      geom_line(data = dat[dat$color==3,], aes(x=x, y=y,colour="3"), size=0.7, alpha=0.8)+
      geom_line(data = dat[dat$color==1,], aes(x=x, y=y,colour="1"), size=0.7, alpha=0.8)+
      theme_bw() +
      scale_color_manual(labels = labs, values=c("#7570b3", "#d95f02", "#1b9e77")[idx])+
      ggtitle(title)+
      # scale_fill_grey() +
      ylab(ylab_)+
      xlab(xlab_)+
      coord_cartesian(ylim=c(0, max(breaksy)), xlim=c(-4, 4)) +
      scale_x_continuous(breaks=breaksx, labels = labelsx)+
      scale_y_continuous(breaks=breaksy, labels = labelsy)+
      theme(text = element_text(size=10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(colour = "black"),
            legend.title = element_blank(),
            legend.spacing.x = unit(3, "pt"),
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
            # legend.position="none",
            legend.position=c(0.19,.82),
            legend.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            legend.key.width = unit(12,"pt"),
            plot.title = element_text(size=10))
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
  
  p <- grid.arrange(grobs=plist, ncol=4, nrow=1#, vp=viewport(width=1, height=1, clip = TRUE),
                    # top=textGrob("First axis", rot = 0, vjust = 0.9,gp=gpar(fontsize=10)),
  )
  
}
distribution.type_3(expressions = T)

