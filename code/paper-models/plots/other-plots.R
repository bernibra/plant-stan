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


