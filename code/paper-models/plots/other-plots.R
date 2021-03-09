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


