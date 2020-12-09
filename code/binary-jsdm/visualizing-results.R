check_results_latest <- function(d, model){
  alphas <- precis(model, pars = "alpha", depth=2)
  betas <- precis(model, pars = "beta", depth=3)
  sigmas <- precis(model, pars = "gamma", depth=3)
  
  # sigmas <- precis(model, pars = "sigma_beta", depth=3)
  N <- length(unique(d$dataset$id))
  alpha_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$alpha[1])
  beta1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta1[1])
  beta2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$beta2[1])
  sigma1_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta1[1])
  sigma2_r <- sapply(1:N, function(x) d$dataset[d$dataset$id==x,]$sigma_beta2[1])
  d_alpha <- data.frame(N=1:N, id= c(rep("real", length(alpha_r)), rep("estimated", length(alphas$mean))),value=c(alpha_r,alphas$mean) , sd=c(rep(0,length(alpha_r)),alphas$sd))
  d_beta1 <- data.frame(N=1:N, id= c(rep("real", length(beta1_r)), rep("estimated", length(betas[1:N,]$mean))),value=c(beta1_r,betas[1:N,]$mean) , sd=c(rep(0,length(beta1_r)),betas[1:N,]$sd))
  d_beta2 <- data.frame(N=1:N, id= c(rep("real", length(beta2_r)), rep("estimated", length(betas[(N+1):(N+N),]$mean))),value=c(beta2_r,betas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(beta2_r)),betas[(N+1):(N+N),]$sd))
  d_sigma1 <- data.frame(N=1:N, id= c(rep("real", length(sigma1_r)), rep("estimated", length(sigmas[1:N,]$mean))),value=c(sigma1_r,sigmas[1:N,]$mean) , sd=c(rep(0,length(sigma1_r)),sigmas[1:N,]$sd))
  d_sigma2 <- data.frame(N=1:N, id= c(rep("real", length(sigma2_r)), rep("estimated", length(sigmas[(N+1):(N+N),]$mean))),value=c(sigma2_r,sigmas[(N+1):(N+N),]$mean) , sd=c(rep(0,length(sigma2_r)),sigmas[(N+1):(N+N),]$sd))
  
  
  p1 <- ggplot(d_alpha, aes(x=N, y=value, group=id, color=id)) + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
  p2 <- ggplot(d_beta1, aes(x=N, y=value, group=id, color=id)) + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
  p3 <- ggplot(d_beta2, aes(x=N, y=value, group=id, color=id)) + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
  p4 <- ggplot(d_sigma1, aes(x=N, y=value, group=id, color=id)) + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
  p5 <- ggplot(d_sigma2, aes(x=N, y=value, group=id, color=id)) + 
    geom_pointrange(aes(ymin=value-sd, ymax=value+sd)) + theme_linedraw()
  figure <- grid.arrange(p1, p2, p3,
                         ncol = 1, nrow = 3)
  print(figure)
  figure <- grid.arrange(p4, p5,
                         ncol = 1, nrow = 2)
  print(figure)
}