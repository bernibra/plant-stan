# SDM with binary data

# Model 1.0 is a simple poisson regression to predict the number of species in a particular place
model1.0 <- "
data{
  int N;
  int K;
  int obs[N];
  real bio[N,K];
}
parameters{
  real alpha;
  vector[K] beta;
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  
  for (i in 1:N) {
     p[i] = alpha;
     for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j];
     }
     p[i] = inv_logit(p[i]);
  }
  
  obs ~ binomial(1, p);
}
"
