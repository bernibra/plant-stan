# SDM with binary data

# Model 1.0 is a simple poisson regression to predict the number of species in a particular place
# without predictors
model1.0 <- "
data{
  int N;
  int K;
  int sp[N];
}
parameters{
  real alpha;
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(3,0.5);

  for (i in 1:N) {
     p[i] = alpha;
     p[i] = exp(p[i]);
  }
  
  sp ~ poisson(p);
}
generated quantities{
    vector[N] log_lik;
    vector[N] p;

    for ( i in 1:N ) {
         p[i] = alpha;
         p[i] = exp(p[i]);
         log_lik[i] = poisson_lpmf(sp[i] | p[i]);
    }
}
"

# Model 1.1 is a simple poisson regression to predict the number of species in a particular place 
# with environmental predictors
model1.1 <- "
data{
  int N;
  int K;
  int sp[N];
  real bio[N,K];
}
parameters{
  real alpha;
  vector[K] beta;
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(3,0.5);
  beta ~ normal(0,0.5);
  
  for (i in 1:N) {
     p[i] = alpha;
     for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j];
     }
     p[i] = exp(p[i]);
  }
  
  sp ~ poisson(p);
}
generated quantities{
    vector[N] log_lik;
    vector[N] p;

    for ( i in 1:N ) {
         p[i] = alpha;
         for (j in 1:K){
             p[i] = p[i] + beta[j] * bio[i, j];
         }
         p[i] = exp(p[i]);
         log_lik[i] = poisson_lpmf(sp[i] | p[i]);
    }
}
"
