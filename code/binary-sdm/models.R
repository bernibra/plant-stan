# SDM with binary data

# Model 1.0 is a simple logistic regression where I do not generate log-odds values so I can't calculate WAIC directly
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

model1.01 <- "
data{
  int N;
  int K;
  int obs[N];
  matrix[N,K] bio;
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
  
  p = alpha + bio * beta;
  p = inv_logit(p);  

  obs ~ binomial(1, p);
}
"

# Model 1.1 is the same as Model 1.0 but I do calculate the log-likelihood so that I can estimate WAIC and LOO scores
model1.1 <- "
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
  alpha ~ normal(0,1.3);
  beta ~ normal(0,1.3);
  
  for (i in 1:N) {
     p[i] = alpha;
     for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j];
     }
     p[i] = inv_logit(p[i]);
  }
  
  obs ~ binomial(1, p);
}
generated quantities{
    vector[N] log_lik;
    vector[N] p;
    
    for (i in 1:N) {
      p[i] = alpha;
      for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j];
      }
      p[i] = inv_logit(p[i]);
      log_lik[i] = binomial_lpmf(obs[i] | 1, p[i]);
    }
}
"
