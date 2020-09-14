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

# Model 1.01 is exactly the same as Model 1.0, I just wanted to write things in matrix form.
# It turns out it is not much faster, so it makes more sense to just use the above notation
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

# Model 2.0 is a logistic regression with quadratic terms for the different predictors
# where I do not generate log-odds values so I can't calculate WAIC directly
model2.0 <- "
data{
  int N;
  int K;
  int obs[N];
  real bio[N,K];
  real bio2[N,K];
}
parameters{
  real alpha;
  vector[K] beta;
  vector[K] beta2;
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  beta2 ~ normal(0,0.5);
  
  for (i in 1:N) {
     p[i] = alpha;
     for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j] + beta2[j] * bio2[i, j];
     }
     p[i] = inv_logit(p[i]);
  }
  
  obs ~ binomial(1, p);
}
"


# Model 2.1 is a logistic regression with quadratic terms for the different predictors
# where I do generate log-odds values so I can calculate WAIC directly
model2.1 <- "
data{
  int N;
  int K;
  int obs[N];
  real bio[N,K];
  real bio2[N,K];
}
parameters{
  real alpha;
  vector[K] beta;
  vector[K] beta2;
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  beta2 ~ normal(0,0.5);
  
  for (i in 1:N) {
     p[i] = alpha;
     for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j] + beta2[j] * bio2[i, j];
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
          p[i] = p[i] + beta[j] * bio[i, j] + beta2[j] * bio2[i, j];
      }
      p[i] = inv_logit(p[i]);
      log_lik[i] = binomial_lpmf(obs[i] | 1, p[i]);
    }
}
"

# Model 3.0 is a zero-inflated logistic regression
# where I do not generate log-odds values so I can't calculate WAIC directly
model3.0 <- "
data{
  int N;
  int K;
  int obs[N];
  real bio[N,K];
}
parameters{
  real alpha1;
  real alpha2;
  vector[K] beta;
}
model{
  vector[N] p1;
  vector[N] p2;

  //priors
  alpha1 ~ normal(0,1);
  alpha2 ~ normal(0,1);
  beta ~ normal(0,1);
  
  for (i in 1:N) {
     p1[i] = alpha1;
     p2[i] = alpha2;
     for (j in 1:K){
          p1[i] = p1[i] + beta[j] * bio[i, j];
     }
     p1[i] = inv_logit(p1[i]);
     p2[i] = inv_logit(p2[i]);
     
     if (obs[i] == 0)
        target += log_sum_exp(bernoulli_lpmf(0 | p1[i]), 
                              bernoulli_lpmf(1 | p1[i])
                               + bernoulli_lpmf(0 | p2[i]));
     else
        target +=  bernoulli_lpmf(1 | p1[i]) + bernoulli_lpmf(1 | p2[i]);
  }
}
"

# Model 4.0 is a zero-inflated logistic regression with quadratic terms
# where I do not generate log-odds values so I can't calculate WAIC directly

model4.0 <- "
data{
  int N;
  int K;
  int obs[N];
  real bio[N,K];
  real bio2[N,K];
}
parameters{
  real alpha1;
  real alpha2;
  vector[K] beta;
  vector[K] beta2;
}
model{
  vector[N] p1;
  vector[N] p2;

  //priors
  alpha1 ~ normal(0,1.3);
  alpha2 ~ normal(-1,0.3);
  beta ~ normal(0,0.5);
  beta2 ~ normal(0,0.25);
  
  for (i in 1:N) {
     p1[i] = alpha1;
     p2[i] = alpha2;
     for (j in 1:K){
          p1[i] = p1[i] + beta[j] * bio[i, j] + beta2[j] * bio2[i, j];
     }
     p1[i] = inv_logit(p1[i]);
     p2[i] = inv_logit(p2[i]);
     
     if (obs[i] == 0)
        target += log_sum_exp(bernoulli_lpmf(0 | p1[i]), 
                              bernoulli_lpmf(1 | p1[i])
                               + bernoulli_lpmf(0 | p2[i]));
     else
        target +=  bernoulli_lpmf(1 | p1[i]) + bernoulli_lpmf(1 | p2[i]);
  }
}
"

