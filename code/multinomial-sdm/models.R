# SDM with binary data

# Model 1.0 is a simple logistic regression where I do not generate log-odds values so I can't calculate WAIC directly
model1.0 <- "
data{
  int N;
  int K;
  int L;
  int obs[N];
  real bio[N,K];
}
parameters{
  ordered[L] alpha;
  vector[K] beta;
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(0,1.5);
  beta ~ normal(0,1);
  
  for (i in 1:N) {
     p[i] = 0;
     for (j in 1:K){
          p[i] = p[i] + beta[j] * bio[i, j];
     }
     obs[i] ~ ordered_logistic( p[i] , alpha );
  }
  
}
"

# Model 2.0 is a zero-inflated logistic regression
# where I do not generate log-odds values so I can't calculate WAIC directly
model2.0 <- "
data{
  int N;
  int K;
  int L;
  int obs[N];
  real bio[N,K];
}
parameters{
  real gamma;
  ordered[L] alpha;
  vector[K] beta;
}
model{
  vector[N] p1;
  vector[N] p2;

  //priors
  gamma ~ normal(0,1.3);
  alpha ~ normal(0,1.5);
  beta ~ normal(0,0.5);
  
  for (i in 1:N) {
     p1[i] = gamma;
     p2[i] = 0;
     for (j in 1:K){
          p2[i] = p2[i] + beta[j] * bio[i, j];
     }
     p1[i] = inv_logit(p1[i]);

     if (obs[i] == 1)
        target += log_sum_exp(bernoulli_lpmf(0 | p1[i]),
                              bernoulli_lpmf(1 | p1[i])
                               + ordered_logistic_lpmf(obs[i] | p2[i] , alpha ));
     else
        target +=  bernoulli_lpmf(1 | p1[i]) + ordered_logistic_lpmf(obs[i] | p2[i] , alpha );
  }
}
"

# Model 3.0 is a ordered categorical regression with gaussian RBF, where I do not generate log-odds values so I can't calculate WAIC directly
model3.0 <- "
functions{
  real radial(real x, real mu, real epsilon, real beta) {
    real K;
    K = x - mu;
    K = K/epsilon;
    K = pow(K, 2);
    K = exp(-0.5*K);
    K = beta * K;
    return K;
  }
}
data{
  int N;
  int K;
  int L;
  int obs[N];
  real bio[N,K];
}
parameters{
  ordered[L] alpha;
  vector[K] beta;
  vector[K] mu_bar;
  real<lower=0> epsilon[K];
}
model{
  vector[N] p;

  //priors
  alpha ~ normal(0,1.5);
  beta ~ normal(0,1);
  mu_bar ~ normal(0,1);
  epsilon ~ exponential( 1 );
  
  for (i in 1:N) {
     p[i] = 0;
     for (j in 1:K){
          p[i] = p[i] + radial(bio[i, j], mu_bar[j], epsilon[j], beta[j]);
     }
     obs[i] ~ ordered_logistic( p[i] , alpha );
  }
}
"

