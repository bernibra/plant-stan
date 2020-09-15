# SDM with binary data

# link1 is the link function for model1
link1 <- "
data{
  int N;
  int N_samples;
  int K;
  real bio[N,K];
  vector[N_samples] alpha;
  real beta[N_samples, K];
}
parameters{
}
model{
}
generated quantities {
  real p[N_samples, N];

  for(i in 1:N){
      for(j in 1:N_samples){
          p[j, i] = alpha[j];
          for(k in 1:K){
              p[j, i] = p[j, i] + beta[j, k] * bio[i, k];
          }
          p[j, i] = inv_logit(p[j, i]);
      }
  }
}
"
# #### Instead of writing link functions in R, I could instead do something along the lines of what I wrote below.
# # This crashes in R and I am not very sure why... I should probably debug it but I find writing the link function in
# # R a bit more intuitive and it can be made quite efficient.
# post <- extract(mfit_1.1)
# pred <- stan(model_code = link1,
#              data = list(N=length(obs),
#                          K=length(variables),
#                          N_samples = length(post$alpha),
#                          bio=bio,
#                          alpha=post$alpha, 
#                          beta=post$beta),
#              chains = 1, iter = 100,
#              algorithm = "Fixed_param")
# 
