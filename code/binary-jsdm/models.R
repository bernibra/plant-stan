# JSDM with binary data

# Model 1.0 is a simple logistic regression where I do not generate log-odds values so I can't calculate WAIC directly
model1.0 <- "
data{
    int N;
    int K;
    int L;
    int obs[N];
    real bio[N,K];
    int id[N];
}
parameters{
    vector[L] zalpha;
    real alpha;
    real<lower=0> sigma_a;
    vector[K] beta[L];
    vector[K] zbeta;
    vector<lower=0>[K] sigma_b;
    corr_matrix[K] Rho;
}
model{
    vector[N] p;

    //priors
    alpha ~ normal(0,1.3);
    zalpha ~ normal(0,1);
    sigma_a ~ exponential( 1 );
    Rho ~ lkj_corr( 2 );
    sigma_b ~ exponential( 1 );
    zbeta ~ normal(0,1.3);
    beta ~ multi_normal( zbeta , quad_form_diag(Rho , sigma_b) );
    
    for (i in 1:N) {
       p[i] = zalpha[id[i]] * sigma_a + alpha;
       for (j in 1:K){
          p[i] = p[i] + beta[id[i], j] * bio[i, j];
       }
       p[i] = inv_logit(p[i]);
    }
  
    obs ~ binomial(1, p);
}
"

# Model 1.1 is a simple logistic regression where I do generate log-odds values so I can calculate WAIC directly
model1.1 <- "
data{
    int N;
    int K;
    int L;
    int obs[N];
    real bio[N,K];
    int id[N];
}
parameters{
    vector[L] zalpha;
    real alpha;
    real<lower=0> sigma_a;
    vector[K] beta[L];
    vector[K] zbeta;
    vector<lower=0>[K] sigma_b;
    corr_matrix[K] Rho;
}
model{
    vector[N] p;

    //priors
    alpha ~ normal(0,1.3);
    zalpha ~ normal(0,1);
    sigma_a ~ exponential( 1 );
    Rho ~ lkj_corr( 2 );
    sigma_b ~ exponential( 1 );
    zbeta ~ normal(0,1.3);
    beta ~ multi_normal( zbeta , quad_form_diag(Rho , sigma_b) );
    
    for (i in 1:N) {
       p[i] = zalpha[id[i]] * sigma_a + alpha;
       for (j in 1:K){
          p[i] = p[i] + beta[id[i], j] * bio[i, j];
       }
       p[i] = inv_logit(p[i]);
    }
  
    obs ~ binomial(1, p);
}
generated quantities{
    vector[N] log_lik;
    vector[N] p;
    
    for (i in 1:N) {
       p[i] = zalpha[id[i]] * sigma_a + alpha;
       for (j in 1:K){
          p[i] = p[i] + beta[id[i], j] * bio[i, j];
       }
       p[i] = inv_logit(p[i]);
       log_lik[i] = binomial_lpmf( obs[i] | 1 , p[i] );
    }
}
"

# Model 1.1 is a simple logistic regression where I do generate log-odds values so I can calculate WAIC directly, and I am
# using non-centered priors
model1.2 <- "
data{
    int N;
    int K;
    int L;
    int obs[N];
    real bio[N,K];
    int id[N];
}
parameters{
    vector[L] zalpha;
    real alpha;
    real<lower=0> sigma_a;
    matrix[2,10] zbeta;
    vector<lower=0>[K] sigma_b;
    cholesky_factor_corr[K] L_Rho_b;
}
transformed parameters{
    matrix[L,K] beta;
    beta = (diag_pre_multiply(sigma_b, L_Rho_b) * zbeta)';
}
model{
    vector[N] p;
    sigma_b ~ exponential( 1 );
    L_Rho_b ~ lkj_corr_cholesky( 2 );
    sigma_a ~ exponential( 1 );
    zalpha ~ normal( 0 , 1 );
    alpha ~ normal( 0 , 1.3 );
    to_vector( z_id ) ~ normal( 0 , 1 );
    for ( i in 1:9120 ) {
        p[i] = zalpha[id[i]] * sigma_a + alpha + beta[id[i], 1] * bio1[i] + beta[id[i], 2] * bio2[i];
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
generated quantities{
    matrix[2,2] Rho_id;
    Rho_id = multiply_lower_tri_self_transpose(L_Rho_id);
}
"


bla <- "
data{
  int obs[9120];
  vector[9120] bio2;
  vector[9120] bio1;
  int id[9120];
}
parameters{
  matrix[2,10] z_id;
  vector[10] alpha;
  real alpha_bar;
  real<lower=0> sigma_alpha;
  cholesky_factor_corr[2] L_Rho_id;
  vector<lower=0>[2] sigma_id;
}
transformed parameters{
  matrix[10,2] beta;
  beta = (diag_pre_multiply(sigma_id, L_Rho_id) * z_id)';
}
model{
    vector[9120] p;
    sigma_id ~ exponential( 1 );
    L_Rho_id ~ lkj_corr_cholesky( 2 );
    sigma_alpha ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1 );
    alpha ~ normal( 0 , 1 );
    to_vector( z_id ) ~ normal( 0 , 1 );
    for ( i in 1:9120 ) {
        p[i] = alpha[id[i]] * sigma_alpha + alpha_bar + beta[id[i], 1] * bio1[i] + beta[id[i], 2] * bio2[i];
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
generated quantities{
    vector[9120] log_lik;
    vector[9120] p;
    matrix[2,2] Rho_id;
    Rho_id = multiply_lower_tri_self_transpose(L_Rho_id);
    for ( i in 1:9120 ) {
        p[i] = alpha[id[i]] * sigma_alpha + alpha_bar + beta[id[i], 1] * bio1[i] + beta[id[i], 2] * bio2[i];
        p[i] = inv_logit(p[i]);
    }
    for ( i in 1:9120 ) log_lik[i] = binomial_lpmf( obs[i] | 1 , p[i] );
}

"