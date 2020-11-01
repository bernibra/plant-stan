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
    real alpha_bar;
    real<lower=0> sigma_a;
    matrix[K,L] zbeta;
    vector[K] beta_bar;
    vector<lower=0>[K] sigma_b;
    cholesky_factor_corr[K] L_Rho_b;
}
transformed parameters{
    matrix[L,K] beta;
    vector[L] alpha;
    beta = (diag_pre_multiply(sigma_b, L_Rho_b) * zbeta)';
    for (i in 1:L){
        alpha[i] = zalpha[i] * sigma_a + alpha_bar;
        for (j in 1:K){
          beta[i, j] = beta[i, j] + beta_bar[j];
        }
    }
}
model{
    vector[N] p;
    sigma_b ~ exponential( 1 );
    sigma_a ~ exponential( 1 );
    L_Rho_b ~ lkj_corr_cholesky( 2 );
    zalpha ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ normal( 0 , 1.3 );
    
    for (i in 1:N) {
       p[i] = alpha[id[i]];
       for (j in 1:K){
          p[i] = p[i] + beta[id[i], j] * bio[i, j];
       }
       p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
generated quantities{
    vector[N] log_lik;
    vector[N] p;
    matrix[2,2] Rho_b;
    Rho_b = multiply_lower_tri_self_transpose(L_Rho_b);
    
    for (i in 1:N) {
       p[i] = alpha[id[i]];
       for (j in 1:K){
          p[i] = p[i] + beta[id[i], j] * bio[i, j];
       }
       p[i] = inv_logit(p[i]);
       log_lik[i] = binomial_lpmf( obs[i] | 1 , p[i] );
    }
}
"

model2.0 <- "
functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data{
    int N;
    int K;
    int L;
    int obs[N];
    real bio[N,K];
    int id[N];
    matrix[L,L] Dmat;
}
parameters{
    vector[L] zalpha;
    real alpha_bar;
    real<lower=0> sigma_a;
    vector[L] beta1;
    vector[L] beta2;
    vector[K] beta_bar;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] etasq;
    vector<lower=0>[K] rhosq;
}
model{
    vector[N] p;
    matrix[L,L] SIGMA1;
    matrix[L,L] SIGMA2;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 1 );
    beta_bar ~ normal( 0 , 1.3 );
    alpha_bar ~ normal( 0 , 1.3 );
    zalpha ~ normal( 0 , 1 );
    SIGMA1 = cov_GPL2(Dmat, etasq[1], rhosq[1], sigma_b[1]);
    SIGMA2 = cov_GPL2(Dmat, etasq[2], rhosq[2], sigma_b[2]);
    beta1 ~ multi_normal( rep_vector(0,L) , SIGMA1 );
    beta2 ~ multi_normal( rep_vector(0,L) , SIGMA2 );
    for ( i in 1:N ) {
        p[i] = zalpha[id[i]] * sigma_a + alpha_bar + (beta1[id[i]] + beta_bar[1]) * bio[i,1] + (beta2[id[i]] + beta_bar[2]) * bio[i,2];
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"


"functions{
    matrix cov_GPL2(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data{
    int L;
    int obs[8000];
    vector[8000] bio2;
    vector[8000] bio1;
    int id[8000];
    matrix[40,40] Dmat;
}
parameters{
    vector[L] zbeta1;
    vector[L] zbeta2;
    vector[40] alpha;
    real alpha_bar;
    real beta_bar2;
    real beta_bar1;
    real<lower=0> etasq2;
    real<lower=0> etasq1;
    real<lower=0> rhosq2;
    real<lower=0> rhosq1;
    real<lower=0> si2;
    real<lower=0> si1;
    real<lower=0> sigma_a;
}
transformed parameters{
    vector[L] beta1;
    vector[L] beta2;
    matrix[L,L] L_SIGMA1;
    matrix[L,L] SIGMA1;
    matrix[L,L] L_SIGMA2;
    matrix[L,L] SIGMA2;
    SIGMA2 = cov_GPL2(Dmat, etasq2, rhosq2, si2);
    L_SIGMA2 = cholesky_decompose(SIGMA2);
    SIGMA1 = cov_GPL2(Dmat, etasq1, rhosq1, si1);
    L_SIGMA1 = cholesky_decompose(SIGMA1);
    beta2 = L_SIGMA2 * zbeta2 + beta_bar2;
    beta1 = L_SIGMA1 * zbeta1 + beta_bar1;
}
model{
    vector[8000] p;
    sigma_a ~ exponential( 1 );
    si1 ~ exponential( 1 );
    si2 ~ exponential( 1 );
    rhosq1 ~ exponential( 0.5 );
    rhosq2 ~ exponential( 0.5 );
    etasq1 ~ exponential( 1 );
    etasq2 ~ exponential( 1 );
    beta_bar1 ~ normal( 0 , 1 );
    beta_bar2 ~ normal( 0 , 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    alpha ~ normal( 0 , 1 );
    zbeta2 ~ normal( 0 , 1 );
    zbeta1 ~ normal( 0 , 1 );
    for ( i in 1:8000 ) {
        p[i] = alpha[id[i]] * sigma_a + alpha_bar + beta1[id[i]] * bio1[i] + beta2[id[i]] * bio2[i];
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"
