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
    matrix[K,L] zbeta;
    vector[K] beta_bar;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] etasq;
    vector<lower=0>[K] rhosq;
}
transformed parameters{
    matrix[K,L] beta;
    vector[L] alpha;
    matrix[L, L] L_SIGMA[K];
    for(i in 1:K){
        L_SIGMA[i] = cholesky_decompose(cov_GPL2(Dmat, etasq[i], rhosq[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA[i]') + beta_bar[i];
    }
    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    vector[N] p;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 1 );
    beta_bar ~ normal( 0 , 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    zalpha ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        p[i] = alpha[id[i]];
        for (j in 1:K){
          p[i] = p[i] + beta[j, id[i]] * bio[i, j];
       }
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"

model3.0 <- "
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
    real radial(real x, real mu) {
        real K;
        K = x - mu;
        K = pow(K, 2);
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
    matrix[K,L] zbeta;
    matrix<lower=0>[K,L] sigma_beta;
    vector[K] beta_bar;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] etasq;
    vector<lower=0>[K] rhosq;
}
transformed parameters{
    matrix[K,L] beta;
    vector[L] alpha;
    matrix[L, L] L_SIGMA[K];
    for(i in 1:K){
        L_SIGMA[i] = cholesky_decompose(cov_GPL2(Dmat, etasq[i], rhosq[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA[i]') + beta_bar[i];
    }
    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    vector[N] p;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 1 );
    beta_bar ~ normal( 0 , 1 );
    to_vector(sigma_beta) ~ normal( 0 , 0.2 );
    alpha_bar ~ normal( 0 , 1.3 );
    zalpha ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        p[i] = alpha[id[i]];
        for (j in 1:K){
          p[i] = p[i] - sigma_beta[j, id[i]] * radial(bio[i, j], beta[j, id[i]]);
       }
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"


model3.1 <- "
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
    real radial(real x, real mu) {
        real K;
        K = x - mu;
        K = pow(K, 2);
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
    matrix[K,L] zbeta;
    matrix[K,L] sigma_beta;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] etasq;
    vector<lower=0>[K] rhosq;
}
transformed parameters{
    matrix[K,L] beta;
    vector[L] alpha;
    matrix[L, L] L_SIGMA[K];
    for(i in 1:K){
        L_SIGMA[i] = cholesky_decompose(cov_GPL2(Dmat, etasq[i], rhosq[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA[i]');
    }
    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    vector[N] p;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    rhosq ~ exponential( 0.5 );
    etasq ~ exponential( 1 );
    beta_bar ~ normal( 0 , 1 );
    to_vector(sigma_beta) ~ lognormal( 0 , 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    zalpha ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        p[i] = alpha[id[i]];
        for (j in 1:K){
          p[i] = p[i] - sigma_beta[j, id[i]] * radial(bio[i, j], beta[j, id[i]]);
       }
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"
