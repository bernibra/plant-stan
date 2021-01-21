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

# Here I model the parameters as Gaussian Processes but keep the linear form
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

# Parameters modelled as Gaussian Processes and using a radial basis function in the linear model and variances are not gaussian processes
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
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real<lower=0> sigma_a;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_b;
    vector<lower=0>[K] rhosq_b;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA[K];
    for(i in 1:K){
        L_SIGMA[i] = cholesky_decompose(cov_GPL2(Dmat, etasq_b[i], rhosq_b[i], sigma_b[i]));
        beta[i] = zbeta[i] * (L_SIGMA[i]') + beta_bar[i];
    }
    for(i in 1:K){
        gamma[i] = zgamma[i] * sigma_g[i] + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }
    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    vector[N] p;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ normal( 0 , 1 );
    gamma_bar ~ normal( 0 , 1 );
    zalpha ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );
    to_vector(zgamma) ~ normal( 0 , 1 );

    for ( i in 1:N ) {
        p[i] = alpha[id[i]];
        for (j in 1:K){
          p[i] = p[i] - gamma[j, id[i]] * radial(bio[i, j], beta[j, id[i]]);
       }
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"

# Same as 3.0 but two of the parameters are modelled as Gaussian Processes.
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
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;

}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real<lower=0> sigma_a;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_b;
    vector<lower=0>[K] rhosq_b;
    vector<lower=0>[K] etasq_g;
    vector<lower=0>[K] rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b[i], rhosq_b[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g[i], rhosq_g[i], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    vector[N] p;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    rhosq_g ~ exponential( 0.5 );
    etasq_g ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ normal( 0 , 1 );
    gamma_bar ~ normal( 0 , 1 );
    zalpha ~ normal( 0 , 1 );
    to_vector(zgamma) ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );

    for ( i in 1:N ) {
        p[i] = alpha[id[i]];
        for (j in 1:K){
          p[i] = p[i] - gamma[j, id[i]] * radial(bio[i, j], beta[j, id[i]]);
       }
        p[i] = inv_logit(p[i]);
    }
    obs ~ binomial( 1 , p );
}
"

# This is the same as model 3.1 but appears to be a bit faster. The structure of the model is very different but it helps me vectorize the operations
model3.2 <- "
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
    matrix[K,N] bio;
    int id[N];
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;

}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real<lower=0> sigma_a;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_b;
    vector<lower=0>[K] rhosq_b;
    vector<lower=0>[K] etasq_g;
    vector<lower=0>[K] rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b[i], rhosq_b[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g[i], rhosq_g[i], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    vector[N] p;
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    rhosq_g ~ exponential( 0.5 );
    etasq_g ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ normal( 0 , 1 );
    gamma_bar ~ normal( 0 , 1 );
    zalpha ~ normal( 0 , 1 );
    to_vector(zgamma) ~ normal( 0 , 1 );
    to_vector( zbeta ) ~ normal( 0 , 1 );

    for ( i in 1:N ) {
        p[i] = alpha[id[i]];
        for (j in 1:K){
          p[i] = p[i] - gamma[j, id[i]] * radial(bio[j,i], beta[j, id[i]]);
       }
        obs[i] ~ binomial( 1 , inv_logit(p[i]) );
    }
}
"

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
model4.1 <- "
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
    int L;
    int K;
    int Y[L,N];
    row_vector[N] X1;
    row_vector[N] X2;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;

}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real<lower=0> sigma_a;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_b;
    vector<lower=0>[K] rhosq_b;
    vector<lower=0>[K] etasq_g;
    vector<lower=0>[K] rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b[i], rhosq_b[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g[i], rhosq_g[i], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    rhosq_g ~ exponential( 0.5 );
    etasq_g ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    to_vector(zgamma) ~ std_normal();
    to_vector( zbeta ) ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ bernoulli_logit(alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i]));
    }
}
"

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
model5.1 <- "
functions{
    matrix cov_GPL2(matrix x, matrix y, real sq_alpha, real sq_alpha_t, real sq_rho, real sq_rho_t, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + sq_alpha_t + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) ) + sq_alpha_t * exp(-sq_rho_t * square(y[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + sq_alpha_t + delta;
        return K;
    }
}
data{
    int N;
    int L;
    int K;
    int Y[L,N];
    row_vector[N] X1;
    row_vector[N] X2;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
    matrix[L,L] Dmat_t;
}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real<lower=0> sigma_a;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_b;
    vector<lower=0>[K] rhosq_b;
    vector<lower=0>[K] etasq_g;
    vector<lower=0>[K] rhosq_g;
    vector<lower=0>[K] etasq_tb;
    vector<lower=0>[K] rhosq_tb;
    vector<lower=0>[K] etasq_tg;
    vector<lower=0>[K] rhosq_tg;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_t, etasq_b[i], etasq_tb[i], rhosq_b[i], rhosq_tb[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_g, Dmat_t, etasq_g[i], etasq_tg[i], rhosq_g[i], rhosq_tg[i], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    rhosq_g ~ exponential( 0.5 );
    etasq_g ~ exponential( 1 );
    rhosq_tb ~ exponential( 0.5 );
    etasq_tb ~ exponential( 1 );
    rhosq_tg ~ exponential( 0.5 );
    etasq_tg ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    to_vector(zgamma) ~ std_normal();
    to_vector( zbeta ) ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ bernoulli_logit(alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i]));
    }
}
"

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
model5.2 <- "
functions{
    matrix cov_GPL2(matrix x, matrix y, matrix z, real a1, real a2, real a3, real b1, real b2, real b3, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;

        for (i in 1:(N-1)) {
          K[i, i] = a1 + a2 + a3 + delta;
          for (j in (i + 1):N) {
            K[i, j] = a1 * exp(- b1 * square(x[i,j]) ) + a2 * exp(-b2 * square(y[i,j]) )  + a3 * exp(-b3 * square(z[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = a1 + a2 + a3 + delta;
        return K;
    }
}
data{
    int N;
    int L;
    int K;
    int Y[L,N];
    row_vector[N] X1;
    row_vector[N] X2;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
    matrix[L,L] Dmat_t;
}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real zsigma_a;
    vector[3] zetasq_a;
    vector[3] zrhosq_a;
    vector[K] zsigma_b;
    matrix[K,3] zetasq_b;
    matrix[K,3] zrhosq_b;
    vector[K] zsigma_g;
    matrix[K,3] zetasq_g;
    matrix[K,3] zrhosq_g;
    real<lower=0> sigma_sigma;
    real sigma_mean;
    real<lower=0> rho_sigma;
    real rho_mean;
    real<lower=0> eta_sigma;
    real eta_mean;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_a;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    real sigma_a;
    vector[3] etasq_a;
    vector[3] rhosq_a;
    vector[K] sigma_b;
    matrix[K,3] etasq_b;
    matrix[K,3] rhosq_b;
    vector[K] sigma_g;
    matrix[K,3] etasq_g;
    matrix[K,3] rhosq_g;
    
    sigma_a = exp(zsigma_a * sigma_sigma + sigma_mean);
    etasq_a = exp(zetasq_a * eta_sigma + eta_mean);    
    rhosq_a = exp(zrhosq_a * rho_sigma + rho_mean);
    sigma_b = exp(zsigma_b * sigma_sigma + sigma_mean);
    etasq_b = exp(zetasq_b * eta_sigma + eta_mean);    
    rhosq_b = exp(zrhosq_b * rho_sigma + rho_mean);
    sigma_g = exp(zsigma_g * sigma_sigma + sigma_mean);
    etasq_g = exp(zetasq_g * eta_sigma + eta_mean);    
    rhosq_g = exp(zrhosq_g * rho_sigma + rho_mean);    

    L_SIGMA_a = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_g, Dmat_t, etasq_a[1], etasq_a[2], etasq_a[3], rhosq_a[1], rhosq_a[2], rhosq_a[3], sigma_a));
    alpha = L_SIGMA_a * zalpha + alpha_bar;

    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_g, Dmat_t, etasq_b[i, 1], etasq_b[i, 2], etasq_b[i, 3], rhosq_b[i, 1], rhosq_b[i, 2], rhosq_b[i, 3], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_g, Dmat_t, etasq_g[i, 1], etasq_g[i, 2], etasq_g[i, 3], rhosq_g[i, 1], rhosq_g[i, 2], rhosq_g[i, 3], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

}
model{
    zsigma_a ~ std_normal();
    zsigma_b ~ std_normal();
    zsigma_g ~ std_normal();
    zrhosq_a ~ std_normal();
    zrhosq_a ~ std_normal();
    to_vector(zrhosq_b) ~ std_normal();
    to_vector(zrhosq_g) ~ std_normal();
    to_vector(zrhosq_b) ~ std_normal();
    to_vector(zrhosq_g) ~ std_normal();
    sigma_sigma ~ exponential( 1 );
    sigma_mean ~ std_normal();
    rho_sigma ~ exponential( 1 );
    rho_mean ~ std_normal();
    eta_sigma ~ exponential( 1 );
    eta_mean ~ std_normal();
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    to_vector( zgamma ) ~ std_normal();
    to_vector( zbeta ) ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ bernoulli_logit(alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i]));
    }
}
"

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
model5.3 <- "
functions{
    matrix cov_GPL2(matrix x, matrix y, matrix z, real a1, real a2, real a3, real b1, real b2, real b3, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;

        for (i in 1:(N-1)) {
          K[i, i] = a1 + a2 + a3 + delta;
          for (j in (i + 1):N) {
            K[i, j] = a1 * exp(- b1 * square(x[i,j]) ) + a2 * exp(-b2 * square(y[i,j]) )  + a3 * exp(-b3 * square(z[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = a1 + a2 + a3 + delta;
        return K;
    }
}
data{
    int N;
    int L;
    int K;
    int Y[L,N];
    row_vector[N] X1;
    row_vector[N] X2;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
    matrix[L,L] Dmat_t;
}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real zsigma_a;
    vector[3] zetasq_a;
    vector[3] zrhosq_a;
    vector[K] zsigma_b;
    matrix[K,3] zetasq_b;
    matrix[K,3] zrhosq_b;
    vector[K] zsigma_g;
    matrix[K,3] zetasq_g;
    matrix[K,3] zrhosq_g;
    real<lower=0> sigma_sigma;
    real sigma_mean;
    real<lower=0> rho_sigma;
    real rho_mean;
    real<lower=0> eta_sigma;
    real eta_mean;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_a;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    real sigma_a;
    vector[3] etasq_a;
    vector[3] rhosq_a;
    vector[K] sigma_b;
    matrix[K,3] etasq_b;
    matrix[K,3] rhosq_b;
    vector[K] sigma_g;
    matrix[K,3] etasq_g;
    matrix[K,3] rhosq_g;
    
    sigma_a = exp(zsigma_a * sigma_sigma + sigma_mean);
    etasq_a = exp(zetasq_a * eta_sigma + eta_mean);    
    rhosq_a = exp(zrhosq_a * rho_sigma + rho_mean);
    sigma_b = exp(zsigma_b * sigma_sigma + sigma_mean);
    etasq_b = exp(zetasq_b * eta_sigma + eta_mean);    
    rhosq_b = exp(zrhosq_b * rho_sigma + rho_mean);
    sigma_g = exp(zsigma_g * sigma_sigma + sigma_mean);
    etasq_g = exp(zetasq_g * eta_sigma + eta_mean);    
    rhosq_g = exp(zrhosq_g * rho_sigma + rho_mean);    

    L_SIGMA_a = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_g, Dmat_t, etasq_a[1], etasq_a[2], etasq_a[3], rhosq_a[1], rhosq_a[2], rhosq_a[3], sigma_a));
    alpha = exp(L_SIGMA_a * zalpha + alpha_bar);

    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_g, Dmat_t, etasq_b[i, 1], etasq_b[i, 2], etasq_b[i, 3], rhosq_b[i, 1], rhosq_b[i, 2], rhosq_b[i, 3], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_g, Dmat_t, etasq_g[i, 1], etasq_g[i, 2], etasq_g[i, 3], rhosq_g[i, 1], rhosq_g[i, 2], rhosq_g[i, 3], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

}
model{
    zsigma_a ~ std_normal();
    zsigma_b ~ std_normal();
    zsigma_g ~ std_normal();
    zrhosq_a ~ std_normal();
    zrhosq_a ~ std_normal();
    to_vector(zrhosq_b) ~ std_normal();
    to_vector(zrhosq_g) ~ std_normal();
    to_vector(zrhosq_b) ~ std_normal();
    to_vector(zrhosq_g) ~ std_normal();
    sigma_sigma ~ exponential( 1 );
    sigma_mean ~ std_normal();
    rho_sigma ~ exponential( 1 );
    rho_mean ~ std_normal();
    eta_sigma ~ exponential( 1 );
    eta_mean ~ std_normal();
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    to_vector( zgamma ) ~ std_normal();
    to_vector( zbeta ) ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ binomial(1, exp(- alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i]))) ;
    }
}
"


########################################
########################################
# The following are just trials...
########################################
########################################

# This is the same as model 5.1, but only one variable and I am adding correlations between the hyperparameters of the gaussian process. 
model6.1 <- "
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
    int L;
    int K;
    int Y[L,N];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
    matrix[L,L] Zer;
}
parameters{
    vector[L] zalpha;
    vector[2*L] zbeta;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_b;
    real<lower=0> sigma_g;
    real<lower=0> etasq_b;
    real<lower=0> rhosq_b;
    real<lower=0> etasq_g;
    real<lower=0> rhosq_g;
    real<lower=-1, upper=1> rho;
}
transformed parameters{
    vector[L] alpha;
    vector[2*L] beta;
    matrix[2*L, 2*L] L_SIGMA;
    
    L_SIGMA[1:L,1:L] = cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b);
    L_SIGMA[(L+1):(2*L), (L+1):(2*L)] = cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g);
    L_SIGMA[(L+1):(2*L),1:L] = Zer + rho * sqrt(etasq_b + sigma_b) * sqrt(etasq_g + sigma_g);
    L_SIGMA[1:L,(L+1):(2*L)] = Zer + rho * sqrt(etasq_b + sigma_b) * sqrt(etasq_g + sigma_g);

    L_SIGMA = cholesky_decompose(L_SIGMA);
    
    beta = L_SIGMA*zbeta;
    
    beta[1:L] = beta[1:L] + beta_bar;
    beta[(L+1):(2*L)] = exp(beta[(L+1):(2*L)] + gamma_bar);

    alpha = zalpha * sigma_a + alpha_bar;
}
model{
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    rhosq_g ~ exponential( 0.5 );
    etasq_g ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    zbeta ~ std_normal();
    rho ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ bernoulli_logit(alpha[i] - beta[(i + L)] * columns_dot_self(X1 - beta[i]));
    }
}
"

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
model7.1 <- "
functions{
    matrix cov_GPL2(matrix x, matrix y, real sq_alpha, real sq_alpha_t, real sq_rho, real sq_rho_t, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + sq_alpha_t + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * square(x[i,j]) ) + sq_alpha_t * exp(-sq_rho_t * square(y[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + sq_alpha_t + delta;
        return K;
    }
}
data{
    int N;
    int L;
    int K;
    int Y[L,N];
    row_vector[N] X1;
    row_vector[N] X2;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
    matrix[L,L] Dmat_t;
}
parameters{
    vector[L] zalpha;
    matrix[K,L] zbeta;
    matrix[K,L] zgamma;
    real alpha_bar;
    vector[K] beta_bar;
    vector[K] gamma_bar;
    real<lower=0> sigma_a;
    
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_b;
    vector<lower=0>[K] rhosq_b;
    vector<lower=0>[K] etasq_g;
    vector<lower=0>[K] rhosq_g;
    vector<lower=0>[K] etasq_tb;
    vector<lower=0>[K] rhosq_tb;
    vector<lower=0>[K] etasq_tg;
    vector<lower=0>[K] rhosq_tg;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];
    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_t, etasq_b[i], etasq_tb[i], rhosq_b[i], rhosq_tb[i], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_g, Dmat_t, etasq_g[i], etasq_tg[i], rhosq_g[i], rhosq_tg[i], sigma_g[i]));
        gamma[i] = zgamma[i]*(L_SIGMA_g[i]') + gamma_bar[i];
        gamma[i] = exp(gamma[i]);
    }

    alpha = zalpha * sigma_a + alpha_bar;
}
model{

    Rho ~ lkj_corr( 2 );
    sigma_cafe ~ exponential( 1 );
    b ~ normal( -1 , 0.5 );
    a ~ normal( 5 , 2 );
    {
    vector[2] YY[K];
    vector[2] MU;
    MU = [ a , b ]';
    for ( j in 1:K ) YY[j] = [ a_cafe[j] , b_cafe[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho , sigma_cafe) );
    }
    
    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    etasq_b ~ exponential( 1 );
    rhosq_g ~ exponential( 0.5 );
    etasq_g ~ exponential( 1 );
    rhosq_tb ~ exponential( 0.5 );
    etasq_tb ~ exponential( 1 );
    rhosq_tg ~ exponential( 0.5 );
    etasq_tg ~ exponential( 1 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    to_vector(zgamma) ~ std_normal();
    to_vector( zbeta ) ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ bernoulli_logit(alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i]));
    }
}
"
