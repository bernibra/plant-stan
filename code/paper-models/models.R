# GAM with binary data

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
base.model <- "
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

    alpha = exp(zalpha * sigma_a + alpha_bar);
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
        Y[i] ~ bernoulli(exp(-alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i])));
    }
}
"

# This is the same as model 3.1 and 3.2 but much faster. 
# The structure of the model is very different but it helps me vectorize the operations and reduce the size of the objects used
base.model.traits <- "
functions{
    matrix cov_GPL2(matrix x, matrix y, real a, real b1, real b2, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;

        for (i in 1:(N-1)) {
          K[i, i] = a + delta;
          for (j in (i + 1):N) {
            K[i, j] = a * exp(- b1 * square(x[i,j]) -b2 * square(y[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = a + delta;
        return K;
    }
    matrix cov_GPL2_alpha(matrix x, real a, real b, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;

        for (i in 1:(N-1)) {
          K[i, i] = a + delta;
          for (j in (i + 1):N) {
            K[i, j] = a * exp(- b * square(x[i,j]) );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = a + delta;
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
    real zetasq_a;
    real zrhosq_a;
    vector[K] zsigma_b;
    vector[K] zetasq_b;
    matrix[K,2] zrhosq_b;
    vector[K] zsigma_g;
    vector[K] zetasq_g;
    matrix[K,2] zrhosq_g;
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
    real etasq_a;
    vector[3] rhosq_a;
    vector[K] sigma_b;
    vector[K] etasq_b;
    matrix[K,2] rhosq_b;
    vector[K] sigma_g;
    vector[K] etasq_g;
    matrix[K,2] rhosq_g;
    
    sigma_a = exp(zsigma_a * sigma_sigma + sigma_mean);
    etasq_a = exp(zetasq_a * eta_sigma + eta_mean);    
    rhosq_a = exp(zrhosq_a * rho_sigma + rho_mean);
    sigma_b = exp(zsigma_b * sigma_sigma + sigma_mean);
    etasq_b = exp(zetasq_b * eta_sigma + eta_mean);    
    rhosq_b = exp(zrhosq_b * rho_sigma + rho_mean);
    sigma_g = exp(zsigma_g * sigma_sigma + sigma_mean);
    etasq_g = exp(zetasq_g * eta_sigma + eta_mean);    
    rhosq_g = exp(zrhosq_g * rho_sigma + rho_mean);    

    L_SIGMA_a = cholesky_decompose(cov_GPL2_alpha( Dmat_t, etasq_a, rhosq_a, sigma_a));
    alpha = L_SIGMA_a * zalpha + alpha_bar;
    alpha = exp(alpha);

    for(i in 1:K){
        L_SIGMA_b[i] = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_t, etasq_b[i], rhosq_b[i, 1], rhosq_b[i, 2], sigma_b[i]));
        beta[i] = zbeta[i]*(L_SIGMA_b[i]') + beta_bar[i];
    }
    for(i in 1:K){
        L_SIGMA_g[i] = cholesky_decompose(cov_GPL2(Dmat_g, Dmat_t, etasq_g[i], rhosq_g[i, 1], rhosq_g[i, 2], sigma_g[i]));
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
    zetasq_b ~ std_normal();
    zetasq_g ~ std_normal();
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
        Y[i] ~ bernoulli(-exp(alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i])));
    }
}
"

