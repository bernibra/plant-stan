# JSDM with ordinal data

# Multinomial regression with RBFs in the linear form
model1.1 <- "
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
    int M;
    int Y[L,N];
    row_vector[N] X1;
    row_vector[N] X2;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;

}
parameters{
    //ordered[L] zalpha[M];
    
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
        Y[i] ~ ordered_logistic( gamma[1,i] * columns_dot_self(X1 - beta[1, i]) + gamma[2,i] * columns_dot_self(X2 - beta[2, i]), alpha[i] );
    }
}
"
