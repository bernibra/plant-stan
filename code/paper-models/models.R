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
    real minp;
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
        Y[i] ~ bernoulli(exp(-alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i])) + minp);
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
    real<lower=0> sigma_a;
    real<lower=0> etasq_a;
    real<lower=0> rhosq_a;
    vector<lower=0>[K] sigma_b;
    vector<lower=0>[K] etasq_b;
    matrix<lower=0>[K,2] rhosq_b;
    vector<lower=0>[K] sigma_g;
    vector<lower=0>[K] etasq_g;
    matrix<lower=0>[K,2] rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    matrix[K,L] beta;
    matrix[K,L] gamma;
    matrix[L, L] L_SIGMA_a;
    matrix[L, L] L_SIGMA_b[K];
    matrix[L, L] L_SIGMA_g[K];

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

    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_a ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_a ~ exponential( 0.5 );
    to_vector(rhosq_b) ~ exponential( 0.5 );
    to_vector(rhosq_g) ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    to_vector( zgamma ) ~ std_normal();
    to_vector( zbeta ) ~ std_normal();

    for ( i in 1:L ){
        Y[i] ~ bernoulli(exp(-alpha[i] - gamma[1,i] * columns_dot_self(X1 - beta[1, i]) - gamma[2,i] * columns_dot_self(X2 - beta[2, i])));
    }
}
"

# Categorical 1d
categorical.model.1d <- "
functions{
    matrix cov_GPL2(matrix x, real a, real b, real delta) {
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
    int M;
    int Y[N,L];
    vector[N] X1;
    real minp;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    ordered[M] zphi;
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_b;
    real<lower=0> etasq_b;
    real<lower=0> rhosq_b;
    real<lower=0> sigma_g;
    real<lower=0> etasq_g;
    real<lower=0> rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;
    vector<lower=0>[M] phi;

    alpha = exp(zalpha * sigma_a + alpha_bar);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
    for(i in 1:M){
        phi[i] = exp(-exp(zphi[i]));
    }
}
model{
    real p;
    vector[M+1] prob;

    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    zphi ~ std_normal();
    
    for ( i in 1:L ){
        for (j in 1:N){
           p = exp(- alpha[i] - gamma[i] * pow(X1[j] - beta[i],2));

           prob[1] = 1 - phi[1]*p - minp*(M+1);
           
           for (k in 2:M){
               prob[k]  =  phi[k-1]*p - phi[k]*p + minp;
           }
           
           prob[M+1]  = phi[M]*p + minp;
           //print(prob);

           Y[j,i] ~ categorical(prob);
        }
    }
}
generated quantities{
     vector[L*N] log_lik;
     int l;
     real p;
     vector[M+1] prob;
     
     l = 1;
     
     for ( i in 1:L ){
        for (j in 1:N){
           p = exp(- alpha[i] - gamma[i] * pow(X1[j] - beta[i],2));

           prob[1] = 1 - phi[1]*p - minp*(M+1);
           
           for (k in 2:M){
               prob[k]  =  phi[k-1]*p - phi[k]*p + minp;
           }
           
           prob[M+1]  = phi[M]*p + minp;
           //print(prob);

           log_lik[l] = categorical_lpmf(Y[j,i] | prob);
           l = l + 1;
        }
    }
}
"

# Categorical 1d
categorical.model.skew.1d <- "
functions{
    matrix cov_GPL2(matrix x, real a, real b, real delta) {
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
    int sgn(real x) {
        return x < 0 ? -1 : x > 0;
    }
    real findmax(real delta, real beta, real gamma, real lambda) {
        real maxy;
        real muz;
        maxy = 0.5 * ( 4.0 - pi() ) * pow((delta * sqrt( 2.0 / pi() )),3) / pow(1.0 - 2.0 * pow(delta,2) / pi() ,1.5);
        muz = sqrt( 2.0 / pi() ) * delta;
        maxy = beta + (1.0 / sqrt( 2.0 * gamma)) * (muz - maxy * sqrt( 1.0 - muz * muz ) * 0.5 - 0.5 * sgn(lambda) * exp(- 2.0 * pi() / fabs(lambda) ));
        maxy = exp(- gamma * pow(maxy - beta,2)) * (1.0 + erf((lambda * ( maxy - beta )) * sqrt( gamma ) ));
        return maxy+0.0001;
    }
}
data{
    int N;
    int L;
    int M;
    int Y[N,L];
    vector[N] X1;
    real minp;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    ordered[M] zphi;
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector<lower=-10, upper=10>[L] zlambda;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real lambda_bar;
    real<lower=0> sigma_l;
    real<lower=0> sigma_a;
    real<lower=0> sigma_b;
    real<lower=0> etasq_b;
    real<lower=0> rhosq_b;
    real<lower=0> sigma_g;
    real<lower=0> etasq_g;
    real<lower=0> rhosq_g;
}
transformed parameters{
    vector[L] lambda;
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;
    vector<lower=0>[M] phi;
    real delta;

    lambda = zlambda * sigma_l + lambda_bar;

    alpha = exp(zalpha * sigma_a + alpha_bar);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
    for(i in 1:M){
        phi[i] = exp(-exp(zphi[i]));
    }
    
    for (i in 1:L){
        delta = lambda[i] / ( sqrt( 1.0 + (pow(lambda[i],2)) ));
        gamma[i] = gamma[i] * (1.0 - (2 * pow(delta,2))/pi());
        beta[i] = beta[i] - sqrt( 1.0 / (2 * gamma[i]) ) * (delta * sqrt(2/pi()));
        alpha[i] = log(findmax(delta, beta[i], gamma[i], lambda[i])) + alpha[i];
    }

}
model{
    real p;
    vector[M+1] prob;

    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    sigma_l ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    lambda_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    zlambda ~ std_normal();
    zphi ~ std_normal();
    
    for ( i in 1:L ){
        for (j in 1:N){
           p = exp(- alpha[i] - gamma[i] * pow(X1[j] - beta[i],2)) * (1 + erf((lambda[i] * (X1[j] - beta[i])) * sqrt(gamma[i]) ));

           prob[1] = 1 - phi[1]*p - minp*(M+1);
           
           for (k in 2:M){
               prob[k]  =  phi[k-1]*p - phi[k]*p + minp;
           }
           
           prob[M+1]  = phi[M]*p + minp;
           //print(prob);

           Y[j,i] ~ categorical(prob);
        }
    }
}
generated quantities{
     vector[L*N] log_lik;
     int l;
     real p;
     vector[M+1] prob;
     
     l = 1;
     
     for ( i in 1:L ){
        for (j in 1:N){
           p = exp(- alpha[i] - gamma[i] * pow(X1[j] - beta[i],2)) * (1 + erf((lambda[i] * (X1[j] - beta[i])) * sqrt(gamma[i]) ));

           prob[1] = 1 - phi[1]*p - minp*(M+1);
           
           for (k in 2:M){
               prob[k]  =  phi[k-1]*p - phi[k]*p + minp;
           }
           
           prob[M+1]  = phi[M]*p + minp;
           //print(prob);

           log_lik[l] = categorical_lpmf(Y[j,i] | prob);
           l = l + 1;
        }
    }
}
"
