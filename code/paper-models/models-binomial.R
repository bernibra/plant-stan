# Binomial 1d
base.model.1d <- "
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
    real minp;
    // int K;
    // int Y[K];
    int Y[L, N];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
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

    alpha = exp(zalpha * sigma_a + alpha_bar);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    

}
model{
    // matrix[L,N] p;

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

    for ( i in 1:L ){
        // p[i] = exp(-alpha[i] - gamma[i] * columns_dot_self(X1 - beta[i]));
        Y[i] ~ binomial(1, exp(-alpha[i] - gamma[i] * columns_dot_self(X1 - beta[i])) + minp);
    }
    // Y ~ binomial(1, to_vector(p));
}
generated quantities{
    vector[L*N] log_lik;
    int k;
    
    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           log_lik[k] = binomial_lpmf(Y[i, j] | 1, exp(-alpha[i] - gamma[i] * pow(X1[j] - beta[i],2)) + minp);
           k = k + 1;
        }
    }
}
"

# Binomial 1d
base.model.generror.1d <- "
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
    real minp;
    // int K;
    // int Y[K];
    int Y[L, N];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector[L] znu;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real nu_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_n;
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
    vector[L] nu;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;

    alpha = exp(zalpha * sigma_a + alpha_bar);
    nu = fabs(znu * sigma_n + nu_bar)+1;

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
    for (i in 1:L){
       gamma[i] = sqrt( (gamma[i] * tgamma(3/nu[i])) / tgamma(1/nu[i]));
    }
}
model{
    // matrix[L,N] p;

    sigma_a ~ exponential( 1 );
    sigma_n ~ exponential( 2 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    nu_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    znu ~ std_normal();

    for ( i in 1:L ){
        for( j in 1:N ){
           Y[i, j] ~ binomial(1, exp(-alpha[i] - pow(gamma[i] * fabs(X1[j] - beta[i]), nu[i])) + minp);
        }
    }
}
generated quantities{
    vector[L*N] log_lik;
    int k;
    
    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           log_lik[k] = binomial_lpmf(Y[i, j] | 1, exp(-alpha[i] - pow(gamma[i] * fabs(X1[j] - beta[i]), nu[i])) + minp);
           k = k + 1;
        }
    }
}
"

base.model.generror.1d.multithread <- "
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
    real partial_sum(int[] y_slice,
                   int start, int end,
                   int N, vector alpha, vector beta, vector gamma, vector nu, int[ , ] Y, row_vector X1, real minp ) {
        real lp = 0.0;
        for( i in start:end){
            for( j in 1:N ){
                lp += binomial_lpmf( Y[i, j] | 1 , exp(-alpha[i] - pow(gamma[i] * fabs(X1[j] - beta[i]), nu[i])) + minp);
            }
        }
        return lp;
    }
}
data{
    int N;
    int L;
    real minp;
    // int K;
    // int Y[K];
    int Y[L, N];
    int indices[L];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector[L] znu;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real nu_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_n;
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
    vector[L] nu;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;

    alpha = exp(zalpha * sigma_a + alpha_bar);
    nu = fabs(znu * sigma_n + nu_bar)+1;

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
    for (i in 1:L){
       gamma[i] = sqrt( (gamma[i] * tgamma(3/nu[i])) / tgamma(1/nu[i]));
    }
}
model{
    // matrix[L,N] p;

    sigma_a ~ exponential( 1 );
    sigma_n ~ exponential( 2 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    nu_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    znu ~ std_normal();
    
    int grainsize = 1;

    target += reduce_sum(partial_sum, indices,
                     grainsize,
                     N, alpha, beta, gamma, nu, Y, X1, minp);

}
generated quantities{
    vector[L*N] log_lik;
    int k;
    
    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           log_lik[k] = binomial_lpmf(Y[i, j] | 1, exp(-alpha[i] - pow(gamma[i] * fabs(X1[j] - beta[i]), nu[i])) + minp);
           k = k + 1;
        }
    }
}
"

# Binomial 1d
base.model.traits.1d <- "
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
    real minp;
    int Y[K];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
    matrix[L,L] Dmat_t;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real<lower=0> sigma_a;
    real<lower=0> etasq_a;
    real<lower=0> rhosq_a;
    real<lower=0> sigma_b;
    real<lower=0> etasq_b;
    vector<lower=0>[2] rhosq_b;
    real<lower=0> sigma_g;
    real<lower=0> etasq_g;
    vector<lower=0>[2] rhosq_g;
}
transformed parameters{
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    matrix[L, L] L_SIGMA_a;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;

    L_SIGMA_a = cholesky_decompose(cov_GPL2_alpha( Dmat_t, etasq_a, rhosq_a, sigma_a));
    alpha = L_SIGMA_a * zalpha + alpha_bar;
    alpha = exp(alpha);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, Dmat_t, etasq_b, rhosq_b[1], rhosq_b[2], sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, Dmat_t, etasq_g, rhosq_g[1], rhosq_g[2], sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    

}
model{
    matrix[L,N] p;

    sigma_a ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_a ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_a ~ exponential( 0.5 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();

    for ( i in 1:L ){
        p[i] = exp(-alpha[i] - gamma[i] * columns_dot_self(X1 - beta[i])) + minp;
    }
    
    Y ~ binomial(1, to_vector(p'));
    
}
"


# Skew 1d
skew.model.traits.1d <- "
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
    int Y[L,N];
    real minp;
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector<lower=-10, upper=10>[L] lambda;
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
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    //vector[L] lambda;
    real delta;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;
    
    //lambda = zlambda * sigma_l + lambda_bar;
    
    alpha = exp(zalpha * sigma_a + alpha_bar);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);

    for (i in 1:L){
        delta = lambda[i] / ( sqrt( 1.0 + (pow(lambda[i],2)) ));
        gamma[i] = gamma[i] * (1.0 - (2 * pow(delta,2))/pi());
        beta[i] = beta[i] - sqrt( 1.0 / (2 * gamma[i]) ) * (delta * sqrt(2/pi()));
        alpha[i] = log(findmax(delta, beta[i], gamma[i], lambda[i])) + alpha[i];
    }
    
}
model{
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
    lambda ~ normal( lambda_bar , sigma_l );

    for ( i in 1:L ){
        Y[i] ~ binomial(1, exp(- alpha[i] - gamma[i] * columns_dot_self(X1 - beta[i])) .* (1.0 + erf((lambda[i] * (X1 - beta[i])) * sqrt(gamma[i]) )) + minp);
    }
}
generated quantities{
     vector[L*N] log_lik;
     int k;
     
     k = 1;
     for ( i in 1:L ){
         for (j in 1:N){
            log_lik[k] = binomial_lpmf(Y[i, j] | 1,  exp(- alpha[i] - gamma[i] * pow(X1[j] - beta[i],2)) * (1 + erf((lambda[i] * (X1[j] - beta[i])) * sqrt(gamma[i]) )) + minp);
            k = k + 1;
         }
     }
}
"

base.model.skew.1d.multithread <- "
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
    real partial_sum(int[] y_slice,
                   int start, int end,
                   int N, vector lambda, vector alpha, vector beta, vector gamma, int[ , ] Y, row_vector X1, real minp ) {
        real lp = 0.0;
        real b;
        real g;
        real a;
        real delta;
        
        for( i in start:end){
        
            delta = lambda[i] / ( sqrt( 1.0 + (pow(lambda[i],2)) ));
            g = gamma[i] * (1.0 - (2 * pow(delta,2))/pi());
            b = beta[i] - sqrt( 1.0 / (2 * g) ) * (delta * sqrt(2/pi()));
            a = log(findmax(delta, b, g, lambda[i])) + alpha[i];
            
            lp += binomial_lpmf( Y[i] | 1 , exp(- a - g * columns_dot_self(X1 - b)) .* (1.0 + erf((lambda[i] * (X1 - b)) * sqrt(g) )) + minp);
            
        }
        return lp;
    }
}
data{
    int N;
    int L;
    int Y[L,N];
    int indices[L];
    real minp;
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector<lower=-10, upper=10>[L] lambda;
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
    vector[L] alpha;
    vector[L] beta;
    vector[L] gamma;
    //vector[L] lambda;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;
    
    //lambda = zlambda * sigma_l + lambda_bar;
    
    alpha = exp(zalpha * sigma_a + alpha_bar);

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
}
model{
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
    lambda ~ normal( lambda_bar , sigma_l );

    int grainsize = 1;

    target += reduce_sum(partial_sum, indices,
                     grainsize,
                     N, lambda, alpha, beta, gamma, Y, X1, minp);

}
generated quantities{
     vector[L*N] log_lik;
     vector[L] alpha_;
     vector[L] beta_;
     vector[L] gamma_;
     int k;
     real delta;
     
     k = 1;
     for ( i  in 1:L ){
        delta = lambda[i] / ( sqrt( 1.0 + (pow(lambda[i],2)) ));
        gamma_[i] = gamma[i] * (1.0 - (2 * pow(delta,2))/pi());
        beta_[i] = beta[i] - sqrt( 1.0 / (2 * gamma_[i]) ) * (delta * sqrt(2/pi()));
        alpha_[i] = log(findmax(delta, beta_[i], gamma_[i], lambda[i])) + alpha[i];
     }
     
     for ( i in 1:L ){
         for (j in 1:N){
            log_lik[k] = binomial_lpmf(Y[i, j] | 1,  exp(- alpha_[i] - gamma_[i] * pow(X1[j] - beta_[i],2)) * (1 + erf((lambda[i] * (X1[j] - beta_[i])) * sqrt(gamma_[i]) )) + minp);
            k = k + 1;
         }
     }
}
"

# Binomial 1d
base.model.skew.1d.multithread2 <- "
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
    real partial_sum(int[] y_slice,
                   int start, int end,
                   int N, vector lambda, vector alpha, vector beta, vector gamma, int[ , ] Y, row_vector X1, real minp ) {
        real lp = 0.0;
        real m;
        real v;
        for( i in start:end){
            v = gamma[i] * sqrt((pi()*(1+3*pow(lambda[i],2))-8.0*pow(lambda[i],2))/(2.0*pi()));
            m = beta[i] - 2.0*lambda[i]/(sqrt(pi())*v);
            for( j in 1:N ){
                lp += binomial_lpmf( Y[i, j] | 1 , exp(-alpha[i] - pow(v * fabs(X1[j] - m)/(1+lambda[i]*sgn(X1[j] - m)), 2.0)) + minp);
            }
        }
        return lp;
    }
}
data{
    int N;
    int L;
    real minp;
    // int K;
    // int Y[K];
    int Y[L, N];
    int indices[L];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector[L] zlambda;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real lambda_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_l;
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

    alpha = exp(zalpha * sigma_a + alpha_bar);
    lambda = inv_logit(zlambda*sigma_l + lambda_bar)*2-1;

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
}
model{
    // matrix[L,N] p;

    sigma_a ~ exponential( 1 );
    sigma_l ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
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
    
    int grainsize = 1;

    target += reduce_sum(partial_sum, indices,
                     grainsize,
                     N, lambda, alpha, beta, gamma, Y, X1, minp);
}
generated quantities{
    vector[L*N] log_lik;
    vector[L] betam;
    vector[L] gammav;
    int k;
    
    for (i in 1:L){
        gammav[i] = gamma[i] * sqrt((pi()*(1+3*pow(lambda[i],2))-8.0*pow(lambda[i],2))/(2.0*pi()));
        betam[i] = beta[i] - 2.0*lambda[i]/(sqrt(pi())*gammav[i]);
    }
    
    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           log_lik[k] = binomial_lpmf(Y[i, j] | 1, exp(-alpha[i] - pow(gammav[i] * fabs(X1[j] - betam[i])/(1+lambda[i]*sgn(X1[j] - betam[i])), 2.0)) + minp);
           k = k + 1;
        }
    }
}
"

# Binomial 1d
base.model.skew.generror.1d <- "
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
}
data{
    int N;
    int L;
    real minp;
    // int K;
    // int Y[K];
    int Y[L, N];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector[L] znu;
    vector[L] zlambda;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real nu_bar;
    real lambda_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_n;
    real<lower=0> sigma_l;
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
    vector[L] nu;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;

    alpha = exp(zalpha * sigma_a + alpha_bar);
    nu = exp(znu * sigma_n + nu_bar)+1;
    lambda = inv_logit(zlambda*sigma_l + lambda_bar)*2-1;

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
    for (i in 1:L){
       gamma[i] = gamma[i] * sqrt((pi()*(1+3*pow(lambda[i],2))*tgamma(3.0/nu[i])-pow(16,(1.0/nu[i]))*pow(lambda[i],2)*tgamma(0.5+1.0/nu[i])*tgamma(0.5+1.0/nu[i])*tgamma(1.0/nu[i]))/(pi()*tgamma(1.0/nu[i])));
       beta[i] = beta[i] - pow(2, 2.0/nu[i])*lambda[i]*tgamma(0.5+1.0/nu[i])/(sqrt(pi())*gamma[i]);
    }
}
model{
    // matrix[L,N] p;

    sigma_a ~ exponential( 1 );
    sigma_n ~ exponential( 1 );
    sigma_l ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    lambda_bar ~ std_normal();
    nu_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    znu ~ std_normal();
    zlambda ~ std_normal();

    for ( i in 1:L ){
        for( j in 1:N ){
           Y[i, j] ~ binomial(1, exp(-alpha[i] - pow(gamma[i] * fabs(X1[j] - beta[i])/(1+lambda[i]*sgn(X1[j] - beta[i])), nu[i])) + minp);
        }
    }
}
generated quantities{
    vector[L*N] log_lik;
    int k;
    
    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           log_lik[k] = binomial_lpmf(Y[i, j] | 1, exp(-alpha[i] - pow(gamma[i] * fabs(X1[j] - beta[i])/(1+lambda[i]*sgn(X1[j] - beta[i])), nu[i])) + minp);
           k = k + 1;
        }
    }
}
"

# Binomial 1d
base.model.skew.generror.1d.multithread <- "
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
    real partial_sum(int[] y_slice,
                   int start, int end,
                   int N, vector lambda, vector alpha, vector beta, vector gamma, vector nu, int[ , ] Y, row_vector X1, real minp ) {
        real lp = 0.0;
        real m;
        real v;
        for( i in start:end){
            v = gamma[i] * sqrt((pi()*(1+3*pow(lambda[i],2))*tgamma(3.0/nu[i])-pow(16,(1.0/nu[i]))*pow(lambda[i],2)*tgamma(0.5+1.0/nu[i])*tgamma(0.5+1.0/nu[i])*tgamma(1.0/nu[i]))/(pi()*tgamma(1.0/nu[i])));
            m = beta[i] - pow(2, 2.0/nu[i])*lambda[i]*tgamma(0.5+1.0/nu[i])/(sqrt(pi())*v);
            for( j in 1:N ){
                lp += binomial_lpmf( Y[i, j] | 1 , exp(-alpha[i] - pow(v * fabs(X1[j] - m)/(1+lambda[i]*sgn(X1[j] - m)), nu[i])) + minp);
            }
        }
        return lp;
    }
}
data{
    int N;
    int L;
    real minp;
    // int K;
    // int Y[K];
    int Y[L, N];
    int indices[L];
    row_vector[N] X1;
    matrix[L,L] Dmat_b;
    matrix[L,L] Dmat_g;
}
parameters{
    vector[L] zalpha;
    vector[L] zbeta;
    vector[L] zgamma;
    vector[L] znu;
    vector[L] zlambda;
    real alpha_bar;
    real beta_bar;
    real gamma_bar;
    real nu_bar;
    real lambda_bar;
    real<lower=0> sigma_a;
    real<lower=0> sigma_n;
    real<lower=0> sigma_l;
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
    vector[L] nu;
    matrix[L, L] L_SIGMA_b;
    matrix[L, L] L_SIGMA_g;

    alpha = exp(zalpha * sigma_a + alpha_bar);
    nu = exp(znu * sigma_n + nu_bar)+1;
    lambda = inv_logit(zlambda*sigma_l + lambda_bar)*2-1;

    L_SIGMA_b = cholesky_decompose(cov_GPL2(Dmat_b, etasq_b, rhosq_b, sigma_b));
    beta = L_SIGMA_b * zbeta + beta_bar;

    L_SIGMA_g = cholesky_decompose(cov_GPL2(Dmat_g, etasq_g, rhosq_g, sigma_g));
    gamma = L_SIGMA_g * zgamma + gamma_bar;
    gamma = exp(gamma);
    
}
model{
    // matrix[L,N] p;

    sigma_a ~ exponential( 1 );
    sigma_n ~ exponential( 1 );
    sigma_l ~ exponential( 1 );
    sigma_b ~ exponential( 1 );
    sigma_g ~ exponential( 1 );
    etasq_b ~ exponential( 1 );
    etasq_g ~ exponential( 1 );
    rhosq_b ~ exponential( 0.5 );
    rhosq_g ~ exponential( 0.5 );
    alpha_bar ~ normal( 0 , 1.3 );
    beta_bar ~ std_normal();
    gamma_bar ~ std_normal();
    lambda_bar ~ std_normal();
    nu_bar ~ std_normal();
    zalpha ~ std_normal();
    zgamma ~ std_normal();
    zbeta ~ std_normal();
    znu ~ std_normal();
    zlambda ~ std_normal();
    
    int grainsize = 1;

    target += reduce_sum(partial_sum, indices,
                     grainsize,
                     N, lambda, alpha, beta, gamma, nu, Y, X1, minp);
}
generated quantities{
    vector[L*N] log_lik;
    vector[L] betam;
    vector[L] gammav;
    int k;
    
    for (i in 1:L){
        gammav[i] = gamma[i] * sqrt((pi()*(1+3*pow(lambda[i],2))*tgamma(3.0/nu[i])-pow(16,(1.0/nu[i]))*pow(lambda[i],2)*tgamma(0.5+1.0/nu[i])*tgamma(0.5+1.0/nu[i])*tgamma(1.0/nu[i]))/(pi()*tgamma(1.0/nu[i])));
        betam[i] = beta[i] - pow(2, 2.0/nu[i])*lambda[i]*tgamma(0.5+1.0/nu[i])/(sqrt(pi())*gammav[i]);
    }
    
    k = 1;
    for ( i in 1:L ){
        for (j in 1:N){
           log_lik[k] = binomial_lpmf(Y[i, j] | 1, exp(-alpha[i] - pow(gammav[i] * fabs(X1[j] - betam[i])/(1+lambda[i]*sgn(X1[j] - betam[i])), nu[i])) + minp);
           k = k + 1;
        }
    }
}
"

