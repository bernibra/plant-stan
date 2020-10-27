source("./prepare-data.R")
# source("./models.R")
library(rethinking)
library(rstan)

####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression.
####

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan <- function(d = NULL, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=F, show.plot=T, ofolder="../../results/models/"){
        
        # Load the data
        if(is.null(d)){
                if(recompile){
                       d <- species_distribution.data(variables=variables)
                       saveRDS(d, file = paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = ""))
                }else{
                       d <- readRDS(file = paste("../../data/processed/jsdm/", paste(variables, collapse = ""), "data.rds", sep = ""))
                }
        }

        # Prepare training data for stan model
        obs <- d$obs
        id <- d$id
        bio <- d[,6:ncol(d)]
        
        # Estimate posterior distributions using rstan
        # Define variables of the model
        dat_1.1 <- list(N=length(obs),
                        K=length(variables),
                        L=length(unique(id)),
                        obs=obs,
                        bio=bio,
                        id=id)
        
        # Set starting values for the parameters
        start_1.1 <- list(
                alpha = rep(0, dat_1.1$L),
                zalpha = 0,
                sigma_a = 0.1,
                beta = matrix(0,dat_1.1$L,dat_1.1$K),
                zbeta = rep(0, dat_1.1$K),
                sigma_b = rep(0.1, dat_1.1$K)
        )
        
        # Initialize data structure
        n_chains_1.1 <- 3
        init_1.1 <- list()
        for ( i in 1:n_chains_1.1 ) init_1.1[[i]] <- start_1.1
        
        # Should I calculate the likelihood?
        if(loglik){
                model_code=model1.1
        }else{
                model_code=model1.0     
        }
        
        # Run stan model
        mfit_1.1 <- stan ( model_code=model_code ,
                           data=dat_1.1 ,
                           chains=n_chains_1.1 ,
                           cores= n_chains_1.1 ,
                           warmup=1000, iter=4000,
                           init=init_1.1 , control = list(adapt_delta = 0.95))
        
        # 
        # dat <- list(obs=obs, id=id, bio1=bio$bio12_p_8110, bio2=bio$bio5_tmaxw_8110)
        # 
        # # Use rethinking package to fit stan model
        # m1 <- ulam(
        #         alist(
        #                 obs ~ dbinom( 1 , p ),
        #                 logit(p) <- alpha[id, 1] + alpha[id, 2] * bio1 + alpha[id, 3] * bio2,
        #                 transpars> matrix[id,3]:alpha <- compose_noncentered( sigma_id , L_Rho_id , z_id ),
        #                 matrix[3,id]:z_id ~ dnorm( 0 , 1),
        #                 cholesky_factor_corr[3]:L_Rho_id ~ lkj_corr_cholesky( 2 ),
        #                 vector[3]:sigma_id ~ dexp(1),
        #                 gq> matrix[3,3]:Rho_id <<- multiply_lower_tri_self_transpose(L_Rho_id)
        #         ) , data=dat , chains=3,  iter = 3000,  log_lik = T, cores = 3)
        # 
        # # Use rethinking package to fit stan model
        # m2 <- ulam(
        #         alist(
        #                 obs ~ dbinom( 1 , p ),
        #                 logit(p) <- a[id] + b[id] * bio1 + c[id] * bio2,
        #                 c(a,b,c)[id] ~ multi_normal(c(za,zb,zc),Rho_id,sigma_a),
        #                 za ~ normal(0,1.3),
        #                 zb ~ normal(0,1.3),
        #                 zc ~ normal(0,1.3),
        #                 sigma_a ~ dexp(1),
        #                 Rho_id ~ lkj_corr(2)
        #         ) , data=dat , chains=3,  iter = 3000,  log_lik = T, cores = 3)
        # 
        # saveRDS(m2, file = paste(ofolder, "model1.rds", sep=""))
        return(mfit_1.1)
}

binomial.stan(recompile = T, ofolder="/cluster/scratch/plant-stan/")

