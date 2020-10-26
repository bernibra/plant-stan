source("./sdm-maxent.R")
source("./models.R")
library(rethinking)
library(rstan)
library(gtools)


####
## Run species distribution model for all species. I wrote the following models:
# - Simple binomial regression.
####

# This first one is just a logistic regression, all "variables" are used as predictors in a linear form.
binomial.stan <- function(d = NULL, idx=128, variables=c("bio5_", "bio6_","bio12_"), recompile = T, loglik=F, show.plot=T, ofolder="../../results/models/"){
        
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
        
        dat <- list(obs=obs, id=id, bio1=bio$bio12_p_8110, bio2=bio$bio5_tmaxw_8110, K=length(unique(id)))
        
        # Use rethinking package to fit stan model
        m1 <- ulam(
                alist(
                        obs ~ dbinom( 1 , p ),
                        logit(p) <- alpha[id, 1] + alpha[id, 2] * bio1 + alpha[id, 3] * bio2,
                        transpars> matrix[K,3]:alpha <- compose_noncentered( rep_vector(sigma_alpha,3) , L_Rho_alpha , z_alpha ),
                        matrix[3,K]:z_alpha ~ dnorm( 0 , 1),
                        cholesky_factor_corr[3]:L_Rho_alpha ~ lkj_corr_cholesky( 2 ),
                        sigma_alpha ~ exponential(1),
                        gq> matrix[3,3]:Rho_alpha <<- Chol_to_Corr( L_Rho_alpha )
                ) , data=dat , chains=3,  iter = 3000,  log_lik = T)
        
        saveRDS(m1, file = paste(ofolder, "model1.rds", sep=""))
        return(m1)
}
