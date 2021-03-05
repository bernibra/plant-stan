library(rethinking)
library(rstan)

m1 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min30-skew-model-traits-1d2.rds", sep=""))
m2 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min30-baseline-model-1d2.rds", sep=""))
m3 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min30-1d-generror.rds", sep=""))
comp <- rethinking::compare(m1, m2, m3)

saveRDS(comp, file = "/cluster/scratch/bemora/plant-stan/comparision.rds")