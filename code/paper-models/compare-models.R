library(rethinking)
library(rstan)

m1 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min20-1d.rds", sep=""))
m2 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min20-1d-generror.rds", sep=""))
m3 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min20-1d-skew.rds", sep=""))
m4 <- readRDS(paste("/cluster/scratch/bemora/plant-stan/min20-1d-skew-generror.rds", sep=""))
comp <- rethinking::compare(m1, m2, m3, m4)

saveRDS(comp, file = "/cluster/scratch/bemora/plant-stan/comparision.rds")
