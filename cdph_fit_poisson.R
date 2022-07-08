library(rstan)

# ----------------------- Constants ------------------------- #

set.seed(123)

# ----------------------- Analysis ------------------------- #

dat <- readRDS(file = 'data/plague_data.rds')

# Subset data for Spatial Poisson model
dat_pois <- list(
  N_obs = dat$N_obs,
  P = 3,
  X = dat$X[dat$I_obs,],
  Y_pos = dat$Y_pos,
  Y_neg = dat$Y_neg
)

print(paste("fraction observed:", round(dat$N_obs/dat$N, 2)))
print(paste("Y (+) total:", sum(dat$Y_pos)))
print(paste("Y (-) total:", sum(dat$Y_neg)))

# Fit STAN model
fit <- stan(
  file = 'poisson.stan', 
  data = dat_pois,
  chains = 4,
  cores = 4)

out_file <- 'cdph_fits/plague_poisson_fit.rds'

saveRDS(fit, file = out_file)
