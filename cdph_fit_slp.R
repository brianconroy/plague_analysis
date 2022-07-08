library(rstan)
library(paws)

# ----------------------- Constants ------------------------- #

OUT_FILE <- 'plague_slp_fit_agg6.rds'

set.seed(123)

Sys.setenv(
  AWS_ACCESS_KEY_ID = AWS_ACCESS_KEY_ID,
  AWS_SECRET_ACCESS_KEY = AWS_SECRET_ACCESS_KEY,
  AWS_REGION = AWS_REGION
)

svc <- paws::s3()

# ----------------------- Analysis ------------------------- #

# retrieve data
svc$download_file(
  Bucket = "plague-analysis",
  Key = 'plague_data.rds',
  Filename = 'plague_data.rds'
)
dat <- readRDS(file = 'plague_data.rds')

print(paste("fraction observed:", round(dat$N_obs/dat$N, 2)))
print(paste("Y (+) total:", sum(dat$Y_pos)))
print(paste("Y (-) total:", sum(dat$Y_neg)))
print(paste("N:", sum(dat$N)))

# Fit STAN model
fit <- stan(
  file = 'slp.stan', 
  data = dat,
  chains = 4,
  cores = 4)

saveRDS(fit, file = OUT_FILE)

svc$put_object(
  Bucket = "plague-analysis",
  Body = OUT_FILE,
  Key = OUT_FILE
)
