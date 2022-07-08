functions{
  
    matrix spatial_cov(matrix x, real phi, real theta, real delta) {
      
        int n = dims(x)[1];
        matrix[n, n] K;
        for (i in 1:(n - 1)) {
          K[i, i] = pow(phi, 2) + delta;
          for (j in (i + 1):n) {
            K[i, j] = pow(phi, 2) * exp(- x[i,j] / theta);
            K[j, i] = K[i, j];
          }
        }
        K[n, n] = pow(phi, 2) + delta;
        return K;
        
    }
    
}
data {
  
  int<lower=0> N_obs;
  int<lower=0> P;
  matrix[N_obs, P] X;
  matrix[N_obs, N_obs] D;
  int<lower=0> Y_pos[N_obs];
  int<lower=0> Y_neg[N_obs];

}
transformed data {
  
  real delta = 1e-9;

}
parameters {
  
  // spatial covariance parameters
  real<lower=0> theta_pos;
  real<lower=0> phi_pos;
  real<lower=0> theta_neg;
  real<lower=0> phi_neg;
  
  // white noise 
  vector[N_obs] xi_pos;
  vector[N_obs] xi_neg;
  
  // covariate slopes and intercepts 
  vector[P] beta_pos;
  vector[P] beta_neg;

}
transformed parameters {
  
  vector[N_obs] eta_pos;
  vector[N_obs] eta_neg;
  
  {
    // residual spatial process for disease positive observations
    matrix[N_obs, N_obs] L_K_pos;
    matrix[N_obs, N_obs] K_pos = spatial_cov(D, phi_pos, theta_pos, delta);
    L_K_pos = cholesky_decompose(K_pos);
    eta_pos = L_K_pos * xi_pos;
  }
  
  {
    // residual spatial process for disease negative observations
    matrix[N_obs, N_obs] L_K_neg;
    matrix[N_obs, N_obs] K_neg = spatial_cov(D, phi_neg, theta_neg, delta);
    L_K_neg = cholesky_decompose(K_neg);
    eta_neg = L_K_neg * xi_neg; 
  }

}
model {
  
  // spatial covariance priors
  theta_pos ~ gamma(100, 2);
  phi_pos ~ std_normal();
  theta_neg ~ gamma(100, 2);
  phi_neg ~ std_normal();
  
  // white noise distributions
  xi_pos ~ std_normal();
  xi_neg ~ std_normal();

  // covariate slope and intercept priors
  beta_pos ~ std_normal();
  beta_neg ~ std_normal();

  Y_pos ~ poisson_log(X * beta_pos + eta_pos);
  
  Y_neg ~ poisson_log(X * beta_neg + eta_neg);

}
