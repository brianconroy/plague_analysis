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
  
  int<lower=1> N;
  int<lower=0> N_obs;
  int<lower=0> I_obs[N_obs];
  int<lower=0> P;
  matrix[N, P] X;
  matrix[N, N] D;
  int<lower=0> Y_pos[N_obs];
  int<lower=0> Y_neg[N_obs];
  int<lower=0> Y_loc[N];
  
}
transformed data {
  
  real delta = 1e-9;

}
parameters {
  
  // spatial covariance parameters
  real<lower=0> theta;
  real<lower=0> phi;
  real<lower=0> theta_pos;
  real<lower=0> phi_pos;
  real<lower=0> theta_neg;
  real<lower=0> phi_neg;
  
  // white noise 
  vector[N] xi;
  vector[N_obs] xi_pos;
  vector[N_obs] xi_neg;
  
  // covariate slopes and intercepts 
  vector[P] beta_pos;
  vector[P] beta_neg;
  vector[P] beta_loc;
  
  // preferential sampling parameters
  real alpha_pos;
  real alpha_neg;

}
transformed parameters {
  
  vector[N] w;
  vector[N_obs] eta_pos;
  vector[N_obs] eta_neg;

  {
    // shared latent process
    matrix[N, N] L_K;
    matrix[N, N] K = spatial_cov(D, phi, theta, delta);
    L_K = cholesky_decompose(K);
    w = L_K * xi;
  }
  
  {
    // residual spatial process for disease positive observations
    matrix[N_obs, N_obs] L_K_pos;
    matrix[N_obs, N_obs] K_pos = spatial_cov(D[I_obs, I_obs], phi_pos, theta_pos, delta);
    L_K_pos = cholesky_decompose(K_pos);
    eta_pos = L_K_pos * xi_pos;
  }
  
  {
    // residual spatial process for disease negative observations
    matrix[N_obs, N_obs] L_K_neg;
    matrix[N_obs, N_obs] K_neg = spatial_cov(D[I_obs, I_obs], phi_neg, theta_neg, delta);
    L_K_neg = cholesky_decompose(K_neg);
    eta_neg = L_K_neg * xi_neg; 
  }

}
model {
  
  // spatial covariance priors
  theta ~ gamma(100, 2);
  phi ~ std_normal();
  theta_pos ~ gamma(100, 2);
  phi_pos ~ std_normal();
  theta_neg ~ gamma(100, 2);
  phi_neg ~ std_normal();
  
  // white noise distributions
  xi ~ std_normal();
  xi_pos ~ std_normal();
  xi_neg ~ std_normal();
  
  // preferential sampling priors
  alpha_pos ~ std_normal();
  alpha_neg ~ std_normal();

  // covariate slope and intercept priors
  beta_pos ~ std_normal();
  beta_neg ~ std_normal();
  beta_loc ~ std_normal();
  
  Y_pos ~ poisson_log(X[I_obs] * beta_pos + eta_pos + alpha_pos * w[I_obs]);
  
  Y_neg ~ poisson_log(X[I_obs] * beta_neg + eta_neg + alpha_neg * w[I_obs]);
  
  Y_loc ~ poisson_log(X * beta_loc + w);
  
} 
