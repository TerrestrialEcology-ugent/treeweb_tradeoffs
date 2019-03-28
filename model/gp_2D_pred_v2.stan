// Fit the hyperparameters of a Gaussian process with an 
// exponentiated quadratic kernel
// allow for fixed effect on the mean

data {
  int<lower=1> N; // sample size
  int<lower=1> D; // dimension of the gaussian process
  int<lower=1> P; // number of predictors
  vector[D] coords[N]; // position
  matrix[N, P] X; //predictor matrix
  vector[N] y; // response
}

parameters {
  real<lower=0> rho; // length-scale parameter
  real<lower=0> alpha; // marginal standard deviation
  real<lower=0> sigma_spat; // spatial error
  vector[P] beta; // regression parameters
  vector[N] gamma; // spatial effect
  real<lower=0> sigma_res; // residual error
}
model {
  matrix[N, N] L_K; 

  matrix[N, N] K; // covariance matrix
  K = cov_exp_quad(coords, alpha, rho);

  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + pow(sigma_spat, 2);
  
  L_K = cholesky_decompose(K);
  
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  sigma_spat ~ normal(0, 1);
  sigma_res ~ normal(0, 1);
  beta ~ cauchy(0, 5);

  gamma ~ multi_normal_cholesky(rep_vector(0,N), L_K);

  y ~ normal(X * beta + gamma, sigma_res);
}

