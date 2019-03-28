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
  real<lower=0> sigma; // residual error
  vector[P] beta; // regression parameters
  
}
transformed parameters {
  matrix[N, N] K; // covariance matrix
  K = cov_exp_quad(coords, alpha, rho);

  // diagonal elements
  for (n in 1:N)
    K[n, n] = K[n, n] + pow(sigma, 2);

}
model {
  matrix[N, N] L_K; 
  
  L_K = cholesky_decompose(K);
  
  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);
  beta ~ cauchy(0, 5);

  y ~ multi_normal_cholesky(X * beta, L_K);
}

