// Fit the hyperparameters of a Gaussian process with an 
// exponentiated quadratic kernel

data {
  int<lower=1> N; // sample size
  int<lower=1> D; // dimension of the gaussian process
  vector[D] x[N]; // position
  vector[N] y; // response
}
transformed data {
  vector[N] mu = rep_vector(0, N); // mean response
}
parameters {
  real<lower=0> rho; // length-scale parameter
  real<lower=0> alpha; // marginal standard deviation
  real<lower=0> sigma; // residual error
  
}
transformed parameters {
  matrix[N, N] K; // covariance matrix
  K = cov_exp_quad(x, alpha, rho);

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

  y ~ multi_normal_cholesky(mu, L_K);
}

