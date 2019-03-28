/** Simple CAR model
* standard no-efficient parametrization
*/

data {
  int<lower = 1> n; //number of locations
  int<lower = 1> p; // number of predictor
  matrix[n, p] X; //predictor matrix
  real y[n]; //the response
  matrix<lower = 0, upper = 1>[n, n] W; // the adjacency matrix
}
transformed data{
  vector[n] zeros;
  matrix<lower = 0>[n, n] D;
  {
    vector[n] W_rowsums;
    for (i in 1:n) {
      W_rowsums[i] = sum(W[i, ]);
    }
    D = diag_matrix(W_rowsums);
  }
  zeros = rep_vector(0, n);
}
parameters {
  vector[p] beta; // regression coefficient
  vector[n] phi; // spatial effect
  real<lower = 0> tau; // spatial precision
  real<lower = 0, upper = 1> alpha; // spatial dependence
  real<lower = 0> sigma;
}
model {
  sigma ~ normal(5,2);
  phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
  beta ~ normal(0, 1);
  tau ~ gamma(2, 2);
  y ~ normal(X * beta + phi, sigma);
}
