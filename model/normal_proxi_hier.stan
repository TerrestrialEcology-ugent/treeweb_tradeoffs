
functions { 
} 
data { 
  int<lower=1> n;  // total number of observations 
  int<lower=1> n_mix; // number of tree compositoin categories
  int<lower=1> k; // number of covariates
  int<lower=1> f; //total number of functions

  vector[n] y;  // response variable 
  int<lower=1, upper=n_mix> mix[n]; // indicator variable of the tree mixture
  int<lower=1, upper=f> f_id; // indicator of the function
  vector[n] proxi; // proximity index
  vector[n] edge; // edge effect
  matrix[n, k] X; // assume that the first column is the intercept, then the categorical indeces and then the conitnuous variable(s)
  
} 
transformed data { 
} 
parameters { 
  vector[k] beta;
  real<lower=0> sigma;  // residual SD 
} 
transformed parameters { 
  vector[n] y_hat; // linear predictor
  for(i in 1:n)
    y_hat[i] = X[i,] * beta;
} 
model { 
  // priors including all constants 
  beta ~ normal(0, 5);
  sigma ~ cauchy(0, 2);
  // likelihood including all constants 
  y ~ normal(y_hat, sigma);
} 
generated quantities { 
  vector[n_mix] Beta_mix;
  real<lower=0> s_mix;
  real<lower=0> s_proxi;
  real<lower=0> s_edge;
  real<lower=0> s_y;

  for(i in 2:n_mix)
     Beta_mix[(i-1)] = beta[1] + beta[i];

  Beta_mix[n_mix] = beta[1] - sum(beta[2:n_mix]);

  s_mix = sd(Beta_mix[mix]);
  s_proxi = sd(proxi) * fabs(beta[(n_mix + 1)]);
  s_edge = sd(edge) * fabs(beta[(n_mix + 2)]);
  s_y = sd(y - y_hat);
} 
