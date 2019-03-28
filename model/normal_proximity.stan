data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables
	matrix[N,K] X; //the model matrix including intercept
	real y[N]; //the response variable
	int<lower=1,upper=id_need> frag_id[N]; //index of the fragment in which each plots are situated
	matrix[N_frag,N_frag] distt; //matrix of distance between all fragments
	vector[N_frag] area; //a vector with the area of all fragment
}

parameters{
	vector[K] beta; //the regression parameters
        real gamma; //isolation effect
	real<lower = 0> sigma; 
	real <lower = 0> alpha;
}

model{
	vector[id_need] isol = isolation(distt, area, id_need, alpha);
	vector[N] mu;

	beta[1] ~ normal(0,10);
	beta[2:K] ~ normal(0,5);
	alpha ~ normal(1, 5);
	sigma ~ normal(0, 5);

	

	for(n in 1:N){
		mu[n] = X[n] * beta + gamma * isol[frag_id[n]];
	}

	y ~ normal_log(mu, sigma); //the likelihood
}

generated quantities{
	real y_gen[N];
	vector[N] mu_gen;
	for(n in 1:N){
		mu_gen[n] = X[n] * beta + alpha;
		y_gen[n] = normal_rng(mu_gen[n], sigma);
	}
}

