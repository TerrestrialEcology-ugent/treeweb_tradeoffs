functions {
	vector isolation(matrix distance, vector area, int id_need, real alpha){
		//int n = rows(distance);
		vector[id_need] isol_tmp = rep_vector(0, id_need);
		vector[id_need] isol_out = rep_vector(0, id_need);
		//int id[n] = 1:n;
		for(i in 1:id_need){
			isol_tmp[i] = sum(area * exp(-alpha * distance[i,])) - area[i]; //minus area[i] to remove area effect of focal patch, based on the fact that exp(0) = 1
		}
		isol_out = (isol_tmp - mean(isol_tmp)) / sd(isol_tmp); //scale
		return isol_out;
	}
}

data{
	int<lower=1> N; //number of observations
	int<lower=1> K; //number of predictor variables
	int<lower=1> id_need; //the number of fragments where data were collected from
	int<lower=1> N_frag; //the total number of fragment in the area, may be different from id_need 
	matrix[N,K] X; //the model matrix including intercept
	real y[N]; //the response variable
	int<lower=1,upper=id_need> frag_id[N]; //index of the fragment in which each plots are situated
	matrix[N_frag,N_frag] distt; //matrix of distance between all fragments
	vector[N_frag] area; //a vector with the area of all fragment
}

parameters{
	vector[K] beta; //the regression parameters
        //real gamma; //isolation effect
	real<lower = 0> sigma; 
	real<lower = 0> alpha;
}

model{
	vector[N] mu;
	vector[id_need] isol  = isolation(distt, area, id_need, alpha);

	beta[1] ~ normal(0,10);
	beta[2:K] ~ normal(0,5);
	alpha ~ normal(1, 5);
	sigma ~ normal(0, 5);

	

	for(n in 1:N){
		mu[n] = X[n] * beta + isol[frag_id[n]];
	}

	y ~ normal_log(mu, sigma); //the likelihood
}

generated quantities{
	real y_gen[N];
	vector[N] mu_gen;
	vector[id_need] isol = isolation(distt, area, id_need, alpha);
	for(n in 1:N){
		mu_gen[n] = X[n] * beta + isol[frag_id[n]];
		y_gen[n] = normal_rng(mu_gen[n], sigma);
	}
}

