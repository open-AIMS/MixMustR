// Code built by MixMustR on 2025-12-02 04:42:00.513561
// User options:
	// sample_tracer = TRUE
	// fix_unsampled = TRUE
	// hierarchical = TRUE
data {
	int<lower=1> N; // number of observations
	int<lower=1> J; // number of sampled (a.k.a considered) sources
	int<lower=1> K; // number of tracers
	matrix[N, K] Y; // data stream 1
	matrix[N, J + 1] ln_rho; // Log mixing proportions, data stream 2
	matrix[N, J + 1] sigma_ln_rho; // Confidence around `ln_rho`
	matrix[J, K] x; // sample means for each source and tracer
	matrix<lower=0>[J, K] s; // sample SDs for each source and tracer
	matrix<lower=1>[J, K] m; // source- and tracer- specific sample size
	int<lower=1> R; // number of levels in the hierarchical structure
	array[N] int<lower=1> YR; // dummy vector of random levels
}
parameters {
	array[N] vector[J + 1] zeta; // zeta from softmax
	array[R] vector[J + 1] randvec; // random effects
	array[J + 1] real<lower=0> randsd; // random effects hyper parameter
	array[J + 1] real<lower=0> ext_sigma; // residual error, likelihood 2
	matrix[J, K] mu; // mean for each source and tracer (J x K)
	matrix<lower=0>[J, K] nu; // degrees of freedom
}
transformed parameters {
	array[N] simplex[J + 1] p; // probability vector for each site
	matrix[J, K] omega; // SD for each source and tracer
	for (n in 1:N) {
		p[n] = softmax(zeta[n]);
	}
	for (j in 1:J) {
		for (k in 1:K) {
			omega[j, k] = sqrt((s[j, k]^2 * (m[j, k] - 1)) / nu[j, k]);
		}
	}
}
model {
	// Sampling sources
	for (j in 1:J) {
		for (k in 1:K) {
			target += chi_square_lpdf(m[j, k] | nu[j, k]);
			target += normal_lpdf(x[j, k] | mu[j, k], s[j, k] / sqrt(m[j, k]));
		}
	}
	// Mixture model
	for (n in 1:N) {
		for (k in 1:K) {
			array[J + 1] real phi;
			// Contribution from sampled sources
			for (j in 1:J) {
				phi[j] = log(p[n][j]) + normal_lpdf(Y[n, k] | mu[j, k], sqrt(p[n][j]^2 * omega[j, k]^2));
			}
			// Contribution from unsampled sources
			real mean_mu_k = mean(mu[, k]);
			real mean_omega_k = sqrt(sum(omega[, k]^2) / J);
			phi[J + 1] = log(p[n][J + 1]) + normal_lpdf(Y[n, k] | mean_mu_k, sqrt(p[n][J + 1]^2 * mean_omega_k^2));
			// Log-sum-exp for numerical stability
			target += log_sum_exp(phi);
		}
	}
	// Independent evaluation of mixing proportion (zeta scale)
	for (n in 1:N) {
		ln_rho[n, ] ~ normal(zeta[n], sigma_ln_rho[n, ]);
		zeta[n] ~ normal(randvec[YR[n]], ext_sigma);
	}
	for (j in 1:(J+1)) {
		randvec[][j] ~ normal(0, randsd[j]);
	}
	randsd ~ normal(0, 1);
	ext_sigma ~ normal(0, 1);
}
