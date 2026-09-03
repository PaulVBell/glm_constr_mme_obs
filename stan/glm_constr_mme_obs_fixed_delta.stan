data {
    int<lower=0> len_a;                             // Length of vector of regression parameters (including intercept).
    int<lower=0> m;                                 // Number of simulations.
    array[m] vector[len_a] b_hat;                   // MLE for regression parameters for each simulation.
    array[m] cholesky_factor_cov[len_a] J;          // Lower Cholesky decomposition of inverse of Fisher information matrix for each simulation.
    vector[len_a] a_hat;                            // MLE for regression parameters for each simulation.
    cholesky_factor_cov[len_a] J_a;                 // Lower Cholesky decomposition of inverse of Fisher information matrix for each simulation.
    int<lower=0> len_hyper_param;                   // 6 or 0 (if no priors)
    vector[len_hyper_param] hyper_param;            // prior hyper params: mu[1,2]_mu, mu[1,2]_sd, eta, tau_gamma (scale).
    vector[len_a] delta;                            // prior hyper params on dissim param: intercept alpha and beta; and slope alpha and beta.
}

parameters {
    array[m] vector[len_a] b_tilde;                 // Regression parameters for each simulation scaled to be standard normal random variables.
    vector[len_a] b0_tilde;                         // Regression parameters for the real world assuming exchangeability scaled to be standard normal random variables.
    vector[len_a] w_tilde;                          // Bias between real world and MME scaled to be standard normal random variables.
    vector[len_a] a_ind;                            // Regression parameters for the real world independent of the MME.
    vector[len_a] mu;                               // Mean of the latent layer.
    cholesky_factor_corr[len_a] L;                  // Cholesky decomposition of the correlation matrix of latent layer.
    vector<lower=0>[len_a] tau;                     // Vector of scaling variance terms for latent layer.
}

transformed parameters {
    // Latent layer (using non central parameterisation)
    array[m] vector[len_a] b;
    for(i in 1:m){
        b[i] = mu + diag_pre_multiply(tau, L) * b_tilde[i];
    }
    vector[len_a] b0 = mu + diag_pre_multiply(tau, L) * b0_tilde;
    vector[len_a] w = tau .* ((1 - delta)^.5) .* w_tilde;
    vector[len_a] a = b0 + w;
}

model {
    // Priors
    if(len_hyper_param > 0){
        target+= normal_lpdf(mu | hyper_param[1:len_a], hyper_param[(len_a+1):(2*len_a)]);
        target+= lkj_corr_cholesky_lpdf(L | hyper_param[2*len_a + 1]);
        target+= cauchy_lpdf(tau |0, hyper_param[2*len_a + 2]);
    }

    for(i in 1:m){
        target+= multi_normal_cholesky_lpdf(b_hat[i] | b[i], J[i]); // "Data" layer
        target+= std_normal_lpdf(b_tilde[i]);                       // Latent layer
    }
    target+= std_normal_lpdf(b0_tilde);
    target+= - 0.5 * delta[1] * w_tilde[1]^2 - log(sqrt(2 * pi())) + 0.5 * log(delta[1]);
    target+= - 0.5 * delta[2] * w_tilde[2]^2 - log(sqrt(2 * pi())) + 0.5 * log(delta[2]);


    // Data layer
    target+= multi_normal_cholesky_lpdf(a_hat | a, J_a);
    target+= multi_normal_cholesky_lpdf(a_hat | a_ind, J_a);
}

generated quantities {
    cov_matrix[len_a] Sigma;
    corr_matrix[len_a] Omega;
    Omega = multiply_lower_tri_self_transpose(L);
    Sigma = quad_form_diag(Omega, tau);
    vector[len_a] b_mme = multi_normal_rng(mu, Sigma);
}


