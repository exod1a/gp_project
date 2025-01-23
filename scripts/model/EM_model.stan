
// implement EM model for drift and diffusion

functions {
  // linear kernel
  matrix lin_cov(real[] x, real[] xs, real var_l, real var_b, real c) {
    
    int R = size(x);
    int C = size(xs);
    
    matrix[R, C] covariance;
    
    for(i in 1:R) {
      for(j in 1:C)
        covariance[i, j] = var_b + var_l * (x[i] - c)*(xs[j] - c);
    }
    
    return covariance;
  }
  
  real generalized_inverse_gaussian_lpdf(real x, int p, real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
  }
  
  // Custom function for prediction generation
  vector gp_pred_rng(real[] x_real, real[] x_pred, int N_real, int N_pred, vector f_r,
                     real stat_var, real length_scale, real sigma, real sigma_b, 
                     real c, int drift_prior) {
    
    vector[N_pred] f_pred;
    {
      matrix[N_real, N_real] K_r;
      matrix[N_pred, N_pred] K_pred;
      matrix[N_real, N_pred] K_r_pred;      
      matrix[N_real, N_real] L_r;
      vector[N_real]         alpha;
      matrix[N_real, N_pred] nu;
      vector[N_pred]         mu;
      matrix[N_pred, N_pred] cov;
      matrix[N_pred, N_pred] nug_pred;

      // kernels for process
      K_r        = cov_exp_quad(x_real, sqrt(stat_var), length_scale);            // known data kernel
      // add nugget to diagonal elements to ensure positive definiteness
      for (i in 1:N_real)
        K_r[i, i] = K_r[i, i] + 1e-9;
      K_r_pred   = cov_exp_quad(x_real, x_pred, sqrt(stat_var), length_scale);    // known + unknown data kernel
      K_pred     = cov_exp_quad(x_pred, sqrt(stat_var), length_scale);            // unknown data kernel
      nug_pred   = diag_matrix(rep_vector(1e-9, N_pred));                         // nugget for positive definiteness

      // only add linear kernel if it is for the drift function
      if (drift_prior) {
        // add linear kernel for drift inference
        K_r      += lin_cov(x_real, x_real, sigma, sigma_b, c); 
        K_r_pred += lin_cov(x_real, x_pred, sigma, sigma_b, c); 
        K_pred   += lin_cov(x_pred, x_pred, sigma, sigma_b, c); 
      }

      // GP transformed mean
      L_r      = cholesky_decompose(K_r);                             // L
      alpha    = mdivide_left_tri_low(L_r, f_r);                      // L^-1 * f
      alpha    = mdivide_right_tri_low(alpha', L_r)';                 // α = L^-T * L^-1 * f
      mu       = K_r_pred' * alpha;                                   // μ = K_r_pred * α
      // GP transformed covariance 
      nu       = mdivide_left_tri_low(L_r, K_r_pred);                 // ν = L^-1 * K_r_pred
      cov      = K_pred - nu' * nu;                                   // K_pred - ν^T * ν
      f_pred   = multi_normal_rng(mu, cov + nug_pred);                // sample from new distribution
    }
    return f_pred;
  }
}

data {
  int<lower=1> N_real;                      // Number of observations
  int<lower=1> N_pred;                      // Number of points where we don't have data

  real x_real[N_real];                      // Data: x-values
  real x_pred[N_pred];                      // x-values where we use the GP for predictive inference
  
  real dx[N_real];                          // Data: change in x
  real dt[N_real];                          // Data: change in t
  
  real c;                                   // linear kernel center
  
  int<lower=0, upper=1> pred_inf;            // do predictive inference or not
}

parameters {
  // GP hyperparameters
  real<lower=0> length_scale_f;             // Length scale of drift
  real<lower=0> stat_var_f;                 // stationary variance of drift 
  
  real<lower=0> length_scale_g;             // Length scale of diffusion
  real<lower=0> stat_var_g;                 // stationary variance of diffusion
  
  real<lower=0> sigma_b;                    // variance of linear kernel
  real<lower=0> sigma;                      // stat var of linear kernel
  
  // latent GP variables
  vector[N_real] eta_f;
  vector[N_real] eta_g;
}

transformed parameters {
  vector[N_real] f;                       // drift GP prior
  vector[N_real] g;                       // diffusion GP prior
  {
    matrix[N_real, N_real] L_f;          // Cholesky decompositions
    matrix[N_real, N_real] L_g;
    
    // drift and diffusion kernels
    matrix[N_real, N_real] K_f = cov_exp_quad(x_real, sqrt(stat_var_f), length_scale_f);
    matrix[N_real, N_real] K_g = cov_exp_quad(x_real, sqrt(stat_var_g), length_scale_g);
    // add linear kernel
    K_f += lin_cov(x_real, x_real, sigma, sigma_b, c);

    // add nugget to diagonal elements to ensure positive definiteness
    for (i in 1:N_real) {
      K_f[i, i] = K_f[i, i] + 1e-9;
      K_g[i, i] = K_g[i, i] + 1e-9;
    }

    L_f = cholesky_decompose(K_f);
    L_g = cholesky_decompose(K_g);
    
    f = L_f * eta_f;
    g = L_g * eta_g;
  }
}

model {
  // Priors
  // RBF kernel hyperpriors
  //target += generalized_inverse_gaussian_lpdf(length_scale_f | -1, 2, 1); 
  length_scale_f ~ inv_gamma(5, 5);
  stat_var_f ~ inv_gamma(2, 2);
  
  length_scale_g ~ inv_gamma(5, 5);
  stat_var_g ~ inv_gamma(2, 2);
  
  // Linear kernel hyperparameters
  sigma ~ inv_gamma(2, 2); 
  sigma_b ~ inv_gamma(2, 2);
  
  // Latent variable priors
  eta_f ~ std_normal();
  eta_g ~ std_normal();
  
  // Likelihood
  for(i in 1:(N_real))
    dx[i] ~ normal(f[i]*dt[i], sqrt(exp(g[i]))*sqrt(dt[i]));
}

generated quantities {
  // Generate predictions
  vector[N_pred] f_pred;
  vector[N_pred] g_pred;
  
  if (pred_inf) {
    // predictive inference
    f_pred = gp_pred_rng(x_real, x_pred, N_real, N_pred, f, stat_var_f, length_scale_f, sigma, sigma_b, c, 1);
    g_pred = gp_pred_rng(x_real, x_pred, N_real, N_pred, g, stat_var_g, length_scale_g, sigma, sigma_b, c, 0);
  }
}
