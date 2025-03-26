if (!requireNamespace("brms", quietly = TRUE)) {
  stop("Package 'brms' is required to run this script.")
}

# Choice-Confidence (CHOCO) distribution to model bimodal analog/Likert scale data (with zeros and ones).
# The idea is to have two separate ordered beta distributions (see Kubinec) modeling the left and the right hand side of the scale.


# Simulation ---------------------------------------------------------


#' Random simulation from the CHOCO distribution:
#' @param n Number of random draws.
#' @param mu Probability of choosing the right side (relative to the left side).
#' @param delta Distance between the two modes (i.e., the separation/discrimination between the two sides)
#' @param phi Shape parameter of the beta distribution.
#' @param k Cutoff parameter for the beta distribution (the lower this value, the higher the likelihood of extreme responses - zeros or ones)
#' @param muleft The center of the left side beta distribution.
#' @param muright The center of the right side beta distribution.
#' @param phileft The shape parameter of the left side beta distribution.
#' @param phiright The shape parameter of the right side beta distribution.
#' @param kleft The cutoff parameter for the left side beta distribution.
#' @param kright The cutoff parameter for the right side beta distribution.
#'
#' @examples
#' hist(rchoco(3000, mu=0.5, delta=0.5, phi=3, k = 0.95), breaks = 100)
#' hist(rchoco(3000, mu=0.6, delta=0.8, phi=5, k = 0.99), breaks = 100)
rchoco <- function(n, mu = 0.5, delta = 0.5, phi = 3, k = 0.95,
                   muleft = delta, muright = delta, phileft = phi, phiright = phi,
                   kleft = k, kright = k) {
  # Overall probabilities for choosing left or right side
  p_left <- 1 - mu    # probability for the left side
  p_right <- mu       # probability for the right side

  # Compute the discrete endpoint probabilities for 0 (left side) and 1 (right side)
  # p0: probability mass at 0 for the left side
  p0 <- 1 - plogis(qlogis(muleft) - qlogis(1 - kleft))
  # p1: probability mass at 1 for the right side
  p1 <- plogis(qlogis(muright) - qlogis(kright))

  # Pre-allocate output vector for efficiency
  y <- numeric(n)

  # Vectorized random assignment of sides for all n draws
  # Each draw is assigned "left" or "right" based on the specified probabilities
  sides <- sample(c("left", "right"), size = n, replace = TRUE, prob = c(p_left, p_right))

  ## Process draws assigned to the left side
  left_idx <- which(sides == "left")
  if (length(left_idx) > 0) {
    # Generate uniform random numbers to decide if each left-side draw is an endpoint (0) or continuous
    left_u <- runif(length(left_idx))
    # For draws with a uniform value less than p0, assign 0 (discrete mass)
    y[left_idx[left_u < p0]] <- 0
    # For the remaining draws, simulate from the continuous left-side beta distribution
    left_cont_idx <- left_idx[left_u >= p0]
    if (length(left_cont_idx) > 0) {
      # Draw from the beta distribution with parameters based on muleft and phileft
      x0 <- rbeta(length(left_cont_idx), shape1 = muleft * phileft, shape2 = (1 - muleft) * phileft)
      # Transform the beta draw from (0,1) to the interval (0,0.5)
      y[left_cont_idx] <- 0.5 - x0 / 2
    }
  }

  ## Process draws assigned to the right side
  right_idx <- which(sides == "right")
  if (length(right_idx) > 0) {
    # Generate uniform random numbers to decide if each right-side draw is an endpoint (1) or continuous
    right_u <- runif(length(right_idx))
    # For draws with a uniform value less than p1, assign 1 (discrete mass)
    y[right_idx[right_u < p1]] <- 1
    # For the remaining draws, simulate from the continuous right-side beta distribution
    right_cont_idx <- right_idx[right_u >= p1]
    if (length(right_cont_idx) > 0) {
      # Draw from the beta distribution with parameters based on muright and phiright
      x1 <- rbeta(length(right_cont_idx), shape1 = muright * phiright, shape2 = (1 - muright) * phiright)
      # Transform the beta draw from (0,1) to the interval (0.5,1)
      y[right_cont_idx] <- 0.5 + x1 / 2
    }
  }

  y
}



# Stan --------------------------------------------------------------------



# Create a stanvars object to pass the custom functions to brms
stanvars_choco <- brms::stanvar(scode = "
real choco_lpdf(real y, real mu, real muleft, real phileft,
                real muright, real phiright, real kleft, real kright) {
  real eps = 1e-8;
  real p_left = 1 - mu;
  real p_right = mu;

  // Handle discrete mass at 0
  if (y < eps) {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    return log(p_left) + log(fmax(p0, eps));
  }
  // Handle discrete mass at 1
  else if (y > (1 - eps)) {
    real p1 = inv_logit(logit(muright) - logit(kright));
    return log(p_right) + log(fmax(p1, eps));
  }
  // Handle continuous part, with special case for y=0.5
  else {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    real p1 = inv_logit(logit(muright) - logit(kright));
    real log1m_p0 = log1m(fmin(p0, 1 - eps));
    real log1m_p1 = log1m(fmin(p1, 1 - eps));

    // Check if y is approximately 0.5 (within eps)
    if (abs(y - 0.5) < eps) {
      // For y=0.5, compute densities approaching 0.5 from left and right
      real x0 = 2 * (0.5 - (0.5 - eps));  // x0 approaches 0+ as y approaches 0.5-
      real x1 = 2 * ((0.5 + eps) - 0.5);  // x1 approaches 0+ as y approaches 0.5+
      real dens_left = p_left * exp(log1m_p0 + log(2) +
                                    beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft));
      real dens_right = p_right * exp(log1m_p1 + log(2) +
                                      beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright));
      // Average the left and right densities
      real avg_dens = (dens_left + dens_right) / 2;
      return log(fmax(avg_dens, eps));
    }
    else if (y < 0.5) {
      real x0 = 2 * (0.5 - y);
      return log(p_left) + log1m_p0 + log(2) +
             beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft);
    }
    else {  // y > 0.5
      real x1 = fmax(2 * (y - 0.5), eps);
      return log(p_right) + log1m_p1 + log(2) +
             beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright);
    }
  }
}
", block = "functions")



stanvars_chocomini <- brms::stanvar(scode = "
real chocomini_lpdf(real y, real mu, real delta, real phi, real k) {
  real eps = 1e-8;
  real d = 0.5 * delta;  // delta in (0,1), d in (0,0.5)
  real muleft = 2 * d;
  real muright = 2 * d;
  real phileft = phi;
  real phiright = phi;
  real kleft = k;
  real kright = k;
  real p_left = 1 - mu;
  real p_right = mu;

  if (y < eps) {
    real p0 = fmax(1 - inv_logit(logit(muleft) - logit(1 - kleft)), eps);
    return log(p_left) + log(fmax(p0, eps));
  }
  else if (y > (1 - eps)) {
    real p1 = fmax(inv_logit(logit(muright) - logit(kright)), eps);
    return log(p_right) + log(fmax(p1, eps));
  }
  else {
    real p0 = fmax(1 - inv_logit(logit(muleft) - logit(1 - kleft)), eps);
    real p1 = fmin(inv_logit(logit(muright) - logit(kright)), 1 - eps);
    real log1m_p0 = log1m(fmax(1 - p0, eps));  // Ensures we never take log(0)
    real log1m_p1 = log1m(fmax(1 - p1, eps));

    if (abs(y - 0.5) < eps) {
      real x0 = eps;
      real x1 = eps;
      real dens_left = p_left * exp(log1m_p0 + log(2) +
                                    beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft));
      real dens_right = p_right * exp(log1m_p1 + log(2) +
                                      beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright));
      real avg_dens = (dens_left + dens_right) / 2;
      return log(fmax(avg_dens, eps));
    }
    else if (y < 0.5) {
      real x0 = 2 * (0.5 - y);
      return log(p_left) + log1m_p0 + log(2) +
             beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft);
    }
    else {
      real x1 = fmax(2 * (y - 0.5), eps);
      return log(p_right) + log1m_p1 + log(2) +
             beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright);
    }
  }
}
", block = "functions")




# Define custom family
choco <- brms::custom_family(
  "choco",
  dpars = c("mu", "muleft", "muright", "phileft", "phiright", "kleft", "kright"),
  links = c("logit", "logit", "logit", "softplus", "softplus", "logit", "logit")
)


chocomini <- brms::custom_family(
  "chocomini",
  dpars = c("mu", "delta", "phi", "k"),
  links = c("logit", "logit", "softplus", "logit")
)




# Predict -----------------------------------------------------------------


# Posterior predict function
posterior_predict_choco <- function(i, prep, ...) {
  # Extract posterior draws of the distributional parameters
  mu       <- brms::get_dpar(prep, "mu", i = i)
  muleft   <- brms::get_dpar(prep, "muleft", i = i)
  muright  <- brms::get_dpar(prep, "muright", i = i)
  phileft  <- brms::get_dpar(prep, "phileft", i = i)
  phiright <- brms::get_dpar(prep, "phiright", i = i)
  kleft    <- brms::get_dpar(prep, "kleft", i = i)
  kright   <- brms::get_dpar(prep, "kright", i = i)

  # Use mapply to vectorize over the posterior draws.
  # For each draw (i.e., each set of parameter values), simulate one response using rchoco.
  yrep <- mapply(function(mu_i, muleft_i, muright_i, phileft_i, phiright_i, kleft_i, kright_i) {
    rchoco(1,
           mu       = mu_i,
           muleft   = muleft_i,
           muright  = muright_i,
           phileft  = phileft_i,
           phiright = phiright_i,
           kleft    = kleft_i,
           kright   = kright_i)
  }, mu, muleft, muright, phileft, phiright, kleft, kright)

  yrep
}


posterior_predict_chocomini <- function(i, prep, ...) {
  # Extract posterior draws for the mini version parameters
  mu    <- brms::get_dpar(prep, "mu", i = i)
  delta <- brms::get_dpar(prep, "delta", i = i)
  phi   <- brms::get_dpar(prep, "phi", i = i)
  k     <- brms::get_dpar(prep, "k", i = i)

  # Vectorize over the draws: For each posterior draw, simulate one response.
  yrep <- mapply(function(mu_i, delta_i, phi_i, k_i) {
    rchoco(1, mu = mu_i, delta = delta_i, phi = phi_i, k = k_i)
  }, mu, delta, phi, k)

  yrep
}

