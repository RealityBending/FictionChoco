# Choice-Confidence (CHOCO) distribution to model bimodal analog/Likert scale data (with zeros and ones).
# The idea is to have two separate ordered beta distributions (see Kubinec, 2023) modeling the left and the right hand side of the scale.

# Regarding Ordered Beta distributions, see:
# - Stan code: https://github.com/saudiwin/ordbetareg/blob/master/beta_logit.stan
# - Python code: https://github.com/saudiwin/ordbetareg_py/blob/main/ordbetareg/model.py

# Simulation ---------------------------------------------------------

#' Random simulation from the CHOCO distribution:
#' @param n Number of random draws.
#' @param mu Probability of choosing the right side (relative to the left side).
#' @param muleft The center of the left side beta distribution.
#' @param phileft The shape parameter of the left side beta distribution.
#' @param kleft The cutoff parameter for the left side beta distribution.
#' @param muright The center of the right side beta distribution (overridden by mud if NULL).
#' @param phiright The shape parameter of the right side beta distribution (overridden by phid if NULL).
#' @param kright The cutoff parameter for the right side beta distribution (overridden by kd if NULL).
#' @param mud Deviation for muright on the logit scale relative to muleft (used if muright is NULL).
#' @param phid Deviation for phiright as a log-multiplier relative to phileft (used if phiright is NULL).
#' @param kd Deviation for kright on the logit scale relative to kleft (used if kright is NULL).
#'
#' @examples
#' hist(rchoco(3000, mu=0.5, muleft=0.5, phileft=3, kleft=0.95), breaks=100)
#' hist(rchoco(3000, mu=0.6, muleft=0.5, phileft=5, kleft=0.99), breaks=100)
#' hist(rchoco(3000, mu=0.4, muleft=0.5, phileft=3, kleft=0.85), breaks=100)
rchoco <- function(n, mu = 0.5, muleft = 0.5, phileft = 3, kleft = 0.95,
                   muright = NULL, phiright = NULL, kright = NULL,
                   mud = 0, phid = 0, kd = 0) {
  # Overall probabilities for choosing left or right side
  p_left <- 1 - mu    # probability for the left side
  p_right <- mu       # probability for the right side

  # Use provided right-side parameters if not NULL, else compute using deviation parameters:
  if (is.null(muright)) muright <- plogis(qlogis(muleft) + mud)
  if (is.null(phiright)) phiright <- phileft * exp(phid)
  if (is.null(kright)) kright <- plogis(qlogis(kleft) + kd)

  # Compute the discrete endpoint probabilities for 0 (left side) and 1 (right side)
  p0 <- 1 - plogis(qlogis(muleft) - qlogis(1 - kleft))      # left-side endpoint probability
  p1 <- plogis(qlogis(muright) - qlogis(kright))              # right-side endpoint probability

  # Pre-allocate output vector for efficiency
  y <- numeric(n)

  # Vectorized random assignment of sides for all n draws
  sides <- sample(c("left", "right"), size = n, replace = TRUE, prob = c(p_left, p_right))

  ## Process draws assigned to the left side
  left_idx <- which(sides == "left")
  if (length(left_idx) > 0) {
    left_u <- runif(length(left_idx))
    y[left_idx[left_u < p0]] <- 0
    left_cont_idx <- left_idx[left_u >= p0]
    if (length(left_cont_idx) > 0) {
      x0 <- rbeta(length(left_cont_idx),
                  shape1 = muleft * phileft,
                  shape2 = (1 - muleft) * phileft)
      y[left_cont_idx] <- 0.5 - x0 / 2
    }
  }

  ## Process draws assigned to the right side
  right_idx <- which(sides == "right")
  if (length(right_idx) > 0) {
    right_u <- runif(length(right_idx))
    y[right_idx[right_u < p1]] <- 1
    right_cont_idx <- right_idx[right_u >= p1]
    if (length(right_cont_idx) > 0) {
      x1 <- rbeta(length(right_cont_idx),
                  shape1 = muright * phiright,
                  shape2 = (1 - muright) * phiright)
      y[right_cont_idx] <- 0.5 + x1 / 2
    }
  }

  y
}

# Stan --------------------------------------------------------------------
# Create a stanvars object to pass the custom functions to brms
choco_stanvars <- function(type = "choco7") {
  if (type == "choco7") { # 7-parameter model with independent right-side parameters
    stancode <- brms::stanvar(scode = "
real choco7_lpdf(real y, real mu, real muleft, real phileft,
                real muright, real phiright, real kleft, real kright) {
  real eps = 1e-8;
  real p_left = 1 - mu;
  real p_right = mu;
  if (y < eps) {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    return log(p_left) + log(fmax(p0, eps));
  } else if (y > (1 - eps)) {
    real p1 = inv_logit(logit(muright) - logit(kright));
    return log(p_right) + log(fmax(p1, eps));
  } else {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    real p1 = inv_logit(logit(muright) - logit(kright));
    real log1m_p0 = log1m(fmin(p0, 1 - eps));
    real log1m_p1 = log1m(fmin(p1, 1 - eps));
    if (abs(y - 0.5) < eps) {
      real x0 = 2 * (0.5 - (0.5 - eps));
      real x1 = 2 * ((0.5 + eps) - 0.5);
      real dens_left = p_left * exp(log1m_p0 + log(2) +
                                    beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft));
      real dens_right = p_right * exp(log1m_p1 + log(2) +
                                      beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright));
      real avg_dens = (dens_left + dens_right) / 2;
      return log(fmax(avg_dens, eps));
    } else if (y < 0.5) {
      real x0 = 2 * (0.5 - y);
      return log(p_left) + log1m_p0 + log(2) +
             beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft);
    } else {
      real x1 = fmax(2 * (y - 0.5), eps);
      return log(p_right) + log1m_p1 + log(2) +
             beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright);
    }
  }
}
", block = "functions")
  } else if (type == "choco7d") { # 7-parameter deviation-based model
    stancode <- brms::stanvar(scode = "
real choco7d_lpdf(real y, real mu, real muleft, real mud, real phileft, real phid, real kleft, real kd) {
  real eps = 1e-8;
  real p_left = 1 - mu;
  real p_right = mu;
  // Derive right-side parameters as deviations from the left-side:
  real muright = inv_logit(logit(muleft) + mud);
  real phiright = phileft * exp(phid);
  real kright = inv_logit(logit(kleft) + kd);
  if (y < eps) {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    return log(p_left) + log(fmax(p0, eps));
  } else if (y > (1 - eps)) {
    real p1 = inv_logit(logit(muright) - logit(kright));
    return log(p_right) + log(fmax(p1, eps));
  } else {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    real p1 = inv_logit(logit(muright) - logit(kright));
    real log1m_p0 = log1m(fmin(p0, 1 - eps));
    real log1m_p1 = log1m(fmin(p1, 1 - eps));
    if (abs(y - 0.5) < eps) {
      real x0 = 2 * (0.5 - (0.5 - eps));
      real x1 = 2 * ((0.5 + eps) - 0.5);
      real dens_left = p_left * exp(log1m_p0 + log(2) +
                                    beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft));
      real dens_right = p_right * exp(log1m_p1 + log(2) +
                                      beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright));
      real avg_dens = (dens_left + dens_right) / 2;
      return log(fmax(avg_dens, eps));
    } else if (y < 0.5) {
      real x0 = 2 * (0.5 - y);
      return log(p_left) + log1m_p0 + log(2) +
             beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft);
    } else {
      real x1 = fmax(2 * (y - 0.5), eps);
      return log(p_right) + log1m_p1 + log(2) +
             beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright);
    }
  }
}
", block = "functions")
  } else {
    stop("Invalid type. Choose either 'choco7' or 'choco7d'.")
  }
  stancode
}

# Define custom families
choco7 <- function(link_mu = "logit", link_muleft = "logit", link_muright = "logit",
                   link_phileft = "softplus", link_phiright = "softplus",
                   link_kleft = "logit", link_kright = "logit") {
  brms::custom_family(
    "choco7",
    dpars = c("mu", "muleft", "muright", "phileft", "phiright", "kleft", "kright"),
    links = c(link_mu, link_muleft, link_muright,
              link_phileft, link_phiright,
              link_kleft, link_kright)
  )
}

choco7d <- function(link_mu = "logit", link_muleft = "logit", link_mud = "identity",
                    link_phileft = "softplus", link_phid = "identity",
                    link_kleft = "logit", link_kd = "identity") {
  brms::custom_family(
    "choco7d",
    dpars = c("mu", "muleft", "mud", "phileft", "phid", "kleft", "kd"),
    links = c(link_mu, link_muleft, link_mud,
              link_phileft, link_phid,
              link_kleft, link_kd)
  )
}
# Predict -----------------------------------------------------------------

# Posterior predict function for choco7
posterior_predict_choco7 <- function(i, prep, ...) {
  mu       <- brms::get_dpar(prep, "mu", i = i)
  muleft   <- brms::get_dpar(prep, "muleft", i = i)
  muright  <- brms::get_dpar(prep, "muright", i = i)
  phileft  <- brms::get_dpar(prep, "phileft", i = i)
  phiright <- brms::get_dpar(prep, "phiright", i = i)
  kleft    <- brms::get_dpar(prep, "kleft", i = i)
  kright   <- brms::get_dpar(prep, "kright", i = i)

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

# Posterior predict function for choco7d (deviation-based)
posterior_predict_choco7d <- function(i, prep, ...) {
  mu      <- brms::get_dpar(prep, "mu", i = i)
  muleft  <- brms::get_dpar(prep, "muleft", i = i)
  mud     <- brms::get_dpar(prep, "mud", i = i)
  phileft <- brms::get_dpar(prep, "phileft", i = i)
  phid    <- brms::get_dpar(prep, "phid", i = i)
  kleft   <- brms::get_dpar(prep, "kleft", i = i)
  kd      <- brms::get_dpar(prep, "kd", i = i)

  yrep <- mapply(function(mu_i, muleft_i, mud_i, phileft_i, phid_i, kleft_i, kd_i) {
    rchoco(1,
           mu      = mu_i,
           muleft  = muleft_i,
           mud     = mud_i,
           phileft = phileft_i,
           phid    = phid_i,
           kleft   = kleft_i,
           kd      = kd_i)
  }, mu, muleft, mud, phileft, phid, kleft, kd)

  yrep
}
