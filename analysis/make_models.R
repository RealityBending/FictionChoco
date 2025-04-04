# Fit model with 16 chains per task

library(brms)
library(cmdstanr)

# Get array task ID from Slurm
task_id <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

# Number of chains per task from Slurm (fixed at 16)
chains_per_task <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Calculate starting chain ID
start_chain <- (task_id - 1) * chains_per_task + 1



# Test --------------------------------------------------------------------

# model <- brm(
#   formula = mpg ~ wt,
#   data = mtcars,
#   family = gaussian(),
#   chains = chains_per_task,    # 16 chains per task
#   cores = chains_per_task,     # Use all 16 CPUs
#   iter = 2000,
#   warmup = 1000,
#   seed = 1234 + start_chain,   # Unique seed per task
#   backend = "cmdstanr",
#   file = paste0("./modeltoy_task_", task_id, ".rds")
# )



# Sample 1 ----------------------------------------------------------------

df <- read.csv("https://raw.githubusercontent.com/RealityBending/FakeFace/refs/heads/main/data/data.csv")
df$Real <- (df$Belief_Answer + 1) / 2  # Rescale
df$Item <- gsub(".jpg", "", df$Stimulus)
df <- df[c("Participant", "Item", "Real")]


# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + (0 + Intercept | Participant),
  muleft ~ 0 + Intercept,
  phileft ~ 0 + Intercept,
  kleft ~ 0 + Intercept + (0 + Intercept | Participant),
  mud = 0,
  phid = 0,
  kd = 0,
  family = choco7d(link_mu = "logit", link_muleft = "logit", link_mud = "identity",
                   link_phileft = "softplus", link_phid = "identity",
                   link_kleft = "logit", link_kd = "identity")
)


# Define priors (brms::default_prior(f, data=df))
priors <- c(
  # prior("normal(0, 1)", class = "sd", lb = 0),
  prior("normal(0, 1)", class = "sd", coef = "", dpar = "", group = "Participant", lb = 0),
  prior("normal(0, 1)", class = "sd", coef = "", dpar = "kleft", group = "Participant", lb = 0),
  prior("normal(0, 0.5)", class = "b", coef = "Intercept"),
  prior("normal(0, 0.3)", class = "b", coef = "Intercept", dpar = "muleft"),
  prior("normal(3, 0.3)", class = "b", coef = "Intercept", dpar = "phileft"),
  prior("normal(1, 0.5)", class = "b", coef = "Intercept", dpar = "kleft")
) |> brms::validate_prior(formula = f, data = df)


model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = choco7d(),
  stanvars = choco_stanvars("choco7d"),
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = 2000,
  warmup = 1000,
  seed = 1234 + start_chain,   # Unique seed per task
  backend = "cmdstanr",
  file = paste0("./sample1_task_", task_id, ".rds")
)