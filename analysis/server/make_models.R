# Fit model with 16 chains per task

library(brms)
library(cmdstanr)
library(cogmod)

options(brms.backend = "cmdstanr")

# Get array task ID from Slurm
task_id <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

# Number of chains per task from Slurm (fixed at 16)
chains_per_task <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

# Calculate starting chain ID
start_chain <- (task_id - 1) * chains_per_task + 1

iter <- 1500 # 50% of that is warmup
path <- "/mnt/lustre/users/psych/dmm56/FictionChoco/"

# Test --------------------------------------------------------------------

# model <- brm(
#   formula = mpg ~ wt,
#   data = mtcars,
#   family = gaussian(),
#   chains = chains_per_task,    # 16 chains per task
#   cores = chains_per_task,     # Use all 16 CPUs
#   iter = 500,
#   seed = 1234 + start_chain,   # Unique seed per task
#   file = paste0(path, "modeltoy_task_", task_id, ".rds")
# )


# Sample 1 ----------------------------------------------------------------


df <- read.csv("https://raw.githubusercontent.com/RealityBending/FictionChoco/refs/heads/main/data/sample1.csv")
df$Sex <- factor(df$Sex, levels = c("Male", "Female"))




# ZOIB
print("===== ZOIB =====")
print(Sys.time())

f <- bf(
  Real ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  phi ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  zoi ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  coi ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item)
)

model <- brm(
  formula = f,
  data = df,
  family = zero_one_inflated_beta(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample1_zoib_task_", task_id, ".rds")
)

# # # # # DO NOT DO XBX
# # #
# # # # # # XBX
# # # # # print("===== XBX =====")
# # # # # f <- bf(
# # # # #   Real ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
# # # # #   phi ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
# # # # #   kappa ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item)
# # # # # )
# # # # #
# # # # # model <- brm(
# # # # #   formula = f,
# # # # #   data = df,
# # # # #   family = xbeta(),
# # # # #   init = 0,
# # # # #   chains = chains_per_task,    # 16 chains per task
# # # # #   cores = chains_per_task,     # Use all 16 CPUs
# # # # #   iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
# # # # #   seed = 1234 + start_chain,   # Unique seed per task
# # # # #   file = paste0(path, "sample1_xbx_task_", task_id, ".rds")
# # # # # )
# # #
# # # # #


# BEXT
print("===== BETAGATE =====")
print(Sys.time())

f <- bf(
  Real ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  phi ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  pex ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  bex ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item)
)

model <- brm(
  formula = f,
  data = df,
  family = cogmod::betagate(),
  stanvars = cogmod::betagate_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample1_betagate_task_", task_id, ".rds")
)



# CHOCO
print("===== CHOCO =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  confright ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  confleft ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  precright ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  precleft ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  pex ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  bex ~ 0 + Intercept + Sex + (1 | Participant) + (1 | Item),
  pmid = 0,
  family = cogmod::choco()
)

model <- brm(
  formula = f,
  data = df,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample1_choco_task_", task_id, ".rds")
)



# CHOCO
print("===== CHOCO - Attractiveness =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex / poly(Attractive, 2) + (poly(Attractive, 2) | Participant) + (1 | Item),
  confright ~ 0 + Intercept + Sex / poly(Attractive, 2) + (poly(Attractive, 2) | Participant),
  confleft ~ 0 + Intercept + Sex / poly(Attractive, 2) + (poly(Attractive, 2) | Participant),
  precright ~ 0 + Intercept + Sex / poly(Attractive, 2) + (poly(Attractive, 2) | Participant),
  precleft ~ 0 + Intercept + Sex / poly(Attractive, 2) + (poly(Attractive, 2) | Participant),
  pex ~ 0 + Intercept + Sex / poly(Attractive, 2) + (poly(Attractive, 2) | Participant),
  bex ~ 0 + Intercept + Sex / poly(Attractive, 2),
  pmid = 0,
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors_base <- c(
  # Variances
  prior("normal(0.5, 0.2)", class = "sd", coef = "Intercept", group = "Item", dpar = ""),
  prior("normal(0.5, 0.2)", class = "sd", coef = "Intercept", group = "Participant", dpar = ""),
  prior("normal(1, 0.2)", class = "sd", coef = "Intercept", group = "Participant", dpar = "confright"),
  prior("normal(1, 0.2)", class = "sd", coef = "Intercept", group = "Participant", dpar = "confleft"),
  prior("normal(2.5, 0.5)", class = "sd", coef = "Intercept", group = "Participant", dpar = "precright"),
  prior("normal(2.5, 0.5)", class = "sd", coef = "Intercept", group = "Participant", dpar = "precleft"),
  prior("normal(0.5, 0.2)", class = "sd", coef = "Intercept", group = "Participant", dpar = "pex"),

  # Intercepts
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = ""),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "confright"),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "confleft"),
  prior("normal(3, 1)", class = "b", coef = "Intercept", dpar = "precright"),
  prior("normal(3, 1)", class = "b", coef = "Intercept", dpar = "precleft"),
  prior("normal(-2, 1)", class = "b", coef = "Intercept", dpar = "pex"),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "bex")
)

# Random
priors <- priors_base
for(p in c("", "confright", "confleft", "precright", "precleft", "pex")) {
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant", coef = "polyAttractive21", dpar = p))
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant", coef = "polyAttractive22", dpar = p))
}

# Coefficients
for(p in c("", "confright", "confleft", "precright", "precleft", "pex", "bex")) {
  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 2, 1), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale", dpar = p))

  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 10, 5), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyAttractive22", dpar = p))
}

priors <- brms::validate_prior(priors, formula = f, data = df)



model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample1_chocoattractiveness_task_", task_id, ".rds")
)



print("===== CHOCO - Beauty =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex / poly(Beauty, 2) + (poly(Beauty, 2) | Participant) + (1 | Item),
  confright ~ 0 + Intercept + Sex / poly(Beauty, 2) + (poly(Beauty, 2) | Participant),
  confleft ~ 0 + Intercept + Sex / poly(Beauty, 2) + (poly(Beauty, 2) | Participant),
  precright ~ 0 + Intercept + Sex / poly(Beauty, 2) + (poly(Beauty, 2) | Participant),
  precleft ~ 0 + Intercept + Sex / poly(Beauty, 2) + (poly(Beauty, 2) | Participant),
  pex ~ 0 + Intercept + Sex / poly(Beauty, 2) + (poly(Beauty, 2) | Participant),
  bex ~ 0 + Intercept + Sex / poly(Beauty, 2),
  pmid = 0,
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- priors_base
# Random
for(p in c("", "confright", "confleft", "precright", "precleft", "pex")) {
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant", coef = "polyBeauty21", dpar = p))
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant", coef = "polyBeauty22", dpar = p))
}

# Coefficients
for(p in c("", "confright", "confleft", "precright", "precleft", "pex", "bex")) {
  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 2, 1), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale", dpar = p))

  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 10, 5), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyBeauty22", dpar = p))
}

priors <- brms::validate_prior(priors, formula = f, data = df)



model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample1_chocobeauty_task_", task_id, ".rds")
)



# Sample 2 ----------------------------------------------------------------

df <- read.csv("https://raw.githubusercontent.com/RealityBending/FictionChoco/refs/heads/main/data/sample2.csv")
df$Sex <- factor(df$Sex, levels = c("Male", "Female"))
df$Condition <- factor(df$Condition, levels = c("Photograph", "AI-Generated"))





print("===== BETAGATE - Sample 2 (attractiveness No-extreme) =====")
print(Sys.time())

f <- bf(
  RealNoex ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant) + (Condition | Item),
  phi ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  pex = 0,
  bex = 0.5
)

model <- brm(
  formula = f,
  data = df,
  family = cogmod::betagate(),
  stanvars = cogmod::betagate_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample2noex_betagate_task_", task_id, ".rds")
)



# CHOCO
print("===== CHOCO - Sample 2 (attractiveness No-extreme) =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  RealNoex ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant) + (Condition | Item),
  confright ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  confleft ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  precright ~ 0 + Intercept + Sex  / poly(Attractive, 2) * Condition + (Condition | Participant),
  precleft ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  pex = 0,
  bex = 0.5,
  pmid ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- priors_base[c(1:6, 8:12),]
priors <- c(priors,
            prior("normal(0, 1)", class = "sd", coef = "Intercept", group ="Participant", dpar = "pmid"),
            prior("normal(-2, 1)", class = "b", coef = "Intercept", dpar = "pmid"))
# Random
for(p in c("", "confright", "confleft", "precright", "precleft")) {
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant",
                                   coef = "ConditionAIMGenerated", dpar = p))
}

# Coefficients
for(p in c("", "confright", "confleft", "precright", "precleft", "pmid")) {
  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 2, 1), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "ConditionAIMGenerated", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated", dpar = p))

  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 10, 5), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyAttractive22", dpar = p))
}


priors <- brms::validate_prior(priors, formula = f, data = df)



model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample2noex_chocoattractiveness_task_", task_id, ".rds")
)




print("===== CHOCO - Sample 2 (Beauty No-extreme) =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  RealNoex ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant) + (Condition | Item),
  confright ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  confleft ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  precright ~ 0 + Intercept + Sex  / poly(Beauty, 2) * Condition + (Condition | Participant),
  precleft ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  pex = 0,
  bex = 0.5,
  pmid ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- priors_base[c(1:6, 8:12),]
priors <- c(priors,
            prior("normal(0, 1)", class = "sd", coef = "Intercept", group ="Participant", dpar = "pmid"),
            prior("normal(-2, 1)", class = "b", coef = "Intercept", dpar = "pmid"))
# Random
for(p in c("", "confright", "confleft", "precright", "precleft")) {
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant",
                                   coef = "ConditionAIMGenerated", dpar = p))
}

# Coefficients
for(p in c("", "confright", "confleft", "precright", "precleft", "pmid")) {
  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 2, 1), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "ConditionAIMGenerated", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated", dpar = p))

  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 10, 5), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyBeauty22", dpar = p))
}


priors <- brms::validate_prior(priors, formula = f, data = df)



model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample2noex_chocobeauty_task_", task_id, ".rds")
)










print("===== BETAGATE - Sample 2 (attractiveness) =====")
print(Sys.time())

f <- bf(
  Real ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant) + (Condition | Item),
  phi ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  pex ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  bex ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition
)

model <- brm(
  formula = f,
  data = df,
  family = cogmod::betagate(),
  stanvars = cogmod::betagate_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample2_betagate_task_", task_id, ".rds")
)


# CHOCO
print("===== CHOCO - Sample 2 (attractiveness) =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant) + (Condition | Item),
  confright ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  confleft ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  precright ~ 0 + Intercept + Sex  / poly(Attractive, 2) * Condition + (Condition | Participant),
  precleft ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  pex ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  bex ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition,
  pmid ~ 0 + Intercept + Sex / poly(Attractive, 2) * Condition + (Condition | Participant),
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- c(priors_base, prior("normal(-2, 1)", class = "b", coef = "Intercept", dpar = "pmid"))
# Random
for(p in c("", "confright", "confleft", "precright", "precleft", "pex")) {
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant",
                                   coef = "ConditionAIMGenerated", dpar = p))
}

# Coefficients
for(p in c("", "confright", "confleft", "precright", "precleft", "pex", "bex", "pmid")) {
  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 2, 1), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "ConditionAIMGenerated", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated", dpar = p))

  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 10, 5), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyAttractive22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyAttractive21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyAttractive22", dpar = p))
}


priors <- brms::validate_prior(priors, formula = f, data = df)



model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample2_chocoattractiveness_task_", task_id, ".rds")
)










# CHOCO
print("===== CHOCO - Sample 2 (Beauty) =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant) + (Condition | Item),
  confright ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  confleft ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  precright ~ 0 + Intercept + Sex  / poly(Beauty, 2) * Condition + (Condition | Participant),
  precleft ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  pex ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  bex ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition,
  pmid ~ 0 + Intercept + Sex / poly(Beauty, 2) * Condition + (Condition | Participant),
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- c(priors_base, prior("normal(-2, 1)", class = "b", coef = "Intercept", dpar = "pmid"))
# Random
for(p in c("", "confright", "confleft", "precright", "precleft", "pex")) {
  priors <- c(priors, prior_string("normal(0, 0.5)", class = "sd", group = "Participant",
                                   coef = "ConditionAIMGenerated", dpar = p))
}

# Coefficients
for(p in c("", "confright", "confleft", "precright", "precleft", "pex", "bex", "pmid")) {
  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 2, 1), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "ConditionAIMGenerated", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated", dpar = p))

  d <- paste0("normal(0, ", ifelse(p %in% c("precright", "precleft"), 10, 5), ")")
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexMale:ConditionAIMGenerated:polyBeauty22", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyBeauty21", dpar = p))
  priors <- c(priors, prior_string(d, class = "b", coef = "SexFemale:ConditionAIMGenerated:polyBeauty22", dpar = p))
}


priors <- brms::validate_prior(priors, formula = f, data = df)



model <- brm(
  formula = f,
  data = df,
  prior = priors,
  family = cogmod::choco(),
  stanvars = cogmod::choco_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0(path, "sample2_chocobeauty_task_", task_id, ".rds")
)


