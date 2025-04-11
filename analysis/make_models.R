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

iter <- 1000 # 50% of that is warmup

# Test --------------------------------------------------------------------

# model <- brm(
#   formula = mpg ~ wt,
#   data = mtcars,
#   family = gaussian(),
#   chains = chains_per_task,    # 16 chains per task
#   cores = chains_per_task,     # Use all 16 CPUs
#   iter = 500,
#   seed = 1234 + start_chain,   # Unique seed per task
#   file = paste0("./modeltoy_task_", task_id, ".rds")
# )



# Sample 1 ----------------------------------------------------------------

df <- read.csv("https://raw.githubusercontent.com/RealityBending/FictionChoco/refs/heads/main/data/sample1.csv")
df$Sex <- factor(df$Sex, levels = c("Male", "Female"))
# df <- df[df$Participant %in% unique(df$Participant)[1:10],]

# ZOIB
print("===== ZOIB =====")
print(Sys.time())

f <- bf(
  Real ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  phi ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  zoi ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  coi ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item)
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
  file = paste0("./sample1_zoib_task_", task_id, ".rds")
)




# # # DO NOT DO XBX
#
# # # # XBX
# # # print("===== XBX =====")
# # # f <- bf(
# # #   Real ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
# # #   phi ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
# # #   kappa ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item)
# # # )
# # #
# # # model <- brm(
# # #   formula = f,
# # #   data = df,
# # #   family = xbeta(),
# # #   init = 0,
# # #   chains = chains_per_task,    # 16 chains per task
# # #   cores = chains_per_task,     # Use all 16 CPUs
# # #   iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
# # #   seed = 1234 + start_chain,   # Unique seed per task
# # #   file = paste0("./sample1_xbx_task_", task_id, ".rds")
# # # )
#
# # #


# BEXT
print("===== BEXT =====")
print(Sys.time())

f <- bf(
  Real ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  phi ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  pex ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  bex ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item)
)

model <- brm(
  formula = f,
  data = df,
  family = cogmod::bext(),
  stanvars = cogmod::bext_stanvars(),
  init = 0,
  chains = chains_per_task,    # 16 chains per task
  cores = chains_per_task,     # Use all 16 CPUs
  iter = iter,  # Number of iterations and warmup must be equal (warmup is 50% iter)
  seed = 1234 + start_chain,   # Unique seed per task
  file = paste0("./sample1_bext_task_", task_id, ".rds")
)



# CHOCO
print("===== CHOCO =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  muleft ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  mudelta ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  phileft ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  phidelta ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  pex ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  bex ~ 0 + Intercept + Sex + (0 + Intercept | Participant) + (0 + Intercept | Item),
  pmid = 0,
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- c(
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "muleft", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "mudelta", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "phileft", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "phidelta", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "pex", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Item", dpar = "bex", lb = 0),
  prior("normal(0, 3)", class = "sd", coef = "", group = "Participant", dpar = "", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "muleft", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "mudelta", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "phileft", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "phidelta", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "pex", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "bex", lb = 0),
  prior("normal(0, 0.5)", class = "b", coef = "Intercept", dpar = ""),
  prior("normal(0, 0.3)", class = "b", coef = "Intercept", dpar = "muleft"),
  prior("normal(0, 0.1)", class = "b", coef = "Intercept", dpar = "mudelta"),
  prior("normal(3, 2)", class = "b", coef = "Intercept", dpar = "phileft"),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "phidelta"),
  prior("normal(-3, 3)", class = "b", coef = "Intercept", dpar = "pex"),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "bex"),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = ""),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = "muleft"),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = "mudelta"),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = "phileft"),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = "phidelta"),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = "pex"),
  prior("normal(0, 1)", class = "b", coef = "SexFemale", dpar = "bex")
) |> brms::validate_prior(formula = f, data = df)



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
  file = paste0("./sample1_choco_task_", task_id, ".rds")
)




# CHOCO
print("===== CHOCO - Attractiveness =====")
print(Sys.time())

# Define formula and priors
f <- brms::bf(
  Real ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant) + (0 + Intercept | Item),
  muleft ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant),
  mudelta ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant),
  phileft ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant),
  phidelta ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant),
  pex ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant),
  bex ~ 0 + Intercept + Sex * poly(Attractive, 2) + (0 + Intercept + poly(Attractive, 2) | Participant),
  pmid = 0,
  family = cogmod::choco()
)


# Define priors (brms::default_prior(f, data=df))
priors <- c(
  prior("normal(0, 3)", class = "sd", coef = "", group = "Item", dpar = "", lb = 0),
  prior("normal(0, 3)", class = "sd", coef = "", group = "Participant", dpar = "", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "muleft", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "mudelta", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "phileft", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "phidelta", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "pex", lb = 0),
  prior("normal(0, 2)", class = "sd", coef = "", group = "Participant", dpar = "bex", lb = 0),
  prior("normal(0, 0.5)", class = "b", coef = "Intercept", dpar = ""),
  prior("normal(0, 0.3)", class = "b", coef = "Intercept", dpar = "muleft"),
  prior("normal(0, 0.1)", class = "b", coef = "Intercept", dpar = "mudelta"),
  prior("normal(3, 2)", class = "b", coef = "Intercept", dpar = "phileft"),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "phidelta"),
  prior("normal(0, 3)", class = "b", coef = "Intercept", dpar = "pex"),
  prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "bex"),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = ""),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = "muleft"),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = "mudelta"),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = "phileft"),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = "phidelta"),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = "pex"),
  prior("normal(0, 1)", class = "b", coef = "polyAttractive21", dpar = "bex"),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = ""),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = "muleft"),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = "mudelta"),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = "phileft"),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = "phidelta"),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = "pex"),
  prior("normal(0, 0.5)", class = "b", coef = "polyAttractive22", dpar = "bex"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = ""),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = "muleft"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = "mudelta"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = "phileft"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = "phidelta"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = "pex"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive21", dpar = "bex"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = ""),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = "muleft"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = "mudelta"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = "phileft"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = "phidelta"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = "pex"),
  prior("normal(0, 0.5)", class = "b", coef = "SexFemale:polyAttractive22", dpar = "bex")
) |> brms::validate_prior(formula = f, data = df)


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
  file = paste0("./sample1_chocoattractiveness_task_", task_id, ".rds")
)







# Sample 2 ----------------------------------------------------------------

# df <- read.csv("https://raw.githubusercontent.com/RealityBending/FakeFace/refs/heads/main/data/data.csv")
# df$Real <- (df$Belief_Answer + 1) / 2  # Rescale
# df$Item <- gsub(".jpg", "", df$Stimulus)
# df <- df[c("Participant", "Item", "Real", "Attractive", "Beauty")]
