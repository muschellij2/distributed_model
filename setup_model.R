#################################################
# Different folder by site
#################################################
library(readr)
synced_folder = "~/Dropbox/Projects/distributed_model"
# this structure is the same on all sites
fols = file.path(synced_folder, 
                 c("formulas", "gradients", "models",
                   "betas"))
sapply(fols, dir.create, showWarnings = FALSE)
model_folder = file.path(synced_folder, "formulas")
gradients_folder = file.path(synced_folder, "gradients")
beta_folder = file.path(synced_folder, "betas")
converged_folder = file.path(synced_folder, "models")
# which model are we running
model_name = "logistic_example"


# setup new model
formula = y ~ x1 + x2
family = binomial()
formula_file = file.path(model_folder, paste0(model_name, ".rds"))
if (!file.exists(formula_file)) {
  L = list(formula = formula,
           family = family)
  readr::write_rds(L, formula_file)
}
