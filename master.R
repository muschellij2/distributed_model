source("package_contents.R")
#################################################
# Different folder by site
#################################################
library(readr)
all_site_names = paste0("site", 1:3)


synced_folder = "~/Dropbox/Projects/distributed_model"
# which model are we running
model_name = "logistic_example"

L = folder_names(synced_folder)
model_folder = L$model_folder
gradients_folder = L$gradients_folder
converged_folder = L$converged_folder
beta_folder = L$beta_folder



final_file = file.path(converged_folder, paste0(model_name, ".rds"))

run = estimate_new_beta(
  model_name, 
  synced_folder, 
  all_site_names = all_site_names,
  tolerance = 1e-8)
  
if (file.exists(final_file)) {
  formula_file = file.path(model_folder, paste0(model_name, ".rds"))
  if (!file.exists(formula_file)) {
    stop("Need to set up new model, see setup_model")
  } else {
    formula = readr::read_rds(formula_file)
  }
  
  
  res = get_current_beta(model_name, synced_folder)
  beta = res$beta
  iteration_number = res$iteration_number

}
