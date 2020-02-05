rm(list=ls())
source("package_contents.R")
# really the only thing you'd ever change after the first time
model_name = "logistic_example"
all_site_names = paste0("site", 1:3)

# Fixed stuff by site
site_name = paste0("site1")

for (site_name in all_site_names) {
  # data source is somewhere else
  data_source = paste0("data/logistic_data_", site_name, ".csv")
  site_data = readr::read_csv(data_source)
  
  #################################################
  # Different folder by site
  #################################################
  library(readr)
  synced_folder = "~/Dropbox/Projects/distributed_model"
  
  gradient_file = estimate_site_gradient(
    model_name, 
    synced_folder, 
    site_name = site_name,
    dataset = site_data,
    all_site_names = all_site_names)
    
  # L = folder_names(synced_folder)
  # model_folder = L$model_folder
  # gradients_folder = L$gradients_folder
  # converged_folder = L$converged_folder
  # beta_folder = L$beta_folder
  # 
  # # which model are we running
  # formula_file = file.path(model_folder, paste0(model_name, ".rds"))
  # 
  # # here we can do checks to see if model is done and stuff
  # all_formula_files = list.files(model_folder, pattern = ".rds")
  # # if (!file.exists(formula_file)) {
  # #   stop(paste0("Formula file: ", formula_file, " doesn't exist!",
  # #               " You may need to contact processing site or check your ", 
  # #               "synced_folder"))
  # # } else {
  # #   formula = readr::read_rds(formula_file)
  # # }
  # # 
  # # # here we can do checks to see if model is done and stuff
  # # res = get_current_beta(model_name, synced_folder)
  # # beta = res$beta
  # # iteration_number = res$iteration_number
  # # 
  # gradient_file = estimate_site_gradient(
  #   model_name = model_name, 
  #   synced_folder = synced_folder, 
  #   site_name = site_name, 
  #   dataset = site_data,
  #   iteration_number = iteration_number,
  #   all_site_names = all_site_names)
  # 
  
  
  
}
