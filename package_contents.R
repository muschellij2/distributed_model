
gradient_value = function(beta = NULL, df, formula, 
                          family = binomial(), iteration_number = 0,
                          shuffle_rows = TRUE) {
  if (shuffle_rows) {
    df = df[sample(nrow(df)), ]
  }
  y = model.frame(formula, data = df)[,1]
  X = model.matrix(formula, data = df)
  cc = complete.cases(X)
  X = X[cc,]
  y = y[cc]
  
  if (is.null(beta)) {
    beta = rep(0, ncol(X))
  }
  linkinv = family$linkinv
  variance <- family$variance
  p = linkinv(X %*% beta)
  # expb = exp(X %*% beta)
  # p = expb / (1 + expb)
  p = c(p)
  gradient = colSums(X * (y - p))
  stopifnot(length(gradient) == length(beta))
  result = list(
    gradient = gradient,
    sample_size = nrow(X),
    iteration_number = iteration_number
  )
  return(result)
}


use_glm_gradient_value = function(
  beta = NULL, df, formula, 
  family = binomial(), iteration_number = 0,
  shuffle_rows = TRUE) {
  
  if (shuffle_rows) {
    df = df[sample(nrow(df)), ]
  }  
  X = model.matrix(formula, data = df)
  if (is.null(beta)) {
    beta = rep(0, ncol(X))
  }  
  start = beta
  print(start)
  mod = glm(
    formula = formula,
    data = df,
    family = family,
    start = start,
    control = list(maxit = 1))
  return(mod)
}



folder_names = function(synced_folder) {
  L = list(
    # this structure is the same on all sites
    model_folder = file.path(synced_folder, "formulas"),
    gradients_folder = file.path(synced_folder, "gradients"),
    beta_folder = file.path(synced_folder, "betas"),
    converged_folder = file.path(synced_folder, "models")
  )
  return(L)
}


master_beta_file = function(model_name, synced_folder) {
  file_list = folder_names(synced_folder)
  beta_folder = file_list$beta_folder
  
  all_beta_files = list.files(
    beta_folder, 
    pattern = paste0("^", model_name, "-iteration.*.rds"),
    full.names = TRUE)
  if (length(all_beta_files) == 0) {
    beta = NULL
    iteration_number = 1
  } else {
    beta_number = sub(".*iteration(.*)[.]rds", "\\1", 
                      basename(all_beta_files))
    beta_number = as.numeric(beta_number)
    beta_list = read_rds(all_beta_files[ which.max(beta_number)])
    beta = beta_list$beta
    iteration_number = beta_list$iteration_number_next
  }
  out_beta_file = file.path(
    beta_folder,
    paste0(
      model_name, 
      sprintf("-iteration%04.0f", iteration_number),
      ".rds")
  )
  return(out_beta_file)
}


aggregate_gradients = function(
  all_gradient_files,
  iteration_number) {
  gradient_list = lapply(all_gradient_files, readr::read_rds)
  names(gradient_list) = all_gradient_files
  
  sum_grads = sapply(gradient_list, function(x) x$gradient)
  if (is.null(beta)) {
    beta = rep(0, nrow(sum_grads))
  }
  stopifnot(is.matrix(sum_grads))
  grad = rowSums(sum_grads)
  
  iter_nums = sapply(gradient_list, function(x) x$iteration_number)
  stopifnot(is.vector(iter_nums))
  stopifnot(all(iter_nums == iteration_number))
  
  ss = sapply(gradient_list, function(x) x$sample_size)
  stopifnot(is.vector(ss))
  n = sum(ss)
  
  grad = grad / n
  result = list(
    gradient = grad,
    total_sample_size = n)
  return(result)
}


get_current_beta = function(model_name, synced_folder) {
  file_list = folder_names(synced_folder)
  beta_folder = file_list$beta_folder
  all_beta_files = list.files(
    beta_folder, 
    pattern = paste0("^", model_name, "-iteration.*.rds"),
    full.names = TRUE)
  
  if (length(all_beta_files) == 0) {
    beta = NULL
    iteration_number = 1
  } else {
    beta_number = sub(".*iteration(.*)[.]rds", "\\1", 
                      basename(all_beta_files))
    beta_number = as.numeric(beta_number)
    beta_list = read_rds(all_beta_files[ which.max(beta_number)])
    beta = beta_list$beta
    iteration_number = beta_list$iteration_number_next
  }
  L = list(
    iteration_number = iteration_number
  )
  L$beta =  beta
  return(L)
}


estimate_site_gradient = function(
  model_name, synced_folder, 
  site_name = "site1", dataset,
  all_site_names = paste0("site", 1:3),
  shuffle_rows = TRUE) {
  
  
  
  site_name = match.arg(site_name, choices = all_site_names)
  file_list = folder_names(synced_folder)
  gradients_folder = file_list$gradients_folder
  model_folder = file_list$model_folder
  
  # which model are we running
  formula_file = file.path(model_folder, 
                           paste0(model_name, ".rds"))
  
  if (!file.exists(formula_file)) {
    stop(paste0("Formula file: ", formula_file, " doesn't exist!",
                " You may need to contact processing site or check your ", 
                "synced_folder"))
  } else {
    formula_list = readr::read_rds(formula_file)
    formula = formula_list$formula
    family = formula_list$family
    if (is.character(family)) {
      family = get(family, envir = .BaseNamespaceEnv)
    }
    if (is.function(family)) {
      family = family()
    }    
    if (!inherits(family, "family")) {
      stop("family specified is not a family object - see setup_model")
    }    
  }
  
  res = get_current_beta(model_name, synced_folder)
  beta = res$beta
  iteration_number = res$iteration_number  
  
  gradient_file = file.path(
    gradients_folder, 
    paste0(model_name, "-", 
           site_name, 
           sprintf("-iteration%04.0f", iteration_number),
           ".rds"))
  all_gradient_files = file.path(
    gradients_folder, 
    paste0(model_name, "-", 
           all_site_names,
           sprintf("-iteration%04.0f", iteration_number),
           ".rds"))
  # here we would simply wait
  # should check if converged
  if (file.exists(gradient_file)) {
    if (!all(file.exists(all_gradient_files))) {
      print("Waiting for other sites to create gradients")
    } else {
      print("Waiting for compute site to create new betas")
    }
  } else {
    print(paste0("Creating Gradient, iteration ", 
                 iteration_number))
    use_glm_gradient_value(beta = beta,
                           df = dataset,
                           formula = formula,
                           family = family,
                           iteration_number = iteration_number,
                           shuffle_rows = shuffle_rows)
    grad = gradient_value(beta = beta, 
                          df = dataset, 
                          formula = formula, 
                          family = family,
                          iteration_number = iteration_number,
                          shuffle_rows = shuffle_rows)
    readr::write_rds(grad, gradient_file)
    rm(grad)
  }
  return(gradient_file)
}


estimate_new_beta = function(
  model_name, synced_folder, 
  all_site_names = paste0("site", 1:3),
  tolerance = 1e-8) {
  
  file_list = folder_names(synced_folder)
  gradients_folder = file_list$gradients_folder  
  beta_folder = file_list$beta_folder  
  converged_folder = file_list$converged_folder
  
  
  final_file = file.path(converged_folder,
                         paste0(model_name, ".rds"))
  
  if (file.exists(final_file)) {
    stop("Model already converged, delete iterations to run again")
  }
  
  
  res = get_current_beta(model_name, synced_folder)
  beta = res$beta
  iteration_number = res$iteration_number
  
  out_beta_file = file.path(
    beta_folder,
    paste0(
      model_name, 
      sprintf("-iteration%04.0f", iteration_number),
      ".rds")
  )
  # list_gradient_files = list.files(
  #   gradients_folder, 
  #   pattern = paste0("^", model_name, ".*",
  #                    sprintf("-iteration%04.0f", iteration_number),
  #                    ".rds"),
  #   full.names = TRUE)
  
  all_gradient_files = file.path(
    gradients_folder, 
    paste0(model_name, "-", 
           all_site_names,
           sprintf("-iteration%04.0f", iteration_number),
           ".rds"))
  
  fe = file.exists(all_gradient_files)
  
  # should check if converged
  if (!file.exists(out_beta_file)) {
    if (!all(fe)) {
      print("Waiting for other sites to create gradients")
      print("Missing files:")
      print(all_gradient_files[!fe])
    } else {
      print(paste0(
        "Reading in gradients, iteration ", iteration_number))
      result = aggregate_gradients(
        all_gradient_files, iteration_number)
      gradient = result$gradient
      total_sample_size = result$total_sample_size
      
      if (is.null(beta)) {
        beta = rep(0, length(gradient))
        epsilon = 10
      } else {
        # see glm.control
        epsilon = max(abs(gradient)/(abs(beta) + 0.1))
      }
      if (epsilon < tolerance) {
        print("Model has converged!")
        final_beta_list = list(
          beta = beta,
          num_iterations = iteration_number,
          gradient = gradient,
          tolerance = tolerance,
          epsilon = epsilon,
          total_sample_size = total_sample_size,
          max_gradient = max(abs(gradient)))
        readr::write_rds(final_beta_list, final_file)
        return(final_file)
      }
      beta = beta + gradient
      beta_list = list(
        beta = beta,
        previous_gradient = gradient,
        total_sample_size = total_sample_size,
        iteration_number_next = iteration_number +  1,
        tolerance = tolerance,
        epsilon = epsilon  
      )
      readr::write_rds(beta_list, out_beta_file)
      rm(beta_list)
      return(out_beta_file)
      
    }
  } else {
    if (!all(fe)) {
      print("Waiting for other sites to create gradients")
      print("Missing files:")
      print(all_gradient_files[!fe])
    } 
  }
  
} 



clear_model = function(
  model_name, synced_folder
) {
  
  file_list = folder_names(synced_folder)
  files = sapply(file_list, function(x) {
    list.files(path = x,
               pattern = paste0("^", model_name, ".*.rds"),
               full.names = TRUE)
  })
  file.remove(unlist(files))
}
