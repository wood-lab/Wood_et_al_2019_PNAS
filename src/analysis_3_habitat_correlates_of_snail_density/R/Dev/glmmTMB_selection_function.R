# These functions iteratively update GAM models following a forward selection approach
# The first function allows the update of a single model with included character term
# The second function iteratively updates a series of models

# The user inputs a base model and a list of character terms present in the data set in the order they would like them added.
# The use can also input the number of updated models they want to create by specifying n_loops


library(glmmTMB)


model.list = list()

# Model Update Function
fun_add = function(mod, x, first = T){
  m = mod
  if(first){
    form = as.formula(paste(".~ 1 + ",x))
  }
  else{
    form = as.formula(paste(".~. + ",x))
  }
  m2 = update(m, form)
  return(m2)
}

# Model term removal Function
fun_remove = function(mod, x){
  m = mod
  
  form = as.formula(paste(".~. - ",x))
  
  m2 = update(m, form)
  return(m2)
}

# Model Selection
forward_model_select = function(base_mod, terms){
  # Check terms
  if(missing(base_mod)){stop("Base model is missing")}
  if(missing(terms)){stop("Update model terms are missing")}
  
  # If n_loops is not specified, set it to the number of terms.
  if(missing(n_loops)){
    n_loops = length(terms) 
  }
  
  # Summary stats and model empty objects
  ind = c()
  remove_vec <- c()
  model.list = list()
  summ_stats = data.frame(matrix(NA, ncol = 6 , nrow = 1))
  colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params","Param")
  
  
  # Forward model selection
  for( i  in 1:(n_loops)){
    
    # Assign base model
    if (i ==1) {
      model.list[[i]] = base_mod
      best_mod <- base_mod
    }else{
      # as.formula will update a base model to '~+' if we don't correct for it
      if(i == 2){
        check = TRUE
      } else {
        check = FALSE
      }
      # Print status
      print(paste("Optimization of model ",i,sep = ""))
      
      # Update model
      model.list[[i]] = fun_add(best_mod, terms[i], first =  check)
      
      # Add model to remove.vec if it does not improve AIC
      if(AIC(model.list[[i]]) <= AIC(best_mod)){
        best_mod <- model.list[[i]]
        print(paste(terms[i],"improves model fit", sep = ""))
      }
      
    }
    # Summary stats
    colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params")
    new_vals <- c(as.character(formula(model.list[[i]])), AIC(model.list[[i]]), BIC(model.list[[i]]), as.numeric(logLik(model.list[[i]])), as.numeric((AIC(model.list[[i]]) - (-logLik(model.list[[i]]) * 2))/2))
    summ_stats <- rbind(summ_stats, new_vals)
  }
  
  
  
  
  # Backward model selection
  for( i  in 1:(n_loops)){
    
      # Print status
      print(paste("Optimization of model ",i,sep = ""))
      
      # Update model
      model.list[[i]] = fun_add(best_mod, terms[i], first =  check)
      
      # Add model to remove.vec if it does not improve AIC
      if(AIC(model.list[[i]]) <= AIC(best_mod)){
        best_mod <- model.list[[i]]
        print(paste(terms[i],"improves model fit", sep = ""))
      }
      
    }
    # Summary stats
    colnames(summ_stats) = c("Model", "AIC", "BIC","Log_Lik","N_Params")
    new_vals <- c(as.character(formula(model.list[[i]])), AIC(model.list[[i]]), BIC(model.list[[i]]), as.numeric(logLik(model.list[[i]])), as.numeric((AIC(model.list[[i]]) - (-logLik(model.list[[i]]) * 2))/2))
    summ_stats <- rbind(summ_stats, new_vals)
  }
  
  
  
  
  
  # Return results
  results_list = list(summ_stats = summ_stats, model.list = model.list)
  
  return(results_list)
  
}