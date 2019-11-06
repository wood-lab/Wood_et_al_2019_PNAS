# Function to expand a design matrix and provide all possible combinations of missing categorical variables.

expand_xq <- function( data_list, map ){

  # Which categorical covariates are turned off? - Set to 0 so no expansion
  x_q <- data_list$x_q
  cat_covar <- grep("cat_covar_", colnames(x_q))
  turned_off <- which(is.na(as.numeric(as.character(map$beta_p))) & is.na(as.numeric(as.character(map$beta_c))))
  turned_off_cat_covar <- turned_off[which(turned_off %in% cat_covar)]
  x_q[,turned_off_cat_covar] <- 0

  data_list$n_na_x_q = rowSums(matrix(is.na( x_q[,cat_covar] ), ncol = ncol( x_q[,cat_covar] )))


  # Expand x_q for missing categorical variables
  x_q_array <- array(NA, dim = c( max(2^data_list$n_na_x_q, na.rm = T) , ncol( x_q ) , nrow(x_q)))

  # Find all combinations
  for(i in 1:nrow(x_q)){


    dat_sub <- x_q[i,]
    na_covar <- which( is.na( dat_sub ) )
    na_cat_covar <- na_covar[which(na_covar %in% cat_covar)]

    # if no missing values
    if(length(na_cat_covar) == 0){
      cat_fill = matrix(dat_sub, nrow = 1)
    }

    # If missing values
    if(length(na_cat_covar) > 0){
    possible_combos <- permutations(2, r = length(na_cat_covar), c(0,1), repeats.allowed = TRUE)
    cat_fill <- matrix(rep(dat_sub, r = 2^length(na_cat_covar)), nrow = 2^length(na_cat_covar), ncol = length(dat_sub), byrow = T)

      for(j in 1:length(na_cat_covar)){
        cat_fill[,na_cat_covar[j]] <- possible_combos[,j]
      }
    }
    # Fill array
    x_q_array[1:nrow(cat_fill),,i] <- cat_fill
  }


  data_list$x_q_array <- x_q_array
  return(data_list)

}
