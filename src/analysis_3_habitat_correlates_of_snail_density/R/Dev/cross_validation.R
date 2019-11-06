# Function to run a k-fold cross validation on a provided data-set and specified delta-regression model. Returns a variety of test statistics inlcuding the mean log-predictive score, root mean squared error (RMSE), median of relative error (MRE), and median of the absolute relative error (MARE).

# User specifies
# -- the number of partitions (k) the cross validation will use, defaults to 10.
# -- the model used:
# --   0 -- Zero-inflated lognormal
# --   1 -- Zero-inflated poisson
# --   2 -- Zero-inflated negative-binomial
# --   3 -- Tweedie distribution
# -- The data used

cross_val <- function(K = 10, model = 1, data){

  library(TMB)

  # Step 0 -- make and compile template file
  setwd("src")
  compile( "spatio_temporal_model.cpp" )
  dyn.load( dynlib("spatio_temporal_model") )
  setwd("../")

  # Step 1 -- divide into partitions
  Partition_q = sample( x=1:K, size=length(data[["c_q"]]), replace=TRUE )
  Partition_i = sample( x=1:K, size=length(data[["y_i"]]), replace=TRUE )
  PredNLL_k = rep(NA, K)

  # Step 2 --Loop through partitions
  for(k in 1:K){
    Params = list("b_p" = rep(0, ncol(data[["X_pi"]])),
                  "b_c" = rep(0, ncol(data[["X_ci"]])),
                  "epsilon_ps" = rep(0, data[["n_s"]]),
                  "epsilon_cs" = rep(0, data[["n_s"]]))

    if(model %in% c(2, 6)){
      Params[["log_sigmas"]] = rep(0, 4)
    } else{
      Params[["log_sigmas"]] = rep(0, 2)
    }
    random = c("epsilon_ps", "epsilon_cs")

    data[["predTF_q"]] = ifelse(Partition_q==k,1,0)
    data[["predTF_i"]] = ifelse(Partition_i==k,1,0)
    data[["model"]] = model

    Obj = MakeADFun( data=Data, parameters=Params, DLL="spatio_temporal_model", random = random)

    # Optimize
    Opt = TMBhelper::Optimize( obj=Obj, newtonsteps=1, getsd=FALSE )
    SD = sdreport( Obj ) # standard errors
    final_gradient = Obj$gr( Opt$par )

    # Check convergence
    if( any(abs(final_gradient)>0.001) | SD$pdHess==FALSE ) stop("Not converged")

    # Report stuff
    Report = Obj$report()
    PredNLL_k[k] = Report$pred_jnll

  }

  # log-Predictive probability per datum
  mean( PredNLL_k / table(Partition_i) )

}





