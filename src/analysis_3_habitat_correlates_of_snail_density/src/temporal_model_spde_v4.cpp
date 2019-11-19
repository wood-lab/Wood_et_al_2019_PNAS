
#include <TMB.hpp>
#include "include/functions.hpp" //  Functions for indexing

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // ------------------------------------------------------------------------------------------------------------------------- //
  // 0. Settings
  DATA_INTEGER( model );
  //  0 -- delta-lognormal
  //  1 -- delta-poisson
  //  2 -- delta-log-normal Poisson
  //  3 -- delta-negative-binomial
  //  4 -- ZI-lognormal
  //  5 -- ZI-poisson
  //  6 -- ZI-log-normal Poisson
  //  7 -- ZI-negative-binomial
  //  8 -- Tweedie distribution
  DATA_INTEGER(debug);


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 1. Data (q = quadrat, s = site, p = presence/absence, c = non-zero density, i = individual snail)
  // -- 1.0. Indices
  DATA_INTEGER( n_s ) ; // Number of sites
  // DATA_INTEGER( n_t ) ; // Number of sampling periods
  DATA_INTEGER( n_v ) ; // Number of unique villages
  // DATA_INTEGER( n_p ) ; // Number of parasite species in snail
  DATA_INTEGER( n_x_q_cont ); // Continuous covariate columns
  DATA_VECTOR( x_q_hat_est ); // Vector of which covariates to impute random effects
  DATA_VECTOR( n_na_x_q );    // Number of missing categorical variables for each observation

  // -- 1.1 Local density data
  DATA_VECTOR( c_q )  ; // Response (count) for each quadrat
  // DATA_IVECTOR( q_q ) ; // Unique quadrat identifiet for each quadrat
  // DATA_IVECTOR( t_q ) ; // Field mission identifier for each quatrat
  DATA_IVECTOR( s_q ) ; // Site identifier for each quatrat
  DATA_IVECTOR( v_q ) ; // Village identifiet for each quadrat: NOT USED
  DATA_MATRIX( x_q ) ; // Covariate design matrix
  DATA_ARRAY( x_q_array ); // Covariate design matrix expanded to include all possible combinations of missing variables
  DATA_IVECTOR( fit_ll ); // Index of wether to include likelihood in cross validation


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 2. Parameters
  // -- 2.0. Design matrix components
  PARAMETER_VECTOR( beta0 );       // 2.1. -- Regression intercepts
  PARAMETER_VECTOR( beta_p );      // 2.1. -- Regression coefficients for probability of presence
  PARAMETER_VECTOR( beta_c );      // 2.1. -- Regression coefficients for local density
  PARAMETER_VECTOR( beta_log );    // 2.3. -- Logistic parameteres for missing categorical data

  // -- 2.1. Random effects components
  PARAMETER_MATRIX( epsilon_mat ); // 2.1. -- Random site effects
  // 0 -- random site effects for density random effects
  // 1 -- random site effects for positive density error
  PARAMETER_MATRIX( omega_mat ); // 2.2. -- Random village effects
  // 0 -- random village effects for density random effects
  // 1 -- random village effects for positive density error
  PARAMETER_VECTOR( gamma_q );     // 2.3. -- Overdispersion parameters of log-normal Poisson

  // -- 2.2. Non-spatial variance components
  PARAMETER_VECTOR( log_sigmas );  // 2.2. -- Variance components
  // 0 -- log SD for random site effects for presence/absence effects
  // 1 -- log SD for random site effects for positive density error
  // 2 -- log SD for random village effects for presence/absence effects
  // 3 -- log SD for random village effects for positive density effects
  // 4 -- log SD for random effects for overdispersion for each quadrat
  // 5 -- log SD for normal/lognormal distribution or p, 1<p<2, for tweedie distribution

  PARAMETER_VECTOR( x_q_cont_missing );// 2.4. -- Missing values of continuous covariates


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 3. Model objects
  // -- Model indices
  int n_proc = 2 ; // number of processes to model (presence/absence, pos density, parasite 1, parasite 2)
  int q, s, proc ; // q = quadrat, i = individual, s = site, t = field mission, p = parasites, proc = process
  int n_q = c_q.size(); // Sample size for number of quadrats
  int overdispersion = 0; // Index of wether to include overdispersion
  if(( model == 2 ) | ( model == 6 )) { overdispersion = 1;}
  int n_cat_vars = x_q.cols() - n_x_q_cont; // Number of categorical predictors
  int n_like_q ; // Number of likelihood statements needed to expand for row q
  // REPORT(n_cat_vars);
  // REPORT( n_q );

  // -- Likelihood components
  matrix<Type> jnll_comp(15, n_q); jnll_comp.setZero();
  // Slot 0 -- presence/absence
  // Slot 1 -- local density
  // Slot 2 -- random site effect for presence/absence
  // Slot 3 -- random site effect for local density
  // Slot 4 -- random village effect for presence/absence
  // Slot 5 -- random village effect for local density
  // Slot 6 -- overdispersion
  // Slot 7 -- empty
  // Slot 8 -- empty
  // Slot 9 -- Probability of covariate being 0
  // Slot 10 -- Probability of covariate being 1
  // Slot 11 -- random effect of missing continuous variables
  Type jnll = 0;

  // ------------------------------------------------------------------------------------------------------------------------- //
  // 4. Impute missing coninuous variables as random effects
  matrix<Type> x_q_hat = x_q; // New dseign matrix with imputed values

  // 4.1. Impute missing values
  int q_ind = 0; // Index for missing continuous variable input
  for(int row = 0; row < x_q.rows(); row++){
    for(int col = 0; col < n_x_q_cont; col++){

      if( isNA( x_q_hat(row, col) ) ){
        x_q_hat(row, col) = x_q_cont_missing( q_ind ); // Input missing N(0,1) random effect

        // Fill expanded dataset with random effects for missing continuous variables
        for( n_like_q = 0; n_like_q < pow(2, n_na_x_q( row ) ); n_like_q++ ){
          x_q_array(n_like_q, col, row ) = x_q_hat(row, col);
        }

        // Missing data likelihood
        if( !isNA( x_q_hat_est( col ) ) ){ // Only estimate the variables of interest

          // Fit or predict?
          if(fit_ll(row) == 1){
            jnll_comp(11, row) -= dnorm( x_q_cont_missing( q_ind ), Type(0.0), Type(1.0), true);
          }

          SIMULATE{
            x_q_cont_missing( q_ind ) = rnorm(Type(0.0), Type(1.0));
          }
        }
        q_ind += 1;
      }
    }
  }



  // ------------------------------------------------------------------------------------------------------------------------- //
  // 5. Probablilty of categorical variables
  vector<Type> pr_x = 1/(1 + exp(-beta_log)); // Probability of factor being 1
  for(int row = 0; row < x_q.rows(); row++){
    for(int col = n_x_q_cont; col < x_q.cols(); col++){

      // Fit available data
      if( !isNA( x_q_hat(row, col)  ) ){
        if(fit_ll(row) == 1){
          //if( !isNA( x_q_hat_est( col ) ) ){
          // Data likelihood
          if(x_q_hat(row, col)  == 0) jnll_comp(9, row) -= log( 1 - pr_x( col ) );  // Probability of covariate being 0
          if(x_q_hat(row, col)  == 1) jnll_comp(10, row) -= log( pr_x( col ) );     // Probability of covariate being 1
          //}
        }
      }
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 5. Get linear predictors
  array<Type> linpred_pq_array( pow(2, n_cat_vars), 2, n_q ); linpred_pq_array.setZero();// array of linear predictor - probability of predictor for each obs
  array<Type> linpred_cq_array( pow(2, n_cat_vars), 2, n_q ); linpred_cq_array.setZero(); // array of linear predictor - probability of predictor for each obs
  vector<Type> zero_prob_q( n_q ); zero_prob_q.setZero(); // Average zero probability for each observations
  vector<Type> pos_count_q( n_q ); pos_count_q.setZero(); // Average positive density for each observation
  vector<Type> pred_cq( n_q ); // predicted abundance

  matrix<Type> x_q_temp = matrix_from_array( x_q_array, 0); // Temporary vector of design matrix row
  x_q_temp.setZero();

  for( q = 0; q < n_q; q++){
    x_q_temp.setZero();
    x_q_temp = matrix_from_array( x_q_array, q); // Extract the qth row of design matrix

    // Set up random effects model and intercept
    for( n_like_q = 0; n_like_q < pow(2, n_cat_vars ); n_like_q++ ){
      linpred_pq_array( n_like_q , 0 , q ) = beta0(0) + epsilon_mat( 0, s_q( q )) + omega_mat( 0, v_q( q )); // Linear predictor of presence/absence
      linpred_cq_array( n_like_q , 0 , q ) = beta0(1) + epsilon_mat( 1, s_q( q )) + omega_mat( 1, v_q( q )); // Linear predictor of abundance

      if( overdispersion == 1 ){
        linpred_cq_array( n_like_q , 0 , q ) += gamma_q( q );
      }
    }

    //
    for( n_like_q = 0; n_like_q < pow(2, n_na_x_q( q ) ); n_like_q++ ){
      // Get linear predictor
      linpred_pq_array( n_like_q , 0 , q ) += ( x_q_temp * beta_p)( n_like_q ) ; // Add linear predictor of presence/absence
      linpred_cq_array( n_like_q , 0 , q ) += ( x_q_temp * beta_c)( n_like_q ) ; // Add linear predictor of abundance

      // If there is NAs - find probability of categorical var = 1
      if( n_na_x_q(q) > 0){ // If NAs are present
        for(int col = n_x_q_cont; col < x_q.cols(); col++){// Loop through categorical columns
          if( isNA( x_q_hat(q, col)  ) ){ // If NA
            if( !isNA( x_q_hat_est( col ) ) ){ // If estimated

              // Add probability of being 1
              if( x_q_temp( n_like_q , col) == 1 ){
                linpred_pq_array( n_like_q , 1 , q ) += log( pr_x( col ) ) ;
                linpred_cq_array( n_like_q , 1 , q ) += log( pr_x( col ) ) ;
              }

              // Add probability of being 0
              if( x_q_temp( n_like_q , col) == 0 ){
                linpred_pq_array( n_like_q , 1 , q ) +=  log( 1 - pr_x( col ) ) ;
                linpred_cq_array( n_like_q , 1 , q ) +=  log( 1 - pr_x( col ) ) ;
              }
            }
          }
        }
      }

      // If NAs are not present
      if( n_na_x_q(q) == 0){
        linpred_pq_array( 0 , 1 , q ) = 0; // Log probability multiplier - since there is no missing data p(x) = 1
        linpred_cq_array( 0 , 1 , q ) = 0; // Log probability multiplier - since there is no missing data p(x) = 1
      }

      // Get average probability of zero
      if(model != 8){
        zero_prob_q( q ) += invlogit(linpred_pq_array(n_like_q, 0, q)) * exp(linpred_pq_array( n_like_q , 1 , q ));
      }
      if(model == 8){ // Tweedie parameters
        zero_prob_q( q ) += exp(linpred_pq_array(n_like_q, 0, q)) * exp(linpred_pq_array( n_like_q , 1 , q ));
      }
      pos_count_q( q ) += exp(linpred_cq_array(n_like_q, 0, q)) * exp(linpred_cq_array( n_like_q , 1 , q ));
    }
  }


  // Expected count across observations
  if(model != 8){
    pred_cq = (1 - zero_prob_q ) * pos_count_q; // Predicted catch
  }


  if( debug == 1 ) {std::cerr << "E1. The first predicted non-zero catch is " <<  linpred_cq_array(0) <<std::endl;}
  if( debug == 1 ) {std::cerr << "E2. The first zero probability is " <<  zero_prob_q(0) <<std::endl;}
  if( debug == 1 ) {std::cerr << "E3. The first predicted catch " <<  pred_cq(0) <<std::endl;}


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 6. Random effects
  for(proc = 0; proc < n_proc; proc++){

    // Random village effects
    for(int v = 0; v < n_v; v++){
      jnll_comp(4 + proc, v) -= dnorm( omega_mat( proc, v ), Type(0.0), exp(log_sigmas( proc + 2 )), true); // Random site effects for presence/absence

      SIMULATE{
        omega_mat( proc, v ) = rnorm( Type(0.0), exp(log_sigmas( proc + 2 ))); // Simulate random effects
      }
    }

    // Random site effects
    for( s = 0; s < n_s; s++){
      jnll_comp(2 + proc, s) -= dnorm( epsilon_mat( proc, s ), Type(0.0), exp(log_sigmas( proc )), true); // Random site effects for presence/absence

      SIMULATE{
        epsilon_mat( proc, s ) = rnorm( Type(0.0), exp(log_sigmas( proc ))); // Simulate random effects
      }

      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J1. The first nll of random site effect for presence/absence  is " <<  jnll_comp(2) << std::endl;}
      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J2. The first nll of random site effects for density random is " <<  jnll_comp(3) << std::endl;}
      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J4. The first nll of random site effects for infection probability is " <<  jnll_comp(4) << " and " << jnll_comp(5) << std::endl;}
    }
  }

  if( overdispersion == 1 ){
    for(q = 0; q < n_q; q ++){
      if(fit_ll(q) == 1){
        jnll_comp(6, q) -= dnorm( gamma_q( q ), Type(0.0), exp(log_sigmas( 4 )), true); // Overdispersion for each quadrat
      }
      SIMULATE{
        gamma_q( q ) = rnorm( Type(0.0), exp(log_sigmas( 4 ))); // Simulate overdispersion
      }
      if(( debug == 1 ) & ( q == 0 )) {std::cerr << "J3. The first nll of overdispersion is " <<  jnll_comp(6) <<std::endl;}
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 7. Probability of abundance conditional on fixed and random effects values
  vector<Type> jnll_tmp(4); jnll_tmp.setZero(); // Temporary likelihood to loop through

  // -- Model 0: Delta-lognormal
  if(model == 0){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;

    for( q = 0; q < n_q; q++){
      jnll_tmp.setZero(); // Reset temporary bit
      for( n_like_q = 0; n_like_q < pow(2, n_na_x_q( q ) ); n_like_q++ ){
        jnll_tmp(0) += (invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ) * exp( linpred_pq_array( n_like_q , 1 , q ) )  ;
        jnll_tmp(2) += ( ( 1 - invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ) * dlognorm( c_q( q ), log( exp( linpred_cq_array( n_like_q , 0 , q ) ) - pow( exp( log_sigmas( 5 ) ) , 2)/2), exp( log_sigmas( 5 ) ), false ) * exp(linpred_cq_array( n_like_q , 1 , q ) ) );
      }
      if(fit_ll(q) == 1){
        if(c_q( q ) == 0) if(jnll_tmp(0) > 0) jnll_comp(0, q) -= log(jnll_tmp(0));
        if(c_q( q ) != 0) if(jnll_tmp(2) > 0) jnll_comp(1, q) -= log(jnll_tmp(2));
      }
    }
  }


  // -- Model 1: Delta-poisson and Model 2: delta-log-normal Poisson
  if((model == 1)| (model == 2)){
    for( q = 0; q < n_q; q++){

      jnll_tmp.setZero(); // Reset temporary bit
      for( n_like_q = 0; n_like_q < pow(2, n_na_x_q( q ) ); n_like_q++ ){
        jnll_tmp(0) += (invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ) * exp( linpred_pq_array( n_like_q , 1 , q ) )  ;
        jnll_tmp(2) += ( ( 1 - invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ) * dpois( c_q( q ), exp( linpred_cq_array( n_like_q , 0 , q ) ), false ) * exp(linpred_cq_array( n_like_q , 1 , q ) ) );

        SIMULATE{
          c_q( q ) =  rpois( exp( linpred_cq_array( n_like_q , 0 , q ) ) ) ; // Poisson process
          c_q( q ) *= rbinom(Type(1), Type(1) - invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ); // Binomial process
          c_q( q ) *= exp( linpred_cq_array( n_like_q , 1 , q ) ) ; // Multiply by probability of occuring
        }
      }
      if(fit_ll(q) == 1){
        if(c_q( q ) == 0) if(jnll_tmp(0) > 0) jnll_comp(0, q) -= log(jnll_tmp(0));
        if(c_q( q ) != 0) if(jnll_tmp(2) > 0) jnll_comp(1, q) -= log(jnll_tmp(2));
      }
    }
  }
  // REPORT( jnll_tmp );
  // REPORT( linpred_cq_array );
  // REPORT( linpred_pq_array );

  // -- delta-negative-binomial
  if(model == 3){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }


  // -- ZI-lognormal
  if(model == 4){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }


  // -- Model 5: ZI poisson and Model 6: ZI-log-normal Poisson
  if((model == 5)|(model == 6)){
    // ZIP-model
    for( q = 0; q < n_q; q++){

      jnll_tmp.setZero(); // Reset temporary bit
      for( n_like_q = 0; n_like_q < pow( n_na_x_q(q), 2); n_like_q++ ){
        jnll_tmp(0) += invlogit( linpred_pq_array( n_like_q , 0 , q ) )  * dpois( Type(0.0), exp( linpred_cq_array( n_like_q , 0 , q ) ), false ) * exp( linpred_pq_array( n_like_q , 1 , q ) )  ;
        jnll_tmp(2) += ( ( 1 - invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ) * dpois( c_q( q ), exp( linpred_cq_array( n_like_q , 0 , q ) ), false ) * exp(linpred_cq_array( n_like_q , 1 , q ) ) );

        SIMULATE{
          c_q( q ) =  rpois( exp( linpred_cq_array( n_like_q , 0 , q ) ) ) ; // Poisson process
          c_q( q ) *= rbinom(Type(1), Type(1) - invlogit( linpred_pq_array( n_like_q , 0 , q ) ) ); // Binomial process
          c_q( q ) *= exp( linpred_cq_array( n_like_q , 1 , q ) ) ; // Multiply by probability of occuring
        }
      }

      if(fit_ll(q) == 1){
        if(c_q( q ) == 0) if(jnll_tmp(0) > 0) jnll_comp(0, q) -= log(jnll_tmp(0));
        if(c_q( q ) != 0) if(jnll_tmp(2) > 0) jnll_comp(1, q) -= log(jnll_tmp(2));
      }

      if(( debug == 1 ) & ( q == 0 )) {std::cerr << "J5. The first nll of 0 abundance is " <<  jnll_comp(0) << " and positive abundance " << jnll_comp(1) << std::endl;}
    }
  }

  // -- ZI-negative-binomial
  if(model == 7){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }


  // -- Model 8 Tweedie
  if(model == 8){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }



  // ------------------------------------------------------------------------------------------------------------------------- //
  // 9. Reporting
  // REPORT( zero_prob_q );
  // REPORT( log_sigmas );
  vector<Type> sigmas = exp(log_sigmas);
  REPORT(sigmas);

  REPORT( pred_cq );
  // REPORT( epsilon_mat );
  // REPORT( gamma_q );
  // REPORT( linpred_pq_array );
  // REPORT( linpred_cq_array );

  ADREPORT( beta0 );
  ADREPORT( beta_p);
  ADREPORT( beta_c);
  // ADREPORT( sigmas);
  // ADREPORT( epsilon_mat );


  if( debug == 1 ){std::cerr << "----------------------------------------------" << std::endl;}

  jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  return jnll;


  // ------------------------------------------------------------------------------------------------------------------------- //
}
