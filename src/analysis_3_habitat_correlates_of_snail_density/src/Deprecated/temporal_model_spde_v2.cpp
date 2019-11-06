
#include <TMB.hpp>
#include "include/functions.hpp" //  Functions for indexing

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // ------------------------------------------------------------------------------------------------------------------------- //
  // 0. Settings
  DATA_INTEGER(debug);


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 1. Data (q = quadrat, s = site, p = presence/absence, c = non-zero density, i = individual snail)
  // -- 1.0. Indices
  DATA_INTEGER( n_s ) ; // Number of sites
  DATA_INTEGER( n_t ) ; // Number of sampling periods
  DATA_INTEGER( n_x_q_cont ); // Continuous covariate columns
  DATA_VECTOR( x_q_hat_est ); // Vector of which covariates to impute random effects

  // -- 1.1 Local density data
  DATA_VECTOR( c_q )  ; // Response (count) for each quadrat
  DATA_IVECTOR( q_q ) ; // Unique quadrat identifiet for each quadrat
  DATA_IVECTOR( t_q ) ; // Field mission identifier for each quatrat
  DATA_IVECTOR( s_q ) ; // Site identifier for each quatrat
  DATA_IVECTOR( v_q ) ; // Village identifiet for each quadrat: NOT USED
  DATA_MATRIX( x_q ) ; // Covariate design matrix 
  // DATA_IVECTOR( predTF_q ); // Index of wether to include likelihood in cross validation
  DATA_VECTOR( x_q_hat_est ); // Do we estimate this variable


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 2. Parameters
  // -- 2.0. Design matrix components
  PARAMETER_VECTOR( beta0 );       // 2.1. -- Regression intercepts
  PARAMETER_VECTOR( beta_p );      // 2.1. -- Regression coefficients for probability of presence
  PARAMETER_VECTOR( beta_c );      // 2.1. -- Regression coefficients for local density
  PARAMETER_VECTOR( beta_log );    // 2.3. -- Logistic parameteres for missing categorical data

  // -- 2.1. Random effects components
  PARAMETER_MATRIX( epsilon_mat ); // 2.2. -- Random site effects
  // 0 -- random site effects for density random effects
  // 1 -- random site effects for positive density error
  // 2 -- random site effects for individual infection of parasite 1
  // 3 -- random site effects for individual infection of parasite 2
  PARAMETER_VECTOR( gamma_q );     // 2.3. -- Overdispersion parameters of log-normal Poisson

  // -- 2.2. Non-spatial variance components
  PARAMETER_VECTOR( log_sigmas );  // 2.2. -- Variance components
  // 0 -- log SD for random site effects for presence/absence effects
  // 1 -- log SD for random site effects for positive density error
  // 2 -- log SD for random site effects for individual infection of parasite 1
  // 3 -- log SD for random site effects for individual infection of parasite 2
  // 4 -- log SD for random effects for overdispersion for each quadrat
  // 5 -- log SD for normal/lognormal distribution or p, 1<p<2, for tweedie distribution

    PARAMETER_MATRIX( x_q_cont_missing );// 2.4. -- Missing values of continuous covariates
  // PARAMETER_MATRIX( x_q_hat_est );// 2.5. -- Missing values of continuous covariates


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 3. Model objects
  // -- Model indices
  int n_proc = 2 ; // number of processes to model (presence/absence, pos density, parasite 1, parasite 2)
  int q, s, t, p, proc ; // q = quadrat, s = site, t = field mission, p = parasites, proc = process
  int n_q = c_q.size(); // Sample size for number of quadrats

  // -- Likelihood components
  vector<Type> jnll_comp(15); jnll_comp.setZero();
  // Slot 0 -- presence/absence
  // Slot 1 -- local density
  // Slot 2 -- random site effect for presence/absence
  // Slot 3 -- random site effect for local density
  // Slot 4 -- deprecated
  // Slot 5 -- deprecated
  // Slot 6 -- overdispersion
  // Slot 7 -- deprecated
  // Slot 8 -- deprecated
  // Slot 9 -- probability of categorical variable being 0
  // Slot 10 -- probability of categorical variable being 1
  // Slot 11 -- probability of missing continuous variables
  Type jnll = 0;

  // ------------------------------------------------------------------------------------------------------------------------- //
  // 4. Impute missing coninuous variables
  matrix<Type> x_q_hat = x_q; // New dseign matrix with imputed values

  // 4.1. Impute missing values
  for(int row = 0; row < x_q.rows(); row++){
    for(int col = 0; col < n_x_q_cont; col++){
    if( isNA( x_q_hat(row, col) ) ){
       x_q_hat(row, col) = x_q_cont_missing( row, col ); // Input missing N(0,1) random effect


      if( !isNA( x_q_hat_est( col ) ) ){ // Only estimate the variables of interest
       jnll_comp(11) -= dnorm( x_q_cont_missing( row, col ), Type(0.0), Type(1.0), true);
        }
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
          if( !isNA( x_q_hat_est( col ) ) ){
          // Data likelihood
          if(x_q_hat(row, col)  == 0) jnll_comp(9) -= log( 1 - pr_x( col - n_x_pq_cont ) ); // Probability of covariate being 0
          if(x_q_hat(row, col)  == 1) jnll_comp(10) -= log( pr_x( col - n_x_pq_cont ) );     // Probability of covariate being 1
        }
      }
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 6. Get linear predictors
  vector<Type> linpred_pq( n_q );
  vector<Type> zero_prob_q( n_q );
  vector<Type> linpred_cq( n_q );
  vector<Type> pred_cq( n_q ); // predicted abundance

  for( q = 0; q < n_q; q++){
    linpred_pq( q ) = beta0(0) + (x_q_hat * beta_p)( q ) + epsilon_mat( 0, s_q( q ))  ; // Linear predictor of presence/absence
    linpred_cq( q ) = beta0(1) + (x_q_hat * beta_c)( q ) + epsilon_mat( 1, s_q( q )) + gamma_q( q ) ; // Linear predictor of abundance
  }


  // Make sure the linear predictor is positive
  if(model != 8){
  zero_prob_q = 1 / (1 + exp( -linpred_pq )); // Probability of absence
  linpred_cq = exp(linpred_cq);
  pred_cq = (1 - zero_prob_q ) * linpred_cq; // Predicted catch
}




  // ------------------------------------------------------------------------------------------------------------------------- //
  // 7. Random effects
  for(proc = 0; proc < n_proc; proc++){
    for( s = 0; s < n_s; s++){
      jnll_comp(2 + proc) -= dnorm( epsilon_mat( proc, s ), Type(0.0), exp(log_sigmas( proc )), true); // Random site effects for presence/absence

      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J1. The first nll of random site effect for presence/absence  is " <<  jnll_comp(2) << std::endl;}
      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J2. The first nll of random site effects for density random is " <<  jnll_comp(3) << std::endl;}
      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J4. The first nll of random site effects for infection probability is " <<  jnll_comp(4) << " and " << jnll_comp(5) << std::endl;}
    }
  }


    for(q = 0; q < n_q; q ++){
      jnll_comp(6) -= dnorm( gamma_q( q ), Type(0.0), exp(log_sigmas( 4 )), true); // Overdispersion for each quadrat
      if(( debug == 1 ) & ( q == 0 )) {std::cerr << "J3. The first nll of overdispersion is " <<  jnll_comp(6) <<std::endl;}
    }
  


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 7. Probability of abundance conditional on fixed and random effects values

    // Delta-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dpois( c_q( q ), linpred_cq( q ), true );
    }
  


  


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 9. Reporting
  REPORT( zero_prob_q );
  REPORT( log_sigmas );
  REPORT( linpred_pq );
  REPORT( linpred_cq );
  REPORT( pred_cq );
  REPORT( c_q );

  ADREPORT( beta_p);
  ADREPORT( beta_c);
  ADREPORT( pr_x );


  if( debug == 1 ){std::cerr << "----------------------------------------------" << std::endl;}

  jnll = jnll_comp.sum();
  REPORT( jnll );
  return jnll;


  // ------------------------------------------------------------------------------------------------------------------------- //
}
