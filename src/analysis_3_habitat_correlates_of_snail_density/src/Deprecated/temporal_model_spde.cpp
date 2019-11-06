
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

  // 0.2. -- Add probability of disease
  DATA_INTEGER( incl_disease );
  // 0 -- Do not include disease model
  // 1 -- Include individual disease model
  DATA_INTEGER(debug);


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 1. Data (q = quadrat, s = site, p = presence/absence, c = non-zero density, i = individual snail)
  // -- 1.0. Indices
  DATA_INTEGER( n_s ) ; // Number of sites
  DATA_INTEGER( n_t ) ; // Number of sampling periods
  DATA_INTEGER( n_p ) ; // Number of parasite species in snail
  DATA_INTEGER( n_x_pq_cont ); // Continuous covariate columns
  DATA_VECTOR( x_pq_hat_est ); // Vector of which covariates to impute random effects

  // -- 1.1 Local density data
  DATA_VECTOR( c_q )  ; // Response (count) for each quadrat
  DATA_IVECTOR( q_q ) ; // Unique quadrat identifiet for each quadrat
  DATA_IVECTOR( t_q ) ; // Field mission identifier for each quatrat
  DATA_IVECTOR( s_q ) ; // Site identifier for each quatrat
  DATA_IVECTOR( v_q ) ; // Village identifiet for each quadrat: NOT USED
  DATA_MATRIX( x_pq ) ; // Covariate design matrix for binomial process
  DATA_MATRIX( x_cq ) ; // Covariate design matrix for abundance process
  // DATA_IVECTOR( predTF_q ); // Index of wether to include likelihood in cross validation

  // -- 1.2. Individual snail data
  DATA_INTEGER( n_i ) ;    // Sample size for number of individual snails
  DATA_MATRIX(y_i);        // Was the snail infected with schisto? 1 = yes, 0 = no, na = if snail was not diagnosable
  DATA_IVECTOR(s_i);       // Site of an individual snail
  DATA_IVECTOR(q_i);       // Unique quadrat of an individual snail relating to row in local density data
  DATA_ARRAY(x_i);         //Covariate design matrix for individual infection (n_i, n_x, n_p)
  // DATA_IVECTOR( predTF_i ); // Index of wether to include likelihood in cross validation

  // ------------------------------------------------------------------------------------------------------------------------- //
  // 2. Parameters
  // -- 2.0. Design matrix components
  PARAMETER_VECTOR( beta0 );       // 2.1. -- Regression intercepts
  PARAMETER_VECTOR( beta_p );      // 2.1. -- Regression coefficients for probability of presence
  PARAMETER_VECTOR( beta_c );      // 2.1. -- Regression coefficients for local density
  PARAMETER_MATRIX( beta_y );      // 2.1. -- Regression coefficients for probability of disease
  PARAMETER_VECTOR( beta_c_hat );  // 2.2 -- Coefficient for conditional dependence of disease probability on density
  PARAMETER_VECTOR( logistic );    // 2.3. -- Logistic parameteres for missing categorical data

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

    PARAMETER_MATRIX( x_pq_cont_missing );// 2.4. -- Missing values of continuous covariates
  // PARAMETER_MATRIX( x_cq_cont_missing );// 2.5. -- Missing values of continuous covariates


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 3. Model objects
  // -- Model indices
  int n_proc = 2 ; // number of processes to model (presence/absence, pos density, parasite 1, parasite 2)
  int q, i, s, t, p, proc ; // q = quadrat, i = individual, s = site, t = field mission, p = parasites, proc = process
  int n_q = c_q.size(); // Sample size for number of quadrats
  int overdispersion = 0; // Index of wether to include overdispersion
  if(( model == 2 ) | ( model == 6 )) { overdispersion = 1;}
  if( incl_disease == 1 ){ n_proc += n_p; }

  // -- Likelihood components
  vector<Type> jnll_comp(15); jnll_comp.setZero();
  // Slot 0 -- presence/absence
  // Slot 1 -- local density
  // Slot 2 -- random site effect for presence/absence
  // Slot 3 -- random site effect for local density
  // Slot 4 -- random site effect for disease probability of parasite 1
  // Slot 5 -- random site effect for disease probability of parasite 2
  // Slot 6 -- overdispersion
  // Slot 7 -- probability of no disease parasite 1
  // Slot 8 -- probability of no disease parasite 2
  // Slot 9 -- probability of disease parasite 1
  // Slot 10 -- probability of disease parasite 2
  // Slot 11 -- probability of missing continuous variables
  Type jnll = 0;

  // ------------------------------------------------------------------------------------------------------------------------- //
  // 4. Impute missing coninuous variables
  matrix<Type> x_pq_hat = x_pq; // New dseign matrix with imputed values
  matrix<Type> x_cq_hat = x_cq; // New dseign matrix with imputed values

  // 4.1. Impute missing values
  for(int row = 0; row < x_pq.rows(); row++){
    for(int col = 0; col < n_x_pq_cont; col++){
    if( isNA( x_pq_hat(row, col) ) ){
       x_pq_hat(row, col) = x_pq_cont_missing( row, col ); // Input missing N(0,1) random effect
       x_cq_hat(row, col) = x_pq_cont_missing( row, col ); //

      if( !isNA( x_pq_hat_est( col ) ) ){ // Only estimate the variables of interest
       jnll_comp(11) -= dnorm( x_pq_cont_missing( row, col ), Type(0.0), Type(1.0), true);
       // jnll_comp(11) -= dnorm( x_cq_cont_missing( row, col ), Type(0.0), Type(1.0), true)
        }
      }
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 5. Impute missing categorical variables

  // 5.1. Impute missing values
  for(int row = 0; row < x_pq.rows(); row++){
    for(int col = n_x_pq_cont; col < x_pq.cols(); col++){
      
    if( isNA( x_pq_hat(row, col) ) ){
       x_pq_hat(row, col) = x_pq_cont_missing( row, col ); // Input missing N(0,1) random effect
       x_cq_hat(row, col) = x_pq_cont_missing( row, col ); //

      if( !isNA( x_pq_hat_est( col ) ) ){ // Only estimate the variables of interest
       jnll_comp(11) -= dnorm( x_pq_cont_missing( row, col ), Type(0.0), Type(1.0), true);
       // jnll_comp(11) -= dnorm( x_cq_cont_missing( row, col ), Type(0.0), Type(1.0), true)
        }
      }
    }
  }

  // ------------------------------------------------------------------------------------------------------------------------- //
  // 5. Get linear predictors
  vector<Type> linpred_pq( n_q );
  vector<Type> zero_prob_q( n_q );
  vector<Type> linpred_cq( n_q );
  vector<Type> pred_cq( n_q ); // predicted abundance

  matrix<Type> linpred_yi( n_i, n_p );
  matrix<Type> zero_prob_yi( n_i, n_p );

  for( q = 0; q < n_q; q++){
    linpred_pq( q ) = beta0(0) + (x_pq_hat * beta_p)( q ) + epsilon_mat( 0, s_q( q )) ; // Linear predictor of presence/absence
    linpred_cq( q ) = beta0(1) + (x_cq_hat * beta_c)( q ) + epsilon_mat( 1, s_q( q )) ; // Linear predictor of abundance
  }

  if( overdispersion == 1 ){
    for(q = 0; q < n_q; q ++){
      linpred_cq( q ) = linpred_cq( q ) + gamma_q( q );
    }
  }

  // Make sure the linear predictor is positive
  if(model != 8){
  zero_prob_q = 1 / (1 + exp( -linpred_pq )); // Probability of absence
  linpred_cq = exp(linpred_cq);
  pred_cq = (1 - zero_prob_q ) * linpred_cq; // Predicted catch
}

if(model == 8){ // Tweedie parameters
  zero_prob_q = exp( linpred_pq ); 
  linpred_cq = exp( linpred_cq );
}


  if(incl_disease == 1){
    for(i = 0; i < n_i; i++){
      for( p = 0; p < n_p; p ++){
        vector<Type> beta_temp = beta_y.row(p); // Temporary storage for multiplication
        linpred_yi( i, p ) = beta0(2 + p) + (matrix_from_array(x_i, p) * beta_temp)(i) + beta_c_hat(p) * linpred_cq( q_i( i ) ) + epsilon_mat( 2 + p, s_i( i ) ) ; // Linear predictor for infection
        zero_prob_yi( i, p  ) = 1/ (1 + exp( - linpred_yi( i, p ) ) ); // Probility of infection
      }
    }
  }

  if( debug == 1 ) {std::cerr << "E1. The first predicted non-zero catch is " <<  linpred_cq(0) <<std::endl;}
  if( debug == 1 ) {std::cerr << "E2. The first zero probability is " <<  zero_prob_q(0) <<std::endl;}
  if( debug == 1 ) {std::cerr << "E3. The first predicted catch " <<  pred_cq(0) <<std::endl;}
  if(( debug == 1 ) & ( incl_disease == 1 )) {std::cerr << "E4. The first infection probability is " <<  zero_prob_yi(0 , 0) <<std::endl;}


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 6. Random effects
  for(proc = 0; proc < n_proc; proc++){
    for( s = 0; s < n_s; s++){
      jnll_comp(2 + proc) -= dnorm( epsilon_mat( proc, s ), Type(0.0), exp(log_sigmas( proc )), true); // Random site effects for presence/absence

      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J1. The first nll of random site effect for presence/absence  is " <<  jnll_comp(2) << std::endl;}
      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J2. The first nll of random site effects for density random is " <<  jnll_comp(3) << std::endl;}
      if(( debug == 1 ) & ( s == 0 )) {std::cerr << "J4. The first nll of random site effects for infection probability is " <<  jnll_comp(4) << " and " << jnll_comp(5) << std::endl;}
    }
  }

  if( overdispersion == 1 ){
    for(q = 0; q < n_q; q ++){
      jnll_comp(6) -= dnorm( gamma_q( q ), Type(0.0), exp(log_sigmas( 4 )), true); // Overdispersion for each quadrat
      if(( debug == 1 ) & ( q == 0 )) {std::cerr << "J3. The first nll of overdispersion is " <<  jnll_comp(6) <<std::endl;}
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 7. Probability of abundance conditional on fixed and random effects values

  // -- Model 1: Delta-lognormal
  if(model == 0){
    // Delta-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dlognorm( c_q(q), log( exp( linpred_cq( q ) ) - pow( exp( log_sigmas( 5 ) ) , 2)/2), exp( log_sigmas( 5 ) ), true );
    }
  }


  // -- Model 1: Delta-poisson and Model 2: delta-log-normal Poisson
  if((model == 1)| (model == 2)){
    // Delta-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dpois( c_q( q ), linpred_cq( q ), true );
    }
  }


  // -- delta-negative-binomial
  if(model == 3){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }


  // -- ZI-lognormal
  if(model == 4){
    // Delta-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) ) + dlognorm( Type(0.0), log( exp( linpred_cq( q ) ) - pow( exp( log_sigmas( 5 ) ) , 2)/2), exp( log_sigmas( 5 ) ), true );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dlognorm( c_q(q), log( exp( linpred_cq( q ) ) - pow( exp( log_sigmas( 5 ) ) , 2)/2), exp( log_sigmas( 5 ) ), true );
    }
  }


  // -- Model 5: ZI poisson and Model 6: ZI-log-normal Poisson
  if((model == 5)|(model == 6)){
    // ZIP-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) ) + dpois( Type(0.0), linpred_cq( q ), true );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dpois( c_q( q ), linpred_cq( q ), true );

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
    for( q = 0; q < n_q; q++){
      jnll_comp(0) -= dtweedie(c_q( q ), linpred_cq( q ), zero_prob_q( q ), invlogit( log_sigmas( 5 ) ) + Type(1), true);
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 8. Probability of infection conditional on abundance, fixed, and random effects values
  // Individual model
  if(incl_disease == 1){
    for( i = 0; i < n_i; i++){
      for( p = 0; p < n_p; p ++){
        if(y_i( i, p ) == 0) jnll_comp(7 + p) -= log( 1 - zero_prob_yi( i, p ) );
        if(y_i( i, p ) == 1) jnll_comp(9 + p) -= log( zero_prob_yi( i, p ) );

        if(( debug == 1 ) & ( i == 0 )) {std::cerr << "J6. The first nll of no infection probability is " <<  jnll_comp(5) << " and " << jnll_comp(6) << std::endl;}
        if(( debug == 1 ) & ( i == 0 )) {std::cerr << "J7. The first nll of infection probability is " <<  jnll_comp(7) << " and " << jnll_comp(8) << std::endl;}
      }
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 9. Reporting
  REPORT( zero_prob_q );
  REPORT( zero_prob_yi );
  REPORT( log_sigmas );
  REPORT( linpred_pq );
  REPORT( linpred_cq );
  REPORT( pred_cq );
  REPORT( c_q );

  ADREPORT( beta_p);
  ADREPORT( beta_c);
  ADREPORT( beta_y);


  if( debug == 1 ){std::cerr << "----------------------------------------------" << std::endl;}

  jnll = jnll_comp.sum();
  REPORT( jnll );
  return jnll;


  // ------------------------------------------------------------------------------------------------------------------------- //
}
