
#include <TMB.hpp>
#include <include/functions.hpp> //  Functions for indexing

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  // ------------------------------------------------------------------------------------------------------------------------- //
  // 0. Setting
  DATA_INTEGER( model );
  //  0 -- delta-lognormal
  //  1 -- delta-poisson
  //  2 -- delta-log-normal Poisson
  //  3 -- delta-negative-binomial
  //  4 -- ZI-lognormal
  //  5 -- ZI-poisson
  //  6 -- ZI-log-normal Poisson
  //  7 -- ZI-negative-binomial

  DATA_INTEGER( spatial_model ); 
  // 0 -- SPDE

  // 0.2. -- Add probability of disease
  DATA_INTEGER( incl_disease );
  // 0 -- Do not include disease model
  // 1 -- Include individual disease model
  DATA_INTEGER(debug);


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 1. Data (q = quadrat, s = site, p = presence/absence, c = non-zero density, i = individual snail)
  // -- 1.1 Local density data
  DATA_INTEGER( n_s ) ; // Number of sites
  DATA_VECTOR( c_q )  ; // Response (count) for each quadrat
  DATA_IVECTOR( q_q ) ; // Unique quadrat identifiet for each quadrat
  DATA_IVECTOR( s_q ) ; // Site identifier for each quatrat
  DATA_IVECTOR( v_q ) ; // Village identifiet for each quadrat: NOT USED
  DATA_MATRIX( x_pq ) ; // Covariate design matrix for binomial process
  DATA_MATRIX( x_cq ) ; // Covariate design matrix for abundance process
  // DATA_IVECTOR( predTF_q ); // Index of wether to include likelihood in cross validation

  // 1.2. -- Individual snail data
  DATA_INTEGER( n_i ) ; // Sample size for number of individual snails
  DATA_INTEGER( n_p ) ; // Number of parasites in snails
  DATA_MATRIX(y_i); // Was the snail infected with schisto? 1 = yes, 0 = no, na = if snail was not diagnosable
  DATA_IVECTOR(s_i); // Site of an individual snail
  DATA_IVECTOR(q_i); // Unique quadrat of an individual snail relating to row in local density data
  DATA_ARRAY(x_i); //Covariate design matrix for individual infection (n_i, n_x, n_p)
  // DATA_IVECTOR( predTF_i ); // Index of wether to include likelihood in cross validation

  // 1.3. -- Spatial model components (identified via q_q or q_i)


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 2. Parameters
  // Used in matrices
  // Row 0 -- presence/absence parameters
  // Row 1 -- Abundance parameters
  // Row 3 -- Individual disease parameters
  PARAMETER_VECTOR( beta0 );       // 2.1. -- Regression intercepts
  PARAMETER_VECTOR( beta_p );    // 2.1. -- Regression coefficients for probability of presence
  PARAMETER_VECTOR( beta_c );    // 2.1. -- Regression coefficients for local density
  PARAMETER_MATRIX( beta_y );    // 2.1. -- Regression coefficients for probability of disease
  PARAMETER_VECTOR( beta_c_hat ); // 2.2 -- Coefficient for conditional dependence of disease probability on density
  PARAMETER_MATRIX( epsilon_mat ); // 2.2. -- Random site effects
  PARAMETER_VECTOR( gamma_q );     // 2.3. -- Overdispersion parameters of log-normal Poisson
  PARAMETER_VECTOR( log_sigmas );  // 2.4. -- Variance components
  // -- 0 = log SD for random site effects for density random effects
  // -- 1 = log SD for random site effects for positive density error
  // -- 2 = log SD for normal/lognormal distribution
  // -- 3 = log SD for random site effects for individual infection of parasite 1
  // -- 4 = log SD for random site effects for individual infection of parasite 2
  // -- 5 = log SD for random effects for overdispersion for each quadrat


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 3. Model objects
  // -- Model indices
  int q, i, s, p; // q = quadrat, i = individual, s = site, p = parasites
  int n_q = c_q.size(); // Sample size for number of quadrats
  int overdispersion = 0; // Index of wether to include overdispersion
  if(model == 2 | model == 6) { overdispersion = 1;}

  // -- Likelihood components
  vector<Type> jnll_comp(11); jnll_comp.setZero();
  // Slot 0 -- presence/absence
  // Slot 1 -- local density
  // Slot 2 -- random site effect for presence/absence
  // Slot 3 -- random site effect for local density
  // Slot 4 -- overdispersion
  // Slot 5 -- probability of no disease parasite 1
  // Slot 6 -- probability of no disease parasite 2
  // Slot 7 -- probability of disease parasite 1
  // Slot 8 -- probability of disease parasite 2
  // Slot 9 -- random site effect for disease probability of parasite 1
  // Slot 10 -- random site effect for disease probability of parasite 2
  Type jnll = 0;



  // ------------------------------------------------------------------------------------------------------------------------- //
  // 4. Get linear predictors
  vector<Type> linpred_pq( n_q );
  vector<Type> zero_prob_q( n_q );
  vector<Type> linpred_cq( n_q );
  vector<Type> pred_cq( n_q ); // predicted abundance

  matrix<Type> linpred_yi( n_i, n_p );
  matrix<Type> zero_prob_yi( n_i, n_p );

  for( q = 0; q < n_q; q++){
    linpred_pq( q ) = beta0(0) + (x_pq * beta_p)( q ) + epsilon_mat( 0, s_q( q )) ; // Linear predictor of presence/absence
    linpred_cq( q ) = beta0(1) + (x_cq * beta_c)( q ) + epsilon_mat( 1, s_q( q )) ; // Linear predictor of abundance
  }

  if( overdispersion == 1 ){
    for(q = 0; q < n_q; q ++){
      linpred_cq( q ) = linpred_cq( q ) + gamma_q( q );
    }
  }

  // Make sure the linear predictor is positive
  linpred_cq = exp(linpred_cq);
  if(debug == 1) {std::cerr << "E1. The first predicted non-zero catch is " <<  linpred_cq(0) <<std::endl;}

  zero_prob_q = 1 / (1 + exp( -linpred_pq )); // Probability of absence
  if(debug == 1) {std::cerr << "E2. The first zero probability is " <<  zero_prob_q(0) <<std::endl;}

  pred_cq = (1 - zero_prob_q ) * linpred_cq; // Predicted catch
  if(debug == 1) {std::cerr << "E3. The first predicted catch " <<  pred_cq(0) <<std::endl;}


  if(incl_disease == 1){
    for(i = 0; i < n_i; i++){
      for( p = 0; p < n_p; p ++){
      vector<Type> beta_temp = beta_y.row(p); // Temporary storage for multiplication
      linpred_yi( i, p ) = beta0(2 + p) + (matrix_from_array(x_i, p) * beta_temp)(i) + beta_c_hat(p) * linpred_cq( q_i( i ) ) + epsilon_mat( 2 + p, s_i( i )); // Linear predictor for infection
      zero_prob_yi( i, p  ) = 1/ (1 + exp( - linpred_yi( i, p ) ) ); // Probility of infection
      }
    }

  }
  if(debug == 1 & incl_disease == 1) {std::cerr << "E4. The first infection probability is " <<  zero_prob_yi(0 , 0) <<std::endl;}

  // ------------------------------------------------------------------------------------------------------------------------- //
  // 5. Random effects
  for( s = 0; s < n_s; s++){
    jnll_comp(2) -= dnorm( epsilon_mat( 0, s ), Type(0.0), exp(log_sigmas( 0 )), true); // Random site effects for presence/absence
    jnll_comp(3) -= dnorm( epsilon_mat( 1, s ), Type(0.0), exp(log_sigmas( 1 )), true); // Random site effects for local density

    if(debug == 1 & s == 0) {std::cerr << "J1. The first nll of random site effect for presence/absence  is " <<  jnll_comp(2) <<std::endl;}
    if(debug == 1 & s == 0) {std::cerr << "J2. The first nll of random site effects for density random is " <<  jnll_comp(3) <<std::endl;}
  }

  if( overdispersion == 1 ){
    for(q = 0; q < n_q; q ++){
      jnll_comp(4) -= dnorm( gamma_q( q ), Type(0.0), exp(log_sigmas( 5 )), true); // Overdispersion for each quadrat
      if(debug == 1 & q == 0) {std::cerr << "J3. The first nll of overdispersion is " <<  jnll_comp(4) <<std::endl;}
    }
  }

  if(incl_disease == 1){
      for( s = 0; s < n_s; s++){
        for( p = 0; p < n_p; p ++){
        jnll_comp(9 + p) -= dnorm( epsilon_mat( 2 + p, s ), Type(0.0), exp(log_sigmas( 3 + p )), true); // Random site effects for individual infection of parasite p
        if(debug == 1 & s == 0) {std::cerr << "J4. The first nll of random site effects for infection probability is " <<  jnll_comp(9) << " and " << jnll_comp(10) << std::endl;}
      }
    }
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 6. Probability of abundance conditional on fixed and random effects values
  // -- Delta-lognormal
  if(model == 0){
    // Delta-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dlognorm( c_q(q), log( exp( linpred_cq( q ) ) - pow( exp( log_sigmas( 2 ) ) , 2)/2), exp( log_sigmas( 2 ) ), true );
    }
  }

  // -- Delta-poisson
  if(model == 1){
    // Delta-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dpois( c_q( q ), linpred_cq( q ), true );
    }
  }

  // -- delta-log-normal Poisson
  if(model == 2){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }

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

  // -- ZI-poisson
  if(model == 5){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }

  // -- ZI-log-normal Poisson
  if(model == 6){
    // ZIP-model
    for( q = 0; q < n_q; q++){
      if(c_q( q ) == 0) jnll_comp(0) -= log( zero_prob_q( q ) ) + dpois( Type(0.0), linpred_cq( q ), true );
      if(c_q( q ) != 0) jnll_comp(1) -= log( 1 - zero_prob_q( q ) ) + dpois( c_q( q ), linpred_cq( q ), true );

      if(debug == 1 & q == 0) {std::cerr << "J5. The first nll of 0 abundance is " <<  jnll_comp(0) << " and positive abundance " << jnll_comp(1) << std::endl;}

    }
  }

  // -- ZI-negative-binomial
  if(model == 7){
    std::cerr<< "Model not yet implemented"<<std::endl;
    return 0;
  }


  // ------------------------------------------------------------------------------------------------------------------------- //
  // 7. Probability of infection conditional on abundance, fixed, and random effects values
  // Individual model
    if(incl_disease == 1){
      for( i = 0; i < n_i; i++){
        for( p = 0; p < n_p; p ++){
          if(y_i( i, p ) == 0) jnll_comp(5 + p) -= log( 1 - zero_prob_yi( i, p ) );
          if(y_i( i, p ) == 1) jnll_comp(7 + p) -= log( zero_prob_yi( i, p ) );

          if(debug == 1 & i == 0) {std::cerr << "J6. The first nll of no infection probability is " <<  jnll_comp(5) << " and " << jnll_comp(6) << std::endl;}
          if(debug == 1 & i == 0) {std::cerr << "J7. The first nll of infection probability is " <<  jnll_comp(7) << " and " << jnll_comp(8) << std::endl;}
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
  REPORT(pred_cq);

  ADREPORT( beta_p);
  ADREPORT( beta_c);
  ADREPORT( beta_y);

  if(debug == 1){std::cerr << "----------------------------------------------" << std::endl;}

  jnll = jnll_comp.sum();
  return jnll;


  // ------------------------------------------------------------------------------------------------------------------------- //
}
