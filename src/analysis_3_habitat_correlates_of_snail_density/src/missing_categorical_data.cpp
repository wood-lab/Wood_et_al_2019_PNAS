
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // Setup
  DATA_VECTOR( x_i ) ; 				// Design matrix of categorical variables with missing values
  DATA_VECTOR( y_i ) ; 				// Response

  PARAMETER_VECTOR( x_i_missing );  // Missing values of categorical covariates to be substituted
  PARAMETER( beta_log );			// Linear predictor for logistic of each categorical variable
  PARAMETER( beta_x );		// Parameter vector for coefficients
  PARAMETER( log_sd ); // Log variance

  // Model objects
  vector<Type> x_i_hat = x_i; // New design matrix with substituted values
  Type jnll = 0; // Joint log-likelihood
  Type pr = 1/(1 + exp(-p_j)); // Probability of factor being 0
  int missing_ind = 0; // Index for missing value

  // Fit model
  for(int i = 0; row < y_i.size() ; row++){

    // Fit available data
    if( !isNA( x_i( i ) ) ){

      // Data likelihood
      if(x_i( i ) == 0) jnll -= log( pr ); // Probability of covariate being 0
      if(x_i( i ) == 1) jnll -= log( 1 - pr ); // Probability of covariate being 1

      // Model likelihood
      jnll -= dnorm( y_i( row ), (x_i * beta_x)( row ), exp( log_sd ), true);
    }

    // Fit missing data
    if( isNA( x_i( i ) ) ){
      x_i_hat( i ) = x_i_missing( missing_ind ); // Input missing N(0,1) random effect
      missing_ind = missing_ind + 1

      // Likelihood if missing value = 0
      jnll -= dnorm( y_i( row ), (x_i * beta_x)( row ), exp( log_sd ), true) + log( pr ) ;

      // Likelihood if missing value = 1
      jnll -= dnorm( y_i( row ), (x_i * beta_x)( row ), exp( log_sd ), true) + log( 1 - pr ) ;
    }
  }


  // Report
  REPORT(x_i_hat);
  return jnll;
}
