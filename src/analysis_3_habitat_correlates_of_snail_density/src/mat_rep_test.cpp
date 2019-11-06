
#include <TMB.hpp>
// Space time
template<class Type>
Type objective_function<Type>::operator() (){

  // Tests how ADREPORT returns matrices
  DATA_MATRIX( x1 );
  DATA_MATRIX( x2 );
  DATA_VECTOR( y );

  PARAMETER_MATRIX( m1);
  PARAMETER( log_sigma );

  Type jnll = 0;

  for(int i = 0; i < y.size(); i ++){
    jnll -= dnorm(y(i), x1(i,0) * m1(0,0) + x1(i,1) * m1(0,1) + x2(i,0) * m1(1,0) + x2(i,1) * m1(1,1), exp(log_sigma), true);
  }

  Type m11 = m1(0,0);
  Type m12 = m1(0,1);
  Type m21 = m1(1,0);
  Type m22 = m1(1,1);

  REPORT( m11 );
  REPORT( m12 );
  REPORT( m21 );
  REPORT( m22 );


  ADREPORT(m1);
  return jnll;
}
