#include <TMB.hpp>
#include "include/functions.hpp" //  Functions for indexing


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(v1);
  DATA_MATRIX(m1);
  DATA_MATRIX(m2);
  DATA_ARRAY(a1);
  DATA_ARRAY(a2);

  // 1. Vectors

  // Trim vector
  vector<Type> v2 = trim_vector(v1, 8);
  REPORT(v2);


  // 2. Matrices

  // Size of n x 1 matrix
  int m1_rows = m1.rows();  // Number of rows
  REPORT(m1_rows);

  int m1_cols = m1.cols();  // Number of columns
  REPORT(m1_cols);

  // Trim matrix
  matrix<Type> m5 = trim_matrix(m2, 50, 2);
  REPORT(m5);


  // 3. Arrays
  vector<int> a2_dim = a1.dim; // Dimensions of matrix
  REPORT(a2_dim);

  // Matrix from n1 x n2 x 1 array
  matrix<Type> m3;
  m3 = matrix_from_array(a1, 0);
  REPORT(m3);

  // Matrix from array
  matrix<Type> m4;
  m4 = matrix_from_array(a2, 0);
  REPORT(m4);

  // 4-D Array
  array<Type> a3(5,5,5,4);
  a3.setZero();
  REPORT(a3);


  return 0;
}
