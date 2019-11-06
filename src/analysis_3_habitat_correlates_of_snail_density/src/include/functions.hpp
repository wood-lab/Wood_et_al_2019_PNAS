// Function for to extract a layer from an array
template<class Type>
matrix<Type> matrix_from_array(array<Type> a1, int sheet){
  vector<int> a1_dim = a1.dim;
  matrix<Type> m1(a1_dim(0), a1_dim(1));

  // If array is array
  if(a1_dim.size() == 3){
    for(int row = 0; row < a1_dim(0); row++){
      for(int col = 0; col < a1_dim(1); col++){
        m1(row, col) = a1(row, col, sheet);
      }
    }
  }

  // If array is matrix
  if(a1_dim.size() == 2){
    for(int row = 0; row < a1_dim(0); row++){
      for(int col = 0; col < a1_dim(1); col++){
        m1(row, col) = a1(row, col);
      }
    }
  }
  return m1;
}

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// Function to trim matrix
template<class Type>
matrix<Type> trim_matrix(matrix<Type> m1, int nrow, int ncol){
  // dim is a DATA_IMATRIX with structure nrow, ncol, ID # of object we wish to interogate
  matrix<Type> m2(nrow, ncol);

  // If array is matrix
  for(int row = 0; row < nrow; row++){
    for(int col = 0; col < ncol; col++){
      m2(row, col) = m1(row, col);
    }
  }
  return m2;
}

// Function to trim vector
template<class Type>
matrix<Type> trim_vector(vector<Type> v1, int length){
  // dim is a DATA_IMATRIX with structure nrow, ncol, ID # of object we wish to interogate
  vector<Type> v2(length);

  // If array is matrix
  for(int i = 0; i < length; i++){
    v2(i) = v1(i);
  }
  return v2;
}


// function to get row vector from array
template<class Type>
vector<Type> col_from_3D_array(array<Type> a1, int col, int sheet){
  vector<int> a1_dim = a1.dim;
  vector<Type> v1(a1_dim(0));

  // If array is array
  if(a1_dim.size() == 3){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col, sheet);
    }
  }

  // If array is matrix
  if(a1_dim.size() == 2){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col);
    }
  }
  return v1;
}

template<class Type>
vector<Type> col_from_4D_array(array<Type> a1, int col, int sheet1, int sheet2){
  vector<int> a1_dim = a1.dim;
  vector<Type> v1(a1_dim(0));

  // If array is 4D array
  if(a1_dim.size() == 4){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col, sheet1, sheet2);
    }
  }

  // If array is 3D array
  if(a1_dim.size() == 3){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col, sheet1);
    }
  }

  // If array is matrix
  if(a1_dim.size() == 2){
    for(int row = 0; row < a1_dim(0); row++){
      v1(row) = a1(row, col);
    }
  }
  return v1;
}

// Function to turn NAs to 0s in matrix in specific columns not specified
// first_col is the first column of categorical variables, save_col are the columns to be left alone - NA otherwise
template<class Type>
vector<Type> matrix_cat_na_to_zero( matrix<Type> m1, int first_col, vector<Type> save_col){
  matrix<Type> m2 = m1;

  for(int row = first_col; row < m1.rows(); row++){ // Start index at first col
    for(int col = 0; col < m1.cols(); col++){
      if( isNA( save_col( col ) ) ){ // If column is not estimated
        if( isNA( m2(row, col) )){ // If value is NA
        m2(row, col) = 0;
      }
    }
  }
}
  return m2;
}


// Function to count NAs
template<class Type>
int n_NA( vector<Type> v1){
  int n_NA = 0;
  for(int i = 0; i < v1.size(); i++){
    if(isNA( v1( i ) )){
      n_NA += 1;
    }
  }
  return n_NA;
}

