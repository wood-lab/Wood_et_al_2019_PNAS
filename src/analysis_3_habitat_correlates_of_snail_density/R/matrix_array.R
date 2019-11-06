# Function to test matrix-vector indexing for debugging


setwd("src")
library(TMB)
library(TMBdebug)
compile( "matrix_array.cpp" )
dyn.load( dynlib("matrix_array") )
setwd("../")

n = 100
v1 <- runif(n, 1, 100)
m1 <- as.matrix(runif(n, 1, 100))
m2 <- matrix(runif(n * 2, 1, 100), nrow = n)
a1 <- matrix(runif(n * 2, 1, 100), nrow = n)
a2 <- array(runif(n * 4, 1, 100), dim = c(n, 2, 2))
data <- list(
  v1 = v1,
  m1 = m1,
  m2 = m2,
  a1 = a1,
  a2 = a2
)

# Step 4 -- Fit model
Obj = MakeADFun( data = data, parameters=list(),type="Fun", DLL="matrix_array")
report <- (Obj$report())  # Note: order of variables NOT the same as .cpp file
