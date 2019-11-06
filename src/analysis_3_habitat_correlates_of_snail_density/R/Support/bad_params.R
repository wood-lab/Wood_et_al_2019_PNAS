# Function to return what parameters are giving issues

find_bad_params <- function(Obj){
  param_check <- TMBhelper::Check_Identifiable(Obj)
  bad_params <- param_check$BadParams[which(param_check$BadParams$Param_check == "Bad"),]
  return(bad_params)
}
