setwd("src")
library(TMB)
library(TMBdebug)
library(TMBhelper)
version = "spatio_temporal_model_spde"
if(data_list$debug == 1){
  #dyn.unload(version)
  #file.remove(paste0(version,".dll"))
  #file.remove(paste0(version,".o"))
}
TMB::compile( paste0(version,".cpp") )
dyn.load(version)
setwd("../")


data_list$incl_disease = 1
params <- build_params(model, data_list, incl_disease = 1)
map <- build_map(model, data_list, incl_disease = 1, params, space = T, space_time = T)

Obj = TMBdebug::MakeADFun( data = data_list, parameters = params, DLL = version, map = map, random = random)
