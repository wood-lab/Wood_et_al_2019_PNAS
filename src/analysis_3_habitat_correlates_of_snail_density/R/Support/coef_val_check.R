

# Step 0 -- Load data and packages
density_data <- read.csv("data/aggregated_data_with_treatments_indicated.csv", na.strings=c(""))
individual_data <- read.csv("data/merged_data.csv", na.strings=c(""))


# Change "sampling_site" to "site"
colnames(density_data)[colnames(density_data) == "sampling_site"] <- "site"

# Step 1 -- variables we want to save
enviro_vars <- c("emergent_stems_cleaned", "floating_veg_g_cleaned", "secchi_depth_cm", "pH", "nitrate_mg.L_NO3", "nitrite_mg.L_NO2", "phosphate_mg.L_PO4", "typha", "phrag", "cyperus", "poa_grass", "ludwigia", "nymphaea", "salvinia", "ceratophyllum", "potamogeton", "pistia", "azolla", "ipomea", "lemna")
site_vars <- c("point_ID", "site", "sampling_date", "field_mission")

# Subset columns
density_data <- density_data[, c(enviro_vars, site_vars)]
individual_data <- individual_data[, c(enviro_vars, site_vars)]

# Merge datasets
merge_vars <- c("point_ID", "site", "field_mission") # variables used to link datasets
combined_dat <- merge(individual_data, density_data, by = merge_vars)

# Check if they are equal
equal_check <- data.frame(variable = enviro_vars, n_mismatch = rep(NA, length(enviro_vars)),  merged_vars = rep(NA, length(enviro_vars)), aggregated_mismatch = rep(NA, length(enviro_vars)), point_ID_site_field_mission = rep(NA, length(enviro_vars)))


for(i in 1:length(enviro_vars)){
  merged_var <- paste0(enviro_vars, ".x")[i] # Variable from merged dataset (individual snails)
  aggregated_var <- paste0(enviro_vars, ".y")[i] # Variable from aggregated count dataset

  equal_test <- as.numeric(as.character(combined_dat[,merged_var])) == as.numeric(as.character(combined_dat[,aggregated_var])) # Test if equal
  mismatch_rows <-  which(equal_test == F) # Which are not equal?
  equal_check$n_mismatch[i] <- length(mismatch_rows)
  equal_check$merged_vars[i] =  paste(combined_dat[mismatch_rows, merged_var], collapse = ", ") # Values that do not match from merged dataset
  equal_check$aggregated_mismatch[i] =  paste(combined_dat[mismatch_rows, aggregated_var], collapse = ", ") # Values that do not match from aggregated dataset
  equal_check$point_ID_site_field_mission[i] = paste(apply( as.matrix(combined_dat[mismatch_rows, merge_vars]) , 1 , paste , collapse = "_" ), collapse = "; ")
}

write.csv(equal_check, file = "data/data_mismatch_id.csv")


