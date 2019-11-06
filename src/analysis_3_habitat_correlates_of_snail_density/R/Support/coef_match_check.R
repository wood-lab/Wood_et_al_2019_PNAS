Mode = function(x){
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = x[1]
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}


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

# Get dataset of unique site and field miessions
unique_sites <- unique(density_data[,c("site", "field_mission")])
for(j in 1:length(enviro_vars)){
  unique_sites$new_col <- NA
  names(unique_sites)[which(names(unique_sites) == "new_col")] <- enviro_vars[j]
}


for(i in 1:nrow(unique_sites)){
  dat_sub <- density_data[which(density_data$site == unique_sites$site[i] & density_data$field_mission == unique_sites$field_mission[i]),]
  for(j in 1:length(enviro_vars)){
    row_sub <- as.numeric(as.character(dat_sub[,enviro_vars[j]]))
    unique_sites[i,enviro_vars[j]] <- sum(row_sub %in% Mode(row_sub))/length(row_sub)
  }
}

write.csv(unique_sites, file = "data/percent_agreement_aggregated.csv")


# Get dataset of unique site and field miessions
unique_sites <- unique(individual_data[,c("site", "field_mission")])
for(j in 1:length(enviro_vars)){
  unique_sites$new_col <- NA
  names(unique_sites)[which(names(unique_sites) == "new_col")] <- enviro_vars[j]
}


for(i in 1:nrow(unique_sites)){
  dat_sub <- individual_data[which(individual_data$site == individual_data$site[i] & individual_data$field_mission == unique_sites$field_mission[i]),]
  for(j in 1:length(enviro_vars)){
    row_sub <- as.numeric(as.character(dat_sub[,enviro_vars[j]]))
    unique_sites[i,enviro_vars[j]] <- sum(row_sub %in% Mode(row_sub))/length(row_sub)
  }
}

write.csv(unique_sites, file = "data/percent_agreement_merged.csv")

