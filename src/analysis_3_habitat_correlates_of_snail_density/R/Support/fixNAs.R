
fixNAs <- function(data, vars.of.import, exclude_rows = c()) {
  # data <- dt %>% filter(field_mission == 2)
  # data <- dt

  newdata <- data
  # vars.of.import <- vars.of.import[1:5]

  # summary(data[, vars.of.import[1:5]])
  # summary(newdata[, vars.of.import[1:5]])

  for (var in vars.of.import) {
    # var <- "secchi_depth_cm"

    print(var)
    any.NA.var <- any(is.na(data[, var]))
    print(any.NA.var)
    # If we have any NAs...
    if(any.NA.var) {
      indices.NA.var <- which(is.na(data[, var]))
      print(indices.NA.var)
      # For each of the NAs...
      for(i in indices.NA.var) {
        # i <- indices.NA.var[1]
        print(i)
        fm.for.i <- data$field_mission[i]
        site.for.i <- data$site[i]
        # print(fm.for.i); print(site.for.i)
        # If not all of the values for that FM for that variable are NA...
        if(!all(is.na(data[data$field_mission == fm.for.i & data$site == site.for.i, var]))) {
          # Replace the NA by the median across the other quadrats in that FM (numeric variables)
          if(grepl("cont_covar", var) ) {# if(is.numeric(data[, var]) ) {
            new.value <-  median(data[data$field_mission == fm.for.i & data$site == site.for.i, var], na.rm=T)
            stopifnot(!is.na(new.value))
            newdata[i, var] <-new.value
            print(paste0("Now the value is ", newdata[i, var]))
          } else if(grepl("cat_covar", var)){ # if(is.factor(data[, var])) {
            new.value <- Mode(data[data$field_mission == fm.for.i & data$site == site.for.i, var])
            stopifnot(!is.na(new.value))
            newdata[i, var] <- new.value
            print(paste0("Now the value is ", newdata[i, var]))
          }
          # If however all the values for that FM for that var are NA
        } else if(all(is.na(data[data$field_mission == fm.for.i & data$site == site.for.i, var]))) {
          if(!all(is.na(data[data$site == site.for.i, var][-exclude_rows]))){
            # Replace by the median/Mode for that site from all the nonNA quadrats in other field missions in non-excluded bits
            if(grepl("cont_covar", var) ) {# if(is.numeric(data[, var]) ) {
              new.value <- median(data[data$site == site.for.i, var][-exclude_rows], na.rm=T)
              stopifnot(!is.na(new.value))
              newdata[i, var] <- new.value
              print(paste0("Now the value is ", newdata[i, var]))
            } else if(grepl("cat_covar", var)){ # if(is.factor(data[, var])) {
              new.value <- Mode(data[data$site == site.for.i, var][-exclude_rows], na.rm=T)
              stopifnot(!is.na(new.value))
              newdata[i, var] <- new.value
              print(paste0("Now the value is ", newdata[i, var]))
            }
          }
          else if(all(is.na(data[data$site == site.for.i, var][-exclude_rows]))){
            # if(!all(is.na(data[data$site == site.for.i, var]))){
            #   # Replace by the median/Mode for for that site across field mission
            #   if(grepl("cont_covar", var) ) {# if(is.numeric(data[, var]) ) {
            #     new.value <- median(data[data$site == site.for.i, var], na.rm=T)
            #     stopifnot(!is.na(new.value))
            #     newdata[i, var] <- new.value
            #     print(paste0("Now the value is ", newdata[i, var]))
            #   } else if(grepl("cat_covar", var)){ # if(is.factor(data[, var])) {
            #     new.value <- Mode(data[data$site == site.for.i, var], na.rm=T)
            #     stopifnot(!is.na(new.value))
            #     newdata[i, var] <- new.value
            #     print(paste0("Now the value is ", newdata[i, var]))
            #   }
            # } else if(all(is.na(data[data$site == site.for.i, var]))){
            # Replace by the median/Mode for all sites that are not excluded
            if(grepl("cont_covar", var) ) {# if(is.numeric(data[, var]) ) {
              new.value <- median(data[-exclude_rows, var], na.rm=T)
              stopifnot(!is.na(new.value))
              newdata[i, var] <- new.value
              print(paste0("Now the value is ", newdata[i, var]))
            } else if(grepl("cat_covar", var)){ # if(is.factor(data[, var])) {
              new.value <- Mode(data[-exclude_rows, var], na.rm=T)
              stopifnot(!is.na(new.value))
              newdata[i, var] <- new.value
              print(paste0("Now the value is ", newdata[i, var]))
              # }
            }
          }
        }
      }
    }
  }

  return(newdata)

  # summary(data[, vars.of.import])
  # summary(newdata[, vars.of.import])

}
