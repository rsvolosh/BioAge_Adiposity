library(stats)
library(dplr)
library(psych)

clock.list <- c('c1', 'c2', 'c3', 'c4', 'c5', 'c6') # clock list (any number of clocks, must be present as column names in dataframe)

n.clocks <- length(clock.list) # number of clocks

# dataframe must contain exact chronological age column (call it c.age)

# make sure all residuals mean the same thing (for example, if you have dnamtl, you should multiply it by -1).

# create standardised residuals function
make.residuals0 = function(dataframe, clock, age, id_col, output) {
  if (is.character(clock)) {
    column_name <- clock
    # Check if the column name exists in the dataframe
    if (!(column_name %in% colnames(dataframe))) {
      stop(paste("Column", column_name, "does not exist in the dataframe."))
    }
    # Convert column name to symbol
    clock <- sym(column_name)
    
  } else {
    # If the input is not a string, assume it's a column object
    column_name <- deparse(substitute(column_input))
    # Check if the column object exists in the dataframe
    if (!(column_name %in% colnames(data))) {
      stop(paste("Column", column_name, "does not exist in the dataframe."))
    }
  }
  dataframe = dataframe %>% filter(!is.na({{clock}})) %>% filter(!is.na({{age}}))
  mod = lm(reformulate(deparse(substitute(age)), deparse(substitute(clock))), data = dataframe)
  dataframe[[paste0("rs.",column_name)]] = rstudent(mod)
  dataframe[[paste0("rst.",column_name)]] = rstandard(mod)
  dataframe[[paste0("r.",column_name)]] = residuals(mod)
  
  if (output == "plot") {
    return(plot(mod))
  } else if (output == "student") {
    return(dplyr::select(dataframe, {{id_col}}, {{clock}}, paste0("rs.",column_name), {{age}}))
  } else if (output == "all") {
    return(dplyr::select(dataframe, {{id_col}}, {{clock}}, paste0("r.",column_name), paste0("rs.",column_name), paste0("rst.",column_name), {{age}}))
  } else if (output == "standard") {
    return(dplyr::select(dataframe, {{id_col}}, {{clock}}, paste0("rst.",column_name), {{age}}))
  } else {
    return(dplyr::select(dataframe, {{id_col}}, {{clock}}, paste0("r.",column_name), {{age}}))
  }
}

# create standardised residuals
for (i in 1:n) {
  residuals1 <- residuals1 %>% left_join(make.residuals0(dataframe = dataframe, clock = clock.list[i], age = c.age, "rs"), by = c("id_col", "c.age")) #id_col is placeholder name for ID column
}

# here, you may check your standardised residuals against z-scores of residuals, correlations, check if removing outliers is warranted, consider data cleaning

# you may want to combine clock estimates by summing them
residuals2 <- residuals1 %>% mutate(
  added.nonpc = (rowSums(dplyr::select(., starts_with("rs."))) - mean(rowSums(dplyr::select(., starts_with("rs."))))) / sd(rowSums(dplyr::select(., starts_with("rs."))))
)

# examine correlations
residuals2 %>% dplyr::select(starts_with("rs.")) %>% cor()

# KMO test
residuals2 %>% dplyr::select(starts_with("rs.")) %>% cor() %>% KMO()

# eigenvalues
residuals2 %>% dplyr::select(starts_with("rs.")) %>% cor() %>% eigen()

# screeplots
residuals2 %>% dplyr::select(starts_with("rs.")) %>% fa.parallel(fm = "ml", fa = "fa")

# single-factor solution
# placeholder name for new combined clock variable factor.age
onefac <- fa(cor(dplyr::select(residuals2, starts_with("rs."))), nfactors = 1)
residuals2[["factor.age"]] <- factor.scores(dplyr::select(residuals2, starts_with("rs.")), onefac$loadings, method = "Bartlett")$scores[,1]

# correlations of all measures
residuals2 %>% dplyr::select(-id_col) %>% cor()
