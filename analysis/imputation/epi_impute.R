devtools::load_all()

epi1 = import.data("IC PC clocks file from Thom - repeats removed - string variables destringed.dta", c(
  uncchdid, pchorvath1, dplyr::all_of(c(clock.list.full, pc.clock.list.full, clock.list.resid)), cell.props
))

# imputing wealth information for the rolling average
imp <- mice::mice(data = epi1, m = 5, method = 'rf', seed = 500)

# completing the data using the iteration closest to the mean
epi.completed <- mice::complete(imp, 5)

# checking imputation
epi_imputed <- epi1 %>%
  mutate(count_NA = rowSums(is.na(.))) %>% # counting missing values
  mutate(test.mean0 = rowMeans(dplyr::select(., (contains("pheno") | contains("grim"))), na.rm=TRUE)) %>% # calculating mean without imputation
  dplyr::select(uncchdid, test.mean0, count_NA) %>%
  left_join(epi.completed, by = 'uncchdid') %>%
  mutate(test.mean1 = rowMeans(dplyr::select(., (contains("pheno") | contains("grim"))), na.rm=TRUE)) %>% # calculating mean after imputation
  dplyr::select(uncchdid, test.mean0, test.mean1, count_NA, (contains("pheno") | contains("grim"))) %>%
  mutate(
    mean_diff = test.mean1 - test.mean0
  )

t.test(epi1$dnamgrimage, epi.completed$dnamgrimage, var.equal = TRUE)
t.test(epi1$pcgrimage, epi.completed$pcgrimage, var.equal = TRUE)

save(epi.completed, file = getfp("memadeR/epicompleted.RData", memade = TRUE))
