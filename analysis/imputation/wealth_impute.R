
devtools::load_all()

# importing data
wealth = import.data("../anthro et al/cebupc.dta", c(survey = year, contains("pc"))) %>%
  filter(survey <= 2005) %>%
  pivot_wider(id_cols = 'uncchdid', names_from = 'survey', values_from = c('pc', 'pc2', 'pc3')) %>%
  dplyr::select(uncchdid, contains("pc"))

# imputing wealth information for the rolling average
imp <- mice::mice(data = wealth, m = 5, method = 'rf', seed = 500)

# completing the data using the iteration closest to the mean
wealth.completed <- mice::complete(imp, 5)

# checking imputation
wealth2 <- wealth %>%
  dplyr::select(-matches("pc[0-9]_")) %>% # only using the first principal component
  mutate(count_NA = rowSums(is.na(.))) %>% # counting missing values
  mutate(pc.mean0 = rowMeans(dplyr::select(., contains("pc_")), na.rm=TRUE)) %>% # calculating mean without imputation
  dplyr::select(uncchdid, pc.mean0, count_NA) %>%
  left_join(dplyr::select(wealth.completed, -matches("pc[0-9]_")), by = 'uncchdid') %>%
  mutate(pc.mean1 = rowMeans(dplyr::select(., contains("pc_")), na.rm=TRUE)) %>% # calculating mean after imputation
  dplyr::select(uncchdid, pc.mean0, pc.mean1, pc_2002, pc_2005, count_NA) %>%
  mutate(
    mean_diff = pc.mean1 - pc.mean0,
    d.pc = pc_2005 - pc_2002
    )

# The imputation does not make a major difference for those present in the 2005 survey
pc.summary <- wealth2 %>% filter(uncchdid %in% unique(epi1$uncchdid)) %>% # getting IDs from those who have outcome variable
  group_by(count_NA) %>% summarise(
    n = n(),
    mean.pc.noimpute = mean(pc.mean0),
    sd.pc.noimpute = sd(pc.mean0),
    mean.pc.imputed = mean(pc.mean1),
    sd.pc.imputed = sd(pc.mean1),
    mean.pc.diff = mean(mean_diff),
    sd.pc.diff = sd(mean_diff),
    mean.pc.prc = mean(mean_diff/pc.mean0) # percent "error"
  )

# exporting wealth2 as RData and Stata file
save(wealth2, file = getfp("memadeR/wealth2.RData", memade = TRUE))

# modifying file for Stata & exporting a Stata file version
wealth2s <- wealth2 %>% rename_with(~ gsub('\\.', '_', .x))
write_dta(wealth2s, getfp("memadeStata/wealth2.dta", memade = TRUE))
