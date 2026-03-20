# this is a test file with R scripts edited from draft1B, using functions
# from the AAA package

# import necessary libraries later
library(AAA)
library(tidyr)

# importing imputed datasets
load(getfp("memadeR/epicompleted.RData", memade = TRUE)) # epigenetic clocks
load(getfp("memadeR/wealth2.RData", memade = TRUE)) # wealth
load(getfp("memadeR/bf.imputed2.RData", memade = TRUE)) # pregnancy and breastfeeding
load(getfp("memadeR/anthro.RData", memade = TRUE)) # anthro file

# import data pre-impute
smoking = import.data(getfp("morhealt.dta"), ofensmok) %>% mutate(
  cursmok = case_when(
    ofensmok == -7 ~ 1, # former smoker
    ofensmok >= -6 ~ 2, # current smoker
    TRUE ~ 0 # non-smoker
  )
) %>% mutate_at(.vars = vars(cursmok), .funs = as_factor)

sex = import.data(getfp("cebuheightlong.dta"), c(sex0 = sex, year)) %>% group_by(uncchdid) %>%
  mutate(sex = min(sex0, na.rm=TRUE)) %>%
  distinct(uncchdid, sex) %>% filter(sex %in% c(0,1)) %>%
  mutate(sex = factor(sex, levels = c(0,1)))

# age at blood collection in 2005
chron.age <- import.data(getfp("BLOODfinal.DTA"), day = "dayblood", month = "monthblood", year = "yearblood", name_date = "date.blood") %>%
  left_join(
    import.data(getfp("birthdate.dta"), birthdate) %>% mutate(birthdate = as.Date("1960-01-01") + birthdate),
    by = "uncchdid"
    ) %>% mutate(
  c.age = as.numeric(date.blood - birthdate)/365.25
) %>% filter(!is.na(c.age)) %>% distinct()

# combining datasets
# epi2 is epi1 merged with demographic vars, anthro, pregnancy, smoking, and wealth
epi2 <- Reduce(function(x, y) left_join(x, y, by = "uncchdid", relationship = "many-to-many"),
       list(sex, chron.age, anthro, smoking, wealth2, bf.imputed2), init = epi.completed) %>%
  arrange(desc(pregord)) %>% arrange(uncchdid) %>%
  distinct(uncchdid, .keep_all = TRUE) %>% mutate(across(
    # for women, replace NAs in relevant variables
    c(cnum.prg, totmonthsbf, ctotmonthsbf, everpreg, everbirth, num.preg.total, num.preg.comp, num.kids, pregnow, mo.bf, mo.bf2, mo.pp, pp2),
    .fns = ~ ifelse(is.na(.), 0, .)
  )) %>% filter(!is.na(sex)) %>% mutate(dnamtladjage.n = -dnamtladjage)

# creating residuals (prior to checking for outliers)
residuals1 <- epi2 %>% dplyr::select(uncchdid, c.age)

for (i in 1:length(c(clock.list, pc.clock.list))) {
  residuals1 <- residuals1 %>% left_join(make.residuals0(epi2, c(clock.list, pc.clock.list)[i], c.age, "all"), by = c("uncchdid", "c.age"))
}

# TO REMOVE LATER

# checking correlations between standardised residuals and z-scores of residuals
for (i in 1:length(c(clock.list, pc.clock.list))) {
  residuals1 <- residuals1 %>% mutate(across(
    (starts_with("r.")),
    .fns = list(z = ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)),
    .names = "{fn}.{col}"
    ))
}

residual.summary <- do.call(rbind, lapply(c(clock.list, pc.clock.list), function(x) cbind(
                                            "clock" = x,
                                            "cor" = cor.res(residuals1, x)$estimate,
                                            "95% CI" = paste0(cor.res(residuals1, x)$conf.int[1], ", ", cor.res(residuals1, x)$conf.int[2])
                                            ))) %>% as.data.frame()

# END TO REMOVE LATER
# I will henceforth use standardised residuals in my analyses instead of z-scores, since it is simpler, easier to explain, and has very high correlation.

# negating pcdnamtl to match the signs of other residuals due to opposite biological meaning
residuals2 <- residuals1 %>% mutate(across(contains("mtl"), .fns = ~ -.)) %>%
  # z-scoring dunedin_po_am to make it comparable with other scocres
  left_join(dplyr::select(epi2, uncchdid, sex, dunedinpoam_45), by = "uncchdid") %>%
    mutate(
      # keeping up with the naming convention to facilitate coding, I am naming the scaled dunedinpoam_45 "rs.dunedin", even though it is not a residual
      rs.dunedin = (dunedinpoam_45 - mean(dunedinpoam_45, na.rm = TRUE))/sd(dunedinpoam_45, na.rm = TRUE),
      # same for PC clocks; dunedinpoam doe snot have a PC version
      rs.pcdunedin = (dunedinpoam_45 - mean(dunedinpoam_45, na.rm = TRUE))/sd(dunedinpoam_45, na.rm = TRUE),
      # ditto
      z.r.pcdunedin = (dunedinpoam_45 - mean(dunedinpoam_45, na.rm = TRUE))/sd(dunedinpoam_45, na.rm = TRUE),
      z.r.dunedin = (dunedinpoam_45 - mean(dunedinpoam_45, na.rm = TRUE))/sd(dunedinpoam_45, na.rm = TRUE)
    )

# creating epi3 (epi2 + residuals) to export for work with Stata
epi3 <- epi2 %>% left_join(dplyr::select(residuals2, uncchdid, starts_with("rs."), starts_with("z.r."), starts_with("r.")), by = "uncchdid")
epi3s <- epi3 %>% rename_with(~ gsub('\\.', '_', .x))
write_dta(epi3s, getfp("memadeStata/epi3.dta", memade = TRUE))

# non-PC residuals
residuals2 %>% dplyr::select(starts_with("rs."), -contains("pc")) %>% cor()
# sex-stratified
residuals2 %>% filter (sex==0) %>% dplyr::select(starts_with("rs."), -contains("pc")) %>% cor() # females
residuals2 %>% filter (sex==1) %>% dplyr::select(starts_with("rs."), -contains("pc")) %>% cor() # males

# PC residuals
residuals2 %>% dplyr::select(starts_with("rs.pc")) %>% cor()
# sex-stratified
residuals2 %>% filter (sex==0) %>% dplyr::select(starts_with("rs.pc")) %>% cor() # females
residuals2 %>% filter (sex==1) %>% dplyr::select(starts_with("rs.pc")) %>% cor() # males

# KMO test
residuals2 %>% dplyr::select(starts_with("rs."), -contains("pc")) %>% cor() %>% KMO() # non-PC clocks
residuals2 %>% dplyr::select(starts_with("rs.pc")) %>% cor() %>% KMO() # PC clocks

# eigenvalues
residuals2 %>% dplyr::select(starts_with("rs."), -contains("pc")) %>% cor() %>% eigen()# non-PC clocks
residuals2 %>% dplyr::select(starts_with("rs.pc")) %>% cor() %>% eigen() # PC clocks

# screeplots
residuals2 %>% dplyr::select(starts_with("rs."), -contains("pc")) %>% fa.parallel(fm = "ml", fa = "fa")
residuals2 %>% dplyr::select(starts_with("rs.pc")) %>% fa.parallel(fm = "ml", fa = "fa")

# factor solution
onefac.nonpc <- fa(cor(dplyr::select(residuals2, starts_with("rs."), -contains("pc"))), nfactors = 1)
residuals2[[factors.nonpc[1]]] <- factor.scores(dplyr::select(residuals2, starts_with("rs."), -contains("pc")), onefac.nonpc$loadings, method = "Bartlett")$scores[,1]

# factor solution
onefac.pc <- fa(cor(dplyr::select(residuals2, starts_with("rs.pc"))), nfactors = 1)
residuals2[[factors.pc[1]]] <- factor.scores(dplyr::select(residuals2, starts_with("rs.pc")), onefac.pc$loadings, method = "Bartlett")$scores[,1]

# looking at correlations
cor(dplyr::select(residuals2, contains("1fac.f1")))
cor(dplyr::select(residuals2, nonpc1fac.f1, matches("^rs.d")))
cor(dplyr::select(residuals2, pc1fac.f1, matches("^rs.pc")))

# merging with epi2
epi3 <- epi2 %>% left_join(
  dplyr::select(residuals2, uncchdid, contains("rs."), starts_with("r."), matches("pc\\d"), contains("added")),
  by = "uncchdid") %>% filter(pregnow != 1) %>% filter(!is.na(bmic)) %>% mutate(across(
    c(bmic, waist, whr, wthr),
    .fns = list(z = ~ (. - mean(., na.rm = TRUE))/sd(., na.rm = TRUE)),
    .names = "{fn}.{col}"
  )) %>% mutate( # categorising
    # Asian-specific BMI cut-points
    bmic3cat = case_when(
      bmic < 23 ~ 0, bmic < 25 ~ 1, !is.na(bmic) ~ 2, TRUE ~ NA),
    bmic4cat = case_when(
      bmic < 18.5 ~ 3, bmic < 23 ~ 0, bmic < 25 ~ 1, !is.na(bmic) ~ 2, TRUE ~ NA),
    # non-Asian BMI cut-points
    bmic3cat.na = case_when(
      bmic < 25 ~ 0, bmic < 30 ~ 1, !is.na(bmic) ~ 2, TRUE ~ NA),
    bmic4cat.na = case_when(
      bmic < 18.5 ~ 3, bmic < 25 ~ 0, bmic < 30 ~ 1, !is.na(bmic) ~ 2, TRUE ~ NA),
    # Asian-specific WC from paper 3b
    waist.high = ifelse(sex == 0, ifelse(waist < 77, 0, ifelse(!is.na(waist), 1, NA)), # this is for rural areas, don't use but mention in discussion
                        ifelse(waist < 85, 0, ifelse(!is.na(waist), 1, NA))),
    # non-Asian WC
    waist.high.na = ifelse(sex == 0, ifelse(waist < 80, 0, ifelse(!is.na(waist), 1, NA)), # use this
                           ifelse(waist < 90, 0, ifelse(!is.na(waist), 1, NA))),
    # Asian-specific WHR from paper 3b
    whr.high = ifelse(sex == 0, ifelse(whr < 0.85, 0, ifelse(!is.na(whr), 1, NA)), # this is for rural areas, don't use but mention in discussion
                      ifelse(whr < 0.91, 0, ifelse(!is.na(whr), 1, NA))),
    # Asian WHR WHO from the big cohort Filipino diabetes paper
    whr.high.na = ifelse(sex == 0, ifelse(whr < 0.80, 0, ifelse(!is.na(whr), 1, NA)),
                         ifelse(whr < 0.90, 0, ifelse(!is.na(whr), 1, NA))),
    # WtHR
    wthr.high = ifelse(wthr < 0.50, 0, ifelse(!is.na(wthr), 1, NA))
  ) %>% group_by(sex) %>% ungroup() %>%
  mutate(
    f.bmic3cat = factor(bmic3cat, levels = c(0,1,2), labels = c('Not obese/overweight', 'Overweight', 'Obese')),
    f.bmic4cat = factor(bmic4cat, levels = c(0,1,2, 3), labels = c('Normal', 'Overweight', 'Obese', 'Underweight')),
    f.bmic3cat.na = factor(bmic3cat.na, levels = c(0,1,2), labels = c('Not obese/overweight', 'Overweight', 'Obese')),
    f.bmic4cat.na = factor(bmic4cat.na, levels = c(0,1,2, 3), labels = c('Normal', 'Overweight', 'Obese', 'Underweight')),
    f.waist.high = factor(waist.high, levels = c(0,1), labels = c("Normal", "High")),
    f.waist.high.na = factor(waist.high.na, levels = c(0,1), labels = c("Normal", "High")),
    f.whr.high = factor(whr.high, levels = c(0,1), labels = c("Normal", "High")),
    f.whr.high.na = factor(whr.high.na, levels = c(0,1), labels = c("Normal", "High")),
    f.wthr.high = factor(wthr.high, levels = c(0,1), labels = c("Normal", "High")),
    f.uw = factor(uw, levels = c(0,1), labels = c("Not underweight", "Underweight"))
  ) %>% mutate(across(c(contains("\\dcat"), contains(".high"), uw), as.factor))

save(epi3, file = getfp("memadeR/epi3.RData", memade = TRUE))

