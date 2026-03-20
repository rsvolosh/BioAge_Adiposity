devtools::load_all()

# importing data
#####
epi1 = import.data(getfp("IC PC clocks file from Thom - repeats removed - string variables destringed.dta"), everything())

load(getfp("memadeR/wealth2.RData", memade = TRUE))

anthro = full_join(
  import.data(getfp("anthdiet2002.dta"), c(lag.weight = weight, lag.height = height, lag.waist = waist, lag.hip = hip)) %>%
    mutate(
      lag.bmic = case_when(lag.weight > 0 & lag.height > 0 ~ lag.weight/(lag.height/100)^2, TRUE ~ NA),
      lag.whr = case_when(lag.waist > 0 & lag.hip > 0 ~ lag.waist/lag.hip, TRUE ~ NA),
      lag.wthr = case_when(lag.waist > 0 & lag.height > 0 ~ lag.waist/lag.height, TRUE ~ NA)
    ),
  import.data(getfp("anthdiet.dta"), c(weight, height, waist, hip)) %>%
    mutate(
      bmic = case_when(weight > 0 & height > 0 ~ weight/(height/100)^2, TRUE ~ NA),
      whr = case_when(waist > 0 & hip > 0 ~ waist/hip, TRUE ~ NA),
      wthr = case_when(waist > 0 & height > 0 ~ waist/height, TRUE ~ NA)
    ),
  by = "uncchdid"
) %>%
  mutate(
    d.bmi = case_when(!is.na(bmic) & !is.na(lag.bmic) ~ bmic - lag.bmic, TRUE ~ NA),
    d.waist = waist - lag.waist,
    d.whr = whr - lag.whr,
    d.wthr = wthr - lag.wthr,
    uw = ifelse(bmic < 18.5, 1, 0)
  )

preg = import.data(getfp("PREGHIST.DTA"), everything(), day = "daypregh", month = "monpregh", year = "yrpreghi", name_date = "dateintvw") %>% distinct(uncchdid, pregord, .keep_all = TRUE)
preg2009 = import.data(getfp("PREGHIST2009.dta"), everything())
blood.dates = import.data(getfp("BLOODfinal.DTA"), day = "dayblood", month = "monthblood", year = "yearblood", name_date = "date.blood")
smoking = import.data(getfp("morhealt.dta"), everything())

sex = import.data(getfp("cebuheightlong.dta"), c(sex, year)) %>% group_by(uncchdid) %>%
  mutate(sex2 = min(sex, na.rm=TRUE)) %>%
  filter(sex %in% c(0,1)) %>%
  mutate(sex = factor(sex2, levels = c(0,1))) %>%
  distinct(uncchdid, sex)

bdays = import.data(getfp("birthdate.dta"), birthdate) %>% mutate(birthdate = as.Date("1960-01-01") + birthdate)
load('/Users/rvolosh/Documents/Papers/Cebu Epigenetics/anthro et al/memadeR/df7.RData') # to get urbanicity

#####
# prep for bf/pregnancy imputation
# age at blood collection in 2005
chron.age <- blood.dates %>% left_join(bdays, by = "uncchdid") %>% mutate(
  c.age = as.numeric(date.blood - birthdate)/365.25
) %>% filter(!is.na(c.age)) %>% distinct()

# compare dates of interview and blood collection
interview.dates.c <- dplyr::select(preg, uncchdid, dateintvw) %>%
  left_join(blood.dates, by = 'uncchdid') %>%
  # different between initial interview and blood collection
  mutate(diff = date.blood - dateintvw) %>% distinct()

# cleaning data (based on conversations with Isabelita Bas at OPS)
preg[preg$uncchdid==20208 & preg$pregord==2,]$yearprg <- 2002
preg[preg$uncchdid==20429 & preg$pregord==2,]$yearprg <- 2004
preg[preg$uncchdid==21001 & preg$pregord==2,]$yearprg <- 2000
preg[preg$uncchdid==21308,]$yearprg <- 2002

preg1 <- preg %>% left_join(blood.dates, by='uncchdid') %>% mutate(
  # who is pregnant now
  pregnow = case_when(yearprg == -9 ~ 1, TRUE ~ 0),
  # create date of pregnancy termination
  dateprgtrm = as.Date(
    paste0(dayprg,"/", monthprg, "/", yearprg), format = '%d/%m/%Y')
) %>% mutate(
  # post-partum (pp)
  # basic 1 - treating miscarriage the same way as birth
  pp.b = case_when(pregnow==1 ~ 5000, TRUE ~ as.numeric(date.blood - dateprgtrm))
) %>% arrange(pp.b) %>% arrange(uncchdid)

# if never stopped breastfeeding kid #1, but kid #2 already born and also breastfed, need different coding
fixed.bf <- preg1 %>% pivot_wider(id_cols = "uncchdid", names_from = "pregord", values_from = c("stillbf", "daysbf", "dateprgtrm")) %>%
  filter(stillbf_1 == 1 & stillbf_2 == 1) %>% mutate(
    stillbf_1 = 0,
    daysbf_1 = as.numeric(dateprgtrm_2 - dateprgtrm_1)
  ) %>% dplyr::select(uncchdid, contains("_1"), contains("_2")) %>%
  pivot_longer(-"uncchdid", cols_vary = "slowest", names_to = c(".value", "pregord"), names_pattern = "(.*)_(.)") %>%
  rename_with(~paste0(., "2"), c("stillbf", "daysbf", "dateprgtrm")) %>% mutate_at(vars(pregord), as.numeric)

preg2 <- preg1 %>% left_join(fixed.bf, by = c('uncchdid', 'pregord')) %>% mutate(
  stillbf = case_when(stillbf2 == 0 ~ 0, TRUE ~ stillbf),
  daysbf = case_when(daysbf2 > 0 ~ daysbf2, TRUE ~ daysbf)
  ) %>% dplyr::select(-stillbf2, -daysbf2, -dateprgtrm2)

# different kinds of post-partum to test
# 1. ignoring miscarriage (including all but miscarried and currently pregnant)
post.partum1 <- preg2 %>% filter(pregterm %in% c(1,2,3,5)) %>%
  mutate(
    pp1 = pp.b/30.4,
    pp1.2 = case_when(pp.b <= 12*30.4 ~ pp.b/30.4, pp.b > 12*30.4 ~ 12, TRUE ~ NA)
  ) %>% dplyr::select(uncchdid, pregord, pp1, pp1.2) %>% distinct()

# 2. ignoring miscarriage less than half the pregnancy (4 months or less)
post.partum2 <- preg2 %>% filter(pregnow != 1) %>%
  mutate(
    pp2 = case_when(pregterm %in% c(1,2,3,5) ~ pp.b/30.4,
                    pregterm == 4 & durmonth > 4 ~ pp.b/30.4,
                    pregterm == 4 & durmonth <= 4 ~ NA, TRUE ~ NA),
    pp2.2 = case_when(pp2 <= 12 ~ pp2, pp2 > 12 ~ 12, TRUE ~ NA)
  ) %>% dplyr::select(uncchdid, pregord, pp2, pp2.2) %>% distinct()

# breastfeeding check
# get IDs
ids.bf <- unique(preg[preg$stillbf==1,]$uncchdid)

summary(as.numeric(interview.dates.c[interview.dates.c$uncchdid %in% ids.bf,]$diff)) # we need info from the future
need.2009 <- (interview.dates.c[interview.dates.c$uncchdid %in% ids.bf & interview.dates.c$diff > 0,]$uncchdid) # number for which we need info from the future
length(need.2009)

# get pregnancy info for women still breastfeeing at time of interview in 2005
bf_length <- preg2009[preg2009$uncchdid %in% ids.bf & preg2009$pregterm %in% c(1,2,5,6),] %>%
  # getting info from before/during 2005 survey
  filter(yearprg < 2006) %>% mutate(
    dateprgtrm = as.Date(
      paste0(dayprg,"/", monthprg, "/", yearprg), format = '%d/%m/%Y')
  ) %>% dplyr::select(uncchdid, pregord, pregterm, dateprgtrm, daysbf) %>%
  arrange(uncchdid, desc(dateprgtrm)) %>% distinct(uncchdid, .keep_all = TRUE)

nrow(bf_length[bf_length$uncchdid %in% need.2009,])

# check against 2005 data
preg.check <- preg2[preg2$stillbf==1,] %>% dplyr::select(uncchdid, pregord, pregterm, dateprgtrm, stillbf, daysbf) %>%
  left_join(
    dplyr::select(bf_length, uncchdid, pregord, daysbf2 = daysbf),
    by=c('uncchdid', 'pregord')) %>%
  # check with post-partum
  left_join(post.partum2, by=c('uncchdid', 'pregord')) %>%
  # for the NAs, check the time difference between blood collection and preg hist interview
  left_join(dplyr::select(interview.dates.c, uncchdid, diff), by='uncchdid') %>%
  # convert diff to months for better comparability
  mutate(diff = as.numeric(diff)/30.4)

# combining datasets
# epi2 is epi1 merged with anthro, pregnancy, and wealth
epi2 <- Reduce(function(x, y) left_join(x, y, by = "uncchdid", relationship = "many-to-many"),
               list(sex, chron.age, anthro, dplyr::select(preg2, -sex, -date.blood), wealth2, dplyr::select(df7[df7$survey==2005,], uncchdid, employed, urbindex)), init = epi1) %>%
  left_join(dplyr::select(preg.check, uncchdid, pp2, diff, daysbf2, pregord), by = c('uncchdid', 'pregord')) %>%
  mutate(
    everpreg = case_when(
      !is.na(pregterm) ~ 1,
      sex == 0 ~ 0, TRUE ~ NA
    )
  )

save(epi2, file = getfp("memadeR/epi2.RData", memade = TRUE)) # all other variables

# dataset for imputation
epi.f0 <- epi2[epi2$sex==0 & !is.na(epi2$pregord),]
# selecting columns of interest
epi.f <- epi.f0 %>% dplyr::select(
  cd8t, cd4t, nk, bcell, mono, gran, diff,
  lag.weight, lag.waist, lag.hip, lag.bmic, lag.whr, lag.wthr,
  weight, height, waist, hip, bmic, whr, wthr,
  brgay05, maristat, gradecom, worknow, childsex, wherborn, pregnow, dateprgtrm, prenatal, pregterm, breastfd, pregord,
  pp.b, pc.mean1, pc_2005, employed, urbindex, pp2, daysbf2, choreprg, choreaft, everpreg, employed, pregwork, workaftr
) %>% mutate_all(as.numeric) %>%
  mutate_at(.vars = vars(worknow, prenatal, pregterm, maristat, breastfd, pregwork, workaftr, choreprg, choreaft, pregnow, everpreg, employed), .funs = as_factor) %>% mutate(childsex = factor(childsex, levels=c(-9, 1, 2)))

#####
# imputation
imp <- mice::mice(data = epi.f, m = 20, method = 'rf', seed = 500)

# check which one deviates the least from the imputed mean
final.imp <- imp$imp$daysbf2 %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans()
final.imp2 <- abs(final.imp - mean(final.imp))

# completing the data using the iteration closest to the mean
epi.f.completed <- mice::complete(imp, which.min(final.imp2)) %>% cbind("uncchdid" = epi.f0$uncchdid)

#####
# exporting completed dataset as RData and Stata file
save(epi.f.completed, file = getfp('memadeR/epi.f.completed.RData', memade = TRUE))

# modifying file for Stata & exporting a Stata file version
epi.f.completeds <- epi.f.completed %>% rename_with(~ gsub('\\.', '_', .x))
write_dta(epi.f.completeds, getfp('memadeStata/epi.f.completed.dta', memade = TRUE))

bf.imputed <- preg2 %>% dplyr::select(uncchdid, pregord, pregterm, dateprgtrm, daysbf, stillbf, pregnow) %>%
  left_join(
    dplyr::select(epi.f.completed, uncchdid, pp2, daysbf2, diff, pregord), by = c('uncchdid', 'pregord')
  ) %>%
  mutate(pp0 = pp2-diff) %>%
  # months breastfeeding variable
  mutate(mo.bf = case_when(
    daysbf >= 0 & daysbf/30.4 > pp2 ~ pp2,
    daysbf >= 0 ~ daysbf/30.4,
    pregterm == 3 | pregterm == 4 | pregterm == 6 ~ 0,
    stillbf == 1 & daysbf2/30.4 >= pp2 ~ pp2,
    stillbf == 1 & daysbf2/30.4 < pp2 ~ daysbf2/30.4,
    daysbf==-9 & stillbf==-9 ~ 0,
    TRUE ~ NA
  )) %>% mutate(mo.bf2 = ifelse(mo.bf > 12, 12, mo.bf)) %>% distinct(uncchdid, pregord, mo.bf, .keep_all = TRUE) %>%
  # months post-partum
  group_by(uncchdid) %>% mutate(mo.pp = min(12, pp2)) %>% ungroup()

total.bf <- bf.imputed %>% group_by(uncchdid) %>% summarise(
  # creating number of pregnancies (includes miscarriage)
  num.preg.total = n(),
  # creating number completed pregnancies
  num.preg.comp = sum(pregterm %in% c(1,2,3,5), na.rm=TRUE),
  # creating number of living kids
  num.kids = sum(pregterm==1, na.rm=TRUE) + 2*sum(pregterm==5, na.rm=TRUE),
  # creating totmonthsbf_2005
  totmonthsbf = sum(mo.bf, na.rm=TRUE)
) %>% mutate(
  # categorical BF variable
  ctotmonthsbf = case_when(totmonthsbf == 0 ~ 0, totmonthsbf <= 6 ~ 1, TRUE ~ 2),
  # categorical pregnancy variable
  cnum.prg = case_when(num.preg.comp > 1 ~ 2, num.preg.comp == 1 ~ 1, TRUE ~ 0),
  # ever gave birth
  everbirth = ifelse(num.preg.comp > 0, 1, 0),
  everpreg = ifelse(num.preg.total > 0, 1, 0)
)

# bf.imputed2 contains imputed breastfeeding info, number of kids had, and categorical number of pregnancies.
bf.imputed2 <- bf.imputed %>% left_join(total.bf, by = 'uncchdid') %>% arrange(desc(pregord)) %>% arrange(uncchdid) %>%
  distinct(uncchdid, .keep_all = TRUE) %>% mutate_at(
    # for women, replace NAs in relevant variables
    .vars = vars(cnum.prg, totmonthsbf, ctotmonthsbf, everpreg, everbirth, num.preg.total, num.preg.comp, num.kids, pregnow, mo.bf, mo.bf2, mo.pp, pp2),
    .funs = ~ifelse(is.na(.), 0, .)
  )

save(bf.imputed2, file = getfp("memadeR/bf.imputed2.RData", memade = TRUE))
