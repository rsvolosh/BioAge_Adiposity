devtools::load_all()

# import data pre-impute
anthro02 <- import.data('/Users/rvolosh/Documents/Papers/Cebu Epigenetics/anthro et al/anthdiet2002.dta', c(lag.weight = weight, lag.height = height, lag.waist = waist, lag.hip = hip, -sexic, -uncmomid))
anthro05 <- import.data('/Users/rvolosh/Documents/Papers/Cebu Epigenetics/anthro et al/anthdiet.dta', c(weight, height, waist, hip, -basehhno, -brgay05, -hhnumb05, -woman05, -momch05, -linenumb))

colnames(anthro02)[colnames(anthro02) %in% colnames(anthro05) & !grepl("uncchdid", colnames(anthro02))] <- sub("$", "02", colnames(anthro02)[colnames(anthro02) %in% colnames(anthro05) & !grepl("uncchdid", colnames(anthro02))])

anthro0 = full_join(
  anthro02, anthro05,
  by = "uncchdid"
) %>% left_join(
    import.data('/Users/rvolosh/Documents/Papers/Cebu Epigenetics/anthro et al/PREGHIST.DTA', pregterm),
    by = "uncchdid"
  ) %>% filter(is.na(pregterm) | pregterm != 6) %>% dplyr::select(-pregterm) %>%
  mutate(across(
    .cols = c(height, waist, hip, weight, lag.waist, lag.height, lag.waist, lag.hip),
    .fns = ~ ifelse(. <= 0, NA, .)
  ))

# fixing typos in the anthro0 file
# reference growth curves to go through individually
growth <- import.data("/Users/rvolosh/Documents/Papers/Cebu Epigenetics/anthro et al/cebuheightlong.dta", c(year, agey, height, weight)) %>% filter(!is.na(uncchdid)) %>% group_by(uncchdid) %>%
  mutate(height = ifelse(
    year > 17 & is.na(height), height[year == 17], ifelse(
      year == 20 & is.na(height) & !is.na(height[year == 16]) & !is.na(height[year == 19]), (height[year == 16] + height[year == 19])/2, height)
  )) %>% ungroup() %>%
  mutate(bmi = weight/(height/100)^2)

# fixing the errors
anthro0[anthro0$uncchdid == 21014,]$lag.height = anthro0[anthro0$uncchdid == 21014,]$height
anthro0[anthro0$uncchdid == 22682,]$lag.height = anthro0[anthro0$uncchdid == 22682,]$height
anthro0[anthro0$uncchdid == 22622,]$lag.height = (anthro0[anthro0$uncchdid == 22622,]$height + growth[growth$uncchdid == 22622 & growth$year == 15,]$height)/2
anthro0[anthro0$uncchdid == 22622,]$lag.weight = (anthro0[anthro0$uncchdid == 22622,]$weight + growth[growth$uncchdid == 22622 & growth$year == 15,]$weight)/2
anthro0[anthro0$uncchdid == 21025,]$lag.height = (anthro0[anthro0$uncchdid == 21025,]$height + growth[growth$uncchdid == 21025 & growth$year == 15,]$height)/2
anthro0[anthro0$uncchdid == 20624,]$lag.height = growth[growth$uncchdid == 20624 & growth$year == 15,]$height

# imputing wealth information for the rolling average
imp <- mice::mice(data = anthro0, m = 5, method = 'rf', seed = 500)

# check which one deviates the most from the imputed mean
final.imp <- c(
  imp$imp$height %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),
  imp$imp$weight %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),
  imp$imp$hip %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),
  imp$imp$waist %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),

  imp$imp$lag.height %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),
  imp$imp$lag.weight %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),
  imp$imp$lag.hip %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans(),
  imp$imp$lag.waist %>% as.data.frame() %>% rename_all(function(x) {paste0('i_', x)}) %>% colMeans()
)
final.imp2 <- c(
  abs(final.imp[1:5] - mean(final.imp[1:5])),
  abs(final.imp[6:10] - mean(final.imp[6:10])),
  abs(final.imp[11:15] - mean(final.imp[11:15])),
  abs(final.imp[16:20] - mean(final.imp[16:20])),

  abs(final.imp[21:25] - mean(final.imp[21:25])),
  abs(final.imp[26:30] - mean(final.imp[26:30])),
  abs(final.imp[31:35] - mean(final.imp[31:35])),
  abs(final.imp[36:40] - mean(final.imp[36:40]))
)

my.min <- c(which.min(final.imp2[1:5]),
            which.min(final.imp2[6:10]),
            which.min(final.imp2[10:15]),
            which.min(final.imp2[16:20]),
            which.min(final.imp2[21:25]),
            which.min(final.imp2[26:30]),
            which.min(final.imp2[31:35]),
            which.min(final.imp2[36:40])) %>%
  median() %>% floor()
# completing the data using the iteration closest to the mean
anthro.completed <- mice::complete(imp, my.min)

anthro <- anthro.completed %>% dplyr::select(
  uncchdid, height, weight, hip, waist, lag.height, lag.weight, lag.hip, lag.waist
) %>% mutate(
  bmic = weight/(height/100)^2,
  whr = waist/hip,
  wthr = waist/height,

  lag.bmic = lag.weight/(lag.height/100)^2,
  lag.whr = lag.waist/lag.hip,
  lag.wthr = lag.waist/lag.height,
) %>% mutate(
  d.bmi = bmic - lag.bmic,
  d.waist = waist - lag.waist,
  d.whr = whr - lag.whr,
  d.wthr = wthr - lag.wthr,
  uw = ifelse(bmic < 18.5, 1, 0)
)

save(anthro, file = getfp("memadeR/anthro.RData", memade = TRUE))
