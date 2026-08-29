
# import data ---------------

load("../data.RData")
load("../sensitivity_data.RData")

# merge ---------------

epi_new3 <- new_data %>% left_join(
  epi_deid2 %>% dplyr::select(
    new_ids, cd8t, cd4t, nk, bcell, mono, gran, sex,
    waist, bmi, whtr, z.bmi, z.waist, z.whtr,
    uw, cursmok, hai, totmonthsbf, num.preg.comp
  )
)

# residualised correlations ---------------

## correlations ---------------

cor.table.ic <- cbind(
  "Clock" = c(
    "Horvath",
    "Hannum",
    "PhenoAge",
    "GrimAge",
    "DNAmTL",
    "DunedinPACE",
    "PCHorvath2",
    "PCHannum",
    "PCPhenoAge",
    "PCGrimAge",
    "PCDNAmTL"
  ),
  "Correlation" = c(
    # non-PC
    sprintf("%.4f", cor(epi_new3$rs.horvath_f, epi_new3$rs.horvath_nf)),
    sprintf("%.4f", cor(epi_new3$rs.hannum_f, epi_new3$rs.hannum_nf)),
    sprintf("%.4f", cor(epi_new3$rs.phenoage_nf, epi_new3$rs.phenoage_f)),
    sprintf("%.4f", cor(epi_new3$rs.grimage_nf, epi_new3$rs.grimage_f, use = "pairwise.complete")),
    sprintf("%.4f", cor(epi_new3$rs.dnamtl_nf, epi_new3$rs.dnamtl_f)),
    sprintf("%.4f", cor(epi_new3$rs.dunedinpace_nf, epi_new3$rs.dunedinpace_f)),
    # PC
    sprintf("%.4f", cor(epi_new3$rs.pchorvath2_f, epi_new3$rs.pchorvath2_nf)),
    sprintf("%.4f", cor(epi_new3$rs.pchannum_f, epi_new3$rs.pchannum_nf)),
    sprintf("%.4f", cor(epi_new3$rs.pcphenoage_nf, epi_new3$rs.pcphenoage_f)),
    sprintf("%.4f", cor(epi_new3$rs.pcgrimage_nf, epi_new3$rs.pcgrimage_f, use = "pairwise.complete")),
    sprintf("%.4f", cor(epi_new3$rs.pcdnamtl_nf, epi_new3$rs.pcdnamtl_f))
  )
) %>%
  # publication-ready table
  kable(
    format='html',
    escape=FALSE
  ) %>%
    kable_styling(bootstrap_options = "striped", full_width = F) %>%
    column_spec(column = 1, bold = TRUE) %>%
    collapse_rows(columns = 1, valign = "top") %>%
    add_header_above(
      header = c(
        setNames(
          object = 2,
          nm = "Correlations between clocks with SNPs and \ncross-hybridizing probes kept vs filtered out"
        )
      )
    )

# recalculating the latent biological age variable -----------------
library(psych)

## new clocks - unfiltered ----------------

unfiltered <- c(
  # non-PC
  "horvath_nf",
  "hannum_nf",
  "phenoage_nf",
  "grimage_nf",
  "dnamtl_nf",
  "dunedinpace_nf",
  # PC
  "pchorvath1_nf",
  "pchorvath2_nf",
  "pchannum_nf",
  "pcphenoage_nf",
  "pcgrimage_nf",
  'pcdnamtl_nf'
)
filtered <- c(
  # non-PC
  "horvath_f",
  'hannum_f',
  'phenoage_f',
  'grimage_f',
  "dnamtl_f",
  "dunedinpace_f",
  # PC
  "pchorvath1_f",
  "pchorvath2_f",
  "pchannum_f",
  "pcphenoage_f",
  "pcgrimage_f",
  "pcdnamtl_f"
)

### non-PC ----------------

cor_nf <- cor(epi_new3 %>% dplyr::select(paste0("rs.", unfiltered)) %>% dplyr::select(-starts_with("rs.pc")))

cor_nf %>% KMO()
cor_nf %>% eigen()
cor_nf %>% fa.parallel(fm = "ml", fa = "fa", n.obs = 1650)

fac_nf <- fa(cor_nf, nfactors = 1, n.obs = 1650)

### PC ----------------

cor_nf_pc <- cor(epi_new3 %>% dplyr::select(
  "rs.pchorvath2_nf",
  "rs.pchannum_nf",
  "rs.pcphenoage_nf",
  "rs.pcgrimage_nf",
  'rs.pcdnamtl_nf',
  "rs.dunedinpace_nf"
))

cor_nf_pc %>% KMO()
cor_nf_pc %>% eigen()
cor_nf_pc %>% fa.parallel(fm = "ml", fa = "fa", n.obs = 1650)

fac_nf_pc <- fa(cor_nf_pc, nfactors = 1, n.obs = 1650)

## new clocks - filtered

### non-PC ----------------

cor_f <- cor(epi_new3 %>% dplyr::select(paste0("rs.", filtered)) %>% dplyr::select(-starts_with("rs.pc")))

cor_f %>% KMO()
cor_f %>% eigen()
cor_f %>% fa.parallel(fm = "ml", fa = "fa", n.obs = 1650)

fac_f <- fa(cor_f, nfactors = 1, n.obs = 1650)

### PC ----------------

cor_f_pc <- cor(epi_new3 %>% dplyr::select(
  "rs.pchorvath2_f",
  "rs.pchannum_f",
  "rs.pcphenoage_f",
  "rs.pcgrimage_f",
  'rs.pcdnamtl_f',
  "rs.dunedinpace_f"
))

cor_f_pc %>% KMO()
cor_f_pc %>% eigen()
cor_f_pc %>% fa.parallel(fm = "ml", fa = "fa", n.obs = 1650)

fac_f_pc <- fa(cor_f_pc, nfactors = 1, n.obs = 1650)

## making the scores ------------------------

epi_new3[["factor_nf"]] <- factor.scores(
  epi_new3 %>% dplyr::select(paste0("rs.", unfiltered)) %>% dplyr::select(-starts_with("rs.pc")),
  fac_nf$loadings, method = "Bartlett")$scores[,1]

epi_new3[["factor_nf_pc"]] <- factor.scores(
  epi_new3 %>% dplyr::select(
    "rs.pchorvath2_nf",
    "rs.pchannum_nf",
    "rs.pcphenoage_nf",
    "rs.pcgrimage_nf",
    'rs.pcdnamtl_nf',
    "rs.dunedinpace_nf"
  ),
  fac_nf_pc$loadings, method = "Bartlett")$scores[,1]

epi_new3[["factor_f"]] <- factor.scores(
  epi_new3 %>% dplyr::select(paste0("rs.", filtered)) %>% dplyr::select(-starts_with("rs.pc")),
  fac_f$loadings, method = "Bartlett")$scores[,1]

epi_new3[["factor_f_pc"]] <- factor.scores(
  epi_new3 %>% dplyr::select(
    "rs.pchorvath2_f",
    "rs.pchannum_f",
    "rs.pcphenoage_f",
    "rs.pcgrimage_f",
    'rs.pcdnamtl_f',
    "rs.dunedinpace_f"
  ),
  fac_f_pc$loadings, method = "Bartlett")$scores[,1]

# correlations among latent variables ---------------

cor.table.fac <- cbind(
  "Factor" = c(
    "FactorAge",
    "PCFactorAge"
  ),
  "Correlation" = c(
    # non-PC
    sprintf("%.4f", cor(epi_new3$factor_nf, epi_new3$factor_f)),
    sprintf("%.4f", cor(epi_new3$factor_nf_pc, epi_new3$factor_f_pc))
  )
) %>%
  # publication-ready table
  kable(
    format='html',
    escape=FALSE
  ) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  column_spec(column = 1, bold = TRUE) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  add_header_above(
    header = c(
      setNames(
        object = 2,
        nm = "Correlations between factor scores made with clocks with SNPs and \ncross-hybridizing probes kept vs filtered out"
      )
    )
  )

# test factor models ------------------

## comparing non-PC factor scores ------------------

### BMI ------------------

m.bmi_nf <- lm(
  reformulate(termlabels = c("z.bmi", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_nf"),
  data = epi_new3
)

m.bmi_f <- lm(
  reformulate(termlabels = c("z.bmi", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_f"),
  data = epi_new3
)

### WC ------------------

m.waist_nf <- lm(
  reformulate(termlabels = c("z.waist", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_nf"),
  data = epi_new3
)

m.waist_f <- lm(
  reformulate(termlabels = c("z.waist", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_f"),
  data = epi_new3
)

### WHtR ------------------

m.whtr_nf <- lm(
  reformulate(termlabels = c("z.whtr", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_nf"),
  data = epi_new3
)

m.whtr_f <- lm(
  reformulate(termlabels = c("z.whtr", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_f"),
  data = epi_new3
)


## comparing PC factor scores ------------------

### BMI ------------------

m.bmi_nf_pc <- lm(
  reformulate(termlabels = c("z.bmi", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_nf_pc"),
  data = epi_new3
)

m.bmi_f_pc <- lm(
  reformulate(termlabels = c("z.bmi", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_f_pc"),
  data = epi_new3
)

### WC ------------------

m.waist_nf_pc <- lm(
  reformulate(termlabels = c("z.waist", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_nf_pc"),
  data = epi_new3
)

m.waist_f_pc <- lm(
  reformulate(termlabels = c("z.waist", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_f_pc"),
  data = epi_new3
)

### WHtR ------------------

m.whtr_nf_pc <- lm(
  reformulate(termlabels = c("z.whtr", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_nf_pc"),
  data = epi_new3
)

m.whtr_f_pc <- lm(
  reformulate(termlabels = c("z.whtr", "sex", "cd8t", "cd4t", "nk", "bcell", "mono", "gran", "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"),
              response = "factor_f_pc"),
  data = epi_new3
)


## make models function ------------------

makemodel <- function(
    response,
    exposure,
    termlabels = c(
      "sex",
      "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
      "cursmok", "hai", "num.preg.comp*totmonthsbf", "uw"
    ),
    data = epi_new3
) {
  model <- lm(
    reformulate(termlabels = c(exposure, termlabels),
                response = response),
    data = data
  )
  return(model)
}

termlist <- c("Unfiltered", "Filtered")

## non-PC clocks ------------------

### BMI ------------------

m.bmi.horvath <- list(
  makemodel("rs.horvath_f", "z.bmi"),
  makemodel("rs.horvath_nf", "z.bmi")
)

m.bmi.hannum <- list(
  makemodel("rs.hannum_f", "z.bmi"),
  makemodel("rs.hannum_nf", "z.bmi")
)

m.bmi.phenoage <- list(
  makemodel("rs.phenoage_f", "z.bmi"),
  makemodel("rs.phenoage_nf", "z.bmi")
)

m.bmi.grimage <- list(
  makemodel("rs.grimage_f", "z.bmi"),
  makemodel("rs.grimage_nf", "z.bmi")
)

m.bmi.dnamtl <- list(
  makemodel("rs.dnamtl_f", "z.bmi"),
  makemodel("rs.dnamtl_nf", "z.bmi")
)

m.bmi.dpace.f <- list(
  makemodel("rs.dunedinpace_f", "z.bmi*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai", "num.preg.comp*totmonthsbf"
  ), data = epi_new3 %>% filter(sex == 0)),
  makemodel("rs.dunedinpace_nf", "z.bmi*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai", "num.preg.comp*totmonthsbf"
  ), data = epi_new3 %>% filter(sex == 0))
)

m.bmi.dpace.m <- list(
  makemodel("rs.dunedinpace_f", "z.bmi*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai"
  ), data = epi_new3 %>% filter(sex == 1)),
  makemodel("rs.dunedinpace_nf", "z.bmi*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai"
  ), data = epi_new3 %>% filter(sex == 1))
)

### WC ------------------

m.waist.horvath <- list(
  makemodel("rs.horvath_f", "z.waist"),
  makemodel("rs.horvath_nf", "z.waist")
)

m.waist.hannum <- list(
  makemodel("rs.hannum_f", "z.waist"),
  makemodel("rs.hannum_nf", "z.waist")
)

m.waist.phenoage <- list(
  makemodel("rs.phenoage_f", "z.waist"),
  makemodel("rs.phenoage_nf", "z.waist")
)

m.waist.grimage <- list(
  makemodel("rs.grimage_f", "z.waist"),
  makemodel("rs.grimage_nf", "z.waist")
)

m.waist.dnamtl <- list(
  makemodel("rs.dnamtl_f", "z.waist"),
  makemodel("rs.dnamtl_nf", "z.waist")
)

m.waist.dpace.f <- list(
  makemodel("rs.dunedinpace_f", "z.waist*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai", "num.preg.comp*totmonthsbf"
  ), data = epi_new3 %>% filter(sex == 0)),
  makemodel("rs.dunedinpace_nf", "z.waist*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai", "num.preg.comp*totmonthsbf"
  ), data = epi_new3 %>% filter(sex == 0))
)

m.waist.dpace.m <- list(
  makemodel("rs.dunedinpace_f", "z.waist*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai"
  ), data = epi_new3 %>% filter(sex == 1)),
  makemodel("rs.dunedinpace_nf", "z.waist*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai"
  ), data = epi_new3 %>% filter(sex == 1))
)

### WHtR ------------------

m.whtr.horvath <- list(
  makemodel("rs.horvath_f", "z.whtr"),
  makemodel("rs.horvath_nf", "z.whtr")
)

m.whtr.hannum <- list(
  makemodel("rs.hannum_f", "z.whtr"),
  makemodel("rs.hannum_nf", "z.whtr")
)

m.whtr.phenoage <- list(
  makemodel("rs.phenoage_f", "z.whtr"),
  makemodel("rs.phenoage_nf", "z.whtr")
)

m.whtr.grimage <- list(
  makemodel("rs.grimage_f", "z.whtr"),
  makemodel("rs.grimage_nf", "z.whtr")
)

m.whtr.dnamtl <- list(
  makemodel("rs.dnamtl_f", "z.whtr"),
  makemodel("rs.dnamtl_nf", "z.whtr")
)

m.whtr.dpace.f <- list(
  makemodel("rs.dunedinpace_f", "z.whtr*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai", "num.preg.comp*totmonthsbf"
  ), data = epi_new3 %>% filter(sex == 0)),
  makemodel("rs.dunedinpace_nf", "z.whtr*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai", "num.preg.comp*totmonthsbf"
  ), data = epi_new3 %>% filter(sex == 0))
)

m.whtr.dpace.m <- list(
  makemodel("rs.dunedinpace_f", "z.whtr*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai"
  ), data = epi_new3 %>% filter(sex == 1)),
  makemodel("rs.dunedinpace_nf", "z.whtr*uw", termlabels = c(
    "cd8t", "cd4t", "nk", "bcell", "mono", "gran",
    "cursmok", "hai"
  ), data = epi_new3 %>% filter(sex == 1))
)

## PC clocks ------------------

### BMI ------------------

m.bmi.horvath_pc <- list(
  makemodel("rs.pchorvath2_f", "z.bmi"),
  makemodel("rs.pchorvath2_nf", "z.bmi")
)

m.bmi.hannum_pc <- list(
  makemodel("rs.pchannum_f", "z.bmi"),
  makemodel("rs.pchannum_nf", "z.bmi")
)

m.bmi.phenoage_pc <- list(
  makemodel("rs.pcphenoage_f", "z.bmi"),
  makemodel("rs.pcphenoage_nf", "z.bmi")
)

m.bmi.grimage_pc <- list(
  makemodel("rs.pcgrimage_f", "z.bmi"),
  makemodel("rs.pcgrimage_nf", "z.bmi")
)

m.bmi.dnamtl_pc <- list(
  makemodel("rs.pcdnamtl_f", "z.bmi"),
  makemodel("rs.pcdnamtl_nf", "z.bmi")
)

### WC ------------------

m.waist.horvath_pc <- list(
  makemodel("rs.pchorvath2_f", "z.waist"),
  makemodel("rs.pchorvath2_nf", "z.waist")
)

m.waist.hannum_pc <- list(
  makemodel("rs.pchannum_f", "z.waist"),
  makemodel("rs.pchannum_nf", "z.waist")
)

m.waist.phenoage_pc <- list(
  makemodel("rs.pcphenoage_f", "z.waist"),
  makemodel("rs.pcphenoage_nf", "z.waist")
)

m.waist.grimage_pc <- list(
  makemodel("rs.pcgrimage_f", "z.waist"),
  makemodel("rs.pcgrimage_nf", "z.waist")
)

m.waist.dnamtl_pc <- list(
  makemodel("rs.pcdnamtl_f", "z.waist"),
  makemodel("rs.pcdnamtl_nf", "z.waist")
)

### WHtR ------------------

m.whtr.horvath_pc <- list(
  makemodel("rs.pchorvath2_f", "z.whtr"),
  makemodel("rs.pchorvath2_nf", "z.whtr")
)

m.whtr.hannum_pc <- list(
  makemodel("rs.pchannum_f", "z.whtr"),
  makemodel("rs.pchannum_nf", "z.whtr")
)

m.whtr.phenoage_pc <- list(
  makemodel("rs.pcphenoage_f", "z.whtr"),
  makemodel("rs.pcphenoage_nf", "z.whtr")
)

m.whtr.grimage_pc <- list(
  makemodel("rs.pcgrimage_f", "z.whtr"),
  makemodel("rs.pcgrimage_nf", "z.whtr")
)

m.whtr.dnamtl_pc <- list(
  makemodel("rs.pcdnamtl_f", "z.whtr"),
  makemodel("rs.pcdnamtl_nf", "z.whtr")
)


# coefficient comparison ------------------

## function to report compared betas ------------------

table_coeff_comp <- function(coeff_comp) {
  t <- paste0(
    sprintf("%.3f", coeff_comp[1]),
    " (",
    sprintf("%.3f", coeff_comp[3]),
    ", ",
    sprintf("%.3f", coeff_comp[4]),
    ")"
  )
  return(t)
}

## differences in factor score betas ------------------

coeff_diffs_facs <- cbind(
  "Factor Score" = c(
    "FactorAge",
    "PCFactorAge"
  ),
  "BMI - Difference (95% CI)" = c(
    table_coeff_comp(coeff_comp(m.bmi_nf, m.bmi_f, "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi_nf_pc, m.bmi_f_pc, "z.bmi", nested = FALSE))
  ),
  "WC - Difference (95% CI)" = c(
    table_coeff_comp(coeff_comp(m.waist_nf, m.waist_f, "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist_nf_pc, m.waist_f_pc, "z.waist", nested = FALSE))
  ),
  "WHtR - Difference (95% CI)" = c(
    table_coeff_comp(coeff_comp(m.whtr_nf, m.whtr_f, "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr_nf_pc, m.whtr_f_pc, "z.whtr", nested = FALSE))
  )
) %>%
  # publication-ready table
  kable(
    format='html',
    escape=FALSE,
    caption = "Adjusted for white blood cell proportions, smoking status, underweight status, HAI, sex, parity and breastfeeding"
  ) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  column_spec(column = 1, bold = TRUE) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  add_header_above(header = c(setNames(object = 4, nm = "Coefficient Differences for Factor Scores")))

## differences in individual clock betas ------------------

coeff_diffs_ics <- cbind(
  "Clock" = c(
    "Horvath",
    "Hannum",
    "PhenoAge",
    "GrimAge",
    "DNAmTL",
    "DunedinPACE - Females",
    "DunedinPACE - Males",
    "PCHorvath2",
    "PCHannum",
    "PCPhenoAge",
    "PCGrimAge",
    "PCDNAmTL"
  ),
  "BMI - Difference (95% CI)" = c(
    table_coeff_comp(coeff_comp(m.bmi.horvath[[1]], m.bmi.horvath[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.hannum[[1]], m.bmi.hannum[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.phenoage[[1]], m.bmi.phenoage[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.grimage[[1]], m.bmi.grimage[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.dnamtl[[1]], m.bmi.dnamtl[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.dpace.f[[1]], m.bmi.dpace.f[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.dpace.m[[1]], m.bmi.dpace.m[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.horvath_pc[[1]], m.bmi.horvath_pc[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.hannum_pc[[1]], m.bmi.hannum_pc[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.phenoage_pc[[1]], m.bmi.phenoage_pc[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.grimage_pc[[1]], m.bmi.grimage_pc[[2]], "z.bmi", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.bmi.dnamtl_pc[[1]], m.bmi.dnamtl_pc[[2]], "z.bmi", nested = FALSE))
  ),
  "WC - Difference (95% CI)" = c(
    table_coeff_comp(coeff_comp(m.waist.horvath[[1]], m.waist.horvath[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.hannum[[1]], m.waist.hannum[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.phenoage[[1]], m.waist.phenoage[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.grimage[[1]], m.waist.grimage[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.dnamtl[[1]], m.waist.dnamtl[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.dpace.f[[1]], m.waist.dpace.f[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.dpace.m[[1]], m.waist.dpace.m[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.horvath_pc[[1]], m.waist.horvath_pc[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.hannum_pc[[1]], m.waist.hannum_pc[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.phenoage_pc[[1]], m.waist.phenoage_pc[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.grimage_pc[[1]], m.waist.grimage_pc[[2]], "z.waist", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.waist.dnamtl_pc[[1]], m.waist.dnamtl_pc[[2]], "z.waist", nested = FALSE))
  ),
  "WHtR - Difference (95% CI)" = c(
    table_coeff_comp(coeff_comp(m.whtr.horvath[[1]], m.whtr.horvath[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.hannum[[1]], m.whtr.hannum[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.phenoage[[1]], m.whtr.phenoage[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.grimage[[1]], m.whtr.grimage[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.dnamtl[[1]], m.whtr.dnamtl[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.dpace.f[[1]], m.whtr.dpace.f[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.dpace.m[[1]], m.whtr.dpace.m[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.horvath_pc[[1]], m.whtr.horvath_pc[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.hannum_pc[[1]], m.whtr.hannum_pc[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.phenoage_pc[[1]], m.whtr.phenoage_pc[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.grimage_pc[[1]], m.whtr.grimage_pc[[2]], "z.whtr", nested = FALSE)),
    table_coeff_comp(coeff_comp(m.whtr.dnamtl_pc[[1]], m.whtr.dnamtl_pc[[2]], "z.whtr", nested = FALSE))
  )
) %>%
  # publication-ready table
  kable(
    format='html',
    escape=FALSE,
    caption = "Adjusted for white blood cell proportions, smoking status, underweight status, HAI, sex, parity and breastfeeding"
  ) %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  column_spec(column = 1, bold = TRUE) %>%
  collapse_rows(columns = 1, valign = "top") %>%
  add_header_above(header = c(setNames(object = 4, nm = "Coefficient Differences for Individual Clocks")))

