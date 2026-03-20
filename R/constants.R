#' @format epigenetic clocks, non-PC versions, excluding DunedinPACE
#' @export
clock.list = c(
  "dnamage",
  "dnamagehannum",
  "dnamphenoage",
  "dnamgrimage",
  "dnamtl"
)

#' @format epigenetic clocks, PC versions, excluding DunedinPACE
#' @export
pc.clock.list = c(
  # "pchorvath1", # Horvath1 multi-tissue predictor, getting rid of it, as lower KMO and less comparable to other clocks
  "pchorvath2", # Horvath2 skin-and-blood clock
  "pchannum",
  "pcphenoage",
  "pcgrimage",
  "pcdnamtl"
)

#' @format epigenetic clocks, non-PC versions, including DunedinPACE
#' @export
clock.list.full = c(
  "dnamage",
  "dnamagehannum",
  "dnamphenoage",
  "dnamgrimage",
  "dnamtl",
  "dunedin"
)

#' @format epigenetic clocks, PC versions, including DunedinPACE
#' @export
pc.clock.list.full = c(
  # "pchorvath1", # Horvath1 multi-tissue predictor, getting rid of it, as lower KMO and less comparable to other clocks
  "pchorvath2", # Horvath2 skin-and-blood clock
  "pchannum",
  "pcphenoage",
  "pcgrimage",
  "pcdnamtl",
  "dunedin"
)

#' @format epigenetic clock residuals, non-PC versions
#' @export
clock.list.resid = c(
  "ageaccelerationresidual",
  "ageaccelerationresidualhannum",
  "ageaccelpheno",
  "ageaccelgrim",
  "dnamtladjage",
  "dunedin"
)

#' @format white blood cell proportion variables
#' @export
cell.props = c("cd8t", "cd4t", "nk", "bcell", "mono", "gran")

#' @format factor names, non-PC versions
#' @export
factors.nonpc = c("factor.age")
#' @format factor names, PC versions
#' @export
factors.pc = c("pcfactor.age")

# for multivariate regression
# outcomes
#' @format all epigenetic clocks and all of their combinations, non-PC versions
#' @export
mvreg_response_ns_nonpc = c(factors.nonpc, gsub("^", "r.", clock.list), "dunedin") # non-standardised, except DunedinPACE is scaled by ten (see CHANGELOG Jun 24, 2025)
#' @export
mvreg_response_nonpc = c(factors.nonpc, gsub("^", "rs.", clock.list), "rs.dunedin")
#' @format all epigenetic clocks and all of their combinations, PC versions
#' @export
mvreg_response_ns_pc = c(factors.pc, gsub("^", "r.", pc.clock.list), "dunedin") # non-standardised, except DunedinPACE is scaled by ten (see CHANGELOG Jun 24, 2025)
#' @export
mvreg_response_pc = c(factors.pc, gsub("^", "rs.", pc.clock.list), "rs.dunedin")
#' @format labels for all epigenetic clocks and all of their combinations
#' @export
mvreg_response_labels = c("Factor", "Horvath", "Hannum", "PhenoAge", "GrimAge", "DNAmTL", "DunedinPACE")
#' @export
mvreg_response_labels_pub = c("PCFactor", "PCHorvath", "PCHannum", "PCPhenoAge", "PCGrimAge", "PCDNAmTL", "DunedinPACE")
#' @export
mvreg_response_labels_pub2 = c("PCFactor", "PCHorvath (mo)", "PCHannum (mo)", "PCPhenoAge (mo)", "PCGrimAge (mo)", "PCDNAmTL (10*kb)", "DunedinPACE (wk/yr)")

# exposures
#' @format all exposure variables
#' @export
exposures = c("bmi", "waist", "whtr")
#' @format all exposure variables, standardised
#' @export
exposures.std = c("z.bmi", "z.waist", "z.whtr")
#' @format all exposure variables
#' @export
exposures_pub = c("BMI", "WC", "WtHR")
#' @format labels of all exposure variables
#' @export
exposure_labels = c("BMI (kg/m2)", "Waist (cm)", "Waist-to-height ratio")

# covariates for fully-adjusted model
#' @format covariates, expluding female reproductive covariates
#' @export
covars.all = c("cursmok", "hai", "uw")
#' @format covariates, including female reproductive covariates
#' @export
covars.fem = c("cursmok", "hai", "num.preg.comp*totmonthsbf", "uw")

# folder where I store my files
data_folder = "~/anthro et al"
output_data = "~/anthro et al"
