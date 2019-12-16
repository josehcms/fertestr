#' Measures of fertility from the Malawi 2008 Census
#'
#' Measures of fertility from the Malawi 2008 Census
#'
#' @format
#'   A data frame with 3 variables:
#'   ages (starting age of five-year age group),
#'   P (average parity of women) and
#'   asfr (age-specific period fertility rates)
#'
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
"data.pf_MWI"

#' Parity data from Cambodia 2008 Census
#'
#' Parity data from Cambodia 2008 Census
#'
#' @format
#'   A data frame with 3 variables:
#'   parity (from 0 to N and including NA response codes 98, 99 or NA for example),
#'   ages (starting age group interval)
#'   women (women counts of respective parity from respective age group)
#'
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
"data.prty_KHM"

#' Parity data from Kenya 1989 Census
#'
#' Parity data from Kenya 1989 Census
#'
#' @format
#'   A data frame with 3 variables:
#'   parity (from 0 to N and including NA response codes 98, 99 or NA for example),
#'   ages (starting age-group interval)
#'   women (women counts of respective parity from respective age group)
#'
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
"data.prty_KEN"

#' Original Brass PF multipliers for linear interpolation
#'
#' Original Brass PF multipliers for linear interpolation
#'
#' @format
#'   A data frame with age coefficients for linear interpolation of cumulate fertility schedules based
#'   on mean age at childbearing (m), parity ratio (P1/P2) and age-specific fertility rates
#'   ratio (f1/f2).
#'
#' @source
#' Brass W. 1975. Methods for Estimating Fertility and Mortality from Limited and Defected Data.
#' North Carolina: Carolina Population Center.
"mult.pf_brass"

#' Coale and Trussel variant for Brass PF multipliers for fertility rates computed
#' from births in a 12-month period by age of mother at end of period
#'
#' Coale and Trussel variant for Brass PF multipliers for fertility rates computed
#' from births in a 12-month period by age of mother at end of period
#'
#' @format
#'   A data frame with age coefficients (a, b and c) for linear interpolation of cumulate fertility
#'   schedules
#'
#' @source
#' United Nations. 1983. Manual X: Indirect techniques for demographic estimation
#' (United Nations publication, Sales No. E.83.XIII.2).
"mult.cltrs_shift"


#' Coale and Trussel variant for Brass PF multipliers for fertility rates computed
#' from births by age of mother at delivery
#'
#' Coale and Trussel variant for Brass PF multipliers for fertility rates computed
#' from births by age of mother at delivery
#'
#' @format
#'   A data frame with age coefficients (a, b and c) for linear interpolation of cumulate fertility
#'   schedules
#'
#' @source
#' United Nations. 1983. Manual X: Indirect techniques for demographic estimation
#' (United Nations publication, Sales No. E.83.XIII.2).
"mult.cltrs_noshift"

