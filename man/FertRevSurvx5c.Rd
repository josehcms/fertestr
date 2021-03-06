% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reverse_survival_clean.R
\name{FertRevSurvx5c}
\alias{FertRevSurvx5c}
\title{Reverse Survival Fertility Estimation Function (x5)}
\source{
Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
  Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
  Population. demographicestimation.iussp.org
}
\usage{
FertRevSurvx5c(
  ages5_c = c(0, 5, 10),
  popx5_c,
  ages5_w = seq(10, 65, 5),
  popx5_w,
  nLx_c,
  l0_c,
  lx5_w,
  asfr = c(0, 0.017, 0.055, 0.057, 0.041, 0.022, 0.007, 0.002),
  asfr_15prior = NULL,
  q15_45f = NULL,
  date_ref
)
}
\arguments{
\item{ages5_c}{children five-year age group ages vector (default 0, 5, 10)}

\item{popx5_c}{children population in five-year age group (x5) matching ages5_c vector}

\item{ages5_w}{women five-year age group ages vector ( default seq( 10, 65, 5 ) )}

\item{popx5_w}{women population in five-year age group format (x5) matching ages_w vector}

\item{nLx_c}{life table person years for ages 0-15 in five year age groups (0L5,5L5,5L10)}

\item{l0_c}{life table radix l0 for children survival estimation}

\item{lx5_w}{women survival function vector in five-year age groups matching ages_w vector}

\item{asfr}{age specific fertility rates for five-year age groups from 10-45 for current period
of estimation}

\item{asfr_15prior}{standardized age specific fertility rates for five-year age groups from 10-45 for the
period of 15 years before the current inquiry period}

\item{q15_45f}{female adult mortality probability for the period of estimation,
period 5 years prior to estimation period, and period 10 years prior to estimation period}

\item{date_ref}{reference date of inquiry given in the following formats:
Y-m-d (4 digit year - 2 digit month - 2 digit day), Y-m (4 digit year - 2 digit month),
Y (4 digit year)}
}
\value{
data.frame with 3 elements:
`year` - reference period of fertility estimation in decimal format;
`TFR`  - estimated total fertility rate;
`births` - estimated total number of births
}
\description{
Reverse Survival Fertility Estimation for five-year age groups of children
}
\examples{

# 1 - User input data
popx5_c <-  c( 1376790, 1472470, 1673110 )

popx5_w <- c(  815930, 780320, 697160, 626430, 361650, 435880,
               393760, 352520, 294280, 230200, 160590, NA )

nLx_c <- c( ( 96270 + 375715 ), 465228, 462406 )

lx5_w <- c( 0.91381, 0.90989, 0.90492, 0.89798, 0.88893, 0.87596,
            0.86029, 0.84188, 0.81791, 0.78472, 0.73735, 0.67316 )

q15_45f <- c( 0.1946, 0.2290, 0.2674 )

asfr <- c( 0.0000, 0.0418,0.1535, 0.1482, 0.1118, 0.0708,
           0.0301, 0.0032 )

asfr_15prior <- c( 0.0000, 0.0533, 0.1974, 0.2144, 0.1836, 0.1332,
                   0.0676, 0.0134 )

FertRevSurvx5c( ages5_c = c( 0, 5, 10 ), popx5_c = popx5_c,
                ages5_w = seq( 10, 65, 5 ), popx5_w = popx5_w,
                nLx_c = nLx_c, lx5_w = lx5_w,
                l0_c = 100000,
                asfr = asfr,
                asfr_15prior = asfr_15prior,
                q15_45f = q15_45f,
                date_ref = '2008-03-03' )


}
