% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reverse_survival_aux.R
\name{SingleAgeLogQuad}
\alias{SingleAgeLogQuad}
\title{Estimate single-age survival functions using Log-Quadratic Model}
\usage{
SingleAgeLogQuad(
  k = NULL,
  e0 = NULL,
  q0_1 = NULL,
  q0_5 = NULL,
  q15_35 = NULL,
  q15_45 = NULL,
  lt = NULL,
  sex = "total"
)
}
\arguments{
\item{k}{log-quad parameter k}

\item{e0}{log-quad parameter e0}

\item{q0_1}{log-quad parameter 0q1}

\item{q0_5}{log-quad parameter 0q5}

\item{q15_35}{log-quad parameter 15q35}

\item{q15_45}{log-quad parameter 15q45}

\item{lt}{reference life table for modeling log quad parameters, default = HMD}

\item{sex}{sex to retrieve HMD sex and life table ('female','male','total' - default)}
}
\value{
single age life table estimated by ungrouping log-quad model estimation
}
\description{
Estimate single-age survival functions using Log-Quadratic Model and ungroup
}
\examples{
SingleAgeLogQuad( e0 = 70, q0_5 = 0.04 )
}
