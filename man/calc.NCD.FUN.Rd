% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/POCD_functions_master.R
\name{calc.NCD.FUN}
\alias{calc.NCD.FUN}
\title{calculate NCD with one dichotom. parameter}
\usage{
calc.NCD.FUN(Input.all, Control.all, impairment, Timepoint.names,
  CTP.names, CTP.number, NCD.names)
}
\arguments{
\item{Input.all}{list fo data frames containing CTPs /RCIs of patients group}

\item{Control.all}{list fo data frames containing CTPs/RCIs of control group}

\item{impairment}{list of data frames containing the information about the impairment of the patients in logicals (TRUE, FALSE NA)}

\item{Timepoint.names}{names of the test timepoints}

\item{CTP.names}{names of CTP tests}

\item{CTP.number}{numer of used CTP test}

\item{NCD.names}{name of column containing the variables for impairment}
}
\value{
list of data frames comtaining NCD diagnosis for each patient
}
\description{
calculate NCD with one dichotom. parameter
}
