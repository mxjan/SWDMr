% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SleepDataManipulation.R
\name{SWdf_AddSD}
\alias{SWdf_AddSD}
\title{Add a sleep deprivation column}
\usage{
SWdf_AddSD(SWdf, Interval)
}
\arguments{
\item{SWdf}{dataframe generated by Read_SW}

\item{Interval}{vector of 2 value for SD interval}
}
\value{
dataframe
}
\description{
Details in description
}
\examples{
SWdf<-AddSD(SWdf,c(24,30))

}
