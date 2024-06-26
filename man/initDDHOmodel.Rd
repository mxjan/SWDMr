% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Exports.R
\name{initDDHOmodel}
\alias{initDDHOmodel}
\title{Initiate DDHO model}
\usage{
initDDHOmodel(object, VarExp, UseDampingRatio = F)
}
\arguments{
\item{object}{A SWDMr object}

\item{VarExp}{The name of the variable to explain by the model}

\item{UseDampingRatio}{The core parameters of the model change for intercept, damping ratio and natural frequency instead of intercept, damping constant and natural frequency}
}
\value{
SWDMr model object
}
\description{
Initiate a driven damped harmonic oscillator model
}
\examples{
model <- initDDHOmodel(swdmr,VarExp = "Arntl")

}
