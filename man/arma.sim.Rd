\name{arma.sim}
\alias{arma.sim}
\title{Generate Random Time Series Sequence.}
\description{
This function is to generate a time series sequence.
}
\usage{
arma.sim(ars, mas, nobs, c = 0, sigma = 1, seed = NA)
}
\arguments{
\item{ars}{AR part.}
\item{mas}{MA part.}
\item{nobs}{Length of time series sequence.}
\item{c}{Constant.}
\item{sigma}{The standrad error.}
\item{seed}{Set seed.}
}
\value{The arma.sim method returns an object of class “Time Series”.}
\examples{
ts <- arma.sim(ars = c(-0.9, -1.4, -0.7, -0.6), mas = c(0.5, -0.4), nobs = 100)
}
