\name{arma.plot}
\alias{arma.plot}
\title{Plot the time series data.}
\description{
This function is to plot the time series figure.
}
\usage{
arma.plot(data)
}
\arguments{
\item{data}{Data containing time series.}
}
\value{The arma.plot method returns a figure.}
\examples{
ts <- arma.sim(ars = c(-0.9, -1.4, -0.7, -0.6), mas = c(0.5, -0.4), nobs = 100)
arma.plot(ts)
}
