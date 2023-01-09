\name{arma.mcb}
\alias{arma.mcb}
\title{Estimate the model confidence bounds of time series model selection method.}
\usage{
arma.mcb(data,
max.p = 5,
max.q = 5,
method = "aic",
stepwise = TRUE,
blength = 15,
btype = "mbb",
bsamples = 100,
seed = NA)
}
\description{
The MCB algorithm can find a upper bound model (UBM) and a lower bound model (LBM) which are nested and where the true model is included between the two as a specific confidence level. And the width of the MCB is the difference between the number of variables in the UBM and LBM. The model selection method with a smaller width under the same confidence level or a larger confidence level under the same width has less uncertainty.
}
\arguments{
\item{data}{A time series sequence.}
\item{max.p}{The max order of AR part.}
\item{max.q}{The max order of MA part.}
\item{method}{Information criterion to be used in model selection.}
\item{stepwise}{If TRUE, will do stepwise selection (faster). Otherwise, it searches over all models. Non-stepwise selection can be very slow, especially for seasonal models.}
\item{blength}{If btype is "mbb" then blength is the fixed block length used in generating the replicate time series. If btype is "sb" then blength is the mean of the geometric distribution used to generate the block lengths. blength should be a positive integer less than the length of data.}
\item{btype}{The type of simulation required to generate the replicate time series. The possible input values are "mbb" (block resampling with fixed block lengths of blength), "sb" (block resampling with block lengths having a geometric distribution with mean blength)}
\item{bsamples}{A positive integer giving the number of bootstrap replicates required.}
\item{seed}{Seed.}
}
\value{The arma.mcb method returns an object of class “dataframe”, which contains the bootstrap coverage rate and the corresponding model confidence bounds of the user-chosen model selection method in the form of lower confidence bound and upper confidence bound.}
\examples{
library(uotm)
ts <- arma.sim(ars = c(-0.9, -1.4, -0.7, -0.6), mas = c(0.5, -0.4), nobs = 100)
\donttest{tmcb <- arma.mcb(ts, max.p = 4, max.q = 4,
stepwise = TRUE, bsamples = 50, blength = 12)}
# the model confidence bounds of time series model selection method
\donttest{print(tmcb)}
}
