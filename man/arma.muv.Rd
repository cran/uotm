\name{arma.muv}
\alias{arma.muv}
\title{Estimate the model uncertainty variance of time series model selection method.}
\usage{
arma.muv(data,
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
By calculating the MUV, we can judge the uncertainty of a model selection method and compare the uncertainty of different model selection methods. The lower the degree of uncertainty in the model selection method, the smaller its MUV.
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
\value{The arma.muv method returns an object of class “double”, which is the model uncertainty variance of the model selection method.}
\examples{
library(uotm)
ts <- arma.sim(ars = c(-0.9, -1.4, -0.7, -0.6), mas = c(0.5, -0.4), nobs = 100)
\donttest{tmuv <- arma.muv(ts, max.p = 4, max.q = 4,
stepwise = TRUE, bsamples = 50, blength = 12)}
# the model uncertainty variance of time series model selection method
\donttest{print(tmuv)}
}
