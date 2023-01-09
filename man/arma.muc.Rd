\name{arma.muc}
\alias{arma.muc}
\title{Plot the model uncertainty curve of time series model selection method.}
\usage{
arma.muc(method, ...)
}
\description{
The MUC of a given model selection method with good performance will tend to arch towards the upper left corner.  Moreover, the area under the MUC, marked as AMUC, can be used as a raw measure of uncertainty of uncertainty for the model selection method under examination. The larger the AMUC, the stabler the corresponding model selection method. Overall, we can decide which method has the best performance according to the shape of the MUC and the AMUC.
}
\arguments{
\item{method}{Model selection method like "aic", "bic" and "aicc".}
\item{...}{Others.}
}
\value{The arma.muc method returns a figure, which is the model uncertainty curve of the model selection method.}
\examples{
library(uotm)
ts <- arma.sim(ars = c(-0.9, -1.4, -0.7, -0.6), mas = c(0.5, -0.4), nobs = 100)

\donttest{tmcb1 <- arma.mcb(ts, method = "aic", max.p = 4, max.q = 4,
bsamples = 50, stepwise = TRUE, blength = 12)}
\donttest{arma.muc(c("aic"), tmcb1)}

# compare the uncertainty of different model selection methods using muc
\donttest{tmcb2 <- arma.mcb(ts, method = "bic", max.p = 4, max.q = 4,
bsamples = 50, stepwise = TRUE, blength = 12)}
\donttest{arma.muc(c("aic", "bic"), tmcb1, tmcb2)}
}
