\name{snp.leg}
\alias{psnp.leg}
\alias{dsnp.leg}
\title{snp.leg}
\description{ snp.leg implement seminonparametric representation of density and distribution function defined on the unit interval}
\arguments{   
   \item{u}{ scalar or vector from the unit interval}
   \item{delta}{parametrization}
}
\value{ Function returns a vector pdf.cdf values }
\author{Konrad Smolinski}
\examples{
n <- 200
delta <- c(1,2,2-3,4) 
u <- seq(from=0.001,to=.999,len=n)
delta <- as.matrix(delta)

pdf.leg <- dsnp.leg(u,delta)
cdf.leg <- psnp.leg(u,delta)

par(mfrow=c(2,1))
plot(pdf.leg, type='l',col='2')
plot(cdf.leg, type='l',col='2')
}   
\keyword{}