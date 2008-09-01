\name{irtMulti}
\alias{irtMulti}
\title{irtMulti}
\description{This is a multiple single item model}  
\usage{} 
\arguments{   
   \item{ r }{is a vector of ordered outcomes }
   \item{ W }{ is a matrix of characteristics }
}
\value{ Function returns an object of class irtSingle  }
\author{}
\examples{
# -----------------------------------------
# Example 1. 
# Do not run
data(elections)
attach(elections)
Wx <- cbind(age,agesq.01,black,bornagain)
rmulti <- cbind(govtakecarewhocant,govguareatsleep)
#res <- irtSingle(r=rmulti,W=Wx )
detach(elections) 
}   
\keyword{}