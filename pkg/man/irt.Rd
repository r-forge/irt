\name{irt}
\alias{irt}
\title{irt}
\description{irt}  
\usage{
 
} 
\arguments{   
   \item{ r }{is a vector of ordered outcomes }
   \item{ W }{ is a matrix of characteristics }
   \item{ thstart }{ starting values }
   \item{ beta }{should beta distribution rather than uniform be taken as base density to be tilted}
   \item{snpdeg}{degree}
   \item{type}{type}   
   \item{nqpoints}{nqpoints}
   \item{qrange}{qrange}
   \item{gradtol}{gradtol}
   \item{iterlim}{iterlim}
}
\value{ Function returns an object ...  }
\author{Richard Spady}
\examples{
# -----------------------------------------
# Example 1. 
# Fast by result may not make sense
data(elections)
attach(elections)
Wx <- cbind(age,agesq.01,black,bornagain)
#res <- irt(r=govtakecarewhocant,W=Wx,iterlim=400)
# summary(res)
detach(elections) 

# -----------------------------------------
# Example 2.
# This computation may take a while (10 min or more... )
data(elections)
attach(elections)

Wx <- cbind(age,agesq.01,black,bornagain,blackbornagain,rel.catholic,rel.nonchr, ed.cat1, ed.cat3,ed.cat4,ed.cat5,income.1,income.3,income.4,income.dk,
parent, hispanic,female)
#res <- irt(r=govtakecarewhocant,W=Wx,iterlim=400)
# summary(res)
detach(elections) 
}   
\keyword{} 