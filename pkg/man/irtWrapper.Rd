\name{irtWrapper}
\alias{irtWrapper}
\title{irt Wrapper}
\description{ Wrapper for IRT}
\usage{}
\arguments{   
   \item{item}{item}
}
\value{ Returns an object of an irt class}
\author{}
\examples{
# do not run
run <- FALSE
if(run){

load('irt1-v2.Rdata')
nm <- names(EITEMS)
# --------------------------------------
# Two steps:
#	a) construct an item
#	b) run irt() on the item

# Data Manipulation

EITEMS <- EITEMS[1:1500,]
W <- W[1:1500,]

# Single Item Example 1
singleItem <- item(EITEMS[,1],W,itemnm=nm[1])
iSet <- irtSet(iterlim=50,pl=1,plt=TRUE,tag="singleItem")
res <- irt(singleItem,iSet=iSet)

# Single Item Example 2. Formula based
# 	Suppose we hava a data set:
dta <- data.frame(govguareatsleep=EITEMS[,1],age=W[,1],agesq.01=W[,2],black=W[,3],bornagain=W[,4])

singleItem <- item(govguareatsleep~age+agesq.01+black,~black, data=dta)
iSet <- irtSet(iterlim=50,pl=1,plt=TRUE,tag="singleItem")
res <- irt(singleItem,iSet=iSet)

print(res)
plot(res)
plot(res,type="grid")
# Multiple Item Example
items <- multipleitem(EITEMS[,1:5],W)
iSet <- irtSet(nb=2,iterlim=10,plt=TRUE,tag="mytag")
res <- irt(items,iSet=iSet)

}
}   
\keyword{}