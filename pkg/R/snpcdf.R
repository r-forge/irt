etFx <- function(ubasis,tvec,scheme,udex,...){

#this is etF() in context of repeated calls with same setup
#thus u, ubasis etc is in environment 

#simplified version of etF(). First stab at efficiency  16 July 2008

#
# We are asked to give etF at points u and t=tvec.
# This does it by quadrature to get M(t) (the normalizing
# constant) and the value of "numerator" integral at each
# point u with a single piecewise integration.
# Each piece in the "scheme" has nq points of quadrature
# between the points in u.

# A scheme INCLUDES an interval that ends at 1.0
# See pwts1()


#nutot <- nrow(scheme)
#nq <- nutot/nu
#utot <- scheme[,1]
#udex <- (1:nu)*nq
#nt <- length(tvec)
#etu <- rep(0,nrow(scheme))

# this is wasteful since utot^j is generally known in advance
# and it is written in Fortran, which is inefficient

etu <- exp(ubasis%*%tvec)*scheme[,2]


qvec <- cumsum(etu)/sum(etu)
ans <- qvec[udex]
ans <- ans[-length(ans)]  #drop the last one; it is 1 by defn
ans
}
# ---------------------------------------------------------
# input all
#etFx2 <- function(tvec,ubasis,scheme,udex){
etFx2 <- function(tvec,ubasis,weights){
	if(is.null(weights)) etu <- exp(ubasis%*%tvec)
	else etu <- exp(ubasis%*%tvec)*weights
	qvec <- cumsum(etu)/sum(etu)
	ans <- qvec[-length(qvec)]  #drop the last one; it is 1 by defn
ans
}
# ---------------------------------------------------------
ietF <- function(u,t){
# a function  to look at iterative construction of exponential distribution functions
# instead of u was: u.1
	nt <- length(t)
	Fnow <- u
	for(i in 1:nt) Fnow <- (exp(t[i]*Fnow)-1)/(exp(t[i])-1)      
return(Fnow)
}
# ---------------------------------------------------------
snpF <- function(u,t,Fmethod){
	# instead of u was: u.1
# a function  to look at iterative construction of exponential distribution functions
if(Fmethod=="snp") Fnow <- psnp.leg(u=u,delta=t)
#if(Fmethod=="snp.cheb") Fnow <- psnp.cheb(u=u,delta=t)
#if(Fmethod=="snp.sin") Fnow <- psnp.sine(u=u,delta=t)
return(Fnow)
}
# ---------------------------------------------------------
mixET <- function(u,t){
# instead of u was: ubasis
# a function  to look at iterative construction of exponential distribution functions
	nt <- length(t)
	Fnow <- matrix(0,nr=nrow(u),nc=1)
	for(i in 1:nt)  Fnow <- Fnow+(exp(t[i]*u[,i])-1)/(exp(t[i])-1)
	#print(Fnow/nt);stop(44)
return(as.vector(Fnow/nt)) 
}
# ---------------------------------------------------------
logisticF <- function(u,t){
# instead of u was: u.log	
	beta <- exp(t)
	expval <- exp(t*u)
	ans <- expval/(1+expval)
	ans
}

#-----------------------------------------------------------
betaFx <- function(u,tvec){
# instead of u was: u.1	
	ans <- pbeta(u,shape1=exp(tvec[1]),shape2=exp(tvec[2]))
return(ans)
}
#-----------------------------------------------------------
