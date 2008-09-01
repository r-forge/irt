# ------------------------------------------------------------------------- 
# ------------------------------------------------------------------------- 
# We are asked to give etF at points u and t=tvec.
# This does it by quadrature to get M(t) (the normalizing
# constant) and the value of "numerator" integral at each
# point u with a single piecewise integration.
# Each piece in the "scheme" has nq points of quadrature
# between the points in u.

# A scheme INCLUDES an interval that ends at 1.0
# See pwts1()
etF <- function(u,tvec,scheme=NULL,nq=8){
	if(is.null(scheme)) scheme <- pwts1(zPoints=u,nqpts=nq)
	nu <- length(u)+1 
	nutot <- nrow(scheme)
	nq <- nutot/nu
	utot <- scheme[,1]
	udex <- (1:nu)*nq
	nt <- length(tvec); js <-1; je<- nt; dubeta <- 1
#	if(beta){ 
#		nt <- length(tvec)-2; js <-3; je <- js+nt-1
#		dubeta <- dbeta(utot,exp(tvec[1]),exp(tvec[2]))
 #  }

	etu <- rep(0,nrow(scheme))
	for(j in js:je){
		if(j==1)etu <- etu+(tvec[j]*(2*utot-1))
		if(j==2)etu <- etu+(tvec[j]*(6*utot^2-6*utot+1))
		if(j==3)etu <- etu+(tvec[j]*(20*utot^3-30*utot^2+12*utot-1))
	}
#	if(beta){ etu <- exp(etu)*dubeta*scheme[,2]}
#	else{etu <- exp(etu)*scheme[,2]}
#	ifelse(beta,etu <- exp(etu)*dubeta*scheme[,2], etu <- exp(etu)*scheme[,2] )
	etu <- exp(etu)*dubeta*scheme[,2]
	qvec <- cumsum(etu)/sum(etu)
	ans <- qvec[udex]
	ans <- ans[-length(ans)]  #drop the last one; it is 1 by defn
ans
}   
#-----------------------------------------------------------               
etF1my <- function(u,tvec,scheme=NULL,nq=8,beta=beta){

	if(is.null(scheme)) scheme <- pwts1(zPoints=u,nqpts=nq)
   nu <- length(u)+1  
	nutot <- nrow(scheme)
#	nq <- nutot/nu
	utot <- scheme[,1]
	udex <- (1:nu)*nq
	nt <- length(tvec); js <-1; je<- nt  

	if(beta){
		nt <- nt-2; js <- js + 2; je <- je + 2
      dubeta <- dbeta(utot,exp(tvec[1]),exp(tvec[2]))
   }
	
	etu <- rep(0,nrow(scheme))

	for(j in js:je){
	   if(j==1)etu <- etu+(tvec[j]*(2*utot-1))
	   if(j==2)etu <- etu+(tvec[j]*(6*utot^2-6*utot+1))
#	   if(j==3)etu <- etu+(tvec[j]*(20*utot^3-30*utot^2+12*utot-1))
	   if(j==3)etu <- etu+(tvec[j]*(20*(utot-1.5)*utot^2+12*utot-1))

   }
	etu <- ifelse(beta,exp(etu)*dubeta*scheme[,2], exp(etu)*scheme[,2])
	qvec <- cumsum(etu)/sum(etu)
	ans <- qvec[udex]
	ans <- ans[-length(ans)]  #drop the last one; it is 1 by defn
return(ans)
}  
# -------------------------------------------------------------------------   