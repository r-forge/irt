# ------------------------------------------------------------------------- 
quadrature <- function(acr,xrange,type="gl"){
	switch(
		gl = GL.YW(acr,xrange),
		sml = quadSmolyak(acr,xrange)
		)	
}

quadSmolyak <- function(acr,xrange=NULL){	
	quad <- smolyak.quad(1,acr)
	pt <- quad$pt; wt <- quad$wt
	if(!is.null(xrange)){
		xL <- diff(xrange)
		wt <- xL*wt
		pt <- xL*(pt-0.5)
	}
	opt <- order(pt)
	pt <- pt[opt]
	wt <- wt[opt]
return(list(pt=pt,wt=wt))
}
# ------------------------------------------------------------------------- 
GL.YW <- function(M,xrange=NULL,epsilon=NULL){
	if(is.null(epsilon))epsilon <- .Machine$double.eps
	if(M%%2==1)stop("M needs to be an even number")
	MM <- (M+1)/2
	Y <- W <- numeric(M)
	for(i in 1:floor(MM)){
		ok <- FALSE
		z <- cos( pi * (i-0.25)/(M+0.5))
			while(ok==FALSE){
				p1 <- 1.0;  p2 <- 0.0
				for(j in 1:M) {
					p3 <- p2;  p2 <- p1
					p1 <- ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
		     }
				pp <- M*(z*p1-p2)/(z^2-1.0)
				z1 <- z
				z <- z-p1/pp
				if(abs(z-z1)<epsilon) ok <- TRUE
			}
			Y[i] <- -z
			Y[M+1-i] <- z
			W[i] <- 2.0/((1-z^2)*pp^2)
			W[M+1-i] <- W[i]
		}
		if(!is.null(xrange)){
			xL <- (xrange[2]-xrange[1])/2.0
			W  <- xL*W
			Y <- xL*Y+(xrange[1]+xrange[2])/2.0

	}
return( cbind(Y,W) )
} 
# -------------------------------------------------------------------------
as1 <- function(nqpoints=50,qrange=4){
	qs1 <- GL.YW(nqpoints,c(-qrange,qrange))
	qs <- list(upts=pnorm(qs1[,1]),zpts=qs1[,1],zwts=qs1[,2])
return(qs)
}
# ------------------------------------------------------------------------- 
as2 <- function(quadlimits=c(-8,-5,-3,-1,1,3,5,8),nzvec=c(6,10,10,20,10,10,6)){
#autocreate piecewise quadrature scheme

	nzPoints <- sum(nzvec)
	zPoints <- NULL
	zWts <- NULL
	nsegs <- length(quadlimits)-1
#	print(nsegs)
#	print(quadlimits)
	for(i in 1:nsegs){
#	 print(nzvec[i])
		uw <- GL.YW(nzvec[i],quadlimits[i:(i+1)])
		zPoints <- c(zPoints,uw[,1])
		zWts <- c(zWts,uw[,2])
	} 
	upts <- pnorm(zPoints)	       
	ans <- list(upts=upts,zpts=zPoints,zwts=zWts)
return(ans)
}
# ------------------------------------------------------------------------- 
	#quadrature scheme for piecewise integration
	#integrating from 0 to 1 and want total integral
	#and each piece--this routine gives concatenation
	#of evaluation points and weights
	#zPoints unevenly spaced over (0,1).   
	
pwts1 <- function(zPoints,nqpts=4){
	zPoints <- c(0.0,zPoints,1.0) 
	nzPoints <- length(zPoints)-1
	YW <- matrix(0,nrow=nqpts*nzPoints,nc=2) 
	for(i in 1:nzPoints) YW[(nqpts*(i-1)+1):(nqpts*i),] <- GL.YW(nqpts,xrange=c(zPoints[i],zPoints[i+1]))   
return(YW)
}
# -------------------------------------------------------------------------


# check note: rescaled smolyak goes beyond 1 !!!
#rg <- c(-4,5)
#x1 <- quadSmolyak(3,xrange=rg)
#res1 <- sum(x1$wt*dnorm(x1$pt))

#a <- length(x1$pt)+1
#x2 <- GL.YW(a,xrange=rg)

#res2 <- sum(x2[,2]*dnorm(x2[,1]))

