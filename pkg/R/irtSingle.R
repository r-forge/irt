#item response model for single item
#W matrix of conditioning variables
#r ordinal responses          
# -------------------------------------------------------------------------
irtSingle <- function(r = stop("item argument is missing"),
		W = stop("design matrix is missing"), 
		W.I=NULL,
		W.sig=NULL,
		thstart=NULL, # is also in settings
		qscheme=as1(nqpoints=90,qrange=4.0),
		scheme=NULL, 
		iSet=irtSet(), 
		dSet=NULL,
		itemnm=NULL){
				
	if(!irtSetValidate(iSet)) stop("iSet object is not of the 'irtSet' class")
	if(is.null(dSet)) dSet <- dtaSet(r=r,W=W)		 
	if(!dtaSetValidate(dSet)) stop("Problems with data settings ...")
	dSet$rnames <- itemnm # because column of a matrix was coerced to a vector without name

	# this is from qscheme()	
	upts <- qscheme$upts; zpts <- qscheme$zpts; zwts <- qscheme$zwts
	# this is from iSet()	
	nb <- iSet$nb
	plt <- iSet$plt
	pl <- iSet$pl
	iterlim <- iSet$iterlim
	gradtol <- iSet$gradtol
	saveres <- iSet$saveres
	tag <- iSet$tag

	# this is from dSet()	
	nwcol <- dSet$nwcol
	nwsigcol <- dSet$nwsigcol
	nwicol  <- dSet$nwicol
	nitems <- dSet$nitems
	nobs <- dSet$nobs
	ncats <- dSet$ncats
	
	ncurves <- ncats -1 
	nicol <- nb #	nicol <- nwicol+nb # this is in irtM1()
	nicoefs <- (ncurves)*nicol; 
	nmucoefs 	<- nwcol
	nsigcoefs <- nwsigcol

	scheme <- pwts1(zPoints=upts,nqpts=8) #nqpts==how many 'in between'
	if(is.null(thstart)){ 
		icoefs <- rep(.05,nicoefs)
		mucoefs <- rep(.02,nmucoefs)
		thstart <- c(icoefs,mucoefs)
	}# validation for thstart would be nice

#	starttime <- proc.time()[3]	

	# function constructor
	irtSingleObjfn <- irtSingleObjfnConstruct(nicoefs,nicol,upts,ncats,ncurves,zwts,zpts,icoefs,scheme,plt,r,W) 

	ans <- nlm(irtSingleObjfn,thstart,print.level=pl,iterlim=iterlim,gradtol=gradtol,fscale=1000,ndigit=12,steptol=1.e-14)
	Ftab <- irtSingleObjfn(ans$estimate,finalcall=TRUE)
	class(Ftab) <- "Ftab"
	JH <- TRUE
	ifelse(JH,JHlist <- JHfun(irtSingleObjfn,ans$estimate),JHlist <- NULL)
	
	result <- list(
				estimate=ans$estimate,value = ans$value,optimans=ans,
				qscheme = qscheme,scheme=scheme, thstart=thstart,
				iSet=iSet,Ftab=Ftab,W=W,r=r,
				nicoefs=nicoefs,ncats=ncats, JHlist=JHlist)
					 	   					
	class(result) <- "irtSingle"
	if(saveres) saveresfun(result,tag) 
return(result)
}      
# -------------------------------------------------------------------------
irtSingleObjfnConstruct <- function(nicoefs,nicol,upts,ncats,ncurves,zwts,zpts,icoefs,scheme,plt,r,W){
	nqpoints <- length(upts)
	utot <- scheme[,1]
	utot1 <- 2*utot-1
	utot2 <- 6*utot^2-6*utot+1
	utot3 <- 20*utot^3-30*utot^2+12*utot-1
	utotl <- cbind(utot1,utot2,utot3)
	
	nu <- length(upts)+1 
	nutot <- nrow(scheme)
	nq <- nutot/nu
	udex <- (1:nu)*nq

	etF <- function(tvec){
		etu <- utotl[,1:length(tvec)]%*%tvec
		etu <- exp(etu)*scheme[,2]
		qvec <- cumsum(etu)/sum(etu)
		ans <- qvec[udex]
		ans <- ans[-length(ans)]  #drop the last one; it is 1 by defn
	}
									
	cP3 <- function(icoefs){
		Ftab <- matrix(0,nc=ncats,nr=nqpoints)  
		Ftab[,ncats] <- 1.0

		for(i in (ncurves:1)) Ftab[,i] <- (1-etF(icoefs[i,]))*Ftab[,i+1]
		Ptab <- rbind( Ftab[,1],apply(Ftab,1,diff) ) # was #			Ptab[,i+1] <- Ftab[,i+1]-Ftab[,i]; Ptab[,1] <- Ftab[,1]
	return(  list( P=Ptab[r,],Ftab=data.frame(x=zpts,Ftab) )  )	
	}
	
	ff <- function(theta,finalcall=FALSE,Jacmode=FALSE){ 
		icoefs <- matrix(theta[1:nicoefs],nc=nicol)
		mucoefs <- theta[-(1:nicoefs)]
		cP2res <- cP3(icoefs)

		P <- cP2res$P ; Ftab <- cP2res$Ftab #; Ftab <- list(cP2res$Ftab)

		muvec <- as.matrix(W)%*%mucoefs	
		dn <- sapply(zpts,dnorm,mean=muvec) # density for muvec f(theta|W) 
		pvec <- (P*dn)%*%zwts
		llfvec <- log(pvec)
		llfnow <- sum(llfvec)
			
		if(finalcall) return(list(Ftab=Ftab))	
		else ifelse(Jacmode, return(-llfvec), return(-llfnow) )
	}
return(ff)
}
# -------------------------------------------------------------------------
JHfun <- function(objFun,theta){

# compute Jacobian
	cat("Compute Jacobian: ")
	jac <- jacobian(objFun,theta,method="Richardson", method.args=list(),Jacmode=TRUE) 
	cat("Jacobian - DONE.")
	cat("\nCompute Hessian : ")
# compute Hessian
	hess <- hessian(objFun,theta,method="Richardson", method.args=list(),Jacmode=FALSE) 
	cat("Hessian - DONE.")
	
# estimates:
	J <- t(jac)%*%jac
	segg <- sqrt(diag(solve(J,tol=1.e-40)))
	hinv <- solve(hess,tol=1.e-40)
	seh <- sqrt(diag(hinv))

	wv <- hinv%*%J%*%hinv
	sew <- sqrt(diag(wv))  	  

return(list(segg=segg,seh=seh,sew=sew))	
}
# -------------------------------------------------------------------------
jacprint <- function(inittime,type=TRUE){
	timenow <- proc.time()[3]
	if( ( (timenow-inittime)%%5 )< 0.2)  {
		if(type) cat("('J') ")
		else cat("('H') ")		
	}
	#ifelse(type, cat(" ('J') "), cat(" ('H') ") )
}
# -------------------------------------------------------------------------
print.irtSingle <- function(x,...){      

	estimate <- x$estimate
	output1 <- data.frame(as.matrix(estimate))
	names(output1) <- "estimates"
	tiltnm <- paste("t",1:x$nicoefs,sep="")
	coefnm <- c(tiltnm,colnames(x$W))   
	rownames(output1)  <- coefnm 
	cat("\n","Coefficients :", "\n") 
	cat("---------------------------------", "\n") 
	print(round(output1,4),...)
}     
# -------------------------------------------------------------------------
summary.irtSingle <- function(object,...){   
	iSet <- object$iSet
	ncats <- object$ncats
	nicoefs <- object$nicoefs
	theta <- object$estimate
	
	nr <- ncats
	W <- object$W
	r <- object$r

	invisible( ifelse(is.null(names(W)), Wnm <- paste("feature",1:ncol(W),sep=" "), Wnm <- names(W)) )
	tiltnm <- paste("tilt",1:nicoefs,sep=" ")
	coefnames <- c(tiltnm,Wnm)

	coeftab <- data.frame(as.matrix(theta),row.names=coefnames)
	names(coeftab) <- "Estimates"

	JHlist <- object$JHlist
	if(!is.null(JHlist)){
		output1 <- cbind(theta,JHlist$segg,JHlist$seh,JHlist$sew)
		cnames <- c(paste("t",1:nicoefs,sep=""),colnames(W))
		colnames(output1) <- c("theta","segg","seh","sew")
		rownames(output1) <- cnames   

		cat("\n","Coefficients Table:", "\n") 
		print(round(output1,4),...)        	
	}
}
# -------------------------------------------------------------------
plot.irtSingle <- function(x,...){
	res <- x
	Ftab <- res$Ftab
	plot(Ftab,...)
}
# -------------------------------------------------------------------