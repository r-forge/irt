# estimation parameters come from iSet
# nitems, nW, nWI, nWsig and corresponding names come from dSet), also ncats
# zpts, upts, zwts come from qscheme; was: zpts <- qscheme$zpts; upts <- qscheme$upts; zwts <- qscheme$zwts
# --------------------------------------------------------------
#was: irtM1 <- function(r = stop("item argument is missing"),
irtMulti <- function(r = stop("item argument is missing"),
		W = stop("design matrix is missing"), 
		W.I=NULL,
		W.sig=NULL,
		thstart=NULL, # is also in settings
		qscheme=as2(),
		scheme=NULL, 
		iSet=irtSet(), 
		dSet=NULL,
		JH=TRUE){

	if(!irtSetValidate(iSet)) stop("iSet object is not of the 'irtSet' class")
	if(is.null(dSet)) dSet <- dtaSet(r=r,W=W,WI=W.I,Wsig=W.sig)		 
	if(!dtaSetValidate(dSet)) stop("Problems with data settings ...")

	upts <- qscheme$upts; zpts <- qscheme$zpts; zwts <- qscheme$zwts
		
	nb <- iSet$nb
	plt <- iSet$plt
	pl <- iSet$pl
	iterlim <- iSet$iterlim
	gradtol <- iSet$gradtol
	reltol <- iSet$reltol
	saveres <- iSet$saveres
	tag <- iSet$tag
	orthobasis <- irtSet$orthobasis
	Fmethod <- iSet$Fmethod
	ctype <- iSet$ctype
	algorithm <- iSet$algorithm
	bestllf <- iSet$bestllf
	nq <- iSet$nq
	orthobasis <- iSet$orthobasis
	bfile <- iSet$bfile
		
	nwcol <- dSet$nwcol
	nwicol <- dSet$nwicol
	nwsigcol <- dSet$nwsigcol
	nitems <- dSet$nitems
	nobs <- dSet$nobs
	ncats <- dSet$ncats

	if(is.null(scheme)) scheme <- pwts1(zPoints=upts,nqpts=nq)        #if we are not given a scheme, make one

	r <- as.matrix(r) # this is important if there is one item
	#now for the 'internal' numerical integration
	#if(Fmethod=="logistic"){u.log <<- log(upts/(1-upts))}

	nu <- length(upts)+1
	nutot <- nrow(scheme)
	nq <- nutot/nu 			# why once again? due to scheme - extended basis
	utot <- scheme[,1]
	udex <- (1:nu)*nq
	etu <- rep(0,nrow(scheme))

	if(orthobasis){ ubasis <- ortbasis(utot,nb)}
	else{ ubasis <- sapply(1:nb,function(i) upts^i) } # note upts -> was u (= upts)
	if(Fmethod=="mixET"){  ubasis <- sapply(1:nb,function(i) upts^i) }# ks added  
	       
	ncurves <- sum(ncats-1)		
	nicol <- nwicol+nb
	nicoefs <- (ncurves)*nicol    
	nmucoefs <- nwcol
	nsigcoefs <- nwsigcol

	if(is.null(thstart)){ 
		icoefs <- rep(.05,nicoefs)
		mucoefs <- rep(.02,nmucoefs)
		sigcoefs <- NULL
		if(nsigcoefs>0 ) sigcoefs <- rep(.01,nsigcoefs)
		thstart <- c(icoefs,mucoefs,sigcoefs)
	}# validation for thstart would be nice	
	
	starttime <- proc.time()[3]

	irtMultiObjfn <- irtMultiObjfnConstruct(ncurves,nicol,nobs,ncats,
								nicoefs,nmucoefs,nsigcoefs,
								zwts,upts,zpts,
								Fmethod, ctype, plt,bestllf,starttime,bfile,
								ubasis,scheme,qscheme,udex,
								r,W,W.I,W.sig,algorithm,
								icoefs,mucoefs,sigcoefs
								)
	
	
	switch( algorithm,
	   	 nlm = ans <- nlm(irtMultiObjfn,thstart,print.level=pl,iterlim=iterlim,steptol=1.e-14,gradtol=gradtol,fscale=1000,ndigit=12),
	   optim = ans <- optim(thstart,irtMultiObjfn,control=list(trace=5,REPORT=1,reltol=reltol,maxit=iterlim),method="BFGS"),
	  nlmraw = ans <- nlm(irtMultiObjfn,thstart,print.level=pl,iterlim=iterlim),						
		  BB = ans <- spg(thstart,irtMultiObjfn,method=1,control=list(triter=1)), 
	  nlminb = ans <- nlminb(thstart,irtMultiObjfn,control=list(trace=1))
	#  	  genoud = ans <- genoud(objfn1,nvars=length(thstart),starting.values=thstart,control=list(trace=5,REPORT=1)),
	#  	  genoud = ans <- genoud(objfn1,nvars=length(thstart),starting.values=thstart,control=list(trace=5,REPORT=1),optSet=objSet,data=dta)),
	)
	#  if() may be necessary here ... check it out
	names(ans)[names(ans)=="par"] <- "estimate"
	names(ans)[names(ans)=="minimum"] <- "value"
	names(ans)[names(ans)=="objective"] <- "value"

#	bestllf <- get("bestllf")

	print("returning from optimization algorithm.")

	FinalCall <- irtMultiObjfn(ans$estimate,finalcall=TRUE)
	print("done")
	Ftab <- FinalCall$Ftab; class(Ftab) <- "Ftab"

	JHlist <- list()
	if(JH){
		JHlist <- JHfun(irtMultiObjfn,ans$estimate)
		}else{
			JHlist$segg <- rep(NA,length(thstart))
			JHlist$sew <- rep(NA,length(thstart))
			JHlist$seh <- rep(NA,length(thstart))
		}

	res <- list(estimate=ans$estimate,value = ans$value,optimans= ans,
				qscheme = qscheme, scheme=scheme, thstart=thstart,
				iSet=iSet,Ftab=Ftab,
				JHlist=JHlist,
				FinalCall=FinalCall,
				nicoefs=nicoefs,ncats=ncats, 
				W=W, r=r,
				irtMultiObjfn=irtMultiObjfn
				)
								
	class(res) <- "irtMulti"
#	if(saveres) saveresfun(res,tag) #well this looks crazy. Intended to produce the result that load(Rfile.name) creates something with that name
return(res)
}
# -------------------------------------------------------------------------
irtMultiObjfnConstruct <- function(ncurves,nicol,nobs,ncats,
							nicoefs,nmucoefs,nsigcoefs,
							zwts,upts,zpts,
							Fmethod, ctype, plt,bestllf,starttime,bfile,
							ubasis,scheme,qscheme,udex,
							r,W,W.I,W.sig,algorithm,
							icoefs,mucoefs,sigcoefs							
							){


    nqpoints <- length(upts)
	Pone <- matrix(1,nr=nobs,nc=nqpoints)
	nitems <- length(ncats)
	nb <- nicoefs/ncurves
	llfnow  <- llfvec <- rep(0,nobs)

	#save repeated index calculations in objfn2 ; to separate parameters
		smudex <- nicoefs+1
		emudex <- nicoefs+nmucoefs
		if(nsigcoefs>0){
		    ssigdex <- emudex+1
		    esigdex <- emudex+nsigcoefs
		}else{
		    ssigdex <- NULL
		    esigdex <- NULL		
		}
		
	itemPMulti <- function(b,r){ # ks update
		nb <- length(b)/(ncats-1)
		curvemat <- matrix(1,nrow=nqpoints,ncol=ncats)

		if(ctype=="1-F"){
			for( j in (ncats-1):1 ){
				bs <- (j-1)*nb+1
				tvec <- b[bs:(bs+nb-1)]
				Fnow <- Fx(tvec=tvec,u=ubasis[udex,],scheme=scheme[udex,2],type=Fmethod)		
				curvemat[,j] <- curvemat[,j+1]*(1-Fnow) 
			}
			pmat <- sapply(2:ncats,function(i) curvemat[,i]-curvemat[,i-1])
			pmat <- t( cbind(curvemat[,1],pmat) ) # transposition: if not here then at the upper level for dn ;( no way out
		}	
		if(ctype=="F"){
			for( j in 1:(ncats-1) ){
				bs <- (j-1)*nb+1
				tvec <- b[bs:(bs+nb-1)]
				Fnow <- Fx(tvec=tvec,u=ubasis[udex,],scheme=scheme[udex,2],type=Fmethod)		
				if(j==1) curvemat[,j] <- Fnow 
				else curvemat[,j] <- curvemat[,j-1]*Fnow
			}
			for(i in 1:(ncats-1)) curvemat[,i] <- 1 - curvemat[,i] #reverse each curve	
			pmat <- sapply(2:ncats,function(i) curvemat[,i]-curvemat[,i-1])
			pmat <- t(cbind(curvemat[,1],pmat) )  # transposition: if not here then at the upper level for dn ;(, no way out
		}
		Ftab <- curvemat[,-ncats]
		P <- pmat[r,]  # note here - dependence on r (items)!!!
	return(list(P=P,Ftab=Ftab))
	}	# end itemPMulti						


	cPMulti <- function(b,r){
		P <- matrix(1,nr=nobs,nc=nqpoints)
		fun <- function(i){
			itemPres <- itemPMulti(b[,i],r[,i])
			P <- itemPres$P*P
	    return( data.frame(x=qscheme$zpts,itemPres$Ftab)	)	# you need this for plotting
		}
		FtabList <- lapply(1:nitems,fun)
		class(FtabList) <- "Ftab"	 	   
	return(list(P=P,FtabList=FtabList))	       #now P is nobs by nqpoints ;integrate it
	} # end cPMulti


	ff <- function(theta,Jacmode=FALSE, finalcall = FALSE,plt=TRUE,...){
		
		icoefs <- matrix(theta[1:nicoefs],nc=nitems) # note a change: was nc=nicol ... ?		
		if(nmucoefs>0)	mucoefs <- theta[smudex:emudex]
		if(nsigcoefs>0)	sigcoefs <- theta[ssigdex:esigdex]
		if(nsigcoefs>0)	sigvec <-exp(W.sig%*%sigcoefs)
		else 			sigvec <- 1
		if(nmucoefs >0)	muvec <- W%*%mucoefs 
		else 			muvec <- 1

    	cPres <- cPMulti(icoefs,r)	
		P <- cPres$P ; FtabList <- cPres$FtabList
	 	dn <- sapply(zpts,dnorm,mean=muvec,sd=sigvec) # density for muvec f(theta|W) 
		Pfunc <- P*dn
		P <- (P*dn)%*%zwts
       
		llfvec <- log(P)
		llfnow <- sum(llfvec)

		if(plt) plotItemsGrdInfo(unclass(FtabList),colnames(r),list(starttime=starttime,llf=llfnow),showinfo=TRUE) # get showinfo from iSet
	
		if(!is.na(llfnow)){
#			if(!Jacmode&(llfnow>bestllf)){ 
#				assign("bestllf",llfnow,envir=environment(irtMulti),inherits=TRUE)
#				wmat(file=bfile,data=matrix(theta,nc=1))
#			}
		}
		
		if(finalcall){
			ftheta <- apply(Pfunc,2,"/",P)
			llfvec <- llfvec						#llfvec=log(P)
			Etheta <- (ftheta)%*%(zpts*zwts)
			Vtheta <- (ftheta)%*%(zpts^2*zwts) 	#moment around zero 
			Vtheta <- Vtheta-Etheta^2
			   #you can get other moments using components in "ans"
			res <- list(Ftab=FtabList,ftheta = ftheta,llfvec=llfvec,Etheta=Etheta, Vtheta=Vtheta,muvec=muvec,sigvec=sigvec)
#			res <- list(Ftab=FtabList)
		return(res)			
		} 
		else ifelse(Jacmode, return(-llfvec), return(-llfnow) )
	} # end ff

return(ff)
}
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
JHfun <- function(objFun,theta){

# compute Jacobian
	cat("Compute Jacobian: ")
	jac <- jacobian(objFun,theta,method="Richardson", method.args=list(),Jacmode=TRUE,plt=FALSE) 
	cat("Jacobian - DONE.")
	cat("\nCompute Hessian : ")
# compute Hessian
	hess <- hessian(objFun,theta,method="Richardson", method.args=list(),Jacmode=FALSE,plt=FALSE) 
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
print.irtMulti <- function(x,...){	
	output1 <- data.frame(x$thstart,as.matrix(x$estimate))
	names(output1) <- c("St. Val.","estimates")
	tiltnm <- paste("t",1:x$nicoefs,sep="")
	coefnm <- c(tiltnm,colnames(x$W))   # substitute with nw from future dataSettings
	rownames(output1)  <- coefnm 

	cat("\n","Coefficients :", "\n") 
	cat("---------------------------------", "\n") 
	print(round(output1,4),...)
		
	out <- list(coeftab=output1)	
invisible(out)	
}
# -------------------------------------------------------------------------
summary.irtMulti <- function(object,startval=FALSE,...){   
	thstart <- object$thstart
	theta <- object$estimate
	W <- object$W
	nicoefs <- object$nicoefs
	JHlist <- object$JHlist

	segg <- JHlist$segg
	seh <- JHlist$seh
	sew <- JHlist$sew	
	
	coeftab <- cbind(theta,segg,seh,sew)
	if(startval) coeftab <- cbind(thstart,coeftab)

	invisible( ifelse(is.null(colnames(W)), Wnm <- paste("feature",1:ncol(W),sep=" "), Wnm <- colnames(W)) )
	tiltnm <- paste("t",1:nicoefs,sep="")
	coefnames <- c(tiltnm,Wnm)	
	rownames(coeftab) <- coefnames   
	
	cat("\n","Coefficients :", "\n") 
	cat("---------------------------------", "\n") 
	print(round(coeftab,4),...)

	out <- list(coeftab=coeftab)	
invisible(out)
}
# -------------------------------------------------------------------
plot.irtMulti <- function(x,...){
	Ftab <- x$Ftab
	plot(Ftab,...)
}


