# July 9 2008
#b=parameters of irt curves
#r=responses
#dvec=description of problem (a vector with nfpar,ncats
# -------------------------------------------------	
cP<- function(b,r,dvec,upts,Fmethod,ctype,ubasis,scheme,udex,qscheme){
	nfpar <- dvec[1]
	ncats <- dvec[-1]
	nitems <- length(ncats)
	ncurves <- sum(ncats)-nitems
	nb <- length(b[1,])/ncurves
	nqpoints <- length(upts)
	nobs <- nrow(r)

	P <- matrix(1,nr=nobs,nc=nqpoints)
	edex <-nb*cumsum(ncats-1)
	sdex <- c(1,nb*cumsum(ncats-1)+1)
	FtabList <- list()
	for(i in 1:nitems){
	       bnow <- b[,sdex[i]:edex[i]]		#temp patch assume every item has 4 categories, use P4
	       itemPres <- itemP(bnow,r[,i],upts,Fmethod,ctype,ubasis,scheme,udex)
		   pnow <- itemPres$P
		
	       P <- pnow*P
#	       Ftab.now[[i]] <- itemPres$thisFtab #$thisFtab # where does this go ...?side effect of P4b() is specimen Ftab in thisFtab, global env
	       FtabList[[i]] <- data.frame(x=qscheme$zpts,itemPres$Ftab)
       }
       class(FtabList) <- "Ftab"	 	   
return(list(P=P,FtabList=FtabList))	       #now P is nobs by nqpoints ;integrate it
}
#-------------------------------------------------------------
itemP <- function(b,r,u,Fmethod,ctype,ubasis,scheme,udex){ # ks update
	nobs <- length(r)
	nqpoints <- length(u)
	ncat <- length(unique(r))
	nb <- length(b)/(ncat-1)
	curvemat <- matrix(1,nrow=nqpoints,ncol=ncat)
	
	if(ctype=="1-F"){
		for( j in (ncat-1):1 ){
			bs <- (j-1)*nb+1
			tvec <- b[bs:(bs+nb-1)]
			Fnow <- Fx(tvec=tvec,u=ubasis[udex,],scheme=scheme[udex,2],type=Fmethod)		
			curvemat[,j] <- curvemat[,j+1]*(1-Fnow) 
		}
		pmat <- sapply(2:ncat,function(i) curvemat[,i]-curvemat[,i-1])
		pmat <- t( cbind(curvemat[,1],pmat) ) # transposition: if not here then at the upper level for dn ;( no way out
	}	
	if(ctype=="F"){
		for( j in 1:(ncat-1) ){
			bs <- (j-1)*nb+1
			tvec <- b[bs:(bs+nb-1)]
			Fnow <- Fx(tvec=tvec,u=ubasis[udex,],scheme=scheme[udex,2],type=Fmethod)		
			if(j==1) curvemat[,j] <- Fnow 
			else curvemat[,j] <- curvemat[,j-1]*Fnow
		}
		for(i in 1:(ncat-1)) curvemat[,i] <- 1 - curvemat[,i] #reverse each curve	
		pmat <- sapply(2:ncat,function(i) curvemat[,i]-curvemat[,i-1])
		pmat <- t(cbind(curvemat[,1],pmat) )  # transposition: if not here then at the upper level for dn ;(, no way out
	}
	Ftab <- curvemat[,-ncat]
	P <- pmat[r,]  # note here - dependence on r (items)!!!
return(list(P=P,Ftab=Ftab))
}
#-------------------------------------------------------------
#  look at the input; maybe has to be change	 
Fx <- function(tvec,u,scheme=NULL,type="ET"){
	switch(type,
			  ET =  etFx2(tvec,u,scheme),
		     iET = ietF(tvec,u),
		    beta = betaFx(tvec,u),
		logistic = logisticF(tvec,u),
#		   mixET = mixET(t=tvec),
#		     snp = snpF(t=tvec),
#		snp.cheb = snpF(t=tvec),
#		 snp.sin = snpF(t=tvec)
	)
}
#--------------------------------

