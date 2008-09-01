# -------------------------------------------------------------------------
irtSet <- function(
			reltol = 1.e-09,
			gradtol = 1.e-11,
			iterlim = 250,
			pl = 2,
			nq = 8,
			nb = 2,
			algorithm = c("nlm","optim","nlmraw","BB","genoud","nlminb"),
			Fmethod = c("ET","logistic","mixET"),
			ctype = c("F","1-F"),
			orthobasis = TRUE,
			bestllf = -99999999999999999.,
			tag = NULL,
			plt=FALSE,
			saveres=TRUE
			){

	if(is.null(tag)){ 	#generate tag
		ifelse(file.exists('runnumber.irt'), runnumber <- rmat('runnumber.irt')+1, runnumber <- 1)
	    wmat(data=matrix(runnumber,1,1),file='runnumber.irt')
#		tag <- as.character(runnumber)
		tag <- paste(Sys.Date(),"-run",as.character(runnumber),sep="")
    }

# default settings
		irtSetObj <- list()
		irtSetObj$reltol <- reltol
		irtSetObj$gradtol <- gradtol
		irtSetObj$iterlim <- iterlim
		irtSetObj$pl <- pl
		irtSetObj$nq <- nq
		irtSetObj$nb <- nb
		irtSetObj$algorithm <- match.arg(algorithm)
		irtSetObj$Fmethod <- match.arg(Fmethod)
		irtSetObj$ctype <- match.arg(ctype)
		irtSetObj$orthobasis <- orthobasis
		irtSetObj$bestllf <- bestllf	
		irtSetObj$tag <- tag
		irtSetObj$plt <- plt
		irtSetObj$saveres <- TRUE
		irtSetObj$bfile <- paste("bbest.",tag,sep="")
		
	class(irtSetObj) <- "irtSet"
return(irtSetObj)	
}
# -------------------------------------------------------------------------
print.irtSet <- function(x,digits = max(5, .Options$digits - 2), ext = FALSE,...){      
	irtSetObj <- x
	cat("\n","Basic settings for irt optimization:")
	cat("\n","========================================","")
	cat("\n","       algorithm :",irtSetObj$algorithm)
	cat("\n","         Fmethod :",irtSetObj$Fmethod)
	cat("\n","           ctype :",irtSetObj$ctype)
	cat("\n","              nb :",irtSetObj$nb)
	cat("\n","         gradtol :",irtSetObj$gradtol)
	cat("\n","          reltol :",irtSetObj$reltol)	
	cat("\n","         iterlim :",irtSetObj$iterlim)
	cat("\n","              pl :",irtSetObj$pl)
	cat("\n","             tag :",irtSetObj$tag)	
	cat("\n","           bfile :",irtSetObj$bfile)
	cat("\n","      orthobasis :",irtSetObj$orthobasis)	
invisible(irtSetObj)
}     
# -------------------------------------------------------------------------
irtSetValidate <- function(irtSetObj){
# validation methods
	valres <- FALSE
	if(class(irtSetObj)=="irtSet") valres <- TRUE
return(valres)
}
# -------------------------------------------------------------------------
# DATA SETTINGS
# -------------------------------------------------------------------------
dtaSet <- function(r,W=NULL,WI=NULL,Wsig=NULL){
	r <- as.matrix(r)
	rnames <- colnames(r)
	nitems <- ncol(r)
	nobs <- nrow(r)

	if(!is.null(W)){ 
		wnames <- colnames(W)
		nwcol <- ncol(W)
	}else{
		wnames <- "There aren't covariates."
		nwcol <- 0				
	}

	if(!is.null(WI)){ 
		winames <- colnames(WI)
		nwicol <- ncol(WI)
	}else{
		winames <- "There aren't individual covariates."
		nwicol <- 0		
	}

	if(!is.null(WI)){ 
		wsignames <- colnames(Wsig)
		nwsigcol <- ncol(Wsig)
	}else{
		wsignames <- "There aren't sigma covariates."
		nwsigcol <- 0		
	}

	ncats <- sapply(1:nitems,function(i) length(unique(r[,i])))
	
	res <-	list(rnames=rnames,nitems=nitems,
				 wnames=wnames, nwcol=nwcol,
				 winames=winames,nwicol =nwicol,
				 wsignames=wsignames,nwsigcol=nwsigcol,
				nobs=nobs,
				ncats = ncats)
	class(res) <- "dtaSet"
return(res)
}
# -------------------------------------------------------------------------
dtaSetValidate <- function(dtaSetObj){
return(TRUE)
}
# -------------------------------------------------------------------------
print.dtaSet <- function(dtaSetObj){
	rnames <- dtaSetObj$rnames
	wnames <- dtaSetObj$wnames

	lnrnames <- length(rnames)
	lnwnames <- length(wnames)

	cat("\n","Data settings:")
	cat("\n","========================================","")
	cat("\n","Items:")
	cat("\n","----------------------------------------","")
	for(i in 1:lnrnames) cat("\n","item",i,"name :",rnames[i])	
		
	if(!is.null(wnames)){
		cat("\n\n","Set of covariates (the same for all items):")
		cat("\n","----------------------------------------","")
		for(i in 1:lnwnames) cat("\n","cov.",i,":",wnames[i])	
	}
}
# -------------------------------------------------------------------------