# -----------------------------------------------------------------
plotItemsStd <- function(...){
	UseMethod("plotItemsStd")
}
plotItemsStd.default <- function(x,...){
	args <- list(...)
	plotFtab(x,...)
}
plotItemsStd.list <- function(x,...){
	plotFtabList(x,...)
}
# ------------------------------------------------------------------------- 
plotFtabList <- function(datalist,ncol=2,xscale=NULL,itemnm=NULL){
	nplots <- length(datalist)
	if(is.null(itemnm)){
		if( is.null(names(datalist)) )	itemnm <-  paste("item",1:nplots,sep=" ")
		else itemnm <-  names(datalist)
	}

	if(iseven(nplots)) nrow <- nplots%/%ncol
	else nrow <- nplots%/%ncol + 1
	if(nplots==1) ncol <- 1 
	parFtab <- par(mfrow=c(nrow,ncol))
	for (i in 1:nplots) plotFtab(datalist[[i]],itemnm=itemnm[i],xscale=xscale)
#invisible(list(parFtab=parFtab))
}      
# -------------------------------------------------------------------------
plotFtab <- function(data,xscale=NULL,itemnm=NULL){
	if(is.null(itemnm)) itemnm <- "item"
	ncol <- ncol(data)
	x <- data[,1]
	plot(x,data[,2],type='l',xlim=xscale,ylim=c(0,1),xlab="theta",ylab= "F(r|theta)")	
 	for(i in 2:ncol) lines(x,data[,i])
	title(main=itemnm)
}
# -------------------------------------------------------------------------
# ------------------------------------------------------------------------- 
plotItemsStdInfo <- function(datalist,ncol=2,xscale=NULL,itemnm=NULL,plotInfoSet=NULL){
	nplots <- length(datalist)
	if(is.null(itemnm)){
		if( is.null(names(datalist)) )	itemnm <-  paste("item",1:nplots,sep=" ")
		else itemnm <-  names(datalist)
	}

	if(iseven(nplots)) nrow <- nplots%/%ncol
	else nrow <- nplots%/%ncol + 1

	if(! is.null(plotInfoSet)){
		timenow <- proc.time()[3]
		if( ( (timenow-plotInfoSet$starttime)%%10 )< 0.3){
	 		parFtab <- par(mfrow=c(nrow,ncol),oma = c(0,0,5,0) )
			for (i in 1:nplots) plotFtab(datalist[[i]],itemnm=itemnm[i],xscale=xscale)
			par(parFtab)
			title <- paste("llf=",plotInfoSet$llf)
			mtext(title,side=3, font=2, cex=1.5, col='red')
		}
	}
	else{	
		parFtab <- par(mfrow=c(nrow,ncol))
		for (i in 1:nplots) plotFtab(datalist[[i]],itemnm=itemnm[i],xscale=xscale)
	}
}      
# -------------------------------------------------------------------------

