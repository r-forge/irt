# some useful function taken from np package
# -------------------------------------------------------------------
ortbasis <- function(utot=10,nb=1){
	if(nb >= 6) stop("'nb' can not exceed 5 for orthonormal basis.")

	ubasis <- matrix(0,nrow=length(utot),ncol=nb)
	utemp <- sapply(1:nb,function(i) utot^i) 
	ubasis[,1] <- (2*utemp[,1]-1)
	if(nb>=2)ubasis[,2] <- (6*utemp[,2]-6*utemp[,1]+1)
	if(nb>=3)ubasis[,3] <- (20*utemp[,3]-30*utemp[,2]+12*utemp[,1]-1)
	if(nb>=4)ubasis[,4] <- (70*utemp[,4]-140*utemp[,3]+90*utemp[,2]-20*utemp[,1]+1)
	if(nb>=5)ubasis[,5] <- (252*utemp[,5]-630*utemp[,4]+560*utemp[,3]-210*utemp[,2]+30*utemp[,1]-1)	
return(ubasis)
}
# -------------------------------------------------------------------
plot.Ftab <- function(x,...,type=c("standard","grid"), plotitems=NULL){

	if(!is.null(plotitems)) Ftab <- x[plotitems]	
	tp <- match.arg(type)
	switch(tp,
		standard = plotItemsStd(unclass(Ftab),...),
		grid = plotItemsGrd(unclass(Ftab),...)
		)
}
# -------------------------------------------------------------------
dlev <- function(x){
  if(is.ordered(x))
    x.dlev <- suppressWarnings(as.numeric(levels(x)))
  if (!is.ordered(x) || any(is.na(x.dlev)))
    x.dlev <- as.numeric(1:nlevels(x))
  x.dlev
}
# -------------------------------------------------------------------
assignList <- function(lst) for(i in 1:length(lst)) assign(names(lst)[i],lst[[i]],pos=+1) # assign list values
# -------------------------------------------------------------------
isodd <- function(nr){	return( ifelse((nr%%2) == 0, FALSE, TRUE) )}
iseven <- function(nr){	return( ifelse(( (nr%%2) == 0 ), TRUE, FALSE) )}
# -------------------------------------------------------------------
# read / write matrix
rmat <- function(file){
	temp <- scan(file)
	dims <- temp[1:2]
	temp <- temp[-c(1,2)]
	matrix(temp,nr=dims[1],byrow=TRUE)
}
# -------------------------------------------------------------------
wmat <- function(file,data){ #inverse of rmat...write a .mat
                             #file for Ox to read
	write(dim(data),file=file)
	write(t(data),file=file,ncol=ncol(data),append=TRUE)
}
# -------------------------------------------------------------------
saveresfun <- function(res,tag=NULL){
	Rfile.name <- paste("irt",tag,".ans",sep='')
	save(res,file=Rfile.name)
}
# -------------------------------------------------------------------
toMatrix <- function(data) {
  tq <- sapply(data, function(y){ if(is.factor(y) ) y <- as.integer(dlev(y)[y]) ; y } )
  dim(tq) <- dim(data) ## cover the corner case of single element d.f.
return(tq)
}
# -------------------------------------------------------------------
toDataFrame <- function(data) {
	tq <- sapply(data, function(y){ if(is.factor(y) ) y <- as.integer(dlev(y)[y]) ; y } )
	dim(tq) <- dim(data) ## cover the corner case of single element d.f.
	tq <- data.frame(tq)
return(tq)
}
# -------------------------------------------------------------------
toFrame <- function(frame) {
  if(!is.data.frame(frame)){
    t.names <- NULL

    if(!(is.vector(frame) || is.factor(frame) || is.matrix(frame)))
      stop(deparse(substitute(frame))," must be a data frame, matrix, vector, or factor")

    if(!is.matrix(frame))
      t.names <- deparse(eval(substitute(substitute(frame)), env = parent.frame()))
    
    frame <- data.frame(frame, check.names=FALSE)
    
    if(!is.null(t.names))
      names(frame) <- t.names
  }
  return(frame)
}