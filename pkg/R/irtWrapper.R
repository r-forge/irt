# This file contains irt Wrapper and Item Constructors.
# NOTE:
# Item Constructors do not support WI and Wsig yet (basic extension soon). 
# However it does support factors (as in the RHS note: item is a matrix that consis of
#	an item, covariates(design) matrix and item determining factors).
# -------------------------------------------------------------------
irt <- function(...){
    args = list(...)
    if( any( c("item","multipleitem") == class(args[[1]]) ) ) UseMethod("irt",args[[1]])
    else stop("First argument is not an item (construct item via item() or multipleitem()")	
}

irt.item <- function(itemObj,...){
	r <- itemObj$item
#	W <- itemObj$dmat
	W <- as.matrix(itemObj$dmat)
	irtSingle(r=r,W=W,...)
}

irt.multipleitem <- function(itemObj,...){
	r <- as.matrix(itemObj$item)
	W <- as.matrix(itemObj$dmat)
	irtMulti(r,W,...)
}

# -------------------------------------------------------------------
# Single Item Control - for irtSingle()
# -------------------------------------------------------------------
item <- function(...){
    args = list(...)
    if( any( c("formula","factor","integer","numeric") == class(args[[1]]) ) ) UseMethod("item",args[[1]])
    else stop(paste("item class '",class(args[[1]]),"' is not allowed."),sep="")
  }
# -------------------------------------------------------------------
item.formula <- function(formula,formula2=NULL, data = parent.frame(),itemnm="item",...){
	
		call <- match.call( expand.dots = FALSE ) 
		m <- match("data", names( call ), 0 ) 				# model
		mf <- call[c(1,m)]						 			# model frame
		mf[[1]] <- as.name( "model.frame" )

		mf$formula <- formula
		emf <- eval.parent(mf)				# list of evaluated model frames of each equation
		terms <- attr( emf, "terms" ) 		# terms and model frames for the individual equations
		item <- model.extract( emf, "response" )
				
		factors <- NULL

		if( length(attr(terms,"term.labels") ) != 0){	# if design matrix exists
			attr(terms,"intercept") <- 0
			W <- model.matrix( terms, emf )
		}else{											# if desig matrix does NOT exist
			attr(terms,"intercept") <- 0
			W <- rep(1,length(item))			
		}


		if(class(formula2) == "formula"){
			mf$formula <- formula2
			emf <- eval.parent(mf)			
			terms <- attr( emf, "terms" ) 
			attr(terms,"intercept") <- 0
			factors <- model.matrix( terms, emf )			
		}
		attr(W,"assign") <- NULL
		attr(factors,"assign") <- NULL

		dmatNames <- colnames(W)
		factorsNames <- colnames(factors)
		#itemnm <- "item"
		itm <- list(item=item,dmat=W,factors=factors,dmatnm=dmatNames,factorsnm=factorsNames,itemnm=itemnm)
		itm <- itemValidate(itm) 

		class(itm) <- "item" 
return(itm)
}
# -------------------------------------------------------------------
item.default <- function(item,designmat=NULL,factors=NULL,itemnm="item",...){
	if(is.null(designmat)) designmat <- rep(1,length(item))
	dmatNames <- colnames(designmat)
	factorsNames <- colnames(factors)	
	itm <- list(item=item,dmat=designmat,factors=factors,dmatnm=dmatNames,factorsnm=factorsNames,itemnm=itemnm) 
	itm <- itemValidate(itm) 

	class(itm) <- "item" 
return(itm)		
}
# -------------------------------------------------------------------
itemValidate <- function(obj){
	item <- as.integer(obj$item) 		# we forse item to be integer! In the future do it for factor.
	lnitem <- length(obj$item)

	if( sum( obj$item == item ) == lnitem ) obj$item <- item
	else 	stop("Item is not of integer type!")

	if(!is.null(obj$designmat)) 
		if(lnitem != nrow(obj$designmat) ) 	stop("length of the item vec different from nrow of Design Mat")

	if(!is.null(obj$factors)) 
		if(lnitem != nrow(obj$factors) ) 	stop("length of the item vec different from nrow of Factors Mat")

return(obj)	
}
# -------------------------------------------------------------------
# Multiple Item Control for irtMulti()
# it returns a matrix of items and matrix of covariates (W). Likewise input in the irtM1() model.
# -------------------------------------------------------------------
multipleitem <- function(...){
    args = list(...)
    if( any( c("matrix","data.frame","numeric") == class(args[[1]]) ) ) UseMethod("multipleitem",args[[1]])
    else stop(paste("multipleitem class '",class(args[[1]]),"' is not allowed."),sep="")
  }
# -------------------------------------------------------------------
multipleitem.default <- function(multipleitem,designmat=NULL,factors=NULL,itemnm=NULL,...){
#	multipleitem <- as.matrix(multipleitem)
	multipleitem <- as.data.frame(multipleitem)

	if(is.null(designmat)) designmat <- rep(1,nrow(multipleitem))
	dmatNames <- colnames(designmat)
	factorsNames <- colnames(factors)	
	if(is.null(itemnm)) itemnm <- colnames(multipleitem)
	names(multipleitem) <- itemnm
	itm <- list(item=multipleitem,dmat=designmat,factors=factors,dmatnm=dmatNames,factorsnm=factorsNames,itemnm=itemnm) 
	itm <- multipleitemValidate(itm) 

	class(itm) <- "multipleitem" 
return(itm)		
}
# -------------------------------------------------------------------
multipleitemValidate <- function(obj){
	nitems <- ncol(obj$item)
	lnitem <- nrow(obj$item)

	for(i in 1:nitems){
		item <- as.integer(obj$item[,i]) 		# we forse item to be integer! In the future do it for factor.
		if( sum( obj$item[,i] == item ) == lnitem ) obj$item[,i] <- item
		else 	stop("Item is not of integer type!")	
	}
	if(!is.null(obj$designmat)) 
		if(lnitem != nrow(obj$designmat) ) 	stop("length of the item vec different from nrow of Design Mat")
	
	if(!is.null(obj$factors)) 
		if(lnitem != nrow(obj$factors) ) 	stop("length of the item vec different from nrow of Factors Mat")
	
return(obj)	
}
# -------------------------------------------------------------------



