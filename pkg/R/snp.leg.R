# -------------------------------------------------------------------------   
# author			: konrad.smolinski@gmail.com
# date			: 15 May 2007     
# current vs	: 15 Jan 2008
# description	: snp-bierenshian pdf & cdf; Legandre polynomials 
#					: delta parametrization implemented; cleaner and shorter
#					: improvement on computation (# of loops reduced);         
# last upgrade : corrected for normalization
# -------------------------------------------------------------------------
# --- eq.(22) Bierens(2006) ----------------------------------------------
psnp.leg <- function(u=c(0.25,0.75),delta=seq(1,len=5)){  
  lu <- length(u) 
  dim <- dim(as.matrix(delta))
  ld <- dim[1]  													 	# global for Pimat 
  if(dim[2] == 1) { delta <- matrix(delta,dim[1],lu) } 	# work on this line

  Pimatfun <- function(i)  Pimat(u[i],ld+1) 	  				# prepare for next line
  Pimatu <- sapply(seq(1,lu), Pimatfun )   					# compute for all u's at once
  L <- Lmat(ld) 

  fun <- function(i){                                   	
		eta <- c(1,delta[,i])
		  P <- matrix(Pimatu[,i],(ld+1),(ld+1))  
		nom <- t(eta)%*%L%*%P%*%t(L)%*%eta   	
  return(nom/sum(eta^2))
  }
  cdf <- sapply(seq(1:lu),fun)
return(cdf)
}    
#---- # eq.(17) Bierens(2006) ------------------------------------------
# delta representation ; uses function rho
dsnp.leg <- function(u=c(0.25,0.75),delta=seq(1,len=5)){
	lu <- length(u)
 	dim <- dim(as.matrix(delta))    
 	ifelse(dim[2] == 1, delta <- matrix(c(1,delta),dim[1]+1,lu), delta <- as.matrix(rbind(1,delta) ))     
	
	rhom <- rho(u,dim[1])      
 	pdf <- apply(rhom*delta,2,sum)^2/apply(as.matrix(delta^2),2,sum)   
return(pdf)
}               
#-----------------------------------------------------------------------
#---- # eq.(7) Bierens(2006) -------------------------------------------
# recursive polynomials for pdf 
rho <- function(u=c(0.2,0.75),n=1){
	rho <- matrix(0,nrow=n+1,ncol=length(u))
	rho[1,] <- 1 ;	rho[2,] <- sqrt(3)*(2*u-1)

	fun <- function(i){  
		k <- i-1
  		a <- sqrt(4*k^2-1)/k
  		b <-  (k-1)*sqrt(2*k+1)/( k*sqrt(2*k-3) )
  		rho[i,] <<- a*(2*u-1)*rho[i-1,] - b*rho[i-2,]  # need to be <<-
  	}
  	if(n>1)	 sapply(seq(3,(n+1)),fun) 
return(rho)
}  
#---- # eq.(21 and forward) Bierens(2006) ------------------------------ 
# L matrix - lower triangular; 
# construct: L.ext - extended by 1st column of zeros
# compute recursively starting from 2nd col of L.ext.  
# at the end: extract important component of L.ext
# constant comp: note the index hift
Lmat <- function(n){
	L.ext <- matrix(0,n+1,n+2)	
	L.ext[1,2] <- 1 ; L.ext[2,2] <- -sqrt(3) ;  L.ext[2,3] <- 2*sqrt(3)	

   rec <- function(m,k){   
#   	a <- sqrt(4*m^2 - 1)/m
#   	b <- ( (m-1)*sqrt(2*m+1) )/( m*sqrt(2*m-3)) 
		
   	a <- sqrt(4*(m-1)^2 - 1)/(m-1)         # note on indexing m
   	b <- ( (m-2)*sqrt(2*(m-1)+1) )/( (m-1)*sqrt(2*(m-1)-3)) 
 		L.ext[m,k] <<-  a*(2*L.ext[m-1,k-1] - L.ext[m-1,k]) - b*L.ext[m-2,k]
# 		assign("L.ext[m,k]", a*(2*L.ext[m-1,k-1] - L.ext[m-1,k]) - b*L.ext[m-2,k], envir=environment(Lmat) )
   }
   if(n>1){    							# otherwise subscript out of bounds                          
		for(j in 2:(n+2)){
			for(i in 3:(n+1)){  if(j<=(i+1)) L.ext[i,j] <- rec(i,j)  }  
		}
	}
return(L <- L.ext[,2:(n+2)])
}     
#---  eq.(one before 22) Bierens(2006) --------------------------------- 
# here is smart trick with matrix indexing - note: computed for u unidimentional
Pimat <- function(u=c(0.25,0.75),n=1){ 
	m <- matrix(0,nrow=n,ncol=n)
	m <- row(m)+col(m) -1    # note: in definition of Pi indexes start from i,j=0
return((u^m)/m)
}    
# ----------------------------------------------------------------------
