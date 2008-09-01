# date: 17 May 2007
# last modified: Nov 4th 2007
#	done: cleaned; multiple delta added; apply function used;
# sugg. changes: 
#		dsnp,psnp: delta creation - speed up
#		psnp: eliminate loops
# ----------------------------------------------------------
dsnp.sine <- function(u=c(0.25,0.75),delta=seq(1,len=4)){
	lu <- length(u) 
	dim <- dim(as.matrix(delta))
	if(dim[2] == 1) { delta <- matrix(delta,dim[1],lu) } # work on this line
	kvec <- seq(1, len=dim[1])
	
	fun <- function(i){
		kvec1 <- cos(kvec*pi*u[i])
		nom <- (1 + sqrt(2)*sum( delta[,i]%*%kvec1) )^2
		return(  nom/(1 + sum(delta[,i]^2) )  )
	}
	pdf <- sapply(seq(1:lu),fun) 
return(pdf)
}
# ----------------------------------------------------------
psnp.sine <- function(u=c(0.25,0.75),delta=seq(1,len=5)){
  lu <- length(u) 
  dim <- dim(as.matrix(delta))
  if(dim[2] == 1) { delta <- matrix(delta,dim[1],lu) } # work on this line
  kvec <- seq(1, len=dim[1])
	
  fun <- function(i){
	comp1 <- (sin(u[i]*pi*kvec))/(pi*kvec)
  	comp2 <- (sin(2*u[i]*pi*kvec))/(2*pi*kvec)
  	comp3 <- 0

	for(m in 1:(dim[1]-1))	for(k in (m+1):dim[1]) 
   	  comp3 <- comp3 + delta[m,i]*delta[k,i]*( (sin((k+m)*pi*u[i]))/(pi*(k+m)) + (sin((k-m)*pi*u[i]))/(pi*(k-m))  )
			
  return( u[i] + (2*sqrt(2)*sum(comp1*delta[,i]) + sum(comp2*(delta[,i]^2) ) + comp3 )/(1 + sum(delta[,i]^2))  )
  }
  cdf <- sapply(seq(1:lu),fun) 
return(cdf)
}
# ----------------------------------------------------------