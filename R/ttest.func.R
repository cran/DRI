ttest.func <-
function(x,y,s0=0, s=NULL){
	n1 <- sum(y==1)
	n2 <- sum(y==2)
	
	p <- nrow(x)
	m1 <- rowMeans(x[,y==1,drop=F])
	m2 <- rowMeans(x[,y==2,drop=F])


	if(is.null(s)){
		s <- sqrt( ((n2-1) * varr(x[, y==2], meanx=m2) + (n1-1) * varr(x[, y==1], meanx=m1) )*(1/n1+1/n2)/(n1+n2-2) ) 
	}

	numer <-  m2 - m1

	dif.obs <- (numer)/(s + s0)
  	return(list(tt=dif.obs,numer=numer, sd=s))
}

