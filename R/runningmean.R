runningmean <-
function (x, k) {
	n <- length(x)
	k2 <- k %/% 2
	y <- vector(mode="numeric", length=n)

	for(j in 1:k2) {
		y[j]=mean(x[1:(j+k2)], na.rm=TRUE)
	}
	for(j in (1+k2):(n-k2)) {
		y[j]=mean(x[(j-k2):(j+k2)], na.rm=TRUE)
	}
	for(j in (n-k2+1):n) {
		y[j]=mean(x[(j-k2):n], na.rm=TRUE)
	}
	return(y)
}

