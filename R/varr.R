`varr` <-
function(x, meanx) {
	n = ncol(x)
	#p = nrow(x)
	z = matrix(1, nrow=n, ncol=1)
	d = x - meanx %*% t(z)
	ans = (d^2) %*% rep(1/(n-1), n)
	ans = drop(ans)
	return (ans)
}

