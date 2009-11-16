scale01 <-
function(x, low=min(x), high=max(x) ) {
	x <- (x-low)/(high - low)
	return(x)
}

