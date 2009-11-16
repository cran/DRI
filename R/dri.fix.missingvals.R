dri.fix.missingvals <-
function(d) {
	d[d == -999999] <- NA
	return (d)
}

