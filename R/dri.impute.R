dri.impute <-
function(d) {
	return(impute.knn(d)$data)
}

