\name{dri.fix.missingvals}
\alias{dri.fix.missingvals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Excel helper function for fixing missing values representation }
\description{
  Helper function for MS Excel DR-Integrator plug-in. Missing values are converted to
	a form interpretable by R
}
\usage{
dri.fix.missingvals(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{data matrix sent from Excel with missing values as -999999}
}
\details{
}
\value{
  \item{d}{data matrix with missing values replaced by NA}
}
\references{ 
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\author{ Keyan Salari, Robert Tibshirani, and Jonathan R. Pollack  }
\note{
}
\seealso{}
\examples{
# Excel helper function
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
