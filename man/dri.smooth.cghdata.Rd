\name{dri.smooth.cghdata}
\alias{dri.smooth.cghdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smooth DNA copy number data over a specific window size}
\description{
 A moving average is calculated by chromosomal location of copy number probes using a specified
 window size.
}
\usage{
dri.smooth.cghdata(DNA.data, Chr, mw = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DNA.data}{ matrix of DNA copy number data (rows ordered by genome position}
  \item{Chr}{ vector of chromosomes on which copy number probes are located (ordered by probes' genome position }
  \item{mw}{ window size to take moving average over}
}
\details{
This function performs one pre-processing method for smoothing copy number data before performing 
further DNA/RNA integrative analyses. To help eliminate missing values, and reduce noise, for each 
sample a moving average is taken over a user-specified window size along each chromosome.
}
\value{
  \item{DNA.smoothed}{ matrix of DNA copy number data smoothed by a moving average window}
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\author{ Keyan Salari, Robert Tibshirani, and Jonathan R. Pollack }
\note{
}
\seealso{
\code{\link{drcorrelate}}, \code{\link{drcorrelate.null}}, \code{\link{drsam}}, 
\code{\link{drsam.null}}, \code{\link{dri.fdrCutoff}}, \code{\link{dri.sig_genes}}, 
\code{\link{dri.heatmap}}, \code{\link{dri.merge.CNbyRNA}}, \code{\link{dri.smooth.cghdata}}, 
\code{\link{runFusedLasso}}
}
\examples{
data(mySampleData)
attach(mySampleData)

DRI.DNA.processed <- dri.smooth.cghdata(DNA.data, Chr, 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
