\name{runFusedLasso}
\alias{runFusedLasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Wrapper function to run cghFLasso }
\description{
 A function to run the Fused Lasso method on DNA copy number data to call copy number alterations.
}
\usage{
runFusedLasso(DNA.data, normal.data = NA, chr, nuc, FDR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DNA.data}{ matrix of disease tissue DNA copy number data}
  \item{normal.data}{ matrix of normal tissue DNA copy number data (optional) }
  \item{chr}{ vector of chromosomes where copy number probes are located (ordered by genome location) }
  \item{nuc}{ numeric vector of nucleotide positions where copy number probes are located (ordered by genome location) }
  \item{FDR}{ false discovery rate to apply to threshold significant gains and losses }
}
\details{
This function performs a pre-processing step to smooth and segment DNA copy number data using the fused
lasso method of Tibshirani & Wang, 2008. The DNA copy number data from the disease tissues are compared 
to that of normal tissues, if supplied, to help call significant copy number alterations. A false discovery
rate is calculated for each copy number alteration and the significance threshold is user-definable. More
details of the Fused Lasso method can be found with the \code{cghFLasso} package and in the reference below.
}
\value{
  \item{data.FL}{matrix of DNA copy number data with significant gains and losses called by Fused Lasso,
			and non-significant changes set to zero}
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/

Tibshirani, R., and Wang, P. (2008) Spatial smoothing and hot spot detection for CGH data using the 
fused lasso, Biostatistics, 9, 18-29.
}
\author{
Keyan Salari, Robert Tibshirani, and Jonathan R. Pollack }
\note{
}
\seealso{
\code{\link{drcorrelate}}, \code{\link{drcorrelate.null}}, \code{\link{drsam}}, 
\code{\link{drsam.null}}, \code{\link{dri.fdrCutoff}}, \code{\link{dri.sig_genes}}, 
\code{\link{dri.heatmap}}, \code{\link{dri.merge.CNbyRNA}}, \code{\link{dri.smooth.cghdata}}, 
\code{\link{runFusedLasso}}
}
\examples{
library(cghFLasso)
data(CGH)
attach(CGH)

DRI.DNA.data.FL <- runFusedLasso(DNA.data=DiseaseArray, 
normal.data=NormalArray, chr=chromosome, nuc=nucposition, FDR=0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
