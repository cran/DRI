\name{dri.merge.CNbyRNA}
\alias{dri.merge.CNbyRNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge DNA copy number and gene expression data sets}
\description{
 A matrix of DNA copy number values and a matrix of gene expression values are merged into
	one matrix using the genome position coordinates of the probes on each microarray platform.
	Copy number probes flanking each gene expression probe are averaged to derive the copy number
	at a given expression probe.
}
\usage{
dri.merge.CNbyRNA(dna.chr, dna.nuc, dna.data, rna.chr, rna.nuc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dna.chr}{ vector of chromosomes on which copy number probes are located }
  \item{dna.nuc}{ numeric vector of nucleotide positions where copy number probes are located }
  \item{dna.data}{ matrix of DNA copy number data }
  \item{rna.chr}{ vector of chromosomes on which expression probes are located }
  \item{rna.nuc}{ numeric vector of nucleotide positions where expression probes are located }
}
\details{
This function performs a pre-processing step to obtain 1-to-1 mappings of DNA copy number 
and gene expression for each measured gene. For each gene expression probe, the closest 5'
and 3' copy number probes are used to calculate an average copy number value. If the same
array platform was used to measure copy number and gene expression, the copy number and 
gene expression measurement will be matched by the probes' genome position coordinates.
This function returns a DNA copy number matrix that has the same number of rows (genes) as
the gene expression matrix.
}
\value{
  \item{average.dna.data }{ a matrix of average DNA copy number values corresponding to each
					gene expression value }
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
require(impute)
data(mySampleData)
attach(mySampleData)

dri.merge.CNbyRNA(dna.chr=Chr, dna.nuc=Nuc, dna.data=DNA.data, rna.chr=Chr, 
rna.nuc=Nuc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
