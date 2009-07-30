\name{dri.sig_genes}
\alias{dri.sig_genes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Determine significant genes for DR-Correlate or DR-SAM analysis }
\description{
Given a cutoff score corresponding to the desired FDR, the list of significant genes is generated, 
along with each gene's FDR.
}
\usage{
dri.sig_genes(cutoff, observed, null_dist, gene_id, gene_name, chr, nuc, bt = TRUE, method = "drcorrelate")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cutoff}{ cutoff score for significance from \code{dri.fdrCutoff} }
  \item{observed}{ vector of observed scores from \code{drcorrelate} or \code{drsam} }
  \item{null_dist}{ matrix of null data from \code{drcorrelate.null} or \code{drsam.null} }
  \item{gene_id}{ vector of gene IDs }
  \item{gene_name}{ vector of gene names }
  \item{chr}{ vector of gene chromosome locations }
  \item{nuc}{ vector of gene nucleotide positions }
  \item{bt}{ either \code{TRUE} or \code{FALSE} indicated whether a 2-tail test was performed }
  \item{method}{ analysis method used, either "drcorrelate" or "drsam" }
}
\details{
}
\value{
  \item{Results.SigGenes}{a list of significant genes, positive and negative, with gene-specific FDRs }
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\author{ Keyan Salari, Robert Tibshirani, and Jonathan R. Pollack }
\note{ 
}
\seealso{
\code{\link{drcorrelate}},\code{\link{drcorrelate.null}},\code{\link{drsam}},
\code{\link{drsam.null}},\code{\link{dri.fdrCutoff}},\code{\link{dri.sig_genes}},\code{\link{dri.heatmap}}
}
\examples{
require(impute)
data(mySampleData)
attach(mySampleData)

# DNA data should contain no missing values - pre-mooth beforehand
# Impute missing values for gene expression data
RNA.data <- dri.impute(RNA.data)

# DR-Correlate analysis to find genes with correlated DNA/RNA measurements
obs <- drcorrelate(DNA.data, RNA.data, method="pearson")
# generate null distribution for FDR calculation (10 permutations)
null <- drcorrelate.null(DNA.data, RNA.data, method="pearson", perm=10)
# identify the correlation cutoff corresponding to your desired FDR
n.cutoff <- dri.fdrCutoff(obs, null, targetFDR=0.05, bt=TRUE)
cutoff <- n.cutoff[2]
# retrieve all genes that are significant at the determined cutoff, and calculate gene-specific FDRs
Results <- dri.sig_genes(cutoff, obs, null, GeneIDs, GeneNames, Chr, Nuc, bt=TRUE, method="drcorrelate") 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
