\name{drcorrelate}
\alias{drcorrelate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute correlations between DNA copy number and gene expression for each gene}
\description{
  A correlation is computed between each gene's DNA copy number and gene expression across
	all the samples. Significant correlations are determined by comparison to a null distribution
	derived from random permutations of the data.
}
\usage{
drcorrelate(DNA, RNA, method = "pearson", tail_p = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DNA}{ matrix of DNA copy number data }
  \item{RNA}{ matrix of gene expression data, samples (columns) in same order as DNA matrix }
  \item{method}{ correlation statistic, either "pearson", "spearman", or "ttest"}
  \item{tail_p}{ top/bottom percent of samples (with respect to the gene's copy number) to use
			 for extremes t-test groups; used only when method = "ttest" }
}
\details{
DR-Correlate aims to identify genes with expression changes explained by underlying CNAs. This tool 
performs an analysis to identify genes with statistically significant correlations between their DNA 
copy number and gene expression levels. Three options for the statistic to measure correlation are 
implemented: (1) Pearson's correlation; (2) Spearman's rank correlation; and (3) an "extremes" t-test. 
For Pearson's and Spearman's correlations, the respective correlation coefficient is computed for each 
gene. For the extremes t-test, a modified Student's t-test (Tusher et al. 2001) is computed for each 
gene, comparing gene expression levels of samples comprising the lowest and highest deciles with respect 
to DNA copy number. That is, for each gene the samples are rank ordered by DNA copy number and samples 
below the 10th percentile and above the 90th percentile form two groups whose means of gene expression 
are compared with a modified t-test. The percentile cutoff defining the two groups is user-adjustable.
}
\value{
  \item{observed}{ a vector of observed correlations for each gene }
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

# DNA data should contain no missing values - pre-smooth beforehand
# Impute missing values for gene expression data
RNA.data <- dri.impute(RNA.data)

# DR-Correlate analysis to find genes with correlated DNA/RNA measurements
obs <- drcorrelate(DNA.data, RNA.data, method="pearson")
# generate null distribution for FDR calculation (10 permutations)
null <- drcorrelate.null(DNA.data, RNA.data, method="pearson", perm=10)
# identify the correlation cutoff corresponding to your desired FDR
n.cutoff <- dri.fdrCutoff(obs, null, targetFDR=0.05, bt=TRUE)
cutoff <- n.cutoff[2]
# retrieve all genes that are significant at the determined cutoff, and
# calculate gene-specific FDRs
Results <- dri.sig_genes(cutoff, obs, null, GeneIDs, GeneNames, Chr, Nuc, 
bt=TRUE, method="drcorrelate") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
