\name{DRI-package}
\alias{DRI-package}
\alias{DRI}
\docType{package}
\title{
DR-Integrator: an analytic tool for integrating DNA copy number and gene expression data
}
\description{
DR-Integrator identifies genes with significant correlations between DNA copy number
alterations and gene expression data, and implements a supervised learning analysis that 
captures genes with significant alterations in both DNA copy number and gene expression
between two sample classes.
}
\details{
\tabular{ll}{
Package: \tab DRI\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2009-11-16\cr
License: \tab GPL-2\cr
}
This package contains two analytic tools: DR-Correlate and DR-SAM.
}
\author{
Keyan Salari, Robert Tibshirani, Jonathan R. Pollack

Maintainer: Keyan Salari <ksalari@stanford.edu>
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}

\keyword{ package }
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

# Optional heatmap plot for significant DR-Correlation genes
sample.names <- colnames(DNA.data)
pdf(file="DRI-Heatmap.pdf", height=8, width=11)
dri.heatmap(Results, DNA.data, RNA.data, sample.names, GeneNames, Chr, Nuc, 
statistic="pearson", color.scheme="RG")
dev.off()

# DR-SAM analysis to find genes with alterations in both DNA and RNA between
# different classes
labels <- c(rep(1,25), rep(2,25)) # 25 samples in class 1 and 25 in class 2
obs <- drsam(DNA.data, RNA.data, labels, transform.type="raw")
# generate null distribution for FDR calculation (10 permutations)
null <- drsam.null(DNA.data, RNA.data, labels, transform.type="raw", 10)
# identify the correlation cutoff corresponding to your desired FDR
n.cutoff <- dri.fdrCutoff(obs$test.summed, null, targetFDR=0.05, bt=TRUE)
cutoff <- n.cutoff[2]
# retrieve all genes that are significant at the determined cutoff, and
# calculate gene-specific FDRs
Results <- dri.sig_genes(cutoff, obs, null, GeneIDs, GeneNames, Chr, Nuc, 
bt=TRUE, method="drsam")
}
