\name{dri.impute}
\alias{dri.impute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Impute missing values for gene expression data}
\description{
  Missing gene expression measurements are imputed by the K-nearest neighbors method
}
\usage{
dri.impute(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ data matrix whose missing values will be imputed }
}
\details{
The K-nearest neighbors method is employed (K=10) to impute missing gene expression values. For any row
with more than 50% missing data, the values are imputed using the overall mean per sample.

Uses the \code{impute.knn} function of the \code{impute} package.
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{data}{ the new imputed data matrix }
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/

Olga Troyanskaya, Michael Cantor, Gavin Sherlock, Pat Brown, Trevor Hastie, Robert Tibshirani,
David Botstein and Russ B. Altman, Missing value estimation methods for DNA microarrays
BIOINFORMATICS Vol. 17 no. 6, 2001 Pages 520-525.
}
\author{ Keyan Salari, Robert Tibshirani, and Jonathan R. Pollack }
\note{
}
\seealso{
\code{\link{drcorrelate}},\code{\link{drcorrelate.null}},\code{\link{drsam}},\code{\link{drsam.null}},
\code{\link{dri.fdrCutoff}},\code{\link{dri.sig_genes}},\code{\link{dri.heatmap}}
}
\examples{
require(impute)
data(mySampleData)
attach(mySampleData)

# DNA data should contain no missing values - pre-mooth beforehand
# Impute missing values for gene expression data
RNA.data <- dri.impute(RNA.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
