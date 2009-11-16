\name{mySampleData}
\alias{mySampleData}
\docType{data}
\title{ DR-Integrator Sample Data}
\description{
  Sample data of 100 genes copy number and gene expression data from 50 cell lines
}
\usage{data(mySampleData)}
\format{
  The format is:
List of 6
 $ GeneIDs  : Factor
 $ GeneNames: Factor
 $ Chr      : int
 $ Nuc      : int
 $ DNA.data : matrix
 $ RNA.data : matrix
}
\details{
}
\source{
Kao, J., Salari, K., Bocanegra, M.C. et al. (2009) Molecular profiling of breast cancer cell lines defines
relevant tumor models and provides a resource for cancer gene discovery. PLoS One.
}
\references{
Salari, K., Tibshirani, R., and Pollack, J.R. (2009) DR-Integrator: a new analytic tool for integrating 
DNA copy number and gene expression data. http://pollacklab.stanford.edu/
}
\examples{
data(mySampleData)
}
\keyword{datasets}
