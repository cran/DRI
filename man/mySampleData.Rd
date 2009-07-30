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
 $ GeneIDs  : Factor w/ 100 levels "_ACTRT2","_AGRN",..: 54 55 50 58 100 44 62 51 52 56 ...
 $ GeneNames: Factor w/ 100 levels "_ACTRT2","_AGRN",..: 54 55 50 58 100 44 62 51 52 56 ...
 $ Chr      : int [1:100] 1 1 1 1 1 1 1 1 1 1 ...
 $ Nuc      : int [1:100] 20463 79326 127660 225560 248946 307582 307582 316075 316090 342167 ...
 $ DNA.data : num [1:100, 1:50] 0 0 0 0 0 0 0 0 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:100] "1" "2" "3" "4" ...
  .. ..$ : chr [1:50] "X184A1" "BT20" "BT474" "BT483" ...
 $ RNA.data : num [1:100, 1:50]     NA -0.791 -1.968     NA  0.319 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:100] "1" "2" "3" "4" ...
  .. ..$ : chr [1:50] "X184A1.1" "BT20.1" "BT474.1" "BT483.1" ...
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
