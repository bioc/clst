\name{clst-package}
\Rdversion{1.1}
\alias{clst-package}
\alias{clst}
\docType{package}
\title{
Classification by local similarity threshold
}

\description{Package for modified nearest-neighbor classification based
on calculation of a similarity threshold distinguishing within-group
from between-group comparisons.}

\details{
\tabular{ll}{
Package: \tab clst\cr
Type: \tab Package\cr
Version: \tab 0.1.5256\cr
Date: \tab 2011-03-14\cr
Depends: \tab R (>= 2.10.0), methods\cr
Imports: \tab ROCR, lattice\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
License: \tab Public domain\cr
Built: \tab R 2.12.2; ; 2011-03-14; unix\cr
}

Index:
\preformatted{
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{clstDemo} \tab clst (source, pdf)\cr
}

TODO: write package overview.

}
\author{
Noah Hoffman

Maintainer: <ngh2@uw.edu>
}
% \references{
% }
\keyword{ package }
\keyword{ classif }
\seealso{
\code{\link[stats]{cmdscale}}
}
\examples{
library(clst)
data(iris)
dmat <- as.matrix(dist(iris[,1:4], method="euclidean"))
groups <- iris$Species
i <- 1
cc <- classify(dmat, groups, dvect=dmat[i,])
cat('query at i =',i,'is species',paste('I.', groups[i]),'\n')
printClst(cc)
i <- 125
cc <- classify(dmat, groups, dvect=dmat[i,])
cat('query at i =',i,'is species',paste('I.', groups[i]),'\n')
printClst(cc)
}
