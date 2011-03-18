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
License: \tab GPL-3\cr
Author: \tab Noah Hoffman <ngh2@uw.edu>\cr 
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
\code{\link[stats]{cmdscale}}, \code{\link{clstutils}}
}
\examples{
library(clst)
packageDescription("clst")
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
