\name{plotDistances}
\alias{plotDistances}
\alias{plotMutinfo}

\title{
Visualize results of \code{link{findThreshold}}
}

\description{

  The functions \code{plotDistances} and \code{plotMutinfo} are used to
visualize the distance threshold calculated by
\code{\link{findThreshold}} in the context of pairwise distances among
objects in the reference set.

}

\usage{
plotDistances(distances, D = NA, interval = NA,
              ylab = "distances", ...)

plotMutinfo(breaks, D = NA, interval = NA,
            xlab = "distance", ylab = "mutual information", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distances}{The \code{$distances} element of the output value of \code{\link{findThreshold}}}

  \item{breaks}{The \code{$breaks} element of the output value of \code{\link{findThreshold}}}
  
  \item{D}{The distance threshold}
  
  \item{interval}{The range of values over which candidiate values of
    PMMI are evaluated.}
  
  \item{xlab}{Label the x axis of the plot.}

  \item{ylab}{Label the y axis of the plot.}

  \item{\dots}{Additional arguments are passed to \code{\link[lattice]{bwplot}}
  (\code{plotDistances}) or \code{\link[lattice]{xyplot}}} (\code{plotMutinfo})
  
}
\details{
\code{plotDistances} produces a box-and-whisker plot contrasting within-
and between-group distances. \code{plotMutinfo} produces a plot of
cutpoints vs mutual information scores.
}
\value{
Returns a lattice grid object.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Noah Hoffman
}
% \note{
% %%  ~~further notes~~
% }

\seealso{
  \code{\link{findThreshold}}
}
\examples{
data(iris)
dmat <- as.matrix(dist(iris[,1:4], method="euclidean"))
groups <- iris$Species
thresh <- findThreshold(dmat, groups)
do.call(plotDistances, thresh)
do.call(plotMutinfo, thresh)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

