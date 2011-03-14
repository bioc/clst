\name{findThreshold}
\alias{findThreshold}
\alias{partition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
findThreshold
}

\description{ Identify a distance threshold predicting whether a
  pairwise distance represents a comparison between objects in the same
  class (within-group comparison) or different classes (between-group
  comparison) given a matrix providing distances between objects and the
  group membership of each object.  }

\usage{

findThreshold(dmat, groups, distances, method = "mutinfo", prob = 0.5,
              na.rm = FALSE, keep.dists = TRUE, roundCuts = 2, minCuts =
              20, maxCuts = 300, targetCuts = 100, verbose = FALSE,
              depth = 1, ...)

partition(dmat, groups, include, verbose = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{Square matrix of pairwise distances.}

  \item{groups}{Object coercible to a factor identifying group
    membership of objects corresponding to either edge of \code{dmat}.}

  \item{include}{vector (numeric or boolean) indicating which elements
  to retain in the output; comparisons including an excluded element
  will have a value of NA}
  
  \item{distances}{Optional output of \code{\link{partition}} provided in the
    place of \code{dmat} and \code{groups}}

  \item{method}{The method for calculating the threshold;  only
    'mutinfo' is currently implemented.}

  \item{prob}{Sets the upper and lower bounds of \var{D} as some quantile of
    the within class distances and between-class differences, respectively.}

  \item{na.rm}{If TRUE, excludes \code{NA} elements in \code{groups} and
    corresponding rows and columns in \code{dmat}. Ignored if
    \code{distances} is provided.}
  
  \item{keep.dists}{If TRUE, the output will contain the \code{distances}
  element (output of \code{\link{partition}}).}
  
  \item{roundCuts}{Number of digits to round cutoff values (see Details)}

  \item{minCuts}{Minimal length of vector of cutoffs (see Details).}
  
  \item{maxCuts}{Maximal length of vector of cutoffs (see Details)}
  
  \item{targetCuts}{Length of vector of cutoffs if conditions met by
    \code{minCuts} and \code{maxCuts} are not met (see Details).}
  
  \item{verbose}{Terminal output is produced if TRUE.}

  \item{depth}{Private argument used to track level of recursion.}

  \item{\dots}{Extra arguments are ignored.}
}

\details{
\code{findThreshold} is used internally in \code{\link{classify}}, but
may also be used to calculate a starting value of $D$.

\code{partition} is used to transform a square (or lower triangular)
   distance matrix into a \code{data.frame} containing a column of
   distances (\code{$vals}) along with a factor (\code{$comparison})
   defining each distance as a within- or between-group
   comparison. Columns \code{$row} and \code{$col} provide indices of
   corresponding rows and columns of \code{dmat}.

 }

\value{

In the case of \code{findThreshold}, output is a list with elements
decsribed below. In the case of \code{partition}, output is the \code{data.frame}
returned as the element named \code{$distances} in the output of
\code{findThreshold}.

\item{D}{The distance threshold (distance cutoff corresponding to the
  PMMI).}
\item{pmmi}{Value of the point of maximal mutual information (PMMI)}
\item{interval}{A vector of length 2 indicating the upper and lower
  bounds over which values for the threshold are evaluated.}
\item{breaks}{A \code{data.frame} with columns \code{x} and \code{y}
  providing candidiate breakpoints and corresponding mutual information
  values, respectively.}
\item{distances}{If \code{keep.distances} is TRUE, a data.frame
  containing pairwise distances identified as within- or between classes.}
\item{method}{Character corresponding to input argument \code{method}.}
\item{params}{Additional input parameters.}

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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{plotDistances}}, \code{\link{plotMutinfo}}}

\examples{
data(iris)
dmat <- as.matrix(dist(iris[,1:4], method="euclidean"))
groups <- iris$Species
thresh <- findThreshold(dmat, groups, type="mutinfo")
str(thresh)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

