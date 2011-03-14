\name{classify}
\alias{classify}
\alias{classifier}
\alias{classifyIter}
\alias{pull}
\alias{pullTab}

\title{
classify
}
\description{
Functions to perform classification by local similarity threshold.
}
\usage{

classify(dmat, groups, dvect, method = "mutinfo", minScore = 0.45,
         doffset = 0.5, dStart = NA, maxDepth = 10, minGroupSize = 2,
         objNames = names(dvect), keep.data = TRUE, ..., verbose =
         FALSE)

classifyIter(dmat, groupTab, dvect, dStart = NA, multiple = FALSE,
             keep.data = TRUE, ..., verbose = FALSE)

classifier(dmat, groups, dvect, method = 'mutinfo', minScore = 0.45,
           doffset = 0.5, dStart = NA, minGroupSize = 2,
           objNames = names(dvect), keep.data = TRUE, ..., verbose = FALSE,
           depth = 1)

pull(dmat, groups, index)

pullTab(dmat, groupTab, index)

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dmat}{Square matrix of pairwise distances.}

  \item{groups}{Object coercible to a factor identifying group
    membership of objects corresponding to either edge of \code{dmat}.}
  
  \item{groupTab}{a data.frame representing a taxonomy, with columns in
 increasing order of specificity from left to right (ie, Kingdom -->
 Species). Column names are used to name taxonomic ranks. Rows
 correspond to margins of dmat.}
  
  \item{dvect}{numeric vector of distance from query sequence to each
  reference corresponding to margins of dmat.}
  
  \item{method}{The method for calculating the threshold; only
    'mutinfo' is currently implemented.}

  \item{minScore}{Threshold value for the match score to define a match.}
  
  \item{doffset}{Offset used in the denominator of the expression to
  calculate match score to penalize very small groups of reference objects.}
  
  \item{dStart}{start with this value of \code{D}.}

  \item{multiple}{if TRUE, stops at the rank that yields at least one
  match; if FALSE, continues to perform classification until exactly one
  match is identified.}
  
  \item{maxDepth}{Maximum number of iterations that will be attempted to
  perform classification.}
  
  \item{minGroupSize}{The minimal number of members comprising at least one
  group required to attempt classification.}
  
  \item{objNames}{Optional character identifiers for objects
  corresponding to margin of \code{dmat}.}
  
  \item{keep.data}{Populates \code{thresh$distances} (see \code{\link{findThreshold}}) if TRUE.}
  
  \item{verbose}{Terminal output is produced if TRUE.}

  \item{index}{an integer specifying an element in \code{dmat}}
  
  \item{\dots}{see Details}

  \item{depth}{specifies iteration number (not meant to be user-defined)}
  
}
\details{

  \code{classify} performs iterative classification. See the vignette
  vignette for package \pkg{clst} for a description of the classification
  algorithm.
  
  \code{classifier} performs non-iterative classification, and is
typically not called directly by the user.

The functions \code{pull} and \code{pullTab} are used to remove a single
element of \code{dmat} for the purpose of performing classification
agains the remaining elements. The value of these two functions (a list)
can be passed directly to \code{classify} or \code{classifyIter}
directly (see examples).

}

\value{ \code{classify} and \code{classifyIter} return \code{x}, a list
  of lists, one for each iteration of the classifier. Each sub-list
  contains the following named elements:

  \item{depth}{An integer indicating the number of the iteration (where
  x[[i]]$depth == i)}
  
  \item{tally}{a \code{data.frame} with one row for each group or
  reference objects. Columns \code{below} and \code{above} contain
  counts of reference objects with distance values greater than or less
  than \emph{D}, respectively; \code{score}, containing match score
  \eqn{S}; \code{match} is 1 if \eqn{S \ge minScore}, 0 otherwise; and
  the minimum, median, and maximum values of distances to all members of
  the indicated group.}

  \item{details}{a list of two matrices, named "below" and "above",
  itemizing each object with index \emph{i} in the reference set with
  distances below or above the distance threshold \emph{D},
  respectively. Columns include \code{index}, the index \emph{i};
  \code{dist}, the distance between the object and the query; and
  \code{group}, indicating the classification of the object.}

  \item{matches}{Character vector naming groups to which query object
    belongs.}
  
  \item{thresh}{object returned by \code{\link{findThreshold}}}

  \item{params}{a list of input arguments and their
    values}
  
  \item{input}{list containing copies of \code{dvect} and \code{groups}}

}

% \references{
% %% ~put references to the literature/web site here ~
% }

\author{
Noah Hoffman
}

% \note{
  
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{findThreshold}}
}
\examples{

## illustrate classification using the Iris data set
data(iris)
dmat <- as.matrix(dist(iris[,1:4], method="euclidean"))
groups <- iris$Species

## remove one element from the data set and perform classification using
## the remaining elements as the reference set
ind <- 1
cat(paste('class of "unknown" sample is Iris',groups[ind]),fill=TRUE)
cc <- classify(dmat[-ind,-ind], groups[-ind], dvect=dmat[ind, -ind])
printClst(cc)

## this operation can be performed conveinetly using the `pull` function
ind <- 51
cat(paste('class of "unknown" sample is Iris',groups[ind]),fill=TRUE)
cc <- do.call(classify, pull(dmat, groups, ind)) 
printClst(cc)
str(cc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

