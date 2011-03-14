\name{scaleDistPlot}
\alias{scaleDistPlot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Annotated multidimensional scaling plots.
}

\description{ Produces annotated representations of two-dimensional
multidimensional scaling plots using \code{\link[stats]{cmdscale}}.  }

\usage{

scaleDistPlot(dmat, groups, fill, X, O, indices = "no",
              include, display, labels,
              shuffleGlyphs = NA, key = "top",
              keyCols = 4, glyphs,
              xflip = FALSE, yflip = FALSE, ...)

}

\arguments{

  \item{dmat}{Square matrix of pairwise distances.}

  \item{groups}{Object coercible to a factor identifying group
    membership of objects corresponding to either edge of \code{dmat}.}
  
  \item{fill}{vector (logical or indices) of points to fill}
  
  \item{X}{vector of points to mark with an X}
  
  \item{O}{vector of points to mark with a circle}
  
  \item{indices}{label points with indices (all points if 'yes', or a
  subset indicated by a vector)}
  
  \item{include}{boolean or numeric vector of elements to include in
  call to cmdscale}
  
  \item{display}{boolean or numeric vector of elements to include in
  call to display}
  
  \item{labels}{list or data frame with parameters $i indicating indices
  and $text containing labels.}
  
  \item{shuffleGlyphs}{modify permutation of shapes and colors given an
  integer to serve as a random seed.}
  
  \item{key}{'right' (single column), 'top' (variable number of
  columns), or NULL for no key}
  
  \item{keyCols}{number of columns in key}
  
  \item{glyphs}{a data.frame with columns named \code{col} and
  \code{pch} corresponding to elements of \code{unique(groups)}}
  
  \item{xflip}{if TRUE, flip orientation of x-axis}
  
  \item{yflip}{if TRUE, flip orientation of y-axis}
  
  \item{\dots}{additional arguments are passed to \code{\link[lattice]{xyplot}}}
  
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[stats]{cmdscale}}, \code{\link[lattice]{xyplot}}
}
\examples{
data(iris)
dmat <- as.matrix(dist(iris[,1:4], method="euclidean"))
groups <- iris$Species

## visualize pairwise euclidean dstances among items in the Iris data set
fig <- scaleDistPlot(dmat, groups)
plot(fig)

## leave-one-out analysis of the classifier
loo <- lapply(seq_along(groups), function(i){
  do.call(classify, pull(dmat, groups, i))
})
matches <- lapply(loo, function(x) rev(x)[[1]]$matches)
result <- sapply(matches, paste, collapse='-')
confusion <- sapply(matches, length) > 1
no_match <- sapply(matches, length) < 1
plot(scaleDistPlot(dmat, groups, fill=confusion, O=confusion, X=no_match))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

