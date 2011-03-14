\name{printClst}
\alias{printClst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a summary of the classifier output.
}
\description{
Prints a description of the output of \code{classify}.
}
\usage{
printClst(cc, rows = 8, nameWidth = 30, groupNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cc}{Output of \code{\link{classify}}}
  
  \item{rows}{Number of rows corresponding to groups of reference
  objects to show.}
  
  \item{nameWidth}{Character width of group names.}
  
  \item{groupNames}{a named vector containing replacement names for
  groups keyed by categories in \code{groups} (\code{\link{classify}})
  or \code{groupTab} (\code{\link{classifyIter}}).}
  
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
Output value is \code{NULL}; output is to stdout.
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
\code{\link{classify}}, \code{\link{classifyIter}}
}
\examples{
data(iris)
dmat <- as.matrix(dist(iris[,1:4], method="euclidean"))
groups <- iris$Species

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
