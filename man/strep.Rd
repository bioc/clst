\name{strep}
\Rdversion{1.1}
\alias{strep}
\docType{data}
\title{
Streptococcus data set.
}
\description{
Square matrices decsribing pairwise distances among 16s rRNA sequences.
}
\usage{data(strep)}
\format{\preformatted{
List of 5
 $ dmat1 : num [1:150, 1:150] 0 5.81 8.38 10.28 10.64 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:150] "197" "199" "207" "208" ...
  .. ..$ : chr [1:150] "197" "199" "207" "208" ...
 $ dmat2 : num [1:150, 1:150] 0 5.09 3.82 7.21 7.59 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:150] "197" "199" "207" "208" ...
  .. ..$ : chr [1:150] "197" "199" "207" "208" ...
 $ dmat3 : num [1:150, 1:150] 0 5.63 5.81 8.77 9.14 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:150] "197" "199" "207" "208" ...
  .. ..$ : chr [1:150] "197" "199" "207" "208" ...
 $ taxa  : Factor w/ 50 levels "Streptococcus acidominimus",..: 31 44 26 4 4 31 32 39 42 31 ...
 $ abbrev: Factor w/ 50 levels "S acidominimus",..: 31 44 26 4 4 31 32 39 42 31 ...
}}
\details{
The matrices \code{$dmat1}, \code{dmat2}, and \code{dmat3} contain
percent nucleotide difference with indels penalized heavily, little, and
somewhat, respectively.


\code{$taxa} is a factor of species names; abbreviations of the same
names are found in \code{$abbrev}.

}

\examples{
data(strep)
}
\keyword{datasets}
