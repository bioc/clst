\name{actino}
\Rdversion{1.1}
\alias{actino}
\docType{data}
\title{
Actinomyces data set
}
\description{
Square matrices decsribing pairwise distances among 16s rRNA sequences. 
}
\usage{data(actino)}
\format{\preformatted{
List of 5
 $ dmat1 : num [1:146, 1:146] 0 0.763 1.25 10.345 12.771 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:146] "200" "201" "202" "203" ...
  .. ..$ : chr [1:146] "200" "201" "202" "203" ...
 $ dmat2 : num [1:146, 1:146] 0 0.574 1.044 5.669 8.409 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:146] "200" "201" "202" "203" ...
  .. ..$ : chr [1:146] "200" "201" "202" "203" ...
 $ dmat3 : num [1:146, 1:146] 0 0.763 1.25 8.571 11.233 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:146] "200" "201" "202" "203" ...
  .. ..$ : chr [1:146] "200" "201" "202" "203" ...
 $ taxa  : Factor w/ 33 levels "Actinomyces bowdenii",..: 12 12 12 23 20 20 8 22 12 20 ...
 $ abbrev: Factor w/ 33 levels "A bowdenii","A canis",..: 12 12 12 23 20 20 8 22 12 20 ...
}}
\details{
The matrices \code{$dmat1}, \code{dmat2}, and \code{dmat3} contain
percent nucleotide difference with indels penalized heavily, little, and
somewhat, respectively.


\code{$taxa} is a factor of species names; abbreviations of the same
names are found in \code{$abbrev}.

}

\examples{
data(actino)
}
\keyword{datasets}
