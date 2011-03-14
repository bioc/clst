source('unit/setup.R')

library(clst)
data(actino)
data(strep)

test_partition_01 <- function(){

  dmat <- strep$dmat1
  groups <- strep$abbrev

  dists <- partition(dmat, groups, verbose=VERBOSE)
  checkTrue(all(c('vals','comparison') %in% colnames(dists)))
  checkTrue(identical(levels(dists$comparison), c('within','between')))
}
