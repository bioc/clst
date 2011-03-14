source('unit/setup.R')

VERBOSE = FALSE

.setUp <- function(){
  if(VERBOSE){
    cat('\n')
  }
}

library(clst)
data(actino)
data(strep)

test_getThresh_01 <- function(){
  thresh <- clst:::getThresh()
  checkTrue(mode(thresh) == 'list')
  checkTrue(is.na(thresh$D))
}

test_getThresh_02 <- function(){
  D <- 1.0
  thresh <- clst:::getThresh(D=D)
  checkTrue(mode(thresh) == 'list')
  checkTrue(thresh$D==D)
  checkTrue(all(c('breaks','distances') %in% names(thresh)))
  checkTrue(is.null(thresh$breaks))
}

test_findPMMI_strep01 <- function(){

  dmat <- strep$dmat1
  groups <- strep$abbrev

  dists <- partition(dmat, groups, verbose=VERBOSE)
  clst:::findPMMI(dists$vals, dists$comparison, verbose=VERBOSE)  
}

test_findPMMI_strep02 <- function(){

  dmat <- strep$dmat1
  groups <- strep$abbrev

  dists <- partition(dmat, groups, verbose=VERBOSE)
  clst:::findPMMI(dists$vals, dists$comparison, prob=0.5, verbose=VERBOSE)  
}

test_findThreshold_strep01 <- function(){

  dmat <- strep$dmat1
  groups <- strep$abbrev

  thresh1 <- findThreshold(dmat, groups, verbose=VERBOSE)  
  dists <- partition(dmat, groups, verbose=VERBOSE)
  thresh2 <- findThreshold(distances=dists, verbose=VERBOSE)  

  checkTrue(identical(thresh1, thresh2))
}

test_plotDistances_strep01 <- function(){
  dmat <- strep$dmat1
  groups <- strep$abbrev
  distances <- partition(dmat, groups, verbose=VERBOSE)

  pdf(print(pdfName()))
  plot(plotDistances(distances))
  dev.off()  
}

test_plotDistances_strep02 <- function(){
  dmat <- strep$dmat1
  groups <- strep$abbrev
  thresh <- findThreshold(dmat, groups, verbose=VERBOSE)

  pdf(print(pdfName()))
  with(thresh,{
    plot(plotDistances(distances, D=D))
    plot(plotDistances(distances, interval=interval))
    plot(plotDistances(distances, D=D, interval=interval))
  })
  dev.off()
}


test_plotMutinfo_strep01 <- function(){
  dmat <- strep$dmat1
  groups <- strep$abbrev
  thresh <- findThreshold(dmat, groups, verbose=VERBOSE)

  pdf(print(pdfName()))
  with(thresh,{
    plot(plotMutinfo(breaks))
    plot(plotMutinfo(breaks, D=D))
    plot(plotMutinfo(breaks, interval=interval))
    plot(plotMutinfo(breaks, D=D, interval=interval))
  })
  dev.off()
}
