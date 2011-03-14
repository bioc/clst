source('unit/setup.R')

VERBOSE = FALSE

.setUp <- function(){
  if(VERBOSE){
    cat('\n')
  }
}

library(clst)
data(bvseqs)
data(actino)
data(strep)

test_strep01 <- function(){
  i <- 1
  pp <- clst:::pull(strep$dmat1, strep$abbrev, index=i)
  cc <- with(pp, classify(dmat, groups, dvect, verbose=VERBOSE))
  printClst(cc)
}

test_bvseqs01 <- function(){
  i <- 1
  pp <- clst:::pullTab(bvseqs$dmat, bvseqs$groupTab, index=i)
  cc <- with(pp, classifyIter(dmat, groupTab, dvect,
                              verbose=VERBOSE))

  printClst(cc)
  printClst(cc,groupNames=bvseqs$taxNames)
  printClst(cc,nameWidth=15,groupNames=bvseqs$taxNames)

}

## test_bvseqs02 <- function(){

##   ## find Dgenus
##   thresh <- with(bvseqs, findThreshold(dmat, groupTab$order, na.rm=TRUE,
##                                        verbose=VERBOSE))


##   i <- 1
##   pp <- clst:::pullTab(bvseqs$dmat, bvseqs$groupTab, index=i)
##   cc <- with(pp, classifyIter(dmat, groupTab, dvect,
##                               dStart=thresh$D,
##                               verbose=VERBOSE))
  
##   pdf(print(pdfName()))
##   plot.cc(cc)
##   dev.off()  
## }
