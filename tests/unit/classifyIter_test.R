source('unit/setup.R')

VERBOSE = TRUE

.setUp <- function(){
  if(VERBOSE){
    cat('\n')
  }
}

library(clst)
data(bvseqs)

plot.cc <- function(cc){
  for(ci in cc){
    with(ci$thresh, {
      plot.new()      
      ff <- plotDistances(distances, D, interval)
      plot(ff, split=c(col=1, row=1, ncol=2, nrow=1), newpage=FALSE)
      if(!is.null(breaks)){
        ff <- plotMutinfo(breaks, D, interval)
      }else{
        ff <- xyplot(x~y,data=data.frame(x=1:10,y=1:10))
      }
      plot(ff, split=c(col=2, row=1, ncol=2, nrow=1), newpage=FALSE)
    })
  }
}

test_bvseqs01 <- function(){
  i <- 1
  pp <- clst:::pullTab(bvseqs$dmat, bvseqs$groupTab, index=i)
  cc <- with(pp, classifyIter(dmat, groupTab, dvect,
                              verbose=VERBOSE))

  printClst(cc, groupNames=bvseqs$taxNames)
  
  pdf(print(pdfName()))
  plot.cc(cc)
  dev.off()  
}

test_bvseqs02 <- function(){

  startingRank <- 'genus'
  thresh <- with(bvseqs, findThreshold(dmat, groupTab[[startingRank]],
                                       na.rm=TRUE,
                                       verbose=VERBOSE))

  i <- 1
  pp <- clst:::pullTab(bvseqs$dmat, bvseqs$groupTab, index=i)
  cc <- with(pp, classifyIter(dmat, groupTab, dvect,
                              dStart=thresh$D,
                              minScore=0.65,
                              verbose=VERBOSE))

  printClst(cc, groupNames=bvseqs$taxNames)
  
  pdf(print(pdfName()))
  plot.cc(cc)
  dev.off()  
}

test_bvseqs03 <- function(){

  startingRank <- 'family'
  thresh <- with(bvseqs, findThreshold(dmat, groupTab[[startingRank]],
                                       na.rm=TRUE,
                                       verbose=VERBOSE))

  i <- 1
  pp <- clst:::pullTab(bvseqs$dmat, bvseqs$groupTab, index=i)
  cc <- with(pp, classifyIter(dmat, groupTab, dvect,
                              dStart=thresh$D,
                              multiple=FALSE,
                              minScore=0.45,                              
                              verbose=VERBOSE))

  printClst(cc, groupNames=bvseqs$taxNames)
  
  pdf(print(pdfName()))
  plot.cc(cc)
  dev.off()  
}
