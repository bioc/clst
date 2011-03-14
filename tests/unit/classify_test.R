source('unit/setup.R')

VERBOSE = TRUE

.setUp <- function(){
  if(VERBOSE){
    cat('\n')
  }
}

library(clst)
data(actino)
data(strep)

plot.cc <- function(cc){
  for(ci in cc){
    with(ci$thresh, {
      plot.new()
      ff <- plotDistances(distances, D, interval)
      plot(ff, split=c(col=1, row=1, ncol=2, nrow=1), newpage=FALSE)
      ff <- plotMutinfo(breaks, D, interval) 
      plot(ff, split=c(col=2, row=1, ncol=2, nrow=1), newpage=FALSE)
    })
  }
}

test_strep01 <- function(){
  i <- 1
  pp <- clst:::pull(strep$dmat1, strep$abbrev, index=i)
  cc <- with(pp, classify(dmat, groups, dvect, verbose=VERBOSE))

  pdf(print(pdfName()))
  plot.cc(cc)
  dev.off()  
}

