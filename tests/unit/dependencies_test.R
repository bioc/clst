source('unit/setup.R')

test_ROC <- function(){

  library(ROC)
  pv <- package.version('ROC')
  if(VERBOSE) print(paste('ROC',pv))

}


