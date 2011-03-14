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

test_strep01 <- function(){
  i <- 1
  pp <- clst:::pull(strep$dmat1, strep$abbrev, index=i)
  cc <- with(pp, clst:::classifier(dmat, groups, dvect, verbose=VERBOSE))
  
  checkTrue(!is.null(cc$details))
  checkTrue(!is.null(cc$input$dvect))
  checkTrue(!is.null(cc$input$groups))
}

test_strep02 <- function(){
  i <- 1
  pp <- clst:::pull(strep$dmat1, strep$abbrev, index=i)
  cc <- with(pp, clst:::classifier(dmat, groups, dvect, keep.data=FALSE,
                                    verbose=VERBOSE))
  checkTrue(is.null(cc$details))
  checkTrue(is.null(cc$input$dvect))
  checkTrue(is.null(cc$input$groups))  
}

test_strep03 <- function(){
  i <- 50
  pp <- clst:::pull(strep$dmat1, strep$abbrev, index=i)
  cc <- with(pp, clst:::classifier(dmat, groups, dvect, verbose=VERBOSE))
  
  checkTrue(!is.null(cc$details))
  checkTrue(!is.null(cc$input$dvect))
  checkTrue(!is.null(cc$input$groups))
}
