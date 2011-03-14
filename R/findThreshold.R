getThresh <- function(D=NA, pmmi=NA, interval=c(NA,NA),
                      breaks=NULL, distances=NULL,
                      method=NA, params=list()){
  as.list(environment())
}

findThreshold <- function(dmat, groups,
                          distances,
                          method='mutinfo',
                          prob=0.5,
                          na.rm=FALSE,
                          keep.dists=TRUE,
                          roundCuts=2,
                          minCuts=20, maxCuts=300,
                          targetCuts=100,
                          verbose=FALSE,
                          depth=1, ## not to be set by user
                          ...
                          ){
  
  stopifnot(method %in% c('mutinfo'))

  if(verbose){
    spaces <- ifelse(depth > 0, paste(rep('   ',depth),collapse=''), '')
    indent <- gettextf('%sfindThreshold (depth %s):', spaces, depth)
    cat(gettextf('%s (in): method = %s\n', indent, method))
  }

  params <- list(prob=prob, roundCuts=roundCuts,
                 minCuts=minCuts, maxCuts=maxCuts,
                 targetCuts=targetCuts)
  
  if(na.rm){
    if(missing(dmat) | missing(groups)){
      stop('Error in findThreshold: `na.rm` can only be used if `dmat` and `groups` are provided.')
    }

    ok <- !is.na(groups)
    groups <- groups[ok]
    dmat <- dmat[ok,ok]
  }
  
  if(missing(distances)){
    if(!identical(class(dmat),'matrix')){
      dmat <- as.matrix(dmat)
    }
    distances <- partition(dmat, groups)
  }
  
  if(verbose){
    categoryTab <- table(distances$comparison)
    cat(gettextf('%s (in): within: %s   between: %s\n',
                 indent, categoryTab['within'], categoryTab['between']))
  }

  if(method=='mutinfo'){
    ## point of maximal mutual information
    thresh <- findPMMI(
                       comparison=distances$comparison,
                       vals=distances$vals,
                       prob=prob,
                       roundCuts=roundCuts,
                       minCuts=minCuts, maxCuts=maxCuts,
                       targetCuts=targetCuts,
                       verbose=verbose
                       )    
    if(verbose){
      cat(gettextf('%s (out): %s = %.2f %s = %.2f\n',
                   indent,
                   "D", thresh$D,
                   "pmmi", thresh$pmmi))
    }
  }

  if(keep.dists){
    thresh$distances <- distances
  }

  thresh$params <- params
  
  return(thresh)
}

findPMMI <- function(vals, comparison, prob=NA, verbose=FALSE, ...){
  
  ## * vals - vector of distances
  ## * comparison - factor with levels c('within','between')
  ## * prob - a numeric between 0.5 and 1 that will be used to define the
  ##   range over which MI is calculated as [quantile(withins,prob), quantile(betweens,1-prob)]
  ## * cutArgs = list() containing arguments for findCutpoints
  
  spl <- split(vals, comparison)
  midpoints <- sapply(spl, median)
  spl <- spl[order(midpoints)]

  ## limit range to values that will generate two categories when
  ## used to cut vals
  interval <- c(min(spl$between), max(spl$within))
  if(interval[1] > interval[2]){
    interval <- range(vals)
  } 
  
  if(!is.na(prob)){
    ## define interval according to quantile=prob of w's and b's
    stopifnot(prob >= 0.5 && prob <= 1)
    lowerBound <- quantile(spl$within, prob)
    interval[1] <- max(lowerBound, interval[1])

    upperBound <- quantile(spl$between, 1-prob)
    interval[2] <- min(upperBound, interval[2])
  }

  cuts <- findCutpoints(vals, interval, verbose=verbose, ...)

  ## define truth so big distances correspond to 1, as that is what
  ## rocdemo.sca wants
  roc <- simpleROC(truth=as.numeric(comparison)-1,
                   data=vals, cutpts = cuts)
  mi <- MIfromROC(roc, length(spl$between), length(spl$within))

  ## identify cutpoint at maximum mutual information;
  ## in case of ties, use min
  d <- cuts[which.max(mi)]
  pmmi <- ifelse(length(d) > 0, max(mi, na.rm=TRUE), NA)
  D <- ifelse(length(d) > 0 & pmmi > 0, d, NA)

  breaks <- data.frame(x=cuts, y=mi)

  getThresh(D=D, pmmi=pmmi, interval=interval, breaks=breaks, method='mutinfo')
}

findCutpoints <- function(vals, interval, roundCuts=2,
                          minCuts=-Inf, maxCuts=Inf,
                          targetCuts=100,
                          verbose=FALSE){

  rounded <- round(
                   subset(vals, vals > interval[1] & vals < interval[2]),
                   digits=roundCuts
                   )
  cuts <- c(interval[1], unique(rounded), interval[2])

  if(length(cuts) > minCuts && length(cuts) < maxCuts){
    sort(cuts) ## TODO: necessary?
  }else{
    seq(interval[1], interval[2], length.out=targetCuts)
  }
}

simpleROC <- function(truth, data, cutpts,
                      markerLabel='marker',caseLabel='case'){

  ## mostly copied from ROC::rocdemo.sca; fixes warning in
  ## "if(is.na(cutpts))" and removes option of providing rule

  if (!all(sort(unique(truth)) == c(0, 1)))
    stop("'truth' variable must take values 0 or 1")

  if(missing(cutpts)) {
    udata <- unique(sort(data))
    delta <- min(diff(udata))/2
    cutpts <- c(udata - delta, udata[length(udata)] + delta)
  }
  np <- length(cutpts)

  ## requires import(ROC) in NAMESPACE
  rocResult <- .C("ROC", as.integer(truth), as.double(data),
                  as.double(cutpts), as.integer(length(truth)),
                  as.integer(length(cutpts)), spec = double(np),
                  sens = double(np), PACKAGE = "ROC")

  new("rocc", spec=rocResult$spec, sens=rocResult$sens, rule=dxrule.sca,
      cuts=cutpts, markerLabel=markerLabel, caseLabel=caseLabel)
}

MIfromROC <- function(roc, ngp1, ngp2, offset=0.5){
  ## add offset to each cell to avoid error caused by log10(0)
  t1 = roc@sens*ngp1
  t3 = ngp1 - t1
  t4 = ngp2 * roc@spec
  t2 = ngp2 - t4

  freqs = (cbind(t1,t2,t3,t4) + offset)/(ngp1+ngp2)
  MI = apply(freqs, 1, function(x) -sum(x*log2(x)))

  fr2 = (cbind(t1+t2, t3+t4) + offset)/(ngp1 + ngp2)
  mx = apply(fr2, 1, function(x) -sum(x*log2(x)))

  orig = c(ngp1, ngp2)/(ngp1+ngp2)
  Hg = -sum(orig*log2(orig))
  Hg+mx-MI
}

plotDistances <- function(distances,
                           D=NA, interval=NA,
                           ylab='distances', ...){
  
  ff <- bwplot(vals~comparison, data=distances,
               panel = function(x, y, ...){
                 panel.bwplot(x, y, ...)
                 panel.abline(h=D, lwd=2, col='red')
                 panel.abline(h=interval, lwd=1, col='red', lty=2)
               },
               ylab=ylab,
               ...)
  return(ff)
}

plotMutinfo <- function(breaks,
                     D=NA, interval=NA,
                     xlab='distance',
                     ylab='mutual information', ...){
  
  ff <- xyplot(y ~ x, data=breaks,
               panel = function(x, y, ...){
                 panel.xyplot(x, y, col='black')
                 panel.abline(v=D, col='red', lwd=2)
                 panel.abline(v=interval, lwd=1, col='red', lty=2)
               },
               xlab=xlab,
               ylab=ylab,
               ...
               )

  return(ff)
}
