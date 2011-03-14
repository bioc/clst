partition <- function(dmat, groups, include, verbose=FALSE){

  ## Transform a square (or lower triangular) distance matrix into a
  ## data.frame containing a column of distances ($vals) along with a
  ## factor defining each distance as a within- or between-group
  ## comparison.
    
  ## * dmat - matrix of distances
  ## * groups - character vector defining categories of objects
  ##   compared in dmat
  ## * include - vector (numeric or boolean) indicating which elements to
  ##   retain in the output; comparisons including an excluded element
  ##   will have a value of NA

  if(any(is.na(groups))){
    stop('groups may not contain NA values')
  }
  
  N <- length(groups)
  stopifnot(ncol(dmat) == N)
  
  if(missing(include)){
    include <- rep(TRUE, N)
  }else if(mode(include)!='logical'){
    include <- include %in% seq(N)
  }

  dmat[!include,!include] <- NA

  ## square matrix defining membership in either same or different groups
  same <- sapply(groups, function(member) member==groups)
  different <- !same

  ## each comparison should only be represented once
  same[upper.tri(same, diag=TRUE)] <- FALSE
  different[upper.tri(different, diag=TRUE)] <- FALSE
  
  ## row and column indices
  rowmat <- col(dmat)
  colmat <- row(dmat)

  distances <- data.frame(
                          vals=c(dmat[same],dmat[different]),
                          comparison=factor(
                              c(
                                rep('within',sum(same)),
                                rep('between',sum(different))
                                ),
                              levels=c('within','between'),
                              ordered=TRUE),
                          row=c(rowmat[same], rowmat[different]),
                          col=c(colmat[same], colmat[different])
                          )
      
  ## sanity check
  stopifnot(nrow(distances) == N*(N-1)/2)

  return(distances)
}
