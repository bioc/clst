source('unit/setup.R')

library(clst)

test_actino <- function(){
  data(actino)
  checkTrue(exists('actino'))
  checkTrue(all(c("dmat1","dmat2","dmat3","taxa","abbrev") %in% names(actino)))
}

test_strep <- function(){
  data(strep)
  checkTrue(exists('strep'))
  checkTrue(all(c("dmat1","dmat2","dmat3","taxa","abbrev") %in% names(strep)))
}

test_bvseqs <- function(){
  data(bvseqs)
  with(bvseqs, {    
    checkTrue(all(rownames(groupTab) == rownames(dmat)))
    ids <- setdiff(unique(do.call(c, as.list(groupTab))),NA)
    checkTrue(all(ids %in% names(taxNames)))
  })
}
