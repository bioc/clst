classify <- function(dmat,
                     groups,
                     dvect,
                     method='mutinfo',
                     minScore=0.45,
                     doffset=0.5,
                     dStart=NA,
                     maxDepth=10,
                     minGroupSize=2,
                     objNames=names(dvect), ## can be a vector
                     keep.data=TRUE,
                     ...,
                     verbose=FALSE
                     ){
  
  ## check for incompatibilities among arguments
  stopifnot(ncol(dmat)==length(groups) && nrow(dmat) == length(groups))
  stopifnot(length(dvect)==length(groups))
  stopifnot(length(objNames)==length(groups))

  ## for use of match.call and eval, see R-lang.pdf sec 6.5
  origCall <- thisCall <- match.call()

  cc <- list()
  keep <- rep(TRUE,length(groups))
  for(depth in seq(maxDepth)){
    ## args unique to classifier
    thisCall[[1]] <- as.name("classifier")
    thisCall$depth <- depth

    ## shared args
    thisCall$dmat <- dmat[keep,keep]
    thisCall$groups <- factor(groups[keep])
    thisCall$dvect <- dvect[keep]
    thisCall$objNames <- objNames[keep]
    thisCall$dStart <- ifelse(depth==1,dStart,NA)

    ## args unique to classify
    thisCall$maxDepth <- NULL

    if(verbose){
      spaces <- ifelse(depth > 0, paste(rep('   ',depth),collapse=''), '')
      indent <- gettextf('%sclassify   (depth %s):', spaces, depth)
      cat(gettextf('%s\n',indent))
    }

    ci <- eval(thisCall, sys.frame(sys.parent()))
    cc[[depth]] <- ci

    keep <- groups %in% ci$matches
    
    ## should we call classify again?
    checks <- list(
                   threshIsValid = !is.na(ci$thresh$D),
                   moreThanOneMatch = length(ci$matches) > 1,
                   remainingGroupsLargeEnough = max(table(groups[keep])) > minGroupSize,
                   remainingNotAllMatches = ifelse(
                       depth == 1, TRUE,
                       !setequal(groups[keep], cc[[depth-1]]$matches)
                      )
                   )
    keepGoing <- all(unlist(checks))
    
    if(verbose){
      for(n in names(checks)){
        cat(gettextf('%s %s: %s \n',indent, n, checks[[n]]))
      }
      if(keepGoing) cat(gettextf('%s running classify from depth %s\n',indent, depth))
      else cat(gettextf('%s stopping at depth %s\n',indent, depth))
    }
    
    if(!keepGoing){
      break
    }

  }

  ## if last iteration has no matches, remove it from the output
  if(depth > 1 & length(ci$matches) == 0){
    if(verbose){
      cat(gettextf('%s no matches on final iteration, returning previous result\n',indent))
    }
    cc <- cc[1:depth-1]
  }

  return(cc)
}

classifyIter <- function(dmat, groupTab, dvect, dStart=NA,
                         multiple=FALSE,
                         keep.data=TRUE,
                         ..., verbose=FALSE){

  ## Input
  ## -----
  ##
  ## * dmat - (see classify)
  ## * groupTab - a data.frame representing a taxonomy, with columns
  ##   in increasing order of specificity from left to right (ie,
  ##   Kingdom --> Species). Column names are used to name taxonomic
  ##   ranks. Rows correspond to margins of dmat.
  ## * dvect - numeric vector of distance from query sequence to each
  ##   reference corresponding to margins of dmat.
  ## * dStart - start with this value of D in place of the initial
  ##   calculation.
  ## * multiple - if TRUE, stops at the rank that yields at
  ##   least one match; if FALSE, continues to perform classification
  ##   until exactly one match is identified.
  ## * ... - additional arguments passed to clst::classify
  ## * verbose - TRUE for verbose screen output.

  stopifnot(nrow(dmat) == nrow(groupTab))
  stopifnot(all(rownames(dmat) == rownames(groupTab)))
  stopifnot(length(dvect) == nrow(dmat))

  lowestRank <- findLowestRank(groupTab, req=any, verbose=verbose)
  groups <- factor(groupTab[,lowestRank])
  nGroups <- length(setdiff(unique(groups),NA))

  if(verbose){
    cat(gettextf('classifyIter: starting rank is %s\n',
                 colnames(groupTab)[lowestRank]))
  }

  ## accumulate results in cc
  cc <- list()
  ii <- !is.na(groups)

  depth <- 1
  continue <- TRUE
  for(depth in seq(ncol(groupTab))){
    if(verbose){
      spaces <- ifelse(depth > 0, paste(rep('   ',depth),collapse=''), '')
      indent <- gettextf('%sclassifyIter (depth %s):', spaces, depth)
      ugroups <- unique(groups[ii])
      cat(gettextf('%s %s reference objects in %s starting groups [%s]\n',
                   indent, sum(ii), length(ugroups), paste(ugroups, collapse=', ')))
    }

    ## reset dStart to original value if provided
    dStart <- ifelse(depth==1, dStart, NA)

    ## address the special case of a single starting group
    if(is.na(dStart) & nGroups == 1){
      dd <- dmat[ii,ii]
      ## TODO: which of these approaches for calculating D in the
      ## setting of a single group is best?

      ## TODO: another approach here would be to maintain a table of
      ## relatedness among groups, and recruit members of the next most
      ## closely related group to assist in the calculation of D
      
      ## dStart <- median(dd[upper.tri(dd)], na.rm=TRUE)
      ## dStart <- quantile(dd[upper.tri(dd)], 0.95)
      dStart <- max(dd[upper.tri(dd)])

      if(verbose){
        cat(gettextf('%s only one group met screening criteria - setting dStart to %.4f\n',
                     indent, dStart))
      }
    }
    
    cci <- do.call(classify,
                   c(list(
                          dmat=dmat[ii,ii],
                          dvect=dvect[ii],
                          groups=factor(groups[ii]),
                          dStart=dStart,
                          keep.data=keep.data,
                          verbose=verbose
                          ),
                     ...)
                   )

    ci <- cci[[length(cci)]]

    ## add some attributes unique to classifyIter (as opposed to clst:::classify)
    ci$rank <- lowestRank
    
    matchCount <- length(ci$matches)

    ##     cat(gettextf('>>>>>>>>>>>>>> depth %s\n',depth))
    ##     show.classify(cci)
    ##     str(ci)
    ##     cat(gettextf('matchCount: %s\n',matchCount))
    ##     cat(gettextf('<<<<<<<<<<<<< depth %s\n',depth))

    cc[[depth]] <- ci
    
    if(matchCount == 1 | (multiple & matchCount > 1)){      
      ## we're done
      if(verbose) cat(gettextf('%s matchCount = %s --> breaking\n',
                               indent, matchCount))
      break
    }else{
      ## no matches - define new groups as the next highest level
      if(nrow(unique(groupTab)) == 1){
        if(verbose) cat(gettextf('%s single group represented at all remaining ranks - stopping',indent))
        break
      }

      lowestRank <- findLowestRank(groupTab[,1:(lowestRank-1), drop=FALSE],
                                   req=any, verbose=FALSE)

      if(is.na(lowestRank)){
        if(verbose) cat(gettextf('%s have exceeded highest rank - stopping',indent))
        break
      }else{
        if(verbose) cat(gettextf('%s no matches - redefining groups at rank %s (%s)\n',
                                 indent, lowestRank, names(lowestRank)))
        groups <- groupTab[,lowestRank]
        groupTab <- groupTab[,1:(lowestRank-1), drop=FALSE]
        ii <- !is.na(groups)
        nGroups <- length(setdiff(unique(groups),NA))
      }
    }

  } ## end loop

  if(verbose){
    printClst(cc)
  }

  return(cc)

}

classifier <- function(dmat,
                       groups,
                       dvect,
                       method='mutinfo',
                       minScore=0.45,
                       doffset=0.5,
                       dStart=NA,
                       minGroupSize=2,
                       objNames=names(dvect),
                       keep.data=TRUE,
                       ...,
                       verbose=FALSE,
                       depth=1
                       ){

  ## initialize a list of lists to contain the output
  params <- as.list(formals(classifier))
  mc <- match.call(expand.dots=TRUE)
  mc[[1]] <- NULL
  params[names(mc)] <- as.list(mc)

  params$dmat <- NULL
  params$groups <- NULL
  params$dvect <- NULL
  params$verbose <- NULL
  params[['...']] <- NULL

  output <- list(
                 depth=depth,
                 tally=NULL,
                 details=NULL,
                 matches=character(),
                 thresh=NULL,
                 params=params,
                 input=list(
                     dvect=NULL,
                     groups=NULL
                     )
                 )

  ## performs non-iterative classification
  if(verbose){
    spaces <- ifelse(depth > 0, paste(rep('   ',depth),collapse=''), '')
    indent <- gettextf('%sclassifier (depth %s):', spaces, depth)
    cat(gettextf('%s %s elements in reference set\n',indent,length(dvect)))
  }

  ## check for adequacy of data; return empty result if
  ## - inadequate minGroupSize; or
  ## - only a single group represented and dStart undefined

  groupTab <- table(factor(groups))

  if(max(groupTab) < minGroupSize | (is.na(dStart) & length(groupTab)==1)){
    ## return empty result showing no matches
    if(verbose){
      cat(gettextf('%s largest group has insufficient number of members or only a single group is represented (see table of groups below)\n', indent))
      print(groupTab)
    }
    output$thresh <- getThresh(method=method)

    if(length(dvect>0)){
      output$tally <- getTally(dvect, groups, bkpt=NA, minScore, doffset)
    }else{
      output$tally <- getTally(empty=TRUE)
    }

  }else{
    ## proceed with classification

    ## drop unused levels
    groups <- factor(groups)

    ## calculate a breakpoint if none is provided
    if(!is.na(dStart)){
      if(verbose) cat(gettextf('%s setting D to dStart=%.3f\n',indent, dStart))
      thresh <- getThresh(method='user', D=dStart)
      if(keep.data){
        thresh$distances <- partition(dmat, groups)
      }
      bkpt <- dStart
    }else{
      if(verbose) cat(gettextf('%s running findThreshold, method="%s"\n', indent, method))
      thresh <- findThreshold(dmat=dmat, groups=groups,
                              method=method,
                              keep.dists=keep.data,
                              verbose=verbose,
                              depth=depth,
                              ...
                              )
      bkpt <- thresh$D
    }

    if(verbose){
      cat(gettextf('%s bkpt = %.3f\n', indent, bkpt))
    }

    ## $tally contains the actual classification result
    tally <- getTally(dvect, groups, bkpt, minScore, doffset)

    matches <- rownames(tally)[tally$match==1]

    if(verbose){
      cat(gettextf('%s matches include %s\n', indent,
                   paste(matches, collapse=', ')))
    }

    ## ###################################################
    ## add results to output
    output$thresh <- thresh
    output$tally <- tally
    output$matches <- matches

    ## $details shows distance score for comparison
    ## of query object to each member of the reference set
    if(length(matches) > 0 & keep.data){
      details <- split(x=data.frame(
                           index=1:length(dvect),
                           dist=dvect,
                           group=factor(groups),
                           row.names=objNames
                           ),
                       f=ifelse(dvect < bkpt, 'below', 'above'))

      for(i in 1:length(details)){
        details[[i]] <- details[[i]][order(details[[i]]$dist),]
      }
      output$details <- details
    }
  }

  if(keep.data){
    output$input <- list(dvect=dvect, groups=groups)
  }

  return(output)
}

pullTab <- function(dmat, groupTab, index){
  
  stopifnot(nrow(groupTab) == nrow(dmat))

  output = list()
  output$dvect <- dmat[index,-index]
  names(output$dvect) <- rownames(dmat)[-index]
  output$dmat <- dmat[-index,-index,drop=FALSE]
  output$groupTab <- groupTab[-index,,drop=FALSE]
  return(output)
}


pull <- function(dmat, groups, index){
  
  stopifnot(length(groups) == nrow(dmat))

  output = list()
  output$dvect <- dmat[index,-index]
  names(output$dvect) <- rownames(dmat)[-index]
  output$qgroup <- as.character(groups[index])
  output$dmat <- dmat[-index,-index,drop=FALSE]
  output$groups <- groups[-index]
  return(output)
}

getTally <- function(dvect, groups, bkpt, minScore, doffset, empty=FALSE){

  if(empty){
    return(
           data.frame(
                      below=NA,
                      above=NA,
                      score=NA,
                      match=NA,
                      min=NA,
                      median=NA,
                      max=NA
                      )
           )
  }

  tt <- table(groups,
              factor(
                     ifelse(dvect <= bkpt, 'below', 'above'),
                     levels=c('below','above'),
                     ordered=TRUE
                     )
              )

  ss <- split(x=dvect, f=groups)

  ## drop=FALSE preserves rownames when length(tt) == 1
  tally <- data.frame(below=unclass(tt)[,'below',drop=FALSE],
                      above=unclass(tt)[,'above',drop=FALSE])

  ## score each match
  tally$score <- apply(tally,
                       MARGIN=1,
                       function(rr){
                         rr['below']/(sum(rr, na.rm=TRUE)+doffset)
                       }
                       )

  tally$match <- ifelse(tally$score < minScore, 0L, 1L)

  tally$min=sapply(ss, function(x) ifelse(length(x)>0, min(x), NA))
  tally$median=sapply(ss, median)
  tally$max=sapply(ss, function(x) ifelse(length(x)>0, max(x), NA))

  colsToRound <- c('score','min','median','max')
  tally[,colsToRound] <- round(tally[,colsToRound],2)

  return(tally)
}

findLowestRank <- function(taxTable,
                           req=all,
                           verbose=FALSE){

  ## Find the lowest rank with either all groups (require='all') or
  ## any group (require='any') defined.
  ##
  ## * taxTable - data.frame containing taxonomic hierarcy
  ## * req - builtin any or all

  defined <- apply(taxTable,
                   MARGIN=2,
                   function(x){
                     req(!is.na(x))
                   })

  ranks <- colnames(taxTable)

  if(any(defined)){
    lowestDefinedRank <- max(which(defined))
    names(lowestDefinedRank) <- ranks[lowestDefinedRank]
  }else{
    lowestDefinedRank <- NA
  }

  if(verbose){
    ## print(unique(taxTable))
    cat(gettextf('findLowestRank: all defined ranks: %s\n',
                 paste(ranks[defined], collapse=', ')))
    cat(gettextf('findLowestRank: lowest defined rank: %s (%s)\n',
                 lowestDefinedRank,
                 ifelse(is.na(lowestDefinedRank), 'NA',
                        names(lowestDefinedRank))))
  }

  return(lowestDefinedRank)
}

printClst <- function(cc, rows=8, nameWidth=30, groupNames){

  ## * cc - a list; output of classify or classifyIter
  ## * rows - maximum number of rows to display (NA for all)
  ## * nameWidth - maximum width of group names
  ## * groupNames - a named vector containing replacement names for groups
  ##   keyed by categories in groups (classify) or groupTab (classifyIter).

  safeStr <- safeStr
    
  matches <- cc[[length(cc)]]$matches
  if(!missing(groupNames)){
    matches <- ifelse(is.na(groupNames[matches]), matches, groupNames[matches])
  }

  msg <- ifelse(length(matches) > 0,paste(matches, collapse=', '),'(no match)')
  cat(gettextf('
============================================================
 matches: %s
============================================================\n',msg))
  
  ## show results of classification
  for(i in 1:length(cc)){
    level <- cc[[i]]

    vals <- as.list(level$params)    
    vals$depth <- level$depth
    vals$D <- safeStr(level$thresh$D, '%.2f')
    vals$prob <- level$thresh$params$prob
    method <- vals$method <- level$thresh$method

    ## classify16s-specific
    vals$rank <- names(level$rank)

    tally <- cc[[i]]$tally

    if(!missing(groupNames)){
      newRowNames <- groupNames[rownames(tally)]
      if(any(is.na(newRowNames))){
        msg <- gettextf('show.classify: groupNames is missing one or more groups at rank=%s: \n[%s]\n',
                        vals$rank,
                        paste(rownames(tally)[is.na(newRowNames)], collapse=', '))       
        ## show the original group label
        newRowNames <- ifelse(is.na(newRowNames), rownames(tally), newRowNames)
        ## stop(msg)
        warning(msg)
      }
      rownames(tally) <- newRowNames
      ## vals$matches <- paste(groupNames[level$matches], collapse=', ')
      vals$matches <- paste(with(level, ifelse(is.na(groupNames[matches]), matches, groupNames[matches])), collapse=', ')

    }else{
      vals$matches <- paste(level$matches, collapse=', ')
    }

    if(method == 'mutinfo'){    
      vals$pmmi <- safeStr(level$thresh$pmmi, '%.3f')
      templateStr <- '
--------------- rank %(rank)s (depth %(depth)s) --------------------
thresh: prob = %(prob)s
classifer: minScore = %(minScore)s  doffset = %(doffset)s dStart = %(dStart).2f
method = %(method)s  D = %(D)s pmmi = %(pmmi)s 
matches: %(matches)s
'
    }

    cat(template(templateStr, vals))

    msg <- ''
    if(!is.na(rows)){
      nr <- seq(min(rows, nrow(tally)))
      if(length(nr) < nrow(tally)){
        msg <- gettextf('(... plus %s groups not shown)\n',
                        nrow(tally)-length(nr))
      }
    }else{
      nr <- seq(nrow(tally))
    }

    ## process rownames
    rr <- rownames(tally)
    nchars <- nchar(rr)
    tooLong <- nchars > nameWidth
    if(any(tooLong)){
      beginning <- 2
      rr[tooLong] <- gettextf('%s...%s',
                              substr(rr[tooLong], 1, beginning),
                              substr(rr[tooLong], nchars[tooLong]-nameWidth+beginning+3,
                                     nchars[tooLong])
                              )
    }

    ## ensure uniqueness
    rrList <- as.list(rep(0,length(unique(rr))))
    names(rrList) <- unique(rr)

    for(i in seq_along(rr)){
      name <- rr[i]
      if(rrList[[name]] > 0){
        rr[[i]] <- gettextf('%s%s',name,rrList[[name]])
      }

      rrList[[name]] <- rrList[[name]] + 1
    }

    rownames(tally) <- rr

    ##print(tally[order(-tally$match, -tally$score, tally$median),][rr,])
    print(tally[order(-tally$score, tally$median),][nr,])
    cat(msg)
  }

  cat('----------------------------------------------\n')
}
