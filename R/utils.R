template <- function(fstr, argv, default='', verbose=FALSE){
  rstr <- '(?<=%)\\([^)]+\\)'

  ## support for vectorized operation
  if(length(fstr) > 1){
    return(sapply(fstr, function(x){
      template(x, argv, default)}, USE.NAMES=FALSE))
  }
  
  if(is.data.frame(argv)){
    return(apply(argv, MARGIN=1, function(x){
      template(fstr, as.list(x), default)}))
  }
  
  starts <- gregexpr(rstr, fstr, perl=TRUE)[[1]]    

  if (starts[1] != -1){        
    words <- substring(fstr, starts+1, starts+attr(starts,'match.length')-2)    

    ## restrict to keys in fstr; ensure proper order; also replace
    ## zero-length elements in repl (eg, NULL, character(0))    
    repl <- ifelse(sapply(argv[words],length) > 0, argv[words], default)

    if(verbose){
      cat('template:\n')
      cat(fstr)
      cat('\n')
      cat(gsub(rstr, '', fstr, perl=TRUE))
      cat('\n')
      str(repl)
    }
    
    do.call(
            gettextf,
            c(list(fmt=gsub(rstr, '', fstr, perl=TRUE)), repl)
            )
  }
  else{
    ## no replacement directives were found
    fstr
  }
}

safeStr <- function(val, floatfmt='%.2f'){
  if(is.null(val)){
    output <- '*NULL*'
  }else if(length(val) == 0){
    output <- 'len(0)'
  }else{
    output <- sapply(val,
                     function(vv){
                       if(is.na(vv)){
                         return('*NA*')
                       }else if(is.numeric(vv) && !is.integer(vv)){
                         return(gettextf(floatfmt, vv))
                       }else{
                         return(gettextf('%s', vv))
                       }
                     }
                     )
  }
  return(output)
}

abbrevTaxName <- function(taxName, genusLength=1){
  
  if(length(taxName) == 1){
    
    spl = unlist(strsplit(taxName," "))
    first = spl[1]
    last = spl[2]
    if (is.na(last) || first == "genus" || first == "family"){
      output = taxName
    }
    else{
      output = paste(substring(first,1,genusLength), last)
    }
  }else{
    output <- sapply(taxName, abbrevTaxName, genusLength)
  }
  
  return(output)
  
}

