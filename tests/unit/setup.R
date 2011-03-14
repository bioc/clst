unit_output <- 'unit_output'

## create destination for output files
dir.create(unit_output, showWarnings=FALSE)

VERBOSE <- TRUE

unitFuncName <- function(){
  unlist(sapply(sys.frames(), function(f)f$funcName))[1]
}

pdfName <- function(){
  gettextf('%s.pdf',file.path(unit_output,unitFuncName()))
}
