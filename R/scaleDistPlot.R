mkrep <- function(x, N){
  rep(x,
      ceiling(N/length(x))
      )[1:N]
}

getGlyphs <- function(N, shuffleGlyphs=NA,
                      colors=c(
                          'grey50',
                          'red',
                          'green',
                          'blue',
                          'orange',
                          'purple',
                          'brown',
                          'olivedrab',
                          'black'
                          ),
                      shapes=c(cir=21, sq=22, diam=23, tri=24, utri=25)
                      ){

  ## * N - number of groups to be identified with a glyph of distinct
  ##   shape and color.
  ## * shuffleGlyphs - modify permutation of shapes and colors given an
  ##   integer to serve as a ransom seed.

  if(N < 1){
    stop('N must be > 0')
  }
  
  if(is.numeric(shuffleGlyphs)){
    ## shuffle order of glyphs and maybe sample fewer
    set.seed(shuffleGlyphs)
    colors <- sample(colors, size=length(colors) - sample(c(0,2),1))
    shapes <- sample(shapes, size=length(shapes))
  }

  glyphs <- data.frame(
                       col=mkrep(colors, N),
                       pch=mkrep(shapes, N)
                       )

  return(glyphs)
}

scaleDistPlot <- function(dmat, groups, fill, X,
                          O, indices='no',
                          include,
                          display,
                          labels,
                          shuffleGlyphs=NA,
                          key='top',
                          keyCols=4,
                          glyphs,
                          xflip=FALSE,
                          yflip=FALSE,
                          ...){
  ## perform multidimensional scaling of dmat using cmdscale(stats)
  ## * dmat - square matrix of distances
  ## * groups - factor defining group membership of objects compared in dmat
  ## * fill - vector (logical or indices) of points to fill
  ## * X - vector of points to mark with an X
  ## * O - vector of points to mark with a circle
  ## * indices - label points with indices (all points if 'yes', or a subset indicated
  ##   by a vector)
  ## * include - boolean or numeric vector of elements to include in call to cmdscale
  ## * display - boolean or numeric vector of elements to include in call to display
  ## * labels - list or data frame with parameters $i indicating indices
  ##   and $text containing labels.
  ## * shuffleGlyphs - NA or an integer (argument to set.seed)
  ## * key - 'right' (single column), 'top' (variable number of columns),
  ##   or NULL for no key
  ## * keyCols - number of columns in key
  ## * glyphs - a data.frame with columns named $col and $pch corresponding
  ##   to elements of unique(groups) - eg, output of getGlyphs
  ## * xflip,yflip - if true, flip orientation of x or y axis, respectively
  ## * ... - additional arguments are passed to xyplot

  if(missing(groups)){
    groups <- factor(rep('group 1',nrow(dmat)))
  }else if(is.null(levels(groups))){
    groups <- factor(groups)
  }

  nPoints <- length(groups)
  nGroups <- length(levels(groups))

  ## ensure include is a logical vector
  if(missing(include)){
    include <- rep(TRUE,nPoints)
  }else if(is.numeric(include)){
    include <- seq(nPoints) %in% include
  }

  plotdata <- data.frame(i=seq_along(groups),
                         groups=groups,
                         include=include,
                         x=NA,
                         y=NA
                         )

  plotdata[include,c('x','y')] <- cmdscale(dmat[include,include])[,1:2]

  ## ensure display is a logical vector
  if(!missing(display)){
    if(is.numeric(display)){
      display <- seq(nPoints) %in% display
    }
    plotdata[!display,c('x','y')] <- NA
  }

  groupTab <- table(groups[include])
  gOK <- levels(groups) %in% names(groupTab[groupTab > 0])

  addO <- !missing(O)
  addX <- !missing(X)
  addLabels <- !missing(labels)

  if(identical(indices,'yes')){
    indices <- seq_along(groups)
    addIndices <- TRUE
  }else if(identical(indices,'no')){
    addIndices <- FALSE
  }else{
    addIndices <- TRUE
  }

  ## cex
  ## cex <- rep(1, nPoints)

  ## col and pch
  if(missing(glyphs)){
    glyphs <- getGlyphs(nGroups, shuffleGlyphs)
  }
  my.col <- as.character(glyphs$col)
  my.pch <- glyphs$pch

  ## attach(expand.grid(my.col=colors, my.pch=shapes)) NOTE: tried
  ## expand.grid - doesn't result in a permutation as effective as
  ## mkrep

  col <- my.col[groups]
  pch <- my.pch[groups]

  bg <- rep('transparent',nPoints)

  if(!(missing(fill))){
    ## change background for filled points
    bg[fill] <- col[fill]
    ## black outline for filled points
    col[fill] <- 'black'
  }

  xlim <- range(plotdata$x, na.rm=TRUE)
  ylim <- range(plotdata$y, na.rm=TRUE)

  if(xflip) plotdata$x <- -plotdata$x
  if(yflip) plotdata$y <- -plotdata$y

  if(!is.null(key)){
    ## adjust cex factor for small values of nGroups
    if(nGroups < 4){
      keycex <- 1
    }else{
      keycex <- 1/log(nGroups)
    }

    key <- list(
                space=key,
                text=list(
                    sapply(
                           ## omit if all representatives of group is excluded
                           which(gOK),
                           function(i){
                             gettextf('%s (%s)', names(groupTab)[i], groupTab[i])
                           }
                           ),
                    cex=keycex
                    ),
                points=list(
                    pch=my.pch[gOK],
                    col=my.col[gOK],
                    cex=keycex
                    ),
                columns=ifelse(key=='top',
                    max(min(length(which(gOK)), keyCols),1), ## ensure > 0
                    1)
                )
  }

  ff <- xyplot(y~x,
               data=plotdata,
               panel=function(x, y, ...){

                 panel.xyplot(x, y,
                              col=col,
                              pch=pch,
                              fill=bg,
                              ...)

                 ## mark points with "X"
                 if(addX){
                   panel.xyplot(x[X], y[X],
                                ##col=col[X],
                                col='black',
                                pch=4,
                                cex=1.5)
                 }

                 ## mark points with circles
                 if(addO){
                   panel.xyplot(x[O],y[O],
                                col='black',
                                pch=21, # circle
                                cex=2)
                 }

                 ## show index
                 if(addIndices){
                   panel.text(x[indices], y[indices],
                              labels=plotdata$i[indices],
                              pos=2, ## to the left
                              cex=1)
                 }


                   if(addLabels){
                     panel.text(x[labels$i], y[labels$i],
                                labels=labels$text,
                                pos=4, ## to the right
                                cex=1)
                   }
               },
               aspect=1,
               key=key,
               ...
               )

}
