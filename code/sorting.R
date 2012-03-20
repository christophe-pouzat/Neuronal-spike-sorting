cstValueSgts <- function ## Finds segements of data where the derivative is null
### Finds segments on a numeric vector where the first derivative is null
(##title<< Finds null derivative segments
 x ##<< a numeric vector
 ) {

  dx <- as.integer(abs(diff(x,2))/2 <= .Machine$double.eps) ## use "order 2" derivative estimates
  ddx <- diff(dx) ## ddx is zero where the derivative did not change,
                  ## it is one if the derivative becomes null and -1
                  ## is the derivative stops being null
  s <- (1:length(ddx))[ddx==1] ## s contains the indexes of the first
                               ## points of segments where the
                               ## derivative is null
  sapply(s,
         function(b) {
           n <- ddx[(b+1):length(ddx)]
           c(b+1,min((1:length(n))[n==-1]))
         }
         )
### an intgeger matrix with as many columns as segments where the
### derivative is null, the first row contains the index of the
### beginning of each segment, the second row contains the sgments'
### lengths.
}
###########################################################################
###########################################################################
###########################################################################
peaks <- function ## Find local maxima
### Finds the local maxima in a vector, or time series, or in each
### column of a matrix.
(##title<< Local maxima
 x, ##<< vector.
 span = 3, ##<< a peak is defined as an element in a sequence which is
           ##greater than all other elements within a window of width
           ##span centered at that element. The default value is 3,
           ##meaning that a peak is bigger than both of its
           ##neighbors. Default: 3.
 strict = TRUE ##<< logical flag: if TRUE, an element must be strictly
               ##greater than all other values in its window to be
               ##considered a peak. Default: TRUE.
 ) {
  stopifnot(is.vector(x))
  z <- embed(rev(x), dim = span)
  z <- z[rev(seq(nrow(z))), ]
  s <- span%/%2
  v <- max.col(z,"first") == 1 + s
  z <- c(rep(FALSE, s), v)
  ans <- c(z[1:(length(z) - s)], 
           rep(FALSE, span - 1))
  res <- ans
  res[is.na(res)] <- FALSE
  res <- (1:length(res))[res]
  attr(res,"call") <- match.call()
  attr(res,"nIDx") <- length(x)
  class(res) <- "eventsPos"
  res
  ##note<< This is adapted from the function with the same name in package splus2R.

### An object of class "eventsPos" is returned. This is essentially a vector of indexes.
}
###########################################################################
###########################################################################
###########################################################################
as.eventsPos <- function ## Transforms vector into eventsPos object
### Transforms vector into eventsPos object
(##title<< as.eventsPos
 x, ##<< an integer vector with strictly increasing elements.
 start, ##<< an integer, the sampling point at which observation started.
 end ##<< an integer, the sampling point at which observation ended.
 ) {
  x <- as.integer(x)
  stopifnot(all(diff(x)>0))
  if (missing(start)) start <- floor(x)
  if (missing(end)) end <- ceiling(x)
  stopifnot(all(x>=start))
  stopifnot(all(x<=end))
  attr(x,"call") <- match.call()
  attr(x,"nIDx") <- end-start+1
  class(x) <- "eventsPos"
  x
### An object of class "eventsPos" is returned. This is essentially a vector of indexes.
}
###########################################################################
###########################################################################
###########################################################################
print.eventsPos <- function ## Print method for eventsPos objects.
### Print method for eventsPos objects.
(##title<< print.evtsPos
 x, ##<< an eventsPos object.
 ... ##<< not used but included for compatibility with generic method.
 ) {
  cat("\neventsPos object with indexes of ", length(x)," events. \n", sep = "")
  cat("  Mean inter event interval: ", round(mean(diff(x)),digits=2), " sampling points, corresponding SD: ", round(sd(diff(x)),digits=2), " sampling points \n", sep = "")
  cat("  Smallest and largest inter event intervals: ", paste(range(diff(x)),collapse=" and "), " sampling points. \n\n",sep= "")
  
}
###########################################################################
###########################################################################
###########################################################################
cutSglEvt <- function ## make a single "cut" on data
### Produces a "cut" on data at a specific location with a specific length
(##title<< Make cuts on data
 evtPos, ##<< a numeric or integer interpretable as an indice: the
         ##posistion at which cuts will be produced
 data, ##<< a numeric vector of matrix containing the data. If
       ##vector the argument is converted as a single column matrix
       ##internally. The matrix rows are indexed by sampling points
       ##and its columns by recording sites / channels.
 before=14, ##<< an integer: the number of sampling points within the
            ##cut before the reference time given by evtPos.
 after=30 ##<< an integer: the number of sampling points after the
          ##reference time.
 ) {
  evtPos <- as.integer(evtPos) ## make sure evtPos is an integer
  before <- as.integer(before) ## make sure before is an integer
  stopifnot(0 <= before) ## make sure before is positive or null
  after <- as.integer(after)
  stopifnot(0 <= after) ## make sure after is positive or null
  if (is.vector(data)) data <- matrix(data,nc=1)
  ns <- dim(data)[2]
  dl <- dim(data)[1]
  stopifnot(0<evtPos, evtPos<=dl) ## make sure evtPos is within range
  sl <- before+after+1 ## the length of the cut
  keep <- -before:after + evtPos
  within <- 1 <= keep & keep <= dl
  kw <- keep[within]
  res <- sapply(1:ns,
                function(idx) {
                  v <- numeric(sl)
                  v[within] <- data[kw,idx]
                  v
                }
                )
  as.vector(res)
### A numeric vector with the cut(s). When several recording sites
### are used the cuts of each individual sites are placed one after
### the other.
}
###########################################################################
###########################################################################
###########################################################################
mkEvents <- function ## Make events' matrix
### Make events matrix out of data and events positions.
(##title<< Make events' matrix
 positions, ##<< an integer vector with events' positions as indices /
            ##sampling points or an 'eventsPos' object.
 data, ##<< a numeric vector of matrix containing the data or a 'ts' or 'mts' object. If
       ##vector the argument is converted as a single column matrix
       ##internally. The matrix rows are indexed by sampling points
       ##and its columns by recording sites / channels.
 before=14, ##<< an integer: the number of sampling points within the
            ##cut before the reference time given by evtPos.
 after=30 ##<< an integer: the number of sampling points after the
          ##reference time.
 ) {
  positions <- unclass(positions)
  data <- unclass(data)
  res <- sapply(positions,
                cutSglEvt,
                data,
                before,
                after)
  the.call <- match.call()
  attr(res,"positions") <- positions
  attr(res,"delta") <- NULL
  attr(res,"data") <- the.call[["data"]]
  attr(res,"before") <- before
  attr(res,"after") <- after
  attr(res,"numberOfSites") <- ifelse(is.matrix(data),dim(data)[2],1)
  attr(res,"call") <- match.call()
  class(res) <- "events"
  res
### A matrix with before + after + 1 rows and as many columns as
### elements in positions. Each column is an "event", that is, a set
### of cuts on the data. Attribute "positions" contains the value of
### the argument with the same name and attribute "data" contains the
### name of the corresponding argument, attribute "before" contains
### the value of the argument with the same name, attribute "after"
### contains the value of the argument with the same name, attribute
### "numberOfSites" contains the number of recording sites. Attribute
### "delta" is used when events are realligned on their mean waveforms
### (during "jitter cancellation"). The returned matrix is given an
### "events" class.
}
###########################################################################
###########################################################################
###########################################################################
summary.events <- function ## Summary method for events objects.
### Summary method for events objects.
(##title<< summary.events
 object, ##<< an events object.
 ... ##<< not used but included for compatibility with generic method.
 ) {
  b <- attr(object,"before")
  a <- attr(object,"after")
  ns <- attr(object,"numberOfSites")
  cat("\nevents object deriving from data set: ",attr(object,"data"),".\n",sep="")
  cat(" Events defined as cuts of ", a+b+1, " sampling points on each of the ",ns, " recording sites.\n",sep="")
  cat(" The 'reference' time of each event is located at point ", b+1, " of the cut.\n",sep="")
  if (!is.null(attr(object,"delta"))) {
    cat(" Events were realigned on median event.\n",sep="")
  }
  cat(" There are ", length(attr(object,"positions")), " events in the object.\n\n",sep="")
### Nothing is returned. A short description of the object is printed.  
}
###########################################################################
###########################################################################
###########################################################################
mkNoise <- function ## Make a "noise" matrix
### Makes a matrix of "noise events" out of "genuine events" positions and data.
(##title<< Make "noise" matrix
 positions, ##<< an integer vector with events' positions as indexes /
            ##sampling points.
 data, ##<< a numeric vector of matrix containing the data. If
       ##vector the argument is converted as a single column matrix
       ##internally. The matrix rows are indexed by sampling points
       ##and its columns by recording sites / channels.
 before=14, ##<< an integer: the number of sampling points within the
            ##cut before the reference time given by evtPos.
 after=30, ##<< an integer: the number of sampling points after the
           ##reference time.
 safetyFactor=2, ##<< a numeric: a factor applied to both before and
                 ##after in order to ensure an extra distance between
                 ##detected events and noise cuts.
 size=2000 ##<< an integer: the maximum number of noise events to cut.
 ) {
  positions <- unclass(positions)
  data <- unclass(data)
  if (!is.matrix(data)) data <- matrix(data,nc=1)
  size <- as.integer(size)
  stopifnot(0 < size) ## make sure size is a positive integer
  sl <- before+after+1
  ns <- dim(data)[2]
  i1 <- diff(positions) ## inter events intervals
  nbI <- (i1-round(safetyFactor*sl))%/%sl ## number of noise sweeps
                                          ## one can cut from each
                                          ## interval
  nbPossible <- min(size,
                    sum((nbI)[nbI>0])
                    )
  ## allocate next the memory for the noise events
  noiseMatrix <- matrix(0,
                        nr=ns*sl,
                        nc=nbPossible
                        )
  
  iV <- (1:length(i1))[nbI>0] ## A vector containing the indexes of
                              ## the (inter event) intervals from
                              ## which at least one noise sweep can be
                              ## cut.
  iIdx <- 1 ## an index running over the inter event intervals from
            ## which noise events can be cut.
  nInI <- nbI[iV[iIdx]] ## the number of noise sweeps that can be cut
                        ## from the "non empty" inter event interval
                        ## iV[iIdx].
  nIdx <- 1 ## An index running over the noise sweeps.
  noisePositions <- integer(nbPossible)
  while (nIdx <= nbPossible) {
    uInI <- 1 ## An index running over the noise sweeps that will be
              ## cut from a given "non empty" inter event interval.
    iPos <- positions[iV[iIdx]] + round(safetyFactor*sl)
    noisePositions[nIdx] <- iPos
    while (uInI <= nInI & 
           nIdx <= nbPossible
           ) {
      ii <- (-before:after) + iPos
      ns <- as.vector(data[ii,])
      noiseMatrix[,nIdx] <- ns
      nIdx <- nIdx + 1
      iPos <- iPos + sl
      uInI <- uInI + 1
    } ## End of while loop on uInI
    iIdx <- iIdx + 1
    nInI <- nbI[iV[iIdx]]
  } ## End of while loop on nIdx

  the.call <- match.call()
  attr(noiseMatrix,"positions") <- noisePositions
  attr(noiseMatrix,"delta") <- NULL
  attr(noiseMatrix,"data") <- the.call[["data"]]
  attr(noiseMatrix,"before") <- before
  attr(noiseMatrix,"after") <- after
  attr(noiseMatrix,"numberOfSites") <- ifelse(is.matrix(data),dim(data)[2],1)
  attr(noiseMatrix,"call") <- match.call()
  class(noiseMatrix) <- "events"
  noiseMatrix
### A matrix with before + after + 1 rows and as many columns as
### "noise events". Each column is a "noise event", that is, a set of
### cuts on the data which does not overlap with an actual event's
### cut. Attribute "positions" contains the location of the "reference
### time" of each cut and attribute "call" contains the function call
### which generated the result.
}
###########################################################################
###########################################################################
###########################################################################
"[.events" <- function ## Subsetting method for "events"
### Subsetting method for columns of events' objects
(##title<< Subsetting method for "events"
 x, ## an "events" object
 i, ## should be missing only subsetting of columns allowed.
 j, ## index vector for the columns
 drop = FALSE
 ) {
  y <- NextMethod("[")
  if (!missing(i)) return(NULL)
  if (is.matrix(y) && dim(y)[2] > 1) {
    attr(y,"positions") <- attr(x,"positions")[j]
    attr(y,"delta") <- attr(x,"delta")
    attr(y,"data") <- attr(x,"data")
    attr(y,"before") <- attr(x,"before")
    attr(y,"after") <- attr(x,"after")
    attr(y,"numberOfSites") <- attr(x,"numberOfSites")
    attr(y,"call") <- match.call()
    class(y) <- "events"
  }
  y
    
### an 'events' object is returned if there is more than one event
### left. A numeric vector is returned otherwise.
}
###########################################################################
###########################################################################
###########################################################################
t.events <- function ## Transpose method for "events"
### Transpose method for "events"
(##title<< t.events
 x ##<< an "events" object
 ) {
  t(unclass(x))
### a matrix, the transpose of the events matrix is returned. The
### 'events' class is lost.
}
###########################################################################
###########################################################################
###########################################################################
mean.events <- function ## Mean method for "events"
### Mean method for "events"
(##title<< mean.events
 x, ##<< an "events" object
 ... ##<< not used but necessary for generic method compatibility 
 ) {
  apply(unclass(x),1,mean,...)
### a numeric vector: the mean at each sampling point across events.
}
###########################################################################
###########################################################################
###########################################################################
median.events <- function ## Median method for "events"
### Median method for "events"
(##title<< median.events
 x, ##<< an "events" object
 na.rm = FALSE ##<< a logical flag, should not available values be
               ##removed? Default set to 'FALSE'.
 ) {
  apply(unclass(x),1,median,na.rm)
### a numeric vector: the median at each sampling point across events.
}
###########################################################################
###########################################################################
###########################################################################
'-.events' <- function ## Minus method for 'events' objects
### Minus method for 'events' object used to substract the mean or the
### median event for instance.
(##title<< '-.events' 
 e1, ##<< an 'events' object
 e2 ##<< a numeric vector with as many elements as e1 as rows.
 ) {
  stopifnot(length(e2) == dim(e1)[1])
  res <- unclass(e1)-e2
  attr(res,"positions") <- attr(e1,"positions")
  attr(res,"delta") <- attr(e1,"delta")
  attr(res,"data") <- attr(e1,"data")
  attr(res,"before") <- attr(e1,"before")
  attr(res,"after") <- attr(e1,"after")
  attr(res,"numberOfSites") <- attr(e1,"numberOfSites")
  attr(res,"call") <- match.call()
  class(res) <- "events"
  res
### An 'events' object.  
}
###########################################################################
###########################################################################
###########################################################################
plot.events <- function ## Plot method for "events" objects.
### Plots events of a data matrix with median event and its MAD.
(##title<< Plot events of a data matrix
 x, ##<< an "events" object.
 y=NULL,
 evts.lwd = 0.1, ##<< the 'lwd' used to plot the events.
 medAndMad = TRUE, ##<< a logical flag. Should the median and mad be
                   ##computed and superposed on the plot? Default
                   ##'TRUE'.
 evts.col = "black", ##<< a character or an integer specifying the
                     ##color used for the events. Default 'black'.
 med.col = "red", ##<< a character or an integer specifying the color
                  ##used for the median. Default 'red'.
 mad.col = "blue", ##<< a character or an integer specifying the color
                   ##used for the mad. Default 'blue'.
 x.bar = NULL, ##<< a numeric scalar, the length of abscissa scale
               ##bar. No bar drawn if set to 'NULL' which the default.
 y.bar = NULL ##<< a numeric scalar, the length of ordinate scale
              ##bar. No bar drawn if set to 'NULL' which the default.
 ) {

  nsites <- attr(x,"numberOfSites")
  ne <- dim(x)[2]
  cl <- dim(x)[1]/nsites
  ylim <- range(x)
  matplot(x,type="n",xlab="",ylab="",axes=FALSE,ylim=ylim)
  if (nsites > 1) {
    ii <- 2*(1:(nsites %/% 2))
    rect((ii-1)*cl,ylim[1],ii*cl,ylim[2],col="grey80",border=NA)
  }
  matlines(x,col=evts.col,lty=1,lwd=evts.lwd)
  if (medAndMad) {
    med <- apply(x,1,median)
    mad <- apply(x,1,mad)
    lines(med,col=med.col)
    lines(mad,col=mad.col)
  }
  if (!is.null(x.bar)) segments(x0=0,y0=ylim[1]+0.1*diff(ylim),x1=x.bar)
  if (!is.null(y.bar)) segments(x0=0,y0=ylim[1]+0.1*diff(ylim),y1=ylim[1]+0.1*diff(ylim)+y.bar)
### Nothing is returned. A plot is creates as a side effect.
}
###########################################################################
###########################################################################
###########################################################################
lines.events <- function ## Lines method for "events" objects.
### Add lines of events to the current device.
(##title<< Add lines with events to an existing plot
 x, ##<< an "events" object.
 evts.lwd = 0.1, ##<< the 'lwd' used to plot the events.
 evts.col = "black", ##<< a character or an integer specifying the
                     ##color used for the events. Default 'black'.
 ...
 ) {
  matlines(x,col=evts.col,lty=1,lwd=evts.lwd,...)
### Nothing is returned. Lines are added to an existing plot as a side effect.
}
###########################################################################
###########################################################################
###########################################################################
print.events <- function ## Print method for events objects
### Print method for events objects
(##title<< print.events
 x, ##<< an 'events' object
 ... ##<< additional parameters passed to the plot.events methods called internally
 ) {
  plot.events(x,...)
### Nothing is returned. A plot is creates as a side effect.
}
###########################################################################
###########################################################################
###########################################################################
sinc <- function ## sinc function value
### Returns the sinc of its argument or one of ist first three derivatives. 
(##title<< Sinc function
 x, ##<< a numeric vector
 deriv=0 ##<< a integer equal to 0, 1, 2 or 3. The order of the
         ##derivative of the sinc function. Default set at 0 (no
         ##derivative).
 ) {
  ##details<< the sinc function (as well as its derivatives) is not
  ##defined at 0 but its limit exists and is returned.
  deriv <- as.integer(deriv)
  if (deriv < 0 || deriv > 3) 
    stop("'deriv' must be between 0 and 3")

  if (deriv == 0)
    ## Maxima code to get the limit at 0:
    ## limit(sin(%pi*x)/%pi/x,x,0);
    res <- ifelse(x==0,1,sin(pi*x)/pi/x)
  if (deriv == 1)
    ## Maxima code to get the first derivative and its limit at 0:
    ## diff(sin(%pi*x)/%pi/x,x,1);
    ## limit(diff(sin(%pi*x)/%pi/x,x,1),x,0);
    res <- ifelse(x==0,0,
                  cos(pi*x)/x-
                  sin(pi*x)/pi/x^2)
  if (deriv == 2)
    ## Maxima code to get the second derivative and its limit at 0:
    ## diff(sin(%pi*x)/%pi/x,x,2);
    ## limit(diff(sin(%pi*x)/%pi/x,x,2),x,0);
    res <- ifelse(x==0,-pi^2/3,
                  -pi*sin(pi*x)/x+
                  2*sin(pi*x)/(pi*x^3)-
                  2*cos(pi*x)/x^2)
  if (deriv == 3)
    ## Maxima code to get the third derivative and its limit at 0:
    ## diff(sin(%pi*x)/%pi/x,x,3);
    ## limit(diff(sin(%pi*x)/%pi/x,x,3),x,0);
    res <- ifelse(x==0,0,
                  3*pi*sin(pi*x)/x^2-
                  6*sin(pi*x)/(pi*x^4)-
                  pi^2*cos(pi*x)/x+
                  6*cos(pi*x)/x^3)

  res
### a numeric vector with the sinc values (or the derivatives).
}
###########################################################################
###########################################################################
###########################################################################
sincfun <- function ## make sinc interpolating function
### Makes a sinc interpolating function given (optional) abscissa and ordinate values.
(##title<< Create sinc interpolating function.
 x, ##<< an integer vector of sampl points at which observations given
    ##in y were made, or a numeric vector with ordinate values. In the
    ##first case, the difference between successive elements of x must
    ##be excatly 1.
 y ##<< a numeric vector of ordinate values if x is given.
 ) {
  if (missing(y)) {
    y <- x
    x <- seq(along=y)
  } else {
    force(y)
  }
  if (any(diff(x) != 1))
    stop("'x' values must be one unit appart")

  function(t,deriv=0) {
    sapply(t,
           function(u) 
           (-1)^deriv*
           sum(y*sinc(x-u,deriv))
           )
  }
### a function of two arguments: t (a numeric vector) and deriv (an
### integer equal to 0, 1, 2 or 3). This functions returns the sinc
### interpolation.
}
###########################################################################
###########################################################################
###########################################################################
shiftEvent <- function ## Shifts an event in time
### Shifts an event by an arbitrary amount using a sinc interpolation
(##title<< Shifts event in time with sinc interpolation
 evtPos, ##<< an integer: the position in sample points of the event of interest.
 delta=0, ##<< a numeric: the shift magnitude. If f(t) is the waveform
          ##of the "centered" event, the function returns f(t+delta).
 data, ##<< a numeric vector of matrix containing the data.
 before, ##<< a positive integer, the number of sampling points to
         ##use before 'evtPos' + 'delta' in the cut.
 after, ##<< a positive integer, the number of sampling points to
        ##use after 'evtPos' + 'delta' in the cut.
 method = c("sinc","spline"), ##<< a character string: the method used
                              ##for the interpolation, '"sinc"'
                              ##performs sinc interpolation (default)
                              ##and '"spline"' performs a cubic spline
                              ##interpolation (faster).
 before4sinc=100, ##<< a positive integer, the number of sampling points to
             ##use before evtPos in order to build the sinc function
             ##interpolation.
 after4sinc=100 ##<< a positive integer, the number of sampling points to
           ##use after evtPos in order to build the sinc function
           ##interpolation.
 ) {
  evtPos <- as.integer(evtPos) ## make sure evtPos is an integer
  before <- as.integer(before) ## make sure before is an integer
  stopifnot(0 <= before) ## make sure before is positive or null
  after <- as.integer(after)
  stopifnot(0 <= after) ## make sure after is positive or null
  if (is.vector(data)) data <- matrix(data,nc=1)
  ns <- dim(data)[2]
  dl <- dim(data)[1]
  ii <- -before:after
  before4sinc <- as.integer(before4sinc) ## make sure before4sinc is an integer
  stopifnot(0 <= before4sinc) ## make sure before4sinc is positive or null
  after4sinc <- as.integer(after4sinc)
  stopifnot(0 <= after4sinc) ## make sure after4sinc is positive or null
  sl <- (before4sinc+after4sinc+1)
  jj <- -before4sinc:after4sinc
  keep <- jj+evtPos
  gg <- 1 <= keep & keep <= dl
  if (method[1] == "sinc") {
    evtF <- lapply(1:ns,
                   function(i) {
                     v <- numeric(sl)
                     v[gg] <- data[keep[gg],i]
                     sincfun(jj,v)
                   }
                   )
  }
  if (method[1] == "spline") {
    evtF <- lapply(1:ns,
                   function(i) {
                     v <- numeric(sl)
                     v[gg] <- data[keep[gg],i]
                     splinefun(jj,v)
                   }
                   )
  }
  unlist(lapply(evtF, function(l) l(ii+delta)))
### A numeric vector with ('before'+'after'+1) x number of recording
### sites elements.
}
###########################################################################
###########################################################################
###########################################################################
shiftValue <- function ## Estimates shift with different methods
### Estimates the shift between an event and a reference waveform with
### three possible methods: first order Taylor expansion, second order
### Taylor expansion and exact "sinc" based interpolation.
(##title<< Shift estimation
 evt, ##<< a numeric vector: the event.
 waveform, ##<< a numeric vector of the same length as 'evt': the
           ##reference waveform.
 method=c("order 2",
   "order 1",
   "Newton-Raphson",
   "optimize"), ##<< a character string with the method used: '"order
                ##2"' for a second order Taylor expansion (default),
                ##'"order 1"' for a first order Taylor expansion,
                ##'"Newton-Raphson"' for a Newton-Raphson algorithm
                ##based optimization and '"optimize"' for an
                ##'optimize' based optimization. The last two methods
                ##require functional estimation (for interpolation)
                ##and that's done with a cubic spline instead of a
                ##sinc (its much faster and nearly as accurate).
 waveformD, ##<< a numeric vector of the same length as 'waveform'
            ##with the first derivative of 'waveform', used only if
            ##'method' is one of the Taylor expansion methods.
 waveformDD, ##<< a numeric vector of the same length as 'waveform'
             ##with the second derivative of 'waveform', used only if
             ##'method' is '"order 2"',
 interval=c(-5,5), ##<< a numeric vector with two elements: the
                   ##smallest and largest values allowed for the
                   ##shift. Default: -5,5. Used if method is not
                   ##'"order 1"'.
 ... ##<< additional arguments passed to 'optimize' if method is not
     ##'"order 1"'.
 ) {
  if (method[1] == "order 1") {
    if (missing(waveformD)) waveformD <- c(0,diff(waveform,2)/2,0)
    result <- sum((evt-waveform)*waveformD)/sum(waveformD^2)
  }
  if (method[1] == "order 2") {
    if (missing(waveformD)) waveformD <- c(0,diff(waveform,2)/2,0)
    if (missing(waveformDD)) waveformDD <- c(0,diff(waveformD,2)/2,0)
    cFct <- function(delta) sum((evt-waveform- delta*waveformD - delta^2*waveformDD/2)^2)
    result <- optimize(cFct,interval,...)$minimum
  }
  if (method[1] == "Newton-Raphson") {
    if (missing(waveformD)) waveformD <- c(0,diff(waveform,2)/2,0)
    if (missing(waveformDD)) waveformDD <- c(0,diff(waveformD,2)/2,0)
    ii <- seq(along=evt)
    waveformF <- splinefun(ii,waveform)
    waveformDF <- splinefun(ii,waveformD)
    waveformDDF <- splinefun(ii,waveformDD)
    cFct <- function(delta) sum((evt-waveformF(ii+delta))^2)
    cFctD <- function(delta) -2*sum(waveformDF(ii+delta)*(evt-waveformF(ii+delta)))
    cFctDD <- function(delta) -2*sum(waveformDDF(ii+delta)*(evt-waveformF(ii+delta))-waveformDF(ii+delta)^2)
    nrFct <- function(delta) delta - cFctD(delta)/cFctDD(delta)
    result <- sum((evt-waveform)*waveformD)/sum(waveformD^2)
    result <- nrFct(result)
    while(abs(cFctD(result) > 1e-6)) result <- nrFct(result) 
  }
  if (method[1] == "optimize") {
    ii <- seq(along=evt)
    waveformF <- splinefun(ii,waveform)
    cFct <- function(delta) sum((evt-waveformF(ii+delta))^2)
    result <- optimize(cFct,interval,...)$minimum
  }
  result
### A numeric, the estimated shift.  
}
###########################################################################
###########################################################################
###########################################################################
alignWithProcrustes <- function ## Recursive alignment on template
### Given a set of events' times recursively estimate the template and
### the shift / jitter of individual events. Individual events are
### assumed to be shifted versions of the template plus a zero mean
### noise with finite variance.
(##title<< Procrustes events alignment
 evtsPos, ##<< an integer vector: the positions in sample points of the events of interest.
 data, ##<< a numeric vector of matrix containing the data.
 before, ##<< a positive integer, the number of sampling points to
         ##use before the reference time in the cuts.
 after, ##<< a positive integer, the number of sampling points to
        ##use after the reference time in the cut.
 method.shift = c("sinc","spline"), ##<< a character string: the method used
                              ##for the interpolation, '"sinc"'
                              ##performs sinc interpolation (default)
                              ##and '"spline"' performs a cubic spline
                              ##interpolation (faster).
 before4sinc=100, ##<< a positive integer, the number of sampling points to
                  ##use before evtPos in order to build the sinc function
                  ##interpolation.
 after4sinc=100, ##<< a positive integer, the number of sampling points to
                 ##use after evtPos in order to build the sinc function
                 ##interpolation.
 method.delta=c("order 2",
   "order 1",
   "Newton-Raphson",
   "optimize"), ##<< a character string with the method used for shift estimation: '"order
                ##2"' for a second order Taylor expansion (default),
                ##'"order 1"' for a first order Taylor expansion,
                ##'"Newton-Raphson"' for a Newton-Raphson algorithm
                ##based optimization and '"optimize"' for an
                ##'optimize' based optimization. The last two methods
                ##require functional estimation (for interpolation)
                ##and that's done with a cubic spline instead of a
                ##sinc (its much faster and nearly as accurate).
 interval=c(-5,5), ##<< a numeric vector with two elements: the
                   ##smallest and largest values allowed for the
                   ##shift. Default: -5,5. Used if method is not
                   ##'"order 1"'.
 tol=1, ##<< a positive number. The maximal allowed pointwise
        ##difference between two successive estimates of the template
        ##waveform normalized by a robust estimate of the
        ##corresponding SE.
 maxIt=10, ##<< a positive integer, the maximal number of allowed iterations.
 plot=TRUE, ##<< a logical flag. Should a plot showing successive
            ##estimates of the ideal waveform be shown. Default set
            ##to 'TRUE'.
 ... ##<< additional arguments passed to 'optimize' if method is not
     ##'"order 1"'.
 ) {
  evtsPos <- unclass(evtsPos)
  data <- unclass(data)
  ne <- length(evtsPos)
  isne <- 1/sqrt(ne)
  ev0 <- mkEvents(positions=evtsPos,data,before,after)
  t0 <- apply(ev0,1,median)
  evD <- apply(ev0,2,function(x) c(0,diff(x,2)/2,0))
  evDD <- apply(evD,2,function(x) c(0,diff(x,2)/2,0))
  evD.med <- apply(evD,1,median)
  evDD.med <- apply(evDD,1,median)
  ev.delta <- apply(ev0,2,
                    function(x)
                    shiftValue(evt=x,
                               waveform=t0,
                               method=method.delta,
                               waveformD=evD.med,
                               waveformDD=evDD.med,
                               interval=interval,
                               ...))
  deltaC <- ev.delta
  ev1 <- sapply(1:ne,
                function(i) shiftEvent(evtsPos[i],
                                       -deltaC[i],
                                       data,
                                       before,
                                       after,
                                       method.shift
                                       )
                )
  t1 <- apply(ev1,1,median)
  tmad <- apply(ev1,1,mad,na.rm=TRUE)*isne
  tDiff <- max(abs(t0-t1)/tmad)
  it <- 1
  if (plot) {
    plot(t0,type="l",lwd=1,
         main=paste("Iteration:",it-1),
         ylab="Amplitude")
    lines(t1,col=2,lwd=1)
    lines((t1-t0)/tmad,col=4)
    abline(h=tol,lty=2)
    abline(h=-tol,lty=2)
  }
  t0 <- t1
  while (tDiff > tol && it < maxIt) {
    evD <- apply(ev1,2,function(x) c(0,diff(x,2)/2,0))
    evDD <- apply(evD,2,function(x) c(0,diff(x,2)/2,0))
    evD.med <- apply(evD,1,median)
    evDD.med <- apply(evDD,1,median)
    ev.delta <- apply(ev1,2,
                      function(x)
                      shiftValue(evt=x,
                                 waveform=t0,
                                 method=method.delta,
                                 waveformD=evD.med,
                                 waveformDD=evDD.med,
                                 interval=interval,
                                 ...))
    deltaC <- deltaC + ev.delta
    ev1 <- sapply(1:ne,
                  function(i) shiftEvent(
                                         evtsPos[i],
                                         -deltaC[i],
                                         data,
                                         before,
                                         after,
                                         method.shift
                                         )
                  )
    t1 <- apply(ev1,1,median)
    tmad <- apply(ev1,1,mad,na.rm=TRUE)*isne
    tDiff <- max(abs(t0-t1)/tmad)
    it <- it+1
    cat(paste("Template difference: ",
              round(tDiff,digits=3),
              ", tolerance: ",
              round(tol,digits=3),
              "\n_______________________\n",
              sep=""
              )
        )
    if (plot) {
      plot(t0,type="l",lwd=1,
           main=paste("Iteration:",it-1),
           ylab="Amplitude")
      lines(t1,col=2,lwd=1)
      lines((t1-t0)/tmad,col=4)
      abline(h=tol,lty=2)
      abline(h=-tol,lty=2)
    }
    t0 <- t1
  }
  attr(ev1,"positions") <- evtsPos
  attr(ev1,"delta") <- deltaC
  attr(ev1,"data") <- match.call()[["data"]]
  attr(ev1,"before") <- before
  attr(ev1,"after") <- after
  attr(ev1,"numberOfSites") <- ifelse(is.matrix(data),dim(data)[2],1)
  attr(ev1,"call") <- match.call()
  class(ev1) <- "events"
  ev1
### a matrix with as many columns as events in evtsPos. Each column
### contains the estimated jitter followed by the shifted even.
}
###########################################################################
###########################################################################
###########################################################################
mkSimpleShift <- function ##Builds a "simple" shift evaluating function
### Builds a "simple" shift evaluating function based on a first order
### Taylor approximation. If g(t) = f(t+delta) + eps is observed,
### where delta is the shift and eps a zero mean noise with finite SD,
### the shift evaluation is based on the following approximation: g(t)
### ~ f(t) + delta * f'(t) + eps.
(##title<< Builds a "simple" shift evaluating function
 waveform, ##<< a numeric vector: the sampled function f above.
 waveformD=c(0,diff(waveform)/2,0), ##<< a numeric
                                                   ##vector: the
                                                   ##sampled function
                                                   ##f' above.
 waveformDD=c(0,diff(waveformD)/2,0)
 ) {

  above <- waveform > 0
  waveform <- waveform[above]
  waveformD <- waveformD[above]
  waveformDD <- waveformDD[above]
  ##squaredNorm <- sum(centeredWaveformD^2)
  function(delta,evt) sum((evt[above]-waveform- delta*waveformD - delta^2*waveformDD/2)^2)
### A function of a single argument, 'evt', a numeric vector returning
### the estimated shift.
}
###########################################################################
###########################################################################
###########################################################################
unjitter <- function ## cancel jitter and (optionally) estimate a
                     ## global scaling factor for a given pair of
                     ## event and template
### Cancel jitter (i.e., performs shift registration) and optionally
### estimate a global scaling factor for a pair of event end template.
(##title<< Cancel jitter and estimate scaling factor
 evtPos, ##<< an integer: the position in sample points of the event of interest.
 tempV, ##<< a numeric vector: the amplitude of the ideal waveform
        ##assumed to have generated the event on each recording site.
 interval=c(-5,5), ##<< a numeric vector with two elements: the
                   ##smallest and largest values allowed for the
                   ##jitter. Default: -5,5.
 pScale=FALSE, ##<< a logical flag. Should a scale factor be
               ##estimated? Default, FALSE.
 data, ##<< a numeric vector of matrix containing the data.
 before=100, ##<< a positive integer, the number of sampling points to
             ##use before evtPos in order to build the sinc function
             ##interpolation.
 after=100, ##<< a positive integer, the number of sampling points to
             ##use after evtPos in order to build the sinc function
             ##interpolation.
 ... ##<< additional (optional) arguments passed to function optim
     ##called internally when pScale is TRUE.
 ) {
  evtPos <- as.integer(evtPos) ## make sure evtPos is an integer
  before <- as.integer(before) ## make sure before is an integer
  stopifnot(0 <= before) ## make sure before is positive or null
  after <- as.integer(after)
  stopifnot(0 <= after) ## make sure after is positive or null
  if (is.vector(data)) data <- matrix(data,nc=1)
  ns <- dim(data)[2]
  dl <- dim(data)[1]
  sl <- (before+after+1)
  jj <- -before:after
  keep <- jj+evtPos
  gg <- 1 <= keep & keep <= dl
  evtF <- lapply(1:ns,
                 function(i) {
                   v <- numeric(sl)
                   v[gg] <- data[keep[gg],i]
                   sincfun(jj,v)
                 }
                 )
  osl <- length(tempV)/ns
  templateM <- matrix(tempV,nc=ns)
  peakPos <- which.max(apply(templateM,1,sum))
  ii <- 1:osl-peakPos+1
  costF <- function(delta) {
    se <- sapply(1:ns,
                 function(i) 
                 evtF[[i]](ii+delta)
                 )
    sum((templateM-se)^2)
  }

  delta <- optimize(costF,interval)$minimum
  se <- sapply(1:ns,
               function(i) 
               evtF[[i]](ii+delta)
               )
  se <- as.vector(se)
  if (!pScale) {
    c(delta,1,se)
  } else {
    p <- unname(coef(lm(se ~ tempV-1)))
    span <- diff(range(interval))
    x.min <- interval[1]
    toR <- function(x) {
      (x-x.min)/span
    }
    toI <- function(y) {
      span*exp(y)/(1+exp(y))+x.min
    }
    costF2 <- function(v) {
      delta <- toI(v[1])
      p <- exp(v[2])
      se <- sapply(1:ns,
                   function(i) 
                   evtF[[i]](ii+delta)
                   )
      sum((p*templateM-se)^2)
    }
    par0 <- c(toR(delta),log(p))
    if (!is.finite(costF2(par0))) {
      fit <- optim(par=par0,fn=costF2,method="Nelder-Mead",...)
    } else {
      fit <- tryCatch(optim(par=par0,fn=costF2,method="BFGS",...),
                      finally=optim(par=par0,fn=costF2,method="Nelder-Mead",...)
                      )
    }
    if (fit$convergence==0) {
      delta <- toI(fit$par[1])
      p <- exp(fit$par[2])
      se <- sapply(1:ns,
                   function(i) 
                   evtF[[i]](ii+delta)
                   )
      se <- as.vector(se)
      c(delta,p,se)
    } else {
      rep(NA,osl*ns+2)
    }
  }
### A numeric vector with the estimated jitter, the scaling factor,
### the unjittered (shift registered) event.
}
###########################################################################
###########################################################################
###########################################################################
shiftRegister <- function
### Recursive shift registration / jitter cancellation.
(##title<< Recursive shift registration.
 evtsPos, ##<< an integer vector with the positions (in sampling
          ##points) of the jittered / shifted events.
 temp0, ##<< a numeric vector: the initial guess for the amplitude of
        ##the ideal waveform assumed to have generated the event on
        ##each recording site.
 interval=c(-5,5), ##<< a numeric vector with two elements: the
                   ##smallest and largest values allowed for the
                   ##jitter. Default: -5,5.
 pScale=FALSE, ##<< a logical flag. Should a scale factor be
               ##estimated? Default, FALSE.
 data, ##<< a numeric vector of matrix containing the data.
 before=100, ##<< a positive integer, the number of sampling points to
             ##use before evtPos in order to build the sinc function
             ##interpolation.
 after=100, ##<< a positive integer, the number of sampling points to
            ##use after evtPos in order to build the sinc function
            ##interpolation.
 cluster, ##<< a snow cluster object.
 tol=1, ##<< a positive number. The maximal allowed pointwise
        ##difference between two successive estimates of the template
        ##waveform normalized by a robust estimate of the
        ##corresponding SE.
 maxIt=10, ##<< a positive integer, the maximal number of allowed iterations.
 plot=FALSE, ##<< a logical flag. Should a plot showing successive
             ##estimates of the ideal waveform be shown. Default set
             ##to 'FALSE'.
 ... ##<< additional (optional) arguments passed to function optim
     ##called internally when pScale is TRUE.
 ) {
  
  tempDiff <- max(abs(temp0))
  temp1 <- temp0
  ne <- length(evtsPos)
  isne <- 1/sqrt(ne)
  it <- 1
  while (tempDiff > tol && it <= maxIt) {
    ujM <- clusterApply(cl,
                        evtsPos,
                        unjitter,
                        temp1,
                        interval,
                        pScale,
                        data,
                        before,
                        after,
                        ...)
    ujM <- sapply(ujM, function(l) l)
    temp0 <- temp1
    temp1 <- apply(ujM[-(1:2),],1,median,na.rm=TRUE)
    tmad <- apply(ujM[-(1:2),],1,mad,na.rm=TRUE)*isne
    tempDiff <- max(abs(temp1-temp0)/tmad)
    it <- it+1
    cat(paste("Template difference: ",
              round(tempDiff,digits=3),
              ", tolerance: ",
              round(tol,digits=3),
              "\n_______________________\n",
              sep=""
              )
        )
    if (plot) {
      plot(temp0,type="l",lwd=1,
           main=paste("Iteration:",it-1),
           ylab="Amplitude")
      lines(temp1,col=2,lwd=1)
      lines((temp1-temp0)/tmad,col=4)
      abline(h=tol,lty=2)
      abline(h=-tol,lty=2)
    }
  }
  ujM
### a matrix with as many columns as events in evtsPos. Each column is
### the result of the call to function unjitter on that event at the
### last iteration.
}
###########################################################################
###########################################################################
###########################################################################
dkw.test <- function ## Two samples Dvorestky-Kiefer-Wolfowitz (DKW) test
### Performs a two samples Dvorestky-Kiefer-Wolfowitz (DKW) test
(##title<< Two samples DKW test
 sample1, ##<< a numeric vector with the first sample.
 sample2, ##<< a numeric vector with the second sample.
 alpha=0.01, ##<< a positive number smaller than 1, the probability of type I error.
 fun=FALSE ##<< a logical flag, should a list of functions returning
           ##the confidence bands be returned? Default FALSE.
 ) {
  
  cf <- sqrt(log(2/sqrt(alpha))/2)
  n1 <- length(sample1)
  n2 <- length(sample2)
  sample1 <- sort(sample1)
  sample2 <- sort(sample2)
  m1 <- (1:n1)/n1
  m1 <- stepfun(sample1, c(0,m1))
  u1 <- (1:n1)/n1 + cf/sqrt(n1)
  u1 <- ifelse(u1 > 1,1,u1)
  u1 <- stepfun(sample1, c(u1[1],u1))
  l1 <- (1:n1)/n1 - cf/sqrt(n1)
  l1 <- ifelse(l1 < 0,0,l1)
  l1 <- stepfun(sample1,c(0,l1))
  m2 <- (1:n2)/n2
  m2 <- stepfun(sample2, c(0,m2))
  u2 <- (1:n2)/n2 + cf/sqrt(n2)
  u2 <- ifelse(u2 > 1,1,u2)
  u2 <- stepfun(sample2, c(u2[1],u2))
  l2 <- (1:n2)/n2 - cf/sqrt(n2)
  l2 <- ifelse(l2 < 0,0,l2)
  l2 <- stepfun(sample2,c(0,l2))

  if (!fun) {
    all(u1(sample1) > l2(sample1) &
        l1(sample1) < u2(sample1)
        )
  } else {
    list(u1=u1,
         m1=m1,
         l1=l1,
         u2=u2,
         m2=m2,
         l2=l2)
  }
### TRUE or FALSE if the test succeeds, respectively fails when
### argument fun is FALSE. When the latter is TRUE a list of 6
### functions is returned. Each of the six component is a step
### function. m1, resp. m2, are the ECDFs of sample 1, resp. sample
### 2. u1, reps. u2, are the upper bounds of sample 1, resp. sample
### 2. l1, resp. l2, are the lower bounds of sample 1, resp. sample
### 2. When the two bands are plotted together, the probability that
### they do not overlap at one point (at least) if the null hypothesis
### (samples 1 and 2 were generated by the same distribution) is
### correct is alpha.
}
###########################################################################
###########################################################################
###########################################################################
mkGenModel <- function ## Makes a generative model
### Makes a generative model from ideal waveforms estimates of "pure"
### events, jitter samples, sacling factors samples and spike trains.
(##title<< Make a generative model
 positions, ##<< an integer vector with events' positions as indexes /
            ##sampling points.
 data, ##<< a numeric vector of matrix containing the data. If
       ##vector the argument is converted as a single column matrix
       ##internally. The matrix rows are indexed by sampling points
       ##and its columns by recording sites / channels.
 shiftRegisterList, ##<< a list returned by function shiftRegister
                    ##(applied to the same data set!).
 goodEvents, ##<< a logical vector with the events in spikeTrain which
             ##where considered as "good" or "pure" (i.e., not
             ##superposed events).
 classification, ##<< an integer vector with the classification for
                 ##each of the "good" events in spikeTrain.
 before=20, ##<< an integer: the number of sampling points within the
            ##cut before the reference time given by evtPos.
 after=70 ##<< an integer: the number of sampling points after the
          ##reference time.
 ) {
  positions <- as.integer(positions) ## make sure positions is an
                                     ## integer
  before <- as.integer(before) ## make sure before is an integer
  stopifnot(0 <= before) ## make sure before is positive or null
  after <- as.integer(after)
  stopifnot(0 <= after) ## make sure after is positive or null
  if (is.vector(data)) data <- matrix(data,nc=1)
  ns <- dim(data)[2] ## number of recording sites
  dl <- dim(data)[1] ## data length
  xxL <- -before:after
  cl <- length(xxL) # the length of the "cuts".
  classification <- as.integer(classification)
  nn <- length(unique(classification)) ## the number of neurons in the model.
  stopifnot( nn == length(shiftRegisterList)) ## shiftRegisterList
                                              ## must have as many
                                              ## elements as there are
                                              ## neurons in the model!

  ## make new cuts on the data using only the "good" events and cancel
  ## the jitter of those guys.
  goodPos <- positions[goodEvents]
  iiV <- integer(length(classification))
  for (i in 1:nn) iiV[classification==i] <- 1:sum(classification==i)
  ii <- 1:length(xxL)
  evtsL <- sapply(1:sum(goodEvents),
                  function(i) {
                  ev <- matrix(0,nr=length(xxL),nc=ns)
                  tt <- xxL+goodPos[i]
                  keep <- 0 < tt & tt <= dl
                  ev[ii[keep],] <- data[tt[keep],]
                  fL <- lapply(1:4,
                               function(j) sincfun(xxL,ev[,j])
                               )
                  res <- sapply(fL,
                                function(f)
                                f(xxL+shiftRegisterList[[classification[i]]][1,iiV[i]]))
                  as.vector(res)
                }
                )
  
  ## get the matrix of ideal waveforms, that is a matrix whose columns
  ## contain the median event of each neuron.
  long.med <- sapply(1:nn, function(i) apply(evtsL[,classification==i],1,median))

  ## Create functional estimates of the ideal waveforms 
  wF <- lapply(1:nn,
               function(j) {
                 res <- lapply(1:ns,
                               function(i) sincfun(xxL,long.med[ii+(i-1)*cl,j])
                               )
                 names(res) <- paste("site_",1:4,sep="")
                 res
               }
               )

  ## Extract jitter samples
  jsL <- lapply(shiftRegisterList, function(l) l[1,])
  ## Extract scaling factor samples
  sfL <- lapply(shiftRegisterList, function(l) l[2,])
  ## Make isi list
  nbs <- length(positions)
  timeAndClass <- integer(nbs*3)
  dim(timeAndClass) <- c(nbs,3)
  timeAndClass[,1] <-  positions
  relF <- sapply(1:nn, function(n) sum(classification==n))/length(classification)
  minISI <- sapply(1:nn, function(n) min(diff(goodPos[classification==n])))
  timeAndClass[goodEvents,2] <- classification
  for (i in (1:nbs)[!goodEvents]) {
    timeAndClass[i,2:3] <- sample(1:nn,2,prob=relF)
  } 
  stL <- lapply(1:nn, function(j) positions[sapply(1:nbs, function(i) any(timeAndClass[i,-1] == j))])
  isiL <- lapply(1:nn, function(n) {res <- diff(stL[[n]]);res[res>=minISI[n]]})
  for (k in 1:10) {
    for (i in (1:nbs)[!goodEvents]) {
      timeAndClass[i,2:3] <- sample(1:nn,2,prob=relF)
    }
    stL <- lapply(1:nn, function(j) positions[sapply(1:nbs, function(i) any(timeAndClass[i,-1] == j))])
    for (n in 1:nn) {
      res <- diff(stL[[n]])
      isiL[[n]] <- c(isiL[[n]],res[res>=minISI[n]])
    }
  }
  ## Prepare result
  res <- lapply(1:nn,
                function(n) {
                  list(fwe=wF[[n]],
                       jS=jsL[[n]],
                       sfS=sfL[[n]],
                       isiS=isiL[[n]],
                       fm=1
                       )
                }
                )
  names(res) <- paste("neuron_",1:nn,sep="")
  class(res) <- "generativeModel"
  res
### A list with as many components as neurons in the model. Each
### component is itself a list with the following components: fwe, a
### list of functions, the functional waveform estimate of the neuron
### on each recording site (each component of the list corresponds to
### one recording site); jS, a numeric vector, the empirical sample of
### jitters estimated for the neuron; sfS, a numeric vector, the
### empirical sample of scaling factors estimated for the neuron;
### isiS, a numeric vector, the empirical inter spike intervals
### sequence USING ONLY THE "GOOD" EVENTS OF THE NEURON; fm, a numeric
### scalar, a factor used latter during the simulations to correct for
### the fact that the isi sample excluded the "bad" events, its value
### is just the length of the positions vector divided by the number
### of "good" events in the data set. The returned object is given a
### "generativeModel" class.
}
###########################################################################
###########################################################################
###########################################################################
simulate.generativeModel <- function ## Simulate data with a generative model
### Simulate a data matrix from a generative model
(##title<< Simulate data with a generativeModel object
 object, ##<< an object of generativeModel class.
 nsim = 1000, ##<< number of response vectors to simulate. Defaults to
              ##'1000'.
 seed = NULL, ##<< an object specifying if and how the random number
              ##generator should be initialized (‘seeded’), either
              ##‘NULL’ or an integer that will be used in a call to
              ##‘set.seed’ before simulating the response. Defaults to
              ##'NULL'
 noiseMatrix, ##<< a numeric matrix containing a noise sample.
 ... ##<< not used put included for compatibility with the generic method.
 ) {
  ## Stuff taken verbatim from simulate.lm
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                     # initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  ns <- length(object[[1]][["fwe"]]) ## the number of recording sites
  nn <- length(object) ## the number of neurons
  ## get mean firing rates of neurons
  nu <- sapply(object, function(l) l$fm/mean(l$isiS))
  nu.total <- sum(nu)
  sim.duration <- 5*nsim/nu.total
  ## simulate individual trains are renewal processes
  trainL <- lapply(1:nn,
                   function(n) {
                     nb <- round(sim.duration*nu[n])
                     tt <- cumsum(jitter(sample(object[[n]][["isiS"]]/object[[n]][["fm"]],nb,replace=TRUE),amount=0.5))
                     while (max(tt) < sim.duration) {
                       tn <- cumsum(jitter(sample(object[[n]][["isiS"]]/object[[n]][["fm"]],nb,replace=TRUE),amount=0.5))
                       tt <- c(tt, max(tt) + tn)
                     }
                     tt[tt <= sim.duration]
                   }
                   )
  trainS <- sort(unlist(trainL))
  last.time <- trainS[nsim]
  trainL <- lapply(trainL, function(l) l[l <= last.time])

  ## Make noise free data
  after <- evalq(length(x),env=environment(object[[1]][["fwe"]][[1]]))
  simD <- matrix(0,nr=last.time+after,nc=ns)
  ## Add events to trace
  ii <- evalq(x,env=environment(object[[1]][["fwe"]][[1]]))
  for (n in 1:nn){
    tt <- trainL[[n]]
    nbs <- length(tt)
    ##jj <- runif(nbs,min=-0.5,max=0.5) ## jitter are uniform on (-0.5,0.5)
    sf <- sample(object[[n]][["sfS"]],nbs,replace=TRUE)
    for (j in 1:nbs) {
      roundPeakPos <- round(tt[j])
      delta <- tt[j]-roundPeakPos
      ee <- sapply(object[[n]][["fwe"]], function(f) f(ii-delta))*sf[j]
      pp <- ii+roundPeakPos
      keep <- 0 < pp & pp <= (last.time+after)
      simD[pp[keep],] <- simD[pp[keep],] + ee[keep,]
    }
  }

  cl <- dim(noiseMatrix)[1]/ns
  nbn <- dim(simD)[1] %/% cl
  nSamp <- sample(1:dim(noiseMatrix)[2],nbn+1,replace=TRUE)
  bigNoiseMatrix <- noiseMatrix[,nSamp]
  noiseD <- sapply(1:ns, function(i) as.vector(bigNoiseMatrix[(1:cl)+(i-1)*cl,]))
  simD <- simD + noiseD[1:dim(simD)[1],]
  trainS <- unlist(trainL)
  classification <- unlist(lapply(1:nn, function(n) rep(n,length(trainL[[n]]))))
  attr(simD,"positions") <- trainS
  attr(simD,"classification") <- as.integer(classification)
  attr(simD,"call") <- match.call()
  simD
### A numeric matrix of simulated raw data. The returned
### matrix has three attributes: positions, a numeric vector
### containing the reference / spike time of each event;
### classification, an integer vector with the neuron to which the
### previous reference time belongs; call, the matched call.
}
###########################################################################
###########################################################################
###########################################################################
superResolution <- function ## Find overlapping events from simulation
### Resolves superpositions from simulated data
(##title<< Superposition resolution in simulated data
 simMatrix ##<< a numeric matrix obtained by a call to
           ##'simulate.generativeModel; that is, a matrix with 3
           ##attributes: 'positions', 'classification' and 'call'.
 ) {
  before <- attr(simMatrix,"call")[["before"]]
  after <- attr(simMatrix,"call")[["after"]]
  classification <- attr(simMatrix,"classification")
  positions <- attr(simMatrix,"positions")
  nn <- length(unique(classification))
  ne <- length(positions)
  sapply(1:ne,
         function(eIdx) {
           et <- positions[eIdx]
           res <- rep(NA,nn)
           res[classification[eIdx]] <- 0
           on <- (1:nn)[-classification[eIdx]]
           for (oIdx in on) {
             st <- positions[classification == oIdx] - et
             st <- st[-after < st & st < before]
             if (length(st) > 0) res[oIdx] <- round(st[which.min(abs(st))])
           }
           res
         }
         )
### An integer matrix with as many rows as neurons in the model and as
### many columns as events in simMatrix. The neuron whose reference
### time corresponds to the position of the event (p_e) gets a 0
### entry. Every other neuron with a position included in the interval
### (-after:before) within the position of the event gets an entry
### corresponding to its position shift with respect to p_e. Negative
### values mean that the neuron fired before p_e. Neurons without an
### event in the above range (-after:before + p_e) get an NA entry.
}
###########################################################################
###########################################################################
###########################################################################
labelEvents <- function ## Generate a label for each event based on
                        ## the output of superResolution
### Generates a label for each event from the output matrix of
### superResolution
(##title<< Make labels for events
 superM ##<< an integer matrix: the result of a call to function
        ##'superResolution'.
 ) {
  ne <- dim(superM)[2] ## number of events
  nn <- dim(superM)[1] ## number of neurons
  nIdx <- 1:nn

  evtLbl <- function(x) {
    gIdx <- nIdx[!is.na(x)]
    lbl <- sapply(gIdx,
                  function(g) {
                    if (x[g] == 0) return(paste("c",g,sep=""))
                    if (x[g] < 0) return(paste("b",g,sep=""))
                    if (x[g] > 0) return(paste("a",g,sep=""))
                  }
                  )
    paste(lbl,collapse="_")
  }
  sapply(1:ne, function(i) evtLbl(superM[,i]))
### A vector of characters where 'cX' means that neuron 'X' is
### centered in the event, 'bY' means that neuron 'Y' fired before the
### event's center and 'aZ' means that neuron 'Z' fired after.
}
###########################################################################
###########################################################################
###########################################################################
writeGGobiXML <- function ## Writes a data matrix to file suitable for GGobi
### Writes a data matrix to a file on disc in an XML format suitable for GGobi
(##title<< Write matrix to XML file for GGobi visualization
 x, ##<< a numeric matrix or a data frame containing the data. If
    ##matrix events form rows and variables columns.
 x.name = "data", ##<< a character vector with default set to
                  ##'data'. The name under which the data will appear
                  ##in GGobi once loaded.
 filename, ##<< a character string: the name of the file where data
           ##will be written.
 x.description = "" ##<< a character vector with a short description
                     ##of the data set. Default set to '""'.
 ) {

  data.name <- "data"
  default.color <- "0"
  default.glyph <- "plus 1"
  catvars1 <- NULL
  x.colors <- NULL
  x.glyphs <- NULL
  x.id <- NULL
 
  ## Check x class
  M <- x
  class(M) <- "matrix"
  M <- data.frame(M)
  
  ## Write the header information
  cat(sep="",
      "<?xml version=\"1.0\"?>\n<!DOCTYPE ggobidata SYSTEM \"ggobi.dtd\">\n",
      file = filename)
  cat(sep="",
      "<ggobidata count=\"",
      1,
      "\">\n",
      file = filename,
      append = TRUE)
  cat(sep="",
      "<data name=\"",
      x.name,
      "\">\n",
      file = filename,
      append = TRUE)
  cat(sep="",
      "<description>\n",
      file = filename,
      append = TRUE)
  cat(sep="",
      x.description,
      "\n",
      file = filename,
      append = TRUE)
  cat(sep="",
      "</description>\n",
      file = filename,
      append = TRUE)

  p1 <- ncol(M)
  n1 <- nrow(M)
  cat(p1,n1,"\n")

  var.name1<-colnames(M)
  if (is.null(var.name1))
    for (i in 1:p1)
      var.name1<-c(var.name1, paste("Var ",i))

  cat(sep="",
      "<variables count=\"",
      p1,
      "\">\n",
      file = filename,
      append = TRUE)
  
  for (i in 1:p1) {
    
    if (is.factor(M[,i])) {
      l1<-length(levels(M[,i]))
      cat(sep="",
          "  <categoricalvariable name=\"",
          var.name1[i],
          "\" >\n",
          file = filename,
          append = TRUE)
      cat("    <levels count=\"",
          l1,
          "\" >\n",
          sep="",
          file = filename,
          append = TRUE)
      
      for (j in 1:l1) { 
        cat("    <level value=\"",
            j,
            "\" >",
            levels(M[,i])[j],
            "</level>\n",
            sep="",
            file = filename,
            append = TRUE)
      } ## End of for loop on j
      
      cat("    </levels>\n",
          file = filename,
          append = TRUE)
      cat("  </categoricalvariable>\n",
          file = filename,
          append = TRUE)
    } else if (i %in% catvars1) {
      cat(sep="",
          "  <categoricalvariable name=\"",
          var.name1[i],
          "\" levels=\"auto\"/>\n",
          file = filename,
          append = TRUE)
    } else {
      cat(sep="",
          "  <realvariable name=\"",
          var.name1[i],
          "\"/>\n",
          file = filename,
          append = TRUE)
    }
  } ## End of conditional on is.factor(M[,i])
  
  cat(sep="",
      "</variables>\n",
      file = filename,
      append = TRUE)
  cat(sep="",
      "<records count=\"",
      n1,
      "\" glyph=\"",
      default.glyph,
      "\" color=\"",
      default.color,
      "\" >\n",
      file = filename,
      append = TRUE)
  
  row.name1 <- rownames(M)

  if (is.null(row.name1))
    row.name1 <- c(1:n1)
  
  if (is.null(x.id)) 
    x.id <- c(1:n1)
  
  if (length(x.colors) != n1) {
    if (!is.null(x.colors))
      cat("Length of data 1 colors vector is not the same as the number of rows.\n")
    x.colors <- rep(default.color,n1)
  } ## End of conditional on length(x.colors) != n1
  if (length(x.glyphs) != n1) {
    if (!is.null(x.glyphs))
      cat("Length of data 1 glyphs vector is not the same as the number of rows.\n")
    x.glyphs <- rep(default.glyph,n1)
  } ## End of conditional on length(x.glyphs) != n1

  for(i in 1:n1) {  

    cat(sep="",
        "<record id=\"",
        x.id[i],
        "\" label=\"",
        row.name1[i],
        "\" ",
        file = filename,
        append = TRUE)
    
    cat(sep="",
        "color=\"",
        x.colors[i],
        "\" ",
        file = filename,
        append = TRUE)
    
    cat(sep="",
        "glyph=\"",
        x.glyphs[i],
        "\"",
        file = filename,
        append = TRUE)
    
    cat(sep="",
        ">\n",
        file = filename,
        append = TRUE)
    
    for (j in 1:p1)
      cat(M[i,j],
          " ",
          file = filename,
          append = TRUE)
    
    cat(sep="",
        "\n</record>\n",
        file = filename,
        append = TRUE)
    
  } ## End of for loop on i
  
  cat(sep="",
      "</records>\n</data>\n",
      file = filename,
      append = TRUE) 

  ## wrap-up file
  cat(sep="",
      "</ggobidata>",
      file = filename,
      append = TRUE)
### Nothing is returned, the function is used for its side effect: a file is written to disk.  
}
###########################################################################
###########################################################################
###########################################################################
explore.prcomp <- function ## Explore PCA results
### Plots the mean +/- one of the principal components times a user
### supplied factor.
(x, ##<< an object of class '"prcomp"'.
 pc=1, ##<< an integer: the pc index to add to the mean.
 factor=2, ##<< a numeric, the scaling factor; that is, the plot shows mean +/- factor * pc.
 m.col="black", ##<< a character string or an integer, the color used for mean.
 u.col="red", ##<< a character string or an integer, the color used for mean + factor * pc.
 l.col="blue", ##<< a character string or an integer, the color used for mean - factor * pc.
 xlab="Index", ##<< a character string with the abscissa label.
 ylab="Amplitude", ##<< a character string with the ordinate label.
 main, ##<< a character string with the title. If 'missing' (default) one is automatically generated.
 ... ##<< additional arguments passed to 'plot'.
 ) {
  if (missing(main)) {
    w <- x$sdev[pc]^2/sum(x$sdev^2)
    main <- paste("PC ",pc," (",round(100*w,digits=1),"%)",sep="")
  }
  u <- x$center + factor * x$rotation[,pc]
  l <- x$center - factor * x$rotation[,pc]
  ylim=range(c(l,u))
  plot(x$center,type="l",xlab=xlab,ylab=ylab,col=m.col,main=main,ylim=ylim,...)
  lines(u,col=u.col,...)
  lines(l,col=l.col,...)
### Nothing is returned, the fucntion is used for its side effect: a plot is created.
}
###########################################################################
###########################################################################
###########################################################################

###########################################################################
### Extra stuff for the "mother waveform plus scaling factors approach ####
###########################################################################
getScaledFit <- function ## estimate mother waveform plus scaling factors
### Given an event recorded from several recording sites a "mother
### waveform" with associated scaling factors are estimated.
(##title<< Mother waveform plus scaling factors estimation
 evt, ##<< a numeric vector: an event cut on several neighboring
      ##recording sites.
 ns=4, ##<< a positive integer: the number of recording sites.
 before=14, ##<< an integer: the number of sampling points before the
            ##"reference time" in the cut.
 w=NULL ##<< a numeric positive vector of weights: the estimated
        ##standard deviation at on each recording site. If NULL (default)
        ##the standard deviation is supposed to be homogenous (at 1).
 ) {

  evt <- matrix(evt,nc=ns)
  sl <- dim(evt)[1] ## the length of the cut on each site
  bigC <- which.max(evt[before+1,]) ## find site with largest amplitude

  ## get initial guess for ideal waveform (i.e., the observed
  ## one on the site where the amplitude is largest
  X <- evt[,bigC]
  ## get initial guesses for the scaling factors with a linear
  ## regression
  sF <- sapply((1:ns)[-bigC],
               function(idx) {
                 Y <- evt[,idx]
                 coef(lm(Y~X-1))[1]
               }
               )

  p0 <- c(sF,X)
  if (is.null(w)) {
    w <- rep(1,ns)
  } else {
    w <- rep(w,length.out=ns)
  }
  ## define a cost function (to be optimised): the residual
  ## sum of squares
  sF <- numeric(ns)
  res <- matrix(0,nr=sl,nc=ns)
  mW <- numeric(sl)
  costF <- function(p) {
    sF[bigC] <- 1
    sF[-bigC] <- p[1:(ns-1)]
    mW <- p[-(1:(ns-1))]
    for (i in 1:ns) res[,i] <- (evt[,i]-sF[i]*mW)/w[i]
    sum(res^2)
  }
  
  ## Get optimal estimation using the
  ## Broyden-Fletcher-Goldfarb-Shanno method
  theFit <- optim(p0,costF,method="BFGS")
  ## Check convergence and if necessary do some extra iterations
  if (theFit$convergence > 0) {
    extraFitIdx <- 1
    while(theFit$convergence > 0 & extraFitIdx < 5) {
      theFit <- optim(theFit$par,costF,method="BFGS")
      extraFitIdx <- extraFitIdx+1
    } ## end of while loop
  } ## end of conditional on theFit$convergence > 0
  pBest <- theFit$par
  scaleF <- rep(1,ns)
  scaleF[-bigC] <- pBest[1:(ns-1)]
  mW <- pBest[-(1:(ns-1))]

  ## return a list
  list(prediction=mW*rep(scaleF,each=sl),
       cost=theFit$value,
       allInOne=c(theFit$value,scaleF),
       fit=theFit)
### a list with four components: prediction is the fitted event; cost
### is the value of the cost function, that is, the residual sum of
### squares; allInOne is a numeric vector whose first element is the
### RSS and whose following elements are the scaling factors; fit is
### the object returned by optim, the optimization routine.
}
###########################################################################
###########################################################################
###########################################################################

###########################################################################
########### Interactive methods and functions #############################
###########################################################################
explore <- function ## generic method definition
(##title<< Generic method definition
 x,
 ...
 ) {
  UseMethod("explore")
}

explore.ts <- function ## intercative time series visualization
### Intercative time series visualization
(##title<< Intercative time series visualization
 x, ##<< a '"ts"' (time series) object
 offsetFactor=0.5, ##<< a numeric controlling the the spacing between
                   ##the recording sites on the plot. Smaller values
                   ##lead to closer spacing.
 ... ##<< additional arguments passed to 'matplot' and 'plot'
     ##functions called inernally.
 ) {

  stopifnot(inherits(x,"ts"))
  yRange <- range(x,na.rm=TRUE)
  plotPara <- list(tlim = start(x)[1] + c(0,0.1),
                   ylim = yRange,
                   yMin = yRange[1],
                   yMax = yRange[2],
                   firstTime = start(x)[1],
                   lastTime = end(x)[1],
                   keepGoing = TRUE)
  
  nFct <- function() {
    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    rightTime <- rightTime + timeRange
    if (rightTime > plotPara$lastTime) {
      cat("Recording end reached.\n ")
      rightTime <- plotPara$lastTime
    }
    plotPara$tlim <- c(rightTime - timeRange, rightTime)
    plotPara
  }
  fFct <- function() {
    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    leftTime <- leftTime - timeRange
    if (leftTime < plotPara$firstTime) {
      cat("Recording end reached.\n ")
      leftTime <- plotPara$firstTime
    }
    plotPara$tlim <- c(leftTime, leftTime + timeRange)
    plotPara
  }
  qFct <- function() {
    plotPara$keepGoing <- FALSE
    plotPara
  }
  ## Function tFct definition
  ## Allows the user to change the recording duration displayed on the window
  ## The user is invited to enter a factor which will be used to multiply the
  ## present duration displayed.
  ## If the resulting duration is too long a warning is given and the whole
  ## recording is shown.
  ## If possible the center of the displayed window is conserved.
  tFct <- function() {

    presentWindowLength <- diff(range(plotPara$tlim))
    tMessage <- paste("Present duration displayed: ", presentWindowLength, " \n", sep = "")
    tMessage <- paste(tMessage,
                      "By what factor do you want to multiply it? \n", sep = "")
    
    theFactor <- as.numeric(readline(tMessage))
    
    if (theFactor <= 0) {
      cat("A negative or null factor does not make sense.\n")
      return(plotPara)
    } ## End of conditional on theFactor <= 0

    ## Check that the new display length is reasonable
    totalLength <- plotPara$lastTime - plotPara$firstTime
    if (theFactor * presentWindowLength >= totalLength) {
      cat("Cannot show more data than available but only the entire record.\n ")
      plotPara$tlim[1] <- plotPara$firstTime
      plotPara$tlim[2] <- plotPara$lastTime
      return(plotPara)
    }

    windowCenter <- plotPara$tlim[1] + presentWindowLength / 2
    newLeft <- windowCenter - theFactor * presentWindowLength / 2
    newRight <- windowCenter + theFactor * presentWindowLength / 2
    
    if (!(newLeft >= plotPara$firstTime & newRight <= plotPara$lastTime)) {
      if (newLeft <= plotPara$firstTime) {
        cat("Cannot show data before the recording started, the displayed center wont be conserved.\n ")
        plotPara$tlim[1] <- plotPara$firstTime
        plotPara$tlim[2] <- plotPara$tlim[1] + theFactor * presentWindowLength
      }
      if (newRight >= plotPara$lastTime) {
        cat("Cannot show data after the recording ended, the displayed center wont be conserved.\n ")
        plotPara$tlim[2] <- plotPara$lastTime
        plotPara$tlim[1] <- plotPara$tlim[2] - theFactor * presentWindowLength
      }
      return(plotPara)
    } ## End of conditional on !(newLeft >= plotPara$firstTime & newRight <= plotPara$lastTime)

    plotPara$tlim[1] <- newLeft
    plotPara$tlim[2] <- newRight
    return(plotPara)
    
  }
  ## End of function tFct definition

  ## Function rFct definition
  ## Allows the user to change the maximal value displayed on the abscissa
  ## The user is invited to enter a value.
  rFct <- function() {

    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    tMessage <- paste("Present latest time displayed: ",
                      rightTime,
                      "\n", sep = "")
    tMessage <- paste(tMessage,
                      "What new latest time do want (return leaves things unchanged)? \n", sep = "")

    theNewTime <- as.numeric(readline(tMessage))
    
    if (is.na(theNewTime)) { ## Nothing entered, leave things unchanged 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)
    
    if (theNewTime <= plotPara$firstTime) {
      ## This choice does not make sense
      cat("Cannot display data before recording started.\n")
      return(plotPara)
    }

    if (theNewTime > plotPara$lastTime) {
      cat("Recording end reached.\n ")
      rightTime <- plotPara$lastTime
    } else {
      if (theNewTime <= leftTime) {
        ## The new latest time entered is smaller that the earliest time displayed
        cat("The new latest time is smaller than the earliest, adjustement will be made.\n")
        leftTime <- theNewTime - timeRange
        if (leftTime < plotPara$firstTime) {
          cat("Adjustment requires a change in displayed duration.\n")
          leftTime <- plotPara$firstTime
        }
      } ## End of conditional on theNewTime <= leftTime 
      rightTime <- theNewTime
    } ## End of conditional on theNewTime > plotPara$lastTime
  
    plotPara$tlim <- c(leftTime, rightTime)
    plotPara
    
  }

  ## Function lFct definition
  ## Allows the user to change the minimal value displayed on the abscissa
  ## The user is invited to enter a value.
  lFct <- function() {

    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    tMessage <- paste("Present earliest time displayed: ",
                      leftTime,
                      "\n", sep = "")
    tMessage <- paste(tMessage,
                      "What new earliest time do want (return leaves things unchanged)? \n", sep = "")

    theNewTime <- as.numeric(readline(tMessage))

    if (is.na(theNewTime)) { ## Nothing entered, leave things unchanged 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)

    if (theNewTime >= plotPara$lastTime) {
      ## This choice does not make sense
      cat("Cannot display data after recording ended.\n")
      return(plotPara)
    }

    if (theNewTime < plotPara$firstTime) {
      cat("Recording start reached.\n ")
      leftTime <- plotPara$firstTime
    } else {
      if (theNewTime >= rightTime) {
        ## The new earliest time entered is larger that the latest time displayed
        cat("The new earliest time is larger than the latest, adjustement will be made.\n")
        rightTime <- theNewTime + timeRange
        if (rightTime > plotPara$lastTime) {
          cat("Adjustment requires a change in displayed duration.\n")
          rightTime <- plotPara$lastTime
        }
      } ## End of conditional on theNewTime <= leftTime 
      leftTime <- theNewTime
    } ## End of conditional on theNewTime > plotPara$lastTime
    
    plotPara$tlim <- c(leftTime, rightTime)
    plotPara
    
  }


  ## Function yMaxFct definition
  ## Allows the user to change the maximal value displayed on the ordinate
  ## The user is invited to enter a value.
  yMaxFct <- function() {

    presentWindowRange <- range(plotPara$ylim)
    tMessage <- paste("Present range displayed: [",
                      paste(presentWindowRange, collapse = ","),
                      "] \n", sep = "")
    tMessage <- paste(tMessage,
                      "What new maximal ordinate value do want (return goes back to maximum)? \n", sep = "")
    
    theFactor <- as.numeric(readline(tMessage))

    if (is.na(theFactor)) {
      plotPara$ylim <- c(presentWindowRange[1],plotPara$yMax) 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)
    
    if (theFactor <= plotPara$ylim[1]) {
      cat("The maximum should be larger than the minimum.\n")
      return(plotPara)
    } ## End of conditional on theFactor <= plotPara$ylim[1]

    plotPara$ylim <- c(presentWindowRange[1],theFactor) 
    return(plotPara)
    
  }
  ## End of function yMaxFct definition

  ## Function yMinFct definition
  ## Allows the user to change the minimal value displayed on the ordinate
  ## The user is invited to enter a value.
  yMinFct <- function() {

    presentWindowRange <- range(plotPara$ylim)
    tMessage <- paste("Present range displayed: [",
                      paste(presentWindowRange, collapse = ","),
                      "] \n", sep = "")
    tMessage <- paste(tMessage,
                      "What new minimal ordinate value do want (return goes back to minimum)? \n", sep = "")
    
    theFactor <- as.numeric(readline(tMessage))

    if (is.na(theFactor)) {
      plotPara$ylim <- c(plotPara$yMin, presentWindowRange[2]) 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)
    
    if (theFactor >= plotPara$ylim[2]) {
      cat("The minimum should be smaller than the maximum.\n")
      return(plotPara)
    } ## End of conditional on theFactor >= plotPara$ylim[2]

    plotPara$ylim <- c(theFactor, presentWindowRange[2]) 
    return(plotPara)
    
  }
  ## End of function yMinFct definition

  show <- function(x,
                   plotPara,
                   ...) {

    s <- plotPara$tlim[1]
    e <- plotPara$tlim[2]
    y.m <- plotPara$ylim[1]
    y.M <- plotPara$ylim[2]
    m <- unclass(window(x,start=s,end=e))
    if (class(m) == "matrix") {
      m <- apply(m,2,function(x) ifelse(x < y.m, y.m,x))
      m <- apply(m,2,function(x) ifelse(x > y.M, y.M,x))
      ns <- dim(m)[2]
      offset <- c(0,-(1:(ns-1))*(y.M-y.m))
      m <- t(t(m)+offset*offsetFactor)
      matplot(m,type="l",lty=1,axes=FALSE,xlab="",ylab="",...)
    } else {
      m[m<y.m] <- y.m
      m[m>y.M] <- y.M
      plot(m,type="l",lty=1,axes=FALSE,xlab="",ylab="",ylim=c(y.m,y.M),...)
    }
  }

  plot.new()
  par(mar=c(0.5,0.5,0.5,0.5))
  show(x,plotPara,...)
  
  myMessage <- "Make a choice:\n n or 'return' (next); f (former); l (lower abscissa limit); r (upper abscissa limit) \n t (time scale); Y (upper ordinate limit); y (lower ordinate limit); q (quit) \n "

  while(plotPara$keepGoing) {
    
    myChoice <- readline(myMessage)

    plotPara <- switch(myChoice,
                       n = nFct(),
                       f = fFct(),
                       l = lFct(),
                       r = rFct(),
                       t = tFct(),
                       Y = yMaxFct(),
                       y = yMinFct(),
                       q = qFct(),
                       nFct()
                       )

    show(x,plotPara,...)
    
  } ## End of while loop on keepGoing

  dev.off()
  invisible()
### Nothing returned. Used for its side effects: interactive data
### exploration through plots.
}
###########################################################################
###########################################################################
###########################################################################
explore.eventsPos <- function ## Explore methods for eventsPos objects
### Explore methods for eventsPos objects
(##title<< explore.eventsPos
 x, ##<< an eventsPos object.
 y, ##<< a time series object.
 offsetFactor=0.5, ##<< a numeric controlling the the spacing between
                   ##the recording sites on the plot. Smaller values
                   ##lead to closer spacing.
 events.pch=16, ##<< an integer of a character: the ploting character
                ##used to indicate events.
 events.col=2, ##<< an integer or a character string coding the color used to indicate the event.
 ... ##<< additional arguments passed to 'matplot' and 'plot'
     ##functions called inernally.
 ) {
  stopifnot(inherits(y,"ts"))
  yRange <- range(y,na.rm=TRUE)
  plotPara <- list(tlim = start(y)[1] + c(0,0.1),
                   ylim = yRange,
                   yMin = yRange[1],
                   yMax = yRange[2],
                   firstTime = start(y)[1],
                   lastTime = end(y)[1],
                   keepGoing = TRUE)
  
  nFct <- function() {
    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    rightTime <- rightTime + timeRange
    if (rightTime > plotPara$lastTime) {
      cat("Recording end reached.\n ")
      rightTime <- plotPara$lastTime
    }
    plotPara$tlim <- c(rightTime - timeRange, rightTime)
    plotPara
  }
  fFct <- function() {
    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    leftTime <- leftTime - timeRange
    if (leftTime < plotPara$firstTime) {
      cat("Recording end reached.\n ")
      leftTime <- plotPara$firstTime
    }
    plotPara$tlim <- c(leftTime, leftTime + timeRange)
    plotPara
  }
  qFct <- function() {
    plotPara$keepGoing <- FALSE
    plotPara
  }
  ## Function tFct definition
  ## Allows the user to change the recording duration displayed on the window
  ## The user is invited to enter a factor which will be used to multiply the
  ## present duration displayed.
  ## If the resulting duration is too long a warning is given and the whole
  ## recording is shown.
  ## If possible the center of the displayed window is conserved.
  tFct <- function() {

    presentWindowLength <- diff(range(plotPara$tlim))
    tMessage <- paste("Present duration displayed: ", presentWindowLength, " \n", sep = "")
    tMessage <- paste(tMessage,
                      "By what factor do you want to multiply it? \n", sep = "")
    
    theFactor <- as.numeric(readline(tMessage))
    
    if (theFactor <= 0) {
      cat("A negative or null factor does not make sense.\n")
      return(plotPara)
    } ## End of conditional on theFactor <= 0

    ## Check that the new display length is reasonable
    totalLength <- plotPara$lastTime - plotPara$firstTime
    if (theFactor * presentWindowLength >= totalLength) {
      cat("Cannot show more data than available but only the entire record.\n ")
      plotPara$tlim[1] <- plotPara$firstTime
      plotPara$tlim[2] <- plotPara$lastTime
      return(plotPara)
    }

    windowCenter <- plotPara$tlim[1] + presentWindowLength / 2
    newLeft <- windowCenter - theFactor * presentWindowLength / 2
    newRight <- windowCenter + theFactor * presentWindowLength / 2
    
    if (!(newLeft >= plotPara$firstTime & newRight <= plotPara$lastTime)) {
      if (newLeft <= plotPara$firstTime) {
        cat("Cannot show data before the recording started, the displayed center wont be conserved.\n ")
        plotPara$tlim[1] <- plotPara$firstTime
        plotPara$tlim[2] <- plotPara$tlim[1] + theFactor * presentWindowLength
      }
      if (newRight >= plotPara$lastTime) {
        cat("Cannot show data after the recording ended, the displayed center wont be conserved.\n ")
        plotPara$tlim[2] <- plotPara$lastTime
        plotPara$tlim[1] <- plotPara$tlim[2] - theFactor * presentWindowLength
      }
      return(plotPara)
    } ## End of conditional on !(newLeft >= plotPara$firstTime & newRight <= plotPara$lastTime)

    plotPara$tlim[1] <- newLeft
    plotPara$tlim[2] <- newRight
    return(plotPara)
    
  }
  ## End of function tFct definition

  ## Function rFct definition
  ## Allows the user to change the maximal value displayed on the abscissa
  ## The user is invited to enter a value.
  rFct <- function() {

    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    tMessage <- paste("Present latest time displayed: ",
                      rightTime,
                      "\n", sep = "")
    tMessage <- paste(tMessage,
                      "What new latest time do want (return leaves things unchanged)? \n", sep = "")

    theNewTime <- as.numeric(readline(tMessage))
    
    if (is.na(theNewTime)) { ## Nothing entered, leave things unchanged 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)
    
    if (theNewTime <= plotPara$firstTime) {
      ## This choice does not make sense
      cat("Cannot display data before recording started.\n")
      return(plotPara)
    }

    if (theNewTime > plotPara$lastTime) {
      cat("Recording end reached.\n ")
      rightTime <- plotPara$lastTime
    } else {
      if (theNewTime <= leftTime) {
        ## The new latest time entered is smaller that the earliest time displayed
        cat("The new latest time is smaller than the earliest, adjustement will be made.\n")
        leftTime <- theNewTime - timeRange
        if (leftTime < plotPara$firstTime) {
          cat("Adjustment requires a change in displayed duration.\n")
          leftTime <- plotPara$firstTime
        }
      } ## End of conditional on theNewTime <= leftTime 
      rightTime <- theNewTime
    } ## End of conditional on theNewTime > plotPara$lastTime
  
    plotPara$tlim <- c(leftTime, rightTime)
    plotPara
    
  }

  ## Function lFct definition
  ## Allows the user to change the minimal value displayed on the abscissa
  ## The user is invited to enter a value.
  lFct <- function() {

    leftTime <- plotPara$tlim[1]
    rightTime <- plotPara$tlim[2]
    timeRange <- rightTime - leftTime
    tMessage <- paste("Present earliest time displayed: ",
                      leftTime,
                      "\n", sep = "")
    tMessage <- paste(tMessage,
                      "What new earliest time do want (return leaves things unchanged)? \n", sep = "")

    theNewTime <- as.numeric(readline(tMessage))

    if (is.na(theNewTime)) { ## Nothing entered, leave things unchanged 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)

    if (theNewTime >= plotPara$lastTime) {
      ## This choice does not make sense
      cat("Cannot display data after recording ended.\n")
      return(plotPara)
    }

    if (theNewTime < plotPara$firstTime) {
      cat("Recording start reached.\n ")
      leftTime <- plotPara$firstTime
    } else {
      if (theNewTime >= rightTime) {
        ## The new earliest time entered is larger that the latest time displayed
        cat("The new earliest time is larger than the latest, adjustement will be made.\n")
        rightTime <- theNewTime + timeRange
        if (rightTime > plotPara$lastTime) {
          cat("Adjustment requires a change in displayed duration.\n")
          rightTime <- plotPara$lastTime
        }
      } ## End of conditional on theNewTime <= leftTime 
      leftTime <- theNewTime
    } ## End of conditional on theNewTime > plotPara$lastTime
    
    plotPara$tlim <- c(leftTime, rightTime)
    plotPara
    
  }


  ## Function yMaxFct definition
  ## Allows the user to change the maximal value displayed on the ordinate
  ## The user is invited to enter a value.
  yMaxFct <- function() {

    presentWindowRange <- range(plotPara$ylim)
    tMessage <- paste("Present range displayed: [",
                      paste(presentWindowRange, collapse = ","),
                      "] \n", sep = "")
    tMessage <- paste(tMessage,
                      "What new maximal ordinate value do want (return goes back to maximum)? \n", sep = "")
    
    theFactor <- as.numeric(readline(tMessage))

    if (is.na(theFactor)) {
      plotPara$ylim <- c(presentWindowRange[1],plotPara$yMax) 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)
    
    if (theFactor <= plotPara$ylim[1]) {
      cat("The maximum should be larger than the minimum.\n")
      return(plotPara)
    } ## End of conditional on theFactor <= plotPara$ylim[1]

    plotPara$ylim <- c(presentWindowRange[1],theFactor) 
    return(plotPara)
    
  }
  ## End of function yMaxFct definition

  ## Function yMinFct definition
  ## Allows the user to change the minimal value displayed on the ordinate
  ## The user is invited to enter a value.
  yMinFct <- function() {

    presentWindowRange <- range(plotPara$ylim)
    tMessage <- paste("Present range displayed: [",
                      paste(presentWindowRange, collapse = ","),
                      "] \n", sep = "")
    tMessage <- paste(tMessage,
                      "What new minimal ordinate value do want (return goes back to minimum)? \n", sep = "")
    
    theFactor <- as.numeric(readline(tMessage))

    if (is.na(theFactor)) {
      plotPara$ylim <- c(plotPara$yMin, presentWindowRange[2]) 
      return(plotPara)
    } ## End of conditional on is.na(theFactor)
    
    if (theFactor >= plotPara$ylim[2]) {
      cat("The minimum should be smaller than the maximum.\n")
      return(plotPara)
    } ## End of conditional on theFactor >= plotPara$ylim[2]

    plotPara$ylim <- c(theFactor, presentWindowRange[2]) 
    return(plotPara)
    
  }
  ## End of function yMinFct definition

  show <- function(x,
                   y,
                   plotPara,
                   ...) {

    s <- plotPara$tlim[1]
    e <- plotPara$tlim[2]
    y.m <- plotPara$ylim[1]
    y.M <- plotPara$ylim[2]
    firstIdx <- round(max(1,s*frequency(y)))
    lastIdx <- round(min(end(y)[1]*frequency(y),e*frequency(y)))
    ii <- firstIdx:lastIdx
    xx <- x[firstIdx <= x & x <= lastIdx]
    if (class(y)[1] == "mts") {
      m <- y[ii,]
      if (length(xx) > 0)
        mAtx <- as.matrix(y)[xx,,drop=FALSE]
    } else {
      m <- y[ii]
      if (length(xx) > 0)
        mAtx <- y[xx]
    }
    if (class(m) == "matrix") {
      m <- apply(m,2,function(x) ifelse(x < y.m, y.m,x))
      m <- apply(m,2,function(x) ifelse(x > y.M, y.M,x))
      ns <- dim(m)[2]
      offset <- c(0,-(1:(ns-1))*(y.M-y.m))
      m <- t(t(m)+offset*offsetFactor)
      matplot(m,type="l",lty=1,axes=FALSE,xlab="",ylab="",...)
      if (length(xx) > 0) {
        mAtx <- t(t(mAtx)+offset*offsetFactor)
        matpoints(xx-ii[1]+1,mAtx,pch=events.pch,col=events.col)
      }
    } else {
      
      m[m<y.m] <- y.m
      m[m>y.M] <- y.M
      plot(m,type="l",lty=1,axes=FALSE,xlab="",ylab="",ylim=c(y.m,y.M),...)
      if (length(xx) > 0)
        points(xx-ii[1]+1,mAtx,pch=events.pch,col=events.col)
    }
  }

  plot.new()
  par(mar=c(0.5,0.5,0.5,0.5))
  show(x,y,plotPara,...)
  
  myMessage <- "Make a choice:\n n or 'return' (next); f (former); l (lower abscissa limit); r (upper abscissa limit) \n t (time scale); Y (upper ordinate limit); y (lower ordinate limit); q (quit) \n "

  while(plotPara$keepGoing) {
    
    myChoice <- readline(myMessage)

    plotPara <- switch(myChoice,
                       n = nFct(),
                       f = fFct(),
                       l = lFct(),
                       r = rFct(),
                       t = tFct(),
                       Y = yMaxFct(),
                       y = yMinFct(),
                       q = qFct(),
                       nFct()
                       )
    show(x,y,plotPara,...)
    
  } ## End of while loop on keepGoing

  dev.off()
  invisible()
### Nothing returned. Used for its side effects: interactive data
### exploration through plots.
}
###########################################################################
###########################################################################
###########################################################################
matchTemplate <- function #Compare cuts to set of templates
### Given an object of class events and a set of template functions, a
### systematic comparison of each event with each template is
### performed with jitter correction.
(evts, ##<< an object of class events
 templateList, ##<< a list of lists template functions. One list of
               ##functions per template and one function per recording
               ##site.
 interval=c(-5,5) ##<< a numeric vector of length 2, the range over
                  ##which the jitter is minimized.
 ) {
  stopifnot("events" %in% class(evts))
  before <- attr(evts,"before")
  after <- attr(evts,"after")
  bb <- (round(interval[1])-0.5):(round(interval[2])+0.5)
  cc <- round(interval[1]):round(interval[2])
  positions <- attr(evts,"positions")
  n <- dim(evts)[2]
  temp <- sapply(templateList,
                 function(l) as.vector(sapply(l, function(f) f(-before:after))))
  res <- matrix(integer(3*n),nr=3,nc=n)
  rownames(res) <- c("clusterID","position","δx1000")
  for (eIdx in 1:n) {
    evt <- evts[,eIdx]
    best <- which.min(apply((temp-evt)^2,2,sum))
    δ <- shiftValue(evt,temp[,best],"optimize",interval=interval)
    δprime <- δ - cc[findInterval(δ,bb)]
    pos <- positions[eIdx] - cc[findInterval(δ,bb)]
    res[,eIdx] <- c(best,pos,round(δprime*1000))
  }
  attr(res,"call") <- match.call()
  attr(res,"data") <- attr(evts,"data")
  attr(res,"templateList") <- templateList
  class(res) <- c("eventsMatched","matrix")
  res
### A integer matrix with 3 rows: clusterID, position and δx1000 and
### as many columns as events in evts is returned. The first row
### contains the number / index of the best template (in the minimal
### Euclidean distance sense), the second row contains the position in
### sampling point of the event's peak and the third row contains the
### jitter value multiplied by 1000 and round to the nearest
### integer. The returned matrix has 3 attributes: call, the matched
### call; data, the symbol of the data object used to create evts and
### templateList, the value of the argument with that name. The
### returned matrix is given an eventsMatched class.
}
###########################################################################
###########################################################################
###########################################################################
"[.eventsMatched" <- function #Subsetting method for eventsMatched objects
### Subsetting method for eventsMatched objects
(x, ##<< an eventsMatchedObject
 i, ##<< should be missing only subsetting of columns allowed.
 j, ##<< index vector for the columns.
 drop = FALSE ##<< logical.
 ) {
  y <- NextMethod("[")
  if (is.matrix(y) && dim(y)[2] >= 1) {
    attr(y,"data") <- attr(x,"data")
    attr(y,"call") <- match.call()
    attr(y,"templateList") <- attr(x,"templateList")
    class(y) <- "eventsMatched"
  }
  y
### An eventsMatched object is returned
}
###########################################################################
###########################################################################
###########################################################################
print.eventsMatched <- function #print method for eventsMatched objects
### Print method for eventsMatched objects.
(x, ##<< an eventsMatchedObject
 ... ##<< not used but required for compatibility with generic method
 ) {
  nc <- dim(x)[2]
  rn <- rownames(x)
  attributes(x) <- NULL
  x <- matrix(x,nc=nc)
  rownames(x) <- rn
  print(x)
### The matrix part of the object is printed.
}
###########################################################################
###########################################################################
###########################################################################
onePerClique <- function #Select one event per clique
### A clique is defined as a sequence of events with inter-event
### intervals smaller than a critical length. This function identifies
### cliques and return the single event giving the best match.
(object, ##<< an eventsMatched object.
 criticalLength=45, ##<< an integer, the number of sampling points
                    ##beyond which events are not considered as
                    ##members of the same clique.
) {
  stopifnot("eventsMatched" %in% class(object))
  rawData <- eval(attr(object,"data"))
  templateFctL <- attr(object,"templateList")
  evtsTimes <- object["position",]
  cliques <- list(list(1))
  nbCliques <- 1
  for (i in 2:length(evtsTimes)) {
    if (evtsTimes[i]-evtsTimes[i-1] <= criticalLength) {
      cliques[[nbCliques]][length(cliques[[nbCliques]])+1] <- i
    } else {
      nbCliques <- nbCliques+1
      cliques[[nbCliques]] <- i
    }
  }
  clusterV <- object["clusterID",]
  δV <- object["δx1000",]/1000
  len <- dim(rawData)[1]
  x <- evalq(x,env=environment(templateFctL[[1]][[1]]))
  selected <- integer(length(cliques))
  for (cIdx in seq(along=cliques)) {
    cliqueMembers <- unlist(cliques[[cIdx]])
    nbMembers <- length(cliqueMembers)
    if (nbMembers > 1) {
      domain <- range(evtsTimes[cliqueMembers])+ c(-1,1)*criticalLength
      refData <- rawData[domain[1]:domain[2],]
      RSS <- sapply(cliqueMembers,
                    function(mIdx) {
                      origin <- clusterV[mIdx]
                      δ <- δV[mIdx]
                      fL <- templateFctL[[origin]]
                      xx <- x+evtsTimes[mIdx]
                      yy <- xx - domain[1] + 1
                      keep <- domain[1] <= xx & xx <= domain[2]
                      pred <- sapply(fL, function(f) f(x[keep]+δ))
                      res <- refData[yy[keep],] - pred
                      sum(res^2)
                    })
      keepMember <- which.min(RSS)
      cliqueMembers <- cliqueMembers[keepMember]
    }
    selected[cIdx] <- cliqueMembers
  }
  selected
### A vector of indexes with the column of object selected as "best" event in each clique.
}
###########################################################################
###########################################################################
###########################################################################
predict.eventsMatched <- function #predict method for eventsMatched objects
### Here prediction is understood as complete data predicition, the
### data being the raw data from which the eventsMatched object was
### derived.
(object, ##<< an eventsMatched object.
 ... ##<< not used but required for compatibility with generic predict method.
 ) {
  rawData <- eval(attr(object,"data"))
  templateFctL <- attr(object,"templateList")
  evtsTimes <- object["position",]
  clusterV <- object["clusterID",]
  δV <- object["δx1000",]/1000
  if ("ts" %in% class(rawData)) {
    freq <- frequency(rawData)
    start <- start(rawData)[1]
    makeTS <- TRUE
  } else {
    makeTS <- FALSE
  }
  len <- dim(rawData)[1]
  idealData <- matrix(0,nr=len,nc=dim(rawData)[2])
  x <- evalq(x,env=environment(templateFctL[[1]][[1]]))
  for (eIdx in seq(along=δV)) {
    origin <- clusterV[eIdx]
    δ <- δV[eIdx]
    fL <- templateFctL[[origin]]
    xx <- x+evtsTimes[eIdx]
    keep <- 1 <= xx & xx <= len
    pred <- sapply(fL, function(f) f(x[keep]+δ))
    idealData[xx[keep],] <- idealData[xx[keep],] + pred
  }
  if (makeTS) ts(idealData,frequency=freq,start=start)
  else idealData
### A vector, matrix, one dimensional or multidimensional time series
### depending on the class of the object to which the data symbol
### attribute of object evaluates. Essentially a each position (second
### row) of object, a template corresponding the cluster given by the
### first row with the jitter specified by the third row is added to
### an ideal "noise free" trace.
}
###########################################################################
###########################################################################
###########################################################################
residuals.eventsMatched <- function #residuals method for eventsMatched object
### Here the residuals are obtained by subtracting from the raw data
### from which the object was derived, the ideal prediction.
(object, ##<< an eventsMatched object.
 ... ##<< not used but required for compatibility with generic predict method.
 ) {
  eval(attr(object,"data")) - predict(object)
### A vector, matrix, one dimensional or multidimensional time series
### depending on the class of the object to which the data symbol
### attribute of object evaluates. Essentially a each position (second
### row) of object, a template corresponding the cluster given by the
### first row with the jitter specified by the third row is subtracted from
### the raw data.
}
###########################################################################
###########################################################################
###########################################################################
fuseMatches <- function #Fuse two eventsMatched objects derived from the same data with the same model
### Given two eventsMatched objects derived from the same data with
### the same model (set of templates), return a new eventsMatched
### object with a single occurrence of each event.
(match1, ##<< an eventsMatched object.
 match2  ##<< an eventsMatched object.
 ) {
  stopifnot("eventsMatched" %in% class(match1))
  stopifnot("eventsMatched" %in% class(match2))
  templateList <- attr(match1,"templateList")
  stopifnot(identical(templateList,attr(match2,"templateList")))
  p1 <- match1["position",]
  p2 <- match2["position",]
  stopifnot(length(p1) >= length(p2))
  data <- attr(match1,"data")
  m1 <- unclass(match1)
  m2 <- unclass(match2)
  res <- cbind(m1,m2)
  positions <- c(p1,p2)
  sIdx <- sort.int(positions,index.return=TRUE)$ix
  positions <- positions[sIdx]
  res <- res[,sIdx]
  attr(res,"call") <- match.call()
  attr(res,"data") <- attr(evts,"data")
  attr(res,"templateList") <- templateList
  class(res) <- c("eventsMatched","matrix")
  res
### An events matched object.
}
###########################################################################
###########################################################################
###########################################################################
