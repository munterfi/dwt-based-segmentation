#' Find peaks
#' 
#' Find peaks (maxima) in a time series.
#'
#' @param x numerical vector taken as a time series
#' @param nups minimum number of increasing steps before a peak is reached
#' @param ndowns minimum number of decreasing steps after the peak
#' @param zero can be `+', `-', or `0'; how to interprete succeeding steps of the same value: increasing, decreasing, or special
#' @param peakpat define a peak as a regular pattern, such as the default pattern ``[+]1,[-]1,''; if a pattern is provided, the parameters nups and ndowns are not taken into account (eg. peakpat = "[+]{2,}[0]*[-]{2,}") 
#' @param minpeakheight the minimum (absolute) height a peak has to have to be recognized as such
#' @param minpeakdistance the minimum distance (in indices) peaks have to have to be counted
#' @param threshold the minimum
#' @param npeaks the number of peaks to return
#' @param sortstr logical; should the peaks be returned sorted in decreasing oreder of their maximum value
#' @param plot logical; should de time series, peaks, their width and pominence be plotted?
#'
#' @return Returns a matrix where each row represents one peak found. The first column gives the height, the second the position/index where the maximum is reached,
#' the third column is the peak width at half prominence and the forth column is the prominence.
#' 
#' @export
#' 
#' @note The original version on this function stems from the \code{pracma} package,
#' but this version provides the same output as the original matlab version.
#'
#' @examples
#' findpeaks(x)
findpeaks <- function(x,nups = 1, ndowns = nups, zero = "0", peakpat = NULL, 
                      minpeakheight = -Inf, minpeakdistance = 1,
                      threshold = 0, npeaks = 0, sortstr = FALSE, plot = FALSE)
{
  stopifnot(is.vector(x, mode="numeric") || length(is.na(x)) == 0)
  if (! zero %in% c('0', '+', '-'))
    stop("Argument 'zero' can only be '0', '+', or '-'.")
  
  # transform x into a "+-+...-+-" character string
  xc <- paste(as.character(sign(diff(x))), collapse="")
  xc <- gsub("1", "+", gsub("-1", "-", xc))
  # transform '0' to zero
  if (zero != '0') xc <- gsub("0", zero, xc)
  
  # generate the peak pattern with no of ups and downs
  if (is.null(peakpat)) {
    peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
  }
  
  # generate and apply the peak pattern
  rc <- gregexpr(peakpat, xc)[[1]]
  if (rc[1] < 0) return(NULL)
  
  # get indices from regular expression parser
  x1 <- rc
  x2 <- rc + attr(rc, "match.length")
  attributes(x1) <- NULL
  attributes(x2) <- NULL
  
  # find index positions and maximum values
  n <- length(x1)
  xv <- xp <- numeric(n)
  for (i in 1:n) {
    xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
    xv[i] <- x[xp[i]]
  }
  
  # eliminate peaks that are too low
  inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >= threshold)
  
  # combine into a matrix format
  X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
  
  # eliminate peaks that are near by
  if (minpeakdistance < 1)
    warning("Handling 'minpeakdistance < 1' is logically not possible.")
  
  # sort according to peak height
  if (sortstr || minpeakdistance > 1) {
    sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
    X <- X[sl, , drop = FALSE]
  }
  
  # return NULL if no peaks
  if (length(X) == 0) return(c())
  
  # find peaks sufficiently distant
  if (minpeakdistance > 1) {
    no_peaks <- nrow(X)
    badpeaks <- rep(FALSE, no_peaks)
    
    # eliminate peaks that are close to bigger peaks
    for (i in 1:no_peaks) {
      ipos <- X[i, 2]
      if (!badpeaks[i]) {
        dpos <- abs(ipos - X[, 2])
        badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
      }
    }
    # select the good peaks
    X <- X[!badpeaks, , drop = FALSE]
  }
  
  # return only the first 'npeaks' peaks
  if (npeaks > 0 && npeaks < nrow(X)) {
    X <- X[1:npeaks, , drop = FALSE]
  }
  
  # sort peaks
  X <- X[order(X[,2]),]
  
  # plot peaks
  if(plot) {
    plot(x, type = "l", col = rgb(0,0.5,1), ylab="")
    points(X[,2], X[,1]+((max(x)-min(x))/100), pch = 25, bg="black")
  }

  # peak prominence and width at half prominence
  for (p in 1:nrow(X)) {
    ind = X[p, 2]
    # peak prominence
    for (l in ind:1) {
      if(x[l] > X[p, 1]) {break}
    }
    for (r in ind:length(x)) {
      if(x[r] > X[p, 1]) {break}
    }
    prom = X[p, 1] - max(min(x[l:ind]), min(x[ind:r]))
    X[p, 4] = prom
    if(plot) {lines(rep(X[p, 2],2), c(X[p, 1]-prom, X[p, 1]), col = rgb(1, 0.6, 0))}
    # peak width at half prominence
    for (l in ind:1) {
      if(x[l] < prom/2) {break}
    }
    for (r in ind:length(x)) {
      if(x[r] < prom/2) {break}
    }
    X[p, 3] = r-l
    if(plot) {lines( c(l, r), rep(prom/2,2), col = rgb(1, 0.6, 0))}
  }
  return(X)
}
