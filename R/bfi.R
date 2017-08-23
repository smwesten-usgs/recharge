#' Baseflow Separation
#'
#' Extract baseflow from a daily streamflow record using the method proposed by
#'the Institute of Hydrology (1980a,b), and as implemented by Wahl and Wahl (1988).
#'
#'The original description of constructing turning points proposed dividing
#'water years into 5-day increments and selecting the minimum flow for each period.
#'Gustard and others (1992) describe using continuous 5-day increments rather than restarting
#'on the water year.
#'
#' @param discharge the daily streamflow to be separated missing values are not permitted
#'within the time specified by \code{Start} and \code{end}.
#' @param date the date for each \code{discharge} value, should be of class "Date." Missing values
#'are not permitted.
#' @param by string describing how to construct turning points: "water year", "calendar_year", or
#'   "continuous".
#' @param f the factor for identifying turning points.
#' @param N the time perod for period for calculating the turning points.
#' @param STAID the station identifier for the data.
#' @references Gustard, A., Bullock, A., and Dixon, J.M., 1992,
#'Low flow estimation in the United Kingdom: Institue of Hydrology
#'Report No. 108, 88 p. and appendixes.
#'
#'Institute of Hydrology, 1980a, Low flow studies:
#'Wallingford, Oxon, United Kingdom, Report No. I, 41 p.
#'
#'Institute of Hydrology, 1980b, Low flow studies: Wallingford, Oxon,
#'United Kingdom, Report No. 3, p. 12- 19.
#'
#'Wahl, K.L., and Wahl, T.L., 1988, BFI—A computer program for determining an index
#' to base flow: US Bureau of Reclamation Water Resources Research Laboratory,
#'  accessed August 25, 2004: U.S. Bureau of Reclamation.
#'
#' @return an object of class "baseflow" and inherits class "data.frame" of the selected data, a data
#'frame of the baseflow information, and other information about the analysis.
#' @keywords baseflow
#' @examples
#'
#'\dontrun{
#'# Process by calendar year as that is the retrieval range
#'ChopBFI <- with(ChoptankFlow, bfi(Flow, datetime, by="calendar year",
#'STAID="01491000"))
#'head(ChopBFI, 20)
#'}
#'@export
bfi <- function(date, discharge, by="water year", f=0.9, N=5L, STAID="Unknown") {

  STAID <- as.character(STAID[1L])
  discharge <- pmax(discharge, 0.0001 ) # Convert 0 to a small number

  if(any(is.na(discharge)))
    stop("Missing discharge values.")
  if(any(diff(as.double(date)) != 1))
    stop("Date data are not continuous.")
  by <- match.arg(by, c("water year", "calendar year", "continuous"))
  if(by == "calendar year") {
    ## Cut pts
    Cut <- c(seq(0, 360, by=N), 366)
    ## Process data by years
    year <- lubridate::year(date)
    dayno <- yday # from lubridate
  } else if(by == "water year") {
    Cut <- c(seq(0, 360, by=N), 366)
    year <- waterYear(date, numeric=TRUE)
    dayno <- function(x) {
      Jul <- as.integer(x)
      lWY <- waterYear(x, numeric=TRUE) - 1L
      baseWY <- as.integer(as.Date(ISOdate(lWY, 9, 30)))
      return(Jul - baseWY)
    }
  } else { # for now continuous
    Cut <- c(seq(0, length(discharge) + N, by=N))
    year <- rep(1L, length(discharge))
    dayno <- function(x) seq(1L, by=1L, length.out=length(x))
  }
  DF <-  data.frame(date=date, discharge=discharge)
  ret <- by(DF, year, function(DF) {
    Jul <- dayno(DF$Dates)
    Grp <- cut(Jul, Cut, labels=FALSE)
    retval <- tapply(DF$Q, Grp, function(x) {
      Min <- min(x)
      Wch <- which(Min == x)[1L]
      return(c(Min, Wch))
    }) # end of lapply
    retval <- do.call('rbind', retval)
    retval[,2L] <- retval[, 2L]+Cut[unique(Grp)]
    return(retval)}
    )
  ## Note: if the first year is a partial year, then need to subtract
  ##  the [starting day number] from the pointer (ret[,2]) for that first
  ##  year. The correction is applied secondarily.
  yeartbl <- cumsum(c(0, table(year))) # The number to add to each pointer
  for(i in seq( length( ret ) ))
    ret[[i]][, 2L] <- ret[[i]][, 2L] + yeartbl[i]
  Jstrt <- dayno(date[1L])
  if(Jstrt > 1L)
    ret[[1L]][, 2L] <- ret[[1L]][, 2L] - Jstrt + 1L
  ## The matrix ret are the minima and the index to the value
  ret <- do.call('rbind', ret)
  ## Now calculate the turning points
  TP <- rep(FALSE, nrow(ret))
  for(i in seq(2L, nrow(ret)-1L)) {
    if(ret[i, 1L] == 0) {
      TP[i] <- TRUE
    } else if(ret[i - 1L, 1L] == 0) {
      TP[i] <- f*ret[i, 1L] <= ret[i+1L, 1L]
    } else if(ret[i + 1L, 1L] == 0) {
      TP[i] <- f*ret[i, 1L] <= ret[i-1L, 1L]
    } else {
      TP[i] <- f*ret[i, 1L] <= min(ret[i-1L, 1L], ret[i+1L, 1L])
    }
  }
  TPdat <- ret[TP,]
  BaseQ <- rep(NA, length=length(Flow))
  for(i in seq(1L, nrow(TPdat)-1L)) {
    Rng <- seq(TPdat[i, 2L], TPdat[i+1L, 2L])
    if(TPdat[i, 1L] == 0 || TPdat[i+1L, 1L] == 0) { # Use linear interpolation
      BaseQ[Rng] <- pmin(Flow[Rng], seq(TPdat[i, 1L], TPdat[i+1L, 1L],
                                          length.out=TPdat[i+1L, 2L] - TPdat[i, 2L]+1L))
    } else
      BaseQ[Rng] <- pmin(Flow[Rng], exp(seq(log(TPdat[i, 1L]), log(TPdat[i+1L, 1L]),
                                            length.out=TPdat[i+1L, 2L] - TPdat[i, 2L]+1L)))
  }
  TurnPt <- rep(" ", length(Flow))
  TurnPt[TPdat[,2]] <- "*"
  retval <- data.frame(Dates=Dates, BaseQ=BaseQ, Flow=Flow, TurnPt=TurnPt)
  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID
  attr(retval, "type") <- "bfi"
  class(retval) <- c("baseflow", "data.frame")
  return(retval)
}

