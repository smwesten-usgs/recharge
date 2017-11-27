#' Estimate recharge by means of hydrograph recession displacement (RORA)
#'
#' Estimate recharge by the method of hydrograph recession displacement.
#'
#' @param date the date of the flow, can be either character or class "Date."
#' @param discharge the mean daily discharge for the corresponding date.
#' @param da the drainage area of the basin in square miles.
#' @param recessIndex the recession index, in days per common log cycle.
#'This is typically estimated using \code{recess}.
#' @param minimum_discharge the value to use for the minimum value in \code{discharge}. Any
#'value in \code{discharge} less than the value for \code{minimum_discharge} is set to the
#'value for \code{minimum_discharge}.
#' @param incAnteRec a value to add to the base antecedent recession time,
#'in days. In general, this should always be 0.
#' @param STAID a character string to be used as the station identifier. This is
#'used only for documentation when printing or plotting.
#' @return An object of class rora, which has these components:\cr
#'iyear, the year of the observed streamflow. \cr
#'imonth, the month of the observed streamflow. \cr
#'iday, the day of the observed streamflow. \cr
#'discharge, the observed streamflow. \cr
#'Nobs, the number of observed values of streamflow. \cr
#'iyearst, the starting year of the recharge analysis. \cr
#'iyearen, the ending year of the recharge analysis. \cr
#'minQ, the value of minQ. \cr
#'idiff, the value of incAnteRec. \cr
#'k, the value of recessIndex. \cr
#'te, the time to end day of the recession following a peak in streamflow, in days. \cr
#'ta, the time at the critical time after the previous peak. \cr
#'qp, the streamflow of the peak. \cr
#'qa, the streamflow the critical time after the previous peak. \cr
#'qb, the streamflow at the critical time that would have occurred in the
#'absence of the current and any subsequent peaks. \cr
#'qc, the streamflow at the critical time that would have occurred
#'in the absence of any subsequent peaks. \cr
#'c, the average value for the current peak calculated from the difference
#'between the flow during recession and flow computed from the recession
#'index for each day between the peak and the critical time. \cr
#'delq, the difference in flow between the hypothetical flow at the
#'critical time after the current peak and the hypothetical flow at
#'the critical time after the previous peak. \cr
#'rech, the estimated recharge for the peak, in inches. \cr
#'year, the year of the peak. \cr
#'mon, the month of the peak. \cr
#'day, the day of the peak. \cr
#'npeaks, the number of peaks. \cr
#'itbase, the antecedent flow base time, in days. This is the minimum time
#'from a peak to when the flow can be considered to be ground-water discharge.\cr
#'ierr, the error code. \cr
#'STAID, the station identifier.
#'
#' @useDynLib recharge roraf
#' @export
bf_rora <- function(date, discharge, da, recessIndex,
                 minimum_discharge=0.01, incAnteRec=0, STAID=NULL) {

  discharge <- discharge[ !is.na( discharge ) ]
  date <- date[ !is.na( discharge ) ]

  iyear  <- lubridate::year( date )
  imonth <- lubridate::month( date )
  iday   <- lubridate::day( date )

  N <- length(iyear)

  retval <- .Fortran("roraf",
                     iyear=as.integer(iyear),
                     imonth=as.integer(imonth),
                     iday=as.integer(iday),
                     discharge=as.double(discharge),
                     nobs=as.integer(N),
                     da=as.double(da),
                     iyearst=as.integer( min(iyear) ),
                     iyearen=as.integer(max( iyear ) ),
                     minQ=as.double( minimum_discharge ),
                     idiff=as.integer(incAnteRec),
                     k=as.double(recessIndex),
                     te=integer(N),
                     ta=numeric(N),
                     qp=numeric(N),
                     qa=numeric(N),
                     qb=numeric(N),
                     qc=numeric(N),
                     c=numeric(N),
                     delq=numeric(N),
                     rech=numeric(N),
                     year=integer(N),
                     month=integer(N),
                     day=integer(N),
                     npeaks=integer(1),
                     itbase=integer(1),
                     ierr=integer(1))
  if(retval$ierr < 0)
    warning(paste("Error code from rora: ", retval$ierr))
  if(retval$ierr > 0)
    stop(paste("Error code from rora: ", retval$ierr))
  Npeaks <- retval$npeaks
  length(retval$te) <- Npeaks
  length(retval$ta) <- Npeaks
  length(retval$qp) <- Npeaks
  length(retval$qa) <- Npeaks
  length(retval$qb) <- Npeaks
  length(retval$qc) <- Npeaks
  length(retval$c) <- Npeaks
  length(retval$delq) <- Npeaks
  length(retval$rech) <- Npeaks
  length(retval$year) <- Npeaks
  length(retval$month) <- Npeaks
  length(retval$day) <- Npeaks

  retval$peakdate <- lubridate::ymd(paste(retval$year,retval$month,retval$day,sep="-") )
  retval$delta_t <- numeric( Npeaks )
  retval$delta_t[ 2:Npeaks ] <- as.numeric( retval$peakdate[2:Npeaks] - retval$peakdate[1:Npeaks-1])
  retval$delta_1[1] <- retval$peakdate[1] - date[1]
  retval$baseflow_cfs <- retval$rech / retval$delta_t / 12. * da * 27878400. / 86400.

  if(!is.null(STAID))
    retval$STAID <- STAID
  oldClass(retval) <- "rora"
  return(retval)
}

na2miss <- function( x, to = -99999 ) {

  if( inherits(x,'factor')) {
    levs <- c(levels(x), as.character(to))
    x <- as.vector(x)
    x[is.na(x)] <- to
    return( factor(x, levels=levs))
  }

  x[ is.na(x) ] <- to
  return(x)
}
