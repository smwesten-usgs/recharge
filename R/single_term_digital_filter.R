#' Baseflow Separation by use of a single-term digital filter.
#'
#' Extract baseflow from a daily streamflow record using the method described by
#'Nathan and MacMahon (1990).
#'
#' @param date vector of dates corresponding to each \code{discharge}, should be of class "Date."
#' Missing values are not permitted.
#' @param discharge the daily streamflow to be separated missing values are not permitted
#'within the time specified by \code{Start} and \code{end}.
#' @param alpha filter parameter.
#' @param STAID the station identifier for the data.
#'
#' @return an object of class "baseflow" and inherits class "data.frame" of the selected data,
#'a data frame of the baseflow information, and other information about the analysis.
#'
#' @keywords baseflow
#' @examples
#'
#'\dontrun{
#'}
#'@export
bf_single_term_filter <- function(date, discharge, alpha, STAID="Unknown") {

    ## Start of code: initial processing
  STAID <- as.character(STAID[1L])
  discharge <- pmax(discharge, 0.000001 ) # Convert 0 to a small number
  if(any(is.na(discharge)))
    stop("Missing values in discharge vector.")
  if(any(diff(as.double(date)) != 1))
    stop("Date data are not continuous.")

  n <- length( discharge )

  baseflow  <- rep(0.0, n)
  quickflow <- rep(0.0, n)

  forward_run <- function( quickflow, discharge, alpha ) {

    n <- length( discharge )


    for ( i in 2:n) {
      quickflow[i] <- alpha * quickflow[i-1] + ( 1. + alpha ) / 2. * ( discharge[i] - discharge[i-1] )
    }


    quickflow <- pmin( quickflow, discharge )
    quickflow <- pmax( quickflow, rep(0.0, n) )

    return(quickflow)

  }

  reverse_run <- function( quickflow, discharge, alpha ) {

    n <- length( discharge )


    for ( i in (n-1):1) {
      quickflow[i] <- alpha * quickflow[i+1] + ( 1. + alpha ) / 2. * ( discharge[i] - discharge[i+1] )
    }


    quickflow <- pmin( quickflow, discharge )
    quickflow <- pmax( quickflow, rep(0.0, n) )

    return(quickflow)

  }

    quickflow <- forward_run( quickflow, discharge, alpha )
    quickflow <- reverse_run( quickflow, discharge, alpha )
    quickflow <- forward_run( quickflow, discharge, alpha )

    baseflow <- discharge - quickflow

  retval <- data.frame( date=date, discharge=discharge, baseflow=round( baseflow, 3),
                        quickflow = round( quickflow, 3 ) )

  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID
  attr(retval, "type") <- "digital"
  class(retval) <- c("single_term_digital", "data.frame")
  return(retval)
}
