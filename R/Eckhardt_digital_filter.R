#' Baseflow Separation by use of an Eckhardt Filter
#'
#' Extract baseflow from a daily streamflow record using the method described by
#'Eckhardt (2005).
#'
#' @param date vector of dates corresponding to each \code{discharge}, should be of class "Date."
#' Missing values are not permitted.
#' @param discharge the daily streamflow to be separated missing values are not permitted
#'within the time specified by \code{Start} and \code{end}.
#' @param BFImax maximum BFI to use in the filtering.
#' @param alpha filter parameter.
#' @param STAID the station identifier for the data.
#' @references Eckhardt, K., 2005, How to construct recursive digital filters for
#'  baseflow separation: Hydrological Processes, v. 19, no. 2, p. 507–515.
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
bf_eckhardt_filter <- function(date, discharge, BFImax, alpha, STAID="Unknown") {

    ## Start of code: initial processing
  STAID <- as.character(STAID[1L])
  discharge <- pmax(discharge, 0.000001 ) # Convert 0 to a small number
  if(any(is.na(discharge)))
    stop("Missing values in discharge vector.")
  if(any(diff(as.double(date)) != 1))
    stop("Date data are not continuous.")

  n <- length( discharge )

  baseflow  <- numeric( n )
  quickflow <- numeric( n )

  baseflow <- rep( discharge[1], n )

  for ( i in 2:n) {
    baseflow[i] <-   ( ( ( 1. - BFImax ) * alpha * baseflow[i-1] )
                      + ( BFImax * ( 1. - alpha ) * discharge[i] ) ) / ( 1. - alpha * BFImax )
  }

  #bf[i] <-(((1 - BFI)* a* bf[i-1]) + ((1-a)* BFI* discharge [i])) /(1- a*BFI)

  baseflow <- pmin( baseflow, discharge )
  baseflow <- pmax( baseflow, 0.0 )


  quickflow <- discharge - baseflow

  retval <- data.frame( date=date, discharge=discharge, baseflow=baseflow, quickflow=quickflow )

  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID
  attr(retval, "type") <- "part"
  class(retval) <- c("Eckhardt", "data.frame")
  return(retval)
}
