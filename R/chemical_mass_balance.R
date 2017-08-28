#' Baseflow Separation by use of a chemical mass balance
#'
#' Extract baseflow from a daily streamflow record using conductance or other continuously
#' monitored conservative chemicals as markers of baseflow and surface flow.
#'
#' @param date vector of dates corresponding to each \code{discharge}, should be of class "Date."
#' Missing values are not permitted.
#' @param discharge the daily streamflow to be separated missing values are not permitted
#' @param C daily chemical concentration values
#' @param C_baseflow assumed chemical concentration associated with water derived from baseflow
#' @param C_runoff assumed chemical concentration associated with water deriver from surface flow
#' @param STAID the station identifier for the data.
#' @references
#' @return an object of class "baseflow" and inherits class "data.frame" of the selected data,
#'a data frame of the baseflow information, and other information about the analysis.
#'
#' @keywords baseflow
#' @examples
#'
#'\dontrun{
#'}
#'@export
bf_chem_mass_balance <- function(date, discharge, C,
                                 C_baseflow=max(C),
                                 C_runoff=min(C),
                                 STAID="Unknown") {

    ## Start of code: initial processing
  STAID <- as.character(STAID[1L])
  discharge <- pmax(discharge, 0.000001 ) # Convert 0 to a small number
  if(any(is.na(discharge)))
    stop("Missing values in discharge vector.")
  if(any(is.na(C)))
    stop("Missing values in chemical concentration vector.")
  if(any(diff(as.double(date)) != 1))
    stop("Date data are not continuous.")

  n <- length( discharge )

  baseflow  <- numeric( n )
  quickflow <- numeric( n )

  q90 <- quantile( discharge, 0.1 )
  C90 <- quantile( C, 0.90 )
  C_baseflow <- min( median( C[ discharge <= q90 ] ), C90 )

  baseflow <- discharge * ( C - C_runoff ) / ( C_baseflow - C_runoff )

  baseflow <- pmin( baseflow, discharge )
  baseflow <- pmax( baseflow, 0.0 )


  quickflow <- discharge - baseflow

  retval <- data.frame( date=date, discharge=discharge, baseflow=round( baseflow, 3),
                        quickflow = round( quickflow, 3 ) )

  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID
  attr(retval, "type") <- "part"
  class(retval) <- c("CMB", "data.frame")
  return(retval)
}
