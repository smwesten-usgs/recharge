#' Baseflow Separation
#'
#' Extract baseflow from a daily streamflow record using the method described by
#'Sloto and Crouse (1996).
#'
#' @param Flow the daily streamflow to be separated missing values are not permitted
#'within the time specified by \code{Start} and \code{end}.
#' @param Dates the date for each \code{x}, should be of class "Date." Missing values
#'are not permitted.
#' @param Start the start date for the analysis, can be either a character string or
#'class "Date."
#' @param End the end date for the analysis, can be either a character string or
#'class "Date."
#' @param da the drainage area of the basin in square miles.
#' @param select a character string indicating which method to use for the
#'baseflow in the output dataset. Must be one of "sliding," "local minimum,"
#'or "fixed." Onle the first letter is required.
#' @param STAID the station identifier for the data.
#' @references Sloto, R.A. and Crouse, M.Y., 1996, HYSEP: 
#'A COMPUTER PROGRAM FOR STREAMFLOW HYDROGRAPH SEPARATION AND ANALYSIS:
#'U.S. geological Survey Water-Resources Investigations
#'Report 96-4040. 46 p.
#'
#' @return an object of class "baseflow" and inherits class "data.frame" of the selected data,
#'a data frame of the baseflow information, and other information about the analysis.
#'
#' @keywords baseflow
#' @examples
#'
#'\dontrun{
#'library(smwrData)
#'data(ChoptankFlow)
#'# Process by calendar year as that is the retrieval range
#'ChopPart <- with(ChoptankFlow, hysep(Flow, datetime, da=113,
#'STAID="01491000"))
#'ChopPart
#'}
#'@export
hysep <- function(Flow, Dates, Start=NULL, End=NULL, da,
                 select="sliding", STAID="Unknown") {
  ## Start of code: initial processing
  STAID <- as.character(STAID[1L])
  if(is.null(Start))
    Start <- Dates[1L]
  else if(is.character(Start))
    Start <- as.Date(Start)
  if(is.null(End))
    End <- Dates[length(Dates)]
  else if(is.character(End))
    End <- as.Date(End)
  sel <- (Dates >= Start) & (Dates <= End)
  Dates <- Dates[sel]
  Flow <- pmax(Flow[sel], 1e-99) # Convert 0 to a small number
  if(any(is.na(Flow)))
    stop("Missing values between ", Start, " and ", End)
  if(any(diff(as.double(Dates)) != 1))
    stop("Date data are not continuous between Start and End")
  select <- match.arg(select, c("sliding", "local minimum", "fixed"))
  Nact <- max(da^0.2, 1)
  N2star <- max((((2*Nact) ) %/% 2)*2 + 1, 3)
  ## Set up for fixed--construct intervals of length N2star
  Nobs <- length(Flow)
  Ngrp <- ceiling(Nobs / N2star)
  Grps <- inverse.rle(list(lengths=rep(N2star, Ngrp), values=seq(Ngrp)))
  length(Grps) <- Nobs # Truncate if necessary
  ## Compute the fixed method
  Mins <- tapply(Flow, Grps, min)
  Fixed <- Mins[as.character(Grps)]
  ## Now the sliding method
  Slide <- sapply(seq(N2star, Nobs), function(i)
    min(Flow[seq(i - N2star + 1L, i)])
  )
  SlB <- Slide[1L]
  SlE <- Slide[length(Slide)]
  Nfil <- (N2star - 1L) / 2L
  Slide <- c(rep(SlB, Nfil), Slide, rep(SlE, Nfil))
  ## And the local minimum
  Mid <- as.integer((N2star) / 2)
  LocMin <- sapply(seq(N2star, Nobs), function(i)
    min(Flow[seq(i - N2star + 1L, i)]) == Flow[i - Mid]
  )
  LocMin <- c(rep(FALSE, Nfil), LocMin, rep(FALSE, Nfil))
  ## Need to trap short periods where only 1 local minimum
  if(sum(LocMin) == 1L) {
    warning("Only one local minimum in calibration period")
    LocMin <- pmax(Flow[LocMin], 0.01)
  } else
    LocMin <- exp(approx(which(LocMin), log(pmax(Flow[LocMin], 0.01)), xout=seq(Nobs), rule=2)$y)
  LocMin <- pmin(Flow, LocMin) # recover 0s and tails
  if(select == "fixed")
    BaseQ <- Fixed
  else if(select == "sliding")
    BaseQ <- Slide
  else 
    BaseQ <- LocMin
  retval <- data.frame(Dates=Dates, BaseQ=round(BaseQ, 3L), 
                       Flow=Flow, Fixed=Fixed, Sliding=Slide, 
                       LocalMin=LocMin)
  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID
  attr(retval, "type") <- "hysep"
  class(retval) <- c("baseflow", "data.frame")
  return(retval)
}
