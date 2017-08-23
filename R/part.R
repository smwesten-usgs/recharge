#' Baseflow Separation
#'
#' Extract baseflow from a daily streamflow record using the method described by
#'Rutledge (1998).
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
#' @param STAID the station identifier for the data.
#' @references Rutledge, A.T., 1998, Computer programs for describing the recession of
#'ground-water discharge and for estimating mean ground-water recharge and discharge
#'from streamflow records---Update: U.S. geological Survey Water-Resources Investigations
#'Report 98-4148. 43 p.
#'
#' @return an object of class "baseflow" and inherits class "data.frame" of the selected data,
#'a data frame of the baseflow information, and other information about the analysis.
#' @note The estimates from this routine will occasionally differ slightly from the original
#'FORTRAN code from Rutledge (1998). It is expected that those small differences are due to
#'rounding differences affecting comparison between numbers that differ by very small amounts.
#'There are also differences in the initial baseflow estimates between this version and the 
#'original version. Those differences are due to slightly different initialization routines.
#'
#'The estimate of baseflow from this routine is computed from a linear interpolation of the
#'\code{Part1} and \code{Part2} estimates rather than the curvilinear interpolation that was 
#'proposed without detials in Rutledge (1998).
#' @keywords baseflow
#' @examples
#'
#'\dontrun{
#'library(smwrData)
#'data(ChoptankFlow)
#'# Process by calendar year as that is the retrieval range
#'ChopPart <- with(ChoptankFlow, part(Flow, datetime, da=113,
#'STAID="01491000"))
#'ChopPart
#'}
#'@export
part <- function(Flow, Dates, Start=NULL, End=NULL, da,
                 STAID="Unknown") {
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
  Nact <- max(da^0.2, 1)
  N <- as.integer(ceiling(Nact))
  NF <- max(N-1L, 1L)
  NC <- max(N, 2L)
  NC1 <- NC + 1L
  ## From the flow chart in Rutledge, with additions for 0 flows
  ## The variable suffixes are F is the floor of Nact, C is the
  ## ceiling of Nact, and C1 is the ceiling plus 1. These correspond
  ## to the three values of N in Rutledge.
  ##
  # Step 1 set up data (Flow already done)
  ## ALLGW is set to logical: TRUE (*) and FALSE (0)
  ALLGWF <- ALLGWC <- ALLGWC1 <- rep(FALSE, length(Flow))
  BaseQF <- BaseQC <- BaseQC1 <- rep(NA_real_, length(Flow))
  # Step 2 Recored all GW flow where antecendent recession OK
  DiffQ <- c(0, diff(Flow))
  AnteF <- na2miss(filter(DiffQ <= 0, rep(1, NF), sides=1), 0)
  AnteC <- na2miss(filter(DiffQ <= 0, rep(1, NC), sides=1), 0)
  AnteC1 <- na2miss(filter(DiffQ <= 0, rep(1, NC1), sides=1), 0)
  ALLGWF <- ifelse(AnteF == NF, TRUE, ALLGWF)
  BaseQF <- ifelse(ALLGWF, Flow, BaseQF)
  ALLGWC <- ifelse(AnteC == NC, TRUE, ALLGWC)
  BaseQC <- ifelse(ALLGWC, Flow, BaseQC)
  ALLGWC1 <- ifelse(AnteC1 == NC1, TRUE, ALLGWC1)
  BaseQC1 <- ifelse(ALLGWC1, Flow, BaseQC1)
  # Step 3 Revise all GW where necessary
  CkQ <- (Flow > 1e-9) & (Flow/shiftData(Flow, k=-1, fill=1) > 1.258925)
  ALLGWF <- ifelse(ALLGWF & CkQ, FALSE, ALLGWF)
  ALLGWC <- ifelse(ALLGWC & CkQ, FALSE, ALLGWC)
  ALLGWC1 <- ifelse(ALLGWC1 & CkQ, FALSE, ALLGWC1)
  # Step 4 Interpolate Baseflows
  Seq <- seq(length(Flow))
  BaseQF <- exp(approx(Seq[ALLGWF], log(Flow[ALLGWF]), xout=Seq, rule=2)$y)
  BaseQC <- exp(approx(Seq[ALLGWC], log(Flow[ALLGWC]), xout=Seq, rule=2)$y)
  BaseQC1 <- exp(approx(Seq[ALLGWC1], log(Flow[ALLGWC1]), xout=Seq, rule=2)$y)
  # Steps 5, 6, and 4 for each F, C, C1
  while(any(CkQ <- (BaseQF > Flow + 0.000001))) { # Avoid rounding problems
    CkQ <- CkQ & !ALLGWF # The trouble makers
    Ck0 <- eventNum(!ALLGWF, reset=TRUE) # Each block of !ALLGW
    CkE <- unique(Ck0[CkQ])
    ## Find the largest ratio (log difference) in each block of !ALLGW
    for(i in CkE) {
      Sel <- which(Ck0 == i) # Select from this group
      MaxR <- BaseQF[Sel]/Flow[Sel]
      Pck <- which.max(MaxR)
      ALLGWF[Sel[Pck]] <- TRUE
      BaseQF[Sel[Pck]] <- Flow[Sel[Pck]]
    }
    ## Redo 4
    BaseQF <- exp(approx(Seq[ALLGWF], log(Flow[ALLGWF]), xout=Seq, rule=2)$y)
    BaseQF <- ifelse(BaseQF < 1e-6, 0, BaseQF) # Clean up
  }
  
  while(any(CkQ <- (BaseQC > Flow + 0.000001))) { # Avoid rounding problems
    CkQ <- CkQ & !ALLGWC # The trouble makers
    Ck0 <- eventNum(!ALLGWC, reset=TRUE) # Each block of !ALLGW
    CkE <- unique(Ck0[CkQ])
    ## Find the largest ratio (log difference) in each block of !ALLGW
    for(i in CkE) {
      Sel <- which(Ck0 == i) # Select from this group
      MaxR <- BaseQC[Sel]/Flow[Sel]
      Pck <- which.max(MaxR)
      ALLGWC[Sel[Pck]] <- TRUE
      BaseQC[Sel[Pck]] <- Flow[Sel[Pck]]
    }
    ## Redo 4
    BaseQC <- exp(approx(Seq[ALLGWC], log(Flow[ALLGWC]), xout=Seq, rule=2)$y)
    BaseQC <- ifelse(BaseQC < 1e-6, 0, BaseQC)
  }
  
  while(any(CkQ <- (BaseQC1 > Flow + 0.000001))) { # Avoid rounding problems
    CkQ <- CkQ & !ALLGWC1 # The trouble makers
    Ck0 <- eventNum(!ALLGWC1, reset=TRUE) # Each block of !ALLGW
    CkE <- unique(Ck0[CkQ])
    ## Find the largest ratio (log difference) in each block of !ALLGW
    for(i in CkE) {
      Sel <- which(Ck0 == i) # Select from this group
      MaxR <- BaseQC1[Sel]/Flow[Sel]
      Pck <- which.max(MaxR)
      ALLGWC1[Sel[Pck]] <- TRUE
      BaseQC1[Sel[Pck]] <- Flow[Sel[Pck]]
    }
    ## Redo 4
    BaseQC1 <- exp(approx(Seq[ALLGWC1], log(Flow[ALLGWC1]), xout=Seq, rule=2)$y)
    BaseQC1 <- ifelse(BaseQC1 < 1e-6, 0, BaseQC1)
  }
  ## Wrap up
  ## Compute the linear interpolation of baseflow
  Ffact <- NC - Nact # Must be between 0 and 1
  BaseQ <- BaseQF*Ffact + BaseQC*(1-Ffact)
  retval <- data.frame(Dates=Dates, BaseQ=round(BaseQ, 3L), 
                       Flow=round(Flow, 3L), Part1=BaseQF, Part2=BaseQC, 
                       Part3=BaseQC1)
  if(!is.null(STAID))
    attr(retval, "STAID") <- STAID
  attr(retval, "type") <- "part"
  class(retval) <- c("baseflow", "data.frame")
  return(retval)
}
  
