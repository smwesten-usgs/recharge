#' Recessions
#'
#' Identify surface-water recessions
#'
#' Note that zero flows are set to a value so the the common log is -2.5 in order to
#'easily process recessions that go to 0 flow.
#'
#' The scaled RMSE of the recession regression is the root-mean-squared error divided
#'by the square root of the recession index. Large values of the scaled RMSE indicate
#'a lack of linearity, so \code{check.srmse} can be used as a filter to remove poor fits,
#'but at some risk of rejecting acceptable recessions. The default value of 0.1, seems
#'to be provide good balance of rejection for poor fits and acceptance of good fits.
#'
#' @param date vector of dates for each \code{x}, should be of class "Date." Missing values
#'are not permitted.
#' @param discharge the streamflow data to be analyzed. Missing values are not permitted.
#' @param by the months to subselect for recessions.
#' @param min.duration the minimum duration for a recession to be selected.
#' @param max.duration the maximum duration for a recession to be selected.
#' @param check.srmse reject the recession if the scaled RMSE of the recesssion
#'regression exceeds this value. See \bold{details}.
#' @param STAID the station identifier for the data.
#' @return an object of class "recess" and inherits class "data.frame" of the selected data, a data
#'frame of the recession information, and other information about the analysis.
#' @keywords recession
#' @importFrom robust lmRob
#' @examples
#'
#'\dontrun{
#'library(smwrData)
#'data(ChoptankFlow)
#'with(ChoptankFlow, recess(Flow, datetime, STAID="0191000"))
#'}
#' @export
recess <- function(date, discharge, by=NULL,
                   min.duration=10, max.duration=300,
									 check.srmse=0.1, STAID="Unknown") {
  ## Coding history:
  ##    2006May18 DLLorenz Initial coding.
  ##    2006Jul24 DLLorenz Added recess function to create a recess object
  ##    2006Jul26 DLLorenz Various enhancements
  ##    2006Jul27 DLLorenz Refinements to better agree with the original
  ##    2006Sep08 DLLorenz Standardization, and rename to recess2, as per
  ##                       suggestion of Al Rutledge.
  ##    2006Sep13 DLLorenz Bug fix
  ##    2008May20 DLLorenz remove NAs from data set
  ##    2013Apr26 DLLorenz Begin conversion to R
	##    2013Jun13 DLLorenz Added max.duration to limit uncharacteristic recessions
  ##
  ## Function needed for this routine:
  eventSel <- function(x, seqno, minseq=10) {
    ## Select the starting and ending value of sequences longer than a
    ## specified minimum
    NR <- length(x)
    if(NR < minseq)
      return(NULL)
    if(seqno[1L] == 0)
      return(NULL)
    return(c(Start=x[1L], End=x[NR], Length=NR))
  }

  recessSel <- function(x, rec, dur) {
    x2 <- seq(along=x)
    x2 <- tapply(x2, rec, function(y, x, rec, dur) eventSel(x[y], rec[y], dur), x=x, rec=rec, dur=dur)
    x2.recno <- as.integer(names(x2)[!sapply(x2, is.null)])
    if(length(x2.recno) == 0)
      stop("No recessions detected")
    x2 <- do.call("rbind", x2)
    x2 <- as.data.frame(x2)
    names(x2) <- c("StartQ", "EndQ", "Length")
    x2$Index <- x2.recno
    return(x2)
  }
  ## Supress warnings
  warn <- options("warn")
  options(warn=-1)
  ## Start of code: initial processing
  STAID <- as.character(STAID[1L])
  if(is.null(by))
    by <- month.abb
  if(is.numeric(by))
    months2Sel <- month.abb[by]
  else
    months2Sel <- by
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
  Flow <- pmax(Flow[sel], 10^(-2.5)) # Convert 0 to a small number
  if(any(is.na(Flow)))
    stop("Missing values between ", Start, " and ", End)
  if(any(diff(as.double(Dates)) != 1))
    stop("Date data are not continuous between Start and End")
  Recess <- eventNum(c(NA, diff(Flow)) <= 0, reset=TRUE)
  Sel <- recessSel(Flow, Recess, dur=min.duration)
  ## Add first estimates of recession indexes
  Index <- Sel$Index
  Len <- Sel$Length
  logMeanQ <- K <- Icept <- DaysFrom1 <- double(length(Index))
  Date <- as.Date(rep(NA, length(Index)))
  SRMSE <- double(length(Index))
  keep <- rep(TRUE, length(Index))

  ## From basic structure of recessions assume that only contiguous points used
  First <- Last <- integer(length(Index))
  for(i in seq(along=Index)) {
    rec <- log10(Flow[Recess == Index[i]])
    Time <- seq(Len[i])
    Date[i] <- Dates[Recess == Index[i]][1L]
    ## Test if month in selected months and flows not all equal and length
    ## less than max.duration
    if(months(Date[i], abbreviate=TRUE) %in% months2Sel &&
    	 	diff(range(Flow)) > 0.001 && Len[i] <= max.duration) {
      ## Use robust regression to identify the most linear part of the recession
      ## Trap flat lines
      if(all(rec == rec[1L])) {
        coefs <- c(0,0)
      } else {
        model <- try(lmRob(Time ~ rec))
        if(class(model) != "try-error")
          temp <- which(model$M.weights > 0)
        else
          temp <- Time # use complete sequence as default
        ## Use lsfit to get results to match the recess program
        tmp <- lsfit(rec[temp], Time[temp])
        coefs <- tmp$coef
        if(coefs[2L] == 0)
        	SRMSE[i] <- Inf
        else
        	SRMSE[i] <- sqrt(mean(tmp$residuals^2)/(-coefs[2L]))
        if(SRMSE[i] > check.srmse)
        	coefs[2L] <- 0 # Rejected immediately below
      }
      if(coefs[2L] == 0) { # Trap bad regression
        keep[i] <- FALSE
      } else { ## The mean log flow of the data used in the regression
        logMeanQ[i] <- mean(rec[temp])
        temp <-  range(temp)
        ## The first point used = START in the x file
        First[i] <- temp[1L]
        ## The last point used = END in the x file
        Last[i] <- temp[2L]
        ## the recesssion index
        K[i] <- -coefs[2L]
        ## An estimate of the time, in days from logMeanQ to 1 flow unit
        DaysFrom1[i] <- coefs[1L] - (First[i] + Last[i]) / 2.0
        ## Compute intercept when x-axis is days from recession
        Icept[i] <- logMeanQ[i] - mean(temp) / coefs[2L]
      }
    } # End if months
    else # Drop recession
      keep[i] <- FALSE
  }
  ## Restore warnings and finish up
  options(warn)
  Recessions <- cbind(Sel, Date=Date, logMeanQ=logMeanQ, K=K, Icept=Icept,
                      DaysFrom1=DaysFrom1, First=First, Last=Last, SRMSE=SRMSE)[keep,]
  if(nrow(Recessions) > 0L)
    rownames(Recessions) <- as.character(seq(nrow(Recessions)))
  retval <- list(Data=na.omit(data.frame(Dates=Dates, Flow=Flow, Recess=Recess)),
                 Recessions=Recessions, Start=Start, End=End, min.duration=min.duration,
                 STAID=STAID, months2Sel=months2Sel)
  attr(retval, "Confirmed") <- FALSE
  oldClass(retval) <- "recess"
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

