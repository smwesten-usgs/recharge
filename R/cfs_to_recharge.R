#' Convert baseflow to a recharge estimate for a basin
#'
#' Convert the daily baseflow in cfs to a net infiltration (recharge) estimate
#' in inches over a basin.
#'
#'@param date the date corresponding to the baseflow (class "Date").
#'@param baseflow the mean daily baseflow for the corresponding date.
#'@param da the drainage area of the basin in square miles.
#'@return A data frame containing the following fields:
#'    \item{date}{date corresponding to the baseflow value}
#'    \item{baseflow}{baseflow value in cfs}
#'    \item{month}{month of year (1-12) of the basteflow value}
#'    \item{day}{day of month (1-31) associated with baseflow value}
#'    \item{year}{calendar year associated with the baseflow value}
#'    \item{wy}{water year associated with the baseflow value}
#'    \item{total_daily_recharge_ft3}{daily recharge in cubic feet}
#'    \item{recharge_in}{daily recharge over basin area, in inches}
#'@export
cfs_to_recharge <- function( date, baseflow, da ) {

  df <- data.frame( date=date, baseflow=baseflow, Date=date )
  df$month <- lubridate::month(df$date)
  df$day   <- lubridate::day(df$date)
  df$year  <- lubridate::year(df$date)
  df$wy    <- dataRetrieval::addWaterYear( df )
  df$Date <- NULL

  df$total_daily_recharge_ft3 <- df$baseflow * 86400

  da_sqft <- ( sqrt( da ) * 5280 )^2

  df$recharge_in <- df$total_daily_recharge_ft3 / da_sqft * 12.

  return( df )

}
