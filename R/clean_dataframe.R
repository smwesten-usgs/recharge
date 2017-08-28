clean_dataframe <- function( date, discharge, conc ) {

  temp_df <- data.frame( date=date, discharge=discharge, conc=conc )
  temp_df <- subset( temp_df, !is.na( temp_df$conc ) )

  firstdate <- min(temp_df$date)
  lastdate <- max(temp_df$date)
  complete_dates <- seq( as.Date(firstdate), as.Date(lastdate), by="day" )

  df <- data.frame( date=complete_dates )
  df <- dplyr::left_join(x=df, y=temp_df, by=c("date" = "date") )

  # eliminate any remaining "-99999"s that are present
  df$discharge[ df$discharge < 0. ] <- NA
  df$conc[ df$conc < 0. ] <- NA

  df$discharge <- imputeTS::na.interpolation( df$discharge )
  df$conc      <- imputeTS::na.interpolation( df$conc )

  return( df )

}
