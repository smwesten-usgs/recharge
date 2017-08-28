# survey the available conductance data at gages within MAP

library(dataRetrieval)
library(dplyr)

download_inventory <- FALSE
download_data <- FALSE

mystates <- c("LA","MO","AR","MS","TN")

if ( download_inventory ) {
  for ( n in 1:length(mystates) ) {
    myConductivitySites  <- dataRetrieval::whatNWISdata( parameterCd=c("00095"), stateCd=mystates[n])
    myDischargeSites <-  dataRetrieval::whatNWISdata( parameterCd=c("00060"), stateCd=mystates[n])
    assign( paste("sites_conductivity_NWIS_",mystates[n],sep=""), myConductivitySites )
    assign( paste("sites_discharge_NWIS_",mystates[n],sep=""), myDischargeSites )
    if ( n == 1 ) {
      combinedConductivitySites <- myConductivitySites
      combinedDischargeSites <- myDischargeSites
    } else {
      combinedConductivitySites <- dplyr::bind_rows( combinedConductivitySites, myConductivitySites )
      combinedDischargeSites <- dplyr::bind_rows( combinedDischargeSites, myDischargeSites )
    }

  }
}

combinedSites <- dplyr::left_join( x=combinedDischargeSites, y=combinedConductivitySites,
                                   by=c("site_no" = "site_no") )

  cut1 <- subset( combinedSites, ( count_nu.x > 1000 ) & ( count_nu.y > 1000 )
                                 & ( site_tp_cd.x == 'ST' )
                                 & ( stat_cd.x=='00003' )
                                 & ( stat_cd.y=='00003' ) )


  gages1 <- dplyr::left_join( x=cut1, y=hydromod_dams, by=c("site_no" = "STAID") )
  gages2 <- dplyr::left_join( x=gages1, y=basin_id, by=c("site_no" = "STAID") )
  gages3 <- dplyr::left_join( x=gages2, y=bound_qa, by=c("site_no" = "STAID") )
  gages4 <- dplyr::left_join( x=gages3, y=soils, by=c("site_no" = "STAID") )
  gages5 <- dplyr::left_join( x=gages4, y=lc06_basin, by=c("site_no" = "STAID") )

  gages6 <- subset( gages5, !is.na( DRAIN_SQKM.x ) )

  #gages_w_dam_info <- subset(gages5, ! is.na( NDAMS_2009 ) )

  gages_w_dam_info <- gages6


if ( download_data ) {
  my_data <- list()
}

for ( n in 1:nrow( gages_w_dam_info ) ) {

  gageid <- gages_w_dam_info$site_no[n]
  da_sq_mi <- as.numeric( gages_w_dam_info$DRAIN_SQKM.x[n] ) * 0.3861
  df_name <- paste("daily_data_",gageid,"_df", sep="" )

  if ( download_data ) {
    day_data <- readNWISdv(siteNumbers=gageid, parameterCd = c("00060","00095") )
    day_data <- renameNWISColumns( day_data )
    day_data <- clean_dataframe( day_data$Date, day_data$Flow, day_data$SpecCond )
    assign( df_name, day_data )
    my_data[[n]] <- get( df_name )
  }

  BFImax <- ( as.numeric(gages_w_dam_info$HGB[n])/100*.85 + as.numeric(gages_w_dam_info$HGA[n])/100*1.25
            + as.numeric(gages_w_dam_info$HGC[n])/100*0.5 + as.numeric(gages_w_dam_info$HGD[n])/100*0.25 )
  df <- my_data[[n]]
  plotit( df, BFImax=BFImax,
          alpha=0.1, da=round(da_sq_mi,1), gageid=gageid,
          extra=paste(" HGA: ", gages_w_dam_info$HGA[n],
                      " HGB: ", gages_w_dam_info$HGB[n],
                      " HGC: ", gages_w_dam_info$HGC[n],
                      " HGD: ", gages_w_dam_info$HGD[n] ) )

}

