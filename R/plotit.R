plotit <- function( df, wy=NA, BFImax=0.8, alpha=0.925, da=NA, gageid=NA, extra=NA ) {

  par( mar=c(5,4,4,5)+0.1)

  df$Date <- df$date
  df <- addWaterYear( df )

  if ( is.na( wy ) ) {
    df_ss <- df
  } else {
    df_ss <- subset(df , waterYear == wy )
  }

  bf_eck <- with(df_ss,bf_eckhardt_filter( date, discharge, BFImax, alpha ) )
  bf_chem <- with(df_ss,bf_chem_mass_balance( date, discharge, conc ) )
  bf_hysep <- with(df_ss, bf_hysep( date=date, discharge=discharge, da=da ) )
  bf_part <- with(df_ss, bf_part( date=date, discharge=discharge, da=da ) )


  with( bf_chem, plot( date, discharge, type="l", col="blue" ) )
#  plot_flow_conc( date, discharge, conc )
  with( bf_chem, points( date, baseflow, col="red", cex=0.24 ) )
  with( bf_eck, lines( date, baseflow, col="green" ) )
  with( bf_hysep, lines( date, baseflow, col="pink" ) )
  with( bf_part, lines( date, baseflow, col="orange" ) )

  par( new=TRUE )
  plot(df_ss$date, df_ss$conc, axes=F, ylim=c(0,max(df_ss$conc, na.rm=TRUE)), xlab="", ylab="",type="l",
       col="green", main="")
  points(df_ss$date,df_ss$conc,pch=21,col="green", cex=0.3)
  axis(4, ylim=c(0,max(df_ss$conc, na.rm=TRUE)),col="black")
  mtext(4,text="Specific conductance, umhos/cm",line=2)

  recharge_chem <- cfs_to_recharge( date=bf_chem$date, baseflow=bf_chem$baseflow, da=da)
  recharge_eck <-  cfs_to_recharge( date=bf_eck$date, baseflow=bf_eck$baseflow, da=da)
  recharge_hysep <-  cfs_to_recharge( date=bf_hysep$date, baseflow=bf_hysep$baseflow, da=da)
  recharge_part <-  cfs_to_recharge( date=bf_part$date, baseflow=bf_part$baseflow, da=da)

  n_years <- nrow( recharge_chem ) / 365.25

  r_chem <- sum( recharge_chem$recharge_in ) / n_years
  r_eck <- sum( recharge_eck$recharge_in ) / n_years
  r_hysep <- sum( recharge_hysep$recharge_in ) / n_years
  r_part <- sum( recharge_part$recharge_in ) / n_years

  title_txt <-paste("Water Year ", wy, ", Eckhart filter parameters: BFImax=",BFImax,
                    "  alpha=", alpha, sep="" )

  if ( !is.na(gageid) ) {
    title_txt <- paste( gageid, "  | ",title_txt, sep="" )
  }

  title( main=title_txt,
         sub=paste("recharge (chem): ",round(r_chem,2),
                   "  recharge (Eckhardt): ", round(r_eck,2),
                   "  recharge (HYSEP): ", round(r_hysep,2),
                   "  recharge (PART): ", round(r_part,2),
                   sub=""))
  mtext(paste("Drainage area: ",da," sq mi",sep=""), side=4, line=3)
  if ( !is.na( extra ) )   mtext(paste(extra), side=4, line=4)

}
