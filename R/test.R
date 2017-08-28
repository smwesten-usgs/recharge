find_optimum_filter_params <- function( date, discharge, conc, BFImax, alpha ) {

  df <- data.frame( date=date, discharge=discharge, conc=conc )
  res <- optim( par=c(0.8 ), fn=update_objective_function, df=df, method="Brent",
                lower=0.1, upper=0.9)

  update_objective_function <- function( params, df ) {

    BFImax <- params[1]
    alpha = 0.925

    bf_chem <- bf_chem_mass_balance( date=df$date, discharge = df$discharge, C = df$conc )
    bf_eckhardt <- bf_eckhardt_filter(date=df$date, discharge=df$discharge,
                                      BFImax=BFImax, alpha=alpha)

    sse <- sum( (bf_chem$baseflow - bf_eckhardt$baseflow )^2 )

    return( sse )

  }



  bf_chem <- bf_chem_mass_balance( date=df$date, discharge = df$discharge, C = df$conc )

  bf_eckhardt <- bf_eckhardt_filter(date=df$date, discharge=df$discharge,
                                    BFImax=0.8, alpha=0.925 )

  res$bf_chem <- bf_chem
  res$bf_eckhardt <- bf_eckhardt
  sse <- sum( (bf_chem$baseflow - bf_eckhardt$baseflow )^2 )

  return( res )

}

