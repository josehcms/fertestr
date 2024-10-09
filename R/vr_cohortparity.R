#' Retrieve Zaba parameters for relational Gompertz model
#'
#' @param x vector of ages to which the parameters must be returned
#' @return a data frame with columns
#' `age` - reference age group;
#' `Fxs` - cumulative fertility for age group for the standard Zaba fertility distribution;
#' `Yxs` - gompit of Fxs;
#' `Fxs_Ratio` - ratio of subsequent `Fxs_Ratio`
#' `Phi`
#' `Phi1`
#' `Phi2`
#' `ex` - parameter for modeling zx - ex = gx
#' `gx` - parameter for modeling zx - ex = gx
#'
#' @keywords internal

get_zaba_params =
  function( x = seq( 15, 50, 5 ) ){

    zaba = DemoToolsData::std.zaba

    std_params =
      data.frame( age = x,
                  Fxs   = NA,
                  Yxs   = NA,
                  Fxs_Ratio = NA,
                  Phi   = NA,
                  Phi1  = NA,
                  Phi2  = NA,
                  ex    = NA,
                  gx    = NA )


    std_params$Fxs = zaba[ zaba$age %in% seq( 15, 50, 5 ), "Fx" ]
    std_params$Yxs = -log( - log( std_params$Fxs ) )

    for( i in 1:7 ){
      std_params$Fxs_Ratio[ i ] = std_params$Fxs[ i ] / std_params$Fxs[ i + 1 ]
    }

    std_params$Phi = -log( -log( std_params$Fxs_Ratio ) )

    for( i in 1:7 ){
      std_params$Phi1[ i ] =
        ( ( std_params$Yxs[ i + 1 ] * exp( - std_params$Yxs[ i + 1 ] ) ) -
            ( std_params$Yxs[ i ] * exp( - std_params$Yxs[ i ] ) ) ) /
        log( std_params$Fxs_Ratio[ i ] )
    }


    for( i in 2:4 ){
      std_params$Phi2[ i ] =
        ( ( std_params$Yxs[ i ] - std_params$Yxs[ i + 1 ] )^2 *
            exp( std_params$Yxs[ i ] + std_params$Yxs[ i + 1 ] ) ) /
        ( exp( std_params$Yxs[ i ] ) - exp( std_params$Yxs[ i + 1 ] ) )^2
    }

    c = mean( std_params$Phi2, na.rm = T )

    std_params$ex = std_params$Phi - std_params$Phi1
    std_params$gx = std_params$Phi1

    return( std_params )

  }

#' Calculate age-specific fertility rates for single ages from a given total fertility level
#'
#' @param TF total fertility level
#' @param alpha parameter estimated from Gompertz relational model
#' @param beta parameter estimated from Gompertz relational model
#' @return a vector of age-specific fertility rates from ages 10-49 by 1-year age group
#'
#' @keywords internal
#'
fx1_from_tf =
  function( TF, alpha, beta ){
    zaba = DemoToolsData::std.zaba
    fx = c()
    for( x in 10:49 ){
      Yxs = zaba[ zaba$age == x, ]$Yx_std
      Yx1s = zaba[ zaba$age == x + 1, ]$Yx_std

      fx =
        c( fx,
           TF * ( exp( - exp( - alpha - beta * Yx1s ) ) -
                    exp( - exp( - alpha - beta * Yxs ) ) ) )

    }
    fx
  }

#' Estimate completeness of VR data using cohort parity from census
#'
#' @param asfr a data.frame with the raw age-specific fertility estimates from VR data for years previous to the census
#' `year` - reference year of age-specific fertility rate estimate;
#' `age_start` - starting age of age group;
#' `age_end`   - ending age of age group;
#' `fx`        - age-specific fertility rate
#'
#' @param parity a data.frame with the information on mean children ever born from census
#' `age_start` - starting age of age group;
#' `age_end`   - ending age of age group;
#' `P`         - age-specific parity rate (mean children ever born by age group)
#'
#' @return a data.frame with the same information of `parity` plus two columns:
#' `Ei` - estimated parity equivalent accumulated from fertility information from VR data
#' `Completeness` - estimated completeness by age group (warning: assess which age groups are more reasonable to use)
#'
cohparitycomp_vr =
  function( asfr, parity ){

    # list years
    years = sort( unique( asfr$year ) )
    # number of years
    m = length( years )
    # max year
    s = max( years )
    # list ages
    ages = sort( unique( asfr$age_start ) )
    # age group diff
    n = unique( diff( ages ) )
    # numer of age groups
    nages = length( ages )

    # make sure all is ordered
    asfr = asfr[ order( asfr$year, asfr$age_start ), ]

    # cumulate fertility within each age-group
    asfr$Fx = NA
    for( t in years ){
      asfr[ asfr$year == t, ]$Fx = cumsum( asfr[ asfr$year == t, ]$fx * n )
    }

    # Ratios between adjacent Fxs
    asfr$Fx_Ratio = NA
    for( t in years ){
      temp = asfr[ asfr$year == t, ]
      asfr[ asfr$year == t, ]$Fx_Ratio = c( temp[ 1 : ( nages - 1 ), ]$Fx / temp[ 2 : nages, ]$Fx, NA )
      rm( temp )
    }

    asfr$zx = - log( - log( asfr$Fx_Ratio ) )

    zaba_std = get_zaba_params()
    c = mean( zaba_std$Phi2, na.rm = T )
    asfr =
      merge(
        asfr,
        zaba_std[ , c( 'age', 'Yxs', 'ex', 'gx' ) ],
        by.x = 'age_end',
        by.y = 'age'
      )

    # reorder
    asfr = asfr[ order( asfr$year, asfr$age_start ), ]

    # estimate alpha and beta from gompertz model
    asfr$alpha = NA
    asfr$beta  = NA
    for( t in years ){
      temp = asfr[ asfr$year == t, ]

      na_index = which( is.na( temp$zx ) )

      yfit = temp$zx - temp$ex
      yfit = yfit[ - na_index ]
      xfit = temp$gx[ - na_index ]
      temp_lm = lm( yfit ~ xfit )
      slope = as.numeric( temp_lm$coefficients[2] )
      intercept = as.numeric( temp_lm$coefficients[1] )
      beta = slope
      alpha = intercept - ( beta - 1 ) ^ 2 * c / 2

      asfr[ asfr$year == t, ]$alpha = alpha
      asfr[ asfr$year == t, ]$beta = beta
    }

    # calculate fertility level equivalent to the cumulant fertility
    asfr$TF =
      asfr$Fx /
      ( exp( - exp( - asfr$alpha - asfr$beta * asfr$Yxs ) ) )

    # get single age fertility schedules
    asfrx1 =
      do.call(
        rbind,
        lapply( years,
                function( t ){
                  TF_t = mean( asfr[ asfr$year == t & asfr$age_start %in% 20:34, ]$TF )
                  alpha_t = unique( asfr[ asfr$year == t & asfr$age_start %in% 20:34, ]$alpha )
                  beta_t  = unique( asfr[ asfr$year == t & asfr$age_start %in% 20:34, ]$beta )

                  asfr_x1 =
                    data.frame(
                      year = t,
                      age_start = 10:49,
                      fx = fx1_from_tf( TF_t, alpha_t, beta_t )
                    )
                  return( asfr_x1 )
                } )
      )

    # calculate the parity equivalent by cumulating birth cohorts
    Ei =
      sapply( 1:nages,
              function( age_index ){
                Ei = 0
                for( j in 0 : ( 5 * age_index + 3 ) ){
                  for( m in ( 5 * age_index + 9 ) : ( 5 * age_index + 13 ) ){
                    t = s - j
                    x  = m - j
                    fi = asfrx1[ asfrx1$year == t & asfrx1$age_start == x, ]$fx
                    if( length( fi ) == 0 ) fi = 0
                    Ei = Ei + fi
                  }
                }
                Ei / 5
              })

    # add to parity data and calculate vr completeness
    parity$Ei = Ei
    parity$Completeness = parity$Ei / parity$P

    return( parity )
}



