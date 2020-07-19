#' Converts a date of multiple formats to decimal
#'
#' Converts a date of multiple formats (year, year-month, year-month-day) to decimal
#'
#' @param date date in the format Y-m-d or string in the formats Y, Y-m, Y-m-d
#' @return decimal of date's year
#' @export

#' @examples
#' date1 <- '2011'
#' # converts date1 for decimal value of '2011-07-31'
#' decimal_anydate(date1)
#'
#' date2 <- '2009-05'
#' # converts date2 for decimal value of '2009-05-15'
#' decimal_anydate(date2)
#'
#' date3 <- '2013-10-31'
#' # converts date3 for decimal value of '2013-10-31'
#' decimal_anydate(date3)

decimal_anydate <-
  function( date ){

    if ( is.na( date ) ){
      return(date)
    }

    datechar <- as.character( date )
    datechar <- gsub('[[:space:]]','', datechar) # remove empty spaces
    # complete format %Y-%m-%d
    if( nchar( datechar) == 10 ){

      if( ! ( as.numeric( substr( datechar, 6, 7 ) ) %in% 1:12 ) ){
        stop( 'Enter a valid month value: 01,02,03,04,05,06,07,08,09,10,11,12')
      }

      if( ! ( as.numeric( substr( datechar, 9, 10 ) ) %in% 1:31 ) ){
        stop( 'Enter a valid day value: 01,02,03,04,05,06,07,08,09,...,29,30,31')
      }

      datedec <- lubridate::decimal_date( as.Date( datechar, '%Y-%m-%d') )
    }

    # year-month format %Y-%m
    if( nchar( datechar) == 7 ){

      if( ! ( as.numeric( substr( datechar, 6, 7 ) ) %in% 1:12 ) ){
        stop( 'Enter a valid month value: 01,02,03,04,05,06,07,08,09,10,11,12')
      }

      datedec <- lubridate::decimal_date( as.Date( paste0( datechar,
                                                           '-15' ), '%Y-%m-%d') )
    }

    # year only format %Y
    if( nchar( datechar) == 4 ){

      datedec <- lubridate::decimal_date( as.Date( paste0( datechar,
                                                           '-07-31' ), '%Y-%m-%d') )
    }

    if( ! ( nchar( datechar ) %in% c( 4, 7, 10 ) ) ){
      stop('Type correct date entry: %Y-%m-%d, Y as 4
           digit value and m and d as 2 digit values.')
    }

    return(datedec)
  }

#' Available locations from WPP 2019 data
#'
#' Provides the list of available locations from WPP 2019 data
#'
#' @return data.frame with two columns name and country_code
#' @export

#' @examples
#' # available countries from wpp2019
#' locs_avail()
#'

locs_avail <- function( ){
  require(wpp2019)
  data('UNlocations')
  return(UNlocations[, c('name','country_code')])
}

#' Get WPP 2019 country codes from country names
#'
#' Provides the list of countries and respective codes available in WPP 2019 or
#' fetch the country code for given country name
#'
#' @param country_name country name or vector of country names
#' @return data.frame with two columns country_name and country_code
#' @export

#' @examples
#' # provides all country codes
#' get_country_code()
#'
#' # provides country codes for a given list of countries
#' names <- c('Brazil','Argentina','Uruguay','Paraguay')
#' get_location_code( names )
#'

get_location_code <- function( country_name ){

  locs_list <- locs_avail()
  country_code_list <- NULL
  invalid_names <- NULL
  for( name in country_name ){
    if( !( name %in% locs_list$name ) ){
      invalid_names <- c( invalid_names, name )
    }
    else{
      country_code_list <-
        c( country_code_list,
           locs_list[locs_list$name == name,]$country_code )
    }
  }

  if( !is.null( invalid_names ) ){
    stop( paste0( invalid_names,
                 ' is not a valid name among wpp 2019 location names.\n' ) )
  }

  return( country_code_list )
}


#' Get WPP 2019 location name from codes
#'
#' @param country_code country code
#' @return country_name
#' @export

#' @keywords internal

get_location_name <- function( country_code ){
  locs_list <- locs_avail()
  country_name <- as.character( locs_list[ locs_list$country_code == country_code, ]$name )
  return( country_name )
}

#' Function to get fertility pattern from Wpp 2019 data
#'
#' @param country_code list of country codes to retrieve fertility pattern data from
#' @param year period of reference to retrieve fertility pattern data
#' @return data.frame with three columns \item{age}: women ages, \item{asfr_std_ref}:
#' age-specific fertility rates for the time-period which contains the reference year and
#' \item{asfr_std_15prior}: age-specific fertility rates for the time-period
#' 15 years prior to the reference period
#'
#' @keywords internal
#'
#'
fetch_FertPattern_Wpp2019 <- function( country_code = NULL, year ){

  require( wpp2019 )
  data('percentASFR')
  year_interv <- findInterval( x = year, vec = seq( 1950, 2020, 5 ) )

  year_sup <- seq( 1950, 2020, 5 )[ year_interv + 1 ]
  year_inf <- seq( 1950, 2020, 5 )[ year_interv ]

  # standardized fertility distribution of selected country
  std_asfr <-
    data.frame(
      age = seq( 10, 45, 5 ),
      asfr_std_ref     = c( 0, percentASFR[ percentASFR$country_code %in% country_code,
                                            c( paste0( year_inf, '-', year_sup) ) ] / ( 5 * 100 ) ),
      asfr_std_15prior = c( 0, percentASFR[ percentASFR$country_code %in% country_code,
                                            c( paste0( year_inf - 15, '-', year_sup - 15) ) ] / ( 5 * 100 ) )
    )

  return(std_asfr)
}

#' Function to get child and adult-women mortality from Wpp 2019 data
#'
#' @param country_code list of country codes to retrieve mortality pattern data from
#' @param year period of reference to retrieve mortality
#'
#' @keywords internal
#'
#'
fetch_MortProb_Wpp2019 <- function( country_code = NULL, year ){

  require( wpp2019 )
  require(MortalityLaws)
  year_interv <- findInterval( x = year, vec = seq( 1950, 2020, 5 ) )

  year_sup <- seq( 1950, 2020, 5 )[ year_interv + 1 ]
  year_inf <- seq( 1950, 2020, 5 )[ year_interv ]
  year_half_interv <- 0.5 * ( year_sup + year_inf )

  # mortality data
  data('mxM')
  data('mxF')

  q_dat <- data.frame()

  for( year_sel in ( year_inf - seq( 0, 20, 5 ) ) ){

    ltM   <- LifeTable( x   = mxM[ mxM$country_code == country_code, ]$age ,
                        mx  = mxM[ mxM$country_code == country_code, c( paste0( year_sel, '-', year_sel + 5 ) ) ] ,
                        lx0 = 1,
                        sex = 'male' )$lt
    ltF   <- LifeTable( x   = mxF[ mxF$country_code == country_code, ]$age ,
                        mx  = mxF[ mxF$country_code == country_code, c( paste0( year_sel, '-', year_sel + 5 ) ) ] ,
                        lx0 = 1,
                        sex = 'female' )$lt

    e0M <- ltM[ ltM$x == 0, ]$ex
    e0F <- ltF[ ltF$x == 0, ]$ex

    # compute infant mortality for both sex using fraction of females at birth 0.4886
    lt_both <- data.frame( x = ltM$x , lx = rep( NA, nrow( ltM ) ) )
    lt_both$lx <- ltF$lx * 0.4886 + ( 1 - 0.4886 )*ltM$lx

    q0_5   <- ( lt_both[ lt_both$x == 0, ]$lx - lt_both[ lt_both$x == 5, ]$lx ) / lt_both[ lt_both$x == 0, ]$lx
    q15_45 <- ( ltF[ ltF$x == 15, ]$lx - ltF[ ltF$x == 60, ]$lx ) / ltF[ ltF$x == 15, ]$lx

    q_dat <-
      rbind(
        q_dat,
        data.frame(
          year_int = paste0( year_sel, '-', year_sel + 5 ),
          year_ref = ( 2 * year_sel + 5) / 2,
          q0_5,
          q15_45,
          e0M,
          e0F
        )
      )
  }

  interpdat_q <- data.frame( )
  years_interp <- rev( unique( q_dat$year_ref ) )

  for ( year_est in ( year - c( 2.5, 7.5, 12.5 ) ) ){

    year1 <- years_interp[ findInterval( year_est, years_interp ) ]
    year2 <- year1 + 5

    q0_5.1 <- q_dat[ q_dat$year_ref == year1, ]$q0_5
    q0_5.2 <- q_dat[ q_dat$year_ref == year2, ]$q0_5

    q15_45.1 <- q_dat[ q_dat$year_ref == year1, ]$q15_45
    q15_45.2 <- q_dat[ q_dat$year_ref == year2, ]$q15_45

    interpdat_q <-
      rbind(
        interpdat_q,
        data.frame(
          years_prior = paste0( year - (year_est) - 2.5, '-', year - (year_est) + 1.5),
          year_est    = year_est,
          q0_5_est    = interpolate( q0_5.1, q0_5.2, year1, year2, year_est),
          q15_45_est    = interpolate( q15_45.1, q15_45.2, year1, year2, year_est)
        )
      )
  }

  output <-
    list(
      lt_data_wpp2019 = q_dat,
      interp_data = interpdat_q
    )

  return( output )

}

#' 2 point interpolation
#'
#' @param y1,y2,x1,x2 Reference known points (x1,y1), (x2,y2)
#' @param x Point for which we want to estimate y based on slope given by points 1 and 2
#'
#' @return y = (x - x1)/(x2-x1) * (y2-y1) + y1
#'
#' @keywords internal
#'
#'
interpolate <- function( y1, y2, x1, x2, x ){
  y <- round( ( ( x - x1 ) / ( x2 - x1 ) ) * ( y2 - y1 ) + y1, 5 )
  return( y )
}

#' Get lx data from model life tables
#'
#' Retrieve survival data for provided model life tables
#'
#' @param lt_family Model Life Table family name (Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West)
#' @param e0 Life expectancy level for life table (lower bound = 20)
#' @param ages age selection of data (single age-interval from 0 to 130)
#' @param sex `Female` or `Male`
#' @return data.frame with selected ages `$age` and survival functions `$lx_std`
#'
#' @keywords internal
#'
#'
get_mlt <- function( lt_family, e0, ages, sex ){

  if( !( lt_family %in% unique( modelLTx1$Family ) ) ){
    stop( 'Enter a model life table family name within the options: Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West' )
  }

  if( e0 < 20 ){
    warning( 'Life expectancy lower than 20, using e0=20 as model life table reference level.')
    e0 <- 20
  }

  MLT       <- modelLTx1[ modelLTx1$Family == lt_family & modelLTx1$Sex == sex & modelLTx1$age %in% ages, ]
  e0_levels <- unique( MLT$E0 )
  e0_inf    <- e0_levels[ findInterval( e0, e0_levels ) ]
  e0_sup    <- e0_levels[ findInterval( e0, e0_levels ) + 1]
  age       <- unique( MLT$age )

  lx_inf    <- MLT[ MLT$E0 == e0_inf, ]$lx / 100000
  lx_sup    <- MLT[ MLT$E0 == e0_sup, ]$lx / 100000
  lx_interp <- interpolate( lx_inf, lx_sup, e0_inf, e0_sup, e0 )

  lx_std <-
    data.frame(
      age    = age,
      lx_std = lx_interp
    )

  return( lx_std )
}


#' Logit lx
#'
#' Cumpute logit of survival function from lx or qx
#'
#' @param qx probability of deaths
#' @param lx survival function
#'
#' @return Yx, survival function logit value
#'
#' @keywords internal
#'
#'
logit <- function( lx = NULL, qx = NULL ){

  if( is.null( qx ) & is.null( lx ) ){
    stop( 'Enter either qx or lx value' )
  }

  if( !is.null( qx ) & is.null( lx ) ){
    lx = 1 - qx
  }

  Yx = ( 1 / 2 ) * log( ( 1 - lx ) / lx )

  return( Yx )

}


#' Estimate alpha parameter for mortality
#'
#' Estimate alpha parameter for mortality of past years basedon provided qx and
#' lx_standard parameters for women and children
#'
#' @param lx_std standard values for mortality estimation
#' @param qx 3 element vector of child mortality (q0_5) or women adult mortality (q15_45) data
#' @param type 'child' if estimating values for children or 'women' for adult women
#'
#' @return alpha, vector of three coefficient elements for 0-4, 5-9 and 10-14 years prior reference date
#'
#' @keywords internal
#'
#'
alphaRevSurv <- function( lx_std, qx, type ){

  Yx_std = logit( lx = lx_std$lx_std )

  if ( type == 'child' ){
    alpha = round( logit( qx = qx ) - Yx_std, 5 )
  }

  if ( type == 'women'){
    Y15_std = logit( lx = lx_std[lx_std$age == 15,]$lx_std )
    Y60_std = logit( lx = lx_std[lx_std$age == 60,]$lx_std )

    alpha = round( ( log( qx ) / 2 ) - ( log( ( 1 - qx ) * exp( 2 * Y60_std ) - exp( 2 * Y15_std ) ) / 2 ) , 4 )
  }

  return( alpha )

}

#' Estimate children cohort survival probabilities Lc
#'
#' Estimate cohort survival probabilities for children aged 0-14 in the reference date
#'
#' @param lx_std standard survival values for children
#' @param age age vector related to standard survival values
#' @param alphaChildren alpha values estimated from alphaRevSurv function
#'
#' @return Lc vector with cohort survival probabilities
#'
#' @keywords internal
#'
#'
childSurvProb <- function( age, lx_std, alphaChildren ){

  cohsr_dat <-
    data.frame(
      age,
      lx_std,
      Yx_std  = logit( lx = lx_std )
    )

  cohsr_dat$Yx0_4   = alphaChildren[1] + cohsr_dat$Yx_std
  cohsr_dat$Yx5_9   = alphaChildren[2] + cohsr_dat$Yx_std
  cohsr_dat$Yx10_14 = alphaChildren[3] + cohsr_dat$Yx_std

  cohsr_dat$Lx0_4   = NA
  cohsr_dat$Lx5_9   = NA
  cohsr_dat$Lx10_14 = NA

  cohsr_dat$Lx0_4[1]   = 0.3 + 0.7 / ( 1 + exp( 2 * cohsr_dat[ cohsr_dat$age == 1, ]$Yx0_4 ) )
  cohsr_dat$Lx5_9[1]   = 0.3 + 0.7 / ( 1 + exp( 2 * cohsr_dat[ cohsr_dat$age == 1, ]$Yx5_9 ) )
  cohsr_dat$Lx10_14[1] = 0.3 + 0.7 / ( 1 + exp( 2 * cohsr_dat[ cohsr_dat$age == 1, ]$Yx10_14 ) )

  for ( i in 2 : ( nrow( cohsr_dat ) - 1 ) ){
    cohsr_dat$Lx0_4[ i ] = 1 / ( 1 + exp( cohsr_dat$Yx0_4[ i ] + cohsr_dat$Yx0_4[ i + 1 ] ) )
    cohsr_dat$Lx5_9[ i ] = 1 / ( 1 + exp( cohsr_dat$Yx5_9[ i ] + cohsr_dat$Yx5_9[ i + 1 ] ) )
    cohsr_dat$Lx10_14[ i ] = 1 / ( 1 + exp( cohsr_dat$Yx10_14[ i ] + cohsr_dat$Yx10_14[ i + 1 ] ) )
  }

  cohsr_dat$Px0_4   = NA
  cohsr_dat$Px5_9   = NA
  cohsr_dat$Px10_14 = NA

  cohsr_dat$Px0_4[1]   = cohsr_dat[ 1, ]$Lx0_4
  cohsr_dat$Px5_9[1]   = cohsr_dat[ 1, ]$Lx5_9
  cohsr_dat$Px10_14[1] = cohsr_dat[ 1, ]$Lx10_14


  for ( i in 2 : ( nrow( cohsr_dat ) - 1 ) ){
    cohsr_dat$Px0_4[ i ]   = cohsr_dat$Lx0_4[ i ] / cohsr_dat$Lx0_4[ i - 1 ]
    cohsr_dat$Px5_9[ i ]   = cohsr_dat$Lx5_9[ i ] / cohsr_dat$Lx5_9[ i - 1 ]
    cohsr_dat$Px10_14[ i ] = cohsr_dat$Lx10_14[ i ] / cohsr_dat$Lx10_14[ i - 1 ]
  }


  Lc = NULL
  P1 = cohsr_dat$Px0_4
  P2 = cohsr_dat$Px5_9
  P3 = cohsr_dat$Px10_14

  for( i in 1 : ( nrow( cohsr_dat ) ) ){
    t = i - 1
    j = i
    S = rep( NA, i )

    while( t >= 0 ){

      if( t %in% 0:2 ){
        S[j] = P1[i - t]
      }

      if ( t %in% 3 : 7 ){
        S[j] = P1[i - t] * ( 1 - ( t - 2 ) / 5 ) + P2[i - t] * ( ( t - 2 ) / 5 )
      }

      if ( t %in% 8 : 12 ){
        S[j] = P2[i - t] * ( 1 - ( t - 5 - 2 ) / 5 ) + P3[i - t] * ( ( t - 5 - 2 ) / 5 )
      }

      if ( t %in% 13 : 14 ){
        S[j] = P3[i - t]
      }

      t = t - 1
      j = j - 1
    }

    Lc = c( Lc, round( prod( S ), 5) )

  }

  return( Lc[1:15] )

}

#' Estimate children cohort survival probabilities Lc
#'
#' Estimate cohort survival probabilities for children aged 0-14 in the reference date
#'
#' @param age age vector related to standard survival values and women population
#' @param lx_std standard survival function for selected women
#' @param women women population for selected ages
#' @param alphaWomen alpha parameters estimated from alphaRevSurv function
#' @param year reference period of estimation
#' @param std_asfr fertility pattern - standardized age-fertility rates
#'
#' @return data.frame with parameters:
#' year: year of reverse survived population
#' AgesWomen: women reproductive period ages
#' popWomen: women population
#' asfr_std: standardized age-specific fetrtility rates for estimation period
#'
#' @keywords internal
#'
womenRevSurv <- function( age, lx_std, women, alphaWomen, year, std_asfr ){

  cohsr_dat <-
    data.frame(
      age,
      lx_std,
      Yx_std  = logit( lx = lx_std ),
      popWomen     = women
    )

  cohsr_dat$Px0_4   <- NA
  cohsr_dat$Px5_9   <- NA
  cohsr_dat$Px10_14 <- NA
  for( x in age[1:10]  ){
    cohsr_dat[cohsr_dat$age == x,]$Px0_4 <-
      ( 1 + exp( 2 * alphaWomen[1] + cohsr_dat[cohsr_dat$age == x,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 5,] $Yx_std) ) /
      ( 1 + exp( 2 * alphaWomen[1] + cohsr_dat[cohsr_dat$age == x + 5,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 10,] $Yx_std) )

    cohsr_dat[cohsr_dat$age == x,]$Px5_9 <-
      ( 1 + exp( 2 * alphaWomen[2] + cohsr_dat[cohsr_dat$age == x,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 5,] $Yx_std) ) /
      ( 1 + exp( 2 * alphaWomen[2] + cohsr_dat[cohsr_dat$age == x + 5,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 10,] $Yx_std) )

    cohsr_dat[cohsr_dat$age == x,]$Px10_14 <-
      ( 1 + exp( 2 * alphaWomen[3] + cohsr_dat[cohsr_dat$age == x,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 5,] $Yx_std) ) /
      ( 1 + exp( 2 * alphaWomen[3] + cohsr_dat[cohsr_dat$age == x + 5,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 10,] $Yx_std) )

  }

  cohsr_dat[cohsr_dat$age == 55,]$Px5_9 <- NA
  cohsr_dat[cohsr_dat$age %in% c( 55, 50 ),]$Px10_14 <- NA

  cohsr_dat$popWomen5  <- NA
  cohsr_dat$popWomen10 <- NA
  cohsr_dat$popWomen15 <- NA

  for( x in age[1:11] ){
    cohsr_dat[cohsr_dat$age == x,]$popWomen5   <-
      cohsr_dat[cohsr_dat$age == x + 5,]$popWomen / cohsr_dat[cohsr_dat$age == x,]$Px0_4
  }

  for( x in age[1:10] ){
    cohsr_dat[cohsr_dat$age == x,]$popWomen10   <-
      cohsr_dat[cohsr_dat$age == x + 5,]$popWomen5 / cohsr_dat[cohsr_dat$age == x,]$Px5_9
  }
  for( x in age[1:9] ){
    cohsr_dat[cohsr_dat$age == x,]$popWomen15   <-
      cohsr_dat[cohsr_dat$age == x + 5,]$popWomen10 / cohsr_dat[cohsr_dat$age == x,]$Px10_14
  }

  pop_fem_new <- cohsr_dat[ , c( 'age', 'popWomen', 'popWomen5', 'popWomen10', 'popWomen15' ) ]

  year_list <- seq( year - 0.5 , year - 0.5 - 14, -1 )

  pop_fem_rs <- data.frame()
  for( t in year_list ){

    if( t <= year & t > year - 5 ){
      pop1 = pop_fem_new$popWomen[pop_fem_new$age %in% seq( 10, 45, 5 )]
      pop2 = pop_fem_new$popWomen5[pop_fem_new$age %in% seq( 10, 45, 5 )]
      t1   = year
      t2   = year - 5
      pop_t = interpolate( pop1, pop2, t1, t2, t )
    }

    if( t <= year - 5 & t > year - 10 ){
      pop1 = pop_fem_new$popWomen5[pop_fem_new$age %in% seq( 10, 45, 5 )]
      pop2 = pop_fem_new$popWomen10[pop_fem_new$age %in% seq( 10, 45, 5 )]
      t1   = year - 5
      t2   = year - 10
      pop_t = interpolate( pop1, pop2, t1, t2, t )
    }

    if( t <= year - 10 & t > year - 15 ){
      pop1 = pop_fem_new$popWomen10[pop_fem_new$age %in% seq( 10, 45, 5 )]
      pop2 = pop_fem_new$popWomen15[pop_fem_new$age %in% seq( 10, 45, 5 )]
      t1   = year - 10
      t2   = year - 15
      pop_t = interpolate( pop1, pop2, t1, t2, t )
    }

    pop_fem_rs <-
      rbind(
        pop_fem_rs,
        data.frame(
          year = t,
          AgesWomen  = pop_fem_new$age[pop_fem_new$age %in% seq( 10, 45, 5 )],
          popWomen   = pop_t,
          asfr_std   = interpolate( std_asfr$asfr_std_ref, std_asfr$asfr_std_15prior, year - 0.5, year - 14 - 0.5, t )
        )
      )
  }

  return( pop_fem_rs )

}


#' Estimate single-age survival functions using Log-Quadratic Model
#'
#' Estimate single-age survival functions using Log-Quadratic Model for children 0-14 and
#' women 10-64
#'
#' @param q0_1 log-quad parameter 0q1
#' @param q0_5 log-quad parameter 0q5
#' @param q15_35 log-quad parameter 15q35
#' @param q15_45 log-quad parameter 15q45
#' @param e0 log-quad parameter e0
#' @param k log-quad parameter k
#' @param lt reference life table for modeling log quad parametes, defaul = HMD
#' @param sex sex to retrieve HMD sex and life table ('female','male','total' - default)
#'
#' @return single age life table estimated by ungrouping log-quad model estimation
#'
#' @keywords internal
#'
#

SingleAgeLogQuad <-
  function( k = NULL,
            e0 = NULL,
            q0_1 = NULL,
            q0_5 = NULL,
            q15_35 = NULL,
            q15_45 = NULL,
            lt = NULL,
            sex = 'total' ){

    require(ungroup)
    require(MortalityEstimate)
    require(MortalityLaws)

    if( is.null(lt) ){
      lt_lqmodel <- HMD719[HMD719$sex == sex, ]
    }

    # fit log-quadratic
    x <- c( 0, 1, seq( 5, 110, by = 5 ) )
    W <- wilmoth(x = x, LT = lt_lqmodel )

    # use available information to retrieve life table
    lt5 <- wilmothLT( W,
                      q0_1 = q0_1, q0_5 = q0_5,
                      q15_35 = q15_35, q15_45 = q15_45,
                      e0 = e0, k = k )

    # ungroup using ages 1 - 110
    lts_model <- pclm( x = lt5$lt$x[ 2:24 ],
                       y = lt5$lt$dx[ 2:24 ],
                       nlast = 20,
                       offset = lt5$lt$Lx[ 2:24 ] )

    # single life table
    lts <-
      LifeTable( x = 0:99,
                 mx = c( lt5$lt$mx[1], fitted( lts_model )[ 1:99 ] ),
                 lx0 = 1,
                 sex = sex )$lt

    return( lts )
  }
