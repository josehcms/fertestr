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

#' Find life table survival function lx data from model life tables
#'
#' Retrieve survival data for provided model life tables
#'
#' @param lt_family Model Life Table family name (Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West)
#' @param e0 Life expectancy level for life table (lower bound = 20)
#' @param ages age selection of data (single age-interval from 0 to 130)
#' @param sex `female` or `male` or `both`
#'
#' @return data.frame with selected ages `age` and survival functions `lx_std`
#'
#' @export
#'
#' @examples
#' find_mlt( lt_family = 'Latin', e0 = 60, ages = 0:5, sex = 'male')
#'
#'
find_mlt <- function( lt_family, e0, ages, sex ){
  modelLTx1 <- DemoToolsData::modelLTx1

  if( !( lt_family %in% unique( modelLTx1$family ) ) ){
    stop( 'Enter a model life table family name within the options: Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West' )
  }

  if( e0 < 20 ){
    warning( 'Life expectancy lower than 20, using e0=20 as model life table reference level.')
    e0 <- 20
  }

  get_mlt <-  function( lt_family, e0, ages, sex ){

    model_lt <- modelLTx1[ tolower( modelLTx1$family ) == tolower( lt_family ) &
                             modelLTx1$sex == tolower(sex) & modelLTx1$age %in% ages, ]

    e0_levels <- unique( model_lt$e0 )
    e0_inf    <- e0_levels[ findInterval( e0, e0_levels ) ]
    e0_sup    <- e0_levels[ findInterval( e0, e0_levels ) + 1]
    age       <- unique( model_lt$age )
    lx_inf    <- model_lt[ model_lt$e0 == e0_inf, ]$lx1 / 100000
    lx_sup    <- model_lt[ model_lt$e0 == e0_sup, ]$lx1 / 100000
    lx_interp <- interpolate( lx_inf, lx_sup, e0_inf, e0_sup, e0 )

    lx_std <-
      data.frame(
        age    = age,
        lx_std = lx_interp
      )

    return( lx_std )
  }

  if( tolower( sex ) == 'both' ){

    lx_m <- get_mlt( lt_family, e0, ages, sex = 'male' )
    lx_f <- get_mlt( lt_family, e0, ages, sex = 'female' )

    lx_std <-
      data.frame(
        age = unique( lx_m$age ),
        lx_std = round( lx_f$lx_std * 0.4886 + ( 1 - 0.4886 ) * lx_m$lx_std, 5 )
      )

  } else{

    lx_std <- get_mlt( lt_family, e0, ages, tolower( sex ) )

  }

  return( lx_std )

}


#' Compute qx probability of death
#'
#' Calculate qx probability of death between age interval age_inf, age_sup for WPP 2019 locations
#'
#' @param location location from WPP 2019 data
#' @param years years to retrieve estimation of qx (numeric)
#' @param age_inf lower bound age of estimation
#' @param age_sup upper_bound age of estimation
#' @param sex female, male or both
#'
#' @return data.frame with 3 columns, years of estimation, age_interval, qx
#'
#' @export
#'
#' @examples
#'
#' # q0_5 for Mexico for periods 0-4, 5-9 and 10-14 before 2010
#' q_calcWpp2019( location = 'Mexico',
#'                years = 2010 - c( 2.5, 7.5, 12.5 ),
#'                sex = 'both',
#'                age_inf = 0,
#'                age_sup = 5 )
#'
#' # q15_45 for Mexican females for periods 0-4, 5-9 and 10-14 before 2010
#' q_calcWpp2019( location = 'Mexico',
#'                years = 2010 - c( 2.5, 7.5, 12.5 ),
#'                sex = 'female',
#'                age_inf = 15,
#'                age_sup = 60 )
#'
q_calcWpp2019 <- function( location, years, sex, age_inf, age_sup ){

  if ( !is.numeric( location ) ){
    location_code <- get_location_code( location )
  } else {
    location_code <- location
  }

  q_df <- data.frame()

  for( year in years ){
    lt <- FetchLifeTableWpp2019( location_code, year, sex = 'both' )

    qx <-
      ( lt$lx[ lt$x == age_inf ] - lt$lx[ lt$x == age_sup ] ) / lt$lx[ lt$x == age_inf ]

    q_df <-
      rbind(
        q_df,
        data.frame(
          year = year,
          age_interval = paste0( age_inf, '-', age_sup ),
          qx = round( qx, 6 ) ) )
  }

  return( q_df )

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
#' @export
#'
#' @examples
#' logit( lx = 0.9 )
#' logit( qx = 0.25)
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

#' Estimate children cohort survival probabilities Lc for x5 version
#'
#' Estimate cohort survival probabilities for children aged 0-4, 5-9, 10-14 in the reference date
#'
#' @param lx_std standard survival values for children
#' @param age age vector related to standard survival values
#' @param alphaChildren alpha values estimated from alphaRevSurv function
#' @param a0_5 a0_5
#' @param a5_5 a5_5
#' @param a10_5 a10_5
#'
#' @return Lc vector with cohort survival probabilities
#'
#' @keywords internal
#'
#'
childSurvProb_x5 <-
  function( age, lx_std, alphaChildren, a0_5, a5_5, a10_5 ){

  lt_est <-
    data.frame(
      age,
      lx_std,
      Yx_std  = logit( lx = lx_std )
    )


  lt_est$Yx0_4   = alphaChildren[1] + lt_est$Yx_std
  lt_est$Yx5_9   = alphaChildren[2] + lt_est$Yx_std
  lt_est$Yx10_14 = alphaChildren[3] + lt_est$Yx_std

  lt_est$lx0_4   = 1 / ( 1 + exp( 2 * lt_est$Yx0_4 ) )
  lt_est$lx5_9   = 1 / ( 1 + exp( 2 * lt_est$Yx5_9 ) )
  lt_est$lx10_14   = 1 / ( 1 + exp( 2 * lt_est$Yx10_14 ) )

  lt_est$Lx0_4   = NA
  lt_est$Lx5_9   = NA
  lt_est$Lx10_14 = NA

  lt_est$Lx0_4[1]   = lt_est$lx0_4[2] * 5 + a0_5 * ( lt_est$lx0_4[1]-lt_est$lx0_4[2] )
  lt_est$Lx5_9[1]   = lt_est$lx5_9[2] * 5 + a0_5 *( lt_est$lx5_9[1]- lt_est$lx5_9[2] )
  lt_est$Lx10_14[1] = lt_est$lx10_14[2] * 5 + a0_5 * ( lt_est$lx10_14[1]-lt_est$lx10_14[2])

  lt_est$Lx0_4[2]   = lt_est$lx0_4[3] * 5 + a5_5 * ( lt_est$lx0_4[2]-lt_est$lx0_4[3] )
  lt_est$Lx5_9[2]   = lt_est$lx5_9[3] * 5 + a5_5 *( lt_est$lx5_9[2]- lt_est$lx5_9[3] )
  lt_est$Lx10_14[2] = lt_est$lx10_14[3] * 5 + a5_5 * ( lt_est$lx10_14[2]-lt_est$lx10_14[3])

  lt_est$Lx0_4[3]   = lt_est$lx0_4[4] * 5 + a10_5 * ( lt_est$lx0_4[3]-lt_est$lx0_4[4] )
  lt_est$Lx5_9[3]   = lt_est$lx5_9[4] * 5 + a10_5 *( lt_est$lx5_9[3]- lt_est$lx5_9[4] )
  lt_est$Lx10_14[3] = lt_est$lx10_14[4] * 5 + a10_5 * ( lt_est$lx10_14[3]-lt_est$lx10_14[4])

  lt_est$Px0_4   = NA
  lt_est$Px5_9   = NA
  lt_est$Px10_14 = NA

  lt_est$Px0_4[1]   = lt_est$Lx0_4[1] / ( 5 )
  lt_est$Px5_9[1]   = lt_est$Lx5_9[1] / ( 5 )
  lt_est$Px10_14[1] = lt_est$Lx10_14[1] / ( 5 )


  for ( i in 2 : 3 ){
    lt_est$Px0_4[ i ]   = lt_est$Lx0_4[ i ] / lt_est$Lx0_4[ i - 1 ]
    lt_est$Px5_9[ i ]   = lt_est$Lx5_9[ i ] / lt_est$Lx5_9[ i - 1 ]
    lt_est$Px10_14[ i ] = lt_est$Lx10_14[ i ] / lt_est$Lx10_14[ i - 1 ]
  }


  Lc = c( ( lt_est$Px0_4[ lt_est$age == 0 ] ),
          ( lt_est$Px5_9[ lt_est$age == 0 ] *
              lt_est$Px0_4[ lt_est$age == 5 ] ),
          ( lt_est$Px10_14[ lt_est$age == 0 ] *
              lt_est$Px5_9[ lt_est$age == 5 ] *
              lt_est$Px0_4[ lt_est$age == 10 ] )
          )

  return( Lc )

}

#' Estimate women survival probabilities
#'
#' Estimate women survival probabilities for females aged 15-49 the reference date
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
#' Estimate single-age survival functions using Log-Quadratic Model and ungroup
#'
#' @param k log-quad parameter k
#' @param e0 log-quad parameter e0
#' @param q0_1 log-quad parameter 0q1
#' @param q0_5 log-quad parameter 0q5
#' @param q15_35 log-quad parameter 15q35
#' @param q15_45 log-quad parameter 15q45
#' @param lt reference life table for modeling log quad parameters, default = HMD
#' @param sex sex to retrieve HMD sex and life table ('female','male','total' - default)
#'
#' @return single age life table estimated by ungrouping log-quad model estimation
#'
#' @export
#'
#' @examples
#' SingleAgeLogQuad( e0 = 70, q0_5 = 0.04 )

SingleAgeLogQuad <-
  function( k = NULL,
            e0 = NULL,
            q0_1 = NULL,
            q0_5 = NULL,
            q15_35 = NULL,
            q15_45 = NULL,
            lt = NULL,
            sex = 'total' ){

    if( is.null(lt) ){
      W <- DemoToolsData::hmd_lqcoeffs[[sex ]]
    } else{
      # fit log-quadratic
      x <- c( 0, 1, seq( 5, 110, by = 5 ) )
      W <- MortalityEstimate::wilmoth(x = x, LT = lt )
    }

    # use available information to retrieve life table
    lt5 <- MortalityEstimate::wilmothLT( W,
                                        q0_1 = q0_1, q0_5 = q0_5,
                                        q15_35 = q15_35, q15_45 = q15_45,
                                        e0 = e0, k = k )

    # ungroup using ages 1 - 110
    lts_model <- ungroup::pclm( x = lt5$lt$x[ 2:24 ],
                               y = lt5$lt$dx[ 2:24 ],
                               nlast = 20,
                               offset = lt5$lt$Lx[ 2:24 ] )

    # single life table
    lts <-
      MortalityLaws::LifeTable( x = 0:99,
                               mx = c( lt5$lt$mx[1], stats::fitted( lts_model )[ 1:99 ] ),
                               lx0 = 1,
                               sex = sex )$lt

    return( lts )
  }


#' Reverse Survival Estimation
#'
#' Estimate TFR levels from processed information of women, children and date
#'
#' @param year reference date of estimation in decimal format
#' @param datWoman women population data.frame
#' @param lxWomen_std female survival functions data.frame for reverse survival of females
#' @param q15_45f female adult mortality probability 3 element vector or single value
#' @param fertPattern female fertility pattern (age-specific standardized rates)
#' @param datChildren children population data.frame
#' @param lxChildren_std children survival functions data.frame for reverse survival of
#' children
#' @param q0_5 children mortality probability 3 element vector or single value
#'
#' @return estimates of TFR by year prior to reference date
#'
#' @keywords internal
#

revSurvMain <-
  function(  year,
             datWomen, lxWomen_std, q15_45f, fertPattern,
             datChildren, lxChildren_std, q0_5  ){

    if( length( q0_5 ) != 3 ){
      if( length( q0_5 ) == 1 ){
        q0_5 <- rep( q0_5, 3)
        warning( 'q0_5 unique value provided - q0_5 set to 3 element vector of same value' )
      } else{
        stop( 'Please provide a 3 element vector for q0_5 or an unique value for all 3 elements')
      }
    }

    if( length( q15_45f ) != 3 ){
      if( length( q15_45f ) == 1 ){
        q15_45f <- rep( q15_45f, 3)
        warning( 'q15_45f unique value provided - q15_45f set to 3 element vector of same value' )
      } else{
        stop( 'Please provide a 3 element vector for q15_45f or an unique value for all 3 elements')
      }
    }

    alphaChildren <- alphaRevSurv( lx_std = lxChildren_std[ lxChildren_std$age == 5, ],
                                   qx = q0_5,
                                   type = 'child' )

    alphaWomen <- alphaRevSurv( lx_std = lxWomen_std,
                                qx =  q15_45f,
                                type = 'women' )

    Lc <- childSurvProb( age = lxChildren_std$age,
                         lx_std = lxChildren_std$lx_std,
                         alphaChildren )

    revSurvWomen <-
      womenRevSurv( age = lxWomen_std$age,
                    lx_std = lxWomen_std$lx_std,
                    women = datWomen$pop_w,
                    alphaWomen,
                    year,
                    fertPattern )

    revSurvBirths <-
      data.frame(
        year = year - seq( 0.5, 14.5, 1 ),
        births = datChildren$pop_c / Lc )

    revSurvTFR <- data.frame()

    for( t in unique( revSurvBirths$year ) ){
      den <- sum( revSurvWomen[ revSurvWomen$year == t, ]$popWomen * revSurvWomen[ revSurvWomen$year == t, ]$asfr_std )
      num <- revSurvBirths[ revSurvBirths$year == t, ]$births

      revSurvTFR <- rbind(
        revSurvTFR,
        data.frame(
          year = t,
          TFR  = num / den,
          births = num
        )
      )
    }

    return( revSurvTFR )
  }


#' Reverse Survival Estimation for 5-year age groups of children
#'
#' Estimate TFR levels from processed information of women, children and date for five year periods before inquiry
#'
#' @param year reference date of estimation in decimal format
#' @param datWoman women population data.frame
#' @param lxWomen_std female survival functions data.frame for reverse survival of females
#' @param q15_45f female adult mortality probability 3 element vector or single value
#' @param fertPattern female fertility pattern (age-specific standardized rates)
#' @param datChildren children population data.frame
#'
#' @return estimates of TFR by year prior to reference date
#'
#' @keywords internal
#

revSurvMainx5 <-
  function(  year,
             datWomen, lxWomen_std, q15_45f, fertPattern,
             datChildren ){

    if( length( q15_45f ) != 3 ){
      if( length( q15_45f ) == 1 ){
        q15_45f <- rep( q15_45f, 3)
        warning( 'q15_45f unique value provided - q15_45f set to 3 element vector of same value' )
      } else{
        stop( 'Please provide a 3 element vector for q15_45f or an unique value for all 3 elements')
      }
    }

    alphaWomen <- alphaRevSurv( lx_std = lxWomen_std,
                                qx =  q15_45f,
                                type = 'women' )

    revSurvWomen <-
      womenRevSurv( age = lxWomen_std$age,
                    lx_std = lxWomen_std$lx_std,
                    women = datWomen$pop_w,
                    alphaWomen,
                    year,
                    fertPattern )

    revSurvBirths <-
      data.frame(
        year = datChildren$year,
        births = datChildren$B )

    revSurvTFR <- data.frame()

    for( t in unique( revSurvBirths$year ) ){
      den <- sum( revSurvWomen[ revSurvWomen$year == t, ]$popWomen * revSurvWomen[ revSurvWomen$year == t, ]$asfr_std )
      num <- revSurvBirths[ revSurvBirths$year == t, ]$births

      revSurvTFR <- rbind(
        revSurvTFR,
        data.frame(
          year = t,
          TFR  = num / den,
          births = num
        )
      )
    }

    return( revSurvTFR )
  }

#' Reverse Survival Estimation for 5-year age group of children (BETA)
#'
#' Estimate TFR levels from processed information of women, children and date
#'
#' @param year reference date of estimation in decimal format
#' @param datWoman women population data.frame
#' @param lxWomen_std female survival functions data.frame for reverse survival of females
#' @param q15_45f female adult mortality probability 3 element vector or single value
#' @param fertPattern female fertility pattern (age-specific standardized rates)
#' @param datChildren children population data.frame
#' @param lxChildren_std children survival functions data.frame for reverse survival of
#' children
#' @param q0_5 children mortality probability 3 element vector or single value
#' @param a0_5 a0_5
#' @param a5_5 a5_5
#' @param a10_5 a10_5
#'
#' @return estimates of TFR by year prior to reference date
#'
#' @keywords internal
#

revSurvMainx5_Beta <-
  function(  year,
             datWomen, lxWomen_std, q15_45f, fertPattern,
             datChildren, lxChildren_std, q0_5, a0_5, a5_5, a10_5  ){

    if( length( q0_5 ) != 3 ){
      if( length( q0_5 ) == 1 ){
        q0_5 <- rep( q0_5, 3)
        warning( 'q0_5 unique value provided - q0_5 set to 3 element vector of same value' )
      } else{
        stop( 'Please provide a 3 element vector for q0_5 or an unique value for all 3 elements')
      }
    }

    if( length( q15_45f ) != 3 ){
      if( length( q15_45f ) == 1 ){
        q15_45f <- rep( q15_45f, 3)
        warning( 'q15_45f unique value provided - q15_45f set to 3 element vector of same value' )
      } else{
        stop( 'Please provide a 3 element vector for q15_45f or an unique value for all 3 elements')
      }
    }

    alphaChildren <- alphaRevSurv( lx_std = lxChildren_std[ lxChildren_std$age == 5, ],
                                   qx = q0_5,
                                   type = 'child' )

    alphaWomen <- alphaRevSurv( lx_std = lxWomen_std,
                                qx =  q15_45f,
                                type = 'women' )

    Lc <- childSurvProb_x5( age = lxChildren_std$age,
                            lx_std = lxChildren_std$lx_std,
                            alphaChildren,
                            a0_5, a5_5, a10_5 )

    revSurvWomen <-
      womenRevSurv( age = lxWomen_std$age,
                    lx_std = lxWomen_std$lx_std,
                    women = datWomen$pop_w,
                    alphaWomen,
                    year,
                    fertPattern )

    revSurvBirths <-
      data.frame(
        year = year - c( 2.5, 7.5, 12.5 ),
        births = ( datChildren$pop_c / Lc ) / 5 )

    revSurvTFR <- data.frame()

    for( t in unique( revSurvBirths$year ) ){
      den <- sum( revSurvWomen[ revSurvWomen$year == t, ]$popWomen * revSurvWomen[ revSurvWomen$year == t, ]$asfr_std )
      num <- revSurvBirths[ revSurvBirths$year == t, ]$births

      revSurvTFR <- rbind(
        revSurvTFR,
        data.frame(
          year = t,
          TFR  = num / den,
          births = num
        )
      )
    }

    return( revSurvTFR )
  }
