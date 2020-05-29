#' Reverse Survival Fertility Estimation for World Population Prospects 2019
#'
#' Reverse Survival Fertility Estimation for countries listed in the World Population Prospects 2019
#'
#' @param country_list list of country codes to perform reverse survival fertility estimation
#' @param year base year for retrospective fertility estimation (1950-2100)
#' @param lt_family family model life table to be used for reverse survival of children and women (Chilean, Far_East_Asian, Latin,
#' General (default), South_Asian, North, South, East, West)
#'
#' @return data.frame with 3 elements: Country (country code), year (reference period of fertility estimation) and
#' TFR (indirect estimated total fertility rate)
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#' ## reverse survival for 5 selected countries in 2010 using UN General mortality profile
#'
#' countries <- c( 32, 76, 380, 508, 752 ) # Argentina, Brazil, Italy, Mozambique, Sweden
#' revSurvWpp( country_list = countries, year = 2010, lt_family = 'General' )
#'
revSurvWpp <- function( country_list = 'all', year, lt_family = 'West' ){

  # if ( country_list == 'all' ){
  #   country_list = unique( popWpp2019x1$Location )
  # }

  if ( ! ( year %in% 1950:2100 ) ){
    stop( 'year value out of bounds, please insert year value within 1950-2100')
  }

  # 1. Set data inputs
  # 1.1 child 0-14 data
  datChildren <-
    popWpp2019x1[ popWpp2019x1$LocID %in% country_list & popWpp2019x1$Time == year & popWpp2019x1$AgeGrp %in% 0:14,
                  c( 'AgeGrp', 'LocID', 'PopTotal' ) ]
  names(datChildren) <- c( 'AgesChildren', 'Country', 'popChildren' )

  # 1.2 adult women data
  popWomen.x1  <-
    popWpp2019x1[ popWpp2019x1$LocID %in% country_list & popWpp2019x1$Time == year & popWpp2019x1$AgeGrp %in% 10:69,
                  c( 'LocID', 'AgeGrp', 'PopFemale' ) ]
  popWomen.x1$age.x5 <-
    cut( popWomen.x1$AgeGrp, breaks = seq( 10, 70, 5 ), labels = seq( 10, 65, 5 ), right = FALSE )

  datWomen <-
    aggregate( popWomen.x1$PopFemale,
               by = list( popWomen.x1$age.x5, popWomen.x1$LocID ),
               FUN = 'sum' )
  names(datWomen) <- c( 'AgesWomen', 'Country', 'popWomen' )

  # 1.3 Fertility pattern data and mortality data

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

  interpolate <- function( y1, y2, x1, x2, x ){
    y <- round( ( ( x - x1 ) / ( x2 - x1 ) ) * ( y2 - y1 ) + y1, 5 )
    return( y )
  }

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
      q0_5   <- ( ltM[ ltM$x == 0, ]$lx - ltM[ ltM$x == 5, ]$lx ) / ltM[ ltM$x == 0, ]$lx
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

  find_MLT <- function( lt_family, e0, ages, sex ){

    if( !( lt_family %in% unique( modelLTx1$Family ) ) ){
      stop( 'Enter a model life table family name within the options: Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West' )
    }

    e0 <- ifelse( e0 < 20, 20, e0 )
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

  estimate_alpha <- function( lx_std, qx, type ){

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

  estimate_Lc <- function( age, lx_std, alphaChildren ){

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

  women_Surv <- function( age, lx_std, women, alphaWomen, year, std_asfr ){

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

  revSurvTFR <- data.frame()

  print( paste0( 'Reverse Survival - WPP2019 Country Data - Reference Year: ', year ) )
  len <- length( country_list )
  i = 1
  for ( country in country_list ){

    print( paste0( country,' - Country number ', i, ' out of ', len  ) )
    i = i + 1

    fertPattern <- data.frame()
    LT.dat <- data.frame()
    q.dat <- data.frame()
    lxChildren_std <- data.frame()
    lxWomen_std    <- data.frame()
    alphaChildren <- NULL
    alphaWomen <- NULL
    Lc <- NULL

    fertPattern <- fetch_FertPattern_Wpp2019( country, year )

    fetch_q_dat <- fetch_MortProb_Wpp2019( country,  year )

    LT.dat <- fetch_q_dat$lt_data_wpp2019
    q.dat <- fetch_q_dat$interp_data


    lxChildren_std <- find_MLT( lt_family, e0 = LT.dat$e0M[1], ages = seq(0,15), sex = 'Male' )

    lxWomen_std <- find_MLT( lt_family, e0 = LT.dat$e0M[1], ages = seq(10,65,5), sex = 'Female' )

    alphaChildren <- estimate_alpha(  lx_std = lxChildren_std[ lxChildren_std$age == 5, ],
                                      qx =  q.dat$q0_5_est,
                                      type = 'child' )

    alphaWomen <- estimate_alpha( lx_std = lxWomen_std,
                                  qx =  q.dat$q15_45_est,
                                  type = 'women' )

    Lc <- estimate_Lc( age = lxChildren_std$age,
                       lx_std = lxChildren_std$lx_std,
                       alphaChildren )


    revSurvWomen <-
        women_Surv( age = lxWomen_std$age,
                    lx_std = lxWomen_std$lx_std,
                    women = datWomen[ datWomen$Country == country, ]$popWomen,
                    alphaWomen,
                    year,
                    fertPattern )

    revSurvBirths <-
      data.frame(
        year = year - seq( 0.5, 14.5, 1 ),
        births = datChildren[ datChildren$Country == country, ]$popChildren / Lc
      )

    for( t in unique( revSurvWomen$year ) ){
      den <- sum( revSurvWomen[ revSurvWomen$year == t, ]$popWomen * revSurvWomen[ revSurvWomen$year == t, ]$asfr_std )
      num <- revSurvBirths[ revSurvBirths$year == t, ]$births

      revSurvTFR <- rbind(
        revSurvTFR,
        data.frame(
          Country = country,
          year = t,
          TFR  = num / den
        )
      )
    }

  }

  return( revSurvTFR )

}



#' Reverse Survival Fertility Estimation Function
#'
#' Reverse Survival Fertility Estimation
#'
#' @param ages_c children ages (default 0:14)
#' @param pop_c children population matching ages_c vector
#' @param lx_c children survival function vector matching ages_c vector
#' @param ages_w women ages (default c( 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65 ))
#' @param pop_women women population matching ages_w vector
#' @param lx_w women survival function matching ages_w vector
#' @param asfr_std standardized age specific fertility rates for five-year age groups from 10-45 for current period
#' of estimation (default set to world 2015-2020 pattern - c( 0, 0.017, 0.055, 0.057, 0.041, 0.022, 0.007, 0.002 ))
#' @param asfr_std_15_prior standardized age specific fertility rates for five-year age groups from 10-45 for the
#' period of 15 years before the current inquiry period
#' @param q0_5 3 element vector for mortality probability between ages 0-4 for the period of estimation,
#' period 5 years prior to estimation period, and period 10 years prior to estimation period
#' @param q15_45 female adult mortality probability for the period of estimation,
#' period 5 years prior to estimation period, and period 10 years prior to estimation period
#' @param year reference year of inquiry
#'
#' @return data.frame with 2 elements: year (reference period of fertility estimation) and
#' TFR (indirect estimated total fertility rate)
#'
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#' ## reverse survival for 2008 Cambodia data (Moultrie et al, 2013)
#' pop_c <-  c( 281260, 261320, 268410, 286810, 278990, 293760, 293490, 302060, 315970, 267190, 326980, 280260, 354120, 356920, 354830 )
#' pop_w <- c(  815930, 780320, 697160, 626430, 361650, 435880, 393760, 352520, 294280, 230200, 160590, NA )
#' lx_c <- c( 1.0000, 0.9320, 0.9275, 0.9228, 0.9165, 0.9125, 0.9110, 0.9094, 0.9079, 0.9063, 0.9048, 0.9032, 0.9017, 0.9001, 0.8986, 0.8970 )
#' lx_w <- c( 0.91381, 0.90989, 0.90492, 0.89798, 0.88893, 0.87596, 0.86029, 0.84188, 0.81791, 0.78472, 0.73735, 0.67316 )
#' q0_5 <-  c( 0.0683, 0.1008, 0.1189)
#' q15_45 <- c( 0.1946, 0.2290, 0.2674)
#' asfr <- c( 0.0000, 0.0418,0.1535, 0.1482, 0.1118, 0.0708, 0.0301, 0.0032 )
#' asfr_std <- asfr/(5 * sum(asfr) )
#' asfr_15prior <- c( 0.0000, 0.0533, 0.1974, 0.2144, 0.1836, 0.1332, 0.0676, 0.0134 )
#' asfr_std_15prior <- asfr_15prior/(5 * sum(asfr_15prior) )
#'
#' FertRevSurv( ages_c = 0:14, pop_c, ages_w = seq(10,65,5), pop_w,  lx_c, lx_w, asfr_std, asfr_std_15prior, q0_5, q15_45, year = 2008)
#'
FertRevSurv <- function( ages_c = 0:14, pop_c,
                         ages_w = seq( 10, 65, 5 ), pop_women,
                         lx_c, lx_w,
                         asfr_std         = c( 0, 0.017, 0.055, 0.057, 0.041, 0.022, 0.007, 0.002 ),
                         asfr_std_15prior = NULL,
                         q0_5 = NULL, q15_45 = NULL,
                         year
                         ){

  # if asfr_std_15prior is null - use asfr_std as unique fertility pattern
  # inputs: 3 elements fpr qx5 and 3 for qx15

  # 1. Set data inputs
  # 1.1 child 0-14 data
  datChildren <-
    data.frame(
      ages_c,
      pop_c
    )

  # 1.2 adult women data
  datWomen <-
    data.frame(
      ages_w,
      pop_w
    )
  # 1.3 Fertility pattern data and mortality data

  fertPattern <-
    data.frame(
      age = seq( 10, 45, 5 ),
      asfr_std_ref = c( asfr_std ),
      asfr_std_15prior = c( asfr_std )
    )

  if( !is.null( asfr_std_15prior) ){
    fertPattern$asfr_std_15prior <-  c( asfr_std_15prior )
  }

  lxChildren_std <-
    data.frame(
      age = c( ages_c, 15 ),
      lx_std = lx_c
    )

  lxWomen_std <-
    data.frame(
      age = ages_w,
      lx_std = lx_w
    )


  interpolate <- function( y1, y2, x1, x2, x ){
    y <- round( ( ( x - x1 ) / ( x2 - x1 ) ) * ( y2 - y1 ) + y1, 5 )
    return( y )
  }

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

  estimate_alpha <- function( lx_std, qx, type ){

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

  estimate_Lc <- function( age, lx_std, alphaChildren ){

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

  women_Surv <- function( age, lx_std, women, alphaWomen, year, std_asfr ){

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

  revSurvTFR <- data.frame()

  print( paste0( 'Reverse Survival Fertility Estimation - Reference Year: ', year ) )

  alphaChildren <- estimate_alpha(  lx_std = lxChildren_std[ lxChildren_std$age == 5, ],
                                    qx = q0_5,
                                    type = 'child' )

  alphaWomen <- estimate_alpha( lx_std = lxWomen_std,
                                qx =  q15_45,
                                type = 'women' )

  Lc <- estimate_Lc( age = lxChildren_std$age,
                     lx_std = lxChildren_std$lx_std,
                     alphaChildren )

  revSurvWomen <-
    women_Surv( age = lxWomen_std$age,
                lx_std = lxWomen_std$lx_std,
                women = datWomen$pop_w,
                alphaWomen,
                year,
                fertPattern )

  revSurvBirths <-
    data.frame(
      year = year - seq( 0.5, 14.5, 1 ),
      births = datChildren$pop_c / Lc
      )

  for( t in unique( revSurvWomen$year ) ){
    den <- sum( revSurvWomen[ revSurvWomen$year == t, ]$popWomen * revSurvWomen[ revSurvWomen$year == t, ]$asfr_std )
    num <- revSurvBirths[ revSurvBirths$year == t, ]$births

    revSurvTFR <- rbind(
      revSurvTFR,
      data.frame(
        year = t,
        TFR  = num / den
        )
      )
    }

  return( revSurvTFR )

}

