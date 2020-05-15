

country_list = c('Argentina')
year = 2010
family = 'Latin'

revSurvWpp <- function( country_list = 'all', year, family = 'West' ){

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

  find_MLT <- function( family, e0, ages, sex ){

    if( !( family %in% unique( modelLTx1$Family ) ) ){
      stop( 'Enter a family name within the options: Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West' )
    }

    e0 <- ifelse( e0 < 20, 20, e0 )
    MLT       <- modelLTx1[ modelLTx1$Family == family & modelLTx1$Sex == sex & modelLTx1$age %in% ages, ]
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


    lxChildren_std <- find_MLT( family, e0 = LT.dat$e0M[1], ages = seq(0,15), sex = 'Male' )

    lxWomen_std <- find_MLT( family, e0 = LT.dat$e0M[1], ages = seq(10,65,5), sex = 'Female' )

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



popWpp2019x1$Location %>% unique

require(data.table)
require(dplyr)

pop_child <-
  data.frame( age = 0 : 14,
              pop =  c( 281260, 261320, 268410, 286810, 278990, 293760,
                        293490, 302060, 315970, 267190, 326980, 280260,
                        354120, 356920, 354830 )
  )

pop_fem <-
  data.frame( age = seq( 10, 65, 5 ),
              pop =  c(  815930, 780320, 697160, 626430, 361650, 435880,
                         393760, 352520, 294280, 230200, 160590, NA )
  )


year = 2010
getpop_wpp2019 <- function( country_code = NULL, country_name = NULL, year ){

  require( wpp2019 )
  data("popF")
  data("popM")
  year_interv <- findInterval( x = year, vec = seq( 1950, 2020, 5 ) )

  year_sup <- seq( 1950, 2020, 5 )[ year_interv + 1 ]
  year_inf <- seq( 1950, 2020, 5 )[ year_interv ]

  # population 0-14
  popM0_14 <- popM[ popM$age %in% c( '0-4', '5-9', '10-14' ) & popM$name == country_name,
                    c( 'age', paste0(year_inf), paste0(year_sup) ) ]

  popF0_14 <- popF[ popM$age %in% c( '0-4', '5-9', '10-14' ) & popM$name == country_name,
                    c( 'age', paste0(year_inf), paste0(year_sup) ) ]

  pop0_14 <- data.frame( age      = popF0_14$age,
                         pop_year_inf = popF0_14[,c( paste0( year_inf ) ) ] + popM0_14[,c( paste0( year_inf ) ) ],
                         pop_year_sup = popF0_14[,c( paste0( year_sup ) ) ] + popM0_14[,c( paste0( year_sup ) ) ]
  )
  pop0_14$exp_rate <- log( pop0_14$pop_year_sup / pop0_14$pop_year_inf ) / ( year_sup - year_inf )
  pop0_14$pop_year <- pop0_14$pop_year_inf * exp( ( year - year_inf ) * pop0_14$exp_rate )

  # female population 10-64
  popF10_64 <- popF[ popF$age %in% c( '10-14', '15-19', '20-24', '25-29', '30-34', '35-39', '40-44',
                                      '45-49', '50-54', '55-59', '60-64', '65-69' ) &
                       popM$name == country_name,
                     c( 'age', paste0(year_inf), paste0(year_sup) ) ]

  popF10_64 <- data.frame( age      = popF10_64$age,
                           pop_year_inf = popF10_64[,c( paste0( year_inf ) ) ],
                           pop_year_sup = popF10_64[,c( paste0( year_sup ) ) ]
  )
  popF10_64$exp_rate <- log( popF10_64$pop_year_sup / popF10_64$pop_year_inf ) / ( year_sup - year_inf )
  popF10_64$pop_year <- popF10_64$pop_year_inf * exp( ( year - year_inf ) * popF10_64$exp_rate )

  output <- list( pop_child = data.frame( age = pop0_14$age, pop = 1000*pop0_14$pop_year ),
                  pop_fem = data.frame( age = popF10_64$age, pop = 1000*popF10_64$pop_year ))

  return(output)
}


pop_child <- getpop_wpp2019(country_name = 'Brazil', year = 2015)$pop_child
pop_fem <- getpop_wpp2019(country_name = 'Brazil', year = 2015)$pop_fem

pop_child <-
  data.frame( age = 0 : 14,
              pop =  c( 2713244, 2694909, 2726957, 2790782, 2870266,
                        2931988, 2894419, 2959192, 2995612, 3188164,
                        3505216, 3352844, 3402242, 3412748, 3493711 )
              )

pop_fem <-
  data.frame( age = seq( 10, 65, 5 ),
              pop =  c(  8441348, 8432004, 8614963, 8643419, 8026854, 7121915,
                         6688796, 6141338, 5305407, 4373877, 3468085, 2616745 )
  )


std_asfr <- getfert_wpp2019(country_name = 'Brazil', year = 2015)


interpolate <- function( y1, y2, x1, x2, x ){
  y_new <- round( ( ( x - x1 ) / ( x2 - x1 ) ) * ( y2 - y1 ) + y1, 5 )
  return( y_new )
}

get_q_wpp2019 <- function( country_code = NULL, country_name = NULL, year_survey ){

  require( wpp2019 )
  require(MortalityLaws)
  year_interv <- findInterval( x = year_survey, vec = seq( 1950, 2020, 5 ) )

  year_sup <- seq( 1950, 2020, 5 )[ year_interv + 1 ]
  year_inf <- seq( 1950, 2020, 5 )[ year_interv ]
  year_half_interv <- 0.5 * ( year_sup + year_inf )

  # mortality data
  data('mxM')
  data('mxF')

  q_dat <- data.frame()

  for( year_sel in ( year_inf - seq( 0, 20, 5 ) ) ){

    ltM   <- LifeTable( x   = mxM[ mxM$name == country_name, ]$age,
                        mx  = mxM[ mxM$name == country_name, c( paste0( year_sel, '-', year_sel + 5 ) ) ],
                        lx0 = 1,
                        sex = 'male' )$lt
    ltF   <- LifeTable( x   = mxF[ mxF$name == country_name, ]$age,
                        mx  = mxF[ mxF$name == country_name, c( paste0( year_sel, '-', year_sel + 5 ) ) ],
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
          name = country_name,
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

  for ( year_est in ( year_survey - c( 2.5, 7.5, 12.5 ) ) ){

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
          name     = country_name,
          years_prior = paste0( year_survey - (year_est) - 2.5, '-', year_survey - (year_est) + 1.5),
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


interp_data <- get_q_wpp2019( country_name = 'Brazil', year_survey = year )$interp_data

family = 'Latin'
e0 = get_q_wpp2019( country_name = 'Brazil', year_survey = 2010 )$lt_data_wpp2019$e0M[1]

find_MLT <- function( family, e0, ages, sex ){

  if( !( family %in% unique( modelLTx1$Family ) ) ){
    stop( 'Enter a family name within the options: Chilean, Far_East_Asian, Latin, General, South_Asian, North, South, East, West' )
  }

  MLT       <- modelLTx1[ modelLTx1$Family == family & modelLTx1$Sex == sex & modelLTx1$age %in% ages, ]
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

lx_std_child <- find_MLT( family, e0, ages = seq(0,15), sex = 'Male' )
lx_std_women <- find_MLT( family, e0, ages = seq(10,65,5), sex = 'Female' )

interp_data$alphaC <- estimate_alpha( lx_std = lx_std_child[lx_std_child$age==5,],
                                      qx = interp_data$q0_5_est,
                                      type = 'child' )

interp_data$alphaW <- estimate_alpha( lx_std = lx_std_women,
                                      qx = interp_data$q15_45_est,
                                      type = 'women' )

child_cohort_sr <- function( age, lx_std, interp_data ){

  cohsr_dat <-
    data.frame(
      age,
      lx_std,
      Yx_std  = logit( lx = lx_std )
    )

  cohsr_dat$Yx0_4   = interp_data[ interp_data$years_prior == '0-4', ]$alphaC   + cohsr_dat$Yx_std
  cohsr_dat$Yx5_9   = interp_data[ interp_data$years_prior == '5-9', ]$alphaC   + cohsr_dat$Yx_std
  cohsr_dat$Yx10_14 = interp_data[ interp_data$years_prior == '10-14', ]$alphaC + cohsr_dat$Yx_std

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

women_sr <- function( age, lx_std, women, interp_data ){

  cohsr_dat <-
    data.frame(
      age,
      lx_std,
      Yx_std  = logit( lx = lx_std ),
      pop     = women
    )

  cohsr_dat$Px0_4   <- NA
  cohsr_dat$Px5_9   <- NA
  cohsr_dat$Px10_14 <- NA
  for( x in age[1:10]  ){
    cohsr_dat[cohsr_dat$age == x,]$Px0_4 <-
      ( 1 + exp( 2 * interp_data[ interp_data$years_prior == '0-4', ]$alphaW + cohsr_dat[cohsr_dat$age == x,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 5,] $Yx_std) ) /
      ( 1 + exp( 2 * interp_data[ interp_data$years_prior == '0-4', ]$alphaW + cohsr_dat[cohsr_dat$age == x + 5,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 10,] $Yx_std) )

    cohsr_dat[cohsr_dat$age == x,]$Px5_9 <-
      ( 1 + exp( 2 * interp_data[ interp_data$years_prior == '5-9', ]$alphaW + cohsr_dat[cohsr_dat$age == x,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 5,] $Yx_std) ) /
      ( 1 + exp( 2 * interp_data[ interp_data$years_prior == '5-9', ]$alphaW + cohsr_dat[cohsr_dat$age == x + 5,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 10,] $Yx_std) )

    cohsr_dat[cohsr_dat$age == x,]$Px10_14 <-
      ( 1 + exp( 2 * interp_data[ interp_data$years_prior == '10-14', ]$alphaW + cohsr_dat[cohsr_dat$age == x,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 5,] $Yx_std) ) /
      ( 1 + exp( 2 * interp_data[ interp_data$years_prior == '10-14', ]$alphaW + cohsr_dat[cohsr_dat$age == x + 5,]$Yx_std
                 + cohsr_dat[cohsr_dat$age == x + 10,] $Yx_std) )

    }

  cohsr_dat[cohsr_dat$age == 55,]$Px5_9 <- NA
  cohsr_dat[cohsr_dat$age %in% c( 55, 50 ),]$Px10_14 <- NA

  cohsr_dat$pop5  <- NA
  cohsr_dat$pop10 <- NA
  cohsr_dat$pop15 <- NA

  for( x in age[1:11] ){
    cohsr_dat[cohsr_dat$age == x,]$pop5   <- cohsr_dat[cohsr_dat$age == x + 5,]$pop / cohsr_dat[cohsr_dat$age == x,]$Px0_4
  }

  for( x in age[1:10] ){
    cohsr_dat[cohsr_dat$age == x,]$pop10   <- cohsr_dat[cohsr_dat$age == x + 5,]$pop5 / cohsr_dat[cohsr_dat$age == x,]$Px5_9
  }
  for( x in age[1:9] ){
    cohsr_dat[cohsr_dat$age == x,]$pop15   <- cohsr_dat[cohsr_dat$age == x + 5,]$pop10 / cohsr_dat[cohsr_dat$age == x,]$Px10_14
  }


  return( cohsr_dat[ , c( 'age', 'pop', 'pop5', 'pop10', 'pop15')] )

}

pop_fem_new <- women_sr( age = lx_std_women$age, lx_std = lx_std_women$lx_std, women = pop_fem$pop, interp_data )

year_list <- seq( year - 0.5 , year - 0.5 - 14, -1 )

pop_fem_rs <- data.frame( year = year, age = pop_fem_new$age[pop_fem_new$age %in% seq( 10, 45, 5 )],
                          pop = pop_fem_new[pop_fem_new$age %in% seq( 10, 45, 5 ),]$pop,
                          asfr_std = std_asfr$asfr_std_ref )

pop_fem_rs <- data.frame()
for( t in year_list ){

  if( t <= year & t > year - 5 ){
    pop1 = pop_fem_new$pop[pop_fem_new$age %in% seq( 10, 45, 5 )]
    pop2 = pop_fem_new$pop5[pop_fem_new$age %in% seq( 10, 45, 5 )]
    t1   = year
    t2   = year - 5
    pop_t = interpolate( pop1, pop2, t1, t2, t )
  }

  if( t <= year - 5 & t > year - 10 ){
    pop1 = pop_fem_new$pop5[pop_fem_new$age %in% seq( 10, 45, 5 )]
    pop2 = pop_fem_new$pop10[pop_fem_new$age %in% seq( 10, 45, 5 )]
    t1   = year - 5
    t2   = year - 10
    pop_t = interpolate( pop1, pop2, t1, t2, t )
  }

  if( t <= year - 10 & t > year - 15 ){
    pop1 = pop_fem_new$pop10[pop_fem_new$age %in% seq( 10, 45, 5 )]
    pop2 = pop_fem_new$pop15[pop_fem_new$age %in% seq( 10, 45, 5 )]
    t1   = year - 10
    t2   = year - 15
    pop_t = interpolate( pop1, pop2, t1, t2, t )
  }

  pop_fem_rs <-
    rbind(
      pop_fem_rs,
      data.frame(
        year = t,
        age  = pop_fem_new$age[pop_fem_new$age %in% seq( 10, 45, 5 )],
        pop  = pop_t,
        asfr_std = interpolate( std_asfr$asfr_std_ref, std_asfr$asfr_std_15prior, year - 0.5, year - 14 - 0.5, t )
      )
    )
}


pop_child$Lc <- child_cohort_sr( age = lx_std_child$age, lx_std = lx_std_child$lx_std, interp_data )
pop_child$births <-  pop_child$pop / pop_child$Lc

pop_child$year_ref = NA
for( i in 0:14 ){
  pop_child$year_ref[i+1] = year - 0.5 - i*1
}

TFT_RS <- data.frame()

for( t in unique( pop_fem_rs$year ) ){
  den <- sum( pop_fem_rs[ pop_fem_rs$year == t, ]$pop * pop_fem_rs[ pop_fem_rs$year == t, ]$asfr_std )
  num <- pop_child[ pop_child$year_ref == t, ]$births

  TFT_RS <- rbind(
    TFT_RS,
    data.frame(
      year = t,
      TFT  = num / den
      )
  )
}

require(ggplot2)
cols <- c("LINE1"="navyblue","LINE2"="tomato3")
lines <- c("LINE1"="solid","LINE2"="dashed")
ggplot() +
#  geom_line( aes( x = TFT_RS$year, y = tft, col = 'LINE1', linetype = 'LINE1' ), size = 1.5 ) +
  geom_line( aes( x = TFT_RS$year, y = TFT_RS$TFT, col = 'LINE2', linetype = 'LINE2'  ), size = 1.5 ) +
  scale_color_manual( values = cols,
                      labels = c( 'Tools for DemEst', 'fertestr' ),
                      name = '') +
  scale_linetype_manual( values = lines,
                         labels = c( 'Tools for DemEst', 'fertestr' ),
                         name = '') +
  theme_classic() +
  theme(
    legend.position = 'top'
  )
ggsave('reverse_surv.png', width = 6, height = 6)

getwd()

devtools::install_github("mpascariu/MortalityEstimate")
require(MortalityEstimate)
HMD719f <- HMD719[HMD719$sex == "female", ]

x <- c(0,1, seq(5, 110, by = 5))

a = (0:110)
wilmoth(a, LT = HMD719f)

W <- wilmoth(x, LT = HMD719f)
L2 <- wilmothLT(W, q0_5 = 0.05, e0 = 65)

lthat.logquad(coefs, a, q = 0.05, k = 0, radix = 1)
?lthat.logquad
?wilmothLT
tmp1 <- read.csv("/home/jose/Downloads/DataProgramsExamples/DataProgramsExamples/Data/coefs.logquad.HMD719.csv")
tmp2 <- array(c(as.matrix(tmp1[, 3:6])), dim=c(24, 3, 4), dimnames=list(ages.5x1, sexes, c("ax", "bx", "cx", "vx")))
coefs <- aperm (tmp2, c(1,3,2))

lthat.any2.logquad(coefs, "Female",  Q5 =0.05, QQa = 0.12 )

?wilmothLT
W$fitted.values
W$coefficients

data('popF')
popF

popWpp2019x1 <- rbind(
  fread('/home/jose/Downloads/WPP2019_PopulationBySingleAgeSex_1950-2019.csv'),
  fread('/home/jose/Downloads/WPP2019_PopulationBySingleAgeSex_2020-2100.csv')) %>%
  .[,.(LocID,Location,Time,MidPeriod,AgeGrp,AgeGrpStart,AgeGrpSpan,PopMale,PopFemale,PopTotal)]%>%
  as.data.frame()


save(popWpp2019x1,file='data/popWpp2019x1.rda')

modelLTx1
pop_child <- data.frame( age = popS[AgeGrp %in% 0:14 & Location == 'Brazil' & Time == 2015 ]$AgeGrp,
                         pop = popS[AgeGrp %in% 0:14 & Location == 'Brazil' & Time == 2015 ]$PopTotal*1000 )

popS$Variant %>% unique
require(data.table)
dat <- as.data.table(modelLTx1) %>% .[ Sex == 'Male', .(Family,age,lx1,qx1,ex1,mx1,E0)]

dat[,q0_5:=(1-lx1[age==5]/lx1[age==0]),.(Family,E0)]

ggplot(data=dat[age %in% 0:4])+
  aes(y = log(qx1), x = log(q0_5), col = as.factor(age)) +
  geom_line()+
  facet_wrap(~Family)+
  theme_bw()

aux <- dat[Family=='North' & age == 1,]

summary(lm(data=aux,log(qx1)~log(q0_5)))

plot(y=unique(dat$q0_5),x=dat[dat$age==0,]$qx1,xlim = c(0,0.65),ylim=c(0,0.65))
lines(y=unique(dat$q0_5),x=dat[dat$age==1,]$qx1,col='red')
lines(y=unique(dat$q0_5),x=dat[dat$age==2,]$qx1,col='blue')
lines(y=unique(dat$q0_5),x=dat[dat$age==3,]$qx1,col='violet')
abline(a = 0, b = 1)


install.packages('MortCast')
require(MortCast)
data(e0Mproj, package = "wpp2019")

country <- "Brazil"
# get target e0
e0m <- as.numeric(subset(e0Mproj, name == country)[-(1:2)])
# project into future
pred <- logquadj(70,70,sex = "male")
?logquadj
logquad(70, sex = "male")
# plot first projection in black and the remaining ones in heat colors
plot(pred$mx[,1], type = "l", log = "y", ylim = range(pred$mx),
     ylab = "male mx", xlab = "Age", main = country)
for(i in 2:ncol(pred$mx)) lines(pred$mx[,i],
                                col = heat.colors(20)[i])
