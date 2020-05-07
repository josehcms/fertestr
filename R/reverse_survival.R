require(wpp2019)
require(data.table)

getpop_wpp2019 <- function( country_code = NULL, country_name = NULL, year ){

  require( wpp2019 )
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
                                      '45-49', '50-54', '55-59', '60-64' ) &
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


pop_child <- getpop_wpp2019(country_name = 'Brazil', year = 2006)

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
          q0_5_est    = interpolate_lx( q0_5.1, q0_5.2, year1, year2, year_est),
          q15_45_est    = interpolate_lx( q15_45.1, q15_45.2, year1, year2, year_est)
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


interp_data <- get_q_wpp2019( country_name = 'Cambodia', year_survey = 2008 )$interp_data

family = 'South_Asian'
e0 = 62.61

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

lx_std_women
age = lx_std_women$age
lx_std = lx_std_women$lx_std
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

  cohsr_dat$pop0_4   <- NA
  cohsr_dat$pop5_9   <- NA
  cohsr_dat$pop10_14 <- NA
  for( x in age[1:11] ){
    cohsr_dat[cohsr_dat$age == x,]$pop0_4   <- cohsr_dat[cohsr_dat$age == x + 5,]$pop / cohsr_dat[cohsr_dat$age == x,]$Px0_4
  }
  for( x in age[1:10] ){
    cohsr_dat[cohsr_dat$age == x,]$pop5_9   <- cohsr_dat[cohsr_dat$age == x + 5,]$pop0_4 / cohsr_dat[cohsr_dat$age == x,]$Px5_9
  }
  for( x in age[1:9] ){
    cohsr_dat[cohsr_dat$age == x,]$pop10_14   <- cohsr_dat[cohsr_dat$age == x + 5,]$pop5_9 / cohsr_dat[cohsr_dat$age == x,]$Px10_14
  }



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


pop_child$Lc <- child_cohort_sr( age = lx_std_child$age, lx_std = lx_std_child$lx_std, interp_data )

pop_child$year_ref = NA
for( i in 1:15 ){
  pop_child$year_ref[i] = 2008 - i*1
}

pop_child$pop_surv <- pop_child$pop / pop_child$Lc

lx_std_w <- find_MLT( family, e0, ages = seq(10,65,5), sex = 'Female' )


female_cohort_sr <-

data("mxF")
data("mxM")
mx_dat <- as.data.table(mxM)




lt$lt

setDT(mxF)
data("popF")

Ages = mxF[name=='Brazil']$age
pop5m_pasex = mxF[name=='Brazil']$`2015-2020`

require(DemoTools)
arr <- smooth_age_5(Value = pop5m_pasex,
                    Age = Ages,
                    method = "Arriaga",
                    OAG = TRUE)
require(beers)
install.packages('demography')
require(demography)
dat = mxF[name %in% c('Brazil', 'Argentina', 'Chile', 'Germany', 'Sweden', 'Uruguay', 'Canada', 'Australia') & age < 15,
          .(name, age, mx = `2010-2015`)]
require(ggplot2)

ggplot( data = dat ) +
  geom_point( aes( x = age, y = log( mx ), color = name ), size = 3 ) +
  geom_line( aes( x = age, y = log( mx ), color = name ), size = 1.25 ) +
  theme_classic()

data("pop")

x = seq(0,100,0.01)
f1 = 1/(x)
f2 = 1/(x + 2)
f3 = 1/(x + 5)
plot(x,f,xlim = c(0,1))
lines(x,f2,col = 'red')
lines(x,f3,col='blue')

setDT(pop)

popFT
