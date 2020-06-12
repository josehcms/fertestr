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
#' @param date_ref reference date of inquiry given in the following formats:
#' Y-m-d (4 digit year - 2 digit month - 2 digit day), Y-m (4 digit year - 2 digit month), Y (4 digit year)
#' @param location_list list of either country names or location IDs from wpp 2019
#' @param lt_family family model life table to be used for reverse survival of children and women (Chilean, Far_East_Asian, Latin,
#' General (default), South_Asian, North, South, East, West)
#'
#' @return data.frame with 2 elements: year (reference period of fertility estimation) and
#' TFR (indirect estimated total fertility rate) plus location name and ID if using wpp 2019 country data
#'
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#' ## reverse survival for 2008 Cambodia census data (Moultrie et al, 2013)
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
#' FertRevSurv( ages_c = 0:14, pop_c, ages_w = seq(10,65,5),
#' pop_w,  lx_c, lx_w, asfr_std, asfr_std_15prior, q0_5, q15_45,
#' date_ref = '2008-03-03')
#'
#' ## reverse survival for 5 selected countries in 2010 using UN General mortality profile
#'
#' countries <- c( 32, 76, 380, 508, 752 ) # Argentina, Brazil, Italy, Mozambique, Sweden
#' FertRevSurv( location_list = countries, date_ref = 2010, lt_family = 'General' )

FertRevSurv <- function( ages_c = 0:14, pop_c,
                         ages_w = seq( 10, 65, 5 ), pop_women,
                         lx_c, lx_w,
                         asfr_std         = c( 0, 0.017, 0.055, 0.057, 0.041, 0.022, 0.007, 0.002 ),
                         asfr_std_15prior = NULL,
                         q0_5 = NULL, q15_45 = NULL,
                         date_ref,
                         location_list = NULL,
                         lt_family = 'West'
                         ){

  if ( !is.null( location_list ) & !is.numeric( location_list ) ){
    location_codes <- get_location_code( location_list )
  } else {
    location_codes <- location_list
  }

  year <- decimal_anydate( date_ref )

  if( is.null( location_list ) ){
    datChildren <-
      data.frame(
        ages_c,
        pop_c
      )

    datWomen <-
      data.frame(
        ages_w,
        pop_w
      )

    fertPattern <-
      data.frame(
        age = seq( 10, 45, 5 ),
        asfr_std_ref = c( asfr_std ),
        asfr_std_15prior = c( asfr_std )
      )

    if( !is.null( asfr_std_15prior) ){
      fertPattern$asfr_std_15prior <-  c( asfr_std_15prior )
    }

    # This is the point to insert log-quadratic mortality estimation function
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

    revSurvTFR <- data.frame()

    print( paste0( 'Reverse Survival Fertility Estimation - Reference date: ',
                    substr( lubridate::date_decimal( year ), 1, 10 ) ) )

    alphaChildren <- alphaRevSurv( lx_std = lxChildren_std[ lxChildren_std$age == 5, ],
                                   qx = q0_5,
                                   type = 'child' )

    alphaWomen <- alphaRevSurv( lx_std = lxWomen_std,
                                qx =  q15_45,
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

  }

  # country estimation:
  else{

    year <- trunc(year)

    if ( ! ( year %in% 1950:2100 ) ){
      stop( 'year value out of bounds, please insert year value within 1950-2100 for using WPP 2019 data for estimation')
    }

    datChildren <-
      popWpp2019x1[ popWpp2019x1$LocID %in% location_codes & popWpp2019x1$Time == year & popWpp2019x1$AgeGrp %in% 0:14,
                    c( 'LocID','AgeGrp', 'PopTotal' ) ]
    names(datChildren) <- c( 'LocID' ,'ages_c', 'pop_c' )

    popWomen.x1  <-
      popWpp2019x1[ popWpp2019x1$LocID %in% location_codes & popWpp2019x1$Time == year & popWpp2019x1$AgeGrp %in% 10:69,
                    c( 'LocID', 'AgeGrp', 'PopFemale' ) ]
    popWomen.x1$age.x5 <-
      cut( popWomen.x1$AgeGrp, breaks = seq( 10, 70, 5 ), labels = seq( 10, 65, 5 ), right = FALSE )

    datWomen <-
      aggregate( popWomen.x1$PopFemale,
                 by = list( LocID = popWomen.x1$LocID, ages_w = popWomen.x1$age.x5 ),
                 FUN = 'sum' )
    names(datWomen) <- c( 'LocID', 'ages_w', 'pop_w' )

    revSurvTFR <- data.frame()

    print( paste0( 'Reverse Survival - WPP 2019 Country Data - Reference Year: ', year ) )
    len <- length( location_codes )
    i = 1

    for ( country in location_codes ){

      loc_name <- get_location_name( country )
      print( paste0( loc_name,' - location number ', i, ' out of ', len  ) )
      i = i + 1

      fertPattern <- data.frame()
      LT.dat <- data.frame()
      q.dat <- data.frame()
      lxChildrenF_std <- data.frame()
      lxChildrenM_std <- data.frame()
      lxChildren_std <- data.frame()
      lxWomen_std    <- data.frame()
      alphaChildren <- NULL
      alphaWomen <- NULL
      Lc <- NULL

      fertPattern <- fetch_FertPattern_Wpp2019( country, year )

      fetch_q_dat <- fetch_MortProb_Wpp2019( country,  year )

      LT.dat <- fetch_q_dat$lt_data_wpp2019
      q.dat <- fetch_q_dat$interp_data


      lxChildrenM_std <- get_mlt( lt_family, e0 = LT.dat$e0M[1], ages = seq(0,15), sex = 'Male' )
      lxChildrenF_std <- get_mlt( lt_family, e0 = LT.dat$e0F[1], ages = seq(0,15), sex = 'Female' )
      # compute lx for both sexes using female at birth factor (pg 69 Watcher - Essential Demographic Methods)
      lxChildren_std <- data.frame( age = lxChildrenM_std$age, lx_std = rep( NA, nrow( lxChildrenM_std ) ) )
      lxChildren_std$lx_std <- lxChildrenF_std$lx_std * 0.4886 + ( 1 - 0.4886 ) * lxChildrenM_std$lx_std

      lxWomen_std <- get_mlt( lt_family, e0 = LT.dat$e0F[1], ages = seq(10,65,5), sex = 'Female' )

      alphaChildren <- alphaRevSurv( lx_std = lxChildren_std[ lxChildren_std$age == 5, ],
                                     qx =  q.dat$q0_5_est,
                                     type = 'child' )

      alphaWomen <- alphaRevSurv( lx_std = lxWomen_std,
                                  qx =  q.dat$q15_45_est,
                                  type = 'women' )

      Lc <- childSurvProb( age = lxChildren_std$age,
                           lx_std = lxChildren_std$lx_std,
                           alphaChildren )

      revSurvWomen <-
        womenRevSurv( age = lxWomen_std$age,
                      lx_std = lxWomen_std$lx_std,
                      women = datWomen[ datWomen$LocID == country, ]$pop_w,
                      alphaWomen,
                      year,
                      fertPattern )

      revSurvBirths <-
        data.frame(
          year = year - seq( 0.5, 14.5, 1 ),
          births = datChildren[ datChildren$LocID == country, ]$pop_c / Lc
        )

      for( t in unique( revSurvWomen$year ) ){
        den <- sum( revSurvWomen[ revSurvWomen$year == t, ]$popWomen * revSurvWomen[ revSurvWomen$year == t, ]$asfr_std )
        num <- revSurvBirths[ revSurvBirths$year == t, ]$births

        revSurvTFR <- rbind(
          revSurvTFR,
          data.frame(
            LocID = country,
            LocName = loc_name,
            year = t,
            TFR  = num / den
          )
        )
      }

    }
  }

  return( revSurvTFR )

}

