#' Reverse Survival Fertility Estimation Function
#'
#' Reverse Survival Fertility Estimation
#'
#' @param ages1_c children single ages vector (default 0:14)
#' @param popx1_c children population in single ages (x1) matching ages1_c vector
#' @param ages5_w women five-year age group ages vector ( default seq( 10, 65, 5 ) )
#' @param popx5_w women population in five-year age group format (x5) matching ages_w vector
#' @param lx1_c children survival function vector in single age groups from 0 to 15
#' @param lx5_w women survival function vector in five-year age groups matching ages_w vector
#' @param asfr age specific fertility rates for five-year age groups from 10-45 for current period
#' of estimation
#' @param asfr_15prior standardized age specific fertility rates for five-year age groups from 10-45 for the
#' period of 15 years before the current inquiry period
#' @param q0_5 3 element vector for mortality probability between ages 0-4 for the period of estimation,
#' period 5 years prior to estimation period, and period 10 years prior to estimation period
#' @param q15_45f female adult mortality probability for the period of estimation,
#' period 5 years prior to estimation period, and period 10 years prior to estimation period
#' @param date_ref reference date of inquiry given in the following formats:
#' Y-m-d (4 digit year - 2 digit month - 2 digit day), Y-m (4 digit year - 2 digit month),
#' Y (4 digit year)
#'
#' @return data.frame with 3 elements:
#' `year` - reference period of fertility estimation in decimal format;
#' `TFR`  - estimated total fertility rate;
#' `births` - estimated total number of births
#'
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#'
#' # 1 - User input data
#' popx1_c <-  c( 281260, 261320, 268410, 286810, 278990, 293760,
#'                293490, 302060, 315970, 267190, 326980, 280260,
#'                354120, 356920, 354830 )
#'
#' popx5_w <- c(  815930, 780320, 697160, 626430, 361650, 435880,
#'                393760, 352520, 294280, 230200, 160590, NA )
#'
#' lx1_c <- c( 1.0000, 0.9320, 0.9275, 0.9228, 0.9165, 0.9125, 0.9110,
#'             0.9094, 0.9079, 0.9063, 0.9048, 0.9032, 0.9017, 0.9001,
#'             0.8986, 0.8970 )
#'
#' lx5_w <- c( 0.91381, 0.90989, 0.90492, 0.89798, 0.88893, 0.87596,
#'             0.86029, 0.84188, 0.81791, 0.78472, 0.73735, 0.67316 )
#'
#' q0_5 <-  c( 0.0683, 0.1008, 0.1189)
#'
#' q15_45f <- c( 0.1946, 0.2290, 0.2674 )
#'
#' asfr <- c( 0.0000, 0.0418,0.1535, 0.1482, 0.1118, 0.0708,
#'            0.0301, 0.0032 )
#'
#' asfr_15prior <- c( 0.0000, 0.0533, 0.1974, 0.2144, 0.1836, 0.1332,
#'                    0.0676, 0.0134 )
#'
#' FertRevSurv( ages1_c = 0:14, popx1_c = popx1_c,
#'              ages5_w = seq( 10, 65, 5 ), popx5_w = popx5_w,
#'              lx1_c = lx1_c, lx5_w = lx5_w,
#'              asfr = asfr,
#'              asfr_15prior = asfr_15prior,
#'              q0_5 = q0_5, q15_45f = q15_45f,
#'              date_ref = '2008-03-03' )
#'
#'
FertRevSurv <- function(ages1_c = 0:14,
                        popx1_c,
                        ages5_w = seq( 10, 65, 5 ),
                        popx5_w,
                        lx1_c,
                        lx5_w,
                        asfr = c( 0, 0.017, 0.055, 0.057, 0.041, 0.022, 0.007, 0.002 ),
                        asfr_15prior = NULL,
                        q0_5 = NULL,
                        q15_45f = NULL,
                        date_ref ){

  year <- decimal_anydate( date_ref )

  datChildren <-
    data.frame( ages = ages1_c, pop_c = popx1_c )

  datWomen <-
    data.frame( ages = ages5_w, pop_w = popx5_w )

  fertPattern <-
    data.frame(
      age = seq( 10, 45, 5 ),
      asfr_std_ref = asfr / ( 5 * sum( asfr ) ),
      asfr_std_15prior = asfr / sum( 5 * asfr )
      )

  if( !is.null( asfr_15prior ) ){
    fertPattern$asfr_std_15prior <- asfr_15prior / sum( 5 * asfr_15prior )
    }

  lxChildren_std <-
    data.frame(
      age = 0:15,
      lx_std = lx1_c
      )

  lxWomen_std <-
    data.frame(
      age = ages5_w,
      lx_std = lx5_w
      )

  print( paste0( 'Reverse Survival Fertility Estimation - Reference date: ',
                 substr( lubridate::date_decimal( year ), 1, 10 ) ) )

  revSurvTFR <- revSurvMain( year,
                             datWomen, lxWomen_std, q15_45f,fertPattern,
                             datChildren, lxChildren_std, q0_5 )

  return( revSurvTFR )

}


#' Reverse Survival Fertility Estimation Function for WPP2019 data
#'
#' Apply reverse survival fertility estimation to WPP2019 country data
#'
#' @param locations list of locations by name or code according to WPP 2019 location list (check fertestr::locs_avail() for list of available locations)
#' @param lt_family model life table family to be used to retrieve single age survival information
#' @param date_ref reference date of inquiry given in the following formats:
#' Y-m-d (4 digit year - 2 digit month - 2 digit day), Y-m (4 digit year - 2 digit month),
#' Y (4 digit year)
#' @param logquad if TRUE estimates lx functions from logquad models instead of using lt_family (default FALSE)
#'
#' @return data.frame with 5 elements:
#' `location_code`: WPP 2019 location code
#' `location_name`: WPP 2019 location name
#' `year` - reference period of fertility estimation in decimal format;
#' `TFR`  - estimated total fertility rate;
#' `births` - estimated total number of births
#'
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#'
#' countries <- c( 32, 76, 380, 508, 752 ) # Argentina, Brazil, Italy, Mozambique, Sweden
#' FertRevSurvWPP( locations = countries, date_ref = 2010, lt_family = 'General' )
#'
#'
FertRevSurvWPP <-
  function( locations, date_ref, lt_family = 'General', logquad = FALSE ){

    if ( !is.numeric( locations ) ){
      location_codes <- get_location_code( locations )
    } else {
      location_codes <- locations
    }

  year <- decimal_anydate( date_ref )

  if ( ! ( ceiling( year ) %in% 1950:2100 ) ){
    stop( 'year value out of bounds, please insert year value within 1950-2100 for using WPP 2019 data for estimation')
  }

  datWomen <-
    FetchPopWpp2019( locations = location_codes, year,
                     ages = 10:69, age_interval = 5, sex = 'female' )
  names( datWomen ) <- c( 'LocID', 'ages_w', 'pop_w' )

  datChildren <-
    FetchPopWpp2019( locations = location_codes, year,
                     ages = 0:14, age_interval = 1, sex = 'total' )
  names( datChildren ) <- c('LocID', 'ages_c', 'pop_c')

  revSurvTFR <- data.frame()

  print( paste0( 'Reverse Survival - WPP 2019 Country Data - Reference Date: ', date_ref ) )
  len <- length( location_codes )
  i = 1

  for ( country in location_codes ){

    loc_name <- get_location_name( country )
    print( paste0( loc_name,' - location number ', i, ' out of ', len  ) )
    i = i + 1

    fertPattern <-
      data.frame(
        age = seq( 10, 45, 5 ),
        asfr_std_ref = FetchFertilityWpp2019( country, year )$asfr_std,
        asfr_std_15prior = FetchFertilityWpp2019( country, ( year - 15 ) )$asfr_std
      )

    q0_5 <-
      q_calcWpp2019 ( location = country,
                      years = year - seq( 2.5, 12.5, 5 ),
                      sex = 'both', age_inf = 0, age_sup = 5 )$qx

    q15_45f <-
      q_calcWpp2019 ( location = country,
                      years = year - seq( 2.5, 12.5, 5 ),
                      sex = 'female', age_inf = 15, age_sup = 60 )$qx

    e0f <- FetchLifeTableWpp2019( country, year, sex = 'female' )$ex[1]
    e0t <- FetchLifeTableWpp2019( country, year, sex = 'both' )$ex[1]

    if( logquad ){

      lq_q15_45 <- q_calcWpp2019 ( location = country,
                                   years = year,
                                   sex = 'both', age_inf = 15, age_sup = 60 )$qx

      lq_q15_45f <- q_calcWpp2019 ( location = country,
                                    years = year,
                                    sex = 'female', age_inf = 15, age_sup = 60 )$qx

      lq_q0_5 <- q_calcWpp2019 ( location = country,
                                 years = year,
                                 sex = 'both', age_inf = 0, age_sup = 5 )$qx

      lq_q0_5f <- q_calcWpp2019 ( location = country,
                                  years = year,
                                  sex = 'female', age_inf = 0, age_sup = 5 )$qx

      lxChildren_std <-
        data.frame( age = 0:15,
                    lx_std = SingleAgeLogQuad( q0_5 = lq_q0_5, q15_45 = lq_q15_45,
                                               sex = 'total' )$lx[ 1:16 ] )

      lxWomen_std <-
        data.frame( age = seq( 10, 65, 5 ),
                    lx_std = SingleAgeLogQuad( q0_5 = lq_q0_5f, q15_45 = lq_q15_45f,
                                               sex = 'female' )$lx[ seq( 11, 66, 5 ) ] )

    } else{

      lxChildren_std <- find_mlt( lt_family,
                                  e0 = e0t,
                                  ages = seq( 0, 15 ), sex = 'both' )

      lxWomen_std <- find_mlt( lt_family,
                               e0 = e0f,
                               ages = seq( 10, 65, 5 ), sex = 'female' )
    }



    revSurvTFR <-
      rbind(
        revSurvTFR,
        data.frame(
          location_code = country,
          location_name = loc_name,
          revSurvMain( year,
                       datWomen[ datWomen$LocID == country, ], lxWomen_std, q15_45f,
                       fertPattern,
                       datChildren[ datChildren$LocID == country, ],
                       lxChildren_std, q0_5 ) ) )

  }

  return( revSurvTFR )
}

