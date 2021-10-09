#' Available locations from WPP 2019 data
#'
#' Provides the list of available locations from WPP 2019 data
#'
#' @return data.frame with four columns `location_name`,`location_code_iso3`,`location_code_iso2`, `location_code`
#' @export

#' @examples
#' # available countries from wpp2019
#' locs_avail()
#'

locs_avail <- function( ){
  locs <- load_named_data( 'WPP_codes', 'DemoToolsData' )
  names( locs ) <-
    c( 'location_code',
       'location_code_iso3',
       'location_code_iso2',
       'location_name' )
  return( locs )
}

#' Get WPP 2019 country codes from country names
#'
#' Provides the list of countries and respective codes available in WPP 2019 or
#' fetch the country code for given country name
#'
#' @param location_names location name or vector of location names
#'
#' @return WPP 2019 `location_code` for given location names
#'
#' @export

#' @examples
#' # provides all location codes
#' get_location_code()
#'
#' # provides location codes for a given list of countries
#' names <- c('Brazil','Argentina','Uruguay','Paraguay')
#' get_location_code( names )
#'

get_location_code <- function( location_names=NULL){

  locs_list <- locs_avail()

  if (is.null(location_names)) return(locs_list$location_code)

  location_code_list <- NULL
  invalid_names <- NULL

  # TR just to give a local binding...
  location_name <- NULL

  for( name in location_names ){
    if( !( name %in% locs_list$location_name ) ){
      invalid_names <- c( invalid_names, name )
    }
    else{
      location_code_list <-
        c( location_code_list,
      subset(locs_list, location_name == name)$location_code)
    }
  }

  if( !is.null( invalid_names ) ){
    stop( paste0( invalid_names,
                  ' is not a valid name among wpp 2019 location names.\n' ) )
  }

  return( location_code_list )
}


#' Get WPP 2019 location name from codes
#'
#' @param location_code WPP 2019 location code
#'
#' @return `location_name` of given `location_code`
#'
#' @export
#'
#' @examples
#' # Get name for Argentina (code 32)
#' get_location_name( 32 )
#'
#' # Get name for list of countries
#' get_location_name( c( 76, 32, 600, 604, 152 ))
#'

get_location_name <- function( location_code ){
  locs_list <- locs_avail()
  location_name <-
    as.character(
      locs_list[ locs_list$location_code %in% location_code, ]$location_name )
  return( location_name )
}

#' Retrieve Fertility Pattern from WPP 2019
#'
#' Retrieve age-specific fertility rates for available WPP 2019 locations
#'
#' @param locations list of locations by name or code according to WPP 2019 location list (check fertestr::locs_avail() for list of available locations)
#' @param year desired reference year to fetch fertility information (numeric)
#'
#' @return data.frame with five columns:
#' `location_code`: WPP 2019 location code
#' `location_name`: WPP 2019 location name
#' `age`: women ages in five-year age groups from 10 to 45;
#' `asfr`: interpolated age-specific fertility rates for desired year
#' `asfr_std`: proportional age-specific fertility rates for desired year (sum(asfr_std*5) = 1)
#'
#' @export
#'
#' @examples
#'
#' # Fertility pattern for Taiwan 2007
#' FetchFertilityWpp2019( locations = 158, year = 2007 )
#'
#' # Fertility pattern for Peru and Paraguay 1994
#' FetchFertilityWpp2019( locations = c( 604, 600 ), year = 1994 )
#'
#' # Fertility pattern for Italy and Australia 2002 (using names)
#' FetchFertilityWpp2019( locations = c( 'Italy', 'Australia'), year = 2002)
#'
#'
FetchFertilityWpp2019 <- function( locations = NULL, year ){

  ASFR <- load_named_data( 'WPP2019_asfr', 'DemoToolsData' )

  asfr_df <- data.frame()

  for( location_code in locations ){

    if ( !is_LocID( location_code ) ){
      location_code <- get_location_code( location_code )
    } else {
      location_code <- as.integer(location_code)
    }

    # Stop if location code is not in lt_wpp dataset
    stopifnot( location_code %in% unique( ASFR$LocID ) )

    year_range <-
      seq( min( ASFR[ ASFR$LocID == location_code, 'Year' ] ),
           max( ASFR[ ASFR$LocID == location_code, 'Year' ] ),
           5 )

    if( year >= min( year_range ) & year < max( year_range ) ){

      year_interv <-
        findInterval( x = year,
                      vec = year_range )

      year_sup <- year_range[ year_interv + 1 ]
      year_inf <- year_range[ year_interv ]

      asfr <-
        interpolate(
          y1 = ASFR[ ASFR$LocID == location_code &
                       ASFR$Year == year_inf, ]$ASFR,
          x1 = year_inf,
          y2 = ASFR[ ASFR$LocID == location_code &
                       ASFR$Year == year_sup, ]$ASFR,
          x2 = year_sup,
          x  = year
        )

    } else if( year < min( year_range ) ){

      # if year < 1953, set ASFR to 1950-55 period
      asfr <- ASFR[ ASFR$LocID == location_code &
                      ASFR$Year == min( year_range ), ]$ASFR

      warning(
        paste0( 'location_code = ',
                location_code,
                ': year is lower than min value for location_code (',
                min( year_range ),
                '). Setting asfr to the period ',
                min( year_range ),
                '-',
                ( min( year_range ) + 5 ) )
      )

    } else if( year >= max( year_range ) ){

      # if year >= 2023, set ASFR to 2018-2023 period
      asfr <- ASFR[ ASFR$LocID == location_code &
                      ASFR$Year == max( year_range ), ]$ASFR
      warning(
        paste0( 'location_code = ',
                location_code,
                ': year is greater or equal the max value for location_code (',
                max( year_range ),
                '). Setting asfr to the period ',
                ( max( year_range ) - 5 ),
                '-',
                max( year_range ) )
      )

    }

    asfr_std <-
      asfr / sum( 5 * asfr )

    asfr_df <-
      rbind( asfr_df,
             data.frame(
               location_code = location_code,
               location_name = get_location_name( location_code ),
               year_ref = year,
               age  = seq( 10, 45, 5 ),
               asfr = c( 0, asfr ),
               asfr_std = c( 0, asfr_std ) ) )
  }

  return( asfr_df )

}


#' Retrieve life table functions from WPP 2019 data
#'
#' Retrieve life tables for selected period and locations from WPP 2019 data
#'
#' @param locations list of location codes or names to retrieve LT information
#' @param year period of reference to generate estimates (numeric format)
#' @param sex sex to retrieve mortality (male, female or both)
#'
#' @return a data.frame with 15 elements:
#' `location_code`: WPP 2019 location code;
#' `location_name`: WPP 2019 location name;
#' `reference_period`: reference year of each period used for interpolating mortality rates;
#' `sex`: both, male or female;
#' `x.int`: age group interval;
#' `x`: age group starting age;
#' `mx`: age-specific mortality rates;
#' `qx`: mortality probability;
#' `ax`: person-years lived by those who die of the age group;
#' `lx`: survival function;
#' `dx`: life table synthetic cohort deaths;
#' `Lx`: person-years lived by age group;
#' `Tx`: cumulate person-years lived by age group;
#' `ex`: life expectancy
#'
#' @export
#'
#' @examples
#' # life table for both sexes Uruguay 1990
#' FetchLifeTableWpp2019( 'Uruguay', 1990 )
#'
#' # female life tables for Japan and Canada 1987
#' FetchLifeTableWpp2019( locations = c( 'Japan', 'Canada' ),
#'                        year = 1987,
#'                        sex = 'female' )
#'
#' # male life tables for Finland in '2014-10-22'
#' FetchLifeTableWpp2019( locations = 'Finland',
#'                        year = lubridate::decimal_date( as.Date( '2014-10-22' ) ),
#'                        sex  = 'male' )
#'
#'
FetchLifeTableWpp2019 <- function( locations = NULL, year, sex = 'both'){

  lt_wpp <-
    load_named_data( 'WPP2019_lt', 'DemoToolsData' )

  sex_code <-
    ifelse( tolower( sex ) == 'both',
            'b',
            ifelse( tolower( sex ) == 'female',
                    'f',
                    ifelse( tolower( sex ) == 'male' ,
                            'm',
                            NA ) ) )

  stopifnot( "Invalid sex name, please set it to 'both', 'male' or 'female'" =
               !is.na( sex_code ) )

  sex_mortlaws <-
    ifelse( sex_code == 'b',
            'total',
            tolower( sex )
            )

  lt_df <- data.frame()

  for( location_code in locations ){

    if ( !is_LocID( location_code ) ){
      location_code <- get_location_code( location_code )
    } else {
      location_code <- as.integer(location_code)
    }

    # Stop if location code is not in lt_wpp dataset
    stopifnot( location_code %in% unique( lt_wpp$LocID ) )

    year_range <-
      seq( min( lt_wpp[ lt_wpp$LocID == location_code, 'Year' ] ),
           max( lt_wpp[ lt_wpp$LocID == location_code, 'Year' ] ),
           5 )

    if( year >= min( year_range ) & year < max( year_range ) ){

      year_interv <-
        findInterval( x = year,
                      vec = year_range )

      year_sup <- year_range[ year_interv + 1 ]
      year_inf <- year_range[ year_interv ]

      mx <-
        interpolate(
          y1 = lt_wpp[ lt_wpp$LocID == location_code &
                         lt_wpp$Year == year_inf &
                         lt_wpp$Sex == sex_code, ]$mx,
          x1 = year_inf,
          y2 = lt_wpp[ lt_wpp$LocID == location_code &
                         lt_wpp$Year == year_sup &
                         lt_wpp$Sex == sex_code, ]$mx,
          x2 = year_sup,
          x  = year
        )

      age <- lt_wpp[ lt_wpp$LocID == location_code &
                       lt_wpp$Year == year_inf &
                       lt_wpp$Sex == sex_code, ]$AgeStart

    } else if( year < min( year_range ) ){

      # if year < 1953, set mx to 1950-55 period
      mx <- lt_wpp[ lt_wpp$LocID == location_code &
                      lt_wpp$Year == min( year_range ) &
                      lt_wpp$Sex == sex_code, ]$mx

      age <- lt_wpp[ lt_wpp$LocID == location_code &
                       lt_wpp$Year == min( year_range ) &
                       lt_wpp$Sex == sex_code, ]$AgeStart

      year_sup <- min( year_range ) + 5
      year_inf <- min( year_range )

      warning(
        paste0( 'location_code = ',
                location_code,
                ': year is lower than min value for location_code (',
                min( year_range ),
                '). Setting mortality rates to the period ',
                min( year_range ),
                '-',
                ( min( year_range ) + 5 ) )
      )

    } else if( year >= max( year_range ) ){

      # if year >= 2023, set mx to 2018-2023 period
      mx <- lt_wpp[ lt_wpp$LocID == location_code &
                      lt_wpp$Year == max( year_range ) &
                      lt_wpp$Sex == sex_code, ]$mx

      age <- lt_wpp[ lt_wpp$LocID == location_code &
                      lt_wpp$Year == max( year_range ) &
                      lt_wpp$Sex == sex_code, ]$AgeStart

      year_sup <- max( year_range )
      year_inf <- max( year_range ) - 5

      warning(
        paste0( 'location_code = ',
                location_code,
                ': year is greater or equal the max value for location_code (',
                max( year_range ),
                '). Setting mortality rates to the period ',
                ( max( year_range ) - 5 ),
                '-',
                max( year_range ) )
      )

    }

    lt <-
      MortalityLaws::LifeTable( x   = age,
                                mx  = mx,
                                lx0 = 1,
                                sex = sex_mortlaws )$lt

    lt_df <-
      rbind( lt_df,
             data.frame( location_code = location_code,
                         location_name = get_location_name( location_code ),
                         year = year,
                         reference_period = paste0( year_inf, '-', year_sup ),
                         sex = tolower( sex ),
                         lt ) )

  }

  return( lt_df )

}

#' Fetch population data from WPP2019
#'
#' @param locations location codes or names from WPP 2019 data
#' @param year period of reference of population data (numeric)
#' @param ages selected ages to retrieve pop data (default - 0:100, all)
#' @param age_interval how to display ages in result - single ages (1 - default)
#' or 5-year age group (5)
#' @param sex sex to retrieve information from (default - 'total', 'male', 'female')
#'
#' @return a data.frame with 3 elements:
#' `LocID`: location code
#' `ages`: age group starting age of population count
#' `pop`: population count (in thousands)
#'
#' @export
#'
#' @examples
#' # Argentina 10:69 females 2010
#' FetchPopWpp2019( 32, 2010, 10:69, 1, 'female')
#'
#' # Ecuador 15:59 males 2004 in five year age groups
#'  FetchPopWpp2019( 'Ecuador', 2004, 15:59, 5, 'male')
#'
#'  # Mexico and Panama 10:15 males 2001 in single year age groups
#'  FetchPopWpp2019( c( 'Mexico', 'Panama'), 2001, 10:15, 1, 'male')
#'
#'
#'
FetchPopWpp2019 <-
  function( locations = NULL,
            year,
            ages = 0:100,
            age_interval = 1,
            sex = 'total' ){

    popWpp2019x1 <-
      load_named_data( 'WPP2019_pop', 'DemoToolsData' )

    sex_code <-
      ifelse( tolower( sex ) == 'total',
              'b',
              ifelse( tolower( sex ) == 'female',
                      'f',
                      ifelse( tolower( sex ) == 'male' ,
                              'm',
                              NA ) ) )

    stopifnot( "Invalid sex name, please set it to 'total', 'male' or 'female'" =
                 !is.na( sex_code ) )

    popx1 <-
      data.frame()

    for( location_code in locations ){

      if ( !is_LocID( location_code ) ){
        location_code <- get_location_code( location_code )
      } else {
        location_code <- as.integer(location_code)
      }

      # Stop if location code is not in lt_wpp dataset
      stopifnot( location_code %in% unique( popWpp2019x1$LocID ) )

      year_range <-
        sort(
          unique(
            popWpp2019x1[ popWpp2019x1$LocID == location_code,]$Year
            )
        )

      if( sex_code == 'm' ){

        pop_aux_df <-
          popWpp2019x1[ popWpp2019x1$LocID == location_code &
                          popWpp2019x1$AgeStart %in% ages,
                        c( 'LocID', 'Year', 'AgeStart', 'PopMale' ) ]

        names( pop_aux_df ) <-
          c( 'LocID', 'Year', 'AgeStart', 'Pop' )

      } else if( sex_code == 'f' ){

        pop_aux_df <-
          popWpp2019x1[ popWpp2019x1$LocID == location_code &
                          popWpp2019x1$AgeStart %in% ages,
                        c( 'LocID', 'Year', 'AgeStart', 'PopFemale' ) ]

        names( pop_aux_df ) <-
          c( 'LocID', 'Year', 'AgeStart', 'Pop' )

      } else{

        pop_aux_df <-
          popWpp2019x1[ popWpp2019x1$LocID == location_code &
                          popWpp2019x1$AgeStart %in% ages,
                        c( 'LocID', 'Year', 'AgeStart', 'PopMale', 'PopFemale' ) ]

        pop_aux_df$Pop <- pop_aux_df$PopMale + pop_aux_df$PopFemale

        pop_aux_df <-
          pop_aux_df[, c( 'LocID', 'Year', 'AgeStart', 'Pop' ) ]

      }

      if( year >= min( year_range ) & year < ( max( year_range ) - 1 ) ){

        year_interv <-
          findInterval( x = year,
                        vec = year_range )

        year_sup <- year_range[ year_interv + 1 ]
        year_inf <- year_range[ year_interv ]

        pop <-
          interpolate(
            y1 = pop_aux_df[ pop_aux_df$Year == year_inf, ]$Pop,
            x1 = year_inf,
            y2 = pop_aux_df[ pop_aux_df$Year == year_sup, ]$Pop,
            x2 = year_sup,
            x  = year
          )

      } else if( year < min( year_range ) ){

        # if year < 1950.5, interpolate pop according to 1950.5-1951.5
        year_sup <- min( year_range ) + 1
        year_inf <- min( year_range )

        pop <-
          interpolate(
            y1 = pop_aux_df[ pop_aux_df$Year == year_inf, ]$Pop,
            x1 = year_inf,
            y2 = pop_aux_df[ pop_aux_df$Year == year_sup, ]$Pop,
            x2 = year_sup,
            x  = year
          )

        warning(
          paste0( 'location_code = ',
                  location_code,
                  ': year is lower than min value for location_code (',
                  min( year_range ),
                  '). Interpolating population using period ',
                  min( year_range ),
                  '-',
                  ( min( year_range ) + 1 ) )
        )

      } else if( year >= max( year_range ) ){

        # if year >= 2020.5, interpolate pop according to 2019.5-2020.5
        year_sup <- max( year_range )
        year_inf <- max( year_range ) - 1

        pop <-
          interpolate(
            y1 = pop_aux_df[ pop_aux_df$Year == year_inf, ]$Pop,
            x1 = year_inf,
            y2 = pop_aux_df[ pop_aux_df$Year == year_sup, ]$Pop,
            x2 = year_sup,
            x  = year
          )

        warning(
          paste0( 'location_code = ',
                  location_code,
                  ': year is lower than min value for location_code (',
                  min( year_range ),
                  '). Interpolating population using period ',
                  ( max( year_range ) - 1 ),
                  '-',
                  max( year_range ) )
        )

      }

      popx1 <-
        rbind(
          popx1,
          data.frame(
            LocID = location_code,
            ages  = ages,
            pop   = pop / 1000
          )
        )

    }

    if( age_interval == 5 ){

      if( ( min( ages )%%10 != 0 & min( ages )%%10 != 5 ) |
          ( max( ages )%%10 != 4 & max( ages )%%10 != 9 ) ){
        stop( 'Please insert min(ages) with final digit either 0 or 5 and max(ages) with final digit either 4 or 9.')
      } else{
        age_inf <- min( ages )
        age_sup <- ifelse( max( ages )%%10 == 4,
                           trunc( max( ages ) / 10 ) * 10,
                           trunc( max( ages ) / 10 ) * 10 + 5 )
      }

      popx1$age.x5 <-
        cut( popx1$ages,
             breaks = seq( age_inf, age_sup + 5, 5 ),
             labels = seq( age_inf, age_sup, 5 ),
             right = FALSE )

      popx5 <-
        stats::aggregate( popx1$pop,
                          by = list( ages = popx1$age.x5, LocID = popx1$LocID ),
                          FUN = 'sum' )

      names( popx5 ) <- c( 'ages', 'LocID', 'pop' )

      return( popx5[ , c( 'LocID', 'ages', 'pop' ) ] )

    } else if( age_interval == 1 ){

      return( popx1 )

    }


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

#' Check if given location code corresponds to a location ID
#'
#' @param location character string of location name, character string of \code{LocID}, or numeric \code{LocID}.
#'
#' @return TRUE if location corresponds to a location ID
#'
#' @keywords internal
#'
#' @export

is_LocID <- function(location){
  location <- as.character(location)
  gsub( "[^\\d]+", "", location, perl = TRUE ) == location
}
