#' Cohort analysis of mean children ever born (CEB)
#'
#' Construct estimates of mean CEB in the reference period of childbearing
#' (mean age of childbearing)
#'
#' @param date_svy survey or census date in the format Y-m-d or string in the formats
#' Y, Y-m, Y-m-d
#' @param ages_w women five-year age group starting-ages vector
#' @param pop_w total women matching ages in ages_w
#' @param ceb total children ever born by women age group
#' @param mean_ceb mean children ever born by women age group
#' @param mac mean age at childbearing matching ages in age_w or unique value for all ages,
#' default = 28
#' @param location_mac retrieve MAC from a given location id or name from WPP2019
#' @param plot_ceb dummy variable to return the plot for mean CEB results
#'
#' @return data.frame with 2 vectors: reference date of mean children ever born (decimal)
#' and mean children ever born
#' @export

#' @examples
#'
#' ## Kenya 1989 census
#' date_svy <- "1989-08-25"
#' ages_w <- seq( 40, 75, 5 )
#' pop_w <- c( 350140, 280920, 230080, 173260, 158140, 111360, 82080, 54220 )
#' ceb <- c( 2532140, 2151920, 1736540, 1314140, 1143740, 816820, 560520, 371060 )
#'
#'ceb_eval( date_svy, ages_w, pop_w, ceb, mac = 28, plot_ceb = TRUE )
#'
#'# Using mac retrieved from WPP2019
#'ceb_eval( date_svy, ages_w, pop_w, ceb, location_mac = 'Kenya', plot_ceb = TRUE )
ceb_eval <- function( date_svy,
                      ages_w,
                      pop_w = NULL,
                      ceb = NULL,
                      mean_ceb = NULL,
                      mac = 28,
                      location_mac = NULL,
                      plot_ceb = TRUE ){

    year_svy <- decimal_anydate( date_svy )


    if( !is.null( location_mac ) ){

      if ( !is.numeric( location_mac ) ){

        location_code <- get_location_code( location_mac )

      } else {

        location_code <- location_mac

      }

      age_interval <- unique( diff( ages_w ) )
      mac <- fetch_mac_Wpp2019( location_code = location_code, year = year_svy )

      cat( paste0( 'WPP 2019 - Mean Age at Childbearing for ',
                   get_location_name( location_code ),
                   '\nReference year ',
                   round( year_svy ),
                   ': ',
                   mac,
                   '\n' ) )
    }

    if( length( mac ) == 1 ){
      mac <- rep( mac, length( ages_w ) )
    }

    age_interval <- unique( diff( ages_w ) )

    if( length( age_interval ) != 1 ){
      stop( 'Please select an unique age interval length.')
    }

    midgroup_ages <- ( ages_w ) + age_interval / 2

    year_ceb <- round( year_svy - ( midgroup_ages - mac ), 3 )

    if( is.null( mean_ceb ) ){

      stopifnot(
        length( ages_w ) == length( pop_w ),
        length( ages_w ) == length( ceb ),
        length( ages_w ) == length( mac )
      )

      print( 'Estimates using mean_ceb values computed from pop_w and ceb.' )

      mean_ceb <- round( ceb / pop_w , 3 )

    } else{

      stopifnot(
        length( ages_w ) == length( mean_ceb ),
        length( ages_w ) == length( mac )
      )

      print( 'Estimates using mean_ceb provided values.' )
    }

    if( plot_ceb ){

      plty_range <- c( 0, ceiling( max( mean_ceb ) ) )
      pltx_range <- c( floor( min( year_ceb ) ), ceiling( max( year_ceb ) ) )

      plot( x = year_ceb, y = mean_ceb,
            type = 'b', col = 'blue',
            xlim = pltx_range, ylim = plty_range,
            xlab = 'Reference year of CEB',
            ylab = 'Children per Women',
            main = 'Cohort analysis of mean CEB' )
      graphics::legend( 'bottomleft',
              legend = paste0( substr( date_svy, 1, 4 ), ' census/survey' ),
              lty = c( 5 ), col = c( 'blue' ), pch = c( 1 ),
              bty = 'n' )
    }

    return( cbind( year_ceb, mean_ceb ) )
  }


#' Function to retrieve mean age of childbearing from Wpp 2019 data
#'
#' @param location_code list of location codes to retrieve fertility pattern data from
#' @param year period of reference to retrieve fertility pattern data
#'
#' @return mean age at childbearing for country and year selected
#'
#' @examples
#' # mac for Honduras 2007
#' fetch_mac_Wpp2019( 340, 2007 )
#'
#'
fetch_mac_Wpp2019 <- function( location_code = NULL, year ){

  percentASFR <- load_named_data('percentASFR', "wpp2019")
  year_interv <- findInterval( x = year, vec = seq( 1950, 2020, 5 ) )

  year_sup <- seq( 1950, 2020, 5 )[ year_interv + 1 ]
  year_inf <- seq( 1950, 2020, 5 )[ year_interv ]

  # standardized fertility distribution of selected country
  asfr_std <- c( percentASFR[ percentASFR$country_code %in% location_code,
                              c( paste0( year_inf, '-', year_sup) ) ] / ( 5 * 100 ) )
  ages <- seq( 15, 45, 5 )

  mac <- round( crossprod( asfr_std, ages )/ sum( asfr_std ), 3 )

  return( mac )
}
