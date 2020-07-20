#' Cohort analysis of mean children ever born (CEB)
#'
#' Construct estimates of mean CEB in the reference period of childbearing
#' (mean age of childbearing)
#'
#' @param date_svy survey or census date in the format Y-m-d or string in the formats
#' Y, Y-m, Y-m-d
#' @param ages_w women starting age of 5 year age group
#' @param pop_w total women matching ages in ages_w
#' @param ceb total children ever born by women age group
#' @param mac mean age at childbearing matching ages in age_w or unique value for all ages,
#' default = 28
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


ceb_eval <-
  function( date_svy,
            ages_w,
            pop_w,
            ceb,
            mac = 28,
            plot_ceb = TRUE ){

    year_svy <- decimal_anydate( date_svy )

    if( length( mac ) == 1 ){
      mac <- rep( mac, length( ages_w ) )
    }

    stopifnot(
      length( ages_w ) == length( pop_w ),
      length( ages_w )  == length( ceb ),
      length( ages_w )  == length( mac )
      )

    midgroup_ages <- ( ages_w ) + 5 / 2

    year_ceb <- round( year_svy - ( midgroup_ages - mac ), 3 )
    mean_ceb <- round( ceb / pop_w , 3 )

    if( plot_ceb ){

      plty_range <- c( 0, ceiling( max( mean_ceb ) ) )
      pltx_range <- c( floor( min( year_ceb ) ), ceiling( max( year_ceb ) ) )

      plot( x = year_ceb, y = mean_ceb,
            type = 'b', col = 'blue',
            xlim = pltx_range, ylim = plty_range,
            xlab = 'Reference year of CEB',
            ylab = 'Children per Women',
            main = 'Cohort analysis of mean CEB' )
      legend( 'bottomleft',
              legend = paste0( substr( date_svy, 1, 4 ), ' census/survey' ),
              lty = c( 5 ), col = c( 'blue' ), pch = c( 1 ),
              bty = 'n' )
    }

    return( cbind(year_ceb,mean_ceb) )
  }


