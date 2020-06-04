#' Converts a date of multiple formats to decimal
#'
#' Converts a date of multiple formats (year, year-month, year-month-day) to decimal
#'
#' @param date date in the format %Y-%m-%d or string in the formats %Y, %Y-%m, %Y-%m-%d
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
