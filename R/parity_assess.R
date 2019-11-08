#' Parity Assessment Function
#'
#' @param data_par A data.frame composed by 3 vectors:
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#' @param na_code A numeric value representing the code for missing parities (default = NA)
#' @param eb_cor TRUE or FALSE for correction of zeros and missing values by El-Badry method (default = FALSE)
#' @param impl_par TRUE or FALSE for correction of implausible parities
#' @param eb_graph TRUE of FALSE for El-Badry diagnose plot output (default = FALSE)
#' @param par_assess TRUE or FALSE for parity assessment of average parity results
#' @param par_assess_graph TRUE or FALSE for parity assessment plot output (default = FALSE)
#' @param age_group Character assuming values 'q' (default) for quinquennial age group or 's' for single age group input
#'
#' @return A data.frame with 2 variables:
#' ages and P for mean parities by age group
#' @export
#' @examples
#' ## Kenya 1989 data:
#' avg.parity(data_par = ken_par_1989)
#' # With El-Badry:
#' avg.parity(data_par = ken_par_1989, eb_cor = TRUE, eb_graph = TRUE)
#' # Correction for implausible parities
#' avg.parity(data_par = ken_par_1989, impl_par = TRUE)


avg.parity <-
  function(
    data_par,
    na_code          = NA,
    eb_cor           = FALSE,
    impl_par         = FALSE,
    eb_graph         = FALSE,
    par_assess       = FALSE,
    par_assess_graph = FALSE,
    age_group        = 'q'){

    # Change na_code for NA
    data_par$parity[data_par$parity==na_code] <- NA

    if ( impl_par ){
      data_par <-
        implausib_parity(data_par,
                         age_group)
    }

    if ( age_group == 's' ){
      # restrict ages to the interval 10 - 54
      data_par <-
        data_par[data_par$ages %in% 10:54,]

      # group ages in five-year age group categories
      data_par$ages <-
        as.numeric(
          paste0(
            cut(data_par$ages,
                breaks = seq(10,55,5),
                labels = seq(10,50,5),
                right  = FALSE)
          )
        )

      data_par <-
        aggregate(
          women ~ ages + parity,
          data = data_par,
          FUN = 'sum'
          )

    }

    # compute maximum unkown women parity percentage
    max_missing_rt <-
      max( 100 *
             tapply( data_par[ is.na(data_par$parity),]$women,
                     data_par[ is.na(data_par$parity),]$ages,
                     sum
                     ) /
             tapply( data_par[!is.na(data_par$parity),]$women,
                     data_par[!is.na(data_par$parity),]$ages,
                     sum
                     )
           )
    if ( eb_cor ){
      if ( max_missing_rt < 2 ){
        warning('Maximum parity missing percentage lower than 2% , verify whether El Badry correction is necessary!')
      }
      data_par <-
        el_badry(data_par,
                 eb_graph)
    }
    if ( eb_cor == FALSE & max_missing_rt >= 2 ){
      warning('Maximum unknown women parity percentage exceeds 2% , El Badry correction is recommended!')
    }

    avg_par <-
      data.frame(
        ages = unique( data_par$ages ),
        P    =
          tapply(data_par[!is.na(data_par$parity), ]$women*data_par[!is.na(data_par$parity), ]$parity,
                 data_par[!is.na(data_par$parity), ]$ages,
                 sum
          ) /
          tapply(data_par$women,
                 data_par$ages,
                 sum
          ),
        row.names = NULL
      )

    if( par_assess ){
      avg_par <-
        parity_assessment(avg_par,
                          par_assess_graph)
    }

    return(avg_par)
  }

#' El-Badry correction Function
#'
#' @param data_par A data.frame composed by 3 vectors:
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#' @param na_code A numeric value representing the code for missing parities (default = NA)
#' @param eb_graph TRUE of FALSE for El-Badry diagnose plot output (default = FALSE)
#'
#' @return The original data_par data.frame with corrected numbers of zero and missing parities
#' estimated by the El-Badry function
#' @export
#' @examples
#' ## Kenya 1989 data:
#' el_badry(data_par = ken_par_1989, eb_graph = T)


el_badry <-
  function(data_par,
           eb_graph = FALSE,
           na_code  = NA){

    # Change na_code for NA
    data_par$parity[data_par$parity==na_code] <- NA

    # Compute women totals by age group
    totals_1 <-
      tapply( data_par$women,
              data_par$ages,
              sum
      )

    # Compute percentage of unknown parities by age group
    Ui <-
      tapply( data_par[is.na(data_par$parity),]$women,
              data_par[is.na(data_par$parity),]$ages,
              sum
      ) /
      totals_1

    # Compute percentage of zero parities by age group
    Zi <-
      tapply( data_par[!is.na(data_par$parity) & data_par$parity == 0,]$women,
              data_par[!is.na(data_par$parity) & data_par$parity == 0,]$ages,
              sum
      ) /
      totals_1

    # Take the intercept of regression Ui vs Zi
    alpha <- lm(as.numeric(Ui)~as.numeric(Zi))$coefficients[1]

    data_par[!is.na(data_par$parity) & data_par$parity == 0,]$women <-
      round(
        ( Ui + Zi - alpha ) * totals_1 ,
        digits=0
      )

    data_par[is.na(data_par$parity),]$women <-
      totals_1 -
      tapply( data_par[!is.na(data_par$parity),]$women,
              data_par[!is.na(data_par$parity),]$ages,
              sum
      )

    if (eb_graph==TRUE){
      plot(
        x = as.numeric(Zi),
        y = as.numeric(Ui),
        main = 'Relation of zero and unknown parities (Zi vs Ui)',
        xlab = 'Zi (zero parity proportions)',
        ylab = 'Ui (unknown parity proportions)',
        pch  = 16,
        cex  = 1.5,
        xlim = c(0,max(as.numeric(Zi),na.rm=T)+0.1),
        ylim = c(0,max(as.numeric(Ui),na.rm=T)+0.05)
      )
      text(x = as.numeric(Zi),
           y = as.numeric(Ui),
           labels = c('15-19','20-24','25-29','30-34','35-39','40-44','45-49'),
           cex=1,
           pos=4
      )
      abline(
        lm( as.numeric(Ui) ~ as.numeric(Zi) ) ,
        col = 'red',
        cex = 1.25 ,
        lty = 3
      )
      legend(
        'topleft',
        legend = paste0( 'Ui = ' ,
                         round(lm( as.numeric(Ui) ~ as.numeric(Zi) )$coefficients[1],3) ,
                         ' + Zi * ',
                         round(lm( as.numeric(Ui) ~ as.numeric(Zi) )$coefficients[2],3)
        ),
        col = 'red',
        lty = 3,
        bty = 'n'
      )

    }

    return(data_par)
  }


#' Implausible Parities correction Function
#'
#' @param data_par A data.frame composed by 3 vectors:
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#' @param age_group Character assuming values 'q' (default) for quinquennial
#' age group or 's' for single age group input
#'
#' @return The original data_par data.frame with implausible parities set to 0
#' @export
#' @examples
#' ## Kenya 1989 data:
#' implausib_parity(data_par = ken_par_1989)

# correct implausible parities
implausib_parity <-
  function(data_par, age_group){
    if (age_group == 'q'){
      #  limit the maximum number of live births that a women may have had to one birth every 18 months from the age of 12
      data_par$women <-
        ifelse ( (data_par$parity > 2 / 3 * ( data_par$ages + 5 ) - 8) & (!is.na(data_par$parity) ), 0, data_par$women )
    }
    else if (age_group == 's'){
      #  limit the maximum number of live births that a women may have had to one birth every 18 months from the age of 12
      data_par$women <-
        ifelse ( (data_par$parity > 2 / 3 * ( data_par$ages ) - 8) & (!is.na(data_par$parity) ), 0, data_par$women )
    }
    return(data_par)
  }


#' El-Badry correction Function
#'
#' @param avg_par A data.frame composed by 2 vectors:
#' @param ages A vector of ages or starting ages of age group intervals
#' @param P A vector of average parities by women age group
#' @param eb_graph TRUE of FALSE for El-Badry diagnose plot output (default = FALSE)
#'
#' @return The original data_par data.frame with corrected numbers of zero and missing parities
#' estimated by the El-Badry function
#' @export
#' @examples
#' ## Kenya 1989 data:
#' el_badry(data_par = ken_par_1989, eb_graph = T)


parity_assessment <-
  function(avg_par,
           par_assess_graph = FALSE){

    avg_par$diff_P <-
      c(NA,diff(avg_par$P))

    # verify if avg_par is monotonically increasing
    verif_monoin <-
      prod(diff(avg_par$P) > 0) == 1

    if( !verif_monoin ){
      warning('Average parity estimates do not increase monotonically. Verify parity data!')
    }

    if ( par_assess_graph ){
      plot(
        x = as.numeric(avg_par$ages),
        y = as.numeric(avg_par$P),
        main = 'Assessment of estimated average parities by age group',
        xlab = 'Ages',
        ylab = 'Average parity (P)',
        pch  = 16,
        cex  = 1.5,
        xlim = c(15,50),
        ylim = c(0,ceiling(max(as.numeric(avg_par$P),na.rm=T)+0.5))
      )
      lines(
        x = as.numeric(avg_par$ages),
        y = as.numeric(avg_par$P),
        type = 'l',
        lty  = 5,
        cex  = 1.25
      )
    }

    return(avg_par)
  }
