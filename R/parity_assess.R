#' Parity Assessment Function
#'
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#'
#' @return A data.frame with 2 variables:
#' ages and P for mean parities by age group
#' @export
#' @examples
#' ## Malawi 2008 Census data:
#' ages_ma = c(15, 20, 25, 30, 35, 40, 45)
#' asfr_ma = c(0.111, 0.245, 0.230, 0.195, 0.147, 0.072, 0.032)
#' P_ma    = c(0.283, 1.532, 2.849, 4.185, 5.214, 6.034, 6.453)
#' brass_pf(P = P_ma, asfr = asfr_ma)
#'
#'

avg.parity <-
  function(
    data_par,
    na_code          = NA,
    el_badry         = FALSE,
    implausib_parity = FALSE,
    eb_graph         = FALSE,
    age_group        = 'q'){



    if ( implausib_parity == TRUE ){

      implausib_parity_function <-
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

      data_par <-
        implausib_parity_function(data_par)

    }

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

    if (el_badry==TRUE){

      if (max_missing_rt < 2){
        warning('Maximum parity missing percentage lower than 2% , verify whether El Badry correction is necessary!')
      }

      el_badry_function <-
        function(data_par,
                 eb_graph=FALSE){

          totals_1 <-
            tapply( data_par$women,
                    data_par$ages,
                    sum
                    )

          Di <-
            tapply( data_par[is.na(data_par$parity),]$women,
                    data_par[is.na(data_par$parity),]$ages,
                    sum
                    ) /
            totals_1

          Zi <-
            tapply( data_par[!is.na(data_par$parity) & data_par$parity == 0,]$women,
                    data_par[!is.na(data_par$parity) & data_par$parity == 0,]$ages,
                    sum
                    ) /
            totals_1

          alpha <- lm(as.numeric(Di)~as.numeric(Zi))$coefficients[1]

          data_par[!is.na(data_par$parity) & data_par$parity == 0,]$women <-
            round(
              ( Di + Zi - alpha ) * totals_1 ,
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
              y = as.numeric(Di),
              main = 'Relation of zero and unknown parities (Zi vs Ui)',
              xlab = 'Zi (zero parity proportions)',
              ylab = 'Ui (unknown parity proportions)',
              pch  = 16,
              cex  = 1.5,
              xlim = c(0,max(as.numeric(Zi),na.rm=T)+0.1),
              ylim = c(0,max(as.numeric(Di),na.rm=T)+0.05)
            )
            text(x = as.numeric(Zi),
                 y = as.numeric(Di),
                 labels = c('15-19','20-24','25-29','30-34','35-39','40-44','45-49'),
                 cex=1,
                 pos=4
            )
            abline(
              lm( as.numeric(Di) ~ as.numeric(Zi) ) ,
              col = 'red',
              cex = 1.25 ,
              lty = 3
            )
            legend(
              'topleft',
              legend = paste0( 'Ui = ' ,
                               round(lm( as.numeric(Di) ~ as.numeric(Zi) )$coefficients[1],3) ,
                               ' + Zi * ',
                               round(lm( as.numeric(Di) ~ as.numeric(Zi) )$coefficients[2],3)
              ),
              col = 'red',
              lty = 3,
              bty = 'n'
            )

          }
          return(data_par)
        }

      data_par <-
        el_badry_function(data_par,
                          eb_graph)


    }

    if (el_badry==FALSE & max_missing_rt >= 2){
      warning('Maximum parity missing percentage exceeds 2% , El Badry correction is recommended!')
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

    return(avg_par)

  }
