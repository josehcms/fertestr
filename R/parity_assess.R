#########################################
### Author  : Jose H C Monteiro da Silva
### Updated : NOV-21-2019
#########################################

#' Average Parity Computation Function
#'
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#' @param na_code A numeric value representing the label for missing parities (default = NA)
#' @param prtyElBadry.set TRUE or FALSE for correction of zeros and missing values by El-Badry method (default = FALSE)
#' @param prtyElBadry.graph TRUE of FALSE for El-Badry diagnose plot output (default = FALSE)
#' @param prtyImp.set TRUE or FALSE for correction of implausible parities (default = FALSE)
#' @param prtyAssess.set TRUE or FALSE for parity assessment of average parity results (default = FALSE)
#' @param prtyAssess.graph TRUE or FALSE for parity assessment plot output (default = FALSE)
#' @param age_group Character assuming values 'q' (default) for quinquennial age group or 's' for single age group input
#'
#' @return A data.frame with 2 variables: ages and P for average parities by age group
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#'
#'   Feeney G. 1991. "Child survivorship estimation: Methods and data analysis",
#'   Asian and Pacific Population Forum 5(2-3):51-55, 76-87. http://hdl.handle.net/10125/3600.
#' @examples
#'
#' library(DemoToolsData)
#' 
#' ### Kenya 1989 data:
#' prtyAverage( ages = data.prty_KEN$ages,
#'              parity = data.prty_KEN$parity,
#'              women = data.prty_KEN$women
#' )
#'
#' # With El-Badry:
#' prtyAverage(  ages = data.prty_KEN$ages,
#'               parity = data.prty_KEN$parity,
#'               women = data.prty_KEN$women,
#'               prtyElBadry.set = TRUE,
#'               prtyElBadry.graph = TRUE )
#'
#' # Correction for implausible parities
#' prtyAverage(  ages = data.prty_KEN$ages,
#'               parity = data.prty_KEN$parity,
#'               women = data.prty_KEN$women,
#'               prtyImp.set = TRUE )
#' ###
#' ### Cambodia 2008 data:
#' prtyAverage( ages = data.prty_KHM$ages,
#'              parity = data.prty_KHM$parity,
#'              women = data.prty_KHM$women )
#'
#'


prtyAverage <-
  function(
           ages,
           parity,
           women,
           na_code           = NA,
           prtyElBadry.set   = FALSE,
           prtyImp.set       = FALSE,
           prtyElBadry.graph = FALSE,
           prtyAssess.graph  = FALSE,
           age_group         = 'q'){

    # stop with lengths are not equal
    stopifnot( all.equal( length(ages), length(parity), length(women) ) )

    # set data frame:
    data_par <-
      data.frame( ages, parity, women )

    # Change na_code for NA
    data_par$parity[data_par$parity==na_code] <- NA

    if ( prtyImp.set ){
      data_par <-
        prtyImplaus( ages   = data_par$ages,
                    parity = data_par$parity,
                    women  = data_par$women,
                    age_group
                    )

      data_par <-
        data.frame(
          ages   = data_par$ages,
          parity = data_par$parity,
          women  = data_par$women.updt
        )
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
        stats::aggregate(
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
    if ( prtyElBadry.set ){
      if ( max_missing_rt < 2 ){
        warning('Maximum parity missing percentage lower than 2% , verify whether El Badry correction is necessary!')
      }

      out.elbadry <-
        prtyElBadry( ages   = data_par$ages,
                    parity = data_par$parity,
                    women  = data_par$women,
                    prtyElBadry.graph
                    )
      data_par <-
        out.elbadry$data_par
    }
    if ( prtyElBadry.set == FALSE & max_missing_rt >= 2 ){
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

    if( prtyAssess.graph ){
      prtyAssess.plot( ages = avg_par$ages,
                      P    = avg_par$P
                      )
    }

    return(avg_par)
  }

#' El-Badry correction Function
#'
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#' @param na_code A numeric value representing the code for missing parities (default = NA)
#' @param prtyElBadry.graph TRUE of FALSE for El-Badry diagnose plot output (default = FALSE)
#'
#' @return A list with 3 elements: $data_par data.frame with corrected numbers of zero and missing parities
#' estimated by the El-Badry function, $prty.Zi data.frame with proportion of zero parities and $prty.Ui with
#' proportion of unkown parities
#' @export
#' @source
#' Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#' Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#' Population. demographicestimation.iussp.org
#'
#' el-Badry MA. 1961.
#' "Failure of enumerators to make entries of zero: errors in recording childless cases in population censuses",
#'  Journal of the American Statistical Association 56(296):909–924.
#'  doi: http://dx.doi.org/10.1080/01621459.1961.10482134
#'
#'  UN Population Division. 1983. Manual X: Indirect Techniques for Demographic Estimation.
#'  New York: United Nations, Department of Economic and Social Affairs, ST/ESA/SER.A/81.
#'  http://www.un.org/esa/population/techcoop/DemEst/manual10/manual10.html
#'
#' @examples
#'
#' library(DemoToolsData)
#'
#' ### Kenya 1989 data:
#' prtyElBadry( ages = data.prty_KEN$ages,
#'              parity = data.prty_KEN$parity,
#'              women = data.prty_KEN$women,
#'              prtyElBadry.graph = TRUE)
#'
#'
#' ### Cambodia 2008 data:
#' prtyElBadry( ages = data.prty_KHM$ages,
#'              parity = data.prty_KHM$parity,
#'              women = data.prty_KHM$women )
#' ###


prtyElBadry <-
  function(ages,
           parity,
           women,
           prtyElBadry.graph = FALSE,
           na_code  = NA){

    # stop with lengths are not equal
    stopifnot( all.equal( length(ages), length(parity), length(women) ) )

    # set data frame:
    data_par <-
      data.frame( ages, parity, women )

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

    # get data.frame of Ui
    Ui.df <-
      data.frame(
        ages = seq(15,45,5),
        Ui = round( Ui, 4 ),
        row.names = NULL
      )

    # Compute percentage of zero parities by age group
    Zi <-
      tapply( data_par[!is.na(data_par$parity) & data_par$parity == 0,]$women,
             data_par[!is.na(data_par$parity) & data_par$parity == 0,]$ages,
             sum
             ) /
      totals_1

    # get data.frame of Zi
    Zi.df <-
      data.frame(
        ages = seq(15,45,5),
        Zi = round( Zi, 4 ),
        row.names = NULL
      )

    # Take the intercept of regression Ui vs Zi
    alpha <- stats::lm(as.numeric(Ui)~as.numeric(Zi))$coefficients[1]

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

    if (prtyElBadry.graph==TRUE){
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
      graphics::text(x = as.numeric(Zi),
                     y = as.numeric(Ui),
                     labels = c('15-19','20-24','25-29','30-34','35-39','40-44','45-49'),
                     cex=1,
                     pos=4
                     )
      graphics::abline(
                  stats::lm( as.numeric(Ui) ~ as.numeric(Zi) ) ,
                  col = 'red',
                  cex = 1.25 ,
                  lty = 3
                )
      graphics::legend(
                  'topleft',
                  legend = paste0( 'Ui = ' ,
                                  round(stats::lm( as.numeric(Ui) ~ as.numeric(Zi) )$coefficients[1],3) ,
                                  ' + Zi * ',
                                  round(stats::lm( as.numeric(Ui) ~ as.numeric(Zi) )$coefficients[2],3)
                                  ),
                  col = 'red',
                  lty = 3,
                  bty = 'n'
                )

    }

    # generate list for outputs of function (zeros proportion, unknown prop and updated parities)
    out.list <-
      list(
        prty.Zi = Zi.df,
        prty.Ui = Ui.df,
        data_par = data_par
      )

    return(out.list)
  }


#' Implausible Parities correction function
#'
#' @param ages A vector of ages or starting ages of age group intervals
#' @param parity A vector of parities for each age group
#' @param women A vector of women counts by parity and age group
#' @param age_group Character assuming values 'q' (default) for quinquennial
#' age group or 's' for single age group input
#'
#' @return The original data_par data.frame with implausible parities set to 0
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#'
#' library(DemoToolsData)
#'
#' ### Kenya 1989 data:
#' prtyImplaus( ages = data.prty_KEN$ages, parity = data.prty_KEN$parity, women = data.prty_KEN$women )
#' 

# correct implausible parities
prtyImplaus <-
  function( ages,
           parity,
           women,
           age_group = 'q'){

    # stop with lengths are not equal
    stopifnot( all.equal( length(ages), length(parity), length(women) ) )

    # set data frame:
    data_par <-
      data.frame( ages, parity, women )

    if (age_group == 'q'){
      #  limit the maximum number of live births that a women may have had to one birth every 18 months from the age of 12
      data_par$women.updt <-
        ifelse ( (data_par$parity > 2 / 3 * ( data_par$ages + 5 ) - 8) & (!is.na(data_par$parity) ), 0, data_par$women )
    }
    else if (age_group == 's'){
      # limit the maximum number of live births that a women may have had to one birth every 18 months from the age of 12
      # women that exceed this patamar is put into missing parity
      data_par$women.updt <-
        ifelse ( (data_par$parity > 2 / 3 * ( data_par$ages ) - 8) & (!is.na(data_par$parity) ), 0, data_par$women )
    }
    # new missing data
    data_par$women.na <-
      ifelse ( data_par$women.updt != data_par$women,  data_par$women, 0)

    # update missing data on original data.frame
    data_par[is.na(data_par$parity),]$women.updt <-
      data_par[is.na(data_par$parity),]$women.updt +
      tapply( data_par$women.na , data_par$ages , sum)

    # set an update flag for age group
    data_par$updt.flag <-
      ifelse ( data_par$women.updt == data_par$women,  0, 1)

    return(data_par[,c('ages','parity','women','women.updt','updt.flag')])
  }


#' Parity Assessment Plot
#'
#' @param ages A vector of ages or starting ages of age group intervals
#' @param P A vector of average parities by women age group
#'
#'
#' @return Average parities vs ages plot
#' @export
#' @source
#'   Moultrie TA, RE Dorrington, AG Hill, K Hill, IM Timæus and B Zaba (eds). 2013.
#'   Tools for Demographic Estimation. Paris: International Union for the Scientific Study of
#'   Population. demographicestimation.iussp.org
#' @examples
#'
#' library(DemoToolsData)
#'
#' ### Malawi 2008 data:
#' prtyAssess.plot( ages = data.pf_MWI$ages, P = data.pf_MWI$P )
#'

prtyAssess.plot <-
  function( ages,
           P
           ){

    # stop with lengths are not equal
    stopifnot( all.equal( length(ages), length(P) ) )

    # set data frame:
    avg_par <-
      data.frame( ages, P )

    # plot data
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
    graphics::lines(
                x = as.numeric(avg_par$ages),
                y = as.numeric(avg_par$P),
                type = 'l',
                lty  = 5,
                cex  = 1.25
              )
  }

#' Parity Assessment function between two surveys
#'
#' Verify if the average number of children ever born of women birth cohorts increase between two surveys/censuses
#'
#' @param ages.t1 vector of ages or starting ages of age group intervals for census/survey year 1
#' @param ages.t2 vector of ages or starting ages of age group intervals for census/survey year 2
#' @param P.t1 vector of average parities by women age group for census/survey year 1
#' @param P.t2 vector of average parities by women age group for census/survey year 2
#' @param time.span census/survey time span (default = 10)
#'
#' @return
#'
#' @export
#' @examples
#'

prtyAssess.cohort <-
  function(ages.t1,
           ages.t2,
           P.t1,
           P.t2,
           time.span = 10){

    # stop with lengths are not equal
    stopifnot( all.equal( length(ages.t1), length(P.t1), length(ages.t2), length(P.t2) ) )

    # set data frame:
    avg_par1 <-
      data.frame( ages.t1, P.t1 )

    avg_par2 <-
      data.frame( ages.t2, P.t2 )

    # create t2 for base t1
    avg_par1$ages.t2 <-
      avg_par1$ages.t1 + time.span

    # merge data for same cohorts
    avg_par.cohort <-
      merge(
        avg_par1,
        avg_par2,
        by  = 'ages.t2'
      )

    # check cohort P diffs
    avg_par.cohort$P.diff <-
      avg_par.cohort$P.t2 - avg_par.cohort$P.t1

    # verify if P.diffs are all positive (average parities of a same cohort should not decrease in between surveys)
    verif.P <-
      prod(avg_par.cohort$P.diff > 0) == 1

    if( !verif.P ){
      warning('Average parity estimates for some cohorts do not increase in between surveys/censuses. Verify parity data!')
    }

    return(avg_par.cohort[,c('ages.t1','ages.t2','P.t1','P.t2','P.diff')])
  }
