#' Gompertz PF Fertility Estimation
#'
#' @param ages A vector of starting ages of five-year age groups ranging from 15 to 45 (default = c(15,20,25,30,35,40,45))
#' @param P A vector of mean parities by five-year age group - same groups as 'ages'
#' @param asfr A vector of age-specific fertility rates by five-year age group - same groups as 'ages'
#' @param level TRUE for correction of fertility level using parity data information, false for correction of fertility shape only
#' @param madef Mother's age definition: '0m' for age at birth of child, '12m' for age at survey for 12 months data (default),
#' '24m' for age at survey for 24 months data, '36m' for age at survey for 36 months data
#'
#' @return A list with 3 elements:
#' pf_data data frame with columns ages, P for mean parities, asfr, Fi for cumulate fertility estimated from Brass coefficients,PF for ratios P/F and adj_asfr for adjusted asfr;
#' tfr_unadj for unadjusted total fertility rate estimate;
#' and tfr_adj for adjusted total fertility rate estimate by applying the selected age-group PF ratio
#' @export
#' @source
#' Brass W, AJ. 1968. Coale Methods  of  analysis  and  estimation.  In:  BRASS,  W.  et  al.  (Ed.).  The demography of tropical Africa. 1. ed. New Jersey: Princeton University Press, p. 88-139.
#' Brass W. 1975. Methods for Estimating Fertility and Mortality from Limited and Defected Data. North Carolina: Carolina Population Center.
#' @examples
#' ## Malawi 2008 Census data:
#' ages_ma = c(15, 20, 25, 30, 35, 40, 45)
#' asfr_ma = c(0.111, 0.245, 0.230, 0.195, 0.147, 0.072, 0.032)
#' P_ma    = c(0.283, 1.532, 2.849, 4.185, 5.214, 6.034, 6.453)
#' fertBrassPF(P = P_ma, asfr = asfr_ma)
#'
#'


# 3) Gompertz function for generating estimates #---------------------
fertGompPF <-
  function( ages            = seq( 15, 45, 5 ),
            asfr,
            P               = NULL,
            level           = FALSE,
            madef           = '12m',
            sel.ages        = c( 20, 25, 30, 35, 40 ),
            plot.diagnostic = TRUE
            ){


    # 1. Adjust the inputs into the correct form for the method
    adjustGompInput <-
      function( ages,
                asfr,
                P,      # default = null
                level,  # default = FALSE
                madef   # default = '12m'
                ){


        # 1.1 Check if level = TRUE and P is provided
        if (level & is.null( P )){
          stop('Please provide vector P for mean number of children ever born by women age group')
        }

        # 1.2 Check if input lengths match
        if (level){

          stopifnot( all.equal( length(ages), length(P), length(asfr)) )

          # 1.3. Adjust data inputs for Gompertz application if asfr and P vectors are given for ages 15+
          if(length(asfr) == 7){
            asfr <-
              c( NA, asfr )
          }

          if(length(P) == 7){
            P <-
              c( NA, P )
          }

        } else{

          stopifnot( all.equal( length(ages), length(asfr)) )

          # 1.3. Adjust data inputs for Gompertz application if asfr vector is given for ages 15+
          if(length(asfr) == 7){
            asfr <-
              c( NA, asfr )
          }

          P <-
            rep( NA, 8 )
        }

        # 1.4. Provide age groups and respective lower bound, upper bound and age.shift
        age.lb <-
          seq( 10, 45, 5)

        age.ub <-
          seq( 15, 50, 5)

        age.group <-
          c( "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49" )

        if ( madef == '0m'){
          age.shift <-
            age.ub
        } else if ( madef == '12m' ){
          age.shift <-
            age.ub - 0.5
        } else if ( madef == '24m' ){
          age.shift <-
            age.ub - 1
        } else if ( madef == '36m' ){
          age.shift <-
            age.ub - 1.5
        } else {
          stop( 'Incorrect madef entry! Try either 0m, 12m, 24m or 36m' )
        }

        # 1.5. return adjusted data for gompertz method
        gomp.dat <-
          data.frame(
            age.group,
            age.lb,
            age.ub,
            age.shift,
            asfr,
            P
          )

        if ( is.null(P) | !level ){
          gomp.dat <-
            gomp.dat[ , - 6 ]
        }

      return(gomp.dat)
      }

    inputGomp.dat <-
      adjustGompInput(
        ages  = ages,
        asfr  = asfr,
        P     = P,
        level = level,
        madef = madef
        )

    # 2. Compute F gompits for estimation of ex, zx and gx points

    F.GompitCalc <-
      function(
        age.group,
        age.shift,
        age.noshift,
        asfr
      ){

        # 2.1 Merge ages and asfr data with correspondent zaba std Fx value for age shift defined in madef
        fgomp.dat <-
          merge(
            data.frame(
              age.group,
              age.shift,
              age.noshift,
              asfr
            ),
            std.zaba[ , c( 'age', 'Fx' ) ],
            by.x = 'age.shift',
            by.y = 'age',
            all.x = TRUE
          )

        names(fgomp.dat)[ncol(fgomp.dat)] <-
          'Fx.std'

        # 2.2 Merge ages and asfr data with correspondent zaba std Fx value for no shift in ages
        fgomp.dat <-
          merge(
            fgomp.dat,
            std.zaba[ , c( 'age', 'Fx' ) ],
            by.x = 'age.noshift',
            by.y = 'age',
            all.x = TRUE
          )

        names(fgomp.dat)[ncol(fgomp.dat)] <-
          'Fx.stdnoshift'

        # 2.3 Generate estimates for gx, ex and zx for F points

        ## A. First Gompit
        fgomp.dat$Yf.std <-
          - log( - log( fgomp.dat$Fx.std ) )

        ## B. Ratio of adjacent groups Fx/Fx+5
        fgomp.dat$Fx_x5.std <-
          NA

        for(i in 1:7){
          fgomp.dat$Fx_x5.std[i] <-
            fgomp.dat$Fx.std[i] / fgomp.dat$Fx.std[i+1]
        }

        ## C. Gompit of ratios equal phi
        fgomp.dat$phi <-
          - log( - log( fgomp.dat$Fx_x5.std ) )

        ## D. Parameter phi' = g(x)
        ## phi'= (exp(Ys.x)*Ys.x+5)-Ys.x*exp(Ys.x+5))/(exp(Ys.x)-exp(Ys.x+5))

        fgomp.dat$phi.1 <-
          NA

        for(i in 1:7){
          fgomp.dat$phi.1[i] <-
            ( ( exp( fgomp.dat$Yf.std[i] ) * fgomp.dat$Yf.std[i+1] ) - fgomp.dat$Yf.std[i] * exp( fgomp.dat$Yf.std[i+1] ) ) /
            ( exp( fgomp.dat$Yf.std[i] ) - exp( fgomp.dat$Yf.std[i+1] ) )
        }

        fgomp.dat$gx <-
         round( fgomp.dat$phi.1, 4 )

        ## E. Parameter phi''
        ## Parameter phi'' - only for 15-30 years old, the mean gives value of parameter c
        ## phi''= (Ys.x-Ys.x+5)?*exp(Ys.x+Ys.x+5)/(exp(Ys.x)-exp(Ys.x+5))

        fgomp.dat$phi.2 <-
          NA

        for(i in 2:4){
          fgomp.dat$phi.2[i] <-
            ( ( fgomp.dat$Yf.std[i] - fgomp.dat$Yf.std[i+1] ) ^ 2 * exp( fgomp.dat$Yf.std[i] + fgomp.dat$Yf.std[i+1] ) ) /
            ( exp( fgomp.dat$Yf.std[i] ) - exp( fgomp.dat$Yf.std[i+1] ) ) ^ 2
        }

        fgomp.dat$c.F <-
         round( mean( fgomp.dat$phi.2, na.rm = T ), 4 )

        ## F. e(x) = difference between ratios gompit (phi) and phi'
        fgomp.dat$ex <-
          round( fgomp.dat$phi - fgomp.dat$phi.1, 4 )

        ## G. Parameter z(x) - based on observed data

        ## cumulate asfr
        fgomp.dat$Fx.obs <-
          c( 0, 5 * cumsum( na.omit( fgomp.dat$asfr ) ) )

        ## Fx ratios
        fgomp.dat$Fx_x5.obs <-
          NA

        for(i in 1:7){
          fgomp.dat$Fx_x5.obs[i] <-
            fgomp.dat$Fx.obs[i] / fgomp.dat$Fx.obs[i+1]
        }

        ## estimating z(x), gompi from cumulated ratios
        fgomp.dat$zx <-
          round( c( NA, ( -log( -log( fgomp.dat$Fx_x5.obs[2:7] ) ) ), NA ), 4 )

        ## H. Return selected arguments
        fgomp.dat <-
          fgomp.dat[ , c( 'age.group', 'age.noshift', 'age.shift', 'asfr', 'Fx.std', 'Fx.stdnoshift', 'gx', 'ex', 'zx', 'c.F')]

        return(fgomp.dat)
      }

    Gomp.Fdat <-
      F.GompitCalc(
        age.group   = inputGomp.dat$age.group,
        age.shift   = inputGomp.dat$age.shift,
        age.noshift = inputGomp.dat$age.ub,
        asfr        = inputGomp.dat$asfr
        )

    # 3. Compute P gompits for estimation of ei, zi and gi points

    P.GompitCalc <-
      function(
        age.group,
        age.noshift,
        P
      ){

        # 3.1 Merge ages and P data with correspondent zaba std Px_x5 value for exact ages
        pgomp.dat <-
          merge(
            data.frame(
              age.group,
              age.noshift,
              P
            ),
            std.zaba[, c( 'age', 'Px_x5' )],
            by.x = 'age.noshift',
            by.y = 'age',
            all.x = TRUE
          )

        names(pgomp.dat)[ncol(pgomp.dat)] <-
          'P.std'

        # 3.2 Generate estimates for gi, ei and zi for P points

        ## A. First Gompi
        pgomp.dat$Yp.std <-
          - log( - log( pgomp.dat$P.std ) )

        ## B. Ratios P(i)/P(i+1)
        pgomp.dat$P.ratio <-
          NA

        for(i in 1:7){
          pgomp.dat$P.ratio[i] <-
            pgomp.dat$P.std[i] / pgomp.dat$P.std[i+1]
        }

        # For last age group, Pi/1
        pgomp.dat$P.ratio[ nrow( pgomp.dat ) ] <-
          pgomp.dat$P.std[ nrow( pgomp.dat ) ] / 1

        ## C. Gompit of ratios equal phi
        pgomp.dat$phi <-
          - log( - log( pgomp.dat$P.ratio ) )

        ## D. parameter phi' = gi
        ## phi'= (exp(Ys.i)*Ys.i+1)-Ys.i*exp(Ys.i+1))/(exp(Ys.i)-exp(Ys.i+1))

        pgomp.dat$phi.1 <-
          NA

        for(i in 1:7){
          pgomp.dat$phi.1[i] <-
            ( ( exp( pgomp.dat$Yp.std[i] ) * pgomp.dat$Yp.std[i+1]) - pgomp.dat$Yp.std[i] * exp( pgomp.dat$Yp.std[i+1] ) ) /
            ( exp( pgomp.dat$Yp.std[i] ) - exp( pgomp.dat$Yp.std[i+1] ) )
        }

        pgomp.dat$phi.1[8] <-
          pgomp.dat$Yp.std[8]

        pgomp.dat$gi <-
          pgomp.dat$phi.1

        ## E. Parameter phi'' - only for 15-30 years old, the mean gives value of parameter c
        ## phi''= (Ys.x-Ys.x+5)?*exp(Ys.x+Ys.x+5)/(exp(Ys.x)-exp(Ys.x+5))?

        pgomp.dat$phi.2 <-
          NA

        for(i in 2:4){
          pgomp.dat$phi.2[i]  <-
            ( ( pgomp.dat$Yp.std[i]- pgomp.dat$Yp.std[i+1]) ^ 2 * exp( pgomp.dat$Yp.std[i] + pgomp.dat$Yp.std[i+1] ) ) /
            ( exp( pgomp.dat$Yp.std[i] ) - exp( pgomp.dat$Yp.std[i+1] ) ) ^ 2
        }

        pgomp.dat$c.P <-
          mean( pgomp.dat$phi.2, na.rm = T )

        ## F. e(i) = difference between ratios gompit (phi) and phi'
        pgomp.dat$ei <-
          pgomp.dat$phi -  pgomp.dat$phi.1

        ## G. Parameter z(i) - based on observed data

        # x / x+5 ratios
        pgomp.dat$Px_x5.obs <-
          NA

        for(i in 1:7){
          pgomp.dat$Px_x5.obs[i] <-
            pgomp.dat$P[i] / pgomp.dat$P[i+1]
        }

        # estimating z(i), gompi from P ratios
        pgomp.dat$zi <-
          c( NA, ( - log( -log( pgomp.dat$Px_x5.obs[2:7] ) ) ), NA )

        ## H. Return selected arguments
        pgomp.dat <-
          pgomp.dat[ , c( 'age.group', 'age.noshift', 'P', 'P.std', 'gi', 'ei', 'zi', 'c.P')]

        return(pgomp.dat)
      }

    if ( level ){
      Gomp.Pdat <-
        P.GompitCalc(
          age.group   = inputGomp.dat$age.group,
          age.noshift = inputGomp.dat$age.ub,
          P           = inputGomp.dat$P
        )
    }

    # 4. Fit Alpha and Beta values

    fitGompPF <-
      function(
        age.group,
        age.ub,
        gi = NULL,
        ei = NULL,
        zi = NULL,
        gx,
        ex,
        zx,
        sel.ages,
        level = FALSE,
        plot.diagnostic,
        c.F,
        c.P = NULL
      ){

        # 4.1. Set data
        fitGomp.dat <-
          data.frame(
            age.group = rep( age.group, 2 ),
            age.ub    = rep( age.ub, 2 ),
            g         = c( gx, gi),
            e         = c( ex, ei),
            z         = c( zx, zi),
            point.lab = c( rep( 'F-Points', length(gx) ), rep( 'P-Points', length(gi) ) )
          )

        fitGomp.dat <-
          fitGomp.dat[ !fitGomp.dat$age.ub %in% c( 15, 50 ) , ] # filter extreme ages

        fitGomp.dat$y <-
          fitGomp.dat$z - fitGomp.dat$e

        fitGomp.dat$x <-
          fitGomp.dat$g

        # 4.2. Fit alpha and beta for F points

        ## A. All points first
        Fall.model <-
          lm(
            y ~ x,
            data = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]
            )

        Fall.beta  <-
          Fall.model$coefficients[2]

        Fall.intercept <-
          Fall.model$coefficients[1]

        ## B. Selected points
        Fsel.model <-
          lm(
            y ~ x,
            data = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]
          )

        Fsel.beta  <-
          Fsel.model$coefficients[2]

        Fsel.intercept <-
          Fsel.model$coefficients[1]

        Fsel.alpha <-
          Fsel.intercept - 0.5 * ( c.F ) * ( Fsel.beta - 1 ) ^ 2

        # 4.3. Fit alpha and beta for P points and joint F and P points if required
        if ( level ){

          ## A. All points first
          Pall.model <-
            lm(
              y ~ x,
              data = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points', ]
            )

          Pall.beta  <-
            Pall.model$coefficients[2]

          Pall.intercept <-
            Pall.model$coefficients[1]

          ## B. Selected points
          Psel.model <-
            lm(
              y ~ x,
              data = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' & fitGomp.dat$age.ub %in% sel.ages, ]
            )

          Psel.beta  <-
            Psel.model$coefficients[2]

          Psel.intercept <-
            Psel.model$coefficients[1]

          Psel.alpha <-
            Psel.intercept - 0.5 * ( c.P ) * ( Psel.beta - 1 ) ^ 2

          ## C. F and P selected points

          FPsel.model <-
            lm(
              y ~ x,
              data = fitGomp.dat[ fitGomp.dat$age.ub %in% sel.ages, ]
            )

          FPsel.beta  <-
            FPsel.model$coefficients[2]

          FPsel.intercept <-
            FPsel.model$coefficients[1]

          FPsel.alpha <-
            FPsel.intercept - 0.5 * mean( c( c.F, c.P ) ) * ( FPsel.beta - 1 ) ^ 2
        }

        # 4.4 Plot diagnostic graphs if required
        if ( plot.diagnostic ){

          if ( level ){

            x11( width = 9, height = 6)
            par( mfrow = c( 1, 2 ) )

            ## A.1 All F points
            plot(
              x    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' ,]$x,
              y    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' ,]$y,
              pch  = 19,
              col  = 'skyblue',
              xlab = 'g()',
              ylab = 'z()-e()'
            )
            abline( reg = Fall.model, col = 'skyblue' , lty = 5)

            ## A.2 All P points
            points(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' ,]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' ,]$y,
              pch = 15,
              col = 'tomato3'
            )
            abline( reg = Pall.model, col = 'tomato3' , lty = 3)
            legend(
              'bottomright',
              legend = c( 'F-points', 'P-points' ),
              lty    = c(  5,  3 ),
              pch    = c( 19, 15 ),
              col    = c( 'skyblue', 'tomato3' ),
              bty    = 'n'
            )
            grid()
            mtext(
              side = 3,
              line = 2.35,
              adj  = 0,
              cex  = 1.15,
              'Figure 1: All F Points and P points'
              )
            mtext(
              side = 3,
              line = 0.60,
              adj  = 0,
              cex  = 0.75,
              paste0( 'F-points linear: y(x) = ', round( Fall.intercept, 4) , ' + ', round( Fall.beta, 4 ) , ' * x\n',
                      'P-points linear: y(x) = ', round( Pall.intercept, 4) , ' + ', round( Pall.beta, 4 ) , ' * x' )
              )
            text(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]$y,
              labels = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]$age.group,
              cex = 0.75,
              pos = 4,
              col = 'skyblue'
            )
            text(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points', ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points', ]$y,
              labels = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points', ]$age.group,
              cex = 0.75,
              pos = 4,
              col = 'tomato3'
            )

            ## B.1 Selected F points
            plot(
              x    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$x,
              y    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$y,
              pch  = 19,
              col  = 'skyblue',
              xlab = 'g()',
              ylab = 'z()-e()'
            )
            abline( reg = Fsel.model, col = 'skyblue' , lty = 5)

            ## B.2 Selected P points
            points(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$y,
              pch = 15,
              col = 'tomato3'
            )
            abline( reg = Psel.model, col = 'tomato3' , lty = 3)

            ## B.3 Combined P and F points
            abline( reg = FPsel.model, col = 'black' , lty = 1, lwd = 0.5)
            legend(
              'bottomright',
              legend = c( 'F-points', 'P-points', 'Combined F and P' ),
              lty    = c(  5,  3 , 1),
              pch    = c( 19, 15 , NA),
              col    = c( 'skyblue', 'tomato3', 'black' ),
              bty    = 'n'
            )
            grid()
            mtext(
              side = 3,
              line = 2.35,
              adj  = 0,
              cex  = 1.15,
              'Figure 2: Selected F Points and P points'
            )
            mtext(
              side = 3,
              line = 0.10,
              adj  = 0,
              cex  = 0.75,
              paste0( 'F-points linear: y(x) = ', round( Fsel.intercept, 4) , ' + ', round( Fsel.beta, 4 ) , ' * x\n',
                      'P-points linear: y(x) = ', round( Psel.intercept, 4) , ' + ', round( Psel.beta, 4 ) , ' * x\n',
                      'Combined F and P linear: y(x) = ', round( FPsel.intercept, 4) , ' + ', round( FPsel.beta, 4 ) , ' * x' )
            )
            text(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$y,
              labels = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$age.group,
              cex = 0.75,
              pos = 4,
              col = 'skyblue'
            )
            text(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'P-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$y,
              labels = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$age.group,
              cex = 0.75,
              pos = 4,
              col = 'tomato3'
            )

          }

          else{
            x11( width = 9, height = 6)
            par( mfrow = c( 1, 2 ) )

            ## A.1 All F points
            plot(
              x    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' ,]$x,
              y    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' ,]$y,
              pch  = 19,
              col  = 'skyblue',
              xlab = 'g()',
              ylab = 'z()-e()'
            )
            abline( reg = Fall.model, col = 'skyblue' , lty = 5)
            grid()
            mtext(
              side = 3,
              line = 2,
              adj  = 0,
              cex  = 1.25,
              'Figure 1: All F Points'
            )
            mtext(
              side = 3,
              line = 0.75,
              adj  = 0,
              cex  = 1,
              paste0( 'F-points linear: y(x) = ', round( Fall.intercept, 4) , ' + ', round( Fall.beta, 4 ) , ' * x')
            )
            text(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]$y,
              labels = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points', ]$age.group,
              cex = 0.75,
              pos = 4,
              col = 'skyblue'
            )

            ## B.1 Selected F points
            plot(
              x    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$x,
              y    = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$y,
              pch  = 19,
              col  = 'skyblue',
              xlab = 'g()',
              ylab = 'z()-e()'
            )
            abline( reg = Fsel.model, col = 'skyblue' , lty = 5)
            grid()
            mtext(
              side = 3,
              line = 2,
              adj  = 0,
              cex  = 1.25,
              'Figure 2: Selected F Points'
            )
            mtext(
              side = 3,
              line = 0.75,
              adj  = 0,
              cex  = 1,
              paste0( 'F-points linear: y(x) = ', round( Fsel.intercept, 4) , ' + ', round( Fsel.beta, 4 ) , ' * x')
            )
            text(
              x = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$x,
              y = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$y,
              labels = fitGomp.dat[ fitGomp.dat$point.lab == 'F-Points' & fitGomp.dat$age.ub %in% sel.ages, ]$age.group,
              cex = 0.75,
              pos = 4,
              col = 'skyblue'
            )
          }

        }
      }
    # 2) Function to select points either from rmse or graphically #-----
    diagnostic_function <- function(data_F,data_P,graph_check=F,rmse_check=T,c_F,c_P){
      require(hydroGOF)
      require(data.table)
      require(ggplot2)
      require(dplyr)

      # Check the values of dummies for correct selection of chosen criteria
      if((graph_check==F & rmse_check==F)|(graph_check==T & rmse_check==T)){
        stop('One option only must be chosen between rmse_check and graph_check')
      }

      # Create data for both points datasets
      data_PF <- rbind(data_F[2:7,.(AGE_GROUP,gx,ex,zx,points='F-Points')],
                       data_P[2:7,.(AGE_GROUP,gx,ex,zx,points='P-Points')]) %>%
        .[,`:=`(y=zx-ex,
                x=gx)] %>%
        .[,group:=paste0(AGE_GROUP,"-",points)]%>%
        .[,.(AGE_GROUP,group,y,x,points)]%>%
        .[y!=Inf,]

      # First testing RMSE criteria #-----------
      if(rmse_check==T){
        # Testing all points first
        model = lm(y~x,data=data_PF)
        beta  = model$coefficients[2]
        alpha = model$coefficients[1]
        data_PF[,pred := x*beta+alpha]
        data_PF[,RMSE:=NA]

        for(i in 1: nrow(data_PF)){
          data_PF$RMSE[i] = round((data_PF$pred[i]-data_PF$y[i])^2,5)
        }

        RMSE_test = sum(trunc(data_PF$RMSE*100)>0)>0 # test if RMSE results has elements higher than 0.00
        data_PF_new <- data_PF
        WITHIN_BOUNDS = (alpha<0.35 & alpha>-0.35) & (beta>0.75 & beta<1.30) # dummy for checking alfa and beta bounds

        # using rmse criteria to update alpha and betas reducing problematic points
        for(rmv_seq in c("40-44-P-Points","40-44-F-Points","35-39-P-Points","35-39-F-Points")){

          if(RMSE_test==T | WITHIN_BOUNDS==F){
            data_PF_new <- data_PF_new[group!=rmv_seq,]
            model = lm(y~x,data=data_PF_new)
            beta  = model$coefficients[2]
            alpha = model$coefficients[1]
            data_PF_new[,pred := x*beta+alpha]
            data_PF_new[,RMSE:=NA]

            for(i in 1: nrow(data_PF_new)){
              data_PF_new$RMSE[i] = round((data_PF_new$pred[i]-data_PF_new$y[i])^2,5)
            }

            # only one element from the list may be higher than 0.05
            RMSE_test     = sum(trunc(data_PF_new$RMSE*100)>5)>1
            # test parameters bounds (alpha and beta)
            WITHIN_BOUNDS = (alpha<0.35 & alpha>-0.35) & (beta>0.75 & beta<1.30)

          } else{
            break
          }
        }

        data_PF <- data_PF_new
        if(WITHIN_BOUNDS == F){
          warning("Alpha and Beta out of bounds, try visual diagnostic for better adjustment")
        }
        ########################

        # adjust alpha and betas with parameters from selected points
        model_PF <- model
        model_F  <- lm(y~x,data_PF_new[points=="F-Points",])
        model_P  <- lm(y~x,data_PF_new[points=="P-Points",])
        beta_F   <- model_F$coefficients[2]
        beta_P   <- model_P$coefficients[2]
        alpha_F  <- model_F$coefficients[1] - (beta_F-1)^(2)*(c_F)/2
        alpha_P  <- model_P$coefficients[1] - (beta_P-1)^(2)*(c_P)/2
        beta_PF  <- model_PF$coefficients[2]
        alpha_PF <- model_PF$coefficients[1] - (beta_PF-1)^(2)*mean(c_F,c_P)/2

        coefficients <- data.table(alpha_PF,beta_PF,alpha_F,beta_F,alpha_P,beta_P)

        output = list(coefficients=coefficients,data_PF=data_PF[,.(AGE_GROUP,y,x,points,pred,RMSE)])

        return(output)
      }
      ##########################################

      # Testing graphical criteria #-----------
      if(graph_check==T){
        # Testing all points first
        model = lm(y~x,data=data_PF)
        beta  = model$coefficients[2]
        alpha = model$coefficients[1]
        data_PF[,pred := x*beta+alpha]
        data_PF[,RMSE:=NA]

        for(i in 1: nrow(data_PF)){
          data_PF$RMSE[i] = round((data_PF$pred[i]-data_PF$y[i])^2,5)
        }

        # print initial parameters and RMSE with all points
        print(alpha)
        print(beta)
        print(data_PF)

        lim_x_F <- c(min(data_PF[points=="F-Points",]$x,na.rm=T) %>% floor,
                     max(data_PF[points=="F-Points",]$x,na.rm=T) %>% ceiling + 0.3)
        lim_y_F <- c(min(data_PF[points=="F-Points",]$y,na.rm=T) %>% floor,
                     max(data_PF[points=="F-Points",]$y,na.rm=T) %>% ceiling + 0.3)
        lim_x_P <- c(min(data_PF[points=="P-Points",]$x,na.rm=T) %>% floor,
                     max(data_PF[points=="P-Points",]$x,na.rm=T) %>% ceiling + 0.3)
        lim_y_P <- c(min(data_PF[points=="P-Points",]$y,na.rm=T) %>% floor,
                     max(data_PF[points=="P-Points",]$y,na.rm=T) %>% ceiling + 0.3)
        lim_x_FP <- range(lim_x_F,lim_x_P)
        lim_y_FP <- range(lim_y_F,lim_y_P)

        # plot all points
        x11(height = 10,width = 10)
        plot(x=data_PF[points=="F-Points",]$x,y=data_PF[points=="F-Points",]$y,
             xlim=lim_x_FP,
             ylim=lim_y_FP,
             col='red',pch=19,
             xlab='g()',ylab = 'z()-e()',
             cex=2)
        points(x=data_PF[points=="P-Points",]$x,y=data_PF[points=="P-Points",]$y, col='blue',pch=15,cex=2)
        abline(a=alpha,b=beta,col="black",lty="dashed")
        grid(col = "gray70", lty = "dotted", equilogs = TRUE)
        legend('topleft',c('F-Points', 'P-Points'),col=c('red','blue'),pch=c(19,15),cex=1.5,bty="n",horiz=T)

        point_sel_incomplete <- T

        while(point_sel_incomplete){
          x11(height = 10,width = 10)
          plot(x=data_PF[points=="F-Points",]$x,y=data_PF[points=="F-Points",]$y,
               xlim=lim_x_FP,
               ylim=lim_y_FP,
               col='red',pch=19,
               xlab='g()',ylab = 'z()-e()',
               cex=2)
          points(x=data_PF[points=="P-Points",]$x,y=data_PF[points=="P-Points",]$y, col='blue',pch=15,cex=2)
          #abline(a=alpha,b=beta,col="black",lty="dashed")
          grid(col = "gray70", lty = "dotted", equilogs = TRUE)
          legend('topleft',c('F-Points', 'P-Points'),col=c('red','blue'),pch=c(19,15),cex=1.5,bty="n",horiz=T)

          # select points to fitt alpha and beta
          sel_points <- identify(x=data_PF$x,y=data_PF$y)

          # new data with selected points
          data_PF_new <- data_PF[sel_points]

          # new model and parameters
          model_new <- lm(y~x,data=data_PF_new)
          beta  <- model_new$coefficients[2]
          alpha <- model_new$coefficients[1]
          data_PF_new[,pred := x*beta+alpha]
          beta_PF  <- beta
          alpha_PF <- alpha - (beta_PF-1)^(2)*mean(c_F,c_P)/2

          # add abline to plot
          abline(a=alpha,
                 b=beta,col="black",lty="dashed")

          # compute RMSE
          data_PF_new[,RMSE:=NA]

          for(i in 1: nrow(data_PF_new)){
            data_PF_new$RMSE[i] = round((data_PF_new$pred[i]-data_PF_new$y[i])^2,5)
          }

          # print results of alpha, beta and rmse
          print(paste0("alpha = ",alpha_PF %>%round(3)))
          print(paste0("beta  = ",beta_PF %>%round(3)))
          print(data_PF_new[,.(AGE_GROUP,points,RMSE)])

          point_sel_check = "x"
          while(point_sel_check!="y"|point_sel_check!="n"){
            # check if customer is satisfied with point selection
            point_sel_check <- readline("Are you done with point selection?(y=yes/n=no)----->\n")
            if(point_sel_check=="y"){
              point_sel_incomplete = FALSE
              break
            }
            if(point_sel_check=="n"){
              point_sel_incomplete = TRUE
              break
            }
            if(point_sel_check!="y"|point_sel_check!="n"){
              cat("Try again! y or n!\n")
            }

          }

        }

        # adjust alpha and betas with parameters from selected points
        model_PF <- lm(y~x,data_PF_new)
        model_F  <- lm(y~x,data_PF_new[points=="F-Points",])
        model_P  <- lm(y~x,data_PF_new[points=="P-Points",])
        beta_F   <- model_F$coefficients[2]
        beta_P   <- model_P$coefficients[2]
        alpha_F  <- model_F$coefficients[1] - (beta_F-1)^(2)*(c_F)/2
        alpha_P  <- model_P$coefficients[1] - (beta_P-1)^(2)*(c_P)/2
        beta_PF  <- model_PF$coefficients[2]
        alpha_PF <- model_PF$coefficients[1] - (beta_PF-1)^(2)*mean(c_F,c_P)/2

        coefficients <- data.table(alpha_PF,beta_PF,alpha_F,beta_F,alpha_P,beta_P)

        output = list(coefficients=coefficients,data_PF=data_PF[,.(AGE_GROUP,y,x,points,pred,RMSE)])

        return(output)
      }

    }
    ######################################################################




  # 3) Run diagnostic function to find alpha and beta parameters for estimation #------
  diagnostic_outputs <- diagnostic_function(data_F,data_P,graph_check = graph_check,rmse_check = rmse_check,c_F,c_P)
  alpha_PF <- diagnostic_outputs$coefficients$alpha_PF
  beta_PF  <- diagnostic_outputs$coefficients$beta_PF
  alpha_F <- diagnostic_outputs$coefficients$alpha_F
  beta_F  <- diagnostic_outputs$coefficients$beta_F
  selected_points <- diagnostic_outputs$data_PF
  ###############################################################################

  # 4) Reestimating P #-----------
  # using alpha and beta estimated to compute Y ajusted of P based on standard gompit
  data_P[,Y_adjust_P:=alpha_PF+beta_PF*Y_std_P]

  # anti-gompit to get the ratio P
  data_P[,Px_x5_adjust:=exp(-exp(-Y_adjust_P))]

  ## Achar o nivel P atraves do relacionamento de informacao do P-medio do periodo com os Ps proporcionais, estimados pelo alfa e beta, entre as idades 15-44 anos (pode testar outras idades)
  adjust_level = mean(data_P$Pi[2:7]/data_P$Px_x5_adjust[2:7])
  data_obs$Pi_adjusted = round(data_P$Px_x5_adjust*adjust_level,3)
  ################################

  # 5) Reestimating F #-------------
  data_F[,Y_std_F_noshift:=-log(-log(Fx_std_noshift))] # exact ages
  data_F[,Y_adjust_F:=alpha_PF+beta_PF*Y_std_F_noshift]

  # anti-gompit
  data_F[,Fx_x5_adjust:=round(exp(-exp(-Y_adjust_F)),5)]

  # adjust by level given by P
  data_F[,Fx_adjust:=round(Fx_x5_adjust*adjust_level,5)]

  # get the ASFR new
  data_obs$ASFR_adjusted=NA
  data_obs$ASFR_adjusted[1]=data_F$Fx_adjust[1]/5

  ## subtract adjacent groups and divide by exposition time (5 years)
  for(i in 2:8){data_obs$ASFR_adjusted[i]=(data_F$Fx_adjust[i]-data_F$Fx_adjust[i-1])/5}
  ##################

  # 6) Generate P/F series #------
  if(age_shift==F){
    data_obs[,AGE_GROUP_SHIFT:=AGE5+4.5]
    data_PF_serie <- merge(data_F[,.(AGE_GROUP_SHIFT)],
                           zaba[,.(AGE_GROUP_SHIFT=age,Fx_std=Fx)],
                           all.x = T)
    data_PF_serie[,Y_std_F:=-log(-log(Fx_std))]
  } else{
    data_PF_serie <- data_F[,.(AGE_GROUP_SHIFT,Fx_std,Y_std_F)]
  }

  data_PF_serie[,Fx_x5_adjust:=exp(-exp(-(alpha_PF+beta_PF*Y_std_F)))]
  data_PF_serie[,Fx_adjust:=Fx_x5_adjust]

  data_PF_serie[1,ASFR_adjusted:=Fx_adjust/5]
  for(i in 2:8){
    data_PF_serie$ASFR_adjusted[i]=(data_PF_serie$Fx_adjust[i]-data_PF_serie$Fx_adjust[i-1])/5
  }

  # vector of cumulant Fs from observed data
  data_PF_serie[,Fi:=NA]
  data_PF_serie$Fi[2:8]=cumsum(data_obs$ASFR[2:8])

  data_PF_serie[,F_cumulant:=NA]
  for (i in 2:7) {
    data_PF_serie$F_cumulant[i] <- 5/(data_PF_serie$Fx_x5_adjust[i])*(data_PF_serie$Fi[i])
  }

  ## correct fertility level defined by F
  adjust_level_F=mean(na.omit(data_PF_serie$F_cumulant))

  ### Get Fs for PF series construction
  data_PF_serie[,age_shift:=seq(12.5,47.5,5)]
  data_PF_serie[,Pi:=data_obs$Pi]

  #Fx*(exp(-exp(-(alfa+beta*Yxp))))
  data_PF_serie <- merge(data_PF_serie,
                         zaba[,.(age_shift=age,Y_std_F_shift=Yx_std)],
                         by="age_shift",
                         all.x=T)

  data_PF_serie[,Fi := adjust_level_F*(exp(-exp(-(alpha_F+beta_F*Y_std_F_shift))))]
  data_PF_serie[,PF:=round(Pi/Fi,3)]


  ###################################

  # 7) Generate results #--------
  results = data.table(AGE           = seq(15,45,5),
                       ASFR          = data_obs$ASFR[2:8],
                       ASFR_adjusted = data_obs$ASFR_adjusted[2:8]) %>%
    .[,TFT_obs:=5*sum(ASFR)]%>%
    .[,TFT_adj:=5*sum(ASFR_adjusted)]

  output <- list(results=results,
                 selected_points=selected_points,
                 parameters = data.table(alpha=alpha_PF,beta=beta_PF),
                 PF_series = data_PF_serie[2:8,.(AGE=seq(15,45,5),age_shift,Fi=round(Fi,3),Pi,PF)])

  ###################################
  return(output)
}
