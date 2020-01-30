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

fertGompPF <-
  function( ages            = seq( 15, 45, 5 ),
            asfr,
            P               = NULL,
            level           = FALSE,
            madef           = '12m',
            sel.ages        = c( 20, 25, 30, 35, 40, 45 ),
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
          fgomp.dat[ , c( 'age.group', 'age.noshift', 'age.shift', 'asfr', 'Fx.std', 'Fx.stdnoshift', 'Fx.obs', 'Yf.std', 'gx', 'ex', 'zx', 'c.F')]

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
          round( pgomp.dat$phi.1, 4 )

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
          round( mean( pgomp.dat$phi.2, na.rm = T ), 4 )

        ## F. e(i) = difference between ratios gompit (phi) and phi'
        pgomp.dat$ei <-
          round( pgomp.dat$phi -  pgomp.dat$phi.1, 4 )

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
          round( c( NA, ( - log( -log( pgomp.dat$Px_x5.obs[2:7] ) ) ), NA ), 4 )

        ## H. Return selected arguments
        pgomp.dat <-
          pgomp.dat[ , c( 'age.group', 'age.noshift', 'P', 'P.std', 'Yp.std', 'gi', 'ei', 'zi', 'c.P')]

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
        if( level ){
          age.group <-
            rep( age.group, 2 )

          age.ub <-
            rep( age.ub, 2 )
        }

        fitGomp.dat <-
          data.frame(
            age.group = age.group,
            age.ub    = age.ub,
            g         = round( c( gx, gi), 4 ),
            e         = round( c( ex, ei), 4 ),
            z         = round( c( zx, zi), 4 ),
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

        # 4.5 Create output data.frame for regression and compute RMSE
        fitGomp.reg <-
          fitGomp.dat[ fitGomp.dat$age.ub %in% sel.ages, ]

        fitGomp.reg$RMSE <-
         round ( ( ( FPsel.intercept + FPsel.beta * fitGomp.reg$x ) - fitGomp.reg$y ) ^ 2, 4 )

        # 4.6 Create coefficient parameters

        if ( level ){
          coeffs.Gomp <-
            data.frame(
              F.beta  = Fsel.beta,
              F.alpha = Fsel.alpha,
              F.intercept = Fsel.intercept,
              P.beta  = Psel.beta,
              P.alpha = Psel.alpha,
              P.intercept = Psel.intercept,
              FP.beta  = FPsel.beta,
              FP.alpha = FPsel.alpha,
              FP.intercept = FPsel.intercept
            )
        } else {
          coeffs.Gomp <-
            data.frame(
              F.beta  = Fsel.beta,
              F.alpha = Fsel.alpha,
              F.intercept = Fsel.intercept
            )
        }

        # 4.7 Create output list with coefficients and point estimates
        out.fitGomp <-
          list(
            coeffs.Gomp = coeffs.Gomp,
            fitGomp.reg = fitGomp.reg
          )

        return(out.fitGomp)
      }

    if ( level ){
      coeffGomp.dat <-
        fitGompPF(
          age.group = inputGomp.dat$age.group,
          age.ub    = inputGomp.dat$age.ub,
          gi        = Gomp.Pdat$gi,
          ei        = Gomp.Pdat$ei,
          zi        = Gomp.Pdat$zi,
          gx        = Gomp.Fdat$gx,
          ex        = Gomp.Fdat$ex,
          zx        = Gomp.Fdat$zx,
          sel.ages  = sel.ages,
          level     = level,
          plot.diagnostic = plot.diagnostic,
          c.F       = unique( Gomp.Fdat$c.F ),
          c.P       = unique( Gomp.Pdat$c.P )
        )
    } else {
      coeffGomp.dat <-
        fitGompPF(
          age.group = inputGomp.dat$age.group,
          age.ub    = inputGomp.dat$age.ub,
          gx        = Gomp.Fdat$gx,
          ex        = Gomp.Fdat$ex,
          zx        = Gomp.Fdat$zx,
          sel.ages  = sel.ages,
          plot.diagnostic = plot.diagnostic,
          c.F       = unique( Gomp.Fdat$c.F )
        )
    }

    # 5. Estimate fmx for true ages (no shift) and correction level from beta and alpha

    AntiGompitCalc <-
      function(
        age.group,
        age.shift,
        age.noshift,
        level = FALSE,
        coeffs.Gomp,
        Fx.std,
        Fx.stdnoshift,
        Fx.obs,
        P.std = NULL,
        P.obs = NULL,
        sel.ages
      ){

        F.alpha <-
          coeffs.Gomp$F.alpha

        F.beta <-
          coeffs.Gomp$F.beta

        P.level <-  NULL

        if ( level & !is.null( P.std ) & !is.null( P.obs ) ){

          F.alpha <-
            coeffs.Gomp$FP.alpha

          F.beta <-
            coeffs.Gomp$FP.beta

          Pmodel.dat <-
            data.frame(
              age.group,
              age.noshift,
              P.obs,
              P.std,
              Yp.std = - log( - log( P.std ) )
            )

          Pmodel.dat$Yp.fit <-
            Pmodel.dat$Yp.std * coeffs.Gomp$FP.beta + coeffs.Gomp$FP.alpha

          Pmodel.dat$P.fit <-
            exp( - exp( - Pmodel.dat$Yp.fit ) )

          Pmodel.dat$actual.cumulant <-
            Pmodel.dat$P.obs / Pmodel.dat$P.fit

          P.level <-
            mean( Pmodel.dat$actual.cumulant[ Pmodel.dat$age.noshift %in% sel.ages ] )

        }

        Fmodel.dat <-
          data.frame(
            age.group,
            age.shift,
            age.noshift,
            Fx.obs,
            Fx.std,
            Fx.stdnoshift,
            Yf.std        = - log( - log( Fx.std ) ),
            Yf.stdnoshift = - log( - log( Fx.stdnoshift ) )
          )

        # shifted ages
        Fmodel.dat$Yf.fit <-
          Fmodel.dat$Yf.std * F.beta + F.alpha

        Fmodel.dat$Fx.fit <-
          exp( - exp( - Fmodel.dat$Yf.fit ) )

        Fmodel.dat$actual.cumulant <-
          Fmodel.dat$Fx.obs / Fmodel.dat$Fx.fit

        F.level <-
          mean( Fmodel.dat$actual.cumulant[ Fmodel.dat$age.noshift %in% sel.ages ] )

        FP.level <-
          ifelse ( is.null( P.level ),
                   F.level,
                   P.level
                   )

        Fmodel.dat$fmx <-
          Fmodel.dat$Fx.fit * FP.level / 5

        for( i in 1 : ( nrow( Fmodel.dat ) - 1 ) ){
          Fmodel.dat$fmx[ i + 1 ] <-
            ( Fmodel.dat$Fx.fit[ i + 1 ] - Fmodel.dat$Fx.fit[ i ] ) * FP.level / 5
        }

        Fmodel.dat$fmx <-
          round( Fmodel.dat$fmx, 4 )

        # conventional ages
        Fmodel.dat$Yf.fitnoshift <-
          Fmodel.dat$Yf.stdnoshift * F.beta + F.alpha

        Fmodel.dat$Fx.fitnoshift <-
          exp( - exp( - Fmodel.dat$Yf.fitnoshift ) )

        Fmodel.dat$fmx.noshift <-
          Fmodel.dat$Fx.fitnoshift * FP.level / 5

        for( i in 1 : ( nrow( Fmodel.dat ) - 1 ) ){
          Fmodel.dat$fmx.noshift[ i + 1 ] <-
            ( Fmodel.dat$Fx.fitnoshift[ i + 1 ] - Fmodel.dat$Fx.fitnoshift[ i ] ) * FP.level / 5
        }

        Fmodel.dat$fmx.noshift <-
          round( Fmodel.dat$fmx.noshift, 4 )

        Fmodel.dat$F.level <-
          F.level

        Fmodel.dat$FP.level <-
          FP.level

        outmodel.dat <-
          Fmodel.dat[, c( 'age.group', 'age.shift', 'fmx', 'age.noshift', 'fmx.noshift', 'F.level', 'FP.level' ) ]

        return( outmodel.dat )

      }

    if ( level ){
      AntiGompit.dat <-
        AntiGompitCalc(
          age.group = Gomp.Fdat$age.group,
          age.shift = Gomp.Fdat$age.shift,
          age.noshift = Gomp.Fdat$age.noshift,
          level = level,
          coeffs.Gomp = coeffGomp.dat$coeffs.Gomp,
          Fx.std = Gomp.Fdat$Fx.std,
          Fx.stdnoshift = Gomp.Fdat$Fx.stdnoshift,
          Fx.obs = Gomp.Fdat$Fx.obs,
          P.std = Gomp.Pdat$P.std,
          P.obs = Gomp.Pdat$P,
          sel.ages = sel.ages
        )
    } else {
      AntiGompit.dat <-
        AntiGompitCalc(
          age.group = Gomp.Fdat$age.group,
          age.shift = Gomp.Fdat$age.shift,
          age.noshift = Gomp.Fdat$age.noshift,
          coeffs.Gomp = coeffGomp.dat$coeffs.Gomp,
          Fx.std = Gomp.Fdat$Fx.std,
          Fx.stdnoshift = Gomp.Fdat$Fx.stdnoshift,
          Fx.obs = Gomp.Fdat$Fx.obs,
          sel.ages = sel.ages
        )
    }


    # 6. Estimate PF series
    PFseries.Calc <-
      function(
        age.group,
        age.lb,
        age.ub,
        P.obs,
        F.level,
        coeffs.Gomp
      ){

        PFseries.dat <-
          data.frame(
            age.group,
            age.mid = (age.ub - age.lb) / 2 + age.lb,
            P       = P.obs
          )

        PFseries.dat <-
          merge(
            PFseries.dat,
            std.zaba[, c( 'age', 'Yx_std' ) ],
            by.x = 'age.mid',
            by.y = 'age'
          )

        PFseries.dat$Fx <-
          round( F.level * exp( - exp( - ( coeffs.Gomp$F.beta * PFseries.dat$Yx_std + coeffs.Gomp$F.intercept ) ) ), 4 )

        PFseries.dat$PF.Ratio <-
          round( PFseries.dat$P / PFseries.dat$Fx, 4 )

        outpf.dat <-
          PFseries.dat[, c('age.group', 'age.mid', 'P', 'Fx', 'PF.Ratio')]

        return(outpf.dat)
      }

    PFseries.dat <-
      PFseries.Calc(
        age.group = inputGomp.dat$age.group,
        age.lb    = inputGomp.dat$age.lb,
        age.ub    = inputGomp.dat$age.ub,
        P.obs     = inputGomp.dat$P,
        F.level   = unique( AntiGompit.dat$F.level ),
        coeffs.Gomp = coeffGomp.dat$coeffs.Gomp
      )

        return(0)
  }




