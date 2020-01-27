#' Gompertz PF Fertility Estimation
#'
#' @param ages A vector of starting ages of five-year age groups ranging from 15 to 45 (default = c(15,20,25,30,35,40,45))
#' @param P A vector of mean parities by five-year age group - same groups as 'ages'
#' @param asfr A vector of age-specific fertility rates by five-year age group - same groups as 'ages'
#' @param adjust_group A vector of age-groups from 'ages' to be used for selection of PF ratios to adjust asfr data
#' (default set to 20 (20-24 five-year age group))
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
  function(ages,
           asfr,
           P,
           age_shift   = TRUE,
           shape       = TRUE,
           level       = TRUE,
           graph_check = FALSE,
           rmse_check  = TRUE){


    adjustGompInput <-
      function(ages, asfr, P){

        # 1. Check if inputs have the correct dimensions
        stopifnot( all.equal( length(ages), length(P), length(asfr)) )

        # 2. Adjust data inputs for Gompertz application if asfr and P vectors are given for ages 15+
        if(length(asfr) == 7){
          asfr <-
            c( NA, asfr )
        }

        if(length(P) == 7){
          P <-
            c( NA, P )
        }

        # 3. Provide lower bound and upper bound
        age.lb <-
          seq( 10, 45, 5)

        age.ub <-
          seq( 15, 50, 5)

        age.group <-
          c( "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49" )

        # 4. return adjusted data for gompertz method
        gomp_data <-
          data.frame(
            age.group,
            age.lb,
            age.ub,
            asfr,
            P
          )

      return(gomp_data)
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

    gomp_data <- adjustGompInput(ages, asfr, P)

    if(age_shift){
      gomp_data$age_shifted <-
        gomp_data$age.lb + 4.5

      gomp_data <-
        merge(
          gomp_data,
          std.zaba[, c( 'age', 'Fx' )],
          by.x = 'age_shifted',
          by.y = 'age',
          all.x = TRUE
        )
    } else{
      gomp_data <-
        merge(gomp_data,
              std.zaba[, c( 'age', 'Fx' )],
              by.x = 'age.lb',
              by.y = 'age',
              all.x = TRUE
        )
    }

    # change name of Fx
    names(gomp_data)[ncol(gomp_data)] <-
      'Fx_std'

  # for P no age shift adjustment is needed
    gomp_data <-
      merge(gomp_data,
            std.zaba[, c( 'age', 'Px_x5' )],
            by.x = 'age.ub',
            by.y = 'age',
            all.x = TRUE
      )

    # change name of Px_x5
    names(gomp_data)[ncol(gomp_data)] <-
      'Px_x5_std'

    # add Fx_std_noshift for exact ages estimation of cumulated fertility Fx_adjust
    gomp_data <-
      merge(gomp_data,
            std.zaba[, c( 'age', 'Fx' )],
            by.x = 'age.ub',
            by.y = 'age',
            all.x = TRUE
      )
    # change name of Fx
    names(gomp_data)[ncol(gomp_data)] <-
      'Fx_std_noshift'

    # 1) Generate estimates for gx, ex and zx for F points #----------
    # First Gompi
    data_F <- gomp_data

    data_F$Y_std_F <-
      -log(-log(data_F$Fx_std))

    # Ratio of adjacent groups Fx/Fx+5
    data_F$Fx_x_5_std <-
      NA
    for(i in 1:7){
      data_F$Fx_x_5_std[i] = data_F$Fx_std[i]/data_F$Fx_std[i+1]
    }

    # Gompi of ratios
    data_F$phi_F <-
      -log(-log(data_F$Fx_x_5_std))

    # Parameter phi' = g(x)
    # phi'= (exp(Ys.x)*Ys.x+5)-Ys.x*exp(Ys.x+5))/(exp(Ys.x)-exp(Ys.x+5))

    data_F$phi_1_F <-
      NA
    for(i in 1:7){
      data_F$phi_1_F[i] = ((exp(data_F$Y_std_F[i])*data_F$Y_std_F[i+1])-data_F$Y_std_F[i]*exp(data_F$Y_std_F[i+1]))/
        (exp(data_F$Y_std_F[i])-exp(data_F$Y_std_F[i+1]))
    }

    data_F$gx <-
      data_F$phi_1_F

  # Parameter phi'' - only for 15-30 years old, the mean gives value of parameter c
  # phi''= (Ys.x-Ys.x+5)?*exp(Ys.x+Ys.x+5)/(exp(Ys.x)-exp(Ys.x+5))?

  data_F$phi_2_F <-
    NA

  for(i in 2:4){
    data_F$phi_2_F[i]=((data_F$Y_std_F[i]-data_F$Y_std_F[i+1])^2*exp(data_F$Y_std_F[i]+data_F$Y_std_F[i+1]))/
      (exp(data_F$Y_std_F[i])-exp(data_F$Y_std_F[i+1]))^2
  }

  c_F <- mean(data_F$phi_2_F, na.rm = T)

  # e(x) = difference between ratios gompit (phi) and phi'
  data_F$ex <-
    data_F$phi_F - data_F$phi_1_F

  ### Parameter z(x) - based on observed data

  # asfr - age specific fertility rates (10-49)
  # Pi - mean number of children ever born from women at age group i

  # cumulate asfr
  data_F$Fx_obs <-
    c(0,5*cumsum(na.omit(data_F$asfr)))

  # Fx ratios
  data_F$Fx_x5_obs <-
    NA

  for(i in 1:7){
    data_F$Fx_x5_obs[i]=data_F$Fx_obs[i]/data_F$Fx_obs[i+1]
  }

  # estimating z(x), gompi from cumulated ratios
  data_F$zx <-
    c( NA, (-log(-log(data_F$Fx_x5_obs[2:7]))), NA )

  ####################

  # 2) Generate estimates for gx, ex and zx for P points #--------
  data_P <-
    gomp_data

  # First Gompi
  data_P$Y_std_P <-
    -log(-log(data_P$Px_x5_std))

  # ratios P(i)/P(i+1)
  data_P$ratio_P <-
    NA
  for(i in 1:7){
    data_P$ratio_P[i] = data_P$Px_x5_std[i]/data_P$Px_x5_std[i+1]
  }

  # For last age group, Pi/1
  data_P$ratio_P[8] <-
    data_P$Px_x5_std[8]/1

  # Gompi of ratios phi_P
  data_P$phi_P <-
    -log(-log(data_P$ratio_P))

  # parameter phi' = gx_P
  # phi'= (exp(Ys.i)*Ys.i+1)-Ys.i*exp(Ys.i+1))/(exp(Ys.i)-exp(Ys.i+1))
  data_P$phi_1_P <-
    NA

  for(i in 1:7){
    data_P$phi_1_P[i] = ((exp(data_P$Y_std_P[i])*data_P$Y_std_P[i+1])-data_P$Y_std_P[i]*exp(data_P$Y_std_P[i+1]))/
      (exp(data_P$Y_std_P[i])-exp(data_P$Y_std_P[i+1]))
  }

  data_P$phi_1_P[8] <-
    data_P$Y_std_P[8]

  data_P$gx <-
    data_P$phi_1_P

  # Parameter phi'' - only for 15-30 years old, the mean gives value of parameter c
  # phi''= (Ys.x-Ys.x+5)?*exp(Ys.x+Ys.x+5)/(exp(Ys.x)-exp(Ys.x+5))?

  data_P$phi_2_P <-
    NA

  for(i in 2:4){
    data_P$phi_2_P[i] = ((data_P$Y_std_P[i]-data_P$Y_std_P[i+1])^2*exp(data_P$Y_std_P[i]+data_P$Y_std_P[i+1]))/
      (exp(data_P$Y_std_P[i])-exp(data_P$Y_std_P[i+1]))^2
  }

  c_P <- mean(data_P$phi_2_P, na.rm = T)

  # e(x) = difference between ratios gompit (phi) and phi'
  data_P$ex <-
    data_P$phi_P - data_P$phi_1_P

  ### Parameter z(x) - based on observed data

  # asfr - age specific fertility rates (10-49)
  # Pi - mean number of children ever born from women at age group i

  # cumulate asfr
  data_P$Px_x5_obs <-
    NA

  for(i in 1:7){
    data_P$Px_x5_obs[i] = gomp_data$P[i]/gomp_data$P[i+1]
  }

  # estimating z(x), gompi from P ratios
  data_P$zx <-
    c( NA, (-log(-log(data_P$Px_x5_obs[2:7]))), NA )
  ####################################################

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
