library(wpp2019)
library(dplyr)
library(fertestr)

pop_c <-  c( 281260, 261320, 268410, 286810, 278990, 293760, 293490, 302060, 315970, 267190, 326980, 280260, 354120, 356920, 354830 )
pop_w <- c(  815930, 780320, 697160, 626430, 361650, 435880, 393760, 352520, 294280, 230200, 160590, NA )
lx_c <- c( 1.0000, 0.9320, 0.9275, 0.9228, 0.9165, 0.9125, 0.9110, 0.9094, 0.9079, 0.9063, 0.9048, 0.9032, 0.9017, 0.9001, 0.8986, 0.8970 )
lx_w <- c( 0.91381, 0.90989, 0.90492, 0.89798, 0.88893, 0.87596, 0.86029, 0.84188, 0.81791, 0.78472, 0.73735, 0.67316 )
q0_5 <-  c( 0.0683, 0.1008, 0.1189)
q15_45 <- c( 0.1946, 0.2290, 0.2674)
asfr <- c( 0.0000, 0.0418,0.1535, 0.1482, 0.1118, 0.0708, 0.0301, 0.0032 )
asfr_std <- asfr/(5 * sum(asfr) )
asfr_15prior <- c( 0.0000, 0.0533, 0.1974, 0.2144, 0.1836, 0.1332, 0.0676, 0.0134 )
asfr_std_15prior <- asfr_15prior/(5 * sum(asfr_15prior) )

FertRevSurv(ages1_c = 0:14,
            popx1_c = pop_c,
            ages5_w = seq(10,65,5),
            popx5_w = pop_w,
            lx1_c = lx_c,
            lx5_w = lx_w,
            asfr = asfr_std,
            asfr_15prior = asfr_std_15prior,
            q0_5 = q0_5,
            q15_45f = q15_45,
            date_ref = 2008)

?FertRevSurv

dic_year <-
  c( '1965-1970' = '1967.5',
    '1970-1975' = '1972.5',
    '1975-1980' = '1977.5',
    '1980-1985' = '1982.5',
    '1985-1990' = '1987.5',
    '1990-1995' = '1992.5',
    '1995-2000' = '1997.5',
    '2000-2005' = '2002.5',
    '2005-2010' = '2007.5',
    '2010-2015' = '2012.5',
    '2015-2020' = '2017.5' )


data("mxF", package = "wpp2019")
data("mxM", package = "wpp2019")
# error analysis
listSample <- mxF[ mxF$country_code > 900, ]$country_code %>% unique

sel_list = NULL
# dados de mortalidade com repeticoes
for( i in mxF$country_code %>% unique ){
  lenF = length(mxF[ mxF$country_code == i, ]$age)
  lenM = length(mxM[ mxM$country_code == i, ]$age)

  if(lenM == 22){
    sel_list = c(sel_list,i)
  }
  if(lenM != 22){
    cat( paste0('Duplicates in mxM: ', i ,' (',mxM[ mxM$country_code == i,]$name %>% unique,')\n') )
    cat( paste0( 'age vector:','\n' ) )
    cat( mxM[ mxM$country_code == i,]$age )
    cat( paste0( '\n' ) )
  }

  if(lenF != 22){
    cat( paste0('Duplicates in mxF: ', i ,' (',mxF[ mxF$country_code == i,]$name %>% unique,')\n') )
    cat( paste0( 'age vector:','\n' ) )
    cat( mxF[ mxF$country_code == i,]$age )
    cat( paste0( '\n' ) )
  }
}


## dat <-
##   rbind(
##     revSurvWpp( country_list = sel_list , year = 2010, family = 'General' ) %>%
##     data.table::as.data.table %>%
##     .[, year_est := 2010],
##     revSurvWpp( country_list = sel_list , year = 2000, family = 'General' ) %>%
##     data.table::as.data.table %>%
##     .[, year_est := 2000],
##     revSurvWpp( country_list = sel_list , year = 1990, family = 'General' ) %>%
##     data.table::as.data.table %>%
##     .[, year_est := 1990],
##     revSurvWpp( country_list = sel_list , year = 1980, family = 'General' ) %>%
##     data.table::as.data.table %>%
##     .[, year_est := 1980]
##   )

## # Issue - countries with mortality level e0 lower than 20 (Rwanda 1990), solved by letting e0<20 = 20

## baselinedat <-
##   tfr[tfr$country_code %in% sel_list,
##       c('country_code', '1965-1970', '1970-1975', '1975-1980', '1980-1985', '1985-1990',
##         '1990-1995', '1995-2000', '2000-2005', '2005-2010', '2010-2015')] %>%
##   data.table::as.data.table %>%
##   melt( id.vars = c('country_code'),
##        measure.vars = c( '1965-1970', '1970-1975', '1975-1980', '1980-1985', '1985-1990',
##                         '1990-1995', '1995-2000', '2000-2005', '2005-2010','2010-2015'),
##        variable.name = 'year',
##        value.name = 'TFR') %>%
##   .[,.( Country = country_code, year = as.numeric( dic_year[as.character(year)] ), TFR.wpp = round(TFR,4) ) ]


## datError <-
##   merge(
##     dat,
##     baselinedat,
##     by = c('Country', 'year')
##   )


## ggplot2::ggplot( data = datError ) +
##   ggplot2::geom_abline( slope = 1, intercept = 0, color = 'tomato3', size = 0.75 ) +
##   ggplot2::geom_point( aes( x = TFR, y = TFR.wpp, color = as.factor(year) ), size = 2 ) +
##   ggplot2::scale_x_continuous( breaks = seq( 0, 15, 1 ), limits = c( 0, 10 ), name = 'TFR - Reverse Survival' ) +
##   ggplot2::scale_y_continuous( breaks = seq( 0, 15, 1 ), limits = c( 0, 10 ), name = 'TFR - WPP 2019\n(mid period)' )  +
##   ggplot2::facet_wrap( ~year_est, ncol = 2 ) +
##   ggplot2::theme_classic() +
##   ggplot2::theme(
##              legend.position = 'top',
##              axis.text = element_text( color = 'black', size = 12 ),
##              panel.grid.major = element_line( color = 'gray80', linetype = 'dashed', size = 0.25 ),
##              strip.text = element_text( size = 12 ),
##              legend.text = element_text( size = 12 )
##            )

## datError[, diff := round( 100*abs(TFR.wpp - TFR)/TFR.wpp, 5 )]

## require(sf)


## codes <-
##   countrycode::codelist %>%
##   data.table::as.data.table %>%
##   .[,.( country_name = country.name.en,
##        country_code_char = iso3c,
##        country_code_num = iso3n)
##     ]

## library(sf)

## world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")


## world$Country = world$un_a3 %>% as.numeric()




## datWorld <-
##   merge(
##     world,
##     datError[year_est == 2010,.(maxdiff_country = max(diff) ), Country] %>%
##     .[,.(Country, maxdiff_country, diff_class = cut( maxdiff_country,
##                                                     breaks = c( 0, 2.5, 5, 7.5, 10, 100 ),
##                                                     labels = c('<2.5%', '[ 2.5; 5.0)','[ 5.0; 7.5)','[ 7.5;10.0)', '>10.0'),
##                                                     right = F ) )],
##     by = 'Country'
##   )


## grDevices::x11()
## ggplot2::ggplot(data = datWorld) +
##   ggplot2::geom_sf(aes(fill = diff_class)) +
##   ggplot2::labs(title = 'Maximum relative difference for TFR estimated by fertestr from TFR of WPP2019',
##                 subtitle = 'Years of estimation: 1980, 1990, 2000, 2010') +
##   ggplot2::scale_fill_viridis_d(option = "plasma", 'Max relative difference (%)\n |TRF_est - TFR_wpp| / TFR_wpp')+
##   ggplot2::theme_classic()


## plot(y = datError$TFR.wpp, x = datError$TFR )
## graphics::abline(a = 0, b = 1)

## dat <-
##   rbind(
##     revSurvWpp( country_list = 'Brazil', year = 1970, family = 'West' ) %>%
##     data.table::as.data.table %>%
##     .[,ref:=1970],
##     revSurvWpp( country_list = 'Brazil', year = 1980, family = 'West' ) %>%
##     data.table::as.data.table %>%
##     .[,ref:=1980],
##     revSurvWpp( country_list = 'Brazil', year = 1991, family = 'West' ) %>%
##     data.table::as.data.table %>%
##     .[,ref:=1991],
##     revSurvWpp( country_list = 'Brazil', year = 2000, family = 'West' ) %>%
##     data.table::as.data.table %>%
##     .[,ref:=2000],
##     revSurvWpp( country_list = 'Brazil', year = 2010, family = 'West' ) %>%
##     data.table::as.data.table %>%
##     .[,ref:=2010]
##   )

## ggplot2::ggplot( data = dat ) +
##   ggplot2::geom_line( aes( x = as.numeric(year), y = TFR, color = as.factor(ref) ), size = 1 ) +
##   ggplot2::scale_x_continuous( breaks = seq(1955,2020,2.5), limits = c( 1955, 2010 ), name = '' ) +
##   ggplot2::scale_y_continuous( breaks = seq( 1, 7, 0.25 ), name = 'TFR' )  +
##   ggplot2::scale_color_manual( values = c( '2010' = 'black', '2000' = 'gray20',
##                                           '1991' = 'gray40', '1980' = 'gray60',
##                                           '1970' = 'gray80' ), name = ''  ) +
##   ggplot2::theme_classic() +
##   ggplot2::theme(
##              legend.position = 'top',
##              axis.text.x = element_text( color = 'black', angle = 90, vjust = 0.5, size = 12 ),
##              axis.text.y = element_text( color = 'black', size = 12 ),
##              panel.grid.major = element_line( color = 'gray90', linetype = 'dashed', size = 0.25 ),
##              strip.text = element_text( size = 12 ),
##              legend.text = element_text( size = 12 )
##            )

## plot( y = datError$diff, x = datError$year )
## dat$Type <- 'fertestr'
## baselinedat$Type <- 'wpp2019'

## ?life.table

## life.table( x = mxM[ mxM$country_code == 1830,]$age,
##            nDx = mxM[ mxM$country_code == 1830, c( paste0( 2010, '-', 2015 ) ) ],
##            nKx = rep(1,23) )
## LifeTable(x = c(0, 1, 5, seq(10,100,5)), mx = unique(mxM[ mxM$country_code == 1830, c( paste0( 2010, '-', 2015 ) ) ]))
