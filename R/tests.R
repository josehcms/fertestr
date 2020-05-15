
require(wpp2019)
require(dplyr)
require(data.table)

lista =  c('Argentina', 'Brazil', 'Chile')

dat <-
  revSurvWpp( country_list = lista, year = 2010, family = 'General' )

data('tfr')

dic_year <-
  c( '1980-1985' = '1982.5',
     '1985-1990' = '1987.5',
     '1990-1995' = '1992.5',
     '1995-2000' = '1997.5',
     '2000-2005' = '2002.5',
     '2005-2010' = '2007.5',
     '2010-2015' = '2012.5',
     '2015-2020' = '2017.5' )


baselinedat <-
  tfr[tfr$name %in% lista,
      c('name', '1990-1995', '1995-2000','2000-2005','2005-2010','2010-2015')] %>%
  as.data.table %>%
  melt( id.vars = c('name'),
        measure.vars = c('1990-1995', '1995-2000','2000-2005','2005-2010','2010-2015'),
        variable.name = 'year',
        value.name = 'TFR') %>%
  .[,.( Country = name, year = as.numeric( dic_year[as.character(year)] ), TFR ) ]


dat$Type <- 'fertestr'
baselinedat$Type <- 'wpp2019'

datplot <- rbind(dat,baselinedat)

require(ggplot2)

ggplot( data = datplot ) +
  geom_line( aes( x = year, y = TFR, linetype = Type, color = Type ), size = 1.25 ) +
  geom_point( aes( x = year, y = TFR, shape = Type, color = Type ), size = 3 ) +
  scale_x_continuous( breaks = seq(1985,2020,2.5), limits = c( 1990, 2020 ), name = '' ) +
  scale_y_continuous( breaks = seq( 1, 9, 0.25 ), name = 'TFR' )  +
  scale_color_manual( values = c( 'fertestr' = 'black', 'wpp2019' = 'tomato3' ), name = ''  ) +
  scale_shape_manual( values = c( 'fertestr' = 16, 'wpp2019' = 19 ), name = ''  ) +
  scale_linetype_manual( values = c( 'fertestr' = 'solid', 'wpp2019' = 'dashed' ), name = '' ) +
  facet_wrap( ~Country, ncol = 3, scales = 'free_y' ) +
  theme_classic() +
  theme(
    legend.position = 'top',
    axis.text.x = element_text( color = 'black', angle = 90, vjust = 0.5, size = 12 ),
    axis.text.y = element_text( color = 'black', size = 12 ),
    panel.grid.major = element_line( color = 'gray90', linetype = 'dashed', size = 0.25 ),
    strip.text = element_text( size = 12 ),
    legend.text = element_text( size = 12 )
  )

ggsave( 'comparison_SouthAmerica.png', width = 8, height = 5 )


# error analysis
listSample <- mxF[ mxF$country_code > 900, ]$country_code %>% unique

# dados de mortalidade com repeticoes
for( i in mxF$country_code %>% unique ){
  lenF = length(mxF[ mxF$country_code == i, ]$age)
  lenM = length(mxM[ mxM$country_code == i, ]$age)

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

mxM[ mxM$country_code == 921, 1:5]

dat <-
  revSurvWpp( country_list = listSample , year = 2010, family = 'General' )

baselinedat <-
  tfr[tfr$name %in% listSample,
      c('name', '1990-1995', '1995-2000','2000-2005','2005-2010','2010-2015')] %>%
  as.data.table %>%
  melt( id.vars = c('name'),
        measure.vars = c('1990-1995', '1995-2000','2000-2005','2005-2010','2010-2015'),
        variable.name = 'year',
        value.name = 'TFR') %>%
  .[,.( Country = name, year = as.numeric( dic_year[as.character(year)] ), TFR.wpp = round(TFR,4) ) ]


datError <-
  merge(
    dat,
    baselinedat,
    by = c('Country', 'year')
    )

datError$relative_diff = round( ( datError$TFR.wpp - datError$TFR )/ datError$TFR.wpp, 4 )
plot(y = datError$relative_diff, x = datError$year )

dat <-
  rbind(
    revSurvWpp( country_list = 'Brazil', year = 1970, family = 'West' ) %>%
      as.data.table %>%
      .[,ref:=1970],
    revSurvWpp( country_list = 'Brazil', year = 1980, family = 'West' ) %>%
      as.data.table %>%
      .[,ref:=1980],
    revSurvWpp( country_list = 'Brazil', year = 1991, family = 'West' ) %>%
      as.data.table %>%
      .[,ref:=1991],
    revSurvWpp( country_list = 'Brazil', year = 2000, family = 'West' ) %>%
      as.data.table %>%
      .[,ref:=2000],
    revSurvWpp( country_list = 'Brazil', year = 2010, family = 'West' ) %>%
      as.data.table %>%
      .[,ref:=2010]
  )

ggplot( data = dat ) +
  geom_line( aes( x = as.numeric(year), y = TFR, color = as.factor(ref) ), size = 1 ) +
  scale_x_continuous( breaks = seq(1955,2020,2.5), limits = c( 1955, 2010 ), name = '' ) +
  scale_y_continuous( breaks = seq( 1, 7, 0.25 ), name = 'TFR' )  +
  scale_color_manual( values = c( '2010' = 'black', '2000' = 'gray20',
                                  '1991' = 'gray40', '1980' = 'gray60',
                                  '1970' = 'gray80' ), name = ''  ) +
  theme_classic() +
  theme(
    legend.position = 'top',
    axis.text.x = element_text( color = 'black', angle = 90, vjust = 0.5, size = 12 ),
    axis.text.y = element_text( color = 'black', size = 12 ),
    panel.grid.major = element_line( color = 'gray90', linetype = 'dashed', size = 0.25 ),
    strip.text = element_text( size = 12 ),
    legend.text = element_text( size = 12 )
  )

plot( y = datError$diff, x = datError$year )
dat$Type <- 'fertestr'
baselinedat$Type <- 'wpp2019'

?life.table

life.table( x = mxM[ mxM$country_code == 1830,]$age,
            nDx = mxM[ mxM$country_code == 1830, c( paste0( 2010, '-', 2015 ) ) ],
            nKx = rep(1,23) )
LifeTable(x = c(0, 1, 5, seq(10,100,5)), mx = unique(mxM[ mxM$country_code == 1830, c( paste0( 2010, '-', 2015 ) ) ]))
