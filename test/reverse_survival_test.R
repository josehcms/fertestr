devtools::install_github('josehcms/fertestr')
require(fertestr)

### Examples:

# 1 - User input data
pop_c <-  c( 281260, 261320, 268410, 286810, 278990, 293760, 293490, 302060, 315970, 267190, 326980, 280260, 354120, 356920, 354830 )
pop_w <- c(  815930, 780320, 697160, 626430, 361650, 435880, 393760, 352520, 294280, 230200, 160590, NA )
lx_c <- c( 1.0000, 0.9320, 0.9275, 0.9228, 0.9165, 0.9125, 0.9110, 0.9094, 0.9079, 0.9063, 0.9048, 0.9032, 0.9017, 0.9001, 0.8986, 0.8970 )
lx_w <- c( 0.91381, 0.90989, 0.90492, 0.89798, 0.88893, 0.87596, 0.86029, 0.84188, 0.81791, 0.78472, 0.73735, 0.67316 )
q0_5 <-  c( 0.0683, 0.1008, 0.1189)
q15_45f <- c( 0.1946, 0.2290, 0.2674 )
asfr <- c( 0.0000, 0.0418,0.1535, 0.1482, 0.1118, 0.0708, 0.0301, 0.0032 )
asfr_15prior <- c( 0.0000, 0.0533, 0.1974, 0.2144, 0.1836, 0.1332, 0.0676, 0.0134 )

FertRevSurv( ages1_c = 0:14, popx1_c = pop_c,
             ages5_w = seq( 10, 65, 5 ), popx5_w = pop_w,
             lx1_c = lx_c, lx5_w = lx_w,
             asfr5 = asfr,
             asfr5_15prior = asfr_15prior,
             q0_5 = q0_5, q15_45f = q15_45f,
             date_ref = '2008-03-03' )

# 2 - using log-quadratic estimated survival functions
ltb <- SingleAgeLogQuad( q15_45 = 0.20, q0_5 = 0.07 )
ltf <- SingleAgeLogQuad( q15_45 = 0.19, q0_5 = 0.05, sex = 'female' )

lx_w <- ltf[ ltf$x %in% seq( 10, 65, 5 ), ]$lx
lx_c <- ltb[ ltb$x %in% seq( 0, 15 ), ]$lx

FertRevSurv( ages1_c = 0:14, popx1_c = pop_c,
             ages5_w = seq( 10, 65, 5 ), popx5_w = pop_w,
             lx1_c = lx_c, lx5_w = lx_w,
             asfr5 = asfr,
             asfr5_15prior = asfr_15prior,
             q0_5 = q0_5, q15_45f = q15_45f,
             date_ref = '2008-03-03' )

# 3 - using model life table estimated survival functions
lx_w <- find_mlt( lt_family = 'General', e0 = 69,
                  ages = seq( 10, 65, 5 ), sex = 'female' )$lx_std
lx_c <- find_mlt( lt_family = 'General', e0 = 67,
                  ages = seq( 0, 15 ), sex = 'both' )$lx_std

FertRevSurv( ages1_c = 0:14, popx1_c = pop_c,
             ages5_w = seq( 10, 65, 5 ), popx5_w = pop_w,
             lx1_c = lx_c, lx5_w = lx_w,
             asfr5 = asfr,
             asfr5_15prior = asfr_15prior,
             q0_5 = q0_5, q15_45f = q15_45f,
             date_ref = '2008-03-03' )

# 4 - country estimation by WPP 2019 data

FertRevSurvWPP( locations = c('Argentina','Italy','Canada','Brazil','Colombia'),
                date_ref = '2015-07-31',
                lt_family = 'West' )

# 5 - using WPP information with user data - survival functions

locs_avail() # use lower-middle-income countries' data

ltb <- FetchLifeTableWpp2019( locations = 1501,
                              year = decimal_anydate('2008-03-03'),
                              sex = 'both' )

ltf <- FetchLifeTableWpp2019( locations = 1501,
                              year = decimal_anydate('2008-03-03'),
                              sex = 'female' )

# use family life tables to retrieve single age lx estimates
lx_c <- find_mlt( lt_family = 'West', e0 = ltb$ex[1], ages = 0:15, sex = 'both' )$lx_std
lx_w <- ltf$lx[ ltf$x %in% seq( 10, 65, 5)]

FertRevSurv( ages1_c = 0:14, popx1_c = pop_c,
             ages5_w = seq( 10, 65, 5 ), popx5_w = pop_w,
             lx1_c = lx_c, lx5_w = lx_w,
             asfr5 = asfr,
             asfr5_15prior = asfr_15prior,
             q0_5 = q0_5, q15_45f = q15_45f,
             date_ref = '2008-03-03' )

# 6 - using WPP information with user data - asfr

locs_avail() # use lower-middle-income countries' data

asfr <- FetchFertilityWpp2019( locations = 'Cambodia',
                               year = decimal_anydate( '2008-03-03' ) )$asfr

asfr_15prior <- FetchFertilityWpp2019( locations = 'Cambodia',
                                       year = decimal_anydate( '1993-03-03' ) )$asfr

FertRevSurv( ages1_c = 0:14, popx1_c = pop_c,
             ages5_w = seq( 10, 65, 5 ), popx5_w = pop_w,
             lx1_c = lx_c, lx5_w = lx_w,
             asfr5 = asfr,
             asfr5_15prior = asfr_15prior,
             q0_5 = q0_5, q15_45f = q15_45f,
             date_ref = '2008-03-03' )

# 7 - using WPP information with user data - qx

locs_avail() # use lower-middle-income countries' data

q0_5 <- q_calcWpp2019( location = 'Cambodia',
                       years = decimal_anydate( '2008-03-03' ) - c( 2.5, 7.5, 12.5),
                       sex = 'both', age_inf = 0, age_sup = 5 )$qx

q15_45f <- q_calcWpp2019( location = 'Cambodia',
                          years = decimal_anydate( '2008-03-03' ) - c( 2.5, 7.5, 12.5),
                          sex = 'female', age_inf = 15, age_sup = 60 )$qx

FertRevSurv( ages1_c = 0:14, popx1_c = pop_c,
             ages5_w = seq( 10, 65, 5 ), popx5_w = pop_w,
             lx1_c = lx_c, lx5_w = lx_w,
             asfr5 = asfr,
             asfr5_15prior = asfr_15prior,
             q0_5 = q0_5, q15_45f = q15_45f,
             date_ref = '2008-03-03' )
