#' Adjust children counts from unidentified mothers
#'
#' Proportionately redistribute children from unidentified mothers (women age index equal to 99)
#'
#' @param owc a data.frame with 5 columns:
#' `source_year` - reference year of data source/data collection;
#' `ages_w`   - women's ages at the interview year;
#' `ages_c`   - chilren's ages at the interview year;
#' `counts_w` - number of women in the age group `ages_w`
#' `counts_c` - number of children in the age group `ages_c` from mothers in the age group `ages_w`
#'
#' @return a data.frame similar to the input data.frame, but with 2 additional columns `k` and `adjcounts_c`:
#' `source_year` - reference year of data source/data collection;
#' `ages_w`   - women's ages at the interview year;
#' `ages_c`   - children's ages at the interview year;
#' `counts_w` - number of women in the age group `ages_w`
#' `counts_c` - number of children in the age group `ages_c` from mothers in the age group `ages_w`
#' `k`        - factor applied for adjusting the number of children
#' `adjcounts_c` - adjusted number of children in the age group `ages_c` from mothers in the age group `ages_w` after proportionately redistributing children from unidentified mothers
#'
#' @keywords internal
#'
#'
owcm_adjust_unmatch =
  function( owc ){

    # check if there are no unmatched children
    flag_match = nrow( owc[ owc$ages_w == 99, ] ) == 0

    if( flag_match ){
      owc_adj = owc[ owc$ages_w != 99, ]
      owc_adj$k = 1
      owc_adj$adjcounts_c = owc_adj$counts_c * owc_adj$k

    } else{
      # total children by age (excluding missing)
      total_c = tapply( owc[ owc$ages_w != 99, 'counts_c'  ], owc[ owc$ages_w != 99, 'ages_c'  ], sum )

      # children with unidentified mom
      unident_c = tapply( owc[ owc$ages_w == 99, 'counts_c'  ], owc[ owc$ages_w == 99, 'ages_c'  ], sum )

      # calculate adjust factor k
      factor_k =
        data.frame(
          ages_c = as.numeric( names( total_c ) ),
          k = 1 + unident_c / total_c
        )

      # adjust children counts
      owc_adj =
        merge(
          owc[ owc$ages_w != 99, ],
          factor_k,
          by = 'ages_c'
        )

      owc_adj$adjcounts_c = owc_adj$counts_c * owc_adj$k
    }

    owc_adj = owc_adj[ , c( 'source_year', 'ages_w', 'ages_c', 'counts_w', 'counts_c', 'k', 'adjcounts_c' ) ]
    owc_adj = owc_adj[ order( owc_adj$source_year, owc_adj$ages_w, owc_adj$ages_c ), ]

    return( owc_adj )
  }

#' Lag/shift function using base R
#'
#' @param x vector
#' @return x shifted up by 1 with NA replacing the last element
#'
#' @keywords internal
#'
lagfunc = function( x ){
  n = length( x )
  c( x[ 2:n ], NA )
}


#' Own-children method for fertility estimation
#'
#' Estimate age-specific fertility rates (single age) using the own-children method
#'
#' @param owc a data.frame with the information on matched children and mothers by age of children and age of mother with 5 columns:
#' `source_year` - reference year of data source/data collection;
#' `ages_w`   - women's ages at the interview year;
#' `ages_c`   - children's ages at the interview year;
#' `counts_w` - number of women in the age group `ages_w`
#' `counts_c` - number of children in the age group `ages_c` from mothers in the age group `ages_w`
#'
#' @param owc_lt a data.frame with the information on survivorship of children and adults. Requires these 3 columns:
#' `age` - age;
#' `Lx_b`   - probability of surviving to age group from x to x+1, both sexes;
#' `Lx_f`   - probability of surviving to age group from x to x+1, females;
#' `year`   - (optional), if there are life tables from multiple years

#' @return a data.frame with 5 columns:
#' `ref_year` - reference year of fertility estimation;
#' `age`      - women's age group;
#' `births`   - number of births from women at each age group;
#' `women`    - number of women in the age group
#' `asfr`     - age-specific fertility rates
#'
#' @export

own_children =
  function( owc, owc_lt ){

    ## 1) redistribution of children with unidentified mother -> proportionately redistribute Children by Child_Age
    owc_adj = owcm_adjust_unmatch( owc )

    # add age when became mother
    owc_adj$ages_m = owc_adj$ages_w - owc_adj$ages_c
    # add reference year of estimates
    owc_adj$ref_year = owc_adj$source_year - owc_adj$ages_c
    owc_adj$year = trunc(owc_adj$ref_year)

    ## 2) estimation of survivorship probabilities for children
    owc_lt$Lx_c = owc_lt$Lx_b # Lx children = Lx both sexes

    if( is.null( owc_lt$year ) ){
      owc_adj =
        merge(
          owc_adj,
          owc_lt[ , c( 'age', 'Lx_c' ) ],
          by.x = 'ages_c',
          by.y = 'age',
          all.x = TRUE,
          all.y = FALSE
        )
    } else{
      owc_adj =
        merge(
          owc_adj,
          owc_lt[ , c( 'year', 'age', 'Lx_c' ) ],
          by.x = c( 'year', 'ages_c' ),
          by.y = c( 'year', 'age' ),
          all.x = TRUE,
          all.y = FALSE
        )
    }


    ## 3) estimation of survivorship probabilities for adult females

    owc_lt$Lx_w = owc_lt$Lx_f # Lx women = Lx adult females
    owc_lt$Lx_m = owc_lt$Lx_f # Lx mothers = Lx adult females

    if( is.null( owc_lt$year ) ){
      # merging information from age at survey
      owc_adj =
        merge(
          owc_adj,
          owc_lt[ , c( 'age', 'Lx_w' ) ],
          by.x = 'ages_w',
          by.y = 'age',
          all.x = TRUE,
          all.y = FALSE
        )

      # merging information from age at motherhood
      owc_adj =
        merge(
          owc_adj,
          owc_lt[ , c( 'age', 'Lx_m' ) ],
          by.x = 'ages_m',
          by.y = 'age',
          all.x = TRUE,
          all.y = FALSE
        )
    } else{
      # merging information from age at survey
      owc_adj =
        merge(
          owc_adj,
          owc_lt[ , c( 'year', 'age', 'Lx_w' ) ],
          by.x = c( 'year', 'ages_w' ),
          by.y = c( 'year', 'age' ),
          all.x = TRUE,
          all.y = FALSE
        )

      # merging information from age at motherhood
      owc_adj =
        merge(
          owc_adj,
          owc_lt[ , c( 'year', 'age', 'Lx_m' ) ],
          by.x = c( 'year', 'ages_m' ),
          by.y = c( 'year', 'age' ),
          all.x = TRUE,
          all.y = FALSE
        )
    }

    # order
    owc_adj = owc_adj[ order( owc_adj$year, owc_adj$ages_m, owc_adj$ages_w, owc_adj$ages_c ), ]

    ## 4) reverse survival of children to estimate back births
    owc_adj$B_temp = owc_adj$adjcounts_c / owc_adj$Lx_c

    # redistribute births according to fraction of the year
    owc_adj$B_temp_lead = NA
    for( this_year in unique( owc_adj$ref_year ) ){
      owc_adj[ owc_adj$ref_year == this_year, ]$B_temp_lead = lagfunc( owc_adj[ owc_adj$ref_year == this_year, ]$B_temp )
    }

    owc_adj$split_factor = owc_adj$source_year - trunc( owc_adj$source_year )
    owc_adj$births = owc_adj$B_temp * owc_adj$split_factor  + ( 1 - owc_adj$split_factor  ) * owc_adj$B_temp_lead

    ## 5) reverse survival of adult females
    owc_adj$W_temp = owc_adj$counts_w * ( owc_adj$Lx_m / owc_adj$Lx_w )

    # redistribute women person-years according to fraction of the year
    owc_adj$W_temp_lead = NA
    for( this_year in unique( owc_adj$ref_year ) ){
      owc_adj[ owc_adj$ref_year == this_year, ]$W_temp_lead = lagfunc( owc_adj[ owc_adj$ref_year == this_year, ]$W_temp )
    }

    owc_adj$women = owc_adj$W_temp * owc_adj$split_factor  + ( 1 - owc_adj$split_factor  ) * owc_adj$W_temp_lead

    ## 6) calculation of ASFRs
    owc_adj$asfr = owc_adj$births / owc_adj$women

    ## 7) prepare output data.frame
    res = owc_adj[ owc_adj$ages_m %in% 15:54, c( 'ref_year', 'ages_m', 'births', 'women', 'asfr' ) ]
    res = res[ order( res$ref_year, res$ages_m ), ]
    res = data.frame( apply( res, 2, function(x) replace( x, is.na( x ), 0 ) ) )

    return( res )

  }
