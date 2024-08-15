
# China 2017 NFS
chn2017_owcm
lt_chn # life table with year info

estim_owcm_chn2017 = own_children( owc = chn2017_owcm, owc_lt = lt_chn )
tapply( estim_owcm_chn2017$asfr, estim_owcm_chn2017$ref_year, sum )

# Colombia 1978 - 1
col1978_owcm
col1978_owcm_lt # life table without year info (West model life table)

estim_owcm_col1978 = own_children( owc = col1978_owcm, owc_lt = col1978_owcm_lt )
tapply( estim_owcm_col1978$asfr, estim_owcm_col1978$ref_year, sum )

# Colombia 1978 - 2

col1978_owcm
lt_col  # life table with year info
estim_owcm_col1978 = own_children( owc = col1978_owcm, owc_lt = lt_col )
tapply( estim_owcm_col1978$asfr, estim_owcm_col1978$ref_year, sum )



