load( 'tests/cohortcomp_vr_testdata.RData' )

cohparitycomp_vr( asfr_chl_vr_1971_1991, parity_chl_1992 )

# for brazil we need to add more data (we need 15-20 preceding years of VR)

asfr_temp =
  asfr_bra_states_2000_2010[ asfr_bra_states_2000_2010$statecode == 31 &
                               asfr_bra_states_2000_2010$age_start %in% 15:39 &
                               asfr_bra_states_2000_2010$year %in% 2000:2009, ]


parity_temp =
  parity_bra_states_2010[ parity_bra_states_2010$statecode == 31 &
                            parity_bra_states_2010$age_start %in% 15:39 , ]

cohparitycomp_vr( asfr_temp, parity_temp )
