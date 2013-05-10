module indices_m
   implicit none
   integer, public, dimension(:), allocatable :: i_n, i_s, i_w, i_e,         &
                                                 i_nn, i_ss, i_ww, i_ee,     &
                                                 i_ne, i_se, i_en, i_wn,     &
                                                 i_wu, i_sv, i_wu2, i_sv2,   &
                                                 i_eu2, i_nv2, i_ev2, i_nu2, &
                                                 i_eu, i_nv 

   integer, public, dimension(:), allocatable :: lwws, lws, lwss, les, lees, &
                                                 less, lwwn, lwnn, leen,     &
                                                 lenn, lsww, lsw, lssw,      &
                                                 lsee, lsse, lnww, lnw,      &
                                                 lnnw, lnee, lnne 

end module indices_m
