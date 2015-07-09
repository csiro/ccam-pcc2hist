! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
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
