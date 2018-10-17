! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module hyblevs_m
   implicit none
   real, save, dimension(:), allocatable :: sig, sigh, anf, bnf, anh, bnh, &
                                            dadnf, dbdnf

   ! For interpolation to pressure levels

   ! Controls whether vertical extrapolation below surface is done
   integer, parameter :: vextrap_default=0, vextrap_lin=1, vextrap_none=2, &
                         vextrap_missing=3, vextrap_t=4
   ! This calculation needs 64 bit precision.
   integer, private, parameter :: r8 = selected_real_kind(p=15, r=307)

end module hyblevs_m
   

