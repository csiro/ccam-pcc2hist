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
    
module kinds_m
   ! Define kinds for single and double precision. Note that on the NEC
   ! single precision may be 64 bit. These parameters are only for 
   ! defining overloadable functions, not for selecting particular values.
   implicit none
   ! Default real kind
   integer, parameter :: sp = selected_real_kind(6)
!   real(kind=sp), parameter :: xxx = 1.0
!  Use default real here so that it works with Intel -r8 -i8
   real, parameter :: xxx = 1.0
   ! This specification works with standard IEEE and with NEC -ew 
   integer, parameter :: dp = selected_real_kind(precision(xxx)+1)
end module kinds_m

