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
    
module xyzinfo_m

   use precis_m
   implicit none
   private
!  Module contains all variables set by setxyz and setaxu. It combines 
!  the xyzinfo.h, map.h. bigxyz4.h vecsuv.h and vecsuva.h include files.

!  xyzinfo
   real(kind=rx), public, dimension(:), allocatable :: x, y, z, rlat, rlong, &
                                                       wts

!  bigxy4
   real(kind=rx), public, dimension(:,:), allocatable :: xx4, yy4

!  map
   real(kind=rx), public, dimension(:), pointer :: emu, emv
   real(kind=rx), public, dimension(:), allocatable :: em,                   &
                                                      f, fu, fv,             &
                                                      dmdx, dmdy, dmdxv, dmdyu 

!  vecsuv
   real(kind=rx), public, dimension(:), pointer :: ax, ay, az, bx, by, bz

!  vecsuva
   real(kind=rx), public, dimension(:), allocatable :: axu, bxv, ayu, byv,   &
                                                       azu, bzv

   real(kind=rx), public, save    :: ds

end module xyzinfo_m

