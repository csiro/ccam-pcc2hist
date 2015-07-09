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
    
module ind_m
   implicit none
   public :: ind, indv
contains
   function ind (i, j, n) result(iq)
     use newmpar_m, only : il  ! This is required when it's not an internal fn.
     integer, intent(in) :: i, j, n
     integer :: iq

     iq = i + (j - 1)*il + n*il*il 

   end function ind
   subroutine indv(iq, i, j, n) 
      use newmpar_m, only : il
      integer , intent(in) :: iq 
      integer , intent(out) :: i 
      integer , intent(out) :: j 
      integer , intent(out) :: n 

      n = (iq - 1)/(il*il) 
      j = 1 + (iq - n*il*il - 1)/il 
      i = iq - (j - 1)*il - n*il*il 

   end subroutine indv

end module ind_m
 
