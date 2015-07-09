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

module jim_utils
   implicit none
   public :: tay
contains

!--------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE  TAY
!  Evaluate the complex function W of Z whose real
!  Taylor series coefficients are RA.
!
!  --> Z    function argument (complex)
!  --> RA   Taylor coefficients (real)
!  --> N    number of coefficients (starting with the linear term)
!  <-- W    Taylor-series approximation of the function (complex)
!--------------------------------------------------------------------------
subroutine tay(z, ra, n, w)
   use precis_m
   integer, intent(in) :: n
   complex(kind=rx), intent(in) :: z
   complex(kind=rx), intent(out) :: w
   real(kind=rx), intent(in), dimension(:) :: ra

   integer :: i

   w = 0.0
   do i = n, 1, -1
      w = (w + ra(i))*z
   end do

end subroutine tay
end module jim_utils
