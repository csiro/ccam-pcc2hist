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
