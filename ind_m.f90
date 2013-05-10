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
 
