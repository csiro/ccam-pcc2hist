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

