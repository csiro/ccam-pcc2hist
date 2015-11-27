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
    
module interp_m

   implicit none

#ifdef parallel_int
   real, pointer, dimension(:,:), public :: xg, yg
   integer, pointer, dimension(:,:), public :: nface
   real, pointer, dimension(:,:,:), public :: xyg
   integer :: xyg_win, nface_win
#else
   real, allocatable, dimension(:,:), public :: xg, yg
   integer, allocatable, dimension(:,:), public :: nface
#endif

!  Type of interpolation
!  Note that parameter int_default = 0 is defined in history module
   integer, parameter, public :: int_normal=0, int_nearest=1, int_lin=2
!  No interplation at all (output on CC grid)
   integer, parameter, public :: int_none=5

contains 

subroutine ints ( s_in, array, int_type )  ! input array (twice), output array

!     can put sx into work array (only 2d)
!     later may wish to save idel etc between array calls
!     this one does linear interp in x on outer y sides
!     doing x-interpolation before y-interpolation

   use newmpar_m
   use indices_m
   use ind_m
#ifdef usenc3
   use netcdf_m, only : NF90_FILL_FLOAT
#else
   use netcdf, only : NF90_FILL_FLOAT
#endif
   use logging_m

   implicit none
   real, intent(in), dimension(:,:)  ::  s_in    ! (il,jl)
   real, intent(out), dimension(:,:) :: array    ! (nxhis,nyhis)
   integer, intent(in) :: int_type

   real, dimension(:,:,:), allocatable, save :: sx
      
   real, dimension(2) :: r
   real, dimension(2:3) :: dmul
   real, dimension(1:4) :: cmul, emul, rmul
   integer :: i, j, n, nn
   integer :: idel, jdel
   integer :: n_n, n_e, n_w, n_s
   real :: xxg, yyg
   real :: cmin, cmax
   integer :: nxhis, nyhis

   if ( int_type == int_none ) then
      array = s_in
      return
   end if
   
   call START_LOG(ints_begin)
   
   nxhis = size(array,1)
   nyhis = size(array,2)

!  This routine is called with a 2D array, (il,jl). Rather than do a reshape
!  in every call do it in the routine

   if (allocated(sx)) then
     if (size(sx,1)/=il+4.or.size(sx,2)/=il+4.or.size(sx,3)/=npanels+1) then
       deallocate(sx)
     end if
   end if
   if (.not.allocated(sx)) allocate ( sx(-1:il+2,-1:il+2,0:npanels) )

!  This is intsb           EW interps done first
!  First extend s arrays into sx - this one -1:il+2 & -1:il+2

   !     first extend s arrays into sx - this one -1:il+2 & -1:il+2
   do n = 0,npanels,2
      sx(1:il,1:il,n) = reshape( s_in(1:il,1+n*il:il+n*il), (/ il, il /) )
      n_w = mod(n+5,6)*il
      n_e = mod(n+2,6)*il
      n_n = mod(n+1,6)*il
      n_s = mod(n+4,6)*il
      do i = 1,il
         sx(0,i,n)    = s_in(il,i+n_w)
         sx(-1,i,n)   = s_in(il-1,i+n_w)
         sx(il+1,i,n) = s_in(il+1-i,1+n_e)
         sx(il+2,i,n) = s_in(il+1-i,2+n_e)
         sx(i,il+1,n) = s_in(i,1+n_n)
         sx(i,il+2,n) = s_in(i,2+n_n)
         sx(i,0,n)    = s_in(il,il+1-i+n_s)
         sx(i,-1,n)   = s_in(il-1,il+1-i+n_s)
      end do
      sx(-1,0,n)      = s_in(il,2+n_w)         ! wws
      sx(0,-1,n)      = s_in(il,il-1+n_s)      ! wss
      sx(0,0,n)       = s_in(il,1+n_w)         ! ws
      sx(il+1,0,n)    = s_in(il,1+n_e)         ! es  
      sx(il+2,0,n)    = s_in(il-1,1+n_e)       ! ees 
      sx(-1,il+1,n)   = s_in(il,il-1+n_w)      ! wwn
      sx(0,il+2,n)    = s_in(il-1,il+n_w)      ! wnn
      sx(il+2,il+1,n) = s_in(2,1+n_e)          ! een  
      sx(il+1,il+2,n) = s_in(1,2+n_e)          ! enn  
      sx(0,il+1,n)    = s_in(il,il+n_w)        ! wn  
      sx(il+1,il+1,n) = s_in(1,1+n_e)          ! en  
      sx(il+1,-1,n)   = s_in(il,2+n_e)         ! ess  
   end do  ! n loop
   do n = 1,npanels,2
      sx(1:il,1:il,n) = reshape( s_in(1:il,1+n*il:il+n*il), (/ il, il /) )
      n_w = mod(n+4,6)*il
      n_e = mod(n+1,6)*il
      n_n = mod(n+2,6)*il
      n_s = mod(n+5,6)*il
      do i = 1,il
         sx(0,i,n)    = s_in(il+1-i,il+n_w)
         sx(-1,i,n)   = s_in(il+1-i,il-1+n_w)
         sx(il+1,i,n) = s_in(1,i+n_e)
         sx(il+2,i,n) = s_in(2,i+n_e)
         sx(i,il+1,n) = s_in(1,il+1-i+n_n)
         sx(i,il+2,n) = s_in(2,il+1-i+n_n)
         sx(i,0,n)    = s_in(i,il+n_s)
         sx(i,-1,n)   = s_in(i,il-1+n_s)
      end do
      sx(-1,0,n)      = s_in(il-1,il+n_w)     ! wws
      sx(0,-1,n)      = s_in(2,il+n_s)        ! wss
      sx(0,0,n)       = s_in(il,il+n_w)       ! ws
      sx(il+1,0,n)    = s_in(1,1+n_e)         ! es
      sx(il+2,0,n)    = s_in(1,2+n_e)         ! ees
      sx(-1,il+1,n)   = s_in(2,il+n_w)        ! wwn   
      sx(0,il+2,n)    = s_in(1,il-1+n_w)      ! wnn  
      sx(il+2,il+1,n) = s_in(1,il-1+n_e)      ! een  
      sx(il+1,il+2,n) = s_in(2,il+n_e)        ! enn  
      sx(0,il+1,n)    = s_in(1,il+n_w)        ! wn  
      sx(il+1,il+1,n) = s_in(1,il+n_e)        ! en  
      sx(il+1,-1,n)   = s_in(2,1+n_e)         ! ess  
   end do  ! n loop      

  select case  ( int_type )
  case ( int_normal )
     do j=1,nyhis
        do i=1,nxhis
           n=nface(i,j)
           idel=int(xg(i,j))
           xxg=xg(i,j)-idel
!       yg here goes from .5 to il +.5
           jdel=int(yg(i,j))
           yyg=yg(i,j)-jdel
           cmul(1) = (1.-xxg)*(2.-xxg)*(-xxg)/6.
           cmul(2) = (1.-xxg)*(2.-xxg)*(1.+xxg)/2.
           cmul(3) = xxg*(1.+xxg)*(2.-xxg)/2.
           cmul(4) = (1.-xxg)*(-xxg)*(1.+xxg)/6.
           dmul(2) = (1.-xxg)
           dmul(3) = xxg
           emul(1) = (1.-yyg)*(2.-yyg)*(-yyg)/6.
           emul(2) = (1.-yyg)*(2.-yyg)*(1.+yyg)/2.
           emul(3) = yyg*(1.+yyg)*(2.-yyg)/2.
           emul(4) = (1.-yyg)*(-yyg)*(1.+yyg)/6.
           cmin = minval(sx(idel:idel+1,jdel:jdel+1,n))
           cmax = maxval(sx(idel:idel+1,jdel:jdel+1,n))  
           rmul(1) = sum(sx(idel:idel+1,  jdel-1,n)*dmul(2:3))
           rmul(2) = sum(sx(idel-1:idel+2,jdel,  n)*cmul(1:4))
           rmul(3) = sum(sx(idel-1:idel+2,jdel+1,n)*cmul(1:4))
           rmul(4) = sum(sx(idel:idel+1,  jdel+2,n)*dmul(2:3))
           array(i,j) = min(max(cmin,sum(rmul(1:4)*emul(1:4))),cmax) ! Bermejo & Staniforth
        enddo ! i loop
     enddo ! j loop
  case (int_nearest)
     ! This case handles the missing values automatically.
     do j=1,nyhis
        do i=1,nxhis
           n=nface(i,j)
           idel=nint(xg(i,j))
           jdel=nint(yg(i,j))
           array(i,j) = sx(idel,jdel,n)
        end do
     end do
  case (int_lin)
     ! Bi-linear. This will probably only be used with vextrap=missing.
     do j=1,nyhis
        do i=1,nxhis
           n = nface(i,j)
           idel = floor(xg(i,j))
           jdel = floor(yg(i,j))
           ! Set missing if any are missing
           if ( sx(idel,jdel,n) == NF90_FILL_FLOAT .or.  &
                sx(idel+1,jdel,n) == NF90_FILL_FLOAT .or. &
                sx(idel+1,jdel+1,n) == NF90_FILL_FLOAT .or. &
                sx(idel,jdel+1,n) == NF90_FILL_FLOAT ) then
              array(i,j) = NF90_FILL_FLOAT
           else
              xxg = xg(i,j)-idel
              yyg = yg(i,j)-jdel
              r(1) = sx(idel,jdel,n) + (sx(idel+1,jdel,n)-sx(idel,jdel,n))*xxg
              r(2) = sx(idel,jdel+1,n) + (sx(idel+1,jdel+1,n)-sx(idel,jdel+1,n))*xxg
              array(i,j) = r(1) + (r(2)-r(1))*yyg
           end if
        end do
     end do
  case default
     print*, "Error, unexpected interpolation option", int_type
  end select

  !deallocate ( sx )
  
  call END_LOG(ints_end)
  
end subroutine ints

end module interp_m
