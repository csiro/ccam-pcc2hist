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
   real, pointer, contiguous, dimension(:,:), public :: xg, yg
   integer, pointer, contiguous, dimension(:,:), public :: nface
   real, pointer, contiguous, dimension(:,:,:), public :: xyg
   integer :: interp_win(2)
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
   use netcdf, only : NF90_FILL_FLOAT

   implicit none
   real, intent(in), dimension(:,:)  ::  s_in    ! (il,jl)
   real, intent(out), dimension(:,:) :: array    ! (nxhis,nyhis)
   integer, intent(in) :: int_type

   real, dimension(ifull) :: s
   real, dimension(:,:,:), allocatable, save :: sx
      
   real, dimension(4) :: r
   integer :: i, j, n, nn
   integer :: idel, jdel
   real :: xxg, yyg, c1, c2, c3, c4, aaa
   integer :: nxhis, nyhis

   if ( int_type == int_none ) then
      array = s_in
      return
   end if
   
   nxhis = size(array,1)
   nyhis = size(array,2)

!  This routine is called with a 2D array, (il,jl). Rather than do a reshape
!  in every call do it in the routine

!  For performance on cherax, do the reshape explicitly
!   s = reshape ( s_in, (/ ifull /) )
   do j=1,jl
      do i=1,il
         s ( i + (j-1)*il ) = s_in(i,j)
      end do
   end do

   if (allocated(sx)) then
     if (size(sx,1)/=il+4.or.size(sx,2)/=il+4.or.size(sx,3)/=npanels+1) then
       deallocate(sx)
     end if
   end if
   if (.not.allocated(sx)) allocate ( sx(-1:il+2,-1:il+2,0:npanels) )

!  This is intsb           EW interps done first
!  First extend s arrays into sx - this one -1:il+2 & -1:il+2
   do n=0,npanels
      do j=1,il
         do i=1,il
            sx(i,j,n)=s(ind(i,j,n))
         enddo  ! i loop
         sx(0,j,n)=s(i_w(ind(1,j,n)))
         sx(-1,j,n)=s(i_ww(ind(1,j,n)))
         sx(il+1,j,n)=s(i_e(ind(il,j,n)))
         sx(il+2,j,n)=s(i_ee(ind(il,j,n)))
      enddo! j loop
      do i=1,il
         sx(i,0,n)=s(i_s(ind(i,1,n)))
         sx(i,-1,n)=s(i_ss(ind(i,1,n)))
         sx(i,il+1,n)=s(i_n(ind(i,il,n)))
         sx(i,il+2,n)=s(i_nn(ind(i,il,n)))
     enddo ! i loop
!    for ew interpolation, sometimes need (different from ns):
!          (-1,0),   (0,0),   (0,-1)   (-1,il+1),   (0,il+1),   (0,il+2)
!         (il+1,0),(il+2,0),(il+1,-1) (il+1,il+1),(il+2,il+1),(il+1,il+2)
     sx(-1,0,n)=s(lwws(n))
     sx(0,0,n)=s(lws(n))
     sx(0,-1,n)=s(lwss(n))
     sx(il+1,0,n)=s(les(n))
     sx(il+2,0,n)=s(lees(n))
     sx(il+1,-1,n)=s(less(n))
     sx(-1,il+1,n)=s(lwwn(n))
     sx(0,il+2,n)=s(lwnn(n))
     sx(il+2,il+1,n)=s(leen(n))
     sx(il+1,il+2,n)=s(lenn(n))
     sx(0,il+1,n)   =s(i_wn(ind(1,il,n)))
     sx(il+1,il+1,n)=s(i_en(ind(il,il,n)))
  enddo ! n loop

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
           do nn=2,3       ! N.B.
              c1=sx(idel-1,jdel+nn-2,n)
              c2=sx(idel  ,jdel+nn-2,n)
              c3=sx(idel+1,jdel+nn-2,n)
              c4=sx(idel+2,jdel+nn-2,n)
              r(nn) = ((1.-xxg)*((2.-xxg)*((1.+xxg)*c2-xxg*c1/3.)           &
                   -xxg*(1.+xxg)*c4/3.) +xxg*(1.+xxg)*(2.-xxg)*c3)/2.
           enddo ! nn loop
!       r       ={(1-x     )*{(2-x     )*[(1+x     )*c2-x     *c1/3]
!         -x     *(1+x     )*c4/3}
!         +x    *(1+x     )*(2-x     )*c3}/2
           do nn=1,4,3       ! N.B.
              c2=sx(idel  ,jdel+nn-2,n)
              c3=sx(idel+1,jdel+nn-2,n)
              r(nn)=(1.-xxg)*c2 +xxg*c3
           enddo! nn loop
!       array(i,j)=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)
!    .             -yyg*(1.+yyg)*r(4)/3.)
!    .             +yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
!      following does Bermejo Staniforth
           aaa=((1.-yyg)*((2.-yyg)*((1.+yyg)*r(2)-yyg*r(1)/3.)                &
                -yyg*(1.+yyg)*r(4)/3.) + yyg*(1.+yyg)*(2.-yyg)*r(3))/2.
           aaa=min( aaa , max( sx(idel,jdel,n),sx(idel+1,jdel,n),             &
                            sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) ) 
           array(i,j)=max( aaa , min( sx(idel,jdel,n),sx(idel+1,jdel,n),      &
                                   sx(idel,jdel+1,n),sx(idel+1,jdel+1,n) ) )
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
end subroutine ints

end module interp_m
