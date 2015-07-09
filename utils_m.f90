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
    
module utils_m
   use kinds_m
   implicit none
   private
   interface str2array
      module procedure str2rarray, str2iarray, str2carray
   end interface
   private  str2rarray, str2iarray, str2carray
   public :: get_new_unit, str2array, assert, lcase
   public :: search_fle
   private :: isrchfle
   interface search_fle
      module procedure isrchfle
   end interface

   public :: search_fgt
   private :: isrchfgt, vsrchfgt
   interface search_fgt
      module procedure isrchfgt, vsrchfgt, vsrchfgt_dp, vxsrchfgt, vxsrchfgt_dp
   end interface

   public :: fpequal
   private :: fpequal_sp, fpequal_dp
   interface fpequal
      module procedure fpequal_sp, fpequal_dp
   end interface

contains

   subroutine get_new_unit ( unum )

!     Get an unused unit number

      integer, intent(out) :: unum
!     nmax is system depedndent, but this should be big enough for practical
!     use
      integer, parameter  :: nmax=100

      integer, parameter :: startno = 7
      integer :: i, ios
      logical :: used, ex

      unum = -1  ! Error return value
      do i = startno, nmax
         inquire(unit=i, iostat = ios, exist = ex, opened = used )
	 if ( ios /= 0 ) then
	    print*, " Error from inquire in get_new_unit "
	    print*, " Last test, unit ", i, " IOSTAT = ", ios
	    stop
	 end if
	 if (ex .and. .not. used ) then
	    unum = i
	    exit
	 end if
      end do
   end subroutine get_new_unit

   subroutine assert (ok, msg)
!  If not ok then display message & stop

      logical, intent(in) :: ok
      character(len=*), intent(in) :: msg
      
      integer, parameter :: lu = 0 ! std. error output unit
      
      if ( .not. ok ) then
         if ( len_trim(msg) == 0 ) then
            write ( lu, '(a)' ) 'assert called with blank message'
         else
            write ( lu, '(a)' )  trim(msg)
         end if
         stop
      end if

   end subroutine assert

   subroutine str2rarray(optarg,arr,n)
      ! Convert a string of values delimited by spaces or commas into 
      ! separate real values.
      character(len=*), intent(in) :: optarg
      real, dimension(:), intent(inout) :: arr
      integer, intent(inout) :: n
      integer :: nvals, i, nmax
      logical :: inword, space, startword

!     Go through the string counting the number of distinct words (delimited
!     by commas or spaces
      inword = .false.
      startword = .true.
      nvals = 0
      do i=1,len_trim(optarg)
         space = ( optarg(i:i) == " " .or. optarg(i:i) == "," )
         ! Only count spaces that end a word
         inword = .not. space
         if ( inword .and. startword ) then
            nvals = nvals + 1
         end if
         startword = space
      end do
      nmax = size(arr)
      if ( nvals + n > nmax ) then
         ! Would overflow array to read these. Return as an error flag
         n = n+nvals
         arr = -huge(1.)
         return
      end if

      read(optarg,*) arr(n+1:n+nvals)
      n = n + nvals

   end subroutine str2rarray

   subroutine str2iarray(optarg,arr,n)
      ! Convert a string of values delimited by spaces or commas into 
      ! separate integer values.
      character(len=*), intent(in) :: optarg
      integer, dimension(:), intent(inout) :: arr
      integer, intent(inout) :: n
      integer :: nvals, i, nmax
      logical :: inword, space, startword

!     Go through the string counting the number of distinct words (delimited
!     by commas or spaces
      inword = .false.
      startword = .true.
      nvals = 0
      do i=1,len_trim(optarg)
         space = ( optarg(i:i) == " " .or. optarg(i:i) == "," )
         ! Only count spaces that end a word
         inword = .not. space
         if ( inword .and. startword ) then
            nvals = nvals + 1
         end if
         startword = space
      end do
      nmax = size(arr)
      if ( nvals + n > nmax ) then
         ! Would overflow array to read these. Return as an error flag
         n = n+nvals
         arr = -huge(1)
         return
      end if

      read(optarg,*) arr(n+1:n+nvals)
      n = n + nvals

   end subroutine str2iarray

   subroutine str2carray(optarg,arr,n)
      ! Convert a string of values delimited by spaces or commas into 
      ! separate words
      character(len=*), intent(in) :: optarg
      character(len=*), dimension(:), intent(inout) :: arr
      integer, intent(inout) :: n
      integer :: i, nmax, iw
      logical :: inword, space, startword

      nmax = size(arr)
!     Go through the string counting the number of distinct words (delimited
!     by commas or spaces
      inword = .false.
      startword = .true.
      do i=1,len_trim(optarg)
         space = ( optarg(i:i) == " " .or. optarg(i:i) == "," )
         ! Only count spaces that end a word
         inword = .not. space
         if ( inword .and. startword ) then
            n = n + 1
            arr(n) = ""
            iw = 0
         end if
         if ( n > nmax ) then
            ! Would overflow array to read these. Return as an error flag
            return
         end if
         if ( inword ) then
            iw = iw + 1
            arr(n)(iw:iw) = optarg(i:i)
         end if
         startword = space
      end do

   end subroutine str2carray

   function isrchfle(array, target) result (ind)
      real, dimension(:), intent(in) :: array
      real, intent(in) :: target
      integer :: ind
      do ind=1,size(array)
         if ( array(ind) <= target ) then
            exit
         end if
      end do
      ! ind = size(array+1) if it gets to the end of the loop without exiting
   end function isrchfle

   function isrchfgt(array, target) result (ind)
      real, dimension(:), intent(in) :: array
      real, intent(in) :: target
      integer :: ind
      do ind=1,size(array)
         if ( array(ind) > target ) then
            exit
         end if
      end do
      ! ind = size(array+1) if it gets to the end of the loop without exiting
   end function isrchfgt

   function vsrchfgt(array, target, nl, nx) result (ind)
      ! Vector version is more efficient with array sizes passed as arguments
      ! This version requires that the array is in increasing order.
      ! If it's decreasing use 
      ! nl+2 - search_fgt(sig(nl:1:-1),target,nl,nx)
      ! Note that this is not checked!
      integer, intent(in) :: nl, nx
      real, dimension(nl), intent(in) :: array
      real, dimension(nx), intent(in) :: target
      integer, dimension(nx) :: ind
      real, dimension(nx) :: rind
      integer :: j

      rind = 1.
      do j=1,nl
         ! Set count to be 1. if it's below the pressure level, 0 if above.
         rind(:) = rind(:) + max( sign(1.0,(target-array(j))), 0.)
      end do
      ind = int ( rind+0.1 )

   end function vsrchfgt

   function vsrchfgt_dp(array, target, nl, nx) result (ind)
      ! Vector version is more efficient with array sizes passed as arguments
      ! This version requires that the array is in increasing order.
      ! If it's decreasing use 
      ! nl+2 - search_fgt(sig(nl:1:-1),target,nl,nx)
      ! Note that this is not checked!
      integer, intent(in) :: nl, nx
      real(kind=dp), dimension(nl), intent(in) :: array
      real(kind=dp), dimension(nx), intent(in) :: target
      integer, dimension(nx) :: ind
      real(kind=dp), dimension(nx) :: rind
      integer :: j

      rind = 1.
      do j=1,nl
         ! Set count to be 1. if it's below the pressure level, 0 if above.
         rind(:) = rind(:) + max( sign(1.0_dp,(target-array(j))), 0._dp)
      end do
      ind = int ( rind+0.1 )

   end function vsrchfgt_dp

   function vxsrchfgt(array, target, nl, nx) result (ind)
      ! Vector version is more efficient with array sizes passed as arguments
      ! This version requires that the array is in increasing order.
      ! If it's decreasing use 
      ! nl+2 - search_fgt(sig(nl:1:-1),target,nl,nx)
      ! Note that this is not checked!

      ! In this version array is two dimensional
      integer, intent(in) :: nl, nx
      real, dimension(nx,nl), intent(in) :: array
      real, dimension(nx), intent(in) :: target
      integer, dimension(nx) :: ind
      real, dimension(nx) :: rind
      integer :: i, j

      rind = 1.
      do j=1,nl
         ! Set count to be 1. if it's below the pressure level, 0 if above.
         do i=1,nx
            rind(i) = rind(i) + max( sign(1.0,(target(i)-array(i,j))), 0.)
         end do
      end do
      ind = int ( rind+0.1 )

   end function vxsrchfgt

   function vxsrchfgt_dp(array, target, nl, nx) result (ind)
      ! Vector version is more efficient with array sizes passed as arguments
      ! This version requires that the array is in increasing order.
      ! If it's decreasing use 
      ! nl+2 - search_fgt(sig(nl:1:-1),target,nl,nx)
      ! Note that this is not checked!

      ! In this version array is two dimensional
      integer, intent(in) :: nl, nx
      real(kind=dp), dimension(nx,nl), intent(in) :: array
      real(kind=dp), dimension(nx), intent(in) :: target
      integer, dimension(nx) :: ind
      real(kind=dp), dimension(nx) :: rind
      integer :: i, j

      rind = 1.
      do j=1,nl
         ! Set count to be 1. if it's below the pressure level, 0 if above.
         do i=1,nx
            rind(i) = rind(i) + max( sign(1.0_dp,(target(i)-array(i,j))), 0._dp)
         end do
      end do
      ind = int ( rind+0.1 )

   end function vxsrchfgt_dp

   elemental function fpequal_sp ( x, y, tol_in ) result (res)
      ! Floating point equality test
      real, intent(in) :: x, y
      real, intent(in), optional :: tol_in
      logical :: res
      real :: tol
      if ( present(tol_in) ) then
         tol = tol_in
      else
         tol = 100.
      end if
      res = abs(x-y) < tol*2*spacing(max(abs(x),abs(y)))
   end function fpequal_sp

   elemental function fpequal_dp ( x, y, tol_in ) result (res)
      ! Floating point equality test
      real(kind=dp), intent(in) :: x, y
      real(kind=dp), intent(in), optional :: tol_in
      logical :: res
      real(kind=dp) :: tol
      if ( present(tol_in) ) then
         tol = tol_in
      else
         tol = 100.
      end if
      res = abs(x-y) < tol*2*spacing(max(abs(x),abs(y)))
   end function fpequal_dp

   function lcase ( str ) result(lstr)
      ! Convert string to lower case
      character(len=*), intent(in) :: str
      character(len=len(str)) :: lstr
      integer :: ia, is
      integer, parameter :: aoffset = iachar('a') - iachar('A')


      do is=1,len(str)
         ia = iachar(str(is:is))
         if (ia >= iachar('A').and.ia <= iachar('Z')) then
            lstr(is:is) = achar(ia+aoffset)
         else
            lstr(is:is) = achar(ia)
         endif
      end do
   end function lcase
   
end module utils_m
