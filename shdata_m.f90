! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
#ifdef usempi3
module shdata_m
   implicit none

   private
   public :: allocshdata, freeshdata

   interface allocshdata
      module procedure allocshdata_r2,allocshdata_r3,allocshdata_i2,allocshdata_i3
   end interface

contains

   subroutine allocshdata_r2(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
      include 'mpif.h'
      real, pointer, dimension(:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_REAL,tsize,ierr)
      lsize = ssize*tsize
      call MPI_Win_allocate_shared( lsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_r2    
    
   subroutine allocshdata_r3(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
      include 'mpif.h'
      real, pointer, dimension(:,:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_REAL,tsize,ierr)
      lsize = ssize*tsize
      call MPI_Win_allocate_shared( lsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_r3

   subroutine allocshdata_i2(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
      include 'mpif.h'
      integer, pointer, dimension(:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_INTEGER,tsize,ierr)
      lsize = ssize*tsize
      call MPI_Win_allocate_shared( lsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_i2

   subroutine allocshdata_i3(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
      include 'mpif.h'
      integer, pointer, dimension(:,:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize, lsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_INTEGER,tsize,ierr)
      lsize = ssize*tsize
      call MPI_Win_allocate_shared( lsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_i3
   
   subroutine freeshdata(win)
      include 'mpif.h'
      integer, intent(in) :: win
      integer :: ierr

      call MPI_Win_Free(win, ierr) 

   end subroutine freeshdata

end module shdata_m
#endif
