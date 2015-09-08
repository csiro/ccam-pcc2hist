#ifdef parallel_int
module shdata_m
   implicit none

   private
   public :: allocshdata
   private :: allocshdata_r3,allocshdata_i2,allocshdata_i3

   interface allocshdata
      module procedure allocshdata_r3,allocshdata_i2,allocshdata_i3
   end interface

contains

   subroutine allocshdata_r3(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
#ifdef usempif
      include 'mpif.h'
#else
      use mpi
#endif
      real, pointer, dimension(:,:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_REAL,tsize,ierr)
      call MPI_Win_allocate_shared( ssize*tsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_r3

   subroutine allocshdata_i2(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
#ifdef usempif
      include 'mpif.h'
#else
      use mpi
#endif
      integer, pointer, dimension(:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_INTEGER,tsize,ierr)
      call MPI_Win_allocate_shared( ssize*tsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_i2

   subroutine allocshdata_i3(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
#ifdef usempif
      include 'mpif.h'
#else
      use mpi
#endif
      integer, pointer, dimension(:,:,:), intent(inout) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer, dimension(:), intent(in) :: sshape
      integer, intent(out) :: win
      integer(kind=MPI_ADDRESS_KIND) :: qsize
      integer :: disp_unit,ierr,tsize
      type(c_ptr) :: baseptr

!     allocted a single shared memory region on each node
      call MPI_Type_size(MPI_INTEGER,tsize,ierr)
      call MPI_Win_allocate_shared( ssize*tsize, 1, MPI_INFO_NULL, node_comm, baseptr, win, ierr)
      if ( node_myid /= 0 ) call MPI_Win_shared_query(win, 0, qsize, disp_unit, baseptr, ierr)
      call c_f_pointer(baseptr, pdata, sshape)

   end subroutine allocshdata_i3

end module shdata_m
#endif
