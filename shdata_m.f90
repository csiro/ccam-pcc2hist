module shdata_m
   implicit none

   private
   public :: cpshdata,allocshdata
   private :: cpshdata_r3,cpshdata_i2,cpshdata_i3
   private :: allocshdata_r3,allocshdata_i2,allocshdata_i3

   interface cpshdata
      module procedure cpshdata_r3,cpshdata_i2,cpshdata_i3
   end interface
   interface allocshdata
      module procedure allocshdata_r3,allocshdata_i2,allocshdata_i3
   end interface

contains

   subroutine cpshdata_r3(pdata, ssize)
      use mpidata_m
#ifdef usenc3
      include 'mpif.h'
#else
      use mpi
#endif
      real, dimension(:,:,:), intent(in) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer(kind=MPI_ADDRESS_KIND), parameter :: disp=0
      integer :: win,ierr,tsize

!     copy the shared data from rank 0 to the first rank on each node
      call MPI_Type_size(MPI_REAL,tsize,ierr)
      call MPI_Win_create(pdata,ssize*tsize,1,MPI_INFO_NULL,MPI_COMM_WORLD,win,ierr)
      call MPI_Win_fence(0,win,ierr)
      if ( node_myid == 0 .and. myid /= 0 ) then
          call MPI_Get(pdata,ssize,MPI_REAL,0,disp,ssize,MPI_REAL,win,ierr)
      end if
      call MPI_Win_fence(0,win,ierr)
      call MPI_Win_free(win,ierr)

   end subroutine cpshdata_r3

   subroutine cpshdata_i2(pdata, ssize)
      use mpidata_m
#ifdef usenc3
      include 'mpif.h'
#else
      use mpi
#endif
      integer, dimension(:,:), intent(in) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer(kind=MPI_ADDRESS_KIND), parameter :: disp=0
      integer :: win,ierr,tsize

!     copy the shared data from rank 0 to the first rank on each node
      call MPI_Type_size(MPI_INTEGER,tsize,ierr)
      call MPI_Win_create(pdata,ssize*tsize,1,MPI_INFO_NULL,MPI_COMM_WORLD,win,ierr)
      call MPI_Win_fence(0,win,ierr)
      if ( node_myid == 0 .and. myid /= 0 ) then
          call MPI_Get(pdata,ssize,MPI_INTEGER,0,disp,ssize,MPI_INTEGER,win,ierr)
      end if
      call MPI_Win_fence(0,win,ierr)
      call MPI_Win_free(win,ierr)

   end subroutine cpshdata_i2

   subroutine cpshdata_i3(pdata, ssize)
      use mpidata_m
#ifdef usenc3
      include 'mpif.h'
#else
      use mpi
#endif
      integer, dimension(:,:,:), intent(in) :: pdata 
      integer(kind=MPI_ADDRESS_KIND), intent(in) :: ssize
      integer(kind=MPI_ADDRESS_KIND), parameter :: disp=0
      integer :: win,ierr,tsize

!     copy the shared data from rank 0 to the first rank on each node
      call MPI_Type_size(MPI_INTEGER,tsize,ierr)
      call MPI_Win_create(pdata,ssize*tsize,1,MPI_INFO_NULL,MPI_COMM_WORLD,win,ierr)
      call MPI_Win_fence(0,win,ierr)
      if ( node_myid == 0 .and. myid /= 0 ) then
          call MPI_Get(pdata,ssize,MPI_INTEGER,0,disp,ssize,MPI_INTEGER,win,ierr)
      end if
      call MPI_Win_fence(0,win,ierr)
      call MPI_Win_free(win,ierr)

   end subroutine cpshdata_i3

   subroutine allocshdata_r3(pdata, ssize, sshape, win)
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
#ifdef usenc3
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
#ifdef usenc3
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
#ifdef usenc3
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
