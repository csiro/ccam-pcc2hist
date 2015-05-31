module mpidata_m

implicit none

private

integer, dimension(:), save, allocatable, public :: ncid_in
integer, dimension(:,:), save, allocatable, public :: ioff, joff
integer, save, public :: myid, nproc
integer, save, public :: pil, pjl, pnpan, pnproc, lproc
integer, save, public :: pil_g, pjl_g

end module mpidata_m