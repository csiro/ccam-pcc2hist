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
    
module mpidata_m

implicit none

private

integer, dimension(:), save, allocatable, public :: ncid_in
#ifdef procformat
integer, dimension(:), save, allocatable, public :: ip_min, ip_max
integer, save, public :: ip_maxcnt
#endif
#ifdef parallel_int
integer, dimension(:,:), save, pointer, contiguous, public :: ioff, joff
integer, dimension(:,:,:), save, pointer, contiguous, public :: ijoff
#else
integer, dimension(:,:), save, allocatable, public :: ioff, joff
#endif
integer, save, public :: myid, nproc
integer, save, public :: pil, pjl, pnpan, pnproc, lproc
integer, save, public :: pil_g, pjl_g
#ifdef parallel_int
integer, save, public :: node_comm, node_myid, node_nproc
integer, save, public :: ijoff_win
#endif
#ifdef usefirstrank
integer, save, public :: node2_comm, node2_myid, node2_nproc
#endif
#ifdef procformat
integer, save, public :: proc_node
#endif

end module mpidata_m
