! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module gldata
   implicit none
!  Arrays for fields that need some extra processing and so 
!  can't be handled by readsave
   real, dimension(:,:),   allocatable, save :: psl, zs, soilt, urban_frac
   real, dimension(:,:),   allocatable, save :: f_cor
   real, dimension(:,:,:), allocatable, save :: u, v, t, q, rh, ql, qf, qs, qg
   real, dimension(:,:,:), allocatable, save :: omega
   real, dimension(:),     allocatable, save :: sig, dsig, zsoil, zse, gosig
   real, dimension(:),     allocatable, save :: shallow_zse
   integer, dimension(:),   allocatable, save :: cable_patch, cable_cohort
   real, dimension(:,:,:), allocatable, save :: tgg
   real, allocatable, dimension(:,:,:) :: zstd, hstd, tstd
   real, allocatable, dimension(:,:,:) :: tmp3d
   real, allocatable, dimension(:,:,:) :: uo_tmp, vo_tmp, thetao_tmp, so_tmp
   real, allocatable, dimension(:,:,:) :: ocn_tmp
   real, allocatable, dimension(:,:,:) :: cp_tmp
   real, allocatable, dimension(:,:,:,:) :: cpc_tmp
   logical, allocatable, dimension(:,:,:) :: ocn_mask
end module gldata
