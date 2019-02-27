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
    
module newmpar_m
   implicit none
   integer, public, save :: il = 20
   integer, public, save :: npanels = 5
   integer, public, save :: jl
   integer, public, save :: kl = 1
   integer, public, save :: ol = 0
   integer, public, save :: cptch = 0
   integer, public, save :: cchrt = 0
   integer, public, save :: ifull
   integer, public, save :: ij
   integer, public, save :: ijk
   integer, public, save :: iquad
end module newmpar_m
   
