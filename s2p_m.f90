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
    
module s2p_m

!  Module for converting from sigma to pressure levels

   use sitop_m
   implicit none
!  Flag to control whether output is on sigma or pressure levels
   logical, save :: use_plevs = .false. 
   logical, save :: use_meters = .false.
   integer, parameter :: npmax = 200
   real, dimension(npmax) :: plevs = 0.0
   real, dimension(npmax) :: mlevs = 9.E9
   integer, save :: nplevs
   integer :: vextrap = vextrap_default
   integer, public, save :: minlev=1, maxlev=HUGE(1)

contains
   subroutine check_plevs
      ! Sanity checks and sorting on the specified output pressure levels.
      use usage_m
      integer :: i, j
      real :: temp

      if ( .not. use_plevs ) then
         return
      end if
      if ( maxval(plevs) == 0.0 ) then
         print*, " Error, use_plevs set but no levels specified "
         call usage()
      end if
!     Find the actual number of pressure levels set.
!     First remove any negative values
      plevs = max ( plevs, 0.0 )
!     Sort into decreasing order. This is only done once on a short list
!     so just use a selection sort. 
      do i=1,npmax
!        Using maxval here is just a trick to reduce the array result of 
!        maxloc to a scalar.
         j = maxval ( maxloc(plevs(i:)) ) + i - 1
         temp = plevs(i)
         plevs(i) = plevs(j)
         plevs(j) = temp
      end do
!     Find the maximum non-zero value
      do i=1,npmax
         if ( plevs(i) == 0 ) then
            exit
         end if
      end do
!     Gives correct value whether the loop exits naturally or not
      nplevs = i-1  
      if  ( nplevs == 0 ) then
         print*,  "Error, no pressure levels set "
      end if
   end subroutine check_plevs

   subroutine check_meters
      ! Sanity checks and sorting on the specified output pressure levels.
      use usage_m
      integer :: i, j
      real :: temp

      if ( .not. use_meters ) then
         return
      end if
      if ( maxval(mlevs) == 0.0 ) then
         print*, " Error, use_mlevs set but no levels specified "
         call usage()
      end if
!     Find the actual number of pressure levels set.
!     First remove any negative values
      mlevs = max ( mlevs, 0.0 )
!     Sort into increasing order. This is only done once on a short list
!     so just use a selection sort. 
      do i=1,npmax
!        Using maxval here is just a trick to reduce the array result of 
!        maxloc to a scalar.
         j = maxval ( minloc(mlevs(i:)) ) + i - 1
         temp = mlevs(i)
         mlevs(i) = mlevs(j)
         mlevs(j) = temp
      end do
!     Find the minimum non-huge value
      do i=1,npmax
         if ( mlevs(i) == 9.E9 ) then
            exit
         end if
      end do
!     Gives correct value whether the loop exits naturally or not
      nplevs = i-1  
      if  ( nplevs == 0 ) then
         print*,  "Error, no height levels set "
         stop
      end if
   end subroutine check_meters   
   
   subroutine vsavehist ( name, array )

!     Version of savehist that optionally does the sigma to pressure 
!     conversion.

      use history, only : savehist, needfld
      use gldata, only : sig, psl, hstd
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: array
!     Make this allocatable rather than automatic because usually it won't be
!     required (check whether this really does save memory. If nplevs=0 for
!     sigma case then it would be zero sized anyway?).
      real, dimension(:,:,:), allocatable, save :: parray

      if ( .not. needfld(name) ) return
      call START_LOG(vsavehist_begin)
      if ( use_plevs ) then
!        sigma to pressure conversion.
         if ( .not. allocated(parray) ) then
           allocate ( parray(size(array,1), size(array,2), nplevs) )
         else if ( size(parray,1)/=size(array,1) .or. size(parray,2)/=size(array,2) .or. size(parray,3)/=nplevs ) then
           deallocate ( parray )
           allocate ( parray(size(array,1), size(array,2), nplevs) )
         end if
         if ( vextrap == vextrap_default ) then
            if ( name == "temp" ) then
               call sitop ( array, parray, sig, plevs(1:nplevs), psl, vextrap_t )
            else
               call sitop ( array, parray, sig, plevs(1:nplevs), psl, vextrap_none )
            end if
         else
            call sitop ( array, parray, sig, plevs(1:nplevs), psl, vextrap )
         end if
         call savehist ( name, parray )
      else if ( use_meters ) then
!        sigma to meters conversion.
         if ( .not. allocated(parray) ) then
           allocate ( parray(size(array,1), size(array,2), nplevs) )
         else if ( size(parray,1)/=size(array,1) .or. size(parray,2)/=size(array,2) .or. size(parray,3)/=nplevs ) then
           deallocate ( parray )
           allocate ( parray(size(array,1), size(array,2), nplevs) )
         end if
         if ( vextrap == vextrap_default ) then
            if (  name == "temp" ) then
               call mitop ( array, parray, vextrap_t )
            else
               call mitop ( array, parray, vextrap_none )
            end if
         else
            call mitop ( array, parray, vextrap )
         end if
         call savehist ( name, parray )
      else
!        Save the sigma level values directly.
         call savehist ( name, array(:,:,minlev:maxlev) )
      end if
      call END_LOG(vsavehist_end)
   end subroutine vsavehist
   
end module s2p_m
