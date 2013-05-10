module s2p_m

!  Module for converting from sigma to pressure levels

   use sitop_m
   implicit none
!  Flag to control whether output is on sigma or pressure levels
   logical, save :: use_plevs = .false. 
   integer, parameter :: npmax = 50
   real, dimension(npmax) :: plevs = 0.0
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

   subroutine vsavehist ( name, array )

!     Version of savehist that optionally does the sigma to pressure 
!     conversion.

      use history, only : savehist, needfld
      use gldata, only : sig, psl
      character(len=*), intent(in) :: name
      real, dimension(:,:,:), intent(in) :: array
!     Make this allocatable rather than automatic because usually it won't be
!     required (check whether this really does save memory. If nplevs=0 for
!     sigma case then it would be zero sized anyway?).
      real, dimension(:,:,:), allocatable :: parray

      if ( .not. needfld(name) ) return
      if ( use_plevs ) then
!        sigma to pressure conversion.
         allocate ( parray(size(array,1), size(array,2), nplevs) )
         if ( vextrap == vextrap_default ) then
            if (  name == "temp" ) then
               call sitop ( array, parray, sig, plevs(1:nplevs), psl, vextrap_t )
            else
               call sitop ( array, parray, sig, plevs(1:nplevs), psl, vextrap_none )
            end if
         else
            call sitop ( array, parray, sig, plevs(1:nplevs), psl, vextrap )
         end if
         call savehist ( name, parray )
         deallocate ( parray )
      else
!        Save the sigma level values directly.
         call savehist ( name, array(:,:,minlev:maxlev) )
      end if
   end subroutine vsavehist
   
end module s2p_m
