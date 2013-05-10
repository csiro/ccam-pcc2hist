module checkver_m

! This module provides a routine to check version numbers of library modules.
! The version number must be a public character parameter in the module
! for this to work.
! The assumption is that the major version number will be increased
! whenever there is a change in the interface.

   implicit none

contains

   subroutine checkver ( name, revision, major_req, minor_req )
      character(len=*), intent(in) :: name
      character(len=*), intent(in) :: revision
      integer, intent(in) :: major_req
      integer, intent(in), optional :: minor_req
      integer :: start, end, major, minor

!     Routine gives an error unless library module major version number 
!     is equal to the requested value. If the minor argument is given, the
!     library minor version number must be greater than or equal to the 
!     requested value.

!     Parse the module revision string      
!     Major number comes after :
      start=index(revision,"Revision:")
      if ( start == 0 ) then
         print*, " Error in checkver. Failed to parse argument ", &
                 trim(revision)
      end if
      end = index(revision,".")
      if ( end == 0 ) then
         print*, " Error in checkver. Failed to parse argument ", &
                 trim(revision)
      end if
      read ( revision(start+9:end-1),* ) major
      start = end+1
      read ( revision(start:len(revision)),* ) minor

      if ( major_req /= major ) then
         print*, " CHECKVER ERROR: ", trim(name), " ", trim(revision), &
                 " Requested major version", major_req
         stop
      end if
      if ( present(minor_req) ) then
         if ( minor_req > minor ) then
            print*, " CHECKVER ERROR: ", trim(name), " ", trim(revision), &
                    " Requested minor version", minor_req
            stop
         end if
      end if

   end subroutine checkver
end module checkver_m

