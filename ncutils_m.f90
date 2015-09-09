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
    
module ncutils_m
   ! Generally useful netcdf stuff
#ifdef usenc3
   use netcdf_m
#else
   use netcdf
#endif
   implicit none
   private
   public ::  check_ncerr, copy_atts, history_append, fix_ncatt

contains

   subroutine check_ncerr(status,mesg)
      integer, intent(in) :: status
      character(len=*), intent(in), optional :: mesg
      if ( status /= nf90_noerr ) then
         if ( present(mesg) ) then
            print*, mesg
         end if
         print*, trim(nf90_strerror(status))
         stop
      end if
   end subroutine check_ncerr

   subroutine copy_atts(ncid_in, ncid_out, vid_in, vid_out)
      ! Copy attributes from variable with vid_in in file ncid_in to
      ! variable vid_out in file ncid_out. 
      ! Note that vid may be NF90_GLOBAL
      integer, intent(in) :: ncid_in, ncid_out, vid_in, vid_out
      integer :: iatt
      character(len=nf90_max_name) :: attname
      integer :: ierr, natts

      if ( vid_in == NF90_GLOBAL ) then
         ierr = nf90_inquire ( ncid_in, nattributes=natts )
      else
         ierr = nf90_inquire_variable ( ncid_in, vid_in, natts=natts )
      end if
      call check_ncerr ( ierr, "copy_atts: Error getting natts")
      do iatt=1,natts
         ierr = nf90_inq_attname(ncid_in, vid_in, iatt, attname)
         call check_ncerr(ierr, "Error getting attribute name in copy_atts")
         ierr = nf90_copy_att(ncid_in, vid_in, attname, ncid_out, vid_out)
         call check_ncerr(ierr, "Error copying attribute")
      end do
   end subroutine copy_atts

   subroutine history_append(ncid, string)
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: string
      integer, parameter :: maxlen=1000
      character(maxlen) :: history
      integer :: attlen, ierr, ierr2, atttype

      ierr = nf90_inquire_attribute(ncid, NF90_GLOBAL, "history", xtype=atttype, len=attlen )
      if ( ierr == NF90_NOERR ) then
         history = ""
         if ( atttype /= NF90_CHAR ) then
            print*, "Error from history_append, history attribute not of type character"
            stop
         end if
         if ( attlen + len_trim(string) + 1 > maxlen ) then
            print*, "Error, increase maxlen in history_append"
            stop
         end if
         ierr2 = nf90_get_att ( ncid, NF90_GLOBAL, "history", history )
         call check_ncerr ( ierr2, "Error reading history attribute ")  
      else if ( ierr == NF90_ENOTATT ) then
         ! Not a fatal error for the attribute to be missing
         history = ""
      else
         call check_ncerr(ierr, "Error getting history attribute")
      end if

      ! Add a newline before the new string
      history = history(1:len_trim(history)) // char(10) // trim(string)
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "history", trim(history))
      call check_ncerr(ierr,"Error defining character attribute")

   end subroutine history_append

   subroutine fix_ncatt ( attval )
      ! Remove any trailing nulls
      character(len=*), intent(inout) :: attval
      integer :: i
      do i=1,len(attval)
         if ( ichar(attval(i:i)) == 0 ) then
            attval(i:i) = " "
         end if
      end do
   end subroutine fix_ncatt
   
end module ncutils_m
