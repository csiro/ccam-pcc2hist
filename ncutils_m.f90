module ncutils_m
   ! Generally useful netcdf stuff
#ifndef PARNETCDF
   use netcdf_m
#else
   use pnetcdf_m
#endif
   implicit none
   private
   public ::  check_ncerr, copy_atts, history_append, fix_ncatt

contains

   subroutine check_ncerr(status,mesg)
      integer, intent(in) :: status
      character(len=*), intent(in), optional :: mesg
      if ( status /= ncf90_noerr ) then
         if ( present(mesg) ) then
            print*, mesg
         end if
         print*, trim(ncf90_strerror(status))
         stop
      end if
   end subroutine check_ncerr

   subroutine copy_atts(ncid_in, ncid_out, vid_in, vid_out)
      ! Copy attributes from variable with vid_in in file ncid_in to
      ! variable vid_out in file ncid_out. 
      ! Note that vid may be NCF90_GLOBAL
      integer, intent(in) :: ncid_in, ncid_out, vid_in, vid_out
      integer :: iatt
      character(len=ncf90_max_name) :: attname
      integer :: ierr, natts

      if ( vid_in == NCF90_GLOBAL ) then
         ierr = ncf90_inquire ( ncid_in, nattributes=natts )
      else
         ierr = ncf90_inquire_variable ( ncid_in, vid_in, natts=natts )
      end if
      call check_ncerr ( ierr, "copy_atts: Error getting natts")
      do iatt=1,natts
         ierr = ncf90_inq_attname(ncid_in, vid_in, iatt, attname)
         call check_ncerr(ierr, "Error getting attribute name in copy_atts")
         ierr = ncf90_copy_att(ncid_in, vid_in, attname, ncid_out, vid_out)
         call check_ncerr(ierr, "Error copying attribute")
      end do
   end subroutine copy_atts

   subroutine history_append(ncid, string)
#ifdef PARNETCDF
      use mpi, only : MPI_OFFSET_KIND
#endif
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: string
      integer, parameter :: maxlen=1000
      character(maxlen) :: history
#ifdef PARNETCDF
      integer(kind=MPI_OFFSET_KIND) :: attlen
      integer ierr, ierr2, atttype
#else
      integer :: attlen, ierr, ierr2, atttype
#endif
      ierr =  ncf90_inquire_attribute(ncid, NCF90_GLOBAL, "history", xtype=atttype, len=attlen )

      if ( ierr == NCF90_NOERR ) then
         history = ""
         if ( atttype /= NCF90_CHAR ) then
            print*, "Error from history_append, history attribute not of type character"
            stop
         end if
         if ( attlen + len_trim(string) + 1 > maxlen ) then
            print*, "Error, increase maxlen in history_append"
            stop
         end if
         ierr2 = ncf90_get_att ( ncid, NCF90_GLOBAL, "history", history )
         call check_ncerr ( ierr2, "Error reading history attribute ")  
      else if ( ierr == NCF90_ENOTATT ) then
         ! Not a fatal error for the attribute to be missing
         history = ""
      else
         call check_ncerr(ierr, "Error getting history attribute")
      end if

      ! Add a newline before the new string
      history = history(1:len_trim(history)) // char(10) // trim(string)
      ierr = ncf90_put_att ( ncid, NCF90_GLOBAL, "history", trim(history))
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
