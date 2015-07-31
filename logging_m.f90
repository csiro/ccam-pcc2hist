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
    
module logging_m

#ifndef usenc3
   use mpi
#endif
   use mpidata_m
implicit none
#ifdef usenc3
   include 'mpif.h'
#endif

private

   public :: start_log, end_log, log_setup
! Timer
   integer, public, save :: model_begin, model_end
   integer, public, save :: paraopen_begin, paraopen_end
   integer, public, save :: paraclose_begin, paraclose_end
   integer, public, save :: openhist_begin, openhist_end
   integer, public, save :: closehist_begin, closehist_end
   integer, public, save :: writehist_begin, writehist_end
   integer, public, save :: ints_begin, ints_end
   integer, public, save :: timeloop_begin, timeloop_end
   integer, public, save :: infile_begin, infile_end
   integer, public, save :: savehist_begin, savehist_end
   integer, public, save :: vsavehist_begin, vsavehist_end
   integer, public, save :: vread_begin, vread_end
   integer, public, save :: paravar2a_begin, paravar2a_end
   integer, public, save :: paravar3a_begin, paravar3a_end
   integer, public, save :: mpiscatter_begin, mpiscatter_end
#ifdef simple_timer
   public :: simple_timer_finalize
#endif
   integer, parameter :: nevents = 15
   real(kind=8), dimension(nevents), save :: tot_time = 0., start_time
   character(len=15), dimension(nevents), save :: event_name

contains

   subroutine start_log ( event )
      integer, intent(in) :: event
      integer :: ierr
#ifdef vampir
      VT_USER_START(event_name(event))
#endif
#ifdef simple_timer
      start_time(event) = MPI_Wtime()
#endif
   end subroutine start_log

   subroutine end_log ( event )
      integer, intent(in) :: event
      integer :: ierr
#ifdef vampir
      VT_USER_END(event_name(event))
#endif

#ifdef simple_timer
      tot_time(event) = tot_time(event) + MPI_Wtime() - start_time(event)
#endif
   end subroutine end_log

   subroutine log_setup()
      integer :: ierr
      integer :: classhandle
#ifdef vampir
#ifdef simple_timer
      write(6,*) "ERROR: vampir and simple_timer should not be compiled together"
      stop
#endif
#endif

      model_begin = 1
      model_end =  model_begin
      event_name(model_begin) = "Main"

      paraopen_begin = 2
      paraopen_end =  paraopen_begin
      event_name(paraopen_begin) = "Paraopen"

      paraclose_begin = 3
      paraclose_end =  paraclose_begin
      event_name(paraclose_begin) = "Paraclose"

      openhist_begin = 4
      openhist_end =  openhist_begin
      event_name(openhist_begin) = "Openhist"

      closehist_begin = 5
      closehist_end =  closehist_begin
      event_name(closehist_begin) = "Closehist"

      writehist_begin = 6
      writehist_end =  writehist_begin
      event_name(writehist_begin) = "Writehist"

      ints_begin = 7
      ints_end =  ints_begin
      event_name(ints_begin) = "Ints"

      timeloop_begin = 8
      timeloop_end =  timeloop_begin
      event_name(timeloop_begin) = "Timeloop"

      infile_begin = 9
      infile_end =  infile_begin
      event_name(infile_begin) = "Infile"

      savehist_begin = 10
      savehist_end =  savehist_begin
      event_name(savehist_begin) = "Savehist"

      vsavehist_begin = 11
      vsavehist_end =  vsavehist_begin
      event_name(vsavehist_begin) = "Vsavehist"

      vread_begin = 12
      vread_end =  vread_begin
      event_name(vread_begin) = "Vread"

      paravar2a_begin = 13
      paravar2a_end =  paravar2a_begin
      event_name(paravar2a_begin) = "Paravar2a"

      paravar3a_begin = 14
      paravar3a_end =  paravar3a_begin
      event_name(paravar3a_begin) = "Paravar3a"

      mpiscatter_begin = 15
      mpiscatter_end =  mpiscatter_begin
      event_name(mpiscatter_begin) = "MPIScatter"

   end subroutine log_setup

#ifdef simple_timer
   subroutine simple_timer_finalize()
      ! Calculate the mean, min and max times for each case
      integer :: i
      integer(kind=4) :: ierr, llen
      real(kind=8), dimension(nevents) :: emean, emax, emin
      llen=nevents
      call MPI_Reduce(tot_time, emean, llen, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0_4, MPI_COMM_WORLD, ierr )
      call MPI_Reduce(tot_time, emax, llen, MPI_DOUBLE_PRECISION, &
                      MPI_MAX, 0_4, MPI_COMM_WORLD, ierr )
      call MPI_Reduce(tot_time, emin, llen, MPI_DOUBLE_PRECISION, &
                      MPI_MIN, 0_4, MPI_COMM_WORLD, ierr )
      if ( myid == 0 ) then
         write(6,*) "==============================================="
         write(6,*) "  Times over all processes"
         write(6,*) "  Routine        Mean time  Min time  Max time"
         do i=1,nevents
            if ( emean(i) > 0. ) then
               ! This stops boundsa, b getting written when they're not used.
               write(*,"(a,3f10.3)") event_name(i), emean(i)/nproc, emin(i), emax(i)
            end if
         end do
      end if

   end subroutine simple_timer_finalize
#endif

end module logging_m
