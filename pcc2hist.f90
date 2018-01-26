! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
program cc2hist

!  This program generates a netcdf history file from the conformal-cubic
!  or conformal-octagon model output file. The interpolation routines are
!  taken from John's plotting program plotg.f.

!  It uses the f90 shallow water model versions of routines setxyz and staguv.
    
!  Modified by MJT to use MPI when reading parallel input files (pcc2hist).
!  Further optimised by Paul Ryan and David Benn.

#ifdef usempi_mod
   use mpi
#endif
    
   use history
   use getopt_m
   use mpidata_m
   use newmpar_m
   use work
   use usage_m
   use height_m
   use s2p_m
   use interp_m
   use checkver_m
   use parm_m, only : rlong0, rlat0, schmidt
   use logging_m

   implicit none

#ifndef usempi_mod
   include 'mpif.h'
#endif

   character(len=MAX_ARGLEN) :: ifile, ofile, cfile

   real :: hres = 0.
   real :: minlon = 0.0, maxlon=360.0, dlon=0.0, &
           minlat=-90.0, maxlat=90.0,  dlat=0.0
   real :: dx=0., dy=0.  ! for TAPM
   integer :: lx=0, ly=0 ! for TAPM

   character(len=MAX_ARGLEN) :: optarg
   integer :: opt, nopt
   type(loption), dimension(4) :: longopts
   integer :: longind
   integer :: kta=0, ktb=999999, ktc=-1, ndate=-1, ntime=-1, k
   integer :: sdate=-1, edate=-1, stime=-1, etime=-1
   character(len=10) :: name
   character(len=80) :: longname
   real :: minsig = 0., maxsig = 1.0

   namelist /input/ kta, ktb, ktc, ndate, ntime,                      &
                    minlon, maxlon, dlon, minlat, maxlat, dlat,       &
                    dx, dy, lx, ly,                                   &
                    minlev, maxlev, minsig, maxsig, use_plevs, plevs, &
                    use_meters, mlevs, sdate, edate, stime, etime,    &
                    hres, debug, ifile, ofile, int_default, vextrap,  &
                    cf_compliant, cordex_compliant,                   &
                    save_ccam_parameters

   integer :: kt, kdate, ktime, ierr, ieof, ntracers
   integer :: mins
   logical :: debug=.false.
   logical :: use_date, use_steps
   include 'revision.h'
   character(len=256) :: source, optionstring=""
   character(len=80)  :: basetime,calendar
   integer :: iyr, imon, iday, ihr, vid
   integer :: base_yr, base_mon, base_day, base_hr
   real :: time
   real, dimension(2) :: time_bnds
   logical :: skip
   type(input_var), dimension(:), pointer, contiguous :: varlist
   integer :: nvars
   type(hist_att), dimension(:), allocatable :: extra_atts, extra_temp
   integer :: veg_int
   real :: time_prev = 0.
   integer :: natts, catts, ival, xtype, attlen, i
   real :: rval
   character(len=80) :: cval, attname

   
#ifndef stacklimit
   ! For linux only
   call setstacklimit(-1)
#endif

   call MPI_Init(ierr)
   comm_world = MPI_COMM_WORLD
   call MPI_Comm_size(comm_world, nproc, ierr) ! Find number of processes
   call MPI_Comm_rank(comm_world, myid, ierr)  ! Find local processor id

   ! Start banner
   if ( myid==0 ) then
      write(6,*) "=============================================================================="
      write(6,*) "CCAM: Starting pcc2hist"
      write(6,*) "=============================================================================="
   end if
   
#ifdef usempi3
   call MPI_Comm_split_type(comm_world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierr) ! Per node communictor
   call MPI_Comm_size(node_comm, node_nproc, ierr) ! Find number of processes on node
   call MPI_Comm_rank(node_comm, node_myid, ierr)  ! Find local processor id on node
#endif

!  Initalise timing logs
   call log_off()
   call log_setup()
   call START_LOG(model_begin)

!  Check versions of library modules.
   call checkver ( "history", history_revision, 7, 4 )

!  Get the command line options
   longopts(1) = loption ( "interp", 1, 0)
   longopts(2) = loption ( "vextrap", 1, 0 )
   longopts(3) = loption ( "cf", 0, 0 )
   longopts(4) = loption ( "cordex", 0, 0 )
   ifile = ""
   ofile = ""
   cfile = ""
   do
      call getopt("c:d:hi:o:r:t:v",nopt,opt,optarg,longopts,longind)
      if ( opt == -1 ) exit  ! End of options
      select case ( char(opt) )
      case ( "c" )
         cfile = optarg
      case ( "d" )
         debug = .true.
      case ( "h" ) 
         call help(cc2hist_revision)
      case ( "i" ) 
         ifile = optarg
      case ( "o" ) 
         ofile = optarg
      case ( "r" ) 
         read(optarg,*) hres
      case ("v")
         if ( myid == 0 ) then
            write(6,*) "pcc2hist ", cc2hist_revision
         end if
      case ( char(0) )
         ! Long options that don't have a short form
         select case ( longind ) 
         case ( 1 ) 
            select case ( optarg )
            case ("nearest")
               int_default = int_nearest
               optionstring = optionstring(:len_trim(optionstring)) // " --interp=nearest"
            case ("linear")
               int_default = int_lin
               optionstring = optionstring(:len_trim(optionstring)) // " --interp=linear"
            case ( "none" )
               int_default = int_none
               optionstring = optionstring(:len_trim(optionstring)) // " --interp=none"
            case ( "tapm" )
               int_default = int_tapm
               optionstring = optionstring(:len_trim(optionstring)) // " --interp=tapm"
            case default
               write(6,*) "Expected nearest, linear or none for interp option"
               call finishbanner
               call MPI_ABORT(comm_world,-1,ierr)
            end select
         case ( 2 )
            select case ( optarg )
            case ("linear")
               vextrap = vextrap_lin
            case ("none")
               vextrap = vextrap_none
               optionstring = optionstring(:len_trim(optionstring)) // " --vextrap=none"
            case ("missing")
               vextrap = vextrap_missing
               optionstring = optionstring(:len_trim(optionstring)) // " --vextrap=missing"
            case default
               print*, "Expected linear, none or missing for vextrap option"
               call finishbanner
               call MPI_ABORT(comm_world,-1,ierr)
            end select
         case ( 3 )
            cf_compliant = .true.
         case ( 4 )
            cordex_compliant = .true.
         case default
            print*, "Unexpected result processing long options", longind
            call finishbanner
            call MPI_ABORT(comm_world,-1,ierr)
         end select
      case default
         if ( myid == 0 ) then
            print*, "Error unknown option "
         end if
         call usage()
      end select
   end do

   ! Read namelist - allows overwriting of command line options
   if ( trim(cfile) == "" ) cfile = "cc.nml"
   if ( myid==0 ) write(6,*) "reading ",trim(cfile)
   open(1,file=trim(cfile))
   read(1,input)   
   
   if ( vextrap == vextrap_missing .and. int_default == int_normal ) then
      print*, "For vextrap=missing option to work, must set interp to linear or nearest"
      call finishbanner
      call MPI_ABORT(comm_world,-1,ierr)
   end if

!  If filenames were not set as options look for them as arguments
   if ( len_trim(ifile) == 0 .and. len_trim(ofile) == 0 ) then
      if ( command_argument_count() /= nopt+1 ) then
         if ( myid == 0 ) then
            print*, "Missing filenames"
         end if
         call usage()
      else
         call getarg(nopt,ifile)
         call getarg(nopt+1,ofile)
      end if
   end if

   if ( len_trim(ifile) == 0 .or. len_trim(ofile) == 0 ) then
      print*, "Error setting input/output filenames"
      call finishbanner
      call MPI_ABORT(comm_world,-1,ierr)
   end if

   if ( use_plevs .and. use_meters ) then
      print *,"Cannot both use_plevs and use_meters together"
      call finishbanner
      call MPI_ABORT(comm_world,-1,ierr)
   end if

!  Check whether ndate and ntime have been set
   use_date = .false.
   if ( .not. ( sdate == -1 .and. edate == -1 .and. &
                stime == -1 .and. etime == -1 ) ) then
      use_date = .true.
   else if ( .not. ( ndate == -1 .and. ntime == -1 ) ) then
      use_date = .true.
      sdate = ndate
      edate = ndate
      stime = ntime
      etime = ntime
   end if
!  Set loop parameters to check all steps
   if ( use_date ) then
      kta=1
      ktb=huge(1)
      ktc=1
   end if
   if ( debug .and. use_date ) then
      print*, "Start", sdate, stime
      print*, "End  ", edate, etime
   end if

   use_steps = .not. use_date 

   call check_plevs
   call check_meters
   
   ! Open parallel input files  
   call paraopen(ifile, NF90_NOWRITE, ncid)
   ! Check that this file wasn't produced by cc2hist itself
   call check_cc2histfile()

   call initialise ( hres,                                               &
                     minlon, maxlon, dlon, minlat, maxlat, dlat,         &
                     dx, dy, lx, ly,                                     &
                     kdate, ktime, ntracers, ksoil, kice, debug,         &
                     nqg )
   call set_hbytes ( 4 )
   call set_htype ( "inst" )
   call set_hfreq ( 1 )

   minlev = max(1, minlev)
   maxlev = min(kl, maxlev)
!  Now find the index of the first sigma level above maxsig
   do k = minlev,maxlev
      if ( sig(k) <= maxsig ) exit
   end do
   minlev = k
!  Sigma level below minsig
   do k = maxlev,minlev,-1
      if ( sig(k) >= minsig ) exit
   end do
   maxlev = k
   nlev = maxlev - minlev + 1
   if ( nlev < 1 ) then
      print*, " Error: No vertical levels selected "
      print*, " Check minlev, maxlev, minsig, maxsig "
      stop
   end if
   if ( use_plevs .or. use_meters ) then
      nlev = nplevs
      minlev = 1
      maxlev = kl
   end if

   call get_var_list(varlist,nvars)

   call histnamelist(1)
   close(1)

!  openhist has to be called before any data is read because it allocates
!  the memory.
   if ( npanels == 5 ) then
      write(source,"(a,a,a,a)" ) "CSIRO conformal-cubic model. Input file: ",&
             trim(ifile(scan(ifile,'/',.TRUE.)+1:)), " Processed by cc2hist ", cc2hist_revision
   else
      write(source,"(a,a,a,a)" ) "CSIRO conformal-octagon model. Input file: ",& 
             trim(ifile(scan(ifile,'/',.TRUE.)+1:)), " Processed by cc2hist ", cc2hist_revision
      print *,"conformal-octagon no longer supported"
      stop
   end if
   if ( len_trim(optionstring) /= 0 ) then
      source = trim(source) // " Options:" // trim(optionstring)
   end if

!     If the input is a netcdf file with a time:units attribute then copy that.
   ierr = nf90_inq_varid (ncid, "time", vid )
   call check_ncerr(ierr, "Error getting time id")
   basetime = ""
   ierr = nf90_get_att(ncid, vid, "units", basetime)
   call check_ncerr(ierr, "Error getting time:units attribute")
   calendar = ""
   ierr = nf90_get_att(ncid, vid, "calendar", calendar)

   if ( cf_compliant .and. basetime(1:13) == "minutes since" ) then
      basetime = "days" // basetime(8:len_trim(basetime))
   end if
   if ( cf_compliant .and. basetime(1:13) == "seconds since" ) then
      basetime = "days" // basetime(8:len_trim(basetime))
   end if

   ierr = nf90_inquire(ncid, nAttributes=natts )
   call check_ncerr(ierr, "Error getting number of global attributes")
   allocate( extra_temp(natts) )
   catts = 0
   do i = 1,natts
      ierr = nf90_inq_attname(ncid, nf90_global, i, attname )
      call check_ncerr(ierr, "Error getting attribute name")
      if ( attname/="il_g"       .and. attname/="jl_g"        .and. &
           attname/="il"         .and. attname/="kl" ) then
         ierr = nf90_inquire_attribute(ncid, nf90_global, attname, xtype=xtype, len=attlen )
         call check_ncerr(ierr, "Error getting attribute type and len")
         select case (xtype)
            case( nf90_int )
               if ( attlen<=1 ) then             
                  catts = catts + 1  
                  ierr = nf90_get_att(ncid, nf90_global, attname, ival )
                  extra_temp(catts) = hist_att(attname, NF90_INT, 0., ival, "")
               end if
            case ( nf90_real )
               if ( attlen<=1 ) then                
                  catts = catts + 1                     
                  ierr = nf90_get_att(ncid, nf90_global, attname, rval )
                  extra_temp(catts) = hist_att(attname, NF90_REAL, rval, 0, "")
               end if
            case ( nf90_char )
               catts = catts + 1                     
               ierr = nf90_get_att(ncid, nf90_global, attname, cval )
               extra_temp(catts) = hist_att(attname, NF90_CHAR, 0., 0, cval)
         end select
      end if
   end do
   catts = catts + 2
   allocate ( extra_atts(catts) )
   extra_atts(1) = hist_att("il", NF90_INT, 0., il, "")
   extra_atts(2) = hist_att("kl", NF90_INT, 0., kl, "")
   extra_atts(3:catts) = extra_temp(1:catts-2)
   deallocate ( extra_temp )
   
   
   if ( calendar /= "" ) then
      if ( cf_compliant ) then
         if ( use_plevs ) then
            call openhist( il, jl, nlev, plevs(1:nplevs), ol, gosig, "_test",      &
                           hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                           source=source, histfilename=ofile, pressure=use_plevs,  &
                           extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil,        &
                           calendar=calendar )
         else if ( use_meters ) then
            call openhist( il, jl, nlev, mlevs(1:nplevs), ol, gosig, "_test",      &
                           hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                           source=source, histfilename=ofile, height=use_meters,   &
                           extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil,        &
                           calendar=calendar )   
         else
            call openhist( il, jl, nlev, sig(minlev:maxlev), ol, gosig, "_test",   &
                          hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis,  &
                          source=source, histfilename=ofile,                       &
                          extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil,         &
                          calendar=calendar )
         end if
      else
         if ( use_plevs ) then
            call openhist ( il, jl, nlev, plevs(1:nplevs), ol, gosig, "_test",      &
                            hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                            source=source, histfilename=ofile, pressure=use_plevs,  &
                            extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil,        &
                            calendar=calendar )
         else if ( use_meters ) then
            call openhist ( il, jl, nlev, mlevs(1:nplevs), ol, gosig, "_test",      &
                            hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                            source=source, histfilename=ofile, height=use_meters,   &
                            extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil,        &
                            calendar=calendar )
         else
            call openhist ( il, jl, nlev, sig(minlev:maxlev), ol, gosig, "_test",      &
                            hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis,    &
                            source=source, histfilename=ofile, extra_atts=extra_atts,  &
                            nsoil=ksoil, zsoil=zsoil, calendar=calendar )
         end if
      end if
   else
       if ( cf_compliant ) then
         if ( use_plevs ) then
            call openhist( il, jl, nlev, plevs(1:nplevs), ol, gosig, "_test",      &
                           hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                           source=source, histfilename=ofile, pressure=use_plevs,  &
                           extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil )
         else if ( use_meters ) then
            call openhist( il, jl, nlev, mlevs(1:nplevs), ol, gosig, "_test",      &
                           hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                           source=source, histfilename=ofile, height=use_meters,   &
                           extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil )
         else
            call openhist( il, jl, nlev, sig(minlev:maxlev), ol, gosig, "_test",   &
                          hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis,  &
                          source=source, histfilename=ofile,                       &
                          extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil )
         end if
      else
         if ( use_plevs ) then
            call openhist ( il, jl, nlev, plevs(1:nplevs), ol ,gosig, "_test",      &
                            hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                            source=source, histfilename=ofile, pressure=use_plevs,  &
                            extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil )
         else if ( use_meters ) then
            call openhist ( il, jl, nlev, mlevs(1:nplevs), ol, gosig, "_test",      &
                            hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis, &
                            source=source, histfilename=ofile, height=use_meters,   &
                            extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil )
         else
            call openhist ( il, jl, nlev, sig(minlev:maxlev), ol, gosig, "_test",       &
                            hlon, hlat, basetime, year=1, nxout=nxhis, nyout=nyhis,     &
                            source=source, histfilename=ofile, extra_atts=extra_atts,   &
                            nsoil=ksoil, zsoil=zsoil )
         end if
      end if  
   end if

!  needfld calls are only valid after openhist
   call final_init(varlist,nvars)

   if ( needfld("zg") .or. use_meters ) then
      call initheight( kl, sig )
   end if

!  If ktc is still -1 and ndate and ntime aren't set then process all fields
   if ( ktc == -1 ) then
      call getstep(ktc)
   end if

   if ( ktc == 0 ) then
      write(6,*) "ERROR: ktc must not equal zero"
      stop
   end if
   
   call log_on()
   call START_LOG(timeloop_begin)
   timeloop: do kt=kta,ktb,ktc
      if ( debug ) then
         print*, "KT", kt
      end if
      do ! Loop through until time found
         call getdate ( kdate, ktime, ieof )
         if ( ieof /= 0 ) then  ! END OF FILE
            exit timeloop
         end if
         if ( debug ) then
            print*, "KTAU, KT", ktau, kt
         end if
         
         if ( use_steps ) then
!           Check whether the correct time has been found
            if ( ktau < kt ) then
               skip = .true.
            else if ( ktau == kt ) then
               if ( debug ) then
                  print*, " Saving (ktau=kt) "
               end if
               skip = .false.
            else
               print *,"WARN: Incorrect time-step ktc" 
               skip = .false. 
            end if

!        Check for matching date and time
         else if ( use_date ) then
!           Check if reached startdate
            if ( kdate < sdate .or. &
                 ( kdate == sdate .and. ktime < stime ) ) then
               if ( debug ) then
                  print*, " before start date"
               end if
!              Go to next timeloop iteration, so kt advances.
               skip = .true. 
!           Check if past enddate
            else if ( kdate > edate .or. &
                      ( kdate == edate .and. ktime > etime ) ) then
               if ( debug ) then
                  print*, " past end date"
               end if
!              Finished
               exit timeloop
            else
               if ( debug ) then
                  print*, " Saving (ndate, ntime) "
               end if
               skip = .false. 
               exit
            end if
         end if

         if ( kt < ktau ) then
           if ( myid == 0 ) then             
              print*, "WARNING: Searching for step", kt
           end if
         end if

         ! calculate mins since start of the year
         call convert_date(mins, kdate, ktime, kt, basetime )
         
         call infile(varlist, nvars, skip, mins)

         if ( myid == 0 ) then
            if ( skip ) then
               write(*,"(a)") " - Skipped"
            else
               write(*,"(a)") " - Saved"
            end if
         end if

         !if ( myid==0 .and. nproc_orig/=nproc ) then
         !   write(6,'(x,a,i0,a,i0,a)') "WARNING: Number of processors(",nproc_orig,&
         !                              ") is not a factor of the number of files(",pnproc,")"
         !   write(6,'(x,a,i0)') "WARNING: Using pcc2hist with the following number of processes: ",nproc
         !end if

         if ( .not. skip ) then
            exit
         end if
         if ( use_date ) then
!           Go to next timeloop iteration, so kt advances.
            cycle timeloop
         end if

      end do

      if ( skip ) then
         cycle
      end if

!     Time in hours
      if  ( use_steps ) then
         time = real(ktau)
      else
         iyr = kdate/10000
         imon = modulo(kdate,10000) / 100
         iday = modulo(kdate,100)
         ihr = ktime/100
!        Need to extend to work over more than one month.
         time = (iday - base_day)*24 + ihr-base_hr
      end if
      if ( cf_compliant ) then 
         time = time/1440. ! Days
         ! History data at time t, really represents values over the preceeding
         ! period in the case of accumulated or average fields.
         ! Should this be the history file time increment or the increment
         ! set for cc2hist?
         time_bnds = (/time_prev,time/)
         call writehist ( ktau, interp=ints, time=time, time_bnds=time_bnds )
         time_prev = time
      else
         call writehist ( ktau, interp=ints, time=time)
      end if

   end do timeloop
   
   call END_LOG(timeloop_end)
   call log_off()

   call writehist(ktau, interp=ints, time=time, endofrun=.true. )
   if ( myid == 0 ) call closehist
   call paraclose

   call END_LOG(model_end)
#ifdef simple_timer
   call simple_timer_finalize
#endif

   ! Complete
   if ( myid==0 ) then
      write(6,*) "------------------------------------------------------------------------------"
      write(6,*) "CCAM: pcc2hist completed successfully"
      call finishbanner
   end if

   call cc2hist_work_close
   
   call MPI_Finalize(ierr)

end program cc2hist
    
subroutine finishbanner

implicit none

! End banner
write(6,*) "=============================================================================="
write(6,*) "CCAM: Finished pcc2hist"
write(6,*) "=============================================================================="

return
end    
