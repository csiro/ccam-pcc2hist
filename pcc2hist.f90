! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2025 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
!  output file.
    
!  Modified by MJT to use MPI when reading parallel input files (pcc2hist).
!  Further optimised by Paul Ryan and David Benn.

! Preprocessor directives:
!   usempimod    - Use fortran90 interface for MPI
!   share_ifullg - Allow shared memory on a node with MPI (usempi3 is always active)

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
#ifdef usempimod
   use mpi
#endif   

   implicit none

#ifndef usempimod   
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
   type(loption), dimension(7) :: longopts
   integer :: longind
   integer :: kta=0, ktb=999999, ktc=-1, ndate=-1, ntime=-1, k
   integer :: sdate=-1, edate=-1, stime=-1, etime=-1
   character(len=10) :: name
   character(len=80) :: longname
   real :: minsig = 0., maxsig = 1.0
   logical :: debug=.false.

   namelist /input/ kta, ktb, ktc, ndate, ntime,                      &
                    minlon, maxlon, dlon, minlat, maxlat, dlat,       &
                    dx, dy, lx, ly,                                   &
                    minlev, maxlev, minsig, maxsig, use_plevs, plevs, &
                    use_meters, mlevs, use_depth, dlevs, use_theta,   &
                    tlevs, use_pvort, vlevs, sdate,                   &
                    edate, stime, etime, hres, debug, ifile, ofile,   &
                    int_default, vextrap, cf_compliant,               &
                    cordex_compliant, save_ccam_parameters,           &
                    ran_compliant, areps_compliant, safe_max,         &
                    fao_potev

   include 'revision.h'
   
   integer :: kt, kdate, ktime, ierr, ieof, ntracers
   integer :: mins, lerr
   integer :: iyr, imon, iday, ihr, vid
   integer :: base_yr, base_mon, base_day, base_hr
   integer :: nvars
   integer :: veg_int, ncount
   integer :: natts, catts, ival, xtype, attlen, i
   integer, dimension(16) :: dumi
   logical :: use_date, use_steps
   logical :: skip
   logical, dimension(12) :: duml
   character(len=256) :: source, optionstring=""
   character(len=80)  :: basetime,calendar
   character(len=120) :: attname
   character(len=256) :: cval
   real :: time
   real :: time_prev = 0.
   real :: rval, dtime
   real, dimension(2) :: time_bnds
   real, dimension(:), allocatable :: xlevs, oxlevs
   real, dimension(11) :: dumr
   type(input_var), dimension(:), pointer, contiguous :: varlist
   type(hist_att), dimension(:), allocatable :: extra_atts, extra_temp

   
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
   
#ifdef share_ifullg
   call MPI_Comm_split_type(comm_world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierr) ! Per node communicator
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
   longopts(5) = loption ( "multioutput", 0, 0 )
   longopts(6) = loption ( "ran", 0, 0 )
   longopts(7) = loption ( "areps", 0, 0 )
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
         case ( 5 )
            single_output = .false.
         case ( 6 )
            ran_compliant = .true. ! depreciated
         case ( 7 )
            areps_compliant = .true.
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
   if ( myid == 0 ) then
      write(6,*) "reading ",trim(cfile)
      open(1,file=trim(cfile))
      read(1,input)
   end if   
   
   call START_LOG(mpibcast_begin)
   dumi(1) = kta
   dumi(2) = ktb
   dumi(3) = ktc
   dumi(4) = ndate
   dumi(5) = ntime
   dumi(6) = lx
   dumi(7) = ly
   dumi(8) = minlev
   dumi(9) = maxlev
   dumi(10) = sdate
   dumi(11) = edate
   dumi(12) = stime
   dumi(13) = etime
   dumi(14) = int_default
   dumi(15) = vextrap
   dumi(16) = safe_max
   call MPI_Bcast( dumi(:), size(dumi), MPI_INTEGER, 0, COMM_WORLD, lerr )
   kta = dumi(1)
   ktb = dumi(2)
   ktc = dumi(3)
   ndate = dumi(4)
   ntime = dumi(5)
   lx = dumi(6)
   ly = dumi(7)
   minlev = dumi(8)
   maxlev = dumi(9)
   sdate = dumi(10)
   edate = dumi(11)
   stime = dumi(12)
   etime = dumi(13)
   int_default = dumi(14)
   vextrap = dumi(15)
   safe_max = dumi(16)
   dumr(1) = minlon
   dumr(2) = maxlon
   dumr(3) = dlon
   dumr(4) = minlat
   dumr(5) = maxlat
   dumr(6) = dlat
   dumr(7) = dx
   dumr(8) = dy
   dumr(9) = minsig
   dumr(10) = maxsig
   dumr(11) = hres
   call MPI_Bcast( dumr(:), size(dumr), MPI_REAL, 0, COMM_WORLD, lerr )
   minlon = dumr(1)
   maxlon = dumr(2)
   dlon = dumr(3)
   minlat = dumr(4)
   maxlat = dumr(5)
   dlat = dumr(6)
   dx = dumr(7)
   dy = dumr(8)
   minsig = dumr(9)
   maxsig = dumr(10)
   hres = dumr(11)
   call MPI_Bcast( plevs, size(plevs), MPI_REAL, 0, COMM_WORLD, lerr )
   call MPI_Bcast( mlevs, size(mlevs), MPI_REAL, 0, COMM_WORLD, lerr )
   call MPI_Bcast( dlevs, size(dlevs), MPI_REAL, 0, COMM_WORLD, lerr )
   call MPI_Bcast( tlevs, size(tlevs), MPI_REAL, 0, COMM_WORLD, lerr )
   call MPI_Bcast( vlevs, size(vlevs), MPI_REAL, 0, COMM_WORLD, lerr )
   duml(1) = use_plevs
   duml(2) = use_meters
   duml(3) = use_depth
   duml(4) = use_theta
   duml(5) = use_pvort
   duml(6) = debug
   duml(7) = cf_compliant
   duml(8) = cordex_compliant
   duml(9) = save_ccam_parameters
   duml(10) = ran_compliant
   duml(11) = areps_compliant
   duml(12) = fao_potev
   call MPI_Bcast( duml(:), size(duml), MPI_LOGICAL, 0, COMM_WORLD, lerr )
   use_plevs = duml(1)
   use_meters = duml(2)
   use_depth = duml(3)
   use_theta = duml(4)
   use_pvort = duml(5)
   debug = duml(6)
   cf_compliant = duml(7)
   cordex_compliant = duml(8)
   save_ccam_parameters = duml(9)
   ran_compliant = duml(10)
   areps_compliant = duml(11)
   fao_potev = duml(12)
   call ccmpi_bcast(ifile,0,COMM_WORLD)
   call ccmpi_bcast(ofile,0,COMM_WORLD)
   call END_LOG(mpibcast_end)
   
   
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

   ncount = 0
   if ( use_plevs ) ncount = ncount + 1
   if ( use_meters ) ncount = ncount + 1
   if ( use_theta ) ncount = ncount + 1
   if ( use_pvort ) ncount = ncount + 1
   if ( ncount>1 ) then
      print *,"Cannot use use_plevs, use_meters, use_theta and use_port together"
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
   call check_depth
   call check_theta
   call check_pvort
   
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
   if ( use_plevs .or. use_meters .or. use_theta .or. use_pvort ) then
      nlev = nplevs
      minlev = 1
      maxlev = kl
   end if
   onlev = ol
   if ( use_depth ) then
      onlev = onplevs
      if ( ol==0 .and. onlev>0 ) then
         print *," Error: No ocean levels avaliable in input file"
         print *," Try using use_depth = F"
         stop
      end if   
   end if

   call get_var_list(varlist,nvars)

   call histnamelist(1)
   if ( myid == 0 ) then
      close(1)
   end if

!  openhist has to be called before any data is read because it allocates
!  the memory.
   write(source,"(a,a,a,a)" ) "CSIRO conformal-cubic model. Input file: ",&
          trim(ifile(scan(ifile,'/',.TRUE.)+1:)), " Processed by cc2hist ", cc2hist_revision
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
   if ( areps_compliant .and. basetime(1:13) == "minutes since" ) then
      basetime = "hours" // basetime(8:len_trim(basetime))
   end if

   ierr = nf90_inquire(ncid, nAttributes=natts )
   call check_ncerr(ierr, "Error getting number of global attributes")
   allocate( extra_temp(natts) )
   catts = 0
   do i = 1,natts
      ierr = nf90_inq_attname(ncid, nf90_global, i, attname )
      call check_ncerr(ierr, "Error getting attribute name")
      if ( attname/="il_g"       .and. attname/="jl_g"        .and. &
           attname/="il"         .and. attname/="kl"          .and. &
           attname/="procmode"   .and. attname/="decomp"      .and. &
           attname/="nproc"      .and. attname/="nrun" ) then
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

   allocate( xlevs(nlev), oxlevs(onlev) ) ! use a pointer here?
   if ( use_plevs ) then
      xlevs(1:nlev) = plevs(1:nlev)
   else if ( use_meters ) then
      xlevs(1:nlev) = mlevs(1:nlev)
   else if ( use_theta ) then
      xlevs(1:nlev) = tlevs(1:nlev)
   else if ( use_pvort ) then
      xlevs(1:nlev) = vlevs(1:nlev)
   else
      xlevs(1:nlev) = sig(minlev:maxlev)
   end if
   if ( use_depth ) then
      oxlevs(1:onlev) = dlevs(1:onlev)
   else
      oxlevs(1:onlev) = gosig(1:onlev) 
   end if
   
   if ( calendar /= "" ) then
      call openhist( il, jl, nlev, xlevs(1:nlev), onlev, cptch, cchrt, oxlevs(1:onlev), &
                     "_test", hlon, hlat, basetime, year=1, nxout=nxhis,     &
                     nyout=nyhis, source=source, histfilename=ofile,         &
                     pressure=use_plevs, height=use_meters, theta=use_theta, &
                     pvort=use_pvort, depth=use_depth,                       &
                     extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil,        &
                     calendar=calendar )
   else
      call openhist( il, jl, nlev, xlevs(1:nlev), onlev, cptch, cchrt, oxlevs(1:onlev), &
                     "_test", hlon, hlat, basetime, year=1, nxout=nxhis,     &
                     nyout=nyhis, source=source, histfilename=ofile,         &
                     pressure=use_plevs, height=use_meters, theta=use_theta, &
                     pvort=use_pvort, depth=use_depth,                       &
                     extra_atts=extra_atts, nsoil=ksoil, zsoil=zsoil )
   end if

   deallocate( xlevs, oxlevs )
   
!  needfld calls are only valid after openhist
   call final_init( varlist, nvars )

   call initheight( kl, sig )

!  If ktc is still -1 and ndate and ntime aren't set then process all fields
   call getstep( kta, ktc )

   ! calculate time interval
   call getdtime( dtime, ktc )
   if ( cf_compliant ) then
      dtime = dtime/1440. ! Days
   else if ( areps_compliant ) then
      dtime = dtime/60.   ! Hours
   end if
   
   call log_on()
   call START_LOG(timeloop_begin)
   timeloop: do kt=kta,ktb,ktc
      if ( debug ) then
         print*, "KT", kt
      end if
      do ! Loop through until time found
         call getdate( kdate, ktime, ieof ) ! writes date to output log/screen
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
              print*, "WARNING: Searching for step kt=", kt
              print*, "         but found step ktau=",ktau
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

         if ( .not. skip ) then
            exit
         end if
         if ( use_date ) then
!           Go to next timeloop iteration, so kt advances.
            cycle timeloop
         end if

      end do ! loop until found

      if ( skip ) then
         cycle
      end if

!     Time in mins (?)
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
      ! History data at time t, really represents values over the preceeding
      ! period in the case of accumulated or average fields.
      ! Should this be the history file time increment or the increment
      ! set for cc2hist?
      !
      ! dtime is the time interval and is used to calculate the time-centre for
      ! CORDEX DRS output.
      if ( cf_compliant ) then 
         time = time/1440. ! Days
         time_bnds = (/time_prev,time/)
         call writehist ( ktau, interp=ints, time=time, time_bnds=time_bnds, dtime=dtime )
         time_prev = time
      else if ( areps_compliant ) then
         time = time/60. ! Hours
         call writehist ( ktau, interp=ints, time=time, dtime=dtime ) 
      else if ( cordex_compliant ) then
         time_bnds = (/time_prev,time/)   
         call writehist ( ktau, interp=ints, time=time, time_bnds=time_bnds, dtime=dtime )   
         time_prev = time   
      else
         call writehist ( ktau, interp=ints, time=time, dtime=dtime )
      end if

   end do timeloop
   
   call END_LOG(timeloop_end)
   call log_off()

   call writehist( ktau, interp=ints, time=time, dtime=dtime, endofrun=.true. )
   call closehist
   call paraclose

   call END_LOG(model_end)
   call simple_timer_finalize

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

end subroutine finishbanner
