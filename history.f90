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
    
! This module provides a new scheme for saving model history data.
! The idea is to hide the arrays used for the accumulation. All 
! manipulation can be done giving just the name of the variable.
! The idea is based on the NCAR CCM3 but it's a lot more convenient in
! Fortran 90.

! Also useful in diagnostic programs, anything that has to write a time 
! series of lat/lon data to a netcdf file.

module history

   use netcdf_m
   use ncutils_m, only : check_ncerr
   use utils_m, only : fpequal
   implicit none

!  Make everything private so that the internal variables of the history 
!  module can only be used by the appropriate routines.
   private   

!  Only these routine names need to be public
   public :: savehist, openhist, closehist, writehist, addfld, &
             histnamelist, needfld, clearhist, hstring

!  Routine that provide access to the namelist control variables.
!  Access is only possible through a routine so that the underlying data
!  structure may be changed without affecting user programs.
   public :: set_htype,  set_hfreq, set_hnames, set_hbytes
   
   public :: spval

   !public :: set_missval

!  Private internal routines
   private :: bindex_hname, sortlist, create_ncvar,                     &
              create_ncfile,                                            &
              hashkey, savehist_work,                                   &
              gsavehist2D, gsavehist3D, gsavehist4D
   private :: fill_cc1
   !private :: create_oldfile
   
   public :: ccmpi_bcast
   interface ccmpi_bcast
      module procedure ccmpi_bcast1s, ccmpi_bcast2s
   end interface   

   character(len=50), public, parameter :: &
        history_revision = "$Revision: 7.10 $"

!  This gets a 64 bit integer
!   integer, private, parameter :: bigint = selected_int_kind(18)
!  On SX-6 range is smaller, though full range is available for addition
!  and subtraction which is all that's required for the key.
   integer, private, parameter :: bigint = selected_int_kind(14)

!  Debugging level. 0 for no debugging, larger values for more detail
   integer, public, save :: hist_debug = 0
!  Debug history at a point
   integer, public, save :: ihdb=1, jhdb=1, khdb=1

!  This variable controls whether to use the old Mark 2 names for history
!  variables or whether to use AMIP2 names. This only affects the variable
!  names in the history file. All calls to savefld must use the Mark 2 name.
!  Note that this hasn't yet been implemented.
   logical, public, save :: amipnames = .false.

!  Maximum length of names of history variables
   integer, parameter :: MAX_NAMELEN = 60

!  Maximum length of string used for key generation. Names must be unique
!  within this length.
   integer, parameter :: MAX_KEYLEN = 10

!  Maximum number of fields allowed
!  Note that this maximum is only used for dimensioning header arrays.
!  The actual data arrays are allocatable, so it doesn't use much space.
!  It can't be made truly dynamic. The best that could be done is to 
!  re-allocate a larger array if required.
   integer, parameter :: nfmax = 20000

!  Initialise hnames so that set_hnames can find the end of the list to append.
   character(len=MAX_NAMELEN), dimension(nfmax), save ::  &
         hnames = "", xnames = ""

!  Integer array of keys corresponding to hnames to make searching faster.
   integer, parameter :: MAX_KEYINDEX = 6 !>=MAX_NAMELEN/MAX_KEYLEN
   integer (bigint), dimension(MAX_KEYINDEX,nfmax), save :: inames

!  Frequency of writing history files. This may be in either steps or minutes
!  depending on the argument of the calling routine. 
   integer, save :: hfreq = 0

!  Number of bytes to use for history variables
   integer, save :: hbytes = 4

!  Default averaging type for file
   integer, save :: ihtype

!  Average time appropriate for time averaged files
   real, save :: avetime
   real, dimension(2), save :: avetime_bnds
   integer, save :: timecount

!  To allow old format files to be written in another directory
   character(len=200), save :: oldprefix

!  Fill value for tsea   
   real, parameter :: spval   = 999.   
   
!  Record for each field in a history file
   type hinfo
      character(len=MAX_NAMELEN)   :: name      ! Mnemonic name
      character(len=MAX_NAMELEN)   :: amip_name ! Mnemonic name (AMIP standard)
      character(len=80)  :: long_name           ! Full name
      character(len=80)  :: std_name            ! IPCC/CF standard name
      character(len=30)  :: units
      real               :: valid_min
      real               :: valid_max
      integer            :: nlevels
!     Controls whether variable is on the "standard" list of monthly means
      logical            :: std
!     Controls whether variable is on the "RAN" list of variables
      logical            :: ran, areps
!     Controls whether variable is on the "tracer" list of variables
      logical            :: tracer
!     Controls whether variable is on the t?_pop list of variables
      integer            :: tn
!     Scale factor for output, e.g. to convert rain from mm/step to mm/day.
      real               :: output_scale 
      logical            :: used
      integer            :: vid  ! netCDF variable identifier
      integer            :: ptr  ! Pointer to position in histarray
      integer            :: count      ! Used for scaling
      integer            :: ave_type   ! Type of averaging
      integer            :: ncid ! netCDF file identifier
      integer            :: procid
      real               :: addoff
      real               :: scalef
      integer            :: int_type   ! Type of interpolation
!     Flag to force variable to be handled as multi-level even if nlevels=1.
      logical            :: multilev
!     For multilevel land surface variables
      logical            :: soil
!     For multilevel ocean variables
      logical            :: water
!     For 2d pop variables
      logical            :: pop2d
!     For 3d pop variables
      logical            :: pop3d
!     For 4d pop variables
      logical            :: pop4d
      ! For CF coordinate attribute
      real               :: coord_height
      character(len=80)  :: coord_name
      character(len=80)  :: coord_stdname
      character(len=80)  :: coord_units
      character(len=80)  :: coord_positive
      ! cell_methods appropriate for variable before any history processing is done
      character(len=30)  :: cell_methods
      ! for daily data
      logical            :: daily
      ! for 6hr data
      logical            :: sixhr
      ! instantaneous
      logical            :: instant
      ! all positive values
      logical            :: all_positive
      ! fill required
      logical            :: fill
   end type hinfo

!  For adding extra global attributes
   type, public :: hist_att
      character(len=nf90_max_name)   :: name
      integer :: type
      real :: rval
      integer :: ival
      ! This length isn't really related to the max_name length but this will do.
      character(len=256)   :: cval
   end type hist_att

!  Type to hold dimension IDs and dimension variable IDs
!  Initialise to -1 to catch errors better
   type, private :: dimarray
      integer :: x = -1
      integer :: y = -1
      integer :: z = -1
      integer :: oz = -1
      integer :: t = -1
      integer :: zsoil = -1
      integer :: cptch = -1
      integer :: cchrt = -1
      ! For bounds 
      integer :: b = -1
      integer :: x_b = -1
      integer :: y_b = -1
      integer :: z_b = -1
      integer :: t_b = -1
      integer :: zsoil_b = -1
   end type dimarray
   
   ! info for output file dimensions and time-steps
   type, private :: hfile
      integer :: id
      integer :: nlevels
      logical :: day
      logical :: sixhr
      logical :: fix
      logical :: inst
      logical :: ocean
      logical :: soil
      logical :: cablepatch
      logical :: cablecohort
      integer :: ncoords
      real, dimension(:), allocatable :: coord_height
      character(len=80), dimension(:), allocatable :: coord_name
   end type hfile

!  Array for accumulating history data.
   real, dimension(:,:,:), allocatable :: histarray

   type (hinfo), dimension(nfmax) :: histinfo 
   type (dimarray), dimension(nfmax) :: histdims
!  Number of history records written
   integer :: histset
   integer :: histset_daily, histset_6hr

!   netCDF file ID of history file
   type(hfile), dimension(:), allocatable, save :: histfile
   type(dimarray), dimension(:), allocatable, save :: histdimvars

!  Total number of fields defined (not necessarily used).
   integer, save :: totflds = 0

   interface savehist
      module procedure gsavehist2D, gsavehist3D, gsavehist4D
   end interface
    
!  Define parameters to identify the processing of history data: ave, max, min
!  or instantaneous value. hist_fixed is for fields like the surface height
!  which don't change but are still useful to have in a history file.

   integer, parameter :: hist_ave=1, hist_max=2, hist_min=3, hist_inst=4, &
                         hist_fixed=5

!  Valid range for 16 bit values
   integer, parameter :: vmin=-32500, vmax=32500

   !real :: missing_value_cordex = 1.00000002004e+20
   real :: missing_value_cordex = nf90_fill_float
   character(len=20), save :: filesuffix

!  Output array size (may be different to model resolution and accumulation
!  array size).
   integer, private, save :: nxhis, nyhis

!  Flag to check whether interpolation is required before output
   logical, private, save :: require_interp = .false.

!  Type of interpolation
   integer, public :: int_default=0

   logical, public :: cf_compliant = .false.
   
   logical, public :: cordex_compliant = .false.
   
   logical, public :: single_output = .true.
   
   logical, public :: ran_compliant = .false. ! depreciated
   
   logical, public :: areps_compliant = .false.
   
! Save CCAM parameters
   logical, public :: save_ccam_parameters = .true.

!  Working arrays
   integer, private, save :: nx_g, ny_g
   real, dimension(:,:,:,:), allocatable, save, private :: hist_a
   
! chunk size
   integer, public, save :: chunk_grid = -1 ! recommend 24
   
! flush output file buffers
   integer, public, save :: safe_max = 3

contains

!-------------------------------------------------------------------

   subroutine hstring(hstr)

      ! Create a string containing the username, machine name and date to 
      ! identify the creation of the history file
   
      character(len=*), intent(out) :: hstr
      !character(len=100) :: logname, hostname
      integer, dimension(8) :: idatetime
      character(len=1000) :: str
   
      !logname = ""
      !call get_environment_variable("LOGNAME",logname)
      !hostname = ""
      !call get_environment_variable("HOSTNAME",hostname)
      call date_and_time(values=idatetime)
      
      ! Avoid overflow in the output string by writing to a long internal string
      str = ""
      write(str,'("Created on ", i4.4,"-", i2.2,"-",i2.2," ",i2.2,":",i2.2,":", i2.2)') &
           idatetime(1:3), idatetime(5:7)
      hstr = str ! Copy truncates if necessary
      
   end subroutine hstring

   subroutine histnamelist ( un, iofatal )
!
!     Read the namelist variables controlling the history.
!     This routine uses local arrays for names etc, so it can control
!     the overwriting of the module variables.
!
#ifdef usempimod
      use mpi
#endif
      use logging_m
      use mpidata_m
#ifndef usempimod
      include 'mpif.h'
#endif      
      integer, intent(in) :: un      ! Unit no to read namelist
!     Controls whether errors in reading the namelist should be fatal.
      logical, intent(in), optional :: iofatal 
      character(len=MAX_NAMELEN) :: htype
      integer :: freq, bytes
      integer :: lerr
      integer, dimension(9) :: dumi
      character(len=MAX_NAMELEN), dimension(nfmax) ::  &
         names, namesx

      namelist / histnl / &
           htype,  &  ! Average type for this file
           hfreq,  &  ! Frequency of writing
           hnames, &  ! Names to include
           xnames, &  ! Names to exclude
           hbytes,  &
           hist_debug, & 
           ihdb, jhdb, khdb,  &
           amipnames, &
           chunk_grid
      integer :: i, j, ierr
      logical :: fatal

      ! Save the current value of the namelist variables and then
      ! reset them
      htype = ""
      freq = hfreq; hfreq = 0
      bytes = hbytes; hbytes = 4
      names = hnames; hnames = ""
      namesx = xnames; xnames = ""
      if ( myid == 0 ) then
        read(un,histnl,iostat=ierr)
      end if

      call START_LOG(mpibcast_begin)
      dumi(1) = ierr
      dumi(2) = hfreq
      dumi(3) = hbytes
      dumi(4) = ihdb
      dumi(5) = jhdb
      dumi(6) = khdb
      dumi(7) = chunk_grid
      dumi(8) = hist_debug
      if ( amipnames ) then
        dumi(9) = 1
      else
        dumi(9) = 0
      end if
      call MPI_BCAST( dumi(:), size(dumi), MPI_INTEGER, 0, COMM_WORLD, lerr )
      ierr = dumi(1)
      hfreq = dumi(2)
      hbytes = dumi(3)
      ihdb = dumi(4)
      jhdb = dumi(5)
      khdb = dumi(6)
      chunk_grid = dumi(7)
      hist_debug = dumi(8)
      amipnames = dumi(9)==1
      call ccmpi_bcast(htype,0,COMM_WORLD)
      call ccmpi_bcast(hnames,0,COMM_WORLD)
      call ccmpi_bcast(xnames,0,COMM_WORLD)
      call END_LOG(mpibcast_end)
      
      if ( ierr /= 0 ) then
         if ( present(iofatal) ) then
            ! Error values > 0 are due to problems with the namelist section
            ! and should be fatal. Errors < 0 might just be that the section
            ! is missing and should not be fatal.
            fatal = iofatal .or. ierr > 0
         else
            fatal = .true.
         end if
         if (fatal) then
            print*, "Error reading history namelist histnl", ierr
            stop
         else
            ! In this case the read failed, so restore the original value
            ! of the variables
            hfreq = freq
            hbytes = bytes
            hnames = names
            xnames = namesx
            return
         end if
      end if

      select case (htype)
      case ("ave")
         ihtype = hist_ave
      case ("max")
         ihtype = hist_max
      case ("min")
         ihtype = hist_min
      case ("inst")
         ihtype = hist_inst
      case ("fixed")
         ihtype = hist_fixed
      case ("")
         ihtype = hist_ave
      case default
         print*, " Error: average type set incorrectly in namelist ", &
                 trim(htype)
         stop
      end select

      if ( hist_debug > 0 ) then
         do i=1,nfmax
            if ( len_trim(hnames(i)) /= 0 ) then
               write(*,"(i3,1x,a)") i, hnames(i)
            else
               exit
            end if
         end do
      end if

   end subroutine histnamelist

   subroutine set_htype ( htype )
!     Routine to allow setting htype directly without namelist.
      character(len=*), intent(in) :: htype

      select case (htype)
      case ("ave")
         ihtype = hist_ave
      case ("max")
         ihtype = hist_max
      case ("min")
         ihtype = hist_min
      case ("inst")
         ihtype = hist_inst
      case ("fixed")
         ihtype = hist_fixed
      case ("")
         ihtype = hist_ave
      case default
         print*, " Error: average type set incorrectly in set_htype ", &
              trim(htype)
         stop
      end select

   end subroutine set_htype
   

   subroutine set_hnames ( name )
!     Routine to allow adding a variable to hnames list directly without 
!     namelist.
      character(len=*), intent(in) :: name
      integer :: i

!     Append this to list of hnames by searching for the first empty entry.
!     This isn't very efficient but this routine is only used at startup.
      do i = 1,nfmax
         if ( len_trim((hnames(i))) == 0 ) then
            hnames(i) = name
            exit
         end if
      end do
      if ( i > nfmax ) then
         print*, " Error, nfmax exceeded in set_hnames "
         stop
      end if

   end subroutine set_hnames

   subroutine set_hbytes ( nbyte )
!     Routine to allow setting hbytes directly without namelist.
      integer, intent(in) :: nbyte
      hbytes = nbyte
   end subroutine set_hbytes

   subroutine set_hfreq ( freq )
!     Routine to allow setting hfreq directly without namelist.
      integer, intent(in) :: freq
      hfreq = freq
   end subroutine set_hfreq

!---------------------------------------------------------------------------
   subroutine addfld(name, long_name, units, valid_min, valid_max,    &
                     nlevels, amip_name, ave_type, std, output_scale, &
                     int_type, multilev, std_name, soil, water,       &
                     pop2d, pop3d, pop4d, coord_height, coord_name,   &
                     coord_stdname, coord_units, coord_positive,      &
                     cell_methods,                                    &
                     ran_type, areps_type, tn_type, daily, sixhr,     &
                     instant, all_positive, tracer_type, fill )
!
!     Add a field to the master list of fields that may be saved.
!
      character(len=*), intent(in)    :: name       ! Mnemonic name
      character(len=*), intent(in)    :: long_name  ! Full name
      character(len=*), intent(in)    :: units
      real, intent(in)                :: valid_min
      real, intent(in)                :: valid_max
      integer, intent(in)             :: nlevels
      character(len=*), intent(in), optional :: amip_name
      character(len=*), intent(in), optional :: ave_type   ! Full name
      logical, intent(in), optional          :: std
      real, intent(in), optional             :: output_scale
      integer, intent(in), optional   :: int_type
      logical, intent(in), optional   :: multilev
      character(len=*), intent(in), optional :: std_name
      logical, intent(in), optional   :: soil
      logical, intent(in), optional   :: water
      logical, intent(in), optional   :: pop2d
      logical, intent(in), optional   :: pop3d
      logical, intent(in), optional   :: pop4d
      real, intent(in), optional :: coord_height
      character(len=*), intent(in), optional :: coord_name
      character(len=*), intent(in), optional :: coord_stdname
      character(len=*), intent(in), optional :: coord_units
      character(len=*), intent(in), optional :: coord_positive
      character(len=*), intent(in), optional :: cell_methods
      logical, intent(in), optional   :: ran_type
      logical, intent(in), optional   :: areps_type
      logical, intent(in), optional   :: tracer_type
      integer, intent(in), optional   :: tn_type
      logical, intent(in), optional   :: daily
      logical, intent(in), optional   :: sixhr
      logical, intent(in), optional   :: instant
      logical, intent(in), optional   :: all_positive
      logical, intent(in), optional   :: fill

!     Local variables corresponding to the optional arguments
      integer :: atype, ltn
      logical :: lstd, lran, lareps, ltracer
      character(len=MAX_NAMELEN) :: aname
      real    :: scale

!     Should check the lengths of the string arguments

      totflds = totflds + 1
      if ( totflds > nfmax ) then
         print*, ' Increase nfmax in history module '
         stop
      end if

      if ( present(ave_type) ) then
         select case (adjustl(ave_type))
         case ("ave")
            atype = hist_ave
         !case ("oave")
         !   atype = hist_oave
         case ("max")
            atype = hist_max
         case ("min")
            atype = hist_min
         case ("inst")
            atype = hist_inst
         case ("fixed")
            atype = hist_fixed
         case default
            print*, " Error: average type set incorrectly"
            print*, " Variable ", name, " ave type ", ave_type
            stop
         end select
      else
!        This is essentially a missing value which will be overridden later.
         atype = 0
      end if
      
      if ( present(amip_name) ) then
         aname = amip_name
      else
         aname = name
      end if
      if ( present(std) ) then
         lstd = std
      else
         lstd = .false.
      end if
      if ( present(ran_type) ) then
         lran = ran_type  
      else
         lran = .false.  
      end if
      if ( present(areps_type) ) then
         lareps = areps_type  
      else
         lareps = .false.  
      end if      
      if ( present(tracer_type) ) then
         ltracer = tracer_type
      else
         ltracer = .false.
      end if
      if ( present(tn_type) ) then
         ltn = tn_type  
      else
         ltn = 0
      end if
      if ( present(output_scale) ) then
         scale = output_scale
      else
         scale = 0.0
      end if

      if ( hist_debug > 0 ) then
         print*, ' Adding ', name, valid_min, valid_max
      end if

!     For completeness, initialise everything here, even variables that
!     must be redefined before use.
      histinfo(totflds)%name         = name
      histinfo(totflds)%amip_name    = aname
      histinfo(totflds)%long_name    = long_name
      if ( present(std_name) ) then
         histinfo(totflds)%std_name    = std_name
      else
         histinfo(totflds)%std_name    = ""
      end if
      histinfo(totflds)%units        = units  
      histinfo(totflds)%valid_min    = valid_min
      histinfo(totflds)%valid_max    = valid_max
      histinfo(totflds)%nlevels      = nlevels   ! number of levels
      histinfo(totflds)%std          = lstd
      histinfo(totflds)%ran          = lran
      histinfo(totflds)%areps        = lareps
      histinfo(totflds)%tracer       = ltracer
      histinfo(totflds)%tn           = ltn
      histinfo(totflds)%output_scale = scale
      histinfo(totflds)%used         = .FALSE.   ! Value for now
      histinfo(totflds)%vid          = 0
      histinfo(totflds)%ptr          = 0         ! start index in history array
      histinfo(totflds)%count        = 0         ! count for averaging
      histinfo(totflds)%ave_type     = atype     ! Default if not over-ridden
      histinfo(totflds)%addoff       = 0.
      histinfo(totflds)%scalef       = 0.
      histinfo(totflds)%ncid         = 0         ! Value for now
      histinfo(totflds)%procid       = 0         ! Output process
      if (present(int_type)) then
         histinfo(totflds)%int_type  = int_type
      else
         histinfo(totflds)%int_type  = int_default
      end if
      if ( present(multilev) ) then
         histinfo(totflds)%multilev = multilev
      else
         histinfo(totflds)%multilev = .false.
      end if
      if ( present(water) ) then
         histinfo(totflds)%water = water
      else
         histinfo(totflds)%water = .false.
      end if
      if ( present(pop2d) ) then
         histinfo(totflds)%pop2d = pop2d
      else
         histinfo(totflds)%pop2d = .false.
      end if
      if ( present(pop3d) ) then
         histinfo(totflds)%pop3d = pop3d
      else
         histinfo(totflds)%pop3d = .false.
      end if
      if ( present(pop4d) ) then
         histinfo(totflds)%pop4d = pop4d
      else
         histinfo(totflds)%pop4d = .false.
      end if
      if ( present(soil) ) then
         histinfo(totflds)%soil = soil
      else
         histinfo(totflds)%soil = .false.
      end if
      if ( present(coord_height) ) then
         histinfo(totflds)%coord_height = coord_height
      else
         histinfo(totflds)%coord_height = -huge(1.)
      end if
      if ( present(coord_name) ) then
         histinfo(totflds)%coord_name = coord_name
      else
         histinfo(totflds)%coord_name = ""
      end if
      if ( present(coord_stdname) ) then
         histinfo(totflds)%coord_stdname = coord_stdname
      else
         histinfo(totflds)%coord_stdname = ""
      end if
      if ( present(coord_units) ) then
         histinfo(totflds)%coord_units = coord_units
      else
         histinfo(totflds)%coord_units = ""
      end if
      if ( present(coord_positive) ) then
         histinfo(totflds)%coord_positive = coord_positive
      else
         histinfo(totflds)%coord_positive = ""
      end if
      if ( present(cell_methods) ) then
         histinfo(totflds)%cell_methods = cell_methods
      else
         histinfo(totflds)%cell_methods = ""
      end if
      if ( present(daily) ) then
         histinfo(totflds)%daily = daily
      else
         histinfo(totflds)%daily = .false.
      end if
      if ( present(sixhr) ) then
         histinfo(totflds)%sixhr = sixhr
      else
         histinfo(totflds)%sixhr = .false.
      end if      
      if ( present(instant) ) then
         histinfo(totflds)%instant = instant
      else
         histinfo(totflds)%instant = .true.
      end if
      if ( present(all_positive) ) then
         histinfo(totflds)%all_positive = all_positive
      else
         histinfo(totflds)%all_positive = .false.
      end if
      if ( present(fill) ) then
         histinfo(totflds)%fill = fill
      else
         histinfo(totflds)%fill = .false.
      end if

   end subroutine addfld
   
!-------------------------------------------------------------------
   subroutine openhist ( nx, ny, nl, sig, ol, cptch, cchrt, gosig, suffix, hlon, hlat, basetime, &
                         year, nxout, nyout, source, histfilename,                               &
                         pressure, height, theta, pvort, depth, extra_atts, hybrid_levels, anf,  &
                         bnf, p0, calendar, nsoil, zsoil )
!
!     Create netCDF history files using information in histinfo array.
!
      use interp_m, only : int_none
      use logging_m
      use mpidata_m

      integer, intent(in) :: nx, ny, nl, ol, cptch, cchrt
      real, intent(in), dimension(:) :: sig
      real, intent(in), dimension(:) :: gosig
      character(len=*), intent(in)   :: suffix  ! Filename suffix
!     Longitudes and latitudes of the output history
      real, intent(in), dimension(:) :: hlon
      real, intent(in), dimension(:) :: hlat
      character(len=*), intent(in) :: basetime
      integer, intent(in), optional  :: year
!     Dimensions of array in output file if different from model resolution
      integer, intent(in), optional  :: nxout, nyout
      character(len=*), intent(in), optional :: source
      character(len=*), intent(in), optional :: histfilename
!     Output uses sigma or pressure as vertical coordinate
      logical, intent(in), optional :: pressure
      logical, intent(in), optional :: height
      logical, intent(in), optional :: depth
      logical, intent(in), optional :: theta
      logical, intent(in), optional :: pvort
      type(hist_att), dimension(:), optional :: extra_atts
      logical, intent(in), optional :: hybrid_levels
      real, dimension(:), intent(in), optional :: anf, bnf
      real, intent(in), optional :: p0
      character(len=*), intent(in), optional :: calendar
      integer, intent(in), optional :: nsoil ! Number of soil levels
      ! Soil depths
      real, dimension(:), intent(in), optional :: zsoil
      real, parameter :: radtodeg=57.29577951 
      real, dimension(:,:), allocatable :: lon_bnds, lat_bnds, zsoil_bnds

      integer :: ierr, hsize, ifld, old_mode
      character(len=200) :: filename
      character(len=300) :: singlefilename
      type(dimarray) :: dims, dimvars
      integer :: ncid, vid, ivar, istart, iend, ilev
      character(len=80) :: longname, units
      character(len=MAX_NAMELEN) :: vname
      logical :: used, use_plevs, use_hyblevs, use_meters, use_depth
      logical :: use_theta, use_pvort
      !logical :: multilev
      logical :: multilev_fld
      !integer :: nlevels
      integer :: nlevels_fld
      real, dimension(totflds) :: coord_height
      character(len=80), dimension(totflds) :: coord_name
      character(len=80), dimension(totflds) :: coord_stdname
      character(len=80), dimension(totflds) :: coord_units
      character(len=80), dimension(totflds) :: coord_positive
      integer :: kc, ncoords, k, pkl
      integer :: ol_fld, nsoil_fld, cptch_fld, cchrt_fld
      logical :: osig_found
      integer :: i, j
      real, dimension(:), allocatable :: cabledata
      real(kind=8) :: dx, dy
      real(kind=8) :: hlonr8, hlatr8
      integer :: ave_type, nlev, rrank
      integer :: cnt, maxcnt, gap, slab

      call START_LOG(openhist_begin)
      
!     Set the module variables
      if ( present (nxout) ) then
         nxhis = nxout
         if ( nxout /= nx ) require_interp = .true.
      else
         nxhis = nx
      end if
      if ( present (nyout) ) then
         nyhis = nyout
         if ( nyout /= ny ) require_interp = .true.
      else
         nyhis = ny
      end if

      if ( size(hlat) /= nyhis ) then
         print*, " Error, mismatch in number of latitudes ", size(hlat), nyhis
         stop
      end if
      if ( size(hlon) /= nxhis ) then
         print*, " Error, mismatch in number of longitudes ", size(hlon), nxhis
         stop
      end if

      if ( hist_debug > 2 ) then
         do ivar = 1,totflds
            print*, ivar, histinfo(ivar)%name, histinfo(ivar)%nlevels
         end do
      end if

      call sortlist

      if ( hist_debug > 2 ) then
         print*, " Names after sorting "
         do ivar = 1,totflds
            print*, ivar, histinfo(ivar)%name, histinfo(ivar)%nlevels
         end do
      end if

!     Save this as a module variable
      filesuffix = suffix

      do ivar = 1,totflds
!        Go until the end of the variable list
         if ( len_trim(hnames(ivar)) == 0 ) then
            exit
         end if
         if ( hist_debug > 2 ) then
            print*, 'Openhist',  ivar, adjustl(hnames(ivar))
         end if 
            
         if ( hnames(ivar) == "all" ) then
            histinfo(1:totflds)%used = .true.
            exit ! All defined so no point checking further
                 ! Except when xnames is used.
         else if ( hnames(ivar) == "std" ) then
!           Use or here so as not to overwrite variables that have been
!           set explicitly.
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%std
         else if ( hnames(ivar) == "areps" ) then
!           Include RAN requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%areps
         else if ( hnames(ivar) == "ran" ) then
!           Include RAN requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%ran
         else if ( hnames(ivar) == "tracer" ) then
!           Include Tracer requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%tracer
         else if ( hnames(ivar) == "pop2d" ) then
!           Include POP 3d requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%pop2d
         else if ( hnames(ivar) == "pop3d" ) then
!           Include POP 3d requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%pop3d
         else if ( hnames(ivar) == "pop4d" ) then
!           Include POP 4d requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%pop4d
         else if ( hnames(ivar) == "pop" ) then
!           Include POP requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%pop2d .or.   &
             histinfo(1:totflds)%pop3d .or. histinfo(1:totflds)%pop4d
         else if ( hnames(ivar) == "t1_pop" ) then
!           Include POP requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. btest(histinfo(1:totflds)%tn,0)
         else if ( hnames(ivar) == "t2_pop" ) then
!           Include POP requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. btest(histinfo(1:totflds)%tn,1)
         else if ( hnames(ivar) == "t3_pop" ) then
!           Include POP requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. btest(histinfo(1:totflds)%tn,2)
         else if ( hnames(ivar) == "t4_pop" ) then
!           Include POP requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. btest(histinfo(1:totflds)%tn,3)
         else if ( hnames(ivar) == "t5_pop" ) then
!           Include POP requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. btest(histinfo(1:totflds)%tn,4)
         else
!           Find the name in histinfo to set the used flag.
            ifld = bindex_hname ( hnames(ivar), &
                                  inames(:,1:totflds), totflds )
            if ( ifld == 0 ) then
               print*, "Error - history variable ", hnames(ivar),  &
                       " is not known. "
               stop
            end if
            histinfo(ifld)%used = .true.
            if ( hist_debug > 4 ) then
               print*,  " Setting used for ", ifld, ivar, &
                    hnames(ivar), histinfo(ifld)%name
            end if
         end if

      end do
         
!     Apply the exclusion list
      do ivar = 1,totflds

!        Go until the end of the variable list
         if ( len_trim(xnames(ivar)) == 0 ) then
            exit
         end if
         if ( hist_debug > 2 ) then
            print*, 'Openhist: exclude',  ivar, adjustl(xnames(ivar))
         end if
            
!        Find the name in histinfo to set the used flag.
         ifld = bindex_hname ( xnames(ivar), inames(:,1:totflds), totflds )
         if ( ifld == 0 ) then
            print*, "Error - excluded history variable ", xnames(ivar)," is not known. "
            stop
         end if
         histinfo(ifld)%used = .false.
         if ( hist_debug > 4 ) then
            print*,  " Overriding used for ", ifld, ivar, xnames(ivar), histinfo(ifld)%name
         end if

      end do

!     Exit on the first file with no variables set.
      if ( .not. any(histinfo(1:totflds)%used) ) then
         return
      end if

!     Check the averaging type
!     For variables for which it isn't set use the file default
      where ( histinfo(1:totflds)%ave_type == 0 )
         histinfo(1:totflds)%ave_type = ihtype
      endwhere

!     The rest of the history files are much simpler with multilevel variables

      if ( present(histfilename) ) then
         filename = histfilename
      else
         write(filename,"(a,i1,a,a)" ) "hist", trim(suffix), ".nc"
      end if

      use_plevs = .false.
      if ( present(pressure) ) then
         use_plevs = pressure
      end if
      use_meters = .false.
      if ( present(height) ) then
         use_meters = height
      end if
      use_hyblevs = .false.
      if ( present(hybrid_levels) ) then
         use_hyblevs = hybrid_levels
      end if
      use_depth = .false.
      if ( present(depth) ) then
         use_depth = depth 
      end if
      use_theta = .false.
      if ( present(theta) ) then
         use_theta = theta
      end if
      use_pvort = .false.
      if ( present(pvort) ) then
         use_pvort = pvort
      end if
      if ( ol > 0 ) then
         osig_found = all(gosig<=1.)
      else 
         osig_found = .true.
      end if    
         
      if ( single_output ) then
         if ( myid == 0 ) then 
            ! check what dimensions need to be written 
            multilev_fld = .false.
            nlevels_fld = 1
            ol_fld = 0
            nsoil_fld = 0
            cptch_fld = 0
            cchrt_fld = 0
            do ifld = 1,totflds
               ! Get a list of the coordinate heights if any
               call coord_list(ifld,ncoords,coord_height,coord_name, &
                               coord_stdname,coord_units,coord_positive)
               ! check dimensons and other file properties
               if ( histinfo(ifld)%used ) then 
                  if ( histinfo(ifld)%nlevels > 1 .or. histinfo(ifld)%multilev ) then
                     multilev_fld = .true.
                  end if   
                  if ( .not.histinfo(ifld)%water .and. .not.histinfo(ifld)%soil .and.  &
                       .not.histinfo(ifld)%pop2d .and. .not.histinfo(ifld)%pop3d .and. &
                       .not.histinfo(ifld)%pop4d ) then
                     nlevels_fld = size(sig) 
                  end if    
                  if ( histinfo(ifld)%water ) then
                     ol_fld = ol
                  end if                
                  if ( histinfo(ifld)%soil ) then
                     nsoil_fld = nsoil
                  end if    
                  if ( histinfo(ifld)%pop2d .or. histinfo(ifld)%pop3d .or. &
                       histinfo(ifld)%pop4d ) then
                     cptch_fld = cptch
                     cchrt_fld = cchrt
                  end if
               end if ! histinfo(ifld)%used
            end do    ! ifld = 1,totflds
            ! patch for AREPS
            if ( areps_compliant ) then
               ! force soil depth output for AREPS 
               nsoil_fld = nsoil 
            end if
            allocate( histfile(1), histdimvars(1) )
            call create_ncfile ( filename, nxhis, nyhis, nlevels_fld, ol_fld, cptch_fld, cchrt_fld,  &
                 multilev_fld, use_plevs, use_meters, use_theta, use_pvort, use_depth, use_hyblevs,  &
                 basetime, coord_name(1:ncoords), coord_stdname(1:ncoords), coord_units(1:ncoords),  &
                 coord_positive(1:ncoords), ncid, dims, dimvars, source, extra_atts,                 &
                 calendar, nsoil_fld, osig_found, instant=.false. )
            histfile(1)%id = ncid
            histfile(1)%day = .false.
            histfile(1)%inst = .true.
            histfile(1)%sixhr = .false.
            histfile(1)%fix = .false.
            histdimvars(1) = dimvars
            histfile(1)%nlevels = nlevels_fld
            histfile(1)%ocean = ol_fld > 0
            histfile(1)%soil = nsoil_fld > 0
            histfile(1)%cablepatch = cptch_fld > 0
            histfile(1)%cablecohort = cchrt_fld > 0
            histfile(1)%ncoords = ncoords
            allocate( histfile(1)%coord_height(1:ncoords) )
            histfile(1)%coord_height(1:ncoords) = coord_height(1:ncoords)
            allocate( histfile(1)%coord_name(1:ncoords) )
            histfile(1)%coord_name(1:ncoords) = coord_name(1:ncoords)
            do ifld = 1,totflds
               histinfo(ifld)%ncid = ncid
               histinfo(ifld)%procid = 0
               histdims(ifld) = dims
            end do   
         else  
            allocate( histfile(0), histdimvars(0) ) 
         end if  
      else
         ! count number of output variables to be interpolated
         maxcnt = 0
         do ifld = 1,totflds
           if ( .not. histinfo(ifld)%used ) then
               cycle
            end if  
            maxcnt = maxcnt + histinfo(ifld)%nlevels
         end do
         ! calculate slab and gap to distribute interpolation across processes
         slab = ceiling(real(maxcnt,8)/real(nproc,8))
         gap = max( nproc/maxcnt, 1 )
         ! assign process id for each output file
         cnt = 0
         do ifld = 1,totflds
            if ( .not. histinfo(ifld)%used ) then
               cycle
            end if  
            histinfo(ifld)%procid = (cnt*gap)/slab
            cnt = cnt + histinfo(ifld)%nlevels
         end do   
         ! create output files
         i = 0
         do ifld = 1,totflds
            if ( histinfo(ifld)%used .and. histinfo(ifld)%procid==myid ) then
               i = i + 1 
            end if    
         end do
         allocate( histfile(i), histdimvars(i) )
         i = 0
         do ifld = 1,totflds
            ! Get a list of the coordinate heights if any
            call coord_list(ifld,ncoords,coord_height,coord_name, &
                            coord_stdname,coord_units,coord_positive)
            ! Determine sizes of dimensions for this file 
            multilev_fld = histinfo(ifld)%nlevels > 1 .or. &
                           histinfo(ifld)%multilev
            if ( histinfo(ifld)%water .or. histinfo(ifld)%soil .or.  &
                 histinfo(ifld)%pop2d .or. histinfo(ifld)%pop3d .or. &
                 histinfo(ifld)%pop4d ) then
               nlevels_fld = 1
            else
               nlevels_fld = histinfo(ifld)%nlevels
            end if    
            if ( histinfo(ifld)%water ) then
               ol_fld = ol
            else
               ol_fld = 0
            end if
            if ( histinfo(ifld)%soil ) then
               nsoil_fld = nsoil 
            else
               nsoil_fld = 0
            end if
            if ( histinfo(ifld)%pop2d .or. histinfo(ifld)%pop3d .or. &
                 histinfo(ifld)%pop4d ) then
               cptch_fld = cptch
               cchrt_fld = cchrt
            else
               cptch_fld = 0
               cchrt_fld = 0
            end if
            if ( histinfo(ifld)%used .and. histinfo(ifld)%procid==myid ) then 
               i = i + 1
               singlefilename = calcfilename(histinfo(ifld)%name,filename)
               if ( histinfo(ifld)%ave_type == hist_fixed ) then
                  call create_ncfile ( singlefilename, nxhis, nyhis, nlevels_fld, ol_fld, cptch_fld, cchrt_fld, &
                       multilev_fld, use_plevs, use_meters, use_theta, use_pvort, use_depth, use_hyblevs,       &
                       "none", coord_name(1:ncoords), coord_stdname(1:ncoords), coord_units(1:ncoords),         &
                       coord_positive(1:ncoords), ncid, dims, dimvars, source, extra_atts, calendar,            &
                       nsoil_fld, osig_found, instant=histinfo(ifld)%instant )
               else    
                  call create_ncfile ( singlefilename, nxhis, nyhis, nlevels_fld, ol_fld, cptch_fld, cchrt_fld, &
                       multilev_fld, use_plevs, use_meters, use_theta, use_pvort, use_depth, use_hyblevs,       &
                       basetime, coord_name(1:ncoords), coord_stdname(1:ncoords), coord_units(1:ncoords),       &
                       coord_positive(1:ncoords), ncid, dims, dimvars, source, extra_atts, calendar,            &
                       nsoil_fld, osig_found, instant=histinfo(ifld)%instant )
               end if
               histfile(i)%id = ncid
               histfile(i)%day = histinfo(ifld)%daily
               histfile(i)%inst = histinfo(ifld)%instant
               histfile(i)%sixhr = histinfo(ifld)%sixhr
               histfile(i)%fix = histinfo(ifld)%ave_type == hist_fixed
               histdimvars(i) = dimvars
               histfile(i)%nlevels = nlevels_fld
               histfile(i)%ocean = ol_fld > 0
               histfile(i)%soil = nsoil_fld > 0
               histfile(i)%cablepatch = cptch_fld > 0
               histfile(i)%cablecohort = cchrt_fld > 0
               histfile(i)%ncoords = ncoords
               allocate( histfile(i)%coord_height(1:ncoords) )
               histfile(i)%coord_height(1:ncoords) = coord_height(1:ncoords)
               allocate( histfile(i)%coord_name(1:ncoords) )
               histfile(i)%coord_name(1:ncoords) = coord_name(1:ncoords)
               histinfo(ifld)%ncid = ncid
               histdims(ifld) = dims
            end if   
         end do
      end if  
      
      do ifld = 1,totflds
         if ( histinfo(ifld)%used .and. myid==histinfo(ifld)%procid ) then
            ncid = histinfo(ifld)%ncid
            dims = histdims(ifld)
            call create_ncvar(histinfo(ifld), ncid, dims)
         end if
      end do

      do i = 1,size(histfile) 
         ncid = histfile(i)%id
         dimvars = histdimvars(i)
         !        Leave define mode 
         ierr = nf90_enddef ( ncid )
         call check_ncerr(ierr, "Error from enddef")
         !        Turn off the data filling to save time.
         ierr = nf90_set_fill ( ncid, NF90_NOFILL, old_mode)
         call check_ncerr(ierr, "Error from set_fill")
         ierr = nf90_put_var ( ncid, dimvars%y, hlat )
         call check_ncerr(ierr,"Error writing latitudes")
         ierr = nf90_put_var ( ncid, dimvars%x, hlon )
         call check_ncerr(ierr,"Error writing longitudes")
      end do
            
      if ( cf_compliant .or. cordex_compliant ) then
         ! Calculate bounds assuming a regular lat-lon grid
         ! Perhaps have optional arguments for the other cases?
         allocate ( lat_bnds(2,size(hlat)), lon_bnds(2,size(hlon)) )
         ! Check if regular grid
         if ( int_default /= int_none ) then
            dx = real(hlon(2),8) - real(hlon(1),8)
            do i = 1,size(hlon)
              hlonr8 = real(hlon(i),8) - 0.5_8*dx
              lon_bnds(1,i) = real(real(nint(hlonr8*1.e5_8),8)*1.e-5_8)
              hlonr8 = real(hlon(i),8) + 0.5_8*dx
              lon_bnds(2,i) = real(real(nint(hlonr8*1.e5_8),8)*1.e-5_8)
            end do  
            ! MJT fix for rounding errors
            do i = 2,size(hlon)
               lon_bnds(1,i) = lon_bnds(2,i-1)
            end do
            do i = 1,size(histfile)
               ncid = histfile(i)%id 
               dimvars = histdimvars(i)
               ierr = nf90_put_var ( ncid, dimvars%x_b, lon_bnds )
               call check_ncerr(ierr,"Error writing longitude bounds")
            end do
         end if
         if ( int_default /= int_none ) then
            dy = real(hlat(2),8) - real(hlat(1),8)
            do j = 1,size(hlat)
              hlatr8 = real(hlat(j),8) - 0.5_8*dy
              lat_bnds(1,j) = real(real(nint(hlatr8*1.e5_8),8)*1.e-5_8)
              hlatr8 = real(hlat(j),8) + 0.5_8*dy
              lat_bnds(2,j) = real(real(nint(hlatr8*1.e5_8),8)*1.e-5_8)
            end do  
            where ( lat_bnds < -90. ) 
               lat_bnds = -90.
            end where
            where ( lat_bnds > 90. ) 
               lat_bnds = 90.
            end where
            ! MJT fix for rounding errors
            do i = 2,size(hlat)
               lat_bnds(1,i) = lat_bnds(2,i-1)
            end do    
            do i = 1,size(histfile)
               ncid = histfile(i)%id 
               dimvars = histdimvars(i)
               ierr = nf90_put_var ( ncid, dimvars%y_b, lat_bnds )
               call check_ncerr(ierr,"Error writing latitude bounds")
            end do   
         end if
      end if
      do i = 1,size(histfile)
         if ( histfile(i)%nlevels > 1 ) then
            ncid = histfile(i)%id 
            dimvars = histdimvars(i)
            ierr = nf90_put_var ( ncid, dimvars%z, sig )
            call check_ncerr(ierr,"Error writing levels")
         end if
      end do   
      if ( use_hyblevs ) then
         if ( .not. present(anf) ) then
            print*, "Error, missing anf argument"
            stop
         end if
         do i = 1,size(histfile)
            if ( histfile(i)%nlevels > 1 ) then
               ncid = histfile(i)%id 
               ierr = nf90_inq_varid(ncid, "anf", vid)
               call check_ncerr(ierr,"Error getting vid for anf")
               ierr = nf90_put_var(ncid, vid, anf)
               call check_ncerr(ierr,"Error writing anf")
               if ( .not. present(bnf) ) then
                  print*, "Error, missing bnf argument"
                  stop
               end if
               ierr = nf90_inq_varid(ncid, "bnf", vid)
               call check_ncerr(ierr,"Error getting vid for bnf")
               ierr = nf90_put_var(ncid, vid, bnf)
               call check_ncerr(ierr,"Error writing bnf")
               if ( .not. present(p0) ) then
                  print*, "Error, missing p0 argument"
                  stop
               end if
               ierr = nf90_inq_varid(ncid, "P0", vid)
               call check_ncerr(ierr,"Error getting vid for p0")
               ierr = nf90_put_var(ncid, vid, p0)
               call check_ncerr(ierr,"Error writing p0")
            end if   
         end do   
      end if
      if ( ol > 0 ) then
         do i = 1,size(histfile)
            if ( histfile(i)%ocean ) then 
               ncid = histfile(i)%id 
               dimvars = histdimvars(i)
               ierr = nf90_put_var ( ncid, dimvars%oz, gosig )
               call check_ncerr(ierr,"Error writing olev")
            end if   
         end do   
      end if
      if ( cptch > 0 ) then
         allocate( cabledata(cptch) )
         do i = 1,cptch
            cabledata(i) = real(i)
         end do
         do i = 1,size(histfile)
            if ( histfile(i)%cablepatch ) then 
               ncid = histfile(i)%id 
               dimvars = histdimvars(i)
               ierr = nf90_put_var ( ncid, dimvars%cptch, cabledata )
               call check_ncerr(ierr,"Error writing cable_patch")
            end if   
         end do   
         deallocate( cabledata )
      end if
      if ( cchrt > 0 ) then
         allocate( cabledata(cchrt) )
         do i = 1,cchrt
            cabledata(i) = real(i)
         end do
         do i = 1,size(histfile)
            if ( histfile(i)%cablecohort ) then 
               ncid = histfile(i)%id 
               dimvars = histdimvars(i)
               ierr = nf90_put_var ( ncid, dimvars%cchrt, cabledata )
               call check_ncerr(ierr,"Error writing cable_cohort")
            end if   
         end do   
         deallocate( cabledata )
      end if
      if ( present(zsoil) .and. present(nsoil) ) then
         if ( nsoil>0 ) then
            do i = 1,size(histfile)
               if ( histfile(i)%soil ) then 
                  ncid = histfile(i)%id 
                  dimvars = histdimvars(i)
                  ierr = nf90_put_var ( ncid, dimvars%zsoil, zsoil )
                  call check_ncerr(ierr,"Error writing depths")
               end if   
            end do  
            if ( cf_compliant .or. cordex_compliant ) then
               ! Soil bounds
               allocate(zsoil_bnds(2, nsoil))
               zsoil_bnds(1,1) = 0.
               zsoil_bnds(2,1) = 2.*zsoil(1)
               do k = 2,nsoil
                  ! Levels are middle of layers
                  zsoil_bnds(2,k) = zsoil_bnds(2,k-1) + 2.*(zsoil(k)-zsoil_bnds(2,k-1))
                  zsoil_bnds(1,k) = zsoil_bnds(2,k-1)
               end do
               do i = 1,size(histfile)
                  if ( histfile(i)%soil ) then 
                     ncid = histfile(i)%id
                     dimvars = histdimvars(i)
                     ierr = nf90_put_var ( ncid, dimvars%zsoil_b, zsoil_bnds )
                     call check_ncerr(ierr,"Error writing depths")
                  end if   
               end do   
               deallocate(zsoil_bnds)
            end if
         end if
      end if

      do i = 1,size(histfile)      
         do kc = 1,histfile(i)%ncoords
            ncid = histfile(i)%id 
            ierr = nf90_inq_varid(ncid, histfile(i)%coord_name(kc), vid)
            call check_ncerr(ierr,"Error getting vid for height coord")
            ierr = nf90_put_var ( ncid, vid, histfile(i)%coord_height(kc) )
            call check_ncerr(ierr,"Error writing coordinate height")
         end do   
      end do
      
!     Allocate the array to hold all the history data
!     Calculate the size by summing the number of fields of each variable

      hsize = 0
      do ifld = 1, totflds
         if ( histinfo(ifld)%used ) then
            histinfo(ifld)%ptr = hsize + 1
            hsize = hsize + histinfo(ifld)%nlevels 
            if ( hist_debug > 0 ) then
               print*, ' Fld ', ifld, 'Ptr', histinfo(ifld)%ptr
            end if
         end if
      end do

      if ( hist_debug > 0 ) then
         print*, ' History array size ', hsize
      end if

!     This is allocated as a single array. It might be cleaner to use
!     individual arrays? A single array has the advantage that it's easier
!     to change the order of the dimensions if that's required for efficiency.

!     In multi-month runs of the climate model it may already be allocated
      if (allocated(histarray)) then
         if ( .not. all ( shape(histarray) == (/pil,pjl*pnpan*lproc,hsize/) ) ) then
            print*, "Error: size of histarray has changed"
            print*, size(histarray), pil, pjl*pnpan*lproc, hsize
            stop
         end if
      else
         allocate ( histarray(pil,pjl*pnpan*lproc,hsize), stat=ierr )
         if ( ierr /= 0 ) then
            print*, "Error: Failed to allocate history array of size", &
                 pil*pjl*pnpan*lproc*hsize, "words"
            stop
         end if
         nx_g = nx
         ny_g = ny
      end if

!     Initialise the history appropriately
      histset = 0
      histset_daily = 0
      histset_6hr = 0
      avetime = 0.0
      ! Initialisation so min/max work
      avetime_bnds(:) = (/huge(1.), -huge(1.)/)
      timecount = 0
      do ifld = 1, totflds
         if ( histinfo(ifld)%used ) then
            istart = histinfo(ifld)%ptr
            iend = istart + histinfo(ifld)%nlevels - 1
            histarray(:,:,istart:iend) = nf90_fill_float
         end if
      end do
      
      call END_LOG(openhist_end)

   end subroutine openhist

   subroutine sortlist

!     Simple insertion sort of the list of names.
!     Assumes all lower case. This should be enforced or checked somewhere

      integer :: i, ipos, j
      type(hinfo) :: temp
      integer (bigint), dimension(MAX_KEYINDEX) :: itemp
      
      if ( MAX_KEYINDEX < MAX_NAMELEN/MAX_KEYLEN ) then
         print*, "Error: MAX_KEYINDEX is too small"
         stop
      end if
          
      do i = 1,totflds
         inames(:,i) = hashkey ( histinfo(i)%name )
      end do

      do i = 1,totflds

!        Find the first element in the rest of the list
         itemp(:) = inames(:,i)
         ipos = i
         do j = i+1,totflds
            if ( hash_lt(inames(:,j),itemp) ) then 
               itemp(:) = inames(:,j)
               ipos = j
            end if
         end do

!        Move the smallest value to position i
         inames(:,ipos) = inames(:,i)
         inames(:,i) = itemp(:)
!        Swap histinfo elements so they keep the same order
         temp = histinfo(ipos)
         histinfo(ipos) = histinfo(i)
         histinfo(i) = temp

      end do

      if ( hist_debug > 1 ) then
         print*, "Sorted NAMES "
         do i = 1,totflds
            print*, histinfo(i)%name, inames(:,i)
         end do
      end if

    end subroutine sortlist
    
    function hash_lt(a,b) result(ans)
       integer(bigint), dimension(:), intent(in) :: a,b
       logical :: ans
       integer :: k
       
       ans = .false.
       do k = 1,MAX_KEYINDEX
          if ( a(k) < b(k) ) then
             ans = .true.
             exit
          else if ( a(k) > b(k) ) then
             ans = .false.
             exit
          end if
       end do
    
    end function hash_lt
    
    function hash_gt(a,b) result(ans)
       integer(bigint), dimension(:), intent(in) :: a,b
       logical :: ans
       integer :: k
       
       ans = .false.
       do k = 1,MAX_KEYINDEX
          if ( a(k) > b(k) ) then
             ans = .true.
             exit
          else if ( a(k) < b(k) ) then
             ans = .false.
             exit
          end if
       end do
    
    end function hash_gt

    subroutine coord_list(ifld,ncoords,coord_height,coord_name, &
                          coord_stdname,coord_units,coord_positive)
       integer, intent(in) :: ifld
       integer, intent(out) :: ncoords
       integer :: kc
       real, dimension(:), intent(inout) :: coord_height
       character(len=*), dimension(:), intent(inout) :: coord_name
       character(len=*), dimension(:), intent(inout) :: coord_stdname
       character(len=*), dimension(:), intent(inout) :: coord_units
       character(len=*), dimension(:), intent(inout) :: coord_positive
    
       ncoords = 0
       if ( histinfo(ifld)%coord_height > -huge(1.) ) then
          ! Check if it's already in list
          do kc = 1,ncoords
             if ( histinfo(ifld)%coord_name == coord_name(kc) ) then 
                exit
             end if
          end do
          if ( kc > ncoords ) then
             ! Value not found
             ncoords = kc
             coord_height(ncoords) = histinfo(ifld)%coord_height
             coord_name(ncoords) = histinfo(ifld)%coord_name
             coord_stdname(ncoords) = histinfo(ifld)%coord_stdname
             coord_units(ncoords) = histinfo(ifld)%coord_units
             coord_positive(ncoords) = histinfo(ifld)%coord_positive
          end if
       end if
    
    end subroutine coord_list   
    
!-------------------------------------------------------------------
   subroutine create_ncvar(vinfo, ncid, dims)

      use mpidata_m
      type(hinfo), intent(inout) :: vinfo
      integer, intent(in) :: ncid
      type(dimarray), intent(in) :: dims

      character(len=MAX_NAMELEN) :: local_name, new_name
      character(len=80) :: cell_methods
      !character(len=80) :: coord_name
      integer :: ierr, vtype, vid, zdim, wdim
      integer :: ndims
      integer, dimension(5) :: chunks
      
      integer(kind=2), parameter :: fill_short = NF90_FILL_SHORT

      local_name = vinfo%name  

      select case ( hbytes )
      case ( 2 )
         vtype = NF90_INT2
      case ( 4 ) 
         vtype = NF90_FLOAT
      case ( 8 )
         vtype = NF90_DOUBLE
      case default
         print*, " Error, unsupported value for hbytes ", hbytes
         stop
      end select
         
      if ( hist_debug > 4 ) then
         print*, "Creating variable ", local_name, vinfo%ave_type, vinfo%nlevels, vinfo%soil
         print*, "DIMID", dims%x, dims%y, dims%z, dims%t, dims%zsoil
      end if
      if ( vinfo%nlevels > 1 .or. vinfo%multilev ) then
         if (vinfo%soil) then
            zdim = dims%zsoil
         else if ( vinfo%water ) then
            zdim = dims%oz
         else if ( vinfo%pop3d ) then
            zdim = dims%cptch
         else if ( vinfo%pop4d ) then
            zdim = dims%cptch
            wdim = dims%cchrt
         else
            zdim = dims%z
         end if
         if ( vinfo%ave_type == hist_fixed ) then
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim /), vid )
            ndims = 3
         else if ( vinfo%pop4d ) then
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim, wdim, dims%t /), vid )
            ndims = 5
         else
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim, dims%t /), vid )
            ndims = 4
         end if
      else
         if ( vinfo%ave_type == hist_fixed ) then
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y /), vid )
            ndims = 2
         else
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, dims%t /), vid )
            ndims = 3
         end if
      end if
      if ( ierr /= 0 ) then
         print*, "Error creating variable ", vinfo
      end if
      call check_ncerr(ierr,"Error creating variable "// local_name)
      vinfo%vid = vid
      
      ! Deflate and chunking
#ifndef usenc3      
      if ( cordex_compliant ) then
         ierr = nf90_def_var_deflate( ncid, vid, 1, 1, 1 ) ! shuffle=1, deflate=1, deflate_level=1
         call check_ncerr(ierr)
      end if
      chunks(:) = 1
      if ( chunk_grid>0 ) then
         chunks(1) = min( nxhis, chunk_grid ) 
         chunks(2) = min( nyhis, chunk_grid )
         ierr = nf90_def_var_chunking( ncid, vid, NF90_CHUNKED, chunks(1:ndims) )
         call check_ncerr(ierr)
      end if  
#endif      

      if ( len_trim(vinfo%long_name) /= 0 ) then
         ierr = nf90_put_att ( ncid, vid, "long_name", vinfo%long_name )
         call check_ncerr(ierr)
      end if
      if ( len_trim(vinfo%std_name) /= 0 ) then
         ierr = nf90_put_att ( ncid, vid, "standard_name", vinfo%std_name )
         call check_ncerr(ierr)
      end if
      if ( len_trim(vinfo%units) /= 0 ) then
         ierr = nf90_put_att ( ncid, vid, "units", vinfo%units )
         call check_ncerr(ierr)
      end if

      if ( hbytes == 2 ) then
         ! 32500/3200 = 1.015625, so with this scaling, end values should
         ! be exactly representable.
         vinfo%scalef = 1.015625 * (vinfo%valid_max - vinfo%valid_min) / real(vmax - vmin)
         vinfo%addoff = 0.5*(vinfo%valid_min+vinfo%valid_max)
         ierr  = nf90_put_att ( ncid, vid, "add_offset", vinfo%addoff )
         call check_ncerr(ierr)
         ierr  = nf90_put_att ( ncid, vid, "scale_factor", vinfo%scalef)
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "valid_min", vmin )
         call check_ncerr(ierr,"Error setting valid min attribute")
         ierr = nf90_put_att ( ncid, vid, "valid_max", vmax )
         call check_ncerr(ierr,"Error setting valid max attribute")
      end if
      
      ! cell_methods from the history averaging process

      ! Possibly time:point should just be left out

      select case (vinfo%ave_type)
      !case (hist_ave, hist_oave)
      case (hist_ave )    
         cell_methods = "time: mean"
      case (hist_max)
         cell_methods = "time: maximum"
      case (hist_min)
         cell_methods = "time: minimum"
      case (hist_inst)
         cell_methods = ""
         ! cell_methods = "time: point"
      case (hist_fixed)
         cell_methods = ""
      case default
         print*, " History internal error: average type set incorrectly"
         stop
      end select

      if ( len_trim(vinfo%cell_methods) > 0 ) then
         ! Doesn't make sense to compose time: mean methods
         if ( .not. (cell_methods == "time: mean" .and. vinfo%cell_methods == "time: mean") ) then
            cell_methods = trim(vinfo%cell_methods) // " " // trim(cell_methods)
         end if
      end if

      if ( cordex_compliant ) then
         if ( len_trim(cell_methods) > 0 ) then
            ierr = nf90_put_att ( ncid, vid, "cell_methods", cell_methods )
            call check_ncerr(ierr,"Error with cell_methods attribute")
         else
            ! Default to time point if otherwise not defined 
            ierr = nf90_put_att ( ncid, vid, "cell_methods", "time: point" )
            call check_ncerr(ierr,"Error with cell_methods attribute")
         end if
      end if   
         
      if ( vtype == NF90_INT2 ) then
         ! Ugly work around to ensure attributes have the correct type on SX6
         if ( vinfo%ave_type /= hist_fixed ) then
            ierr = nf90_put_att ( ncid, vid, "_FillValue", fill_short )
            call check_ncerr(ierr,"Error with INT2 fill value attribute")
            ierr = nf90_put_att ( ncid, vid, "missing_value", fill_short )
            call check_ncerr(ierr,"Error with missing value attribute")
         end if
      else
         ierr = nf90_put_att ( ncid, vid, "_FillValue", missing_value_cordex )
         call check_ncerr(ierr,"Error with FLOAT/DOUBLE fill value attribute")
         ierr = nf90_put_att ( ncid, vid, "missing_value", missing_value_cordex )
         call check_ncerr(ierr,"Error with missing value attribute")
      end if

      if ( vinfo%coord_height > -huge(1.) ) then
         ! This will cover 2m and 10m. CORDEX also has 1000hPa, 925hPa, 800hPa, etc
         ierr = nf90_put_att ( ncid, vid, "coordinates", vinfo%coord_name )
         call check_ncerr(ierr,"Error with coordinates attribute")
      end if

      if ( cordex_compliant ) then
         ierr = nf90_put_att( ncid, vid, "grid_mapping", "crs" )      
      end if   
      
   end subroutine create_ncvar
  
!---------------------------------------------------------------------------
   subroutine create_ncfile ( filename, nxhis, nyhis, nlev, ol, cptch, cchrt, multilev,         &
                 use_plevs, use_meters, use_theta, use_pvort, use_depth, use_hyblevs, basetime, &
                 coord_name, coord_stdname, coord_units, coord_positive, ncid, dims, dimvars,   &
                 source, extra_atts, calendar, nsoil, osig_found, instant )

      use mpidata_m
      use newmpar_m, only : il
      use parm_m, only : schmidt

      character(len=*), intent(in) :: filename
      integer, intent(in) :: nxhis, nyhis, nlev, ol, cptch, cchrt
      logical, intent(in) :: multilev, use_plevs, use_meters, use_theta, use_pvort, use_depth, use_hyblevs
      logical, intent(in) :: osig_found
      character(len=*), intent(in) :: basetime
      character(len=*), dimension(:), intent(in) :: coord_name
      character(len=*), dimension(:), intent(in) :: coord_stdname
      character(len=*), dimension(:), intent(in) :: coord_units
      character(len=*), dimension(:), intent(in) :: coord_positive
      integer, intent(out) :: ncid
      type(dimarray), intent(out) :: dims, dimvars
      character(len=*), intent(in), optional :: source
      type(hist_att), dimension(:), optional :: extra_atts
      character(len=*), intent(in), optional :: calendar
      integer, intent(in), optional :: nsoil ! Number of soil levels
      logical, intent(in), optional :: instant

      integer :: ierr
      integer :: i, k
      character(len=80) :: histstr
      character(len=20) :: tmpname
      character(len=160) :: griddes
      character(len=20) :: gridsize, gridschmidt
      integer :: vid, cmode
      logical :: do_time_bnds

      do_time_bnds = .true. 
      if ( present(instant) .and. ihtype==hist_inst ) then
         do_time_bnds = .not.instant
      else if ( ihtype == hist_fixed ) then   
         do_time_bnds = .false. 
      end if
      
      if ( hist_debug > 0 ) then
         print*, "Creating file ", filename
      end if
      if ( areps_compliant ) then
         ! AREPS must use classic netcdf3 
         ierr = nf90_create(filename, nf90_clobber, ncid) 
      else    
#ifdef usenc3
         ierr = nf90_create(filename, nf90_64bit_offset, ncid)
#else
         !cmode = nf90_netcdf4 .or. nf90_classic_model
         ierr = nf90_create(filename, nf90_netcdf4, ncid)
#endif
      end if
      call check_ncerr ( ierr, "Error in creating history file "//trim(filename) )
               
      ! Create dimensions, lon, lat and rec
      if ( areps_compliant ) then
         ierr = nf90_def_dim ( ncid, "longitude", nxhis, dims%x ) 
      else    
         ierr = nf90_def_dim ( ncid, "lon", nxhis, dims%x )
      end if   
      call check_ncerr(ierr,"Error creating lon dimension")
      if ( hist_debug > 5 ) print*, "Created lon dimension, id",  dims%x
      if ( areps_compliant ) then
         ierr = nf90_def_dim ( ncid, "latitude", nyhis, dims%y ) 
      else    
         ierr = nf90_def_dim ( ncid, "lat", nyhis, dims%y )
      end if   
      call check_ncerr(ierr,"Error creating lat dimension")
      if ( hist_debug > 5 ) print*, "Created lon dimension, id",  dims%y
      if ( cf_compliant .or. cordex_compliant ) then
         ierr = nf90_def_dim ( ncid, "bnds", 2, dims%b )
         call check_ncerr(ierr,"Error creating bnds dimension")
         if ( hist_debug > 5 ) print*, "Created bnds dimension, id",  dims%b
      end if

      ! Only create the lev dimension if one of the variables actually uses it
      if ( nlev > 1 ) then
         if ( areps_compliant ) then
            ierr = nf90_def_dim ( ncid, "level", nlev, dims%z ) 
         else
            if ( use_meters ) then
               ierr = nf90_def_dim ( ncid, "alt", nlev, dims%z )
            else
               ierr = nf90_def_dim ( ncid, "lev", nlev, dims%z )
            end if
         end if   
         call check_ncerr(ierr,"Error creating lev dimension")
         if ( hist_debug > 5 ) print*, "Created lev dimension, id",  dims%z
      end if
      if ( ol > 0 ) then
         ierr = nf90_def_dim ( ncid, "olev", ol, dims%oz )
         call check_ncerr(ierr,"Error creating olev dimension")
         if ( hist_debug > 5 ) print*, "Created olev dimension, id", dims%oz
      end if
      if ( present(nsoil) ) then
         if ( nsoil > 0 ) then
            ierr = nf90_def_dim ( ncid, "depth", nsoil, dims%zsoil )
            call check_ncerr(ierr,"Error creating soil depth dimension")
            if ( hist_debug > 5 ) print*, "Created soil dimension, id",  dims%zsoil
         end if
      end if
      if ( cptch > 0 ) then
         ierr = nf90_def_dim ( ncid, "cable_patch", cptch, dims%cptch )
         call check_ncerr(ierr,"Error creating cable_patch dimension")
         if ( hist_debug > 5 ) print*, "Created cable_patch dimension, id", dims%cptch
      end if
      if ( cchrt > 0 ) then
         ierr = nf90_def_dim ( ncid, "cable_cohort", cchrt, dims%cchrt )
         call check_ncerr(ierr,"Error creating cable_cohort dimension")
         if ( hist_debug > 5 ) print*, "Created cable_cohort dimension, id", dims%cchrt
      end if
      if ( basetime /= "none" ) then
         ierr = nf90_def_dim ( ncid, "time", NF90_UNLIMITED, dims%t )
         call check_ncerr(ierr,"Error creating time dimension")
         if ( hist_debug > 5 ) print*, "Created time dimension, id",  dims%t
      end if   
      
      if ( present(source) ) then
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "source", source )
         call check_ncerr(ierr)
      end if
      if ( cordex_compliant ) then
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "Conventions", "CF-1.11" )
         call check_ncerr(ierr)
      else
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "Conventions", "CF-1.7" )
         call check_ncerr(ierr)
      end if
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "title", "CCAM simulation data" )
      call check_ncerr(ierr)
      !ierr = nf90_put_att ( ncid, NF90_GLOBAL, "contact", "ccam@csiro.au" )
      !call check_ncerr(ierr)
      !ierr = nf90_put_att ( ncid, NF90_GLOBAL, "project", "Undefined CCAM project" )
      !call check_ncerr(ierr)
      ! Make sure it's a null string, otherwise len_trim won't find the end 
      ! properly.
      histstr = "" 
      call hstring(histstr)
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "history", histstr )
      if ( present(extra_atts) ) then
         do i=1,size(extra_atts)
            if ( save_ccam_parameters .or.                                &
                 extra_atts(i)%name=="driving_model_id" .or.              &
                 extra_atts(i)%name=="driving_model_ensemble_number" .or. &
                 extra_atts(i)%name=="driving_experiment_name" ) then
               select case(extra_atts(i)%type)
               case ( NF90_FLOAT )
                  ierr = nf90_put_att ( ncid, NF90_GLOBAL, extra_atts(i)%name, extra_atts(i)%rval )
               case ( NF90_INT )
                  ierr = nf90_put_att ( ncid, NF90_GLOBAL, extra_atts(i)%name, extra_atts(i)%ival )
               case ( NF90_CHAR )
                  ierr = nf90_put_att ( ncid, NF90_GLOBAL, extra_atts(i)%name, trim(extra_atts(i)%cval) )
               case default
                  print*, "Error, unexpected type in extra_atts", extra_atts(i)%type
                  stop
               end select
               call check_ncerr(ierr, "Error writing extra attribute")
            end if   
         end do
      end if
      
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "modeling_realm", "atmos" )
      
      ! describe grid for CORDEX
      write(gridsize,'(I0)') il
      write(gridschmidt,'(G0.6)') real(nint(1.e5_8/real(schmidt,8))*1.e-5_8)
      do i = len_trim(gridschmidt),1,-1
         if (gridschmidt(i:i) /= "0") exit
      end do
      griddes = "Unrotated latitude/longitude grid interpolated from a variable resolution "// &
                "conformal cubic C"//trim(gridsize)//" grid with Schmidt="//trim(gridschmidt(1:i))
      ierr = nf90_put_att( ncid, NF90_GLOBAL, "grid", griddes )

!     Define attributes for the dimensions
!     Is there any value in keeping long_name for latitude and longitude?
!     Possibly add bounds here??
      if ( areps_compliant ) then
         ierr = nf90_def_var ( ncid, "longitude", NF90_FLOAT, dims%x, dimvars%x) 
      else    
         ierr = nf90_def_var ( ncid, "lon", NF90_FLOAT, dims%x, dimvars%x)
      end if   
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "long_name", "longitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "standard_name", "longitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "axis", "X" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "units", "degrees_east" )
      call check_ncerr(ierr)
      if ( cf_compliant .or. cordex_compliant ) then
         ierr = nf90_put_att ( ncid, dimvars%x, "bounds", "lon_bnds" )
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "lon_bnds", NF90_FLOAT, (/dims%b, dims%x/), dimvars%x_b)
         call check_ncerr(ierr)
      end if

      if ( areps_compliant ) then
         ierr = nf90_def_var ( ncid, "latitude", NF90_FLOAT, dims%y, dimvars%y ) 
      else    
         ierr = nf90_def_var ( ncid, "lat", NF90_FLOAT, dims%y, dimvars%y )
      end if  
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "long_name", "latitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "standard_name", "latitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "axis", "Y" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "units", "degrees_north" )
      call check_ncerr(ierr)
      if ( cf_compliant .or. cordex_compliant ) then
         ierr = nf90_put_att ( ncid, dimvars%y, "bounds", "lat_bnds" )
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "lat_bnds", NF90_FLOAT, (/dims%b, dims%y/), dimvars%y_b)
         call check_ncerr(ierr)
      end if

      if ( nlev > 1 ) then
         if ( use_plevs ) then
            if ( areps_compliant ) then 
               ierr = nf90_def_var ( ncid, "level", NF90_FLOAT, dims%z, dimvars%z ) 
            else    
               ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
            end if   
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "pressure_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "hPa" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "down" )
            call check_ncerr(ierr)
         else if ( use_meters ) then
            if ( areps_compliant ) then  
               ierr = nf90_def_var ( ncid, "level", NF90_FLOAT, dims%z, dimvars%z ) 
            else    
               ierr = nf90_def_var ( ncid, "alt", NF90_FLOAT, dims%z, dimvars%z )
            end if   
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "standard_name", "height" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "vertical distance above the surface" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "m" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "up" )
            call check_ncerr(ierr)
         else if ( use_theta ) then
            if ( areps_compliant ) then  
               ierr = nf90_def_var ( ncid, "level", NF90_FLOAT, dims%z, dimvars%z )
            else
               ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
            end if   
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "theta_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "K" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "up" )
            call check_ncerr(ierr)
         else if ( use_pvort ) then
            if ( areps_compliant ) then 
               ierr = nf90_def_var ( ncid, "level", NF90_FLOAT, dims%z, dimvars%z ) 
            else    
               ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
            end if   
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "potential_vorticity_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "PVU" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "up" )
            call check_ncerr(ierr)             
         else if ( use_hyblevs ) then
            ! This is required for grads to recognise it's a vertical dimension
            ! Allowed by CF convention
            if ( areps_compliant ) then 
               ierr = nf90_def_var ( ncid, "level", NF90_FLOAT, dims%z, dimvars%z ) 
            else    
               ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
            end if   
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "standard_name", "atmosphere_hybrid_sigma_pressure_coordinate" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "formula_terms", "a: anf b: bnf ps: psf p0: P0" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "comment",  "p(k) = anf(k)*P0 + bnf(k)*ps, Coordinate value is anf(k)+ bnf(k)")
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "down" )
            call check_ncerr(ierr)
         else
            if ( areps_compliant ) then 
               ierr = nf90_def_var ( ncid, "level", NF90_FLOAT, dims%z, dimvars%z ) 
            else    
               ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
            end if   
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "sigma_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "sigma_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "down" )
            call check_ncerr(ierr)
         end if
         ierr = nf90_put_att ( ncid, dimvars%z, "axis", "Z" )
         call check_ncerr(ierr)
      end if
      
      if ( ol > 0 ) then
         if ( use_depth ) then 
            ierr = nf90_def_var ( ncid, "olev", NF90_FLOAT, dims%oz, dimvars%oz )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%oz, "long_name", "ocean depth" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%oz, "units", "m" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%oz, "positive", "down" )
            call check_ncerr(ierr)
         else    
            ierr = nf90_def_var ( ncid, "olev", NF90_FLOAT, dims%oz, dimvars%oz )
            call check_ncerr(ierr)
            if ( osig_found ) then
               ierr = nf90_put_att ( ncid, dimvars%oz, "long_name", "ocean sigma_level" )
               call check_ncerr(ierr)
               ierr = nf90_put_att ( ncid, dimvars%oz, "units", "sigma_level" )
               call check_ncerr(ierr)
            else
               ierr = nf90_put_att ( ncid, dimvars%oz, "long_name", "ocean zstar_level" )
               call check_ncerr(ierr)
               ierr = nf90_put_att ( ncid, dimvars%oz, "units", "m" )
               call check_ncerr(ierr)                
            end if   
            ierr = nf90_put_att ( ncid, dimvars%oz, "positive", "down" )
            call check_ncerr(ierr)
         end if   
      end if

      if ( present(nsoil) ) then
         if ( nsoil>0 ) then
            ierr = nf90_def_var ( ncid, "depth", NF90_FLOAT, dims%zsoil, dimvars%zsoil)
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%zsoil, "long_name", "depth of soil layers" )
            ierr = nf90_put_att ( ncid, dimvars%zsoil, "standard_name", "depth" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%zsoil, "units", "m" )
            call check_ncerr(ierr)
            if ( cf_compliant .or. cordex_compliant ) then
               ierr = nf90_put_att ( ncid, dimvars%zsoil, "bounds", "depth_bnds" )
               call check_ncerr(ierr)
               ierr = nf90_def_var ( ncid, "depth_bnds", NF90_FLOAT, (/dims%b, dims%zsoil/), dimvars%zsoil_b)
               call check_ncerr(ierr)
            end if
         end if
      end if

      if ( cptch > 0 ) then
         ierr = nf90_def_var ( ncid, "cable_patch", NF90_FLOAT, dims%cptch, dimvars%cptch )
         call check_ncerr(ierr)
      end if
      if ( cchrt > 0 ) then
         ierr = nf90_def_var ( ncid, "cable_cohort", NF90_FLOAT, dims%cchrt, dimvars%cchrt )
         call check_ncerr(ierr)
      end if

      if ( basetime /= "none" ) then
         ierr = nf90_def_var ( ncid, "time", NF90_FLOAT, dims%t, dimvars%t )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, dimvars%t, "units", basetime )
         call check_ncerr(ierr)
         if ( present(calendar) ) then
            ierr = nf90_put_att ( ncid, dimvars%t, "calendar", calendar )
            call check_ncerr(ierr)
         end if
         ierr = nf90_put_att ( ncid, dimvars%t, "standard_name", "time" )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, dimvars%t, "axis", "T" )
         call check_ncerr(ierr)
         if ( (cf_compliant.or.cordex_compliant) .and. do_time_bnds ) then
            ierr = nf90_put_att ( ncid, dimvars%t, "bounds", "time_bnds" )
            call check_ncerr(ierr)
            ierr = nf90_def_var ( ncid, "time_bnds", NF90_FLOAT, (/dims%b, dims%t/), dimvars%t_b)
            call check_ncerr(ierr)
         end if
      end if   

      if ( nlev>1 .and. use_hyblevs ) then
         ierr = nf90_def_var ( ncid, "anf", NF90_FLOAT, dims%z, vid )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "long_name",  "A coefficient for hybrid level calculation")
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "bnf", NF90_FLOAT, dims%z, vid )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "long_name",  "B coefficient for hybrid level calculation")
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "P0", NF90_FLOAT, vid )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "long_name",  "Reference pressure for hybrid level calculation")
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "units",  "hPa")
         call check_ncerr(ierr)
      end if

      ! Define variables for coordinate heights
      do k = 1,size(coord_name)
         if ( coord_name(k) /= "" ) then 
            ierr = nf90_def_var ( ncid,  coord_name(k), NF90_FLOAT, vid )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, vid, "standard_name",  coord_stdname(k))
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, vid, "long_name",  coord_stdname(k))
            call check_ncerr(ierr)            
            ierr = nf90_put_att ( ncid, vid, "units",  coord_units(k))
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, vid, "axis",  "Z")
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, vid, "positive",  coord_positive(k))
            call check_ncerr(ierr)
         end if   
      end do
      
      if ( cordex_compliant ) then
         ! Define grid mapping for CF-1.11
         ierr = nf90_def_var ( ncid, "crs", NF90_INT, vid )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "grid_mapping_name", "latitude_longitude" )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "semi_major_axis", 6371000. )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "inverse_flattening", 0. )
         call check_ncerr(ierr)      
      end if   

   end subroutine create_ncfile
   
!-----------------------------------------------------------------------------
   subroutine closehist
      
      use logging_m
      
      integer :: ierr, i
     
      call START_LOG(closehist_begin)
      
!     The hist_oave files are closed individually in writehist.
      do i = 1,size(histfile)
         ierr = nf90_close ( histfile(i)%id )
         call check_ncerr(ierr,"Error closing history file")
      end do
      deallocate( histfile, histdimvars )
      
      call END_LOG(closehist_end)
      
   end subroutine closehist
!-------------------------------------------------------------------

   subroutine gsavehist2D ( name, array )
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:,:) :: array
      real, dimension(size(array,1),size(array,2),1) :: temp

      call START_LOG(savehist_begin)      

!     Copy to the 3D array that savehist_work expects. The extra level of
!     subroutines is necessary to avoid multiple versions.
      temp(:,:,1) = array
      call savehist_work  ( name, temp, size(array,2) )

      call END_LOG(savehist_end)

   end subroutine gsavehist2D
   
   subroutine gsavehist3D( name, array )
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:,:,:) :: array
      real, dimension(size(array,1),size(array,2),size(array,3)) :: temp

      call START_LOG(savehist_begin)

      call savehist_work ( name, array, size(array,2) )

      call END_LOG(savehist_end)

   end subroutine gsavehist3D
   
   subroutine gsavehist4D( name, array )
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:,:,:,:) :: array
      real, dimension(size(array,1),size(array,2),size(array,3)*size(array,4)) :: temp

      call START_LOG(savehist_begin)

      temp(:,:,:) = reshape(array,(/ size(array,1), size(array,2), size(array,3)*size(array,4) /) )
      call savehist_work ( name, temp, size(array,2) )

      call END_LOG(savehist_end)

   end subroutine gsavehist4D
   
   subroutine savehist_work ( name, array, jlat2 )

!     This routine actually saves the data to the history array.
      character(len=*), intent(in) :: name
      real, dimension(:,:,:) :: array
      integer, intent(in) :: jlat2
      integer :: ifld, js, jn, nl, istart, iend
      integer :: nx
      
!     Check that openhist has been called to allocate the array
      if ( .not. allocated ( histarray ) ) then
         print*, " History error - must call openhist before savehist "
         stop
      end if

!     Loop over history files because variable could occur in several.

      ifld = bindex_hname(name,inames(:,1:totflds),totflds)
      if ( ifld == 0 ) then
         print*, "Error - history variable ", name, " is not known. "
         stop
      end if
      if ( hist_debug > 5 ) then
         print*, "SAVEHIST2D ", name, ifld, &
                 histinfo(ifld)%used
      end if
      if ( .not. histinfo(ifld)%used ) then
         return
      end if

!     Optional check on array sizes
      if ( hist_debug > 0 ) then
         if ( size(histarray,1) /= size(array,1) ) then
            print*, "ERROR: History array has nlon = ", &
                    size(histarray,1), &
                    "Input array has nlon = ", size(array,1)
            print*, "Variable ", trim(name)
            stop
         end if
      end if

      if ( jlat2 > size(histarray,2) ) then
         print*, "jlat out of bounds ", jlat2
         print*, "Variable ", trim(name)
         stop
      end if

      nl = histinfo(ifld)%nlevels
      if ( nl /= size(array,3) ) then
         print*, "History error, savehist called with array ", name
         print*, "Input has ", size(array,3) , &
                 " levels, history array has ", nl
         print*, "Variable ", trim(name)
         stop
      end if

      istart = histinfo(ifld)%ptr
      iend = istart + nl - 1

      select case ( histinfo(ifld)%ave_type)
      case ( hist_ave )
         where ( array /= nf90_fill_float .and. histarray(:,1:jlat2,istart:iend) /= nf90_fill_float ) 
            histarray(:,1:jlat2,istart:iend) =  &
                 histarray(:,1:jlat2,istart:iend) + array
         elsewhere ( array /= nf90_fill_float )
            histarray(:,1:jlat2,istart:iend) = array
         end where    
      case ( hist_max ) 
         where ( array /= nf90_fill_float .and. histarray(:,1:jlat2,istart:iend) /= nf90_fill_float )  
            histarray(:,1:jlat2,istart:iend) =  &
                 max ( histarray(:,1:jlat2,istart:iend), array )
         elsewhere ( array /= nf90_fill_float )
            histarray(:,1:jlat2,istart:iend) = array
         end where 
      case ( hist_min )
         where ( array /= nf90_fill_float .and. histarray(:,1:jlat2,istart:iend) /= nf90_fill_float )     
            histarray(:,1:jlat2,istart:iend) =  &
                 min ( histarray(:,1:jlat2,istart:iend), array )
         elsewhere ( array /= nf90_fill_float )
            histarray(:,1:jlat2,istart:iend) = array
         end where 
      case ( hist_inst, hist_fixed ) 
         histarray(:,1:jlat2,istart:iend) = array
      case default
         print*, "Internal error in history, unknown ave_type", &
              histinfo(ifld)%ave_type
         stop
      end select
     
      if ( any( array /= nf90_fill_float ) ) then
         histinfo(ifld)%count = histinfo(ifld)%count + 1
      end if   

      if ( hist_debug >= 4 ) then
         if ( jhdb==1 ) then
            if ( khdb <= nl ) then
               print*, "History at point ", name,  &
                 histarray(ihdb,1,istart+khdb-1), array(ihdb,1,khdb)
            else
               print*, "History at point ", name,  &
                 histarray(ihdb,1,istart), array(ihdb,1,1)
            end if
         else if ( jlat2 >= jhdb ) then
            ! Global case
            if ( khdb <= nl ) then
               print*, "History at point ", name,  &
                 histarray(ihdb,jhdb,istart+khdb-1), array(ihdb,jhdb,khdb)
            else
               print*, "History at point ", name,  &
                 histarray(ihdb,jhdb,istart), array(ihdb,jhdb,1)
            end if
         end if
      end if

   end subroutine savehist_work

!-------------------------------------------------------------------
   subroutine writehist( istep, time, dtime, endofrun, year, month, interp, time_bnds )

      use mpidata_m
      use newmpar_m, only : cptch, cchrt
      use logging_m
#ifdef usempimod
      use mpi
#else
      include 'mpif.h'
#endif
      
      integer, intent(in) :: istep
      logical, intent(in), optional :: endofrun
      integer, intent(in), optional :: year, month !  For old format files
      real, intent(in) :: time, dtime
      real, dimension(2), intent(in), optional :: time_bnds
      interface
         subroutine interp ( a, b, flag )
            real, intent(in), dimension(:,:) :: a
            real, intent(out), dimension(:,:) :: b
            integer, intent(in) :: flag
         end subroutine interp
      end interface
      optional :: interp

      integer ierr, vid, ifld
      integer ip, n, i
      integer, dimension(5) :: start4D, count4D
      integer, dimension(4) :: start3D, count3D
      integer, dimension(3) :: start2D, count2D
      integer :: istart, iend
      integer :: nlev, ncount, ncid
      integer :: k
      real :: umin, umax, addoff, sf, timeout
!     Temporary for interpolated output
      integer, dimension(:), allocatable, save :: k_indx
      real, dimension ( nxhis, nyhis ) :: htemp
      real, dimension(:,:), allocatable, save :: hist_g
      real, dimension(2) :: time_bnds_local
      logical :: doinc
      logical :: endofday, endof6hr
      integer, save :: safe_count = 0
      integer :: cnt, maxcnt, slab, gap, offset
      integer :: rrank, sizehis, itag
      integer :: nreq, dreq, rreq, maxdreq, maxrreq
      integer, dimension(:), allocatable, save :: ireq, ireq_r, ireq_d
      real, dimension(:,:,:), allocatable, save :: htemp_buff
      logical :: do_time_bnds, do_adj_avetime
      real :: adj_day, adj_6hr
      
      if ( require_interp .and. .not. present(interp) ) then
         print*, " Error, interp argument required for writehist "
         stop
      end if

!     Check whether this file should be written on this step.
      if ( hist_debug > 1 ) then
         print*, "Checking frequency ", istep, hfreq
      end if

!     For files containing time averages, calculate the appropriate
!     average time value. Don't increment this on the final call when
!     endofrun is set.
      doinc = .true.
      if ( present(endofrun) ) then
         doinc = .not. endofrun
      end if
      if ( doinc .and. ihtype == hist_ave ) then
         avetime = avetime + time
         timecount = timecount + 1
         if ( cordex_compliant ) then
            avetime = avetime - 0.5*dtime
         end if
         if ( hist_debug > 0 ) then
            print*, "AVETIME ", time, avetime, timecount
         end if
      end if
      if ( doinc .and. ihtype == hist_ave .and. present(time_bnds) ) then
         ! Bounds for average should be the extreme bounds
         avetime_bnds(1) = min(avetime_bnds(1), time_bnds(1))
         avetime_bnds(2) = max(avetime_bnds(2), time_bnds(2))
         if ( hist_debug > 0 ) then
            print*, "AVETIME_BNDS ", time_bnds, avetime_bnds(:)
         end if
      end if

!     Files with hfreq=0 are written once at end of run
      if ( hfreq == 0 ) then
         if ( present(endofrun) ) then
            if ( .not. endofrun ) then
               return
            end if
         else ! No endofrun argument
            return
         end if
      else   ! hfreq /= 0
         if ( modulo(istep, hfreq ) /= 0 ) then
            return
         end if
!        Files for which hfreq is non-zero shouldn't be written at the end.
         if ( present(endofrun) ) then
            if ( endofrun ) then
               return
            end if
         end if
      end if
     
      call START_LOG(writehist_begin)
      
      if ( hist_debug > 0 ) then
         print*, "Writing history file "
      end if

      histset = histset + 1
      if ( cf_compliant ) then
         endofday = modulo(nint(time*1440.),1440) == 0  
         endof6hr = modulo(nint(time*1440.),360) == 0
         adj_day = 1.         
         adj_6hr = 0.25
      else if ( areps_compliant ) then
         endofday = modulo(nint(time*60.),1440) == 0  
         endof6hr = modulo(nint(time*60.),360) == 0
         adj_day = 24.         
         adj_6hr = 6.
      else    
         endofday = modulo(nint(time),1440) == 0  
         endof6hr = modulo(nint(time),360) == 0
         adj_day = 1440.         
         adj_6hr = 360.
      end if   
      if ( endofday ) histset_daily = histset_daily + 1  
      if ( endof6hr ) histset_6hr = histset_6hr + 1

      ! write time
      call START_LOG(putvar_begin)
      do i = 1,size(histfile)
          
         do_time_bnds = (cf_compliant.or.cordex_compliant) .and. present(time_bnds) .and. &
                        .not.(histfile(i)%inst.and.ihtype==hist_inst) .and.               &
                        .not.(ihtype==hist_fixed)
         do_adj_avetime = cordex_compliant .and. .not.histfile(i)%inst .and. &
                          ihtype==hist_inst
          
         if ( histfile(i)%fix ) then
            ! do nothing 
         else if ( histfile(i)%day ) then
            ! output only contains daily data 
            if ( endofday ) then
               ncid = histfile(i)%id
               ierr = nf90_inq_varid (ncid, "time", vid )
               call check_ncerr(ierr, "Error getting time id")       
               if ( ihtype == hist_ave) then
                  timeout = avetime/timecount
               else
                  timeout = time 
               end if
               if ( do_adj_avetime ) then
                  timeout = timeout - 0.5*adj_day
               end if   
               ierr = nf90_put_var ( ncid, vid, timeout, start=(/histset_daily/) )
               call check_ncerr(ierr, "Error writing time")
               if ( do_time_bnds ) then
                  ierr = nf90_inq_varid (ncid, "time_bnds", vid )
                  call check_ncerr(ierr, "Error getting time_bnds id")
                  if ( ihtype == hist_ave) then
                     ierr = nf90_put_var ( ncid, vid, avetime_bnds(:), start=(/1,histset_daily/) )
                  else
                     time_bnds_local(2) = time_bnds(2)
                     time_bnds_local(1) = time_bnds(2) - adj_day
                     ierr = nf90_put_var ( ncid, vid, time_bnds_local, start=(/1,histset_daily/) )
                  end if
                  call check_ncerr(ierr, "Error writing time_bnds")
               end if   
            end if 
         else if ( histfile(i)%sixhr ) then
            ! output only contains 6-hourly data 
            if ( endof6hr ) then
               ncid = histfile(i)%id
               ierr = nf90_inq_varid (ncid, "time", vid )
               call check_ncerr(ierr, "Error getting time id")       
               if ( ihtype == hist_ave) then
                  timeout = avetime/timecount
               else
                  timeout = time
               end if
               if ( do_adj_avetime ) then
                  timeout = timeout - 0.5*adj_6hr
               end if
               ierr = nf90_put_var ( ncid, vid, timeout, start=(/histset_6hr/) )
               call check_ncerr(ierr, "Error writing time")
               if ( do_time_bnds ) then
                  ierr = nf90_inq_varid (ncid, "time_bnds", vid )
                  call check_ncerr(ierr, "Error getting time_bnds id")
                  if ( ihtype == hist_ave) then
                     ierr = nf90_put_var ( ncid, vid, avetime_bnds(:), start=(/1,histset_6hr/) )
                  else
                     time_bnds_local(2) = time_bnds(2)
                     time_bnds_local(1) = time_bnds(2) - adj_6hr
                     ierr = nf90_put_var ( ncid, vid, time_bnds_local, start=(/1,histset_6hr/) )
                  end if
                  call check_ncerr(ierr, "Error writing time_bnds")
               end if   
            end if  
         else
            ! output contains instantaneous or mixed data 
            ncid = histfile(i)%id
            ierr = nf90_inq_varid (ncid, "time", vid )
            call check_ncerr(ierr, "Error getting time id")
            if ( ihtype == hist_ave) then
               timeout = avetime/timecount
            else
               timeout = time
            end if
            if ( do_adj_avetime ) then
               timeout = timeout - 0.5*dtime 
            end if
            ierr = nf90_put_var ( ncid, vid, timeout, start=(/histset/) )
            call check_ncerr(ierr, "Error writing time")
            if ( do_time_bnds ) then
               ierr = nf90_inq_varid (ncid, "time_bnds", vid )
               call check_ncerr(ierr, "Error getting time_bnds id")
               if ( ihtype == hist_ave) then
                  ierr = nf90_put_var ( ncid, vid, avetime_bnds(:), start=(/1,histset/))
               else
                  ierr = nf90_put_var ( ncid, vid, time_bnds, start=(/1,histset/))
               end if
               call check_ncerr(ierr, "Error writing time_bnds")
            end if   
         end if
      end do
      call END_LOG(putvar_end)

      cnt = 0
!     first pass
!     find total number of levels to store
      do ifld = 1,totflds

         if ( .not.writeout_test(ifld,endof6hr,endofday) ) then
         
            nlev = histinfo(ifld)%nlevels
            ncount = histinfo(ifld)%count
            istart = histinfo(ifld)%ptr
            iend = istart + nlev - 1
            
            if ( histinfo(ifld)%ave_type == hist_ave ) then    
               if ( hist_debug >= 4 ) then
                  print*, "Raw history at point ", histinfo(ifld)%name, &
                    histarray(ihdb,jhdb,istart+khdb-1), ncount
               end if
               where ( histarray(:,:,istart:iend) /= NF90_FILL_FLOAT )
                  histarray(:,:,istart:iend) =   &
                    histarray(:,:,istart:iend) / max( ncount, 1 )
               end where  
            end if
            if ( histinfo(ifld)%output_scale /= 0 ) then
               where ( histarray(:,:,istart:iend) /= NF90_FILL_FLOAT ) 
                  histarray(:,:,istart:iend) = histarray(:,:,istart:iend) * &
                    histinfo(ifld)%output_scale
               end where   
            end if
            if ( hist_debug >= 4 ) then
               print*, "History written at point ", histinfo(ifld)%name, &
                    histarray(ihdb,jhdb,istart+khdb-1)
            end if

!           Even multilevel variables are written one level at a time
            cnt = cnt + iend - istart + 1
            
         end if   
               
      end do ! Loop over fields
      
      ! allocate memory for maxcnt number of interpolation arrays 
      maxcnt = cnt
      if ( maxcnt > 0 ) then
      
         slab = ceiling(real(maxcnt,8)/real(nproc,8))
         gap = max( nproc/maxcnt, 1 )
         if ( slab > 1 .and. gap > 1 ) then
            print *, "Internal error defining slab and gap ",slab,gap
            stop
         end if
         if ( interpolate_test(myid,slab,gap,maxcnt) ) then
            allocate( hist_a(pil,pjl*pnpan,pnproc,1+slab*myid/gap:slab*(myid/gap+1)) ) 
            if ( allocated( hist_g ) ) then
               if ( size(hist_g,1) /= nx_g .or. size(hist_g,2) /= ny_g ) then
                  deallocate( hist_g )
               end if
            end if
            if ( .not.allocated( hist_g ) ) then
               allocate( hist_g(nx_g,ny_g) )
            end if   
         else
            allocate( hist_a(0,0,0,0) )  
         end if
         if ( allocated( k_indx ) ) then
            if ( size(k_indx) < maxcnt ) then
               deallocate( k_indx, ireq_r, ireq_d )
            end if
         end if
         if ( .not.allocated( k_indx ) ) then
            allocate( k_indx(maxcnt), ireq_r(maxcnt), ireq_d(maxcnt) )
         end if

         dreq = 0
         cnt = 0
!        second pass
!        create the k_indx array used to index histarray
!        count number of dreq messages to be sent from interpolation process to
!        the file output process
         do ifld = 1,totflds
             
            if ( .not.writeout_test(ifld,endof6hr,endofday) ) then 
               nlev = histinfo(ifld)%nlevels
               istart = histinfo(ifld)%ptr
               iend = istart + nlev - 1

!              Even multilevel variables are written one level at a time
               do k = istart,iend
                  cnt = cnt + 1
                  k_indx(cnt) = k
                  rrank = ((cnt-1)*gap)/slab
                  if ( rrank==histinfo(ifld)%procid .or. &
                       myid==histinfo(ifld)%procid .or.  &
                       myid==rrank ) then
                     dreq = dreq + 1 
                  end if
               end do   ! k loop
               
            end if   

         end do ! Loop over fields
      
         ! reallocate MPI request and buffers if more memory is needed
         maxdreq = dreq
         if ( allocated( ireq ) ) then
            if ( size(ireq) < maxdreq ) then
               deallocate( htemp_buff, ireq )
            end if
         end if
         if ( .not.allocated( ireq ) ) then
            allocate( htemp_buff(nxhis,nyhis,maxdreq), ireq(maxdreq) )
         end if      

!        now gather data on interpolation processes
         if ( slab > 0 ) then
            call gather_wrap(histarray,hist_a,slab,gap,maxcnt,k_indx)
         end if
 
         cnt = 0
         nreq = 0
         dreq = 0
         sizehis = nxhis*nyhis
!        third pass
!        perform the interpolation and send data to output process
         do ifld = 1,totflds
             
            if ( .not.writeout_test(ifld,endof6hr,endofday) ) then

               nlev = histinfo(ifld)%nlevels
               istart = histinfo(ifld)%ptr
               iend = istart + nlev - 1

!              Even multilevel variables are written one level at a time
               do k = istart, iend

                  cnt = cnt + 1
                  rrank = ((cnt-1)*gap)/slab
                  if ( myid == rrank ) then
                     do ip = 0,pnproc-1
                        do n = 0,pnpan-1
                           hist_g(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g) = &
                              hist_a(1:pil,1+n*pjl:(n+1)*pjl,ip+1,cnt)
                        end do
                     end do
!                    Update fill if required         
                     if ( histinfo(ifld)%fill ) then
                        call fill_cc1(hist_g)
                     end if
                     if ( present(interp) ) then
                        call interp ( hist_g, htemp, histinfo(ifld)%int_type )
                     else
                        htemp = hist_g(:,:)
                     end if
                     if ( histinfo(ifld)%all_positive ) then
                        htemp = max( htemp, 0. )
                     end if   
                  end if

                  call START_LOG(mpisendrecv_begin)
                  itag = 1000 + cnt
                  rrank = ((cnt-1)*gap)/slab
                  if ( rrank /= histinfo(ifld)%procid ) then
                     if ( myid == histinfo(ifld)%procid ) then
                        dreq = dreq + 1 
                        nreq = nreq + 1
                        call MPI_IRecv(htemp_buff(:,:,dreq), sizehis, MPI_REAL, rrank, itag, comm_world, ireq(nreq), ierr)
                        ireq_r(cnt) = nreq
                        ireq_d(cnt) = dreq
                     else if ( myid == rrank ) then
                        dreq = dreq + 1 
                        nreq = nreq + 1 
                        htemp_buff(:,:,dreq) = htemp
                        call MPI_ISend(htemp_buff(:,:,dreq), sizehis, MPI_REAL, histinfo(ifld)%procid, itag, comm_world, &
                                       ireq(nreq), ierr)
                     end if
                  else
                     dreq = dreq + 1 
                     htemp_buff(:,:,dreq) = htemp
                     ireq_d(cnt) = dreq
                  end if
                  call END_LOG(mpisendrecv_end)
            
               end do   ! k loop
               
            end if   

         end do ! Loop over fields
      
         cnt = 0
!        fourth pass
!        write the data from output process
         do ifld = 1,totflds
             
            if ( .not.writeout_test(ifld,endof6hr,endofday) ) then 

               nlev = histinfo(ifld)%nlevels
               vid = histinfo(ifld)%vid
               istart = histinfo(ifld)%ptr
               iend = istart + nlev - 1

!              Even multilevel variables are written one level at a time
               do k = istart, iend

                  cnt = cnt + 1

                  if ( myid == histinfo(ifld)%procid ) then
            
                     call START_LOG(mpisendrecv_begin) 
                     rrank = ((cnt-1)*gap)/slab
                     if ( rrank /= histinfo(ifld)%procid ) then
                        rreq = ireq_r(cnt) 
                        call MPI_Wait(ireq(rreq), MPI_STATUS_IGNORE, ierr)
                     end if   
                     dreq = ireq_d(cnt) 
                     htemp = htemp_buff(:,:,dreq)
                     call END_LOG(mpisendrecv_end)

                     if ( hbytes == 2 ) then
                        addoff = histinfo(ifld)%addoff
                        sf = histinfo(ifld)%scalef
                        umin = sf * vmin + addoff
                        umax = sf * vmax + addoff
                        where ( fpequal(htemp, nf90_fill_float) )
                           htemp = NF90_FILL_SHORT
                        elsewhere
                           ! Put the scaled array back in the original and let
                           ! netcdf take care of conversion to int2
                           htemp = nint((max(umin,min(umax,htemp))-addoff)/sf)
                        end where
                     else
                        where ( fpequal(htemp, nf90_fill_float) )
                           htemp = missing_value_cordex 
                        end where    
                     end if    
                   
                     call START_LOG(putvar_begin)
                     ncid = histinfo(ifld)%ncid
                     if ( histinfo(ifld)%daily .and. .not.single_output ) then
                        if ( endofday ) then
                           if ( histinfo(ifld)%pop4d ) then
                              start4D = (/ 1, 1, mod(k+1-istart-1,cptch)+1, 1+(k+1-istart-1)/cptch,  histset_daily /)
                              count4D = (/ nxhis, nyhis, 1, 1, 1 /)
                              ierr = nf90_put_var ( ncid, vid, htemp, start=start4D, count=count4D )
                           else if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                              start3D = (/ 1, 1, k+1-istart, histset_daily /)
                              count3D = (/ nxhis, nyhis, 1, 1 /)
                              ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
                           else
                              start2D = (/ 1, 1, histset_daily /)
                              count2D = (/ nxhis, nyhis, 1 /)
                              ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
                           end if
                        end if    
                     else if ( histinfo(ifld)%sixhr .and. .not.single_output ) then
                        if ( endof6hr ) then
                           if ( histinfo(ifld)%pop4d ) then
                              start4D = (/ 1, 1, mod(k+1-istart-1,cptch)+1, 1+(k+1-istart-1)/cptch,  histset_6hr /)
                              count4D = (/ nxhis, nyhis, 1, 1, 1 /)
                              ierr = nf90_put_var ( ncid, vid, htemp, start=start4D, count=count4D )
                           else if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                              start3D = (/ 1, 1, k+1-istart, histset_6hr /)
                              count3D = (/ nxhis, nyhis, 1, 1 /)
                              ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
                           else
                              start2D = (/ 1, 1, histset_6hr /)
                              count2D = (/ nxhis, nyhis, 1 /)
                              ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
                           end if
                        end if  
                     else    
                        if ( histinfo(ifld)%pop4d ) then
                           start4D = (/ 1, 1, mod(k+1-istart-1,cptch)+1, 1+(k+1-istart-1)/cptch,  histset /)
                           count4D = (/ nxhis, nyhis, 1, 1, 1 /)
                           ierr = nf90_put_var ( ncid, vid, htemp, start=start4D, count=count4D )
                        else if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                           start3D = (/ 1, 1, k+1-istart, histset /)
                           count3D = (/ nxhis, nyhis, 1, 1 /)
                           ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
                        else
                           start2D = (/ 1, 1, histset /)
                           count2D = (/ nxhis, nyhis, 1 /)
                           ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
                        end if
                     end if   
                     call check_ncerr(ierr, "Error writing history variable "//histinfo(ifld)%name )
                     call END_LOG(putvar_end)

                  end if
               end do   ! k loop
               
!              Zero array ready for next set
               histarray(:,:,istart:iend) = nf90_fill_float
!              Reset the count variable
               histinfo(ifld)%count = 0

            end if   
               
         end do ! Loop over fields
      
         call START_LOG(mpisendrecv_begin)
         if ( nreq > 0 ) then
            call MPI_Waitall(nreq, ireq, MPI_STATUSES_IGNORE, ierr )
         end if   
         call END_LOG(mpisendrecv_end)       

         deallocate(hist_a)
         
      end if ! maxcnt>0   

      avetime = 0.0
      ! Initialisation so min/max work
      avetime_bnds(:) = (/huge(1.), -huge(1.)/)
      timecount = 0

#ifdef safe
      safe_count = safe_count + 1
      if ( safe_count >= safe_max ) then
         if ( myid == 0 ) then
            write(*,"(a)") "Clearing buffers"
         end if    
!        Sync the file so that if the program crashes for some reason 
!        there will still be useful output.
         call START_LOG(putvar_begin)
         do i = 1,size(histfile) 
            ncid = histfile(i)%id 
            ierr = nf90_sync ( ncid )
            call check_ncerr(ierr, "Error syncing history file")
         end do 
         call END_LOG(putvar_end)
         call START_LOG(mpibarrier_begin)
         call mpi_barrier(comm_world,ierr)
         call END_LOG(mpibarrier_end)
         safe_count = 0
      end if
#endif

      call END_LOG(writehist_end)

   end subroutine writehist

!-------------------------------------------------------------------

   subroutine clearhist()

      ! Reset the history outside the normal process

      integer :: ifld
      integer :: istart, iend
      integer :: nlev

      if ( hist_debug > 0 ) then
         print*, "Resetting history"
      end if

      do ifld = 1,totflds
         if ( .not. histinfo(ifld)%used ) then
            cycle
         end if
         nlev = histinfo(ifld)%nlevels
         istart = histinfo(ifld)%ptr
         iend = istart + nlev - 1

!        Zero ready for next set
         histarray(:,:,istart:iend) = nf90_fill_float
!        Reset the count variable
         histinfo(ifld)%count = 0

      end do ! Loop over fields

   end subroutine clearhist

   subroutine gather_wrap(histarray,hist_a,slab,gap,maxcnt,k_indx)
      use mpidata_m, only : nproc, lproc, myid, pil, pjl, pnpan, comm_world
      use logging_m
#ifdef usempimod
      use mpi
#endif
#ifndef usempimod
      include 'mpif.h'
#endif
      integer, intent(in) :: slab, gap, maxcnt
      integer :: istart, iend, ip, k, n, ierr, lsize, lp, iq
      integer, dimension(maxcnt), intent(in) :: k_indx
      real, dimension(:,:,:), intent(in) :: histarray
      real, dimension(:,:,:,:), intent(inout) :: hist_a
      real, dimension(pil*pjl*pnpan*lproc*slab*nproc) :: hist_a_tmp
      integer :: nreq, rreq, ipr, sreq
      integer :: ldone, rcount, jproc
      integer, parameter :: itag = 1
      integer, dimension(2*nproc) :: ireq
      integer, dimension(2*nproc) :: donelist
      real, dimension(pil,pjl*pnpan,lproc,slab,nproc) :: histarray_tmp      
      
      call START_LOG(gatherwrap_begin)

      nreq = 0
      rreq = 0
      
      istart = 1 + slab*(myid/gap) 
      iend = min( istart + slab - 1, maxcnt )
      if ( interpolate_test(myid,slab,gap,maxcnt) ) then
         lsize = pil*pjl*pnpan*lproc*(iend-istart+1)
         do ipr = 0,nproc-1
            nreq = nreq + 1
            call START_LOG(mpigather_begin)
            call MPI_IRecv( hist_a_tmp(lsize*ipr+1:lsize*(ipr+1)), lsize, MPI_REAL, ipr, &
                            itag, comm_world, ireq(nreq), ierr )
            call END_LOG(mpigather_end)
         end do    
      end if
      rreq = nreq
      
      do ip = 0,nproc-1
         istart = 1 + slab*(ip/gap) 
         iend = min( istart + slab - 1, maxcnt )
         if ( interpolate_test(ip,slab,gap,maxcnt) ) then
            do k = istart,iend
               histarray_tmp(:,:,:,k-istart+1,ip+1) = reshape( histarray(:,:,k_indx(k)), (/ pil,pjl*pnpan,lproc /) )
            end do   
            lsize = pil*pjl*pnpan*lproc*(iend-istart+1)
            nreq = nreq + 1
            call START_LOG(mpigather_begin)
            call MPI_ISend( histarray_tmp(:,:,:,:,ip+1), lsize, MPI_REAL, ip, &
                            itag, comm_world, ireq(nreq), ierr )
            call END_LOG(mpigather_end)
         end if
      end do
      
      istart = 1 + slab*(myid/gap)
      iend = min( slab*(myid/gap+1), maxcnt )          
      if ( interpolate_test(myid,slab,gap,maxcnt) ) then
         rcount = rreq
         do while ( rcount > 0 )
            call START_LOG(mpigather_begin)
            call MPI_Waitsome( rreq, ireq, ldone, donelist, MPI_STATUSES_IGNORE, ierr )
            call END_LOG(mpigather_end)          
            rcount = rcount - ldone
            do jproc = 1,ldone
               n = donelist(jproc)
               do k = istart,iend 
                  do lp = 0,lproc-1 
                     iq = lp*pil*pjl*pnpan + (k-istart)*pil*pjl*pnpan*lproc + (n-1)*pil*pjl*pnpan*lproc*(iend-istart+1)
                     hist_a(:,:,lp+(n-1)*lproc+1,k-istart+1) = reshape( hist_a_tmp(iq+1:iq+pil*pjl*pnpan), (/ pil, pjl*pnpan /) )
                  end do
               end do
            end do   
         end do
      end if 
      
      call START_LOG(mpigather_begin)
      sreq = nreq - rreq
      if ( sreq > 0 ) then
         call MPI_Waitall( sreq, ireq(rreq+1:nreq), MPI_STATUSES_IGNORE, ierr )
      end if   
      call END_LOG(mpigather_end)
      
      call END_LOG(gatherwrap_end)
   
   end subroutine gather_wrap

   pure function writeout_test(ifld,endof6hr,endofday) result(ans)

      ! determine if output should be written for this time-step

      integer, intent(in) :: ifld
      logical, intent(in) :: endof6hr, endofday

      integer :: ave_type, nlev
      logical :: ans

      ave_type = histinfo(ifld)%ave_type
      nlev = histinfo(ifld)%nlevels
      ans = .false.

      if ( .not. histinfo(ifld)%used ) then
         ans = .true.
      end if
!     Only write fixed variables in the first history set
      if ( histset > 1 .and. ave_type == hist_fixed ) then
         ans = .true.
      end if
      if ( histinfo(ifld)%daily .and. .not.single_output .and. .not.endofday ) then
         ans = .true.
      end if   
      if ( histinfo(ifld)%sixhr .and. .not.single_output .and. .not.endof6hr ) then
         ans = .true.
      end if 

   end function writeout_test
   
   pure function interpolate_test(myid,slab,gap,maxcnt) result(ans)
   
      integer, intent(in) :: myid, gap, slab, maxcnt
      logical :: ans
      
      ans = mod(myid,gap)==0 .and. myid<gap*ceiling(real(maxcnt,8)/real(slab,8))
   
   end function interpolate_test
   
!-------------------------------------------------------------------
   function needfld ( name ) result (needed)
      character(len=*), intent(in) :: name
      logical :: needed
      integer :: ifld

!     Check if name is in the list for any of the history files
      needed = .false.
      ifld = bindex_hname ( name, inames(:,1:totflds), totflds)
      if ( ifld == 0 ) then
         ! Name not known at all in this case
         needed = .false.
      else if ( histinfo(ifld)%used ) then
         needed = .true.
      end if
   end function needfld
   
!-------------------------------------------------------------------
   function bindex_hname(name, table,nflds) result(ifld)
      character(len=*), intent(in)    :: name
      integer, intent(in) :: nflds
      integer ( bigint ), intent(in), dimension(MAX_KEYINDEX,nflds) :: table
      integer :: ifld
      integer :: i, lower, upper, k
      integer (bigint), dimension(MAX_KEYINDEX) :: key
      integer (bigint) :: fac

!  Lookup "key" in "table" of integers using a binary search.
!  This assumes that the list has been sorted but it doesn't test for this.

      key = hashkey ( name )

      ifld = 0
      lower = 1
      upper = nflds
      do
         i = (lower+upper)/2
!         print*, "Search", lower, upper, i, table(lower), table(upper), table(i), key
         if ( hash_lt(key,table(:,i)) ) then
            upper = i-1
         else if ( hash_gt(key,table(:,i)) ) then
            lower = i + 1
         else
            ifld = i
            exit
         end if
         if ( upper < lower ) then
            exit
         end if
      end do

   end function bindex_hname
!-------------------------------------------------------------------
   function hashkey ( name ) result (key)
!  Generate an integer from a string
      character(len=*), intent(in) :: name
      integer (bigint), dimension(MAX_KEYINDEX) :: key
      integer :: i, k
      integer (bigint) :: fac

      key(:) = int(0,bigint)
      do k = 1,MAX_KEYINDEX
         fac = int(1,bigint) 
         ! This makes numerical order the same as alphabetical order
         do i=min(k*MAX_KEYLEN,len_trim(name)),(k-1)*MAX_KEYLEN+1,-1
            ! Netcdf allowed characters for variables names are in 
            ! range 48 to 122 (alphanumeric and _)
            key(k) = key(k) + fac*int(ichar(name(i:i))-48,bigint)
            fac = fac * int(75,bigint) ! 122-48+1 = 75
         end do   
      end do
      if (hist_debug > 2 ) then
         print*, "  HASHKEY ", name, key
      end if
   end function hashkey
   
   function calcfilename( vname, fname ) result (oname)
! Generate output filename
      character(len=*), intent(in) :: vname, fname
      character(len=300) :: oname
      integer :: dirpos
   
      dirpos = scan( fname, "\/", back=.true. )
      if ( dirpos > 0 ) then
        if ( dirpos == len_trim(fname) ) then
           print *,"Error, last character in output filename is a directory delimiter"
           stop
        end if   
        oname = fname(1:dirpos)//trim(vname)//"_"//trim(fname(dirpos+1:))  
      else
        oname = trim(vname)//"_"//trim(fname)  
      end if
      
      
   end function calcfilename
   
   subroutine fill_cc1(b_io)
!     routine fills in interior of an array which has undefined points
      use newmpar_m
      use indices_m
      use logging_m      
      real, dimension(:,:), intent(inout) :: b_io ! input and output array
      real, dimension(ifull) :: b, a
      integer, dimension(0:5) :: imin, imax, jmin, jmax
      integer :: iminb, imaxb, jminb, jmaxb
      integer :: nrem, iq, neighb, i, j
      integer :: n
      real :: av, avx
      
      call START_LOG(fillcc0_begin)
      
      a(1:ifull) = reshape( b_io(1:il,1:jl), (/ ifull /) )
      
      imin(:) = 1
      imax(:) = il
      jmin(:) = 1
      jmax(:) = il
      
      nrem = 1    ! Just for first iteration
!     nrem_gmin used to avoid infinite loops, e.g. for no sice
      do while ( nrem > 0 )
         nrem = 0
         b(:) = a(:)
         do n = 0,5
            iminb = il
            imaxb = 1
            jminb = il
            jmaxb = 1
            do j = jmin(n),jmax(n)
               do i = imin(n),imax(n)
                  iq = i + (j-1)*il + n*il*il
                  if ( a(iq) == spval ) then
                     neighb = 0
                     av = 0.
                     if ( a(i_n(iq)) /= spval ) then
                        neighb = neighb + 1
                        av = av + a(i_n(iq))
                     end if
                     if ( a(i_e(iq)) /= spval ) then
                        neighb = neighb + 1
                        av = av + a(i_e(iq))
                     end if
                     if ( a(i_w(iq)) /= spval ) then
                        neighb = neighb + 1
                        av = av + a(i_w(iq))
                     end if
                     if ( a(i_s(iq)) /= spval ) then
                        neighb = neighb + 1
                        av = av + a(i_s(iq))
                     end if
                     if ( neighb > 0 ) then
                        b(iq) = av/neighb
                        avx = av
                     else
                        nrem = nrem + 1   ! current number of points without a neighbour
                        iminb = min( i, iminb )
                        imaxb = max( i, imaxb )
                        jminb = min( j, jminb )
                        jmaxb = max( j, jmaxb )
                     end if
                  end if
               end do
            end do
            imin(n) = iminb
            imax(n) = imaxb
            jmin(n) = jminb
            jmax(n) = jmaxb
         end do
         a(:) = b(:)
         if ( nrem == ifull ) then
            print*, "Warn in fill_cc - no points defined"
            a(:) = nf90_fill_float
            nrem = 0
         end if
      end do
 
      b_io(1:il,1:jl) = reshape( a(1:ifull), (/ il, jl /) )
      
      call END_LOG(fillcc0_end)
      
   end subroutine fill_cc1
   
   subroutine ccmpi_bcast1s(ldat,host,comm)
#ifdef usempimod
      use mpi
#endif
#ifndef usempimod
      include 'mpif.h'
#endif   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, lsize
      character(len=*), intent(inout) :: ldat
      integer i
      integer(kind=1), dimension(:), allocatable :: dummy

      ! MJT notes - MS Windows MPI_CHARACTER seems broken

      lhost = host
      lcomm = comm
      lsize = len(ldat)
      allocate( dummy(lsize) )
      do i = 1,lsize
         dummy(i) = int(iachar(ldat(i:i)),1)
      end do
      call MPI_Bcast(dummy,lsize,MPI_BYTE,lhost,lcomm,lerr)
      do i = 1,lsize
         ldat(i:i) = achar(dummy(i))
      end do
      deallocate( dummy )
      
   end subroutine ccmpi_bcast1s     
   
   subroutine ccmpi_bcast2s(ldat,host,comm)
#ifdef usempimod
      use mpi
#endif
#ifndef usempimod
      include 'mpif.h'
#endif   
      integer, intent(in) :: host, comm
      integer(kind=4) :: lcomm, lhost, lerr, llen, lsize, lnum
      character(len=*), dimension(:), intent(inout) :: ldat
      integer i, j
      integer(kind=1), dimension(:,:), allocatable :: dummy

      ! MJT notes - MS Windows MPI_CHARACTER seems broken

      lhost = host
      lcomm = comm
      llen = len(ldat(1))
      lnum = size(ldat)
      lsize = llen*lnum
      allocate( dummy(llen,lnum) )
      do j = 1,lnum
         do i = 1,llen
            dummy(i,j) = int(iachar(ldat(j)(i:i)),1)
         end do   
      end do
      call MPI_Bcast(dummy,lsize,MPI_BYTE,lhost,lcomm,lerr)
      do j = 1,lnum
         do i = 1,llen
            ldat(j)(i:i) = achar(dummy(i,j))
         end do   
      end do
      deallocate( dummy )
      
   end subroutine ccmpi_bcast2s
   
end module history