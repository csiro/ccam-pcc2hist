! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2021 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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

!**********************************************************************
!***   WARNING: No user serviceable parts inside.
!***   Opening the case will invalidate warranty.
!***   Electric shock may result.
!***
!***   Seriously, you should not have to change anything in this module.
!***   To add a new diagnostic just call addfld and savefld. If it's 
!***   expensive to calculate you can add a call to needfld to check
!***   if it's required.
!***   If you really do need to change anything here then it's a failure
!***   in my modular design and I'd rather fix it properly.
!**********************************************************************

! Things to do.

! At the moment there's no way to specify an integer field like a land
! sea mask.

! One disadvantage of this scheme is that extra arrays are required to
! hold the instantaneous winds. In the old scheme these are just taken 
! from the standard model arrays required for the SLT.

! Should it be an error to call savehist on a field which is not in the 
! variable list or should it just be a warning. A warning might allow more
! flexible ways of extending cc2hist?

! Have added a nsoil argument to openhist. Might be cleaner to have a routine
! to set up vertical axes. This would be more like how one would do in in
! cdat for example.

module history

#ifdef usenc_mod
   use netcdf
#else
   use netcdf_m
#endif
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

   !public :: set_missval

!  Private internal routines
   private :: bindex_hname, initval, sortlist, create_ncvar,            &
              create_ncfile,                                            &
              hashkey, qindex_hname, savehist_work,                     &
              gsavehist2D, gsavehist3D, gsavehist4D
   !private :: create_oldfile

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
   integer, parameter :: nfmax = 1000

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
      logical            :: ran
!     Controls whether variable is on the t?_pop list of variables
      integer            :: tn
!     Scale factor for output, e.g. to convert rain from mm/step to mm/day.
      real               :: output_scale 
      logical            :: used
      integer            :: vid  ! netCDF variable identifier
      integer            :: ptr  ! Pointer to position in histarray
      integer            :: count      ! Used for scaling
      integer            :: ave_type   ! Type of averaging
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
      ! cell_methods appropriate for variable before any history processing is done
      character(len=30)  :: cell_methods
   end type hinfo

!  For adding extra global attributes
   type, public :: hist_att
      character(len=nf90_max_name)   :: name
      integer :: type
      real :: rval
      integer :: ival
      ! This length isn't really related to the max_name length but this will do.
      character(len=nf90_max_name)   :: cval
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

!  Array for accumulating history data.
   real, dimension(:,:,:), allocatable :: histarray

   type (hinfo), dimension(nfmax) :: histinfo 
!  Number of history records written
   integer :: histset

!   netCDF file ID of history file
   integer :: histid      

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
   
! Save CCAM parameters
   logical, public :: save_ccam_parameters = .true.

!  MPI working arrays
   integer, private, save :: nx_g, ny_g
#ifdef usempi3
   real, dimension(:,:,:,:), pointer, contiguous :: hist_a
#else
   real, dimension(:,:,:,:), allocatable, save, private :: hist_a
#endif

!  Maximum number of CABLE tiles
   integer, parameter :: maxtile=5

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
      integer, intent(in) :: un      ! Unit no to read namelist
!     Controls whether errors in reading the namelist should be fatal.
      logical, intent(in), optional :: iofatal 
      character(len=MAX_NAMELEN) :: htype
      integer :: freq
      integer :: bytes
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
           amipnames
      integer :: i, j, ierr
      logical :: fatal

      ! Save the current value of the namelist variables and then
      ! reset them
      htype = ""
      freq = hfreq; hfreq = 0
      bytes = hbytes; hbytes = 4
      names = hnames; hnames = ""
      namesx = xnames; xnames = ""
      read(un,histnl,iostat=ierr)
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
                     pop2d, pop3d, pop4d,                             &
                     coord_height, cell_methods, ran_type, tn_type )
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
      character(len=*), intent(in), optional :: cell_methods
      logical, intent(in), optional   :: ran_type
      integer, intent(in), optional   :: tn_type

!     Local variables corresponding to the optional arguments
      integer :: atype, ltn
      logical :: lstd, lran
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
      histinfo(totflds)%nlevels      = nlevels
      histinfo(totflds)%std          = lstd
      histinfo(totflds)%ran          = lran
      histinfo(totflds)%tn           = ltn
      histinfo(totflds)%output_scale = scale
      histinfo(totflds)%used         = .FALSE.   ! Value for now
      histinfo(totflds)%vid          = 0
      histinfo(totflds)%ptr          = 0
      histinfo(totflds)%count        = 0
      histinfo(totflds)%ave_type     = atype     ! Default if not over-ridden
      histinfo(totflds)%addoff       = 0.0
      histinfo(totflds)%scalef       = 0.0
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
      if ( present(cell_methods) ) then
         histinfo(totflds)%cell_methods = cell_methods
      else
         histinfo(totflds)%cell_methods = ""
      end if

   end subroutine addfld
   
!-------------------------------------------------------------------
   subroutine openhist ( nx, ny, nl, sig, ol, cptch, cchrt, gosig, suffix, hlon, hlat, basetime, &
                         year, nxout, nyout, source, histfilename,      &
                         pressure, height, depth, extra_atts, hybrid_levels, anf,  &
                         bnf, p0, calendar, nsoil, zsoil )
!
!     Create netCDF history files using information in histinfo array.
!
      use mpidata_m
      use logging_m

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
      type(dimarray) :: dims, dimvars
      integer :: ncid, vid, ivar, istart, iend, ilev
      character(len=80) :: longname, units
      character(len=MAX_NAMELEN) :: vname
      logical :: used, multilev, use_plevs, use_hyblevs, use_meters, use_depth
      integer, dimension(totflds) :: coord_heights
      integer :: kc, ncoords, k, pkl
      logical :: soil_used, water_used, osig_found
      real :: dx, dy
      integer :: i
      real, dimension(:), allocatable :: cabledata

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

      if ( myid == 0 ) then

         if ( size(hlat) /= nyhis ) then
            print*, " Error, mismatch in number of latitudes ", size(hlat), nyhis
            stop
         end if
         if ( size(hlon) /= nxhis ) then
            print*, " Error, mismatch in number of longitudes ", size(hlon), nxhis
            stop
         end if
      
      end if

      if ( hist_debug > 2 ) then
         do ivar=1,totflds
            print*, ivar, histinfo(ivar)%name, histinfo(ivar)%nlevels
         end do
      end if

      call sortlist

      if ( hist_debug > 2 ) then
         print*, " Names after sorting "
         do ivar=1,totflds
            print*, ivar, histinfo(ivar)%name, histinfo(ivar)%nlevels
         end do
      end if

!     Save this as a module variable
      filesuffix = suffix

      do ivar=1,totflds
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
         else if ( hnames(ivar) == "ran" ) then
!           Include RAN requested variables
            histinfo(1:totflds)%used = &
             histinfo(1:totflds)%used .or. histinfo(1:totflds)%ran
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

      if ( myid == 0 ) then

!        The rest of the history files are much simpler with multilevel variables

         soil_used = .false.
         water_used = .false.
         if ( present(histfilename) ) then
            filename = histfilename
         else
            write(filename,"(a,i1,a,a)" ) "hist", trim(suffix), ".nc"
         end if
!        Is it possible to do a masked maxval?
         multilev = .false.
         ncoords = 0
         do ifld=1,totflds
            if ( hist_debug > 4 ) then
               print*, "Checking variable properties", ifld, &
                    histinfo(ifld)%name, histinfo(ifld)%used, &
                    histinfo(ifld)%nlevels, histinfo(ifld)%soil
            end if
            if ( .not. histinfo(ifld)%used ) cycle

            ! From here only considering variables that are used in this file

            multilev = multilev .or. histinfo(ifld)%nlevels > 1 .or. &
                       histinfo(ifld)%multilev

            ! Check if the file has any soil variables
            soil_used = soil_used .or. histinfo(ifld)%soil

             ! Check if the file has any water variables
            water_used = water_used .or. histinfo(ifld)%water
               
            ! Get a list of the coordinate heights if any
            if ( histinfo(ifld)%coord_height > -huge(1.) ) then
               ! Check if it's already in list
               do kc=1,ncoords
                  if ( nint(histinfo(ifld)%coord_height) == coord_heights(kc) ) then
                     exit
                  end if
               end do
               if ( kc > ncoords ) then
                  ! Value not found
                  ncoords = kc
                  coord_heights(ncoords) = histinfo(ifld)%coord_height
               end if
            end if
         end do

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
         if ( ol > 0 ) then
            osig_found = all(gosig<=1.)
         else 
            osig_found = .true.
         end if    
         call create_ncfile ( filename, nxhis, nyhis, size(sig), ol, cptch, cchrt, multilev,               &
              use_plevs, use_meters, use_depth, use_hyblevs, basetime,                       &
              coord_heights(1:ncoords), ncid, dims, dimvars, source, extra_atts, calendar,   &
              nsoil, zsoil, osig_found )
         histid = ncid

         do ifld = 1,totflds
            if ( histinfo(ifld)%used ) then
               call create_ncvar(histinfo(ifld), ncid, dims)
            end if
         end do

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
         if ( cf_compliant ) then
            ! Calculate bounds assuming a regular lat-lon grid
            ! Perhaps have optional arguments for the other cases?
            allocate ( lat_bnds(2,size(hlat)), lon_bnds(2,size(hlon)) )
            ! Check if regular grid
            if ( maxval(hlon(2:)-hlon(:nx-1)) - minval(hlon(2:)-hlon(:nx-1)) < 1e-4*maxval(hlon(2:)-hlon(:nx-1)) ) then
               dx = hlon(2) - hlon(1)
               lon_bnds(1,:) = hlon - 0.5*dx
               lon_bnds(2,:) = hlon + 0.5*dx
               ierr = nf90_put_var ( ncid, dimvars%x_b, lon_bnds )
               call check_ncerr(ierr,"Error writing longitude bounds")
            end if
            if ( maxval(hlat(2:)-hlat(:ny-1)) - minval(hlat(2:)-hlat(:ny-1)) < 1e-4*maxval(hlat(2:)-hlat(:ny-1)) ) then
               dy = hlat(2) - hlat(1)
               lat_bnds(1,:) = hlat - 0.5*dy
               lat_bnds(2,:) = hlat + 0.5*dy
               where ( lat_bnds < -90. ) 
                  lat_bnds = -90.
               end where
               where ( lat_bnds > 90. ) 
                  lat_bnds = 90.
               end where
               ierr = nf90_put_var ( ncid, dimvars%y_b, lat_bnds )
               call check_ncerr(ierr,"Error writing latitude bounds")
            end if
         end if
         if ( multilev ) then
            ierr = nf90_put_var ( ncid, dimvars%z, sig )
            call check_ncerr(ierr,"Error writing levels")
            if ( use_hyblevs ) then
               if ( .not. present(anf) ) then
                  print*, "Error, missing anf argument"
                  stop
               end if
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
         end if
         if ( ol > 0 ) then
            ierr = nf90_put_var ( ncid, dimvars%oz, gosig )
            call check_ncerr(ierr,"Error writing olev")
         end if
         if ( cptch > 0 ) then
            allocate( cabledata(cptch) )
            do i=1,cptch
               cabledata(i) = real(i)
            end do
            ierr = nf90_put_var ( ncid, dimvars%cptch, cabledata )
            call check_ncerr(ierr,"Error writing cable_patch")
            deallocate( cabledata )
         end if
         if ( cchrt > 0 ) then
            allocate( cabledata(cchrt) )
            do i=1,cchrt
               cabledata(i) = real(i)
            end do
            ierr = nf90_put_var ( ncid, dimvars%cchrt, cabledata )
            call check_ncerr(ierr,"Error writing cable_cohort")
            deallocate( cabledata )
         end if
         if ( present(zsoil) .and. present(nsoil) ) then
            if ( nsoil>0 ) then
               ierr = nf90_put_var ( ncid, dimvars%zsoil, zsoil )
               call check_ncerr(ierr,"Error writing depths")
               if ( cf_compliant ) then
                  ! Soil bounds
                  allocate(zsoil_bnds(2, nsoil))
                  zsoil_bnds(1,1) = 0.
                  zsoil_bnds(2,1) = 2.*zsoil(1)
                  do k=1,nsoil
                     ! Levels are middle of layers
                     zsoil_bnds(2,k) = zsoil_bnds(2,k-1) + 2*(zsoil(k)-zsoil_bnds(2,k-1))
                     zsoil_bnds(1,k) = zsoil_bnds(2,k-1)
                  end do
                  ierr = nf90_put_var ( ncid, dimvars%zsoil_b, zsoil_bnds )
                  call check_ncerr(ierr,"Error writing depths")
                  deallocate(zsoil_bnds)
               end if
            end if
          end if

         do kc = 1,ncoords
            if ( coord_heights(kc) < 10 ) then
               write(vname, "(a,i1.1)") "height", coord_heights(kc)
            else
               write(vname, "(a,i2.2)") "height", coord_heights(kc)
            end if
            ierr = nf90_inq_varid(ncid, vname, vid)
            call check_ncerr(ierr,"Error getting vid for height coord")
            ierr = nf90_put_var ( ncid, vid, real(coord_heights(kc)))
            call check_ncerr(ierr,"Error writing coordinate height")
         end do

#ifdef outsync
!        Sync the file so that if the program crashes for some reason 
!        there will still be useful output.
         ierr = nf90_sync ( ncid )
         call check_ncerr(ierr, "Error syncing history file")
#endif
      
      end if ! myid==0

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
      avetime = 0.0
      ! Initialisation so min/max work
      avetime_bnds(:) = (/huge(1.), -huge(1.)/)
      timecount = 0
      do ifld = 1, totflds
         if ( histinfo(ifld)%used ) then
            istart = histinfo(ifld)%ptr
            iend = istart + histinfo(ifld)%nlevels - 1
            histarray(:,:,istart:iend) = &
               initval(histinfo(ifld)%ave_type)
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
          
      do i=1,totflds
         inames(:,i) = hashkey ( histinfo(i)%name )
      end do

      do i=1,totflds

!        Find the first element in the rest of the list
         itemp(:) = inames(:,i)
         ipos = i
         do j=i+1,totflds
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
         do i=1,totflds
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
    
!-------------------------------------------------------------------
   subroutine create_ncvar(vinfo, ncid, dims)

      use mpidata_m
      type(hinfo), intent(inout) :: vinfo
      integer, intent(in) :: ncid
      type(dimarray), intent(in) :: dims

      character(len=MAX_NAMELEN) :: local_name, new_name
      character(len=80) :: cell_methods, coord_name
      integer :: ierr, vtype, vid, zdim, wdim
      
      integer(kind=2), parameter :: fill_short = NF90_FILL_SHORT

      if ( myid /= 0 ) return
      
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
         else if ( vinfo%pop4d ) then
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim, wdim, dims%t /), vid )
         else
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim, dims%t /), vid )
         end if
      else
         if ( vinfo%ave_type == hist_fixed ) then
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y /), vid )
         else
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, dims%t /), vid )
         end if
      end if
      if ( ierr /= 0 ) then
         print*, "Error creating variable ", vinfo
      end if
      call check_ncerr(ierr,"Error creating variable "// local_name)
      vinfo%vid = vid
      
      if ( cordex_compliant ) then
        ierr = nf90_def_var_deflate( ncid, vid, 1, 1, 1 ) ! shuffle=1, deflate=1, deflate_level=1
      end if

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
         vinfo%scalef = 1.015625 * (vinfo%valid_max - vinfo%valid_min) / float(vmax - vmin)
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

      ! If the variable has cell_methods set, then compose this with the 
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

      if ( len_trim(cell_methods) > 0 ) then
         ierr = nf90_put_att ( ncid, vid, "cell_methods", cell_methods )
         call check_ncerr(ierr,"Error with cell_methods attribute")
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
         ! This will cover 2m and 10m. No other likely values?
         if ( vinfo%coord_height < 10. ) then
            write(coord_name, "(a,i1.1)") "height", nint(vinfo%coord_height)
         else
            write(coord_name, "(a,i2.2)") "height", nint(vinfo%coord_height)
         end if
         ierr = nf90_put_att ( ncid, vid, "coordinates", coord_name )
         call check_ncerr(ierr,"Error with coordinates attribute")
      end if

   end subroutine create_ncvar
  
!---------------------------------------------------------------------------
   subroutine create_ncfile ( filename, nxhis, nyhis, nlev, ol, cptch, cchrt, multilev,            &
                 use_plevs, use_meters, use_depth, use_hyblevs, basetime,            &
                 coord_heights, ncid, dims, dimvars, source, extra_atts, calendar,   &
                 nsoil, zsoil, osig_found )

      use mpidata_m
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nxhis, nyhis, nlev, ol, cptch, cchrt
      logical, intent(in) :: multilev, use_plevs, use_meters, use_depth, use_hyblevs
      logical, intent(in) :: osig_found
      character(len=*), intent(in) :: basetime
      integer, dimension(:), intent(in) :: coord_heights
      integer, intent(out) :: ncid
      type(dimarray), intent(out) :: dims, dimvars
      character(len=*), intent(in), optional :: source
      type(hist_att), dimension(:), optional :: extra_atts
      character(len=*), intent(in), optional :: calendar
      integer, intent(in), optional :: nsoil ! Number of soil levels
      ! Soil depths
      real, dimension(:), intent(in), optional :: zsoil

      integer :: ierr
      integer :: i, k
      character(len=80) :: histstr
      character(len=20) :: tmpname
      integer :: vid

      if ( myid /=0 ) return

      if ( hist_debug > 0 ) then
         print*, "Creating file ", filename
      end if
#ifdef usenc3
      ierr = nf90_create(filename, nf90_64bit_offset, ncid)
#else
      ierr = nf90_create(filename, nf90_netcdf4, ncid)
#endif
      call check_ncerr ( ierr, "Error in creating history file" )
               
!     Create dimensions, lon, lat and rec
      ierr = nf90_def_dim ( ncid, "lon", nxhis, dims%x )
      call check_ncerr(ierr,"Error creating lon dimension")
      if ( hist_debug > 5 ) print*, "Created lon dimension, id",  dims%x
      ierr = nf90_def_dim ( ncid, "lat", nyhis, dims%y )
      call check_ncerr(ierr,"Error creating lat dimension")
      if ( hist_debug > 5 ) print*, "Created lon dimension, id",  dims%y
      if ( cf_compliant ) then
         ierr = nf90_def_dim ( ncid, "bnds", 2, dims%b )
         call check_ncerr(ierr,"Error creating bnds dimension")
         if ( hist_debug > 5 ) print*, "Created bnds dimension, id",  dims%b

      end if

!     Only create the lev dimension if one of the variables actually uses it
      if ( multilev ) then
         if ( use_meters ) then
            ierr = nf90_def_dim ( ncid, "alt", nlev, dims%z )
         else
            ierr = nf90_def_dim ( ncid, "lev", nlev, dims%z )
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
         if ( nsoil > 1 ) then
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
      ierr = nf90_def_dim ( ncid, "time", NF90_UNLIMITED, dims%t )
      call check_ncerr(ierr,"Error creating time dimension")
      if ( hist_debug > 5 ) print*, "Created time dimension, id",  dims%t
      
      if ( present(source) ) then
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "source", source )
         call check_ncerr(ierr)
      end if
      !if ( cf_compliant ) then
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "Conventions", "CF-1.7" )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "title", "CCAM simulation data" )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "contact", "ccam@csiro.au" )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "project", "Undefined CCAM project" )
         call check_ncerr(ierr)
      !end if
!     Make sure it's a null string, otherwise len_trim won't find the end 
!     properly.
      histstr = "" 
      call hstring(histstr)
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "history", histstr )
      if ( present(extra_atts) .and. save_ccam_parameters ) then
         do i=1,size(extra_atts)
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
         end do
      end if
      
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "modeling_realm", "atmos" )

!     Define attributes for the dimensions
!     Is there any value in keeping long_name for latitude and longitude?
!     Possibly add bounds here??
      ierr = nf90_def_var ( ncid, "lon", NF90_FLOAT, dims%x, dimvars%x)
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "long_name", "longitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "standard_name", "longitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "axis", "X" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%x, "units", "degrees_east" )
      call check_ncerr(ierr)
      if ( cf_compliant ) then
         ierr = nf90_put_att ( ncid, dimvars%x, "bounds", "lon_bnds" )
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "lon_bnds", NF90_FLOAT, (/dims%b, dims%x/), dimvars%x_b)
         call check_ncerr(ierr)
      end if

      ierr = nf90_def_var ( ncid, "lat", NF90_FLOAT, dims%y, dimvars%y )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "long_name", "latitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "standard_name", "latitude" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "axis", "Y" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, dimvars%y, "units", "degrees_north" )
      call check_ncerr(ierr)
      if ( cf_compliant ) then
         ierr = nf90_put_att ( ncid, dimvars%y, "bounds", "lat_bnds" )
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "lat_bnds", NF90_FLOAT, (/dims%b, dims%y/), dimvars%y_b)
         call check_ncerr(ierr)
      end if

      if ( multilev ) then
         if ( use_plevs ) then
            ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "pressure_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "hPa" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "down" )
            call check_ncerr(ierr)
         else if ( use_meters ) then
            ierr = nf90_def_var ( ncid, "alt", NF90_FLOAT, dims%z, dimvars%z )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "standard_name", "height" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "long_name", "vertical distance above the surface" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "units", "m" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%z, "positive", "up" )
            call check_ncerr(ierr)
         else if ( use_hyblevs ) then
            ! This is required for grads to recognise it's a vertical dimension
            ! Allowed by CF convention
            ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
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
            ierr = nf90_def_var ( ncid, "lev", NF90_FLOAT, dims%z, dimvars%z )
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
         if ( .not. present(zsoil) ) then
            print*, "Error - missing zsoil argument"
            stop
         end if
         if ( nsoil>0 ) then
            ierr = nf90_def_var ( ncid, "depth", NF90_FLOAT, dims%zsoil, dimvars%zsoil)
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%zsoil, "long_name", "depth of soil layers" )
            ierr = nf90_put_att ( ncid, dimvars%zsoil, "standard_name", "depth" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%zsoil, "units", "m" )
            call check_ncerr(ierr)
            if ( cf_compliant ) then
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
      if ( cf_compliant ) then
         ierr = nf90_put_att ( ncid, dimvars%t, "bounds", "time_bnds" )
         call check_ncerr(ierr)
         ierr = nf90_def_var ( ncid, "time_bnds", NF90_FLOAT, (/dims%b, dims%t/), dimvars%t_b)
         call check_ncerr(ierr)
      end if


      if ( multilev .and. use_hyblevs ) then
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
      do k=1,size(coord_heights)
         ! This code is repeated several times. Could do better
         if ( coord_heights(k) < 10. ) then
            write(tmpname, "(a,i1.1)") "height", coord_heights(k)
         else
            write(tmpname, "(a,i2.2)") "height", coord_heights(k)
         end if
         ierr = nf90_def_var ( ncid, tmpname, NF90_FLOAT, vid )
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "standard_name",  "height")
         ierr = nf90_put_att ( ncid, vid, "units",  "m")
         call check_ncerr(ierr)
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "axis",  "Z")
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "positive",  "up")
         call check_ncerr(ierr)
      end do

   end subroutine create_ncfile
   
!-----------------------------------------------------------------------------
   subroutine closehist
      
      use logging_m
      
      integer :: ierr
     
      call START_LOG(closehist_begin)
      
!     The hist_oave files are closed individually in writehist.
      ierr = nf90_close ( histid )
      call check_ncerr(ierr,"Error closing history file")
      
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
      
      if ( all( array == NF90_FILL_FLOAT ) ) return

!     Check that openhist has been called to allocate the array
      if ( .not. allocated ( histarray ) ) then
         print*, " History error - must call openhist before savehist "
         stop
      end if

!     Loop over history files because variable could occur in several.

      ifld = qindex_hname(name,inames(:,1:totflds),totflds)
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
         where ( abs(array) /= nf90_fill_float .and. histarray(:,1:jlat2,istart:iend) /= nf90_fill_float ) 
            histarray(:,1:jlat2,istart:iend) =  &
                 histarray(:,1:jlat2,istart:iend) + array
         elsewhere
            histarray(:,1:jlat2,istart:iend) = nf90_fill_float 
         end where    
      case ( hist_max ) 
         where ( abs(array) /= nf90_fill_float .and. histarray(:,1:jlat2,istart:iend) /= nf90_fill_float )  
            histarray(:,1:jlat2,istart:iend) =  &
                 max ( histarray(:,1:jlat2,istart:iend), array )
         elsewhere
            histarray(:,1:jlat2,istart:iend) = nf90_fill_float 
         end where 
      case ( hist_min )
         where ( abs(array) /= nf90_fill_float .and. histarray(:,1:jlat2,istart:iend) /= nf90_fill_float )     
            histarray(:,1:jlat2,istart:iend) =  &
                 min ( histarray(:,1:jlat2,istart:iend), array )
         elsewhere
            histarray(:,1:jlat2,istart:iend) = nf90_fill_float 
         end where 
      case ( hist_inst, hist_fixed ) 
         where ( abs(array) /= nf90_fill_float ) 
            histarray(:,1:jlat2,istart:iend) = array
         elsewhere
            histarray(:,1:jlat2,istart:iend) = nf90_fill_float 
         end where    
      case default
         print*, "Internal error in history, unknown ave_type", &
              histinfo(ifld)%ave_type
         stop
      end select
     
      if ( any( abs(array) /= nf90_fill_float ) ) then
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
   subroutine writehist ( istep, endofrun, year, month, interp, time, time_bnds )

      use mpidata_m
      use newmpar_m, only : cptch, cchrt
      use logging_m
      
      integer, intent(in) :: istep
      logical, intent(in), optional :: endofrun
      integer, intent(in), optional :: year, month !  For old format files
      real, intent(in), optional :: time
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
      integer ip, n
      integer, dimension(5) :: start4D, count4D
      integer, dimension(4) :: start3D, count3D
      integer, dimension(3) :: start2D, count2D
      integer :: istart, iend
      integer :: ave_type, nlev, count, ncid
      integer :: k
      real :: umin, umax, addoff, sf
!     Temporary for interpolated output
      integer, dimension(:), allocatable, save :: k_indx
      real, dimension ( nxhis, nyhis ) :: htemp
      real, dimension(:,:), allocatable, save :: hist_g
      logical :: doinc
#ifdef usempi3
      integer :: cnt,maxcnt,interp_nproc
      integer :: slab,offset
#endif
      
      if ( require_interp .and. .not. present(interp) ) then
         print*, " Error, interp argument required for writehist "
         stop
      end if

      call START_LOG(writehist_begin)

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
      if ( doinc .and. ihtype == hist_ave .and. present(time) ) then
         avetime = avetime + time
         timecount = timecount + 1
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

      if ( hist_debug > 0 ) then
         print*, "Writing history file "
      end if

      histset = histset + 1

      if ( myid == 0 ) then

         ncid = histid

         ierr = nf90_inq_varid (ncid, "time", vid )
         call check_ncerr(ierr, "Error getting time id")
         if ( present(time) ) then
            if ( ihtype == hist_ave) then
               ierr = nf90_put_var ( ncid, vid,   &
                    avetime/timecount, start=(/histset/) )
            else
               ierr = nf90_put_var ( ncid, vid, time, start=(/histset/))
            end if
         else
            ierr = nf90_put_var ( ncid, vid, &
                                  real(histset*hfreq), &
                                  start=(/histset/) )
         end if
         call check_ncerr(ierr, "Error writing time")
         if ( cf_compliant .and. present(time_bnds) ) then
            ierr = nf90_inq_varid (ncid, "time_bnds", vid )
            call check_ncerr(ierr, "Error getting time_bnds id")
            if ( ihtype == hist_ave) then
               ierr = nf90_put_var ( ncid, vid, avetime_bnds(:), start=(/1,histset/))
            else
               ierr = nf90_put_var ( ncid, vid, time_bnds, start=(/1,histset/))
            end if
            call check_ncerr(ierr, "Error writing time_bnds")
         end if

         start2D = (/ 1, 1, histset /)
         count2D = (/ nxhis, nyhis, 1 /)
         count3D = (/ nxhis, nyhis, 1, 1 /)
         count4D = (/ nxhis, nyhis, 1, 1, 1 /)
         
      end if

#ifdef usempi3
      cnt = 0
!     first pass
!     find total number of levels to store
      do ifld = 1,totflds
         if ( .not. histinfo(ifld)%used ) then
            cycle
         end if
         ave_type = histinfo(ifld)%ave_type
         nlev = histinfo(ifld)%nlevels
         count = histinfo(ifld)%count

!        Only write fixed variables in the first history set
         if ( histset > 1 .and. ave_type == hist_fixed ) then
            cycle
         end if

         istart = histinfo(ifld)%ptr
         iend = istart + nlev - 1
         if ( ave_type == hist_ave ) then    
            if ( hist_debug >= 4 ) then
               print*, "Raw history at point ", histinfo(ifld)%name,&
                 histarray(ihdb,jhdb,istart+khdb-1), count
            end if
            where ( histarray(:,:,istart:iend) /= NF90_FILL_FLOAT )
               histarray(:,:,istart:iend) =   &
                 histarray(:,:,istart:iend) / max( count, 1 )
            end where  
         end if
         if ( histinfo(ifld)%output_scale /= 0 ) then
            histarray(:,:,istart:iend) = histarray(:,:,istart:iend) * &
              histinfo(ifld)%output_scale
         end if
         if ( hist_debug >= 4 ) then
            print*, "History written at point ", histinfo(ifld)%name,&
                 histarray(ihdb,jhdb,istart+khdb-1)
         end if

!        Even multilevel variables are written one level at a time
         if ( count /= 0 ) then
           cnt = cnt + iend - istart + 1
         end if
               
      end do ! Loop over fields

      slab = ceiling(real(cnt,8)/real(nproc,8))
      maxcnt = cnt
      interp_nproc = ceiling(real(maxcnt,8)/real(max(slab,1),8))
      offset = nproc - interp_nproc
      if ( myid >= offset ) then
         allocate( hist_g(nx_g,ny_g) )
         allocate( hist_a(pil,pjl*pnpan,pnproc,1+slab*(myid-offset):slab*(myid-offset+1)) ) 
      end if
      allocate( k_indx(maxcnt) )

      cnt = 0
!     second pass
!     create the array used to index histarray
      do ifld = 1,totflds
         if ( .not. histinfo(ifld)%used ) then
            cycle
         end if
         ave_type = histinfo(ifld)%ave_type
         nlev = histinfo(ifld)%nlevels
         count = histinfo(ifld)%count

!        Only write fixed variables in the first history set
         if ( histset > 1 .and. ave_type == hist_fixed ) then
            cycle
         end if

         istart = histinfo(ifld)%ptr
         iend = istart + nlev - 1

!        Even multilevel variables are written one level at a time
         if ( count /= 0 ) then
            do k = istart,iend
               cnt = cnt + 1
               k_indx(cnt) = k
            end do   ! k loop
         end if  

      end do ! Loop over fields

!     now do the gather wrap
      if ( slab > 0 ) then
         call gather_wrap(histarray,hist_a,slab,offset,maxcnt,k_indx)
      end if
 
      cnt = 0
!     third pass
!     perform the interpolation and write the data
      do ifld = 1,totflds
         if ( .not. histinfo(ifld)%used ) then
            cycle
         end if
         ave_type = histinfo(ifld)%ave_type
         nlev = histinfo(ifld)%nlevels
         count = histinfo(ifld)%count
         vid = histinfo(ifld)%vid

!        Only write fixed variables in the first history set
         if ( histset > 1 .and. ave_type == hist_fixed ) then
            cycle
         end if

         istart = histinfo(ifld)%ptr
         iend = istart + nlev - 1

!        Even multilevel variables are written one level at a time
         do k = istart, iend

            if ( count == 0 ) then
               if ( hbytes == 2 ) then
                  htemp = NF90_FILL_SHORT
               else
                  htemp = NF90_FILL_FLOAT
               end if
        
            else

               cnt = cnt + 1
               if ( (cnt>=(1+slab*(myid-offset))).and.(cnt<=(slab*(myid-offset+1))) ) then

                  do ip = 0,pnproc-1
                     do n = 0,pnpan-1
                        hist_g(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g) = &
                           hist_a(1:pil,1+n*pjl:(n+1)*pjl,ip+1,cnt)
                     end do
                  end do
                
                  if ( present(interp) ) then
                     call interp ( hist_g, htemp, histinfo(ifld)%int_type )
                  else
                     htemp = hist_g(:,:)
                  end if
                  
               end if

            end if

            call sendrecv_wrap(htemp,cnt,slab,offset)
               
            if ( myid == 0 ) then

               if ( count /= 0 ) then
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
               end if
                   
               call START_LOG(putvar_begin)
               if ( histinfo(ifld)%pop4d ) then
                  start4D = (/ 1, 1, mod(k+1-istart-1,cptch)+1, 1+(k+1-istart-1)/cptch,  histset /)
                  ierr = nf90_put_var ( ncid, vid, htemp, start=start4D, count=count4D )
               else if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                  start3D = (/ 1, 1, k+1-istart, histset /)
                  ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
               else
                  ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
               end if
               call check_ncerr(ierr, "Error writing history variable "//histinfo(ifld)%name )
               call END_LOG(putvar_end)

            end if
                 
         end do   ! k loop
               
!        Zero ready for next set
         histarray(:,:,istart:iend) = initval(ave_type)
!        Reset the count variable
         histinfo(ifld)%count = 0

      end do ! Loop over fields

      if ( myid >= offset ) then
         deallocate(hist_a, hist_g)
      end if
      deallocate(k_indx)
#else
      if ( myid == 0 ) then
         allocate( hist_g(nx_g,ny_g) )
      else
         allocate( hist_g(0,0) )
         allocate( hist_a(0,0,0,0) )
      end if

      do ifld = 1,totflds
         if ( .not. histinfo(ifld)%used ) then
            cycle
         end if
         ave_type = histinfo(ifld)%ave_type
         nlev = histinfo(ifld)%nlevels
         count = histinfo(ifld)%count
         vid = histinfo(ifld)%vid

!        Only write fixed variables in the first history set
         if ( histset > 1 .and. ave_type == hist_fixed ) then
            cycle
         end if

         istart = histinfo(ifld)%ptr
         iend = istart + nlev - 1
         if ( ave_type == hist_ave ) then    
            if ( hist_debug >= 4 ) then
               print*, "Raw history at point ", histinfo(ifld)%name,&
                 histarray(ihdb,jhdb,istart+khdb-1), count
            end if
            where ( histarray(:,:,istart:iend) /= NF90_FILL_FLOAT )
               histarray(:,:,istart:iend) =   &
                 histarray(:,:,istart:iend) / max( count, 1 )
            end where  
         end if
         if ( histinfo(ifld)%output_scale /= 0 ) then
            histarray(:,:,istart:iend) = histarray(:,:,istart:iend) * &
              histinfo(ifld)%output_scale
         end if
         if ( hist_debug >= 4 ) then
            print*, "History written at point ", histinfo(ifld)%name,&
                 histarray(ihdb,jhdb,istart+khdb-1)
         end if
            
         if ( myid == 0 ) then
            if ( allocated(hist_a) ) then 
               if ( size(hist_a,3) /= nlev ) then
                  deallocate( hist_a )
                  allocate( hist_a(pil,pjl*pnpan,nlev,pnproc) )
               end if
            else
               allocate( hist_a(pil,pjl*pnpan,nlev,pnproc) )
            end if    
         end if    

         call gather_wrap(histarray(:,:,istart:iend),hist_a)
            
         if ( myid == 0 ) then
                
!           Even multilevel variables are written one level at a time
            do k = istart, iend

               do ip = 0,pnproc-1   
                  do n = 0,pnpan-1
                     hist_g(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g) = &
                        hist_a(1:pil,1+n*pjl:(n+1)*pjl,k-istart+1,ip+1)
                  end do
               end do
  
               if ( count == 0 ) then
                  if ( hbytes == 2 ) then
                     htemp = NF90_FILL_SHORT
                  else
                     htemp = NF90_FILL_FLOAT
                  end if
               else

                  if ( present(interp) ) then
                     call interp ( hist_g, htemp, histinfo(ifld)%int_type )
                  else
                     htemp = hist_g
                  end if

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
                  
               end if

               call START_LOG(putvar_begin)
               if ( histinfo(ifld)%pop4d ) then
                  start4D = (/ 1, 1, mod(k+1-istart-1,cptch)+1, 1+(k+1-istart-1)/cptch,  histset /)
                  ierr = nf90_put_var ( ncid, vid, htemp, start=start4D, count=count4D )
               else if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                  start3D = (/ 1, 1, k+1-istart, histset /)
                  ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
               else
                  ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
               end if
               call check_ncerr( ierr, "Error writing history variable "//histinfo(ifld)%name )
               call END_LOG(putvar_begin)
                 
            end do   ! k loop
               
         end if ! myid == 0
           
!        Zero ready for next set
         histarray(:,:,istart:iend) = initval(ave_type)
!        Reset the count variable
         histinfo(ifld)%count = 0

      end do ! Loop over fields
         
      deallocate( hist_a, hist_g )
#endif

      avetime = 0.0
      ! Initialisation so min/max work
      avetime_bnds(:) = (/huge(1.), -huge(1.)/)
      timecount = 0

!     Sync the file so that if the program crashes for some reason 
!     there will still be useful output.
#ifdef outsync
      if ( myid == 0 ) then
         ierr = nf90_sync ( ncid )
         call check_ncerr(ierr, "Error syncing history file")
      end if
#endif

      call END_LOG(writehist_end)

   end subroutine writehist

!-------------------------------------------------------------------

   function initval ( ave_type ) result (val)
      integer, intent(in) :: ave_type
      real :: val
!     Set a suitable initialisation value depending on the ave_type
      select case ( ave_type)
      case ( hist_inst, hist_fixed )
         val = NF90_FILL_FLOAT
      case ( hist_ave )    
         val = 0.0
      case (hist_max)
         val = -huge(1.0)
      case (hist_min)
         val = huge(1.0)
      case default
         print*, " History internal error: average type set incorrectly in initval"
         stop
      end select
   end function initval

!-------------------------------------------------------------------

   subroutine clearhist ()

      ! Reset the history outside the normal process

      integer :: ifld
      integer :: istart, iend
      integer :: ave_type, nlev

      if ( hist_debug > 0 ) then
         print*, "Resetting history"
      end if

      do ifld = 1,totflds
         if ( .not. histinfo(ifld)%used ) then
            cycle
         end if
         ave_type = histinfo(ifld)%ave_type
         nlev = histinfo(ifld)%nlevels
         istart = histinfo(ifld)%ptr
         iend = istart + nlev - 1

!        Zero ready for next set
         histarray(:,:,istart:iend) = initval(ave_type)
!        Reset the count variable
         histinfo(ifld)%count = 0

      end do ! Loop over fields

   end subroutine clearhist

#ifdef usempi3
   subroutine gather_wrap(histarray,hist_a,slab,offset,maxcnt,k_indx)
      use mpidata_m, only : nproc, lproc, myid, pil, pjl, pnpan, comm_world
      use logging_m
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif 
      integer, intent(in) :: slab, offset, maxcnt
      integer :: istart, iend, ip, k, n, ierr, lsize, lp, iq
      integer, dimension(maxcnt), intent(in) :: k_indx
      real, dimension(:,:,:), intent(in) :: histarray
      real, dimension(:,:,:,:), pointer, contiguous, intent(inout) :: hist_a
      !real, dimension(pil*pjl*pnpan*lproc*slab*nproc), target :: hist_a_tmp
      real, dimension(pil*pjl*pnpan*lproc*slab*nproc) :: hist_a_tmp
      !real, dimension(:,:,:,:), pointer, contiguous :: hist_a_remap, hist_a_tmp_remap
#ifdef isendrecv      
      integer :: nreq, rreq, ipr
      integer, save :: itag = 0
      integer, dimension(2*nproc) :: ireq
      integer, dimension(MPI_STATUS_SIZE,2*nproc) :: status
      real, dimension(pil,pjl*pnpan,lproc,slab,nproc) :: histarray_tmp      
#else
      real, dimension(pil,pjl*pnpan,lproc,slab) :: histarray_tmp
#endif
      
      call START_LOG(gatherwrap_begin)

#ifdef isendrecv      
      nreq = 0
      rreq = 0
      itag = itag + 1
      istart = 1 + slab*(myid-offset) 
      if ( istart > 0 ) then
         iend = istart + slab - 1
         iend = min( iend, maxcnt )
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
#endif
      
      do ip = 0,nproc-1
         istart = 1 + slab*(ip-offset) 
         if ( istart > 0 ) then
            iend = istart + slab - 1
            iend = min( iend, maxcnt )
#ifdef isendrecv            
            do k = istart,iend
               histarray_tmp(:,:,:,k-istart+1,ip+1) = reshape( histarray(:,:,k_indx(k)), (/ pil,pjl*pnpan,lproc /) )
            end do   
            call START_LOG(mpigather_begin)
            lsize = pil*pjl*pnpan*lproc*(iend-istart+1)
            nreq = nreq + 1
            call MPI_ISend( histarray_tmp(:,:,:,:,ip+1), lsize, MPI_REAL, ip, &
                 itag, comm_world, ireq(nreq), ierr )
            call END_LOG(mpigather_end)
#else
            do k = istart,iend
               histarray_tmp(:,:,:,k-istart+1) = reshape( histarray(:,:,k_indx(k)), (/ pil,pjl*pnpan,lproc /) )
            end do   
            call START_LOG(mpigather_begin)
            lsize = pil*pjl*pnpan*lproc*(iend-istart+1)
            call MPI_Gather(histarray_tmp(:,:,:,:), lsize, MPI_REAL, hist_a_tmp, lsize, MPI_REAL,  &
                            ip, comm_world, ierr)
            call MPI_Barrier(comm_world, ierr) ! avoids crashes on some systems
            call END_LOG(mpigather_end)
#endif
         end if
      end do
      
      istart = 1 + slab*(myid-offset)
      if ( istart > 0 ) then
#ifdef isendrecv
         call START_LOG(mpigather_begin)
         call MPI_Waitall( rreq, ireq, status, ierr )
         call END_LOG(mpigather_end)          
#endif
         iend = slab*(myid-offset+1)
         iend = min( iend, maxcnt )          
         !hist_a_remap(1:pil,1:pjl*pnpan*lproc,1:nproc,istart:iend) => hist_a
         !hist_a_tmp_remap(1:pil,1:pjl*pnpan*lproc,istart:iend,1:nproc) =>    &
         !    hist_a_tmp(1:pil*pjl*pnpan*lproc*(iend-istart+1)*nproc)
         do n = 1,nproc
            do k = istart,iend 
               !hist_a_remap(:,:,n,k) = hist_a_tmp_remap(:,:,k,n)
               do lp = 0,lproc-1 
                  iq = lp*pil*pjl*pnpan + (k-istart)*pil*pjl*pnpan*lproc + (n-1)*pil*pjl*pnpan*lproc*(iend-istart+1)
                  hist_a(:,:,lp+(n-1)*lproc+1,k) = reshape( hist_a_tmp(iq+1:iq+pil*pjl*pnpan), (/ pil, pjl*pnpan /) )
               end do
            end do
         end do
      end if 
      
#ifdef isendrecv
      call START_LOG(mpigather_begin)
      call MPI_Waitall( nreq, ireq, status, ierr )
      call END_LOG(mpigather_end)
      call MPI_Barrier(comm_world, ierr) ! avoids crashes on some systems
#endif
      
      call END_LOG(gatherwrap_end)
   
   end subroutine gather_wrap
   
   subroutine sendrecv_wrap(htemp,cnt,slab,offset)
      use mpidata_m, only : nproc, lproc, myid, comm_world
      use logging_m
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif  
      real, dimension(:,:), intent(inout) :: htemp
      integer, intent(in) :: cnt, slab, offset
      integer :: rrank, sizehis, ierr
   
      call START_LOG(mpisendrecv_begin)
      
      rrank = ceiling(real(cnt,8)/real(max(slab,1),8)) - 1 + offset
      if ( rrank /= 0 ) then
         if ( myid == 0 ) then
            sizehis = size(htemp, 1)*size(htemp, 2) 
            call MPI_Recv(htemp,sizehis,MPI_REAL,rrank,1,comm_world,MPI_STATUS_IGNORE,ierr)
         else if ( myid == rrank ) then
            sizehis = size(htemp, 1)*size(htemp, 2) 
            call MPI_Send(htemp,sizehis,MPI_REAL,0,1,comm_world,ierr)
         end if
      end if
      
      call END_LOG(mpisendrecv_end)
       
   end subroutine sendrecv_wrap
#else
   subroutine gather_wrap(array_in,array_out)
      use mpidata_m, only : nproc, lproc, comm_world, myid
      use logging_m
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif   
      real, dimension(:,:,:), intent(in) :: array_in
      real, dimension(:,:,:,:), intent(out) :: array_out
      real, dimension(size(array_out,1),size(array_out,2),lproc,size(array_in,3),nproc) :: array_temp
      integer :: lsize, ierr, np, lp, k
      
      call START_LOG(gatherwrap_begin)
      lsize = size(array_in)
      call START_LOG(mpigather_begin)
      call MPI_Gather(array_in,lsize,MPI_REAL,array_temp,lsize,MPI_REAL,0,comm_world,ierr)
      call END_LOG(mpigather_end)
      if ( myid == 0 ) then
         do np = 0,nproc-1
            do k = 1,size(array_in,3)
               do lp = 0,lproc-1
                  array_out(:,:,k,lp+np*lproc+1) = array_temp(:,:,lp+1,k,np+1)
               end do
            end do
         end do
      end if
      !call MPI_Barrier(comm_world,ierr) ! avoids crashes on some systems
      call END_LOG(gatherwrap_end)
      
   end subroutine gather_wrap
#endif


!-------------------------------------------------------------------
   function needfld ( name ) result (needed)
      character(len=*), intent(in) :: name
      logical :: needed
      integer :: ifld

!     Check if name is in the list for any of the history files
      needed = .false.
      ifld = qindex_hname ( name, inames(:,1:totflds), totflds)
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

!      key = hashkey ( name )
!     hashkey inlined by hand      
      key(:) = 0
!     This makes numerical order the same as alphabetical order
      do k = 1,MAX_KEYINDEX
         fac = 1 
         do i=min(k*MAX_KEYLEN,len_trim(name)),(k-1)*MAX_KEYLEN+1,-1
            ! Netcdf allowed characters for variables names are in 
            ! range 48 to 122 (alphanumeric and _)
            key(k) = key(k) + fac*(ichar(name(i:i))-48)
            fac = fac * 75
         end do   
      end do

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
   function qindex_hname ( name, table, nflds ) result(ifld)
      character(len=*), intent(in)    :: name
      integer, intent(in) :: nflds
      integer ( bigint ), intent(in), dimension(MAX_KEYINDEX,nflds) :: table
      integer :: ifld
!     Both of these are initially set to one because this is a safe value.
!     It doesn't matter that this points to the wrong thing because it's
!     always checked.
!     Index of the previously found name
      integer, save :: prev = 1
!     Index of the successor of a given index
      integer, save, dimension(nfmax) :: successor = 1

!     First check the successor of the previous name. This will generally 
!     work because things are done in a fixed order. It will only break on
!     the first latitude of the radiation step.

      ifld = successor(prev)
      if ( histinfo(ifld)%name /= name ) then
         ifld = bindex_hname ( name, table, nflds )
         if ( ifld == 0 ) then
            ! Name is not known, return without updating internal fields
            return
         end if
!        Redefine the successor.
         successor(prev) = ifld
      end if

!     Set this for next time
      prev = ifld
         

   end function qindex_hname
!-------------------------------------------------------------------
   function hashkey ( name ) result (key)
!  Generate an integer from a string
      character(len=*), intent(in) :: name
      integer (bigint), dimension(MAX_KEYINDEX) :: key
      integer :: i, k
      integer (bigint) :: fac

      key(:) = 0
      do k = 1,MAX_KEYINDEX
         fac = 1 
         ! This makes numerical order the same as alphabetical order
         do i=min(k*MAX_KEYLEN,len_trim(name)),(k-1)*MAX_KEYLEN+1,-1
            ! Netcdf allowed characters for variables names are in 
            ! range 48 to 122 (alphanumeric and _)
            key(k) = key(k) + fac*(ichar(name(i:i))-48)
            fac = fac * 75
         end do   
      end do
      if (hist_debug > 2 ) then
         print*, "  HASHKEY ", name, key
      end if
   end function hashkey
   
end module history
