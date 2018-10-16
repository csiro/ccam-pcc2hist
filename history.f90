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
    
! This module provides a new scheme for saving model history data.
! The idea is to hide the arrays used for the accumulation. All 
! manipulation can be done giving just the name of the variable.
! The idea is based on the NCAR CCM3 but it's a lot more convenient in
! Fortran 90.

! Also useful in diagnostic programs, anything that has to write a time 
! series of lat/lon data to a netcdf file.

! GrADS limitations. Single vertical dimension.

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

! Each variable has average_type. They should also have interp_type as 
! an attribute, at least where it differs from the file default.
! The file default should be a global attribute. Is this appropriate for the
! cell_methods attribute?

! Have an option in the interpolation to also apply a land-sea mask?

! The routines to save full global field at once would be better if they
! had an option to write instantaneously without requiring the histarray.
! This would be useful for 3D winds in the model which are available globally
! from the sL arrays.

! Individual variables could have double row option via addfld if that was
! useful. Could be if anything was saved from the dynamics.

! At the moment it will reuse existing monthly mean files. This could be
! incorrect if hbytes changes? Does it use the value in the file and does
! the writing check for this properly?

! Need to be able to cope with the monthly mean daily maximum temperature
! and the extreme monthly temperature. Max and min need both an averaging 
! period and a reset period. This might not require any special handling 
! if the daily max is calculated in the rest of the model and savefld called
! just once per day.

! One disadvantage of this scheme is that extra arrays are required to
! hold the instantaneous winds. In the old scheme these are just taken 
! from the standard model arrays required for the SLT.

! Should values outside the valid range be trucated to the range or set
! missing. Should there be separate scale ranges and valid ranges?

! groice is added to twice per step, in seaice and surfupl. Add a check on
! whether a variable has been updated already this step and don't increment
! count.
! With the ice routines need to be careful of the jlat==1 test for updating
! count. There may not be ice there. In this case savefld will only be called
! at some latitudes. If ice retreats over a month might it be called at 
! a different number of latitudes? This has to be avoided (must be avoided 
! at the moment?? )

! Add an accumulated type for these, groice and redice???

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
              create_ncfile, oldwrite, savehist2D,                      &
              savehist3D, hashkey, qindex_hname, savehist_work,         &
              gsavehist2D, gsavehist3D
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

!  Flag for use of double row variables in savehist calls.
   logical, private, save :: double = .false.

!  Maximum length of names of history variables
   integer, parameter :: MAX_NAMELEN = 20

!  Maximum length of string used for key generation. Names must be unique
!  within this length.
   integer, parameter :: MAX_KEYLEN = 20

!  Maximum number of history files
   integer, parameter :: MAX_HFILES = 1
!  Actual number of history files used
!  Set huge value as a trap if openhist not called.
   integer, save :: nhfiles=-huge(1)

!  Maximum number of fields allowed
!  Note that this maximum is only used for dimensioning header arrays.
!  The actual data arrays are allocatable, so it doesn't use much space.
!  It can't be made truly dynamic. The best that could be done is to 
!  re-allocate a larger array if required.
   integer, parameter :: nfmax = 1000

!  Initialise hnames so that set_hnames can find the end of the list to append.
   character(len=MAX_NAMELEN), dimension(nfmax,MAX_HFILES), save ::  &
         hnames = "", xnames = ""

!  Integer array of keys corresponding to hnames to make searching faster.
   integer (bigint), dimension(nfmax), save :: inames

!  Frequency of writing history files. This may be in either steps or minutes
!  depending on the argument of the calling routine. 
   integer, dimension(MAX_HFILES), save :: hfreq = 0

!  Number of bytes to use for history variables
   integer, dimension(MAX_HFILES), save :: hbytes = 4

!  Default averaging type for file
   integer, dimension(MAX_HFILES), save :: ihtype

!  Average time appropriate for time averaged files
   real, dimension(MAX_HFILES), save :: avetime
   real, dimension(2,MAX_HFILES), save :: avetime_bnds
   integer, dimension(MAX_HFILES), save :: timecount

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
!     Scale factor for output, e.g. to convert rain from mm/step to mm/day.
      real               :: output_scale 
      logical, dimension(MAX_HFILES) :: used
      integer, dimension(MAX_HFILES) :: vid  ! netCDF variable identifier
      integer, dimension(MAX_HFILES) :: ptr  ! Pointer to position in histarray
      integer, dimension(MAX_HFILES) :: count      ! Used for scaling
      integer, dimension(MAX_HFILES) :: ave_type   ! Type of averaging
      real, dimension(MAX_HFILES)    :: addoff
      real, dimension(MAX_HFILES)    :: scalef
      integer                        :: int_type   ! Type of interpolation
!     Flag to force variable to be handled as multi-level even if nlevels=1.
      logical                        :: multilev
!     For multilevel land surface variables
      logical                        :: soil
!     For multilevel ocean variables
      logical                        :: water
      ! For CF coordinate attribute
      real                           :: coord_height
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
   integer, dimension(MAX_HFILES) :: histset

!   netCDF file ID of history file
   integer, dimension(MAX_HFILES) :: histid      

!  Total number of fields defined (not necessarily used).
   integer, save :: totflds = 0

   interface savehist
      module procedure savehist2D, savehist3D, gsavehist2D, gsavehist3D
   end interface
    
!  Define parameters to identify the processing of history data: ave, max, min
!  or instantaneous value. hist_fixed is for fields like the surface height
!  which don't change but are still useful to have in a history file.

   integer, parameter :: hist_ave=1, hist_max=2, hist_min=3, hist_inst=4, &
                         hist_fixed=5

!  Valid range for 16 bit values
   integer, parameter :: vmin=-32500, vmax=32500

   real :: missing_value = NF90_FILL_FLOAT
   real :: missing_value_cordex = 1.00000002004e+20
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
   logical, public :: no_underscore = .false.
   
! Save CCAM parameters
   logical, public :: save_ccam_parameters = .true.

!  MPI working arrays
   integer, private, save :: nx_g, ny_g
#ifdef usempi3
   real, dimension(:,:,:,:), pointer, contiguous :: hist_a
#else
   real, dimension(:,:,:,:), allocatable, save, private :: hist_a
#endif

contains

!-------------------------------------------------------------------

   subroutine hstring(hstr)

      ! Create a string containing the username, machine name and date to 
      ! identify the creation of the history file
   
      character(len=*), intent(out) :: hstr
      character(len=100) :: logname, hostname
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
      character(len=MAX_NAMELEN), dimension(MAX_HFILES) :: htype
      integer, dimension(MAX_HFILES) :: freq
      integer, dimension(MAX_HFILES) :: bytes
      character(len=MAX_NAMELEN), dimension(nfmax,MAX_HFILES) ::  &
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

      do j=1,size(htype)
         select case (htype(j))
         case ("ave")
            ihtype(j) = hist_ave
         !case ("oave")
         !   ihtype(j) = hist_oave
         case ("max")
            ihtype(j) = hist_max
         case ("min")
            ihtype(j) = hist_min
         case ("inst")
            ihtype(j) = hist_inst
         case ("fixed")
            ihtype(j) = hist_fixed
         case ("")
            ihtype(j) = hist_ave
         case default
            print*, " Error: average type set incorrectly in namelist ", &
                    trim(htype(j))
            stop
         end select
      end do

      if ( hist_debug > 0 ) then
         do j=1,MAX_HFILES
            do i=1,nfmax
               if ( len_trim(hnames(i,j)) /= 0 ) then
                  write(*,"(2i3,1x,a)") i, j, hnames(i,j)
               else
                  exit
               end if
            end do
         end do
      end if

   end subroutine histnamelist

   subroutine set_htype ( htype, ifile )
!     Routine to allow setting htype directly without namelist.
      character(len=*), intent(in) :: htype
      integer, intent(in), optional :: ifile
      integer :: j

      if ( present(ifile) ) then 
         j = ifile
      else
         j = 1
      end if
      select case (htype)
      case ("ave")
         ihtype(j) = hist_ave
      !case ("oave")
      !   ihtype(j) = hist_oave
      case ("max")
         ihtype(j) = hist_max
      case ("min")
         ihtype(j) = hist_min
      case ("inst")
         ihtype(j) = hist_inst
      case ("fixed")
         ihtype(j) = hist_fixed
      case ("")
         ihtype(j) = hist_ave
      case default
         print*, " Error: average type set incorrectly in set_htype ", &
              trim(htype)
         stop
      end select

   end subroutine set_htype
   

   subroutine set_hnames ( name, ifile )
!     Routine to allow adding a variable to hnames list directly without 
!     namelist.
      character(len=*), intent(in) :: name
      integer, intent(in), optional :: ifile
      integer :: i, j

      if ( present(ifile) ) then 
         j = ifile
      else
         j = 1
      end if

!     Append this to list of hnames by searching for the first empty entry.
!     This isn't very efficient but this routine is only used at startup.
      do i=1,nfmax
         if ( len_trim((hnames(i,j))) == 0 ) then
            hnames(i,j) = name
            exit
         end if
      end do
      if ( i > nfmax ) then
         print*, " Error, nfmax exceeded in set_hnames "
         stop
      end if

   end subroutine set_hnames

   subroutine set_hbytes ( nbyte, ifile )
!     Routine to allow setting hbytes directly without namelist.
      integer, intent(in) :: nbyte
      integer, intent(in), optional :: ifile
      integer :: j

      if ( present(ifile) ) then 
         j = ifile
      else
         j = 1
      end if
      hbytes(j) = nbyte
   end subroutine set_hbytes

   subroutine set_hfreq ( freq, ifile )
!     Routine to allow setting hfreq directly without namelist.
      integer, intent(in) :: freq
      integer, intent(in), optional :: ifile
      integer :: j

      if ( present(ifile) ) then 
         j = ifile
      else
         j = 1
      end if
      hfreq(j) = freq
   end subroutine set_hfreq

!---------------------------------------------------------------------------
   subroutine addfld(name, long_name, units, valid_min, valid_max,    &
                     nlevels, amip_name, ave_type, std, output_scale, &
                     int_type, multilev, std_name, soil, water,       &
                     coord_height, cell_methods, ran_type )
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
      real, intent(in), optional :: coord_height
      character(len=*), intent(in), optional :: cell_methods
      logical, intent(in), optional   :: ran_type

!     Local variables corresponding to the optional arguments
      integer :: atype
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
      histinfo(totflds)%output_scale = scale
      histinfo(totflds)%used(:)      = .FALSE.   ! Value for now
      histinfo(totflds)%vid(:)       = 0
      histinfo(totflds)%ptr(:)       = 0
      histinfo(totflds)%count(:)     = 0
      histinfo(totflds)%ave_type(:)  = atype     ! Default if not over-ridden
      histinfo(totflds)%addoff(:)    = 0.0
      histinfo(totflds)%scalef(:)    = 0.0
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
   subroutine openhist ( nx, ny, nl, sig, ol, gosig, suffix, hlon, hlat, basetime, &
                         doublerow, year, nxout, nyout, source, histfilename,      &
                         pressure, height, depth, extra_atts, hybrid_levels, anf,  &
                         bnf, p0, calendar, nsoil, zsoil )
!
!     Create netCDF history files using information in histinfo array.
!
      use mpidata_m
      use logging_m

      integer, intent(in) :: nx, ny, nl, ol
      real, intent(in), dimension(:) :: sig
      real, intent(in), dimension(:) :: gosig
      character(len=*), intent(in)   :: suffix  ! Filename suffix
!     Longitudes and latitudes of the output history
      real, intent(in), dimension(:) :: hlon
      real, intent(in), dimension(:) :: hlat
      character(len=*), intent(in) :: basetime
      logical, intent(in), optional  :: doublerow 
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
      integer :: ncid, vid, ifile, ivar, istart, iend, ilev
      character(len=80) :: longname, units
      character(len=MAX_NAMELEN) :: vname
      logical :: used, multilev, use_plevs, use_hyblevs, use_meters, use_depth
      integer, dimension(totflds) :: coord_heights
      integer :: kc, ncoords, k, pkl
      logical :: soil_used, water_used
      real :: dx, dy
      
      integer :: i, upos
      integer, parameter :: n_underscore_names = 120
      character(len=MAX_NAMELEN) :: new_name, local_name
      character(len=MAX_NAMELEN), dimension(n_underscore_names) :: underscore_names
      
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

      if ( present ( doublerow ) ) then
         double = doublerow
      end if

!     Flag for double row physics (all savehist calls use double row 
!     variables).
      if ( hist_debug > 0 .and. double ) then
         print*, "Using double row option in history"
      end if

      if ( hist_debug > 2 ) then
         do ivar=1,totflds
            print*, ivar, histinfo(ivar)%name, histinfo(ivar)%nlevels
         end do
      end if

      call sortlist ( histinfo(1:totflds), inames(1:totflds) )

      if ( hist_debug > 2 ) then
         print*, " Names after sorting "
         do ivar=1,totflds
            print*, ivar, histinfo(ivar)%name, histinfo(ivar)%nlevels
         end do
      end if
      
      ! define names that have an underscore
      underscore_names(:) = ""
      underscore_names(1)   = "aero_vis"
      underscore_names(2)   = "alb_ave"
      underscore_names(3)   = "anth_ave"
      underscore_names(4)   = "bc_vis"
      underscore_names(5)   = "bcb_ave"
      underscore_names(6)   = "bcdd_ave"
      underscore_names(7)   = "bce_ave"
      underscore_names(8)   = "bci_s"
      underscore_names(9)   = "bcwd_ave"
      underscore_names(10)  = "bco_s"
      underscore_names(11)  = "cape_ave"
      underscore_names(12)  = "cape_max"
      underscore_names(13)  = "cbas_ave"
      underscore_names(14)  = "climate_agdd5"
      underscore_names(15)  = "climate_alpha20"
      underscore_names(16)  = "climate_biome"
      underscore_names(17)  = "climate_gmd"
      underscore_names(18)  = "climate_dmoist_max20"
      underscore_names(19)  = "climate_dmoist_min20"
      underscore_names(20)  = "climate_ivegt"
      underscore_names(21)  = "climate_max20"
      underscore_names(22)  = "climate_min20"
      underscore_names(23)  = "cnbp_ave"
      underscore_names(24)  = "cnpp_ave"
      underscore_names(25)  = "convh_ave"
      underscore_names(26)  = "cos_zen"
      underscore_names(27)  = "ctop_ave"
      underscore_names(28)  = "dew_ave"
      underscore_names(29)  = "dms_s"
      underscore_names(30)  = "dmsb_ave"
      underscore_names(31)  = "dmse_ave"
      underscore_names(32)  = "dmsso2_ave"
      underscore_names(33)  = "dust1_s"
      underscore_names(34)  = "dust1b_ave"
      underscore_names(35)  = "dust1dd_ave"
      underscore_names(36)  = "dust1e_ave"
      underscore_names(37)  = "dust1wd_ave"
      underscore_names(38)  = "dust2_s"
      underscore_names(39)  = "dust2b_ave"
      underscore_names(40)  = "dust2dd_ave"
      underscore_names(41)  = "dust2e_ave"
      underscore_names(42)  = "dust2wd_ave"
      underscore_names(43)  = "dust3_s"
      underscore_names(44)  = "dust3b_ave"
      underscore_names(45)  = "dust3dd_ave"
      underscore_names(46)  = "dust3e_ave"
      underscore_names(47)  = "dust3dd_ave"
      underscore_names(48)  = "dust3wd_ave"
      underscore_names(49)  = "dust4_s"
      underscore_names(50)  = "dust4b_ave"
      underscore_names(51)  = "dust4dd_ave"
      underscore_names(52)  = "dust4e_ave"
      underscore_names(53)  = "dust4wd_ave"
      underscore_names(54)  = "eg_ave"
      underscore_names(55)  = "epan_ave"
      underscore_names(56)  = "epot_ave"
      underscore_names(57)  = "fbeam_ave"
      underscore_names(58)  = "fg_ave"
      underscore_names(59)  = "fnee_ave"
      underscore_names(60)  = "fpn_ave"
      underscore_names(61)  = "frday_ave"
      underscore_names(62)  = "frp_ave"
      underscore_names(63)  = "frpr_ave"
      underscore_names(64)  = "frpw_ave"
      underscore_names(65)  = "frs_ave"
      underscore_names(66)  = "ga_ave"
      underscore_names(67)  = "iwp_ave"
      underscore_names(68)  = "ldust_vis"
      underscore_names(69)  = "lwp_ave"
      underscore_names(70)  = "mixd_ave"
      underscore_names(71)  = "oc_vis"
      underscore_names(72)  = "ocb_ave"
      underscore_names(73)  = "ocdd_ave"
      underscore_names(74)  = "oce_ave"
      underscore_names(75)  = "oci_s"
      underscore_names(76)  = "oco_s"
      underscore_names(77)  = "ocwd_ave"
      underscore_names(78)  = "pmsl_ave"
      underscore_names(79)  = "rgc_ave"
      underscore_names(80)  = "rgdn_ave"
      underscore_names(81)  = "rgn_ave"
      underscore_names(82)  = "rnet_ave"
      underscore_names(83)  = "rtc_ave"
      underscore_names(84)  = "rtu_ave"
      underscore_names(85)  = "sdust_vis"
      underscore_names(86)  = "sgc_ave"
      underscore_names(87)  = "sgdn_ave"
      underscore_names(88)  = "sgn_ave"
      underscore_names(89)  = "sint_ave"
      underscore_names(90)  = "so2_s"
      underscore_names(91)  = "so2dd_ave"
      underscore_names(92)  = "so2b_ave"
      underscore_names(93)  = "so2e_ave"
      underscore_names(94)  = "so2so4_ave"
      underscore_names(95)  = "so2wd_ave"
      underscore_names(96)  = "so4_s"
      underscore_names(97)  = "so4_vis"
      underscore_names(98)  = "so4b_ave"
      underscore_names(99)  = "so4dd_ave"
      underscore_names(100) = "so4e_ave"
      underscore_names(101) = "so4wd_ave"
      underscore_names(102) = "soc_ave"
      underscore_names(103) = "sot_ave"
      underscore_names(104) = "ssalt_vis"
      underscore_names(105) = "strat_nt"
      underscore_names(106) = "tscrn_ave"
      underscore_names(107) = "tsu_ave"
      underscore_names(108) = "urbantas_ave"
      underscore_names(109) = "wb1_ave"
      underscore_names(110) = "wb2_ave"
      underscore_names(111) = "wb3_ave"
      underscore_names(112) = "wb4_ave"
      underscore_names(113) = "wb5_ave"
      underscore_names(114) = "wb6_ave"
      underscore_names(115) = "wbice1_ave"
      underscore_names(116) = "wbice2_ave"
      underscore_names(117) = "wbice3_ave"
      underscore_names(118) = "wbice4_ave"
      underscore_names(119) = "wbice5_ave"
      underscore_names(120) = "wbice6_ave"

!     Save this as a module variable
      filesuffix = suffix

      nhfiles = 0
      do ifile=1,MAX_HFILES

         do ivar=1,totflds
!           Go until the end of the variable list
            if ( len_trim(hnames(ivar,ifile)) == 0 ) then
               exit
            end if
            if ( hist_debug > 2 ) then
               print*, 'Openhist', ifile,  ivar, adjustl(hnames(ivar,ifile))
            end if
            
            ! add underscore to variable names so they match
            if ( no_underscore ) then
               do i = 1,n_underscore_names
                  local_name = ""
                  local_name = underscore_names(i) 
                  upos = index(local_name,"_")
                  do while ( upos > 0 )  
                     new_name = "" 
                     new_name = local_name(1:upos-1)//local_name(upos+1:)
                     local_name = new_name
                     upos = index(local_name,"_")
                  end do
                  if ( hnames(ivar,ifile) == local_name ) then
                     hnames(ivar,ifile) = underscore_names(i)
                  end if
               end do
            end if   
            
            if ( hnames(ivar,ifile) == "all" ) then
               histinfo(1:totflds)%used(ifile) = .true.
               exit ! All defined so no point checking further
                    ! Except when xnames is used.
            else if ( hnames(ivar,ifile) == "std" ) then
!              Use or here so as not to overwrite variables that have been
!              set explicitly.
               histinfo(1:totflds)%used(ifile) = &
                histinfo(1:totflds)%used(ifile) .or. histinfo(1:totflds)%std
            else if ( hnames(ivar,ifile) == "ran" ) then
!              Include RAN requested variables
               histinfo(1:totflds)%used(ifile) = &
                histinfo(1:totflds)%used(ifile) .or. histinfo(1:totflds)%ran
            else
!              Find the name in histinfo to set the used flag.
               ifld = bindex_hname ( hnames(ivar,ifile), &
                                     inames(1:totflds), totflds )
               if ( ifld == 0 ) then
                  print*, "Error - history variable ", hnames(ivar,ifile),  &
                          " is not known. "
                  stop
               end if
               histinfo(ifld)%used(ifile) = .true.
               if ( hist_debug > 4 ) then
                  print*,  " Setting used for ", ifld, ivar, ifile, &
                       hnames(ivar,ifile), histinfo(ifld)%name
               end if
            end if

         end do
         
!        Apply the exclusion list
         do ivar = 1,totflds
!           Go until the end of the variable list
            if ( len_trim(xnames(ivar,ifile)) == 0 ) then
               exit
            end if
            if ( hist_debug > 2 ) then
               print*, 'Openhist: exclude', ifile,  ivar, adjustl(xnames(ivar,ifile))
            end if
            
!           Find the name in histinfo to set the used flag.
            ifld = bindex_hname ( xnames(ivar,ifile), inames(1:totflds), totflds )
            if ( ifld == 0 ) then
               print*, "Error - excluded history variable ", xnames(ivar,ifile)," is not known. "
               stop
            end if
            histinfo(ifld)%used(ifile) = .false.
            if ( hist_debug > 4 ) then
               print*,  " Overriding used for ", ifld, ivar, ifile, xnames(ivar,ifile), histinfo(ifld)%name
            end if

         end do

!        Exit on the first file with no variables set.
         if ( .not. any(histinfo(1:totflds)%used(ifile)) ) then
            exit
         end if

         nhfiles = nhfiles + 1

!        Check the averaging type
!        For variables for which it isn't set use the file default
         where ( histinfo(1:totflds)%ave_type(ifile) == 0 )
            histinfo(1:totflds)%ave_type(ifile) = ihtype(ifile)
         endwhere

      end do

      if ( nhfiles == 0 ) then
         return
      end if

      if ( nhfiles > 1 .and. present(histfilename) ) then
         print*, "Warning, more than one history file used so requested name"
         print*, "is ignored"
      end if

      if ( myid == 0 ) then

!!        First file may be old average format. If this is the case it has
!!        to be handled differently
!         if ( ihtype(1) == hist_oave ) then
!
!            ifile = 1
!            if ( present(histfilename) ) then
!               ! Use histfile as a path in this case
!               oldprefix = trim(histfilename) // "/s"
!            else
!               oldprefix = "s"
!            end if
!
!            do ifld = 1, totflds
!
!               if ( .not. histinfo(ifld)%used(ifile) ) then
!                  cycle
!               end if
!
!!              lookup again to get long name and units
!               vname = histinfo(ifld)%name
!               longname = histinfo(ifld)%long_name
!               units = histinfo(ifld)%units
!
!               do ilev=1,histinfo(ifld)%nlevels
!                  if ( histinfo(ifld)%nlevels == 1 ) then
!                     vname = histinfo(ifld)%name
!                  else
!                     write(vname,"(a,i2.2)") trim(histinfo(ifld)%name), ilev
!                  end if
!                  write(filename,"(a,a,a,a)" ) trim(oldprefix), trim(vname), &
!                                               trim(suffix), ".nc"
!!                 Check if this file exists. If it does there's no more to
!!                 do, if not go on to create it.
!                  inquire(file=filename, exist=used)
!                  if ( used ) then
!                     if ( hist_debug > 0 ) then
!                        print*, "Using existing file ", filename
!                     end if
!                     cycle
!                  else
!
!                     if ( .not. present(year) ) then
!                        print*, " Year argument to openhist is required for old format files"
!                        stop
!                     end if
!                     call create_oldfile ( filename, nxhis, nyhis, hbytes(ifile),&
!                                           vname, longname, units, &
!                                           histinfo(ifld)%valid_min, &
!                                           histinfo(ifld)%valid_max, &
!                                           hlat, hlon, year )
!
!                  end if
!
!               end do
!            end do
!            istart = 2
!         else
            istart = 1
!         end if ! ihtype(1) == hist_oave

!        The rest of the history files are much simpler with multilevel variables
         do ifile = istart, nhfiles

            soil_used = .false.
            water_used = .false.
            if ( present(histfilename) .and. nhfiles == 1 ) then
               filename = histfilename
            else
               write(filename,"(a,i1,a,a)" ) "hist", ifile, trim(suffix), ".nc"
            end if
!           Is it possible to do a masked maxval?
            multilev = .false.
            ncoords = 0
            do ifld=1,totflds
               if ( hist_debug > 4 ) then
                  print*, "Checking variable properties", ifld, &
                       histinfo(ifld)%name, histinfo(ifld)%used(ifile), &
                       histinfo(ifld)%nlevels, histinfo(ifld)%soil
               end if
               if ( .not. histinfo(ifld)%used(ifile) ) cycle

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
            !if ( soil_used ) then
            !   ! Better to define a new local nsoil variable?
            !   call create_ncfile ( filename, nxhis, nyhis, size(sig), size(gosig), multilev,      &
            !        use_plevs, use_meters, use_depth, use_hyblevs, basetime,                       &
            !        coord_heights(1:ncoords), ncid, dims, dimvars, source, extra_atts, calendar,   &
            !        nsoil, zsoil )
            !else
               call create_ncfile ( filename, nxhis, nyhis, size(sig), size(gosig), multilev,      &
                    use_plevs, use_meters, use_depth, use_hyblevs, basetime,                       &
                    coord_heights(1:ncoords), ncid, dims, dimvars, source, extra_atts, calendar,   &
                    nsoil, zsoil )
            !end if
            histid(ifile) = ncid

            do ifld=1,totflds
               if ( histinfo(ifld)%used(ifile) ) then
                  call create_ncvar(histinfo(ifld), ncid, ifile, dims)
               end if
            end do

!           Leave define mode
         
            ierr = nf90_enddef ( ncid )
            call check_ncerr(ierr, "Error from enddef")

!           Turn off the data filling to save time.
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
            !if ( soil_used .and. present(zsoil) ) then
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
!           Sync the file so that if the program crashes for some reason 
!           there will still be useful output.
            ierr = nf90_sync ( ncid )
            call check_ncerr(ierr, "Error syncing history file")
#endif

         end do  ! Loop over files
      
      end if ! myid==0

!     Allocate the array to hold all the history data
!     Calculate the size by summing the number of fields of each variable

      hsize = 0
      do ifile = 1,nhfiles
         do ifld = 1, totflds
            if ( histinfo(ifld)%used(ifile) ) then
               histinfo(ifld)%ptr(ifile) = hsize + 1
               hsize = hsize + histinfo(ifld)%nlevels 
               if ( hist_debug > 0 ) then
                  print*, ' Fld ', ifld, 'Ptr', histinfo(ifld)%ptr(ifile)
               end if
            end if
         end do
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
      do ifile = 1,nhfiles
         avetime(ifile) = 0.0
         ! Initialisation so min/max work
         avetime_bnds(:,ifile) = (/huge(1.), -huge(1.)/)
         timecount(ifile) = 0
         do ifld = 1, totflds
            if ( histinfo(ifld)%used(ifile) ) then
               istart = histinfo(ifld)%ptr(ifile)
               iend = istart + histinfo(ifld)%nlevels - 1
               histarray(:,:,istart:iend) = &
                  initval(histinfo(ifld)%ave_type(ifile))
            end if
         end do
      end do
      
      call END_LOG(openhist_end)

   end subroutine openhist

   subroutine sortlist ( histinfo, inames )

!     Simple insertion sort of the list of names.
!     Assumes all lower case. This should be enforced or checked somewhere

      type(hinfo), dimension(:), intent(inout) :: histinfo
      integer (bigint), dimension(:), intent(out) :: inames
      integer :: i, ipos, j
      type(hinfo) :: temp
      integer (bigint) :: itemp

      do i=1,size(histinfo)
         inames(i) = hashkey ( histinfo(i)%name )
      end do

      do i=1,size(histinfo)

!        Find the first element in the rest of the list
         itemp = inames(i)
         ipos = i
         do j=i+1,size(histinfo)
            if ( inames(j) < itemp ) then
               itemp = inames(j)
               ipos = j
            end if
         end do

!        Move the smallest value to position i
         inames(ipos) = inames(i)
         inames(i) = itemp
!        Swap histinfo elements so they keep the same order
         temp = histinfo(ipos)
         histinfo(ipos) = histinfo(i)
         histinfo(i) = temp

      end do

      if ( hist_debug > 1 ) then
         print*, "Sorted NAMES "
         do i=1,size(histinfo)
            print*, histinfo(i)%name, inames(i)
         end do
      end if

    end subroutine sortlist
    
!-------------------------------------------------------------------
   subroutine create_ncvar(vinfo, ncid, ifile, dims)

      use mpidata_m
      type(hinfo), intent(inout) :: vinfo
      integer, intent(in) :: ncid
      integer, intent(in) :: ifile
      type(dimarray), intent(in) :: dims

      character(len=MAX_NAMELEN) :: local_name, new_name
      character(len=80) :: cell_methods, coord_name
      integer :: ierr, vtype, vid, zdim
      integer :: upos
      
      integer(kind=2), parameter :: fill_short = NF90_FILL_SHORT

      if ( myid /=0 ) return
      
      local_name = vinfo%name
      if ( no_underscore ) then
         upos = index(local_name,"_")
         do while ( upos > 0 )  
            new_name = "" 
            new_name = local_name(1:upos-1)//local_name(upos+1:)
            local_name = new_name
            upos = index(local_name,"_")
         end do   
      end if    

      select case ( hbytes(ifile) )
      case ( 2 )
         vtype = NF90_INT2
      case ( 4 ) 
         vtype = NF90_FLOAT
      case ( 8 )
         vtype = NF90_DOUBLE
      case default
         print*, " Error, unsupported value for hbytes ", hbytes(ifile)
         stop
      end select
         
      if ( hist_debug > 4 ) then
         print*, "Creating variable ", local_name, vinfo%ave_type(ifile), vinfo%nlevels, vinfo%soil
         print*, "DIMID", dims%x, dims%y, dims%z, dims%t, dims%zsoil
      end if
      if ( vinfo%nlevels > 1 .or. vinfo%multilev ) then
         if (vinfo%soil) then
            zdim = dims%zsoil
         else if ( vinfo%water ) then
            zdim = dims%oz
         else
            zdim = dims%z
         end if
         if ( vinfo%ave_type(ifile) == hist_fixed ) then
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim /), vid )
         else
            ierr = nf90_def_var ( ncid, local_name, vtype, &
                                (/ dims%x, dims%y, zdim, dims%t /), vid )
         end if
      else
         if ( vinfo%ave_type(ifile) == hist_fixed ) then
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
      vinfo%vid(ifile) = vid
      
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

      if ( hbytes(ifile) == 2 ) then
         ! 32500/3200 = 1.015625, so with this scaling, end values should
         ! be exactly representable.
         vinfo%scalef(ifile) = 1.015625 * (vinfo%valid_max - vinfo%valid_min) / float(vmax - vmin)
         vinfo%addoff(ifile) = 0.5*(vinfo%valid_min+vinfo%valid_max)
         ierr  = nf90_put_att ( ncid, vid, "add_offset", vinfo%addoff(ifile) )
         call check_ncerr(ierr)
         ierr  = nf90_put_att ( ncid, vid, "scale_factor", vinfo%scalef(ifile))
         call check_ncerr(ierr)
         ierr = nf90_put_att ( ncid, vid, "valid_min", vmin )
         call check_ncerr(ierr,"Error setting valid min attribute")
         ierr = nf90_put_att ( ncid, vid, "valid_max", vmax )
         call check_ncerr(ierr,"Error setting valid max attribute")
      end if

      ! If the variable has cell_methods set, then compose this with the 
      ! cell_methods from the history averaging process

      ! Possibly time:point should just be left out

      select case (vinfo%ave_type(ifile))
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
         if ( vinfo%ave_type(ifile) /= hist_fixed ) then
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
   subroutine create_ncfile ( filename, nxhis, nyhis, nlev, ol, multilev,            &
                 use_plevs, use_meters, use_depth, use_hyblevs, basetime,            &
                 coord_heights, ncid, dims, dimvars, source, extra_atts, calendar,   &
                 nsoil, zsoil )

      use mpidata_m
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nxhis, nyhis, nlev, ol
      logical, intent(in) :: multilev, use_plevs, use_meters, use_depth, use_hyblevs
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
            ierr = nf90_put_att ( ncid, dimvars%oz, "long_name", "ocean sigma_level" )
            call check_ncerr(ierr)
            ierr = nf90_put_att ( ncid, dimvars%oz, "units", "sigma_level" )
            call check_ncerr(ierr)
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
      
      integer :: ierr, ifile
     
      call START_LOG(closehist_begin)
      
      do ifile = 1,nhfiles
!        The hist_oave files are closed individually in writehist.
         !if ( ihtype(ifile) /= hist_oave ) then
            ierr = nf90_close ( histid(ifile) )
            call check_ncerr(ierr,"Error closing history file")
         !end if
      end do
      
      call END_LOG(closehist_end)
      
   end subroutine closehist
!-------------------------------------------------------------------

   subroutine savehist2D ( name, array, jlat )
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:) :: array
      integer, intent(in) :: jlat
      real, dimension(size(array),1,1) :: temp

      call START_LOG(savehist_begin)

!     Copy to the 3D array that savehist_work expects. The extra level of
!     subroutines is necessary to avoid multiple versions.
!     This copying is extra work if the variable isn't actually saved.
!     Check name first?
      temp(:,1,1) = array
      call savehist_work  ( name, temp, jlat, jlat )

      call END_LOG(savehist_end)

   end subroutine savehist2D
   
   subroutine gsavehist2D ( name, array )
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:,:) :: array
      real, dimension(size(array,1),size(array,2),1) :: temp

      call START_LOG(savehist_begin)      

!     Copy to the 3D array that savehist_work expects. The extra level of
!     subroutines is necessary to avoid multiple versions.
      temp(:,:,1) = array
      call savehist_work  ( name, temp, 1, size(array,2) )

      call END_LOG(savehist_end)

   end subroutine gsavehist2D
   
   subroutine savehist3D( name, array, jlat )
      use logging_m       
      character(len=*), intent(in) :: name
      real, dimension(:,:) :: array
      integer, intent(in) :: jlat
      real, dimension(size(array,1),1,size(array,2)) :: temp

      call START_LOG(savehist_begin)

      temp(:,1,:) = array
      call savehist_work  ( name, temp, jlat, jlat )
      
      call END_LOG(savehist_end)
      
   end subroutine savehist3D
   
   subroutine gsavehist3D( name, array )
      use logging_m
      character(len=*), intent(in) :: name
      real, dimension(:,:,:) :: array
      call START_LOG(savehist_begin)
      call savehist_work ( name, array, 1, size(array,2) )
      call END_LOG(savehist_end)
   end subroutine gsavehist3D
   
   subroutine savehist_work ( name, array, jlat1, jlat2 )

!     This routine actually saves the data to the history array.
      character(len=*), intent(in) :: name
      real, dimension(:,:,:) :: array
      integer, intent(in) :: jlat1, jlat2
      integer :: ifld, ifile, js, jn, nl, istart, iend
      integer :: nx
      
      if ( all( array == NF90_FILL_FLOAT ) ) return

!     Check that openhist has been called to allocate the array
!     nhfiles check is required because when no history is saved the array
!     isn't allocated.
      if ( nhfiles /= 0 .and. .not. allocated ( histarray ) ) then
         print*, " History error - must call openhist before savehist "
         stop
      end if

!     Loop over history files because variable could occur in several.
      do ifile = 1,nhfiles

         ifld = qindex_hname(name,inames(1:totflds),totflds,ifile)
         if ( ifld == 0 ) then
            print*, "Error - history variable ", name, " is not known. "
            stop
         end if
         if ( hist_debug > 5 ) then
            print*, "SAVEHIST2D ", name, ifld, ifile, &
                    histinfo(ifld)%used(ifile)
         end if
         if ( .not. histinfo(ifld)%used(ifile) ) then
            cycle
         end if

!        Optional check on array sizes
         if ( hist_debug > 0 ) then
            if ( double ) then
               if ( 2*size(histarray,1) /= size(array,1) ) then
                  print*, "ERROR (doublerow): History array has nlon = ", &
                          size(histarray,1), &
                          "Input array has nlon = ", size(array,1)
                  print*, "Variable ", trim(name)
                  stop
               end if
            else
               if ( size(histarray,1) /= size(array,1) ) then
                  print*, "ERROR: History array has nlon = ", &
                          size(histarray,1), &
                          "Input array has nlon = ", size(array,1)
                  print*, "Variable ", trim(name)
                  stop
               end if
            end if
         end if

!        Need a check that if jlat1 /= jlat2 they cover the whole dimension.
!        Allowing partial coverage has the risk of partial overlaps.

         if ( jlat1 > jlat2 ) then
            print*, "jlat out of order", jlat1, jlat2
            stop
         end if
         if ( jlat1 < 1  .or. jlat2 > size(histarray,2) ) then
            print*, "jlat out of bounds ", jlat1, jlat2
            print*, "Variable ", trim(name)
            stop
         end if

         if ( double .and. jlat1 /= jlat2 ) then
            print*, "For double option must have jlat1 = jlat2", jlat1, jlat2
            stop
         end if

         if ( double .and. size(array,2) /= 1 ) then
            print*, "Error: double option and size(array,2) = ", size(array,2)
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

         istart = histinfo(ifld)%ptr(ifile)
         iend = istart+nl-1

         if ( double ) then
!           In this case the size of the second dimension must be 1
            js = jlat1
            jn = nyhis+1-jlat1
            nx = size(histarray,1)
            select case ( histinfo(ifld)%ave_type(ifile))
            case ( hist_ave )
               histarray(:,jn,istart:iend) = &
                    histarray(:,jn,istart:iend) + array(1:nx,1,:)
               histarray(:,js,istart:iend) = &
                    histarray(:,js,istart:iend) + array(nx+1:,1,:)
            case ( hist_max )
               histarray(:,jn,istart:iend) = &
                  max ( histarray(:,jn,istart:iend), array(1:nx,1,:) )
               histarray(:,js,istart:iend) = &
                  max ( histarray(:,js,istart:iend), array(nx+1:,1,:) )
            case ( hist_min ) 
               histarray(:,jn,istart:iend) = &
                  min ( histarray(:,jn,istart:iend), array(1:nx,1,:) )
               histarray(:,js,istart:iend) = &
                  min ( histarray(:,js,istart:iend), array(nx+1:,1,:) )
            case ( hist_inst, hist_fixed ) 
               histarray(:,jn,istart:iend) = array(1:nx,1,:)
               histarray(:,js,istart:iend) = array(nx+1:,1,:)
            case default
               print*, "Internal error in history, unknown ave_type", &
                    histinfo(ifld)%ave_type(ifile)
               stop
            end select
         else
            select case ( histinfo(ifld)%ave_type(ifile))
            case ( hist_ave )
               histarray(:,jlat1:jlat2,istart:iend) =  &
                    histarray(:,jlat1:jlat2,istart:iend) + array
            case ( hist_max ) 
               histarray(:,jlat1:jlat2,istart:iend) =  &
                    max ( histarray(:,jlat1:jlat2,istart:iend), array )
            case ( hist_min )
               histarray(:,jlat1:jlat2,istart:iend) =  &
                    min ( histarray(:,jlat1:jlat2,istart:iend), array )
            case ( hist_inst, hist_fixed ) 
               histarray(:,jlat1:jlat2,istart:iend) =  array
            case default
               print*, "Internal error in history, unknown ave_type", &
                    histinfo(ifld)%ave_type(ifile)
               stop
            end select
         end if

         if ( jlat1 == 1 ) then
            histinfo(ifld)%count(ifile) = histinfo(ifld)%count(ifile) + 1
         end if

         if ( hist_debug >= 4 ) then
            if ( jhdb==jlat1 ) then
               if ( khdb <= nl ) then
                  print*, "History at point ", name,  &
                    histarray(ihdb,jlat1,istart+khdb-1), array(ihdb,1,khdb)
               else
                  print*, "History at point ", name,  &
                    histarray(ihdb,jlat1,istart), array(ihdb,1,1)
               end if
            else if ( jlat1 == 1 .and. jlat2 >= jhdb ) then
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

      end do
   end subroutine savehist_work

!-------------------------------------------------------------------
   subroutine writehist ( istep, endofrun, year, month, interp, time, time_bnds )

      use mpidata_m
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
      integer, dimension(4) :: start3D, count3D
      integer, dimension(3) :: start2D, count2D
      integer :: istart, iend
      integer :: ave_type, nlev, ifile, count, ncid
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

      do ifile = 1,nhfiles

!        Check whether this file should be written on this step.
         if ( hist_debug > 1 ) then
            print*, "Checking frequency ", ifile, istep, hfreq(ifile)
         end if

!        For files containing time averages, calculate the appropriate
!        average time value. Don't increment this on the final call when
!        endofrun is set.
         doinc = .true.
         if ( present(endofrun) ) then
            doinc = .not. endofrun
         end if
         if ( doinc .and. ihtype(ifile) == hist_ave .and. present(time) ) then
            avetime(ifile) = avetime(ifile) + time
            timecount(ifile) = timecount(ifile) + 1
            if ( hist_debug > 0 ) then
               print*, "AVETIME ", ifile, time, avetime(ifile), timecount(ifile)
            end if
         end if
         if ( doinc .and. ihtype(ifile) == hist_ave .and. present(time_bnds) ) then
            
            ! Bounds for average should be the extreme bounds
            avetime_bnds(1,ifile) = min(avetime_bnds(1,ifile), time_bnds(1))
            avetime_bnds(2,ifile) = max(avetime_bnds(2,ifile), time_bnds(2))
            if ( hist_debug > 0 ) then
               print*, "AVETIME_BNDS ", ifile, time_bnds, avetime_bnds(:,ifile)
            end if
         end if

!        Files with hfreq=0 are written once at end of run
         if ( hfreq(ifile) == 0 ) then
            if ( present(endofrun) ) then
               if ( .not. endofrun ) then
                  cycle
               end if
            else ! No endofrun argument
               cycle
            end if
         else   ! hfreq(ifile) /= 0
            if ( modulo(istep, hfreq(ifile) ) /= 0 ) then
               cycle
            end if
!           Files for which hfreq is non-zero shouldn't be written at the end.
            if ( present(endofrun) ) then
               if ( endofrun ) then
                  cycle
               end if
            end if
         end if

         if ( hist_debug > 0 ) then
            print*, "Writing history file ", ifile
         end if

         histset(ifile) = histset(ifile) + 1

         if ( myid == 0 ) then

            !if ( ifile == 1 .and. ihtype(ifile) == hist_oave ) then
            !   if ( .not. ( present(year) .and. present(month) ) ) then
            !      print*, "Error year and month arguments to writehist are "
            !      print*, "required for old format files"
            !      stop
            !   end if
!           !   Add a check. Interpolation not supported with old format.
            !   call oldwrite(nxhis,nyhis,year,month)
            !   cycle
            !end if

            ncid = histid(ifile)

            ierr = nf90_inq_varid (ncid, "time", vid )
            call check_ncerr(ierr, "Error getting time id")
            if ( present(time) ) then
               if ( ihtype(ifile) == hist_ave) then
                  ierr = nf90_put_var ( ncid, vid,   &
                       avetime(ifile)/timecount(ifile), start=(/histset(ifile)/) )
               else
                  ierr = nf90_put_var ( ncid, vid, time, start=(/histset(ifile)/))
               end if
            else
               ierr = nf90_put_var ( ncid, vid, &
                                     real(histset(ifile)*hfreq(ifile)), &
                                     start=(/histset(ifile)/) )
            end if
            call check_ncerr(ierr, "Error writing time")
            if ( cf_compliant .and. present(time_bnds) ) then
               ierr = nf90_inq_varid (ncid, "time_bnds", vid )
               call check_ncerr(ierr, "Error getting time_bnds id")
               if ( ihtype(ifile) == hist_ave) then
                  ierr = nf90_put_var ( ncid, vid, avetime_bnds(:,ifile), start=(/1,histset(ifile)/))
               else
                  ierr = nf90_put_var ( ncid, vid, time_bnds, start=(/1,histset(ifile)/))
               end if
               call check_ncerr(ierr, "Error writing time_bnds")
            end if

            start2D = (/ 1, 1, histset(ifile) /)
            count2D = (/ nxhis, nyhis, 1 /)
            count3D = (/ nxhis, nyhis, 1, 1 /)
         
         end if

#ifdef usempi3
         cnt = 0
!        first pass
!        find total number of levels to store
         do ifld = 1,totflds
            if ( .not. histinfo(ifld)%used(ifile) ) then
               cycle
            end if
            ave_type = histinfo(ifld)%ave_type(ifile)
            nlev = histinfo(ifld)%nlevels
            count = histinfo(ifld)%count(ifile)

!           Only write fixed variables in the first history set
            if ( histset(ifile) > 1 .and. ave_type == hist_fixed ) then
               cycle
            end if

            istart = histinfo(ifld)%ptr(ifile)
            iend = istart + nlev - 1
            !if ( ave_type == hist_ave .or. ave_type == hist_oave ) then
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

!           Even multilevel variables are written one level at a time
            if ( count /= 0 ) then
              cnt = cnt + iend - istart + 1
            end if
               
         end do ! Loop over fields

         slab = ceiling(1.0d0*cnt/nproc)
         maxcnt = cnt
         interp_nproc = ceiling(1.0d0*maxcnt/slab)
         offset = nproc - interp_nproc
         if ( myid >= offset ) then
            allocate( hist_a(pil,pjl*pnpan,pnproc,1+slab*(myid-offset):slab*(myid-offset+1)) )
            allocate( hist_g(nx_g,ny_g) )
         end if
         allocate( k_indx(maxcnt) )

         cnt = 0
!        second pass
!        create the array used to index histarray
         do ifld = 1,totflds
            if ( .not. histinfo(ifld)%used(ifile) ) then
               cycle
            end if
            ave_type = histinfo(ifld)%ave_type(ifile)
            nlev = histinfo(ifld)%nlevels
            count = histinfo(ifld)%count(ifile)

!           Only write fixed variables in the first history set
            if ( histset(ifile) > 1 .and. ave_type == hist_fixed ) then
               cycle
            end if

            istart = histinfo(ifld)%ptr(ifile)
            iend = istart + nlev - 1

!           Even multilevel variables are written one level at a time
            if ( count /= 0 ) then
               do k = istart,iend
                  cnt = cnt + 1
                  k_indx(cnt) = k
               end do   ! k loop
            end if  

         end do ! Loop over fields

!        now do the gather wrap
         call gather_wrap(histarray,hist_a,slab,offset,maxcnt,k_indx)

         cnt = 0
!        third pass
!        perform the interpolation and write the data
         do ifld = 1,totflds
            if ( .not. histinfo(ifld)%used(ifile) ) then
               cycle
            end if
            ave_type = histinfo(ifld)%ave_type(ifile)
            nlev = histinfo(ifld)%nlevels
            count = histinfo(ifld)%count(ifile)
            vid = histinfo(ifld)%vid(ifile)

!           Only write fixed variables in the first history set
            if ( histset(ifile) > 1 .and. ave_type == hist_fixed ) then
               cycle
            end if

            istart = histinfo(ifld)%ptr(ifile)
            iend = istart + nlev - 1

!           Even multilevel variables are written one level at a time
            do k = istart, iend

               if ( count == 0 ) then
                  if ( hbytes(ifile) == 2 ) then
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
                     if ( hbytes(ifile) == 2 ) then
                       addoff = histinfo(ifld)%addoff(ifile)
                       sf = histinfo(ifld)%scalef(ifile)
                       umin = sf * vmin + addoff
                       umax = sf * vmax + addoff
                       where ( fpequal(htemp, missing_value) )
                          htemp = NF90_FILL_SHORT
                       elsewhere
                          ! Put the scaled array back in the original and let
                          ! netcdf take care of conversion to int2
                          htemp = nint((max(umin,min(umax,htemp))-addoff)/sf)
                       end where
                     else
                       where ( fpequal(htemp, missing_value) )
                          htemp = missing_value_cordex 
                       end where    
                     end if    
                  end if
                   
                  if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                     start3D = (/ 1, 1, k+1-istart, histset(ifile) /)
                     ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
                  else
                     ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
                  end if
                  call check_ncerr(ierr, "Error writing history variable "//histinfo(ifld)%name )

               end if
                 
            end do   ! k loop
               
!           Zero ready for next set
            histarray(:,:,istart:iend) = initval(ave_type)
!           Reset the count variable
            histinfo(ifld)%count(ifile) = 0

         end do ! Loop over fields

         if ( myid >= offset ) then
            deallocate(hist_a, hist_g)
         end if
         deallocate(k_indx)
#else
         do ifld = 1,totflds
            if ( .not. histinfo(ifld)%used(ifile) ) then
               cycle
            end if
            ave_type = histinfo(ifld)%ave_type(ifile)
            nlev = histinfo(ifld)%nlevels
            count = histinfo(ifld)%count(ifile)
            vid = histinfo(ifld)%vid(ifile)

!           Only write fixed variables in the first history set
            if ( histset(ifile) > 1 .and. ave_type == hist_fixed ) then
               cycle
            end if

            istart = histinfo(ifld)%ptr(ifile)
            iend = istart + nlev - 1
            !if ( ave_type == hist_ave .or. ave_type == hist_oave ) then
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
               allocate( hist_a(pil,pjl*pnpan,nlev,pnproc) )
               allocate( hist_g(nx_g,ny_g) )
            else
               allocate( hist_a(0,0,0,0) )
               allocate( hist_g(0,0) )
            end if    

            call gather_wrap(histarray(:,:,istart:iend),hist_a)
            
            if ( myid == 0 ) then
                
!              Even multilevel variables are written one level at a time
               do k = istart, iend

                  do ip = 0,pnproc-1   
                     do n = 0,pnpan-1
                        hist_g(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g) = &
                           hist_a(1:pil,1+n*pjl:(n+1)*pjl,k-istart+1,ip+1)
                     end do
                  end do
  
                  if ( count == 0 ) then
                     if ( hbytes(ifile) == 2 ) then
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

                     if ( hbytes(ifile) == 2 ) then
                        addoff = histinfo(ifld)%addoff(ifile)
                        sf = histinfo(ifld)%scalef(ifile)
                        umin = sf * vmin + addoff
                        umax = sf * vmax + addoff
                        where ( fpequal(htemp, missing_value) )
                           htemp = NF90_FILL_SHORT
                        elsewhere
                           ! Put the scaled array back in the original and let
                           ! netcdf take care of conversion to int2
                           htemp = nint((max(umin,min(umax,htemp))-addoff)/sf)
                        end where
                     else
                        where ( fpequal(htemp, missing_value) )
                           htemp = missing_value_cordex
                        end where  
                     end if
                  
                  end if

                  if ( nlev > 1 .or. histinfo(ifld)%multilev ) then
                     start3D = (/ 1, 1, k+1-istart, histset(ifile) /)
                     ierr = nf90_put_var ( ncid, vid, htemp, start=start3D, count=count3D )
                  else
                     ierr = nf90_put_var ( ncid, vid, htemp, start=start2D, count=count2D )
                  end if
                  call check_ncerr( ierr, "Error writing history variable "//histinfo(ifld)%name )
                 
               end do   ! k loop
               
            end if ! myid == 0

            deallocate( hist_a, hist_g )
            
!           Zero ready for next set
            histarray(:,:,istart:iend) = initval(ave_type)
!           Reset the count variable
            histinfo(ifld)%count(ifile) = 0

         end do ! Loop over fields
#endif

         avetime(ifile) = 0.0
         ! Initialisation so min/max work
         avetime_bnds(:,ifile) = (/huge(1.), -huge(1.)/)
         timecount(ifile) = 0

!        Sync the file so that if the program crashes for some reason 
!        there will still be useful output.
#ifdef outsync
         if ( myid == 0 ) then
            ierr = nf90_sync ( ncid )
            call check_ncerr(ierr, "Error syncing history file")
         end if
#endif

      end do ! Loop over files

      call END_LOG(writehist_end)

   end subroutine writehist
!-------------------------------------------------------------------
   subroutine oldwrite(nxhis, nyhis, year, month )

      use mpidata_m

!  Write old format files containing only a single field at a single level
      integer, intent(in) :: nxhis, nyhis, year, month
      integer :: ierr, vid, ifld, ilev, yrid
      integer, dimension(4) :: astart, acount
      integer :: istart
      integer :: ave_type, nlev, ifile, count, ncid
      character(len=200) :: filename
      character(len=MAX_NAMELEN) :: vname
      integer, dimension(nxhis,nyhis) :: imat
      real :: addoff, sf, umin, umax
      real, dimension(nxhis,nyhis) :: rmat
      integer :: year1

      if ( myid /= 0 ) return

      ifile = 1
      do ifld = 1, totflds

         if ( .not. histinfo(ifld)%used(ifile) ) then
            cycle
         end if

         ave_type = histinfo(ifld)%ave_type(ifile)
         nlev = histinfo(ifld)%nlevels
         count = histinfo(ifld)%count(ifile)
         if ( count == 0 ) then
            print*, "Error: history variable ", &
                    trim(histinfo(ifld)%name), &
                    " has never been set, skipping"
            cycle
         end if
         do ilev = 1,nlev
            if ( nlev == 1 ) then
               vname = histinfo(ifld)%name
            else
               write(vname,"(a,i2.2)") trim(histinfo(ifld)%name), ilev
            end if
            write(filename,"(a,a,a,a)" ) trim(oldprefix), trim(vname), &
                                         trim(filesuffix), ".nc"

            if ( hist_debug > 0 ) then
               print*, "Opening ", filename
            end if
            ierr = nf90_open(filename, NF90_WRITE, ncid)
            call check_ncerr ( ierr, "Error opening history file "//filename )
            ierr = nf90_inq_varid (ncid, vname, vid )
            call check_ncerr(ierr, "Error getting variable name "//vname)

            istart = histinfo(ifld)%ptr(ifile) + ilev - 1

            if ( ave_type == hist_ave ) then
               where( histarray(:,:,istart) /= NF90_FILL_FLOAT ) 
                  histarray(:,:,istart) = histarray(:,:,istart) / count
               end where   
            end if
            if ( histinfo(ifld)%output_scale /= 0 ) then
               histarray(:,:,istart) = histarray(:,:,istart) * &
                   histinfo(ifld)%output_scale
            end if
            if ( hist_debug >= 4 ) then
               print*, "History written at point ", vname, &
                    histarray(ihdb,jhdb,istart)
            end if

            ierr = nf90_inq_varid (ncid, "year", yrid )
            call check_ncerr(ierr, "Error getting year id")
            ierr = nf90_get_att (ncid, yrid, "valid_min", year1 )
            call check_ncerr(ierr,"Error getting year valid min attribute")
            astart = (/ 1, 1, month, year+1-year1 /)
            acount = (/ nxhis, nyhis, 1, 1 /)

!           Note that this routine can only be called on the first file
            if ( hbytes(1) == 2 ) then
               ierr = nf90_get_att (ncid, vid, "add_offset", addoff )
               call check_ncerr(ierr,"Error getting add_offset attribute")
               ierr = nf90_get_att (ncid, vid, "scale_factor", sf )
               call check_ncerr(ierr,"Error getting scale_factor attribute")
               umin = sf * vmin + addoff
               umax = sf * vmax + addoff
               where ( fpequal(histarray(:,:,istart), missing_value) )
                  imat = NF90_FILL_SHORT
               elsewhere
                  imat = &
                   nint((max(umin,min(umax,histarray(:,:,istart)))-addoff)/sf)
               endwhere
!              Use nf_put_vara_int and let netcdf itself convert to int2.
!              This should be more portable than using integer*2 in f90.
               ierr = nf90_put_var ( ncid, vid, imat, start=astart, &
                                     count=acount )
            else
               where ( fpequal(histarray(:,:,istart), missing_value) )
                  rmat = missing_value_cordex
               elsewhere
                  rmat = histarray(:,:,istart)
               end where   
               ierr = nf90_put_var ( ncid, vid, rmat, start=astart, &
                                     count=acount )
            end if
            call check_ncerr(ierr, &
                 "Error writing history variable "//vname )

!           Zero ready for next set
            histarray(:,:,istart) = initval(ave_type)
!           Reset the count variable
            histinfo(ifld)%count(ifile) = 0

!           Write the year 
            ierr = nf90_put_var ( ncid, yrid, year, start=(/year+1-year1/) )
            call check_ncerr(ierr, "Error writing year")

            ierr = nf90_close(ncid)
            call check_ncerr(ierr,"Error closing file "//filename)

         end do ! Loop over levels
      end do ! Loop over fields

   end subroutine oldwrite

!-------------------------------------------------------------------
   function initval ( ave_type ) result (val)
      integer, intent(in) :: ave_type
      real :: val
!     Set a suitable initialisation value depending on the ave_type
      select case ( ave_type)
      case ( hist_inst, hist_fixed )
         val = NF90_FILL_FLOAT
      !case ( hist_ave, hist_oave )
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
      integer :: ave_type, nlev, ifile

      if ( hist_debug > 0 ) then
         print*, "Resetting history"
      end if
      do ifile=1,nhfiles

         do ifld = 1,totflds
            if ( .not. histinfo(ifld)%used(ifile) ) then
               cycle
            end if
            ave_type = histinfo(ifld)%ave_type(ifile)
            nlev = histinfo(ifld)%nlevels
            istart = histinfo(ifld)%ptr(ifile)
            iend = istart+nlev-1

!           Zero ready for next set
            histarray(:,:,istart:iend) = initval(ave_type)
!           Reset the count variable
            histinfo(ifld)%count(ifile) = 0

         end do ! Loop over fields

      end do ! Loop over files

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
      integer :: istart, iend, ip, k, n, ierr
      integer, dimension(maxcnt), intent(in) :: k_indx
      real, dimension(:,:,:), intent(in) :: histarray
      real, dimension(size(histarray,1),size(histarray,2),slab) :: histarray_tmp
      real, dimension(:,:,:,:), pointer, contiguous, intent(out) :: hist_a
      real, dimension(pil*pjl*pnpan*lproc*slab*nproc), target :: hist_a_tmp
      real, dimension(:,:,:,:), pointer, contiguous :: hist_a_remap, hist_a_tmp_remap
   
      call START_LOG(gatherwrap_begin)
      
      do ip = 0,nproc-1
         istart = 1 + slab*(ip-offset) 
         if ( istart > 0 ) then
            iend = istart + slab - 1
            iend = min( iend, maxcnt )
            do k = istart,iend
               histarray_tmp(:,:,k-istart+1) = histarray(:,:,k_indx(k))
            end do   
            call START_LOG(mpigather_begin)
            call MPI_Gather(histarray_tmp(:,:,:), pil*pjl*pnpan*lproc*(iend-istart+1), MPI_REAL,    &
                            hist_a_tmp, pil*pjl*pnpan*lproc*(iend-istart+1), MPI_REAL,              &
                            ip, comm_world, ierr)
            call END_LOG(mpigather_end)
         end if
      end do
      
      istart = 1 + slab*(myid-offset)
      if ( istart > 0 ) then
         iend = slab*(myid-offset+1)
         iend = min( iend, maxcnt )          
         hist_a_remap(1:pil,1:pjl*pnpan*lproc,1:nproc,istart:iend) => hist_a
         hist_a_tmp_remap(1:pil,1:pjl*pnpan*lproc,istart:iend,1:nproc) =>    &
             hist_a_tmp(1:pil*pjl*pnpan*lproc*(iend-istart+1)*nproc)
         do n = 1,nproc
            do k = istart,iend 
               hist_a_remap(:,:,n,k) = hist_a_tmp_remap(:,:,k,n)
            end do
         end do
      end if          
     
      call END_LOG(gatherwrap_end)
   
   end subroutine gather_wrap
   
   subroutine sendrecv_wrap(htemp,cnt,slab,offset)
      use mpidata_m, only : nproc, lproc, myid, comm_world
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif  
      real, dimension(:,:), intent(inout) :: htemp
      integer, intent(in) :: cnt, slab, offset
      integer :: rrank, sizehis, ierr
   
      rrank = ceiling(1.0d0*cnt/slab) - 1 + offset
      if ( rrank /= 0 ) then
         if ( myid == 0 ) then
            sizehis = size(htemp, 1)*size(htemp, 2) 
            call MPI_Recv(htemp,sizehis,MPI_REAL,rrank,1,comm_world,MPI_STATUS_IGNORE,ierr)
         else if ( myid == rrank ) then
            sizehis = size(htemp, 1)*size(htemp, 2) 
            call MPI_Send(htemp,sizehis,MPI_REAL,0,1,comm_world,ierr)
         end if
      end if
       
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
      call END_LOG(gatherwrap_end)
      
   end subroutine gather_wrap
#endif


!-------------------------------------------------------------------
   function needfld ( name ) result (needed)
      character(len=*), intent(in) :: name
      logical :: needed
      integer :: ifld, ifile

!     Check if name is in the list for any of the history files
      needed = .false.
      do ifile=1,nhfiles
         ifld = qindex_hname ( name, inames(1:totflds), totflds, ifile)
         if ( ifld == 0 ) then
            ! Name not known at all in this case
            needed = .false.
            exit
         else if ( histinfo(ifld)%used(ifile) ) then
            needed = .true.
            exit   ! No need to check any further.
         end if
      end do
   end function needfld
   
!-------------------------------------------------------------------
   function bindex_hname(name, table,nflds) result(ifld)
      character(len=*), intent(in)    :: name
      integer, intent(in) :: nflds
      integer ( bigint ), intent(in), dimension(nflds) :: table
      integer :: ifld
      integer :: i, lower, upper
      integer (bigint) :: key, fac

!  Lookup "key" in "table" of integers using a binary search.
!  This assumes that the list has been sorted but it doesn't test for this.

!      key = hashkey ( name )
!     hashkey inlined by hand      
      key = 0
      fac = 1
!     This makes numerical order the same as alphabetical order
      do i=min(MAX_KEYLEN,len_trim(name)),1,-1
         ! Netcdf allowed characters for variables names are in 
         ! range 48 to 122 (alphanumeric and _)
         key = key + fac*(ichar(name(i:i))-48)
         fac = fac * 75
      end do

      ifld = 0
      lower = 1
      upper = nflds
      do
         i = (lower+upper)/2
!         print*, "Search", lower, upper, i, table(lower), table(upper), table(i), key
         if ( key < table(i) ) then
            upper = i-1
         else if ( key > table(i) ) then
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
   function qindex_hname ( name, table, nflds, ifile ) result(ifld)
      character(len=*), intent(in)    :: name
      integer, intent(in) :: nflds
      integer ( bigint ), intent(in), dimension(nflds) :: table
      integer, intent(in) :: ifile
      integer :: ifld
!     Both of these are initially set to one because this is a safe value.
!     It doesn't matter that this points to the wrong thing because it's
!     always checked.
!     Index of the previously found name
      integer, save, dimension(MAX_HFILES) :: prev = 1
!     Index of the successor of a given index
      integer, save, dimension(nfmax,MAX_HFILES) :: successor = 1

!     First check the successor of the previous name. This will generally 
!     work because things are done in a fixed order. It will only break on
!     the first latitude of the radiation step.

      ifld = successor(prev(ifile),ifile)
      if ( histinfo(ifld)%name /= name ) then
         ifld = bindex_hname ( name, table, nflds )
         if ( ifld == 0 ) then
            ! Name is not known, return without updating internal fields
            return
         end if
!        Redefine the successor.
         successor(prev(ifile),ifile) = ifld
      end if

!     Set this for next time
      prev(ifile) = ifld
         

   end function qindex_hname
!-------------------------------------------------------------------
   function hashkey ( name ) result (key)
!  Generate an integer from a string
      character(len=*), intent(in) :: name
      integer (bigint) :: key
      integer :: i
      integer (bigint) :: fac

      key = 0
      fac = 1
!     This makes numerical order the same as alphabetical order
      do i=min(MAX_KEYLEN,len_trim(name)),1,-1
         ! Netcdf allowed characters for variables names are in 
         ! range 48 to 122 (alphanumeric and _)
         key = key + fac*(ichar(name(i:i))-48)
         fac = fac * 75
      end do
      if (hist_debug > 2 ) then
         print*, "  HASHKEY ", name, key
      end if
   end function hashkey
   
end module history
