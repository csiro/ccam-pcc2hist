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

   use netcdf
   use ncutils_m, only : check_ncerr
   use utils_m, only : fpequal
   implicit none

!  Make everything private so that the internal variables of the history 
!  module can only be used by the appropriate routines.
   private   

!  Only these routine names need to be public
   public :: savehist, openhist, closehist, writehist, addfld, &
             inithist, histnamelist, needfld, clearhist, hstring

!  Routine that provide access to the namelist control variables.
!  Access is only possible through a routine so that the underlying data
!  structure may be changed without affecting user programs.
   public :: set_htype,  set_hfreq, set_hnames, set_hbytes

   public :: set_missval

!  Private internal routines
   private :: bindex_hname, initval, sortlist, create_ncvar,            &
              create_ncfile, create_oldfile, oldwrite, savehist2D,      &
              savehist3D, hashkey, qindex_hname, savehist_work,         &
              gsavehist2D, gsavehist3D

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
   integer, parameter :: MAX_NAMELEN = 10

!  Maximum length of string used for key generation. Names must be unique
!  within this length.
   integer, parameter :: MAX_KEYLEN = 10

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
   integer, parameter :: nfmax = 500

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
                         hist_fixed=5, hist_oave=6

!  Valid range for 16 bit values
   integer, parameter :: vmin=-32500, vmax=32500

   real :: missing_value = NF90_FILL_FLOAT
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

!  MPI working arrays
   real, dimension(:,:,:,:), allocatable, save, private :: hist_a
   real, dimension(:,:,:), allocatable, save, private :: hist_g

contains

!-------------------------------------------------------------------

   subroutine hstring(hstr)

      ! Create a string containing the username, machine name and date to 
      ! identify the creation of the history file
   
      character(len=*), intent(out) :: hstr
      character(len=100) :: logname, hostname
      integer, dimension(8) :: idatetime
      character(len=1000) :: str
   
      logname = ""
      call get_environment_variable("LOGNAME",logname)
      hostname = ""
      call get_environment_variable("HOSTNAME",hostname)
      call date_and_time(values=idatetime)
      
      ! Avoid overflow in the output string by writing to a long internal string
      str = ""
      write(str,'("Created by ",a, " on ", a, " ", i4.4,"-", i2.2,"-",i2.2," ",i2.2,":",i2.2,":", i2.2)') &
           trim(logname), trim(hostname), idatetime(1:3), idatetime(5:7)
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
         case ("oave")
            ihtype(j) = hist_oave
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
      case ("oave")
         ihtype(j) = hist_oave
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

   subroutine set_missval ( missval )
      real, intent(in) :: missval
      missing_value = missval
   end subroutine set_missval
   
!-------------------------------------------------------------------
   subroutine inithist ( nl, mstep )
      integer, intent(in) :: nl    ! Number of levels
      integer, intent(in) :: mstep ! Time step (minutes)
!
!     Initialise the standard list of history variables used in the Mark 3
!     model. This routine is just a convenience for this model and is not
!     required generally.
!
!     3D fields
      call addfld ( "t", "Air temperature", "K", 100., 400., nl, &
                     amip_name="ta", multilev=.true. )
      call addfld ( "u", "Eastward wind", "m/s", -100., 100.,nl, &
                     amip_name="ua", multilev=.true.  ) 
      call addfld ( "v", "Northward wind", "m/s", -100., 100.,nl, &
                     amip_name="va", multilev=.true.  )
      call addfld ( "q", "Specific humidity", "kg/kg", 0., 0.1,nl, multilev=.true.  )
      call addfld('als','surface albedo',' ', 0., 1.,1, std=.true.)
      call addfld('cld','total cloud', ' ', 0., 1.,1, std=.true.)
      call addfld('clh','high cloud',' ', 0., 1.,1, std=.true.)
      call addfld('cll','low cloud',' ', 0., 1.,1, std=.true.)
      call addfld('clm','middle cloud',' ', 0., 1.,1, std=.true.)
      call addfld('evp','evaporation','mm/day',0.,25.,1, std=.true., &
                  output_scale=1440./mstep)
      call addfld('rnd','precipitation','mm/day',0.,100.,1, std=.true., &
                  output_scale=1440./mstep)
      call addfld('hfl','sensible heat flux','W/m2',-200.,1000.,1, std=.true.)
      call addfld('tsu','surface temperature','K',100.,350.,1, std=.true.)
      call addfld('tsc','screen temperature','K',100.,350.,1, std=.true.)
      call addfld('wfg','soil moisture upper',' ',0.,1.,1, std=.true.)
      call addfld('wfb','soil moisture lower',' ',0.,1.,1, std=.true.)
      call addfld('snd','snow depth','cm', 0., 1000.,1, std=.true.)
      call addfld('sid','sea-ice depth', 'm', 0., 10.,1, std=.true.)
      call addfld('tb2','soil temp level 2','K', 100., 350.,1, std=.true.)
      call addfld('tb3','soil temp level 3','K', 100., 350.,1, std=.true.)
      call addfld('vmo','surface wind speed','m/s', 0., 50.,1, std=.true.)
      call addfld('tax','surface stress east','N/m2',-10.,10.,1, std=.true.)
      call addfld('tay','surface stress nth','N/m2',-10.,10.,1, std=.true.)
      call addfld('rnc','convective rainfall','mm/day',0.,100.,1, std=.true., &
                  output_scale=1440./mstep)
      call addfld('run','runoff','mm/day', 0., 25.,1, std=.true., &
                  output_scale=1440./mstep)
      call addfld('inr','canopy interception','mm/day', 0., 25.,1, std=.true.,&
                  output_scale=1440./mstep)
      call addfld('tgg','bare ground temp','K', 100., 350.,1 )
      call addfld('tgf','veg ground temp','K', 100., 350.,1 )
      call addfld('pev','potential evap','mm/day', 0., 250.,1, &
                  output_scale=1440./mstep)
      call addfld('sev','scaling evap','mm/day', 0., 50.,1, &
                  output_scale=1440./mstep)
!     Should this have additional scale factor too??
      call addfld('per','Soil percolation','mm/day', 0., 10.,1)
      call addfld('ico','ice concentration',' ', 0., 1.,1)
      call addfld('itf','ice-ocean heat flux','W/m2',-100.,100.,1)
      call addfld('isf','ice-ocean salt flux','m',-0.01,0.01,1)
      call addfld('psl','mean sea-level pressure', 'hPa', 925., 1075.,1, std=.true.)
      call addfld('psf','surface pressure', 'hPa', 0., 1200.,1 )
      call addfld('thd','mean daily max temp', 'K', 100., 350.,1)
      call addfld('tld','mean daily min temp', 'K', 100., 350.,1)
      call addfld('thg','bare ground Tmax', 'K', 100., 350.,1)
      call addfld('tlg','bare ground Tmin', 'K', 100., 350.,1)
      call addfld('thf','Veg ground Tmax', 'K', 100., 350.,1)
      call addfld('tlf','Veg ground Tmin', 'K', 100., 350.,1)
      call addfld('thm','extreme max temp','K', 100., 350.,1)
      call addfld('tlm','extreme min temp', 'K', 100., 350.,1)
      call addfld('icu','ice zonal velocity', 'm/s', -5., 5.,1)
      call addfld('icv','ice merid velocity','m/s', -5., 5.,1)
      call addfld('sno','Snowfall','mm/day', 0., 25.,1)
      call addfld('rev','Rain evaporation', 'mm/day', 0., 100.,1)
      call addfld('ssb','Snow sublimation', 'mm/day', 0., 50.,1)
      call addfld('clc','convective cloud', ' ', 0., 1.,1)
      call addfld('lwp','liquid water path', 'kg/m3', 0., 3.,1)
      call addfld('pwc','precipitable water', 'mm', 0., 200.,1)
      call addfld('ref','Effective radius for liquid clouds', 'um', 0., 20.,1)
      call addfld('cli','Liquid cloud fraction', ' ', 0., 1.,1)
      call addfld('dtm','Tmlo error','K', -10., 10.,1)
      call addfld('gro','Monthly ice growth', 'm', -10., 10.,1)
      call addfld('ire','Ice redistribution', 'm', -1., 1.,1)
      call addfld('ich','Ice advection', 'm', -1., 1.,1)
      call addfld('acd','Average drag', 'm', -10., 10.,1)
      call addfld("rhs", "Screen level relative humidity","percent",0., 100.,1)
      call addfld("v10m","10m windspeed","m/s",0.,100.,1)

!     SW radiation
      call addfld("sit","SW insolation at TOA","W/m^2",0.,1400.,1)
      call addfld('sot','SW out at top','W/m2', 0., 1000.,1, std=.true.)
      call addfld('soc','SW out clear sky','W/m2',0.,1000.,1, std=.true.)
      call addfld('sgn','Net SW at ground','W/m2',0.,1000.,1, std=.true.)
      call addfld('sgc','Net SW ground clear','W/m2',0.,1000.,1, std=.true.)
      call addfld('sgd','Downward SW ground','W/m2',0.,1000.,1, std=.true.)
!     LW radiation
      call addfld('rtu','LW out at top','W/m2', 0., 1000.,1, std=.true.)
      call addfld('rtc','LW out clear sky','W/m2',0.,1000.,1, std=.true.)
      call addfld('rgn','Net LW at ground','W/m2',-200.,500.,1, std=.true.)
      call addfld('rgc','Net LW ground clear','W/m2',-200.,500.,1, std=.true.)
      call addfld('rgd','Downward LW ground','W/m2',0.,1000.,1, std=.true.)

!     Ice scheme
      call addfld("wls","ice residual divergence","/s",0.,100.,1,std=.true.)
      call addfld("wdf","ice divergence removed by rheology","/s",-100.,100.,&
                   1, std=.true.)

!     Other useful fields to have
      call addfld('omega','Vertical velocity', 'Pa/s', -10., 10.,nl)

   end subroutine inithist
   
!---------------------------------------------------------------------------
   subroutine addfld(name, long_name, units, valid_min, valid_max,    &
                     nlevels, amip_name, ave_type, std, output_scale, &
                     int_type, multilev, std_name, soil, coord_height, cell_methods )
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
      real, intent(in), optional :: coord_height
      character(len=*), intent(in), optional :: cell_methods

!     Local variables corresponding to the optional arguments
      integer :: atype
      logical :: lstd
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
         case ("oave")
            atype = hist_oave
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
   subroutine openhist ( nx, ny, nl, sig, suffix, hlon, hlat, basetime,       &
                         doublerow, year, nxout, nyout, source, histfilename, &
                         pressure, height, extra_atts, hybrid_levels, anf,    &
                         bnf, p0, calendar, nsoil, zsoil )
!
!     Create netCDF history files using information in histinfo array.
!
      use mpidata_m

      integer, intent(in) :: nx, ny, nl
      real, intent(in), dimension(:) :: sig
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
      logical :: used, multilev, use_plevs, use_hyblevs, use_meters
      integer, dimension(totflds) :: coord_heights
      integer :: kc, ncoords, k, pkl
      logical :: soil_used
      real :: dx, dy
      
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
            
            if ( hnames(ivar,ifile) == "all" ) then
               histinfo(1:totflds)%used(ifile) = .true.
               exit ! All defined so no point checking further
                    ! Except when xnames is used.
            else if ( hnames(ivar,ifile) == "std" ) then
!              Use or here so as not to overwrite variables that have been
!              set explicitly.
               histinfo(1:totflds)%used(ifile) = &
                histinfo(1:totflds)%used(ifile) .or. histinfo(1:totflds)%std
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
         do ivar=1,totflds
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

!        First file may be old average format. If this is the case it has
!        to be handled differently
         if ( ihtype(1) == hist_oave ) then

            ifile = 1
            if ( present(histfilename) ) then
               ! Use histfile as a path in this case
               oldprefix = trim(histfilename) // "/s"
            else
               oldprefix = "s"
            end if

            do ifld = 1, totflds

               if ( .not. histinfo(ifld)%used(ifile) ) then
                  cycle
               end if

!              lookup again to get long name and units
               vname = histinfo(ifld)%name
               longname = histinfo(ifld)%long_name
               units = histinfo(ifld)%units

               do ilev=1,histinfo(ifld)%nlevels
                  if ( histinfo(ifld)%nlevels == 1 ) then
                     vname = histinfo(ifld)%name
                  else
                     write(vname,"(a,i2.2)") trim(histinfo(ifld)%name), ilev
                  end if
                  write(filename,"(a,a,a,a)" ) trim(oldprefix), trim(vname), &
                                               trim(suffix), ".nc"
!                 Check if this file exists. If it does there's no more to
!                 do, if not go on to create it.
                  inquire(file=filename, exist=used)
                  if ( used ) then
                     if ( hist_debug > 0 ) then
                        print*, "Using existing file ", filename
                     end if
                     cycle
                  else

                     if ( .not. present(year) ) then
                        print*, " Year argument to openhist is required for old format files"
                        stop
                     end if
                     call create_oldfile ( filename, nxhis, nyhis, hbytes(ifile),&
                                           vname, longname, units, &
                                           histinfo(ifld)%valid_min, &
                                           histinfo(ifld)%valid_max, &
                                           hlat, hlon, year )

                  end if

               end do
            end do
            istart = 2
         else
            istart = 1
         end if ! ihtype(1) == hist_oave

!        The rest of the history files are much simpler with multilevel variables
         do ifile = istart, nhfiles

            soil_used = .false.
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
            if ( soil_used ) then
               ! Better to define a new local nsoil variable?
               call create_ncfile ( filename, nxhis, nyhis, size(sig), multilev,                   &
                    use_plevs, use_meters, use_hyblevs, basetime, coord_heights(1:ncoords), ncid,  &
                    dims, dimvars, source, extra_atts, calendar, nsoil, zsoil )
            else
               call create_ncfile ( filename, nxhis, nyhis, size(sig), multilev,                   &
                    use_plevs, use_meters, use_hyblevs, basetime, coord_heights(1:ncoords), ncid,  &
                    dims, dimvars, source, extra_atts, calendar)
            end if
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
            if ( soil_used .and. present(zsoil) ) then
               ierr = nf90_put_var ( ncid, dimvars%zsoil, zsoil )
               call check_ncerr(ierr,"Error writing depths")
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
            end if

            do kc=1,ncoords
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

!           Sync the file so that if the program crashes for some reason 
!           there will still be useful output.
#ifdef outsync
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
         if ( myid == 0 ) then
            pkl = size(sig)
            allocate( hist_a(pil,pjl*pnpan,pkl,pnproc) )
            allocate( hist_g(nx,ny,pkl) )
         else
            allocate( hist_a(0,0,0,0) )
            allocate( hist_g(0,0,0) )
         end if
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
               iend = istart + histinfo(ifld)%nlevels-1
               histarray(:,:,istart:iend) = &
                  initval(histinfo(ifld)%ave_type(ifile))
            end if
         end do
      end do

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

      character(len=80) :: cell_methods, coord_name
      integer :: ierr, vtype, vid, zdim

      if ( myid /=0 ) return

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
         print*, "Creating variable ", vinfo%name, vinfo%ave_type(ifile), vinfo%nlevels, vinfo%soil
         print*, "DIMID", dims%x, dims%y, dims%z, dims%t, dims%zsoil
      end if
      if ( vinfo%nlevels > 1 .or. vinfo%multilev ) then
         if (vinfo%soil) then
            zdim = dims%zsoil
         else
            zdim = dims%z
         end if
         if ( vinfo%ave_type(ifile) == hist_fixed ) then
!#ifdef usenc3
            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
                                (/ dims%x, dims%y, zdim /), vid )
         else
            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
                                (/ dims%x, dims%y, zdim, dims%t /), vid )
         end if
      else
         if ( vinfo%ave_type(ifile) == hist_fixed ) then
            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
                                (/ dims%x, dims%y /), vid )
         else
            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
                                (/ dims%x, dims%y, dims%t /), vid )
!#else
!            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
!                                (/ dims%x, dims%y, zdim /), &
!                                 vid, deflate_level=1 )
!         else
!            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
!                                (/ dims%x, dims%y, zdim, dims%t /), &
!                                 vid, deflate_level=1 )
!         end if
!      else
!         if ( vinfo%ave_type(ifile) == hist_fixed ) then
!            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
!                                (/ dims%x, dims%y /), vid, &
!                                deflate_level=1 )
!         else
!            ierr = nf90_def_var ( ncid, vinfo%name, vtype, &
!                                (/ dims%x, dims%y, dims%t /), &
!                                 vid, deflate_level=1 )
!#endif
         end if
      end if
      if ( ierr /= 0 ) then
         print*, "Error creating variable", vinfo
      end if
      call check_ncerr(ierr,"Error creating variable "// vinfo%name)
      vinfo%vid(ifile) = vid

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
      case (hist_ave, hist_oave)
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
!         ierr = nf_put_att_int ( ncid, vid, "_FillValue", NF_INT2, 1, int(NF90_FILL_SHORT) )
         ierr = nf90_put_att ( ncid, vid, "_FillValue", NF90_FILL_SHORT )
         call check_ncerr(ierr,"Error with fill value attribute")
!         ierr = nf_put_att_int ( ncid, vid, "missing_value", NF_INT2, 1, int(NF90_FILL_SHORT) )
         ierr = nf90_put_att ( ncid, vid, "missing_value", NF90_FILL_SHORT )
         call check_ncerr(ierr,"Error with missing value attribute")
      else
         ierr = nf90_put_att ( ncid, vid, "_FillValue", missing_value )
         call check_ncerr(ierr,"Error with fill value attribute")
         ierr = nf90_put_att ( ncid, vid, "missing_value", missing_value )
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
         call check_ncerr(ierr,"Error with fill value attribute")
      end if

   end subroutine create_ncvar
  
!---------------------------------------------------------------------------
   subroutine create_ncfile ( filename, nxhis, nyhis, nlev, multilev,                &
                 use_plevs, use_meters, use_hyblevs, basetime, coord_heights, ncid,  &
                 dims, dimvars, source, extra_atts, calendar, nsoil, zsoil )

      use mpidata_m
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nxhis, nyhis, nlev
      logical, intent(in) :: multilev, use_plevs, use_meters, use_hyblevs
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
      ierr = nf90_create(filename, NF90_CLOBBER, ncid)
      !ierr = nf90_create(filename, NF90_64BIT_OFFSET, ncid)
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
      if ( cf_compliant ) then
         ierr = nf90_put_att ( ncid, NF90_GLOBAL, "Conventions", "CF-1.0" )
         call check_ncerr(ierr)
      end if
!     Make sure it's a null string, otherwise len_trim won't find the end 
!     properly.
      histstr = "" 
      !call hstring(histstr)
      ierr = nf90_put_att ( ncid, NF90_GLOBAL, "history", histstr )
      if ( present(extra_atts)) then
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

      if ( present(nsoil) ) then
         if ( .not. present(zsoil) ) then
            print*, "Error - missing zsoil argument"
            stop
         end if
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
   
!---------------------------------------------------------------------------
   subroutine create_oldfile ( filename, nxhis, nyhis, hbytes,        &
                               vname, longname, units, xmin, xmax,    &
                               ylat, xlon, year )

      use mpidata_m
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nxhis, nyhis
      integer, intent(in) :: hbytes 
      character(len=*), intent(in) :: vname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real, intent(in) :: xmin, xmax
      real, dimension(:), intent(in) :: ylat
      real, dimension(:), intent(in) :: xlon
      integer, intent(in) :: year

      integer  ncid, lonid, latid, monid, yrid, vid, old_mode
      integer :: ierr
      integer :: londim, latdim, mondim, yrdim
      integer :: m
      real :: scalef, addoff

      if ( myid /=0 ) return

      if ( hist_debug > 0 ) then
         print*, "Creating file ", filename, xmin, xmax
      end if
      !ierr = nf90_create(filename, NF90_CLOBBER, ncid)
      ierr = nf90_create(filename, NF90_64BIT_OFFSET, ncid)
      call check_ncerr ( ierr, "Error in creating history file" )
               
!     Create dimensions, lon, lat, month and year
      ierr = nf90_def_dim ( ncid, "longitude", nxhis, londim )
      call check_ncerr(ierr,"Error creating lon dimension")
      ierr = nf90_def_dim ( ncid, "latitude", nyhis, latdim )
      call check_ncerr(ierr,"Error creating lat dimension")
      ierr = nf90_def_dim ( ncid, "month", 12, mondim )
      call check_ncerr(ierr,"Error creating month dimension")
      ierr = nf90_def_dim ( ncid, "year", NF90_UNLIMITED, yrdim )
      call check_ncerr(ierr,"Error creating time dimension")
      
!     Define attributes for the dimensions
      ierr = nf90_def_var ( ncid, "longitude", NF90_FLOAT, londim, lonid)
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, lonid, "units", "degrees_east" )
      call check_ncerr(ierr)

      ierr = nf90_def_var ( ncid, "latitude", NF90_FLOAT, latdim, latid )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, latid, "units", "degrees_north" )
      call check_ncerr(ierr)

      ierr = nf90_def_var ( ncid, "month", NF90_INT, mondim, monid )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, monid, "units", "months" )
      call check_ncerr(ierr)

      ierr = nf90_def_var ( ncid, "year", NF90_INT, yrdim, yrid )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, yrid, "units", "years" )
      call check_ncerr(ierr)
      ierr = nf90_put_att ( ncid, yrid, "valid_min", year )
      call check_ncerr(ierr,"Error setting valid min for year")

      if ( hbytes == 2 ) then
         ierr = nf90_def_var ( ncid, vname, NF90_INT2,  &
                             (/ lonid, latid, monid, yrid /), vid )
         call check_ncerr(ierr,"Error creating variable "// vname)
!     In the case of specific humidity we need to cheat and use a 
!     vertically varying scale factor.
!            qmax = 0.05 * sig(k)**2

         ierr = nf90_put_att ( ncid, vid, "valid_min",  vmin )
         call check_ncerr(ierr,"Error setting valid min attribute")
         ierr = nf90_put_att ( ncid, vid, "valid_max", vmax )
         call check_ncerr(ierr,"Error setting valid max attribute")

!        Use 1.01 factor to ensure that xmax and xmin fit within range despite
!        any roundoff.

         scalef = 1.01 * (xmax - xmin) / float(vmax - vmin)
         addoff = 0.5*(xmin+xmax)
         ierr  = nf90_put_att ( ncid, vid, "add_offset", addoff )
         call check_ncerr(ierr)
         ierr  = nf90_put_att ( ncid, vid, "scale_factor", scalef)
         call check_ncerr(ierr)

      else if ( hbytes == 4 ) then
         ierr = nf90_def_var ( ncid, vname, NF90_FLOAT, &
                             (/ lonid, latid, monid, yrid /), vid )
         call check_ncerr(ierr,"Error creating variable "// vname)
      else
         print*, " Error, impossible value for hbytes "
         stop
      end if

      ierr = nf90_put_att ( ncid, vid, "long_name", longname )
      call check_ncerr(ierr)
      if ( len_trim(units) /= 0 ) then
         ierr = nf90_put_att ( ncid, vid, "units", units )
         call check_ncerr(ierr)
      end if

      if ( hbytes == 2 ) then
         ierr = nf90_put_att ( ncid, vid, "missing_value", NF90_FILL_SHORT )
         call check_ncerr(ierr,"Error with missing value attribute")
      else
         ierr = nf90_put_att ( ncid, vid, "missing_value", missing_value )
         call check_ncerr(ierr,"Error with missing value attribute")
      end if

!     Leave define mode
      ierr = nf90_enddef ( ncid )
      call check_ncerr(ierr)

!     Turn off the data filling to save time.
      ierr = nf90_set_fill ( ncid, NF90_NOFILL, old_mode)
      call check_ncerr(ierr)

!     Write the months
      ierr = nf90_put_var ( ncid, monid, (/ (m,m=1,12) /) )
      call check_ncerr(ierr,"Error writing months")
      ierr = nf90_put_var ( ncid, latid, ylat )
      call check_ncerr(ierr,"Error writing latitudes")
      ierr = nf90_put_var ( ncid, lonid, xlon )
      call check_ncerr(ierr,"Error writing longitudes")

      ierr = nf90_close ( ncid )
      call check_ncerr(ierr)

   end subroutine create_oldfile

!-----------------------------------------------------------------------------
   subroutine closehist
      
      integer :: ierr, ifile
     
      do ifile = 1,nhfiles
!        The hist_oave files are closed individually in writehist.
         if ( ihtype(ifile) /= hist_oave ) then
            ierr = nf90_close ( histid(ifile) )
            call check_ncerr(ierr,"Error closing history file")
         end if
      end do
   end subroutine closehist
!-------------------------------------------------------------------

   subroutine savehist2D ( name, array, jlat )
      character(len=*), intent(in) :: name
      real, dimension(:) :: array
      integer, intent(in) :: jlat
      real, dimension(size(array),1,1) :: temp

!     Copy to the 3D array that savehist_work expects. The extra level of
!     subroutines is necessary to avoid multiple versions.
!     This copying is extra work if the variable isn't actually saved.
!     Check name first?
      temp(:,1,1) = array
      call savehist_work  ( name, temp, jlat, jlat )

   end subroutine savehist2D
   
   subroutine gsavehist2D ( name, array )
      character(len=*), intent(in) :: name
      real, dimension(:,:) :: array
      real, dimension(size(array,1),size(array,2),1) :: temp

!     Copy to the 3D array that savehist_work expects. The extra level of
!     subroutines is necessary to avoid multiple versions.
      temp(:,:,1) = array
      call savehist_work  ( name, temp, 1, size(array,2) )

   end subroutine gsavehist2D
   
   subroutine savehist3D( name, array, jlat )
      character(len=*), intent(in) :: name
      real, dimension(:,:) :: array
      integer, intent(in) :: jlat
      real, dimension(size(array,1),1,size(array,2)) :: temp

      temp(:,1,:) = array
      call savehist_work  ( name, temp, jlat, jlat )
   end subroutine savehist3D
   
   subroutine gsavehist3D( name, array )
      character(len=*), intent(in) :: name
      real, dimension(:,:,:) :: array
      call savehist_work  ( name, array, 1, size(array,2) )
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
            case ( hist_ave, hist_oave ) 
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
            case ( hist_ave, hist_oave ) 
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
      real, dimension ( nxhis, nyhis ) :: htemp
      logical :: doinc
      
      if ( require_interp .and. .not. present(interp) ) then
         print*, " Error, interp argument required for writehist "
         stop
      end if

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

            if ( ifile == 1 .and. ihtype(ifile) == hist_oave ) then
               if ( .not. ( present(year) .and. present(month) ) ) then
                  print*, "Error year and month arguments to writehist are "
                  print*, "required for old format files"
                  stop
               end if
!              Add a check. Interpolation not supported with old format.
               call oldwrite(nxhis,nyhis,year,month)
               cycle
            end if

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
            iend = istart+nlev-1
            if ( ave_type == hist_ave .or. ave_type == hist_oave ) then
               if ( hist_debug >= 4 ) then
                  print*, "Raw history at point ", histinfo(ifld)%name,&
                    histarray(ihdb,jhdb,istart+khdb-1), count
               end if
               histarray(:,:,istart:iend) =   &
                 histarray(:,:,istart:iend) / max( count, 1 )
            end if
            if ( histinfo(ifld)%output_scale /= 0 ) then
               histarray(:,:,istart:iend) = histarray(:,:,istart:iend) * &
                 histinfo(ifld)%output_scale
            end if
            if ( hist_debug >= 4 ) then
               print*, "History written at point ", histinfo(ifld)%name,&
                    histarray(ihdb,jhdb,istart+khdb-1)
            end if

            call gatherwrap(histarray(:,:,istart:iend),hist_a)
            
            if ( myid == 0 ) then

               do ip = 0,pnproc-1   
                  do n = 0,pnpan-1
                     hist_g(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g,1:iend-istart+1) = &
                        hist_a(1:pil,1+n*pjl:(n+1)*pjl,1:iend-istart+1,ip+1)
                  end do
               end do
                
!              Even multilevel variables are written one level at a time
               do k=istart, iend

                  if ( count == 0 ) then
                     if ( hbytes(ifile) == 2 ) then
                        htemp = NF90_FILL_SHORT
                     else
                        htemp = NF90_FILL_FLOAT
                     end if
                  else

                     if ( present(interp) ) then
                        call interp ( hist_g(:,:,k+1-istart), htemp, histinfo(ifld)%int_type )
                     else
                        htemp = hist_g(:,:,k+1-istart)
                     end if

                     if ( hbytes(ifile) == 2 ) then
                        addoff = histinfo(ifld)%addoff(ifile)
                        sf = histinfo(ifld)%scalef(ifile)
                        umin = sf * vmin + addoff
                        umax = sf * vmax + addoff
                        where ( fpequal(htemp, missing_value) )
                           htemp = NF90_FILL_SHORT
                        elsewhere
!                       Put the scaled array back in the original and let
!                       netcdf take care of conversion to int2
                           htemp = nint((max(umin,min(umax,htemp))-addoff)/sf)
                        endwhere
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

!           Zero ready for next set
            histarray(:,:,istart:iend) = initval(ave_type)
!           Reset the count variable
            histinfo(ifld)%count(ifile) = 0

         end do ! Loop over fields

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
         do ilev=1,nlev
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

            istart = histinfo(ifld)%ptr(ifile) + ilev-1

            if ( ave_type == hist_ave .or. ave_type == hist_oave ) then
               histarray(:,:,istart) = histarray(:,:,istart) / count
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
               ierr = nf90_put_var ( ncid, vid, histarray(:,:,istart), &
                                     start=astart, count=acount )
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
      case ( hist_ave, hist_oave )
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

   subroutine gatherwrap(array_in,array_out)
#ifndef usenc3
      use mpi
#else
      include 'mpif.h'   
#endif
      real, dimension(:,:,:), intent(in) :: array_in
      real, dimension(:,:,:,:), intent(out) :: array_out
      real, dimension(size(array_out,1),size(array_out,2),size(array_in,3),size(array_out,4)) :: array_temp
      integer :: lsize, ierr
      
      lsize = size(array_in,1)*size(array_in,2)*size(array_in,3)
      call MPI_Gather(array_in,lsize,MPI_REAL,array_temp,lsize,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      array_out(:,:,1:size(array_in,3),:) = array_temp(:,:,:,:)
      
   end subroutine gatherwrap

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
