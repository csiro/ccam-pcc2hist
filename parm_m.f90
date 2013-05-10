module parm_m

! Model control parameters. This module combines parm.h, dates.h and sigs.h
! values set here are defaults that may be overridden by the namelist input.
! Some physical constants are also set.
   use precis_m
   implicit none
   private
   integer, public, save :: ktau = 0
   integer, public, save :: mex = 4
   integer, public, save :: ntau = 48
   integer, public, save :: m = 5
   integer, public, save :: kstart = 1
   integer, public, save :: meso = 6
   integer, public, save :: khdif = 20
   integer, public, save :: nhor = -1
   integer, public, save :: nhort = -1
   integer, public, save :: nhorps = -1
   integer, public, save :: khor = 0
   integer, public, save :: nonl = 0
   integer, public, save :: nxmap = 0
   integer, public, save :: mfix = 0
   real(kind=rx), public, save    :: rlong0 = 0.0_rx
   real(kind=rx), public, save    :: rlat0 = 90.0_rx
   real(kind=rx), public, save    :: schmidt = 1.0_rx
   real(kind=rx), public, save    :: epsp = 0.0_rx
   real(kind=rx), public, parameter :: schm13 = 0.1_rx
   real(kind=rx), public, save    :: hdiff = 0.0_rx
   real(kind=rx), public, save    :: restol =  5.0e-6_rx
   real(kind=rx), public, save    :: tbar = 0.0_rx
   integer, public, save :: norder = 4 ! Order of FD approx. in nonlin
   logical, public, save :: agrid = .true.  ! Use A grid
   integer, public, save :: nrot = 1
   integer, public, save :: nupg = 0  ! Old value, JLM now uses 1
   integer, public, save :: nuv = 0
!  Options for staggering/unstaggering of winds
   integer, public, save :: nstag  = 0
   integer, public, save :: nstagu = 0
   integer, public, save :: mstagpt = 4
   integer, public, save :: ntang = 2   ! For setxyz

   integer, public, save :: ktime
   integer, public, save :: kdate 
   real(kind=rx), public, save    :: dt = 1800
   real(kind=rx), public, save    :: timer
   real(kind=rx), public, save    :: timeg

   real(kind=rx), public, parameter :: rearth = 6.371220e6_rx
   real(kind=rx), public, parameter :: pi = 3.141592653589793238_rx
   real(kind=rx), public, parameter :: rtod=180.0_rx/pi
   real(kind=rx), public, parameter :: dtor=pi/180.0_rx
   real(kind=rx), public, parameter :: grav = 9.806_rx
   real(kind=rx), public, parameter :: omega = 7.292e-5_rx

   real(kind=rx), public, save :: alpha = 0.0_rx ! Wind rotation angle 
                                                 ! (degrees)

   real(kind=rx), public, parameter :: eps = 0.0_rx ! uncentering parameter 
                                                    ! (not used)

   logical, public, save :: newdepts = .false.

!  Just do advection
   logical, public, save :: advect_only = .false.
end module parm_m
