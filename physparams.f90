module physparams
  ! F90 version of model PHYSPARAMS.f include file
   use precis_m, only : rx

   real, parameter :: stefbo = 5.67e-8  ! Stefan-Boltzmann constant
   real(kind=rx), parameter :: erad = 6.37122e6, eradsq = erad*erad !Radius of earth
   real, parameter :: cp = 1004.64  ! Specific heat of dry air at const P
   real, parameter :: rdry = 287.04 ! Specific gas const for dry air
   real, parameter :: epsil = 0.622 ! Ratio molec wt of H2O vapour to dry air
   real, parameter :: rvap = 461.   ! Gas constant for water vapour
   real, parameter :: hlf = 3.35e5  ! Latent heat of fusion :: at 0 deg C
   real, parameter :: hl = 2.5e6    !Latent heat of vaporization :: at 0 deg. C
   real, parameter :: hls = hl+hlf  !  "      "   " sublimation
   real, parameter :: hlcp = hl/cp, hlfcp = hlf/cp, hlscp = hlcp+hlfcp
   real, parameter :: grav = 9.80616 ! Acceleration of gravity
   
   real, parameter :: sq2 = 1.414213562373092 ! Square root of 2
   real, parameter :: cappa = rdry/cp 
   real, parameter :: tomg = 2*7.2921233e-5 ! 2*omega
   real, parameter :: pi = 3.14159265358979 !Good ol' pi
   real, parameter :: degtorad = pi/180.

end module physparams
