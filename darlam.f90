module darlam
   ! DARLAM module for cc2hist
   implicit none
   private
   logical, public :: darlam_grid = .false. 
   integer, public :: il_darlam, jl_darlam

   ! DARLAM grid parameters
   ! Make these private so ds doesn't interfere with the CCAM ds
   real :: ds, du, tanl, rnml, stl1, stl2
   ! Shared between lconset and lconll
   real :: sign, rotate, refj, rj, rmult, rk, pole
   real, parameter :: pi=3.1415926536, rearth=6371.22e3, cfact=1.45444e-4

   public :: read_darlam_namelist, darlam_atts, lconset, lconij, lconll
contains

   subroutine read_darlam_namelist(un)
      integer, intent(in) :: un
      integer :: il, jl
      namelist /darlam_grid/  il, jl, ds, du, tanl, rnml, stl1, stl2
      ! Set huge default values so that anything that's not set will cause
      ! an obvious problem.
      il = -huge(1)
      jl = -huge(1)
      ds = -huge(1.)
      du = -huge(1.)
      tanl = -huge(1.)
      rnml = -huge(1.)
      stl1 = -huge(1.)
      stl2 = -huge(1.)
      read(un,darlam_grid)
      il_darlam = il
      jl_darlam = jl
   end subroutine read_darlam_namelist

   function darlam_atts() result ( res )
      use history, only : hist_att
      use netcdf, only : nf90_real
      type(hist_att), dimension(6) :: res
      res = (/ hist_att("darlam_ds", NF90_REAL, ds, 0, ""),     &
               hist_att("darlam_du", NF90_REAL, du, 0, ""),     &
               hist_att("darlam_tanl", NF90_REAL, tanl, 0, ""),     &
               hist_att("darlam_rnml", NF90_REAL, rnml, 0, ""),     &
               hist_att("darlam_stl1", NF90_REAL, stl1, 0, ""),     &
               hist_att("darlam_stl2", NF90_REAL, stl2, 0, "") /)
   end function darlam_atts
   
   subroutine lconset
!     rotate added to do S. Pole properly Wed  11-30-1994
!     this is because usually increasing j goes northwards, but not so for
!     S. Pole polar stereographic
!     rewritten completely by jlm August 1993 & May 1995 for SH & NH  (lamb.for)
!     it is cleaned up, & expects -ve latitudes for S. Hem.
!     tanl with proper sign is latitude where j=1
!     du is i value along standard longitude rnml
!     all colatitudes from pi/2, so need sign=-1. for SH

!     this subroutine converts the grid coordinates ri,rj of a point to its
!     latitude (rlat) and longitude (rlong) or vice versa.
!     Lambert conformal projection is used.
!     stl1,stl2  standard latitudes   -ve for S. Hem.
!     ds         length of grid unit(m) at standard latitudes
!     rearth     radius of earth(m)

      real :: chi1, chi2, chiref2
!     e.g. gu=137.057, du=19.76, tanl=-60.2, stl1=-40., stl2=-10., rnml=130.
!jlm  use chi for colatitudes in radians
      pole=.5*pi
      sign=1.
      if(.5*(stl1+stl2).lt.0.)sign=-1.       ! decides which hemisphere
      rotate=1.
      if(stl1.lt.0..and.stl1.eq.stl2)rotate=-1.  ! for grid over S. pole
      chi1=pole-sign*stl1*pi/180.
      chi2=pole-sign*stl2*pi/180.
      if(stl1.eq.stl2)then
        rk=1.
        rmult=(1.+cos(chi1))  *rearth/ds
      else
        rk=log(sin(chi1)/sin(chi2))/log(tan(.5*chi1)/tan(.5*chi2))    ! i.e. K
        rmult=sin(chi1)*tan(.5*chi1)**(-rk)  *rearth/(ds*rk)
      endif
      chiref2=.5*(pole-sign*tanl*pi/180.)
      refj=1.+rotate*sign*rmult*tan(chiref2)**rk
   end subroutine lconset

   subroutine lconij(rlong,rlat,ri,rj,theta)
      real, intent(in) :: rlong, rlat
      real, intent(out) :: ri, rj, theta
      real :: ther, chion2, r
      theta=rk*(rlong-rnml)
      ther=theta*pi/180.
      theta=-sign*theta
!     theta gives the angle needed to rotate the lat-long coordinates
!        anticlockwise to align with the l_c coordinates
!     For usual Australian grids, theta is +ve east of longitude rnml
!     N.B. theta is used in setup routines to get u,v from analyses as:
!              u = (wmer*sin(theta) + wzon*cos(theta))
!              v = (wmer*cos(theta) - wzon*sin(theta))
!     so want theta ~ -(rlong-rnml) in NH; and theta ~ (rlong-rnml) in SH
      chion2=.5*(pole-sign*rlat*pi/180.)
!     print *,'chion2,tan(chion2) ',chion2,tan(chion2)
      r=rmult*tan(chion2)**rk
      rj=refj-sign*r*cos(ther)*rotate   ! Wed  11-30-1994
!     rj=refj-sign*r*cos(ther)
      ri=du+r*sin(ther)*rotate          ! Wed  11-30-1994
   end subroutine lconij
   
   subroutine lconll(rlong,rlat,ri,rj)
      real, intent(in) :: ri, rj
      real, intent(out) :: rlong, rlat
      real :: alpha, beta, chion2, ther
      alpha=ri-du
      beta=rotate*sign*(refj-rj)
!     print *,'alpha,beta',alpha,beta
      if(alpha.ne.0..or.beta.ne.0.)then
         ther=atan2( alpha,beta)
      else
         print *,'At pole: alpha, beta in lconll =',alpha,beta
         stop
      endif
      rlong=rnml+rotate*ther*180./(pi*rk)
      if (rlong > 360.) rlong = rlong-360.
      if ( rk==1. ) then
         chion2 = sign*atan(sqrt(alpha**2+beta**2)/rmult)
      else
         chion2 = sign*atan( (beta/(rmult*cos(ther)))**(1./rk) )
      endif
      rlat = (sign*pole-2.*chion2)*180./pi  ! Sun  05-07-1995
   end subroutine lconll

end module darlam
