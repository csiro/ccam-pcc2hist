module vertutils_m

!  Routines involving the vertical coordinate

   implicit none
   private
   public :: setsig, sig2sigh, sig2ds, setsig_lbl, initial, height, &
             initial_full, check_plevs, height_std

   logical, public :: vertutils_debug = .false.
   integer, public :: vertutils_idebug=1, vertutils_jdebug=1

contains

   subroutine setsig ( sig, sigh, n )
      real, intent(out), dimension(:) :: sig, sigh
      integer, intent(in) :: n
      integer :: k
!     Calculate sigma levels using CSIRO model formula

!     Set half levels
      do k=1,n+1
         sigh(k) = (n+1-k)**2*(n-2+2.*k)/n**3
      end do
!     Set full levels
      do k=1,n
         sig(k) = 0.5*(sigh(k)+sigh(k+1))
      end do
 
   end subroutine setsig
   
   subroutine sig2sigh(sig, sigh)

!     Calculate half level values of sigma from the full levels.
!     For the CSIRO model the level is at the mid-point of the layer. 

      real, intent(in),  dimension(:) :: sig
      real, intent(out), dimension(:) :: sigh
      integer :: k

      sigh(1) = 1.0
      do k=1,size(sig)-1
         sigh(k+1) = 2*sig(k) - sigh(k)
      end do
      sigh(size(sig)+1) = 0.0

   end subroutine sig2sigh

   subroutine sig2ds(sig, ds)

!     Calculate layer thickness from the full level sigma values.
!     For the CSIRO model the level is at the mid-point of the layer. 

      real, intent(in),  dimension(:) :: sig
      real, intent(out), dimension(:) :: ds
      real, dimension(size(sig)+1) :: sigh

      call sig2sigh(sig,sigh)
      ds = sigh(1:size(ds)) - sigh(2:)

   end subroutine sig2ds

   subroutine setsig_lbl ( sig, sigh, n )
      real, intent(out), dimension(:) :: sig, sigh
      integer, intent(in) :: n
      integer :: k
      real, dimension(n+1) :: plev
!     Calculate sigma levels as used in the LBL calculations described in 
!     Schwarzkopf and Fels 1991. These are equally spaced in pressure below 
!     100 mb and log spaced above that. The surface pressure is assumed to
!     be 1013 and there are flux levels at both 1013 and 1000.
!     Use n=122 to match their levels with a top at 0.001 mb.

!     Equally spaced from 1 to 0.1
      plev(1) = 1013.
      do k=2,min(n+1,47)
         plev(k) = 1000 - 20.0*(k-2)
      end do

!     Log spaced above this with 15 levels per decade
      do k=48,n
         plev(k) = exp ( log(100.) - (k-47)*log(10.0)/15.0 )
      end do

      plev(n+1) = 0.0

      print*, "PLEV", plev
      sigh = plev / plev(1)

!     Set full levels
      do k=1,n
         sig(k) = 0.5*(sigh(k)+sigh(k+1))
      end do
 
   end subroutine setsig_lbl

subroutine initial(nl,sig,sigh,am,algf,alogg)

!  Calculate the A matrix needed to integrate hydrostatic equation

   use physparams
   implicit none
   integer, intent(in) :: nl
   real, intent(in), dimension(nl) :: sig
   real, intent(in), dimension(nl+1) :: sigh
   real, intent(out), dimension(nl,nl) :: am
   real, intent(out), dimension(nl) :: algf, alogg

   real, dimension(nl,nl) :: amh
   real, dimension(nl) :: algh, atm, atp, btm, btp
   integer :: i, j, k

!  Sigma level arrays
   do k=1,nl
      algf(k)=alog(sig(k))
      algh(k)=alog(sigh(k))
   end do
   do k=2,nl
      alogg(k) = 1.0 / (algf(k) - algf(k-1))
   end do

   do k=1,nl-1
      btm(k)=1.0/(algf(k)-algf(k+1))
      btp(k)=-btm(k)
      atm(k)=-algf(k+1)*btm(k)
      atp(k)=algf(k)*btm(k)
   end do

!****  temp/geopotl coupling for phi=am*t+phist
!***  based upon t=a+b.log(sigma) : energy conserving -
!****  changes made to dynamical terms (in DYNMNL).
   am = 0.
   amh = 0.

!  Define temporary full level geopotential matrix values
   am(1,1)=rdry*(atm(1)*(algh(1)-algf(1))+btm(1)*(algh(1)**2-algf(1)**2)*0.5)
   am(1,2)=rdry*(atp(1)*(algh(1)-algf(1))+btp(1)*(algh(1)**2-algf(1)**2)*0.5)
   am(2,1)=rdry*(atm(1)*(algh(1)-algf(2))+btm(1)*(algh(1)**2-algf(2)**2)*0.5)
   am(2,2)=rdry*(atp(1)*(algh(1)-algf(2))+btp(1)*(algh(1)**2-algf(2)**2)*0.5)
   do i=3,nl
      do j=1,i
         if ( j == i-1 ) then
            am(i,j) = am(i-1,j)    &
              + rdry*(atm(i-1)*(algf(i-1)-algf(i))+btm(i-1)*(algf(i-1)**2  &
              - algf(i)**2)*0.5)
         else if ( j == i) then
            am(i,j) = am(i-1,j)    &
              + rdry*(atp(i-1)*(algf(i-1)-algf(i))+btp(i-1)*(algf(i-1)**2  &
              - algf(i)**2)*0.5)
         else 
            am(i,j) = am(i-1,j)
         end if
      end do
   end do

!  The model redefines the full level am from a half level matrix. This
!  is needed for energy conservation. However for diagnostic purposes it's
!  more accurate to skip this part.

end subroutine initial

subroutine initial_full(nl,sig,sigh,am,amh,algf,alogg)

!  Calculate the A matrix needed to integrate hydrostatic equation
!  Exactly as in model, calculate AM and AMH

   use physparams
   implicit none
   integer, intent(in) :: nl
   real, intent(in), dimension(nl) :: sig
   real, intent(in), dimension(nl+1) :: sigh
   real, intent(out), dimension(nl,nl) :: am, amh
   real, intent(out), dimension(nl) :: algf, alogg

   real, dimension(nl) :: algh, atm, atp, btm, btp
   integer :: i, j, k

   call initial(nl,sig,sigh,am,algf,alogg)

!  Sigma level arrays
   do k=1,nl
      algh(k)=alog(sigh(k))
   end do
   do k=1,nl-1
      btm(k)=1.0/(algf(k)-algf(k+1))
      btp(k)=-btm(k)
      atm(k)=-algf(k+1)*btm(k)
      atp(k)=algf(k)*btm(k)
   end do

!  Create half level geopotential matrix values
   amh(1,1) = rdry*(atm(1)*(algh(1)-algh(2))+btm(1)*(algh(1)**2-algh(2)**2)*0.5)
   amh(1,2)=rdry*(atp(1)*(algh(1)-algh(2))+btp(1)*(algh(1)**2-algh(2)**2)*0.5)
   do i=2,nl-1
      amh(i,:)=am(i,:)
      ! j=i
      amh(i,i) = amh(i,i) + &
         rdry*(atm(i)*(algf(i)-algh(i+1))+btm(i)*(algf(i)**2-algh(i+1)**2)*0.5)
      ! j=i+1
      amh(i,i+1) = amh(i,i+1) + &
         rdry*(atp(i)*(algf(i)-algh(i+1))+btp(i)*(algf(i)**2-algh(i+1)**2)*0.5)
   end do
   amh(nl,:)=2.0*am(nl,:)-amh(nl-1,:)
!  Create final full level geopotential matrix values
   am(1,:) = 0.5*amh(1,:)
   do i=2,nl
      am(i,:)=0.5*(amh(i-1,:)+amh(i,:))
   end do

end subroutine initial_full


!-----------------------------------------------------------------
subroutine height(lon,lat2,nl,tg,qg,zg,pg,am,phi)

!  Calculate geopotential height on sigma levels using virtual temperature
   use physparams
   implicit none
   integer, intent(in) :: lon, lat2, nl
   real, intent(in), dimension(lon,lat2,nl) :: tg, qg
   real, intent(in), dimension(lon,lat2)    :: zg, pg
   real, intent(in), dimension(nl,nl) :: am
   real, intent(out), dimension(lon,lat2,nl) :: phi

   real, dimension(lon,nl) :: tv
   integer :: j, k, lg 

   do lg=1,lat2

!     Calculate virtual temperature.
      tv = tg(:,lg,:) * (epsil+qg(:,lg,:))/(epsil*(1.+qg(:,lg,:)))

!     Calculate height on sigma levels
      do k=1,nl
         phi(:,lg,k) = zg(:,lg) * grav
         do j=1,nl
            phi(:,lg,k) = phi(:,lg,k)+am(k,j)*tv(:,j)
         end do
      end do
   end do       ! Latitude loop

   phi = phi / grav

   return
end subroutine height

subroutine height_std(use_hyb,lon,lat,nl,tg,qg,zg,pg,am,sigin,phistd,pstd,nstd)

!  Calculate geopotential height on standard pressure levels using 
!  virtual temperature
   use physparams
   use hyblevs_m
   implicit none
   logical, intent(in) :: use_hyb
   integer, intent(in) :: lon, lat, nl
   real, intent(in), dimension(lon,lat,nl) :: tg, qg
   real, intent(in), dimension(lon,lat)    :: zg, pg
   real, intent(in), dimension(nl,nl) :: am
   real, intent(in), dimension(nl) :: sigin
   integer, intent(in) :: nstd
   real, intent(in), dimension(nstd) :: pstd
   real, intent(out), dimension(lon,lat,nstd) :: phistd

   real, dimension(lon,nl) :: tv, phi
   integer :: j, k, lg, mg, kstd, kk, jlev
   real :: plev, plog, temp, qemp1, qemp2
   real, dimension(lon,nl) :: muf, prf, fff1

   real, parameter :: lapse = 6.5e-3
!  Extrapolation factor for temperature. Set to zero for no extrapolation.
   real, parameter :: exfac = rdry * lapse / grav 
!  Crtical lapse rate to decide where to extrapolate from
   real, parameter :: critfac = rdry * 1.e-3 / grav 

   do lg=1,lat

      if ( use_hyb ) then
         do k=1,nl
            ! muf = dp/dn
            muf(:,k) = dadnf(k) + dbdnf(k)*pg(:,lg)
            prf(:,k) = anf(k)+(bnf(k)*pg(:,lg))
            fff1(:,k) = muf(:,k)/prf(:,k)*sig(k)
         end do
      else
         do k=1,nl
            prf(:,k) = pg(:,lg)*sigin(k)
         end do
         fff1 = 1.0
      end if

!     Calculate virtual temperature.
      tv = tg(:,lg,:) * (epsil+qg(:,lg,:))/(epsil*(1.+qg(:,lg,:)))

!     Calculate height on model levels
      do k=1,nl
         phi(:,k) = zg(:,lg) * grav
         do j=1,nl
            phi(:,k) = phi(:,k)+am(k,j)*tv(:,j)*fff1(:,j)
         end do
      end do

      if ( vertutils_debug .and. lg==vertutils_jdebug ) then
         print*, "Temperature", tg(vertutils_idebug,lg,:)
         print*, "Virtual temperature", tv(vertutils_idebug,:)
         print*, "Surface pressure", pg(vertutils_idebug,lg)
         print*, "Surface height", zg(vertutils_idebug,lg)
         print*, "Height on model levels", phi(vertutils_idebug,:)
      end if

!     Now calculate the height at the standard pressure levels.
!     Everything from here done with pressure so it works for both
!     sigma and hybrid levels.
      do kstd = 1,nstd
         do mg=1,lon
            plev = pstd(kstd)
            plog = alog(plev)
            if( plev < prf(mg,nl) ) then
               jlev=nl
               temp=tv(mg,nl)
            else if ( plev .gt. prf(mg,1) ) then
               jlev=1
               ! Extrapolate temperature using standard lapse rate
               ! Avoid surface inversions by starting from the lowest level
               ! with a positive lapse rate. This is important to get a smooth
               ! 1000 hPa height over the Himalayas. Start from second level
               ! to avoid extreme temperature gradients near the surface
               do k=2,4 ! Level 5 is a reasonable BL top
                  ! Extrapolate upwards using critical lapse rate
                  temp=tv(mg,k)*(prf(mg,k+1)/prf(mg,k))**critfac
                  if ( tv(mg,k+1) < temp ) exit
               end do
               temp=tv(mg,k)*(plev/prf(mg,k))**exfac
               temp = 0.5 * ( temp + tv(mg,k) )
               if ( vertutils_debug .and. lg==vertutils_jdebug .and. mg==vertutils_idebug) then
                  print*, "Extrapolating temp from level", k
               end if
            else
               do kk=1,nl
                  jlev = kk
!                 Find first level above the standard level
                  if ( prf(mg,kk) <= plev ) exit
               end do
!              Interpolate the temperature in log pressure. Target level is 
!              between jlev-1 and jlev
               qemp1 = ( plog - alog(prf(mg,jlev-1)) )*tv(mg,jlev)
               qemp2 = ( alog(prf(mg,jlev)) - plog )*tv(mg,jlev-1)
               temp = (qemp1+qemp2)/(alog(prf(mg,jlev))-alog(prf(mg,jlev-1)))
!              Use average of interpolated and level temperature
               temp = 0.5 * ( temp + tv(mg,jlev) )
            end if
!
            phistd(mg,lg,kstd) = phi(mg,jlev) - &
                                 rdry*temp*(plog - alog(prf(mg,jlev)))
         end do
      end do    ! Std level loop
   end do       ! Latitude loop

   phistd = phistd / grav

   if ( vertutils_debug ) then
      print*, "Height on pressure levels", phistd(vertutils_idebug,vertutils_jdebug,:)
   end if

end subroutine height_std

   subroutine check_plevs(plevs)
      ! Sanity checks and sorting on the specified output pressure levels.
      real, dimension(:), intent(inout) :: plevs
      integer :: i, j
      real :: temp

      if ( size(plevs) == 0 ) then
         print*, "Error, no pressure levels specified"
         stop
      end if
      if ( minval(plevs) < 0. ) then
         print*, "Error, negative pressure set"
         stop
      end if
      if ( maxval(plevs) > 1100. ) then
         print*, "Error, pressure value is too large"
         stop
      end if
      ! Sort into decreasing order
      do i=1,size(plevs)
!        Using maxval here is just a trick to reduce the array result of 
!        maxloc to a scalar.
         j = maxval ( maxloc(plevs(i:)) ) + i - 1
         temp = plevs(i)
         plevs(i) = plevs(j)
         plevs(j) = temp
      end do
   end subroutine check_plevs

end module vertutils_m
