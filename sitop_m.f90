! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module sitop_m
   implicit none
   ! Controls whether vertical extrapolation below surface is done
   integer, parameter :: vextrap_default=0, vextrap_lin=1, vextrap_none=2, &
                         vextrap_missing=3, vextrap_t=4

   ! Within the model, all sitop calls will have the same dimensions
   ! ix = 0 indicates that it hasn't been initialised yet
   integer, private, save :: ix=0, iy, nsglvs, nprlvs, nsgm1
   integer, private, save :: oix=0, oiy, onsglvs, onprlvs
   ! Will be dimensioned (ix,iy)
   integer, private, allocatable, save, dimension(:,:) :: &
            mexdn1, mexdn2, mexup1, mexup2, min1, min2
   ! Will be dimensioned (ix,nprlvs,iy)
   integer,  private, allocatable, save, dimension(:,:,:) :: jaa, ojaa
   real, private, allocatable, save, dimension(:,:,:) :: sig, osig, ow1
   ! (nsglvs)
   real, private, allocatable, save, dimension(:) :: siglvs, h, x1, x2, a, c
   real, private, allocatable, save, dimension(:) :: osiglvs

contains
   subroutine sitop_setup ( sigr, prelvs, press, maxlev, minlev )

      ! Setup for sitop. Calculate things that don't depend on the field
      ! being interpolated.
      use utils_m, only : assert, search_fgt
      real, dimension(:), intent(in)       :: sigr   ! Sigma levels
      real, dimension(:), intent(in)       :: prelvs ! Pressure levels
      real, dimension(:,:), intent(in)     :: press  ! Surface pressure
      real, parameter :: cmu1=0.5, clmdam=0.5

      integer, intent(inout) :: maxlev, minlev
      integer :: i, k, ii, ns, iii, kp, ks
      integer :: j, j1, jj 
      integer :: maxlev_g, minlev_g, ierr
      real :: v1, w1, v2, w2, h1, h2, h3, z

!----------------------------------------------------------------------
!
      if ( ix == 0 ) then
         ! Not initialised yet
         nsglvs = size(sigr)
         nprlvs = size(prelvs)
         ix = size(press,dim=1)
         iy = size(press,dim=2)
         allocate ( jaa(ix,nprlvs,iy), sig(ix,nprlvs,iy) )
         allocate ( mexdn1(ix,iy), mexdn2(ix,iy), mexup1(ix,iy), &
                    mexup2(ix,iy), min1(ix,iy), min2(ix,iy) )
         allocate ( siglvs(nsglvs), h(nsglvs), x1(nsglvs), x2(nsglvs), &
                    a(nsglvs), c(nsglvs) )

         c(1) = cmu1
         a(nsglvs) = clmdam
!        Sigma levels in increasing order
         siglvs = sigr(nsglvs:1:-1)
         nsgm1 = nsglvs-1
         do i = 1,nsgm1
            h(i) = siglvs(i+1)-siglvs(i)
         end do
         do i = 2,nsgm1
            x1(i) = 0.5/(h(i)+h(i-1))
            x2(i) = h(i)/h(i-1)
            a(i) = h(i)*x1(i)
            c(i) = h(i-1)*x1(i)
         end do
         c(nsglvs) = 0.
      end if

!     Reverse data and pressure levels.
      do k = 1,nprlvs
         sig(:,nprlvs+1-k,:) = prelvs(k)/press
      end do

      do j = 1,iy
!
!     mexup1, mexup2    mexdn1,mexdn2    min1,min2 set
!     the upper and lower limits on upward and downward
!     extrapolation and interpolation.
!     if there are no levels in a given class,the limits are
!     set to zero.
!
         do i = 1,ix
            mexup1(i,j) = 0
            mexup2(i,j) = 0
            do ii = 1,nprlvs
               if ( sig(i,ii,j)<siglvs(1) ) mexup2(i,j)=ii
            end do
            if ( mexup2(i,j)>0 ) mexup1(i,j)=1
            mexdn1(i,j)=0
            mexdn2(i,j)=0
            do ii=1,nprlvs
               iii=nprlvs+1-ii
               if ( sig(i,iii,j)>siglvs(nsglvs) ) mexdn1(i,j)=iii
            end do
            if ( mexdn1(i,j)>0 ) mexdn2(i,j)=nprlvs

            if ( mexup2(i,j)==0 ) then
               !             no upward extrapolation
               if ( mexdn1(i,j)==0 ) then
                  min1(i,j)=1
                  min2(i,j)=nprlvs
               else
                  min1(i,j)=1
                  min2(i,j)=mexdn1(i,j)-1
                  if ( mexdn1(i,j)==1 ) min1(i,j)=0
               end if
            else 
!
!         upward extrapolation
!
               if ( mexdn1(i,j)==0 ) then
                  min1(i,j)=mexup2(i,j)+1
                  min2(i,j)=nprlvs
                  if ( mexup2(i,j)>=nprlvs ) min1(i,j)=0
                  if ( mexup2(i,j)>=nprlvs ) min2(i,j)=0
               else
                  min1(i,j)=mexup2(i,j)+1
                  min2(i,j)=mexdn1(i,j)-1
                  if ( mexdn1(i,j)==(mexup2(i,j)+1) ) min1(i,j)=0
                  if ( mexdn1(i,j)==(mexup2(i,j)+1) ) min2(i,j)=0
               end if
            end if
         end do

!        Calculate jaa, the index of the lowest sigma level above the
!        given pressure level.
         do kp = 1,nprlvs 
            jaa(:,kp,j) = search_fgt(siglvs(:nsgm1),sig(:,kp,j),nsgm1,ix)
         end do
         
      end do

      minlev = 1
      maxlev = nsglvs
      
   end subroutine sitop_setup


   subroutine sitop ( grids, gridp, vextrap )


!     This routine interpolates a field "grids" on sigma surfaces
!     to "gridp" on pressure surfaces
!     The interpolation is performed using cubic splines
!     prelvs(nprlvs) are the pressure levels and siglvs(nsglvs) those
!     of sigma
!
!     All features of the interpolating spline are determined except
!     for the impostion of one condition at each end
!     These are prescribed through the quantities
!      "cmu1","clmdam","c1" and "cm"
!
!     For specified slopes at each end - gsd(1) and gsd(nsglvs) - have
!      cmu1=clmdam=0.0 , c1=gsd(1) , cm=gsd(nsglvs)
!
!     For specified second derivative - gsdc(1) and gsdd(nsglvs) - have
!      cmu1=clmdam=0.5 ,
!     c1=1,5*(gs(2     )-gs(1       ))/h(1       )-
!        h(1       )*gsdd(1     )*0.25
!     cm=1.5*(gs(nsglvs)-gs(nsglvs-1))/h(nsglvs-1)+
!        h(nsglvs-1)*gsdd(nsglvs)*0.25
!
!     Note the case gsdd( )=gsdd(nsglvs)=0.0 is a particular case
!     if the above and is refered to as the "natural" spline
!
!     The theory upon which the routine is based may be found in
!     Ahlberg,J.A.,E.N.Nilson and J.L.Wals",1967 : The Theory of Spline
!     and Their Applications. New York , Academic Press , 284 pp.
!     (pages 9-15)
!
!     Note that this uses the natural spline, with the extrapolation
!     imposed separately. This means that the derivatives are not continuous
!     across the boundary between interpolation and extrapolation.

#ifdef usenc_mod
      use netcdf, only : NF90_FILL_FLOAT
#else
      use netcdf_m, only : NF90_FILL_FLOAT
#endif      
      use utils_m, only : assert, search_fgt
      use physparams
      ! Lapse rate used in temperature extrapolation
      real, parameter :: extrap_lapse = 6.5e-3
      real, parameter :: konst = extrap_lapse*rdry/grav
      real, dimension(:,:,:), intent(in)   :: grids  ! Sigma field
      real, dimension(:,:,:), intent(out)  :: gridp  ! Pressure field
      integer, intent(in) :: vextrap       ! Controls vertical extrapolation

      real, dimension(size(grids,dim=1),size(grids,dim=3)) :: d, gd, gs
      real,    dimension(size(gridp,dim=1),size(gridp,dim=3)) :: gp
      real, dimension(size(grids,dim=1)) :: x3
      real, parameter :: cmu1=0.5, clmdam=0.5
      integer :: i, k, ii, ns, iii, kp, ks
      integer :: minmin1, maxmin1, maxmin2, j, j1, jj 
      real :: v1, w1, v2, w2, h1, h2, h3, z

!----------------------------------------------------------------------
!
      call assert ( ix /= 0, "Error, sitop called before sitop_setup")
!     Check all the array sizes match
      call assert ( nsglvs == size ( grids, dim=3 ), &
                  "Error, number of sigma levels doesn't match in sitop" )
      call assert ( nprlvs == size ( gridp, dim=3 ), &
                  "Error, number of pressure levels doesn't match in sitop" )

      do j=1,iy
     
         do ns=1,nsglvs
            gs(:,nsglvs+1-ns) = grids(:,j,ns)
         end do
         d(:,1) = 1.5*(gs(:,2)-gs(:,1)) / h(1)
         d(:,nsglvs) = 1.5*(gs(:,nsglvs)-gs(:,nsglvs-1)) / h(nsglvs-1)
         do ns=2,nsgm1
            d(:,ns) = 3.0*x1(ns)*( (gs(:,ns)-gs(:,ns-1))*x2(ns) + &
                                   (gs(:,ns+1)-gs(:,ns))/x2(ns) )
         end do
!
!     calculate spline derivatives at orginal points
!
         gd(:,1) = -c(1)
         do ns = 2,nsglvs
            x3 = 1. / (1.+a(ns)*gd(:,ns-1))
            gd(:,ns) = -c(ns)*x3
            d(:,ns) = (d(:,ns)-a(ns)*d(:,ns-1)) * x3
         end do
         gd(:,nsglvs) = d(:,nsglvs)
         do ns = 1,nsgm1
            k = nsglvs-ns
            gd(:,k) = gd(:,k)*gd(:,k+1)+d(:,k)
         end do
!
!     Find the extreme values of the interpolation. The assumption
!     used here is that most of the values will be interpolated and there
!     will be just a few odd extrapolations. Therefore it is worth 
!     vectorising the interpolation section, even if a few values are
!     calculated unnecessarily. These will be fixed up in the extrapolation
!     section
         minmin1 = minval(min1(:,j))
         maxmin1 = maxval(min1(:,j))
         maxmin2 = maxval(min2(:,j))
!
!     Do interpolation
!
         if ( maxmin1/=0 ) then

            do ii = max(1,minmin1),maxmin2
               do i = 1,ix
!              Use max to ensure stay in array bounds. Any points where this
!              is used will be fixed in the extrapolation section.
                  j1=max(jaa(i,ii,j)-1,1)
                  jj=jaa(i,ii,j)
                  v1=siglvs(jaa(i,ii,j))-sig(i,ii,j)
                  w1=sig(i,ii,j)-siglvs(j1)
                  v2=v1*v1
                  w2=w1*w1
                  h1=h(j1)
                  h2=1.0/(h1*h1)
                  h3=h2/h1
                  z=(gd(i,j1)*v2*w1-gd(i,jaa(i,ii,j))*w2*v1)*h2
                  gp(i,ii) = z + (gs(i,j1)*v2*(w1+w1+h1)+  &
                                gs(i,jaa(i,ii,j))*w2*(v1+v1+h1))*h3
               end do
            end do
            
         end if
!
!     linear extrapolation
!
!     here we use the options of extrapolation as explained above
!
!
         select case ( vextrap )
         case ( vextrap_lin ) 
            do i = 1,ix
               if ( mexup1(i,j)/=0 ) then
                  do ii = 1,mexup2(i,j)
                     gp(i,ii) = gs(i,1)+gd(i,1)*(sig(i,ii,j)-siglvs(1))
                  end do
               end if
               if ( mexdn1(i,j)/=0 ) then
                  do ii = mexdn1(i,j),nprlvs
                     gp(i,ii) = gs(i,nsglvs)+gd(i,nsglvs)*(sig(i,ii,j)-siglvs(nsglvs))
                  end do
               end if
            end do
         case ( vextrap_none )
            do i = 1,ix
               do ii = 1,mexup2(i,j)
                  gp(i,ii) = gs(i,1)
               end do
               if ( mexdn1(i,j)/=0 ) then
                  do ii = mexdn1(i,j),nprlvs
                     gp(i,ii) = gs(i,nsglvs)
                  end do
               end if
            end do
         case ( vextrap_missing )
            do i = 1,ix
               do ii = 1,mexup2(i,j)
                  gp(i,ii) = NF90_FILL_FLOAT ! local missing value 
               end do
               if ( mexdn1(i,j)/=0 ) then
                  do ii = mexdn1(i,j),nprlvs
                     gp(i,ii) = NF90_FILL_FLOAT ! local missing value
                  end do
               end if
            end do
         case ( vextrap_t )
            ! Extrapolate using the specified lapse rate
            ! temp = ts[0,0] * sig ** k
            do i = 1,ix
               ! Linear extrapolation at the top.
               do ii = 1,mexup2(i,j)
                  gp(i,ii) = gs(i,1)+gd(i,1)*(sig(i,ii,j)-siglvs(1))
               end do
               if ( mexdn1(i,j)/=0 ) then
                  do ii = mexdn1(i,j),nprlvs
                     gp(i,ii) = gs(i,nsglvs)*(sig(i,ii,j)/siglvs(nsglvs))**konst
                  end do
               end if
            end do
         case default
            print*, "Error in sitop, unexpected value of vextrap", vextrap
            stop
         end select

         do ii = 1,nprlvs
            gridp(:,j,nprlvs+1-ii) = gp(:,ii)
         end do

      end do

   end subroutine sitop

   
   subroutine mitop_setup ( sigr, mtrlvs, height, zs, tg, qg, maxlev, minlev )

      ! Setup for mitop. Calculate things that don't depend on the field
      ! being interpolated.
      use physparams
      use utils_m, only : assert, search_fgt
      real, dimension(:), intent(in)       :: sigr      ! Sigma levels
      real, dimension(:), intent(in)       :: mtrlvs    ! Height levels
      real, dimension(:,:,:), intent(in)   :: height    ! Output levels
      real, dimension(:,:,:), intent(in)   :: tg        ! Temperature
      real, dimension(:,:,:), intent(in)   :: qg        ! Water vapor
      real, dimension(:,:), intent(in)     :: zs        ! surface height
      real, dimension(size(tg,dim=1),size(tg,dim=3)) :: tv
      real, parameter :: cmu1=0.5, clmdam=0.5

      integer, intent(inout) :: maxlev, minlev
      integer :: i, k, ii, ns, iii, kp, ks
      integer :: j, j1, jj 
      integer :: maxlev_g, minlev_g, ierr
      real :: v1, w1, v2, w2, h1, h2, h3, z
      real :: mtrphys

!----------------------------------------------------------------------
!
      if ( ix == 0 ) then
         ! Not initialised yet
         nsglvs = size(sigr)
         nprlvs = size(mtrlvs)
         ix = size(height,dim=1)
         iy = size(height,dim=2)
         allocate ( jaa(ix,nprlvs,iy), sig(ix,nprlvs,iy) )
         allocate ( mexdn1(ix,iy), mexdn2(ix,iy), mexup1(ix,iy), &
                    mexup2(ix,iy), min1(ix,iy), min2(ix,iy) )
         allocate ( siglvs(nsglvs), h(nsglvs), x1(nsglvs), x2(nsglvs), &
                    a(nsglvs), c(nsglvs) )

         c(1) = cmu1
         a(nsglvs) = clmdam
!        Sigma levels in increasing order
         siglvs = sigr(nsglvs:1:-1)
         nsgm1 = nsglvs-1
         do i = 1,nsgm1
            h(i) = siglvs(i+1)-siglvs(i)
         end do
         do i = 2,nsgm1
            x1(i) = 0.5/(h(i)+h(i-1))
            x2(i) = h(i)/h(i-1)
            a(i) = h(i)*x1(i)
            c(i) = h(i-1)*x1(i)
         end do
         c(nsglvs) = 0.
      end if

!     Reverse data and heights.
      z = grav/stdlapse
      do j = 1,iy
         tv = tg(:,j,:) * (epsil+qg(:,j,:))/(epsil*(1.+qg(:,j,:)))    
         do i = 1,ix
            ii = 1
            do k = 1,nprlvs
               if ( mtrlvs(k) < height(i,j,1) ) then
                  !mtrphys = mtrlvs(k)/mtrlvs(nprlvs)*(mtrlvs(nprlvs)-zs(i,j)) + zs(i,j)
                  mtrphys = mtrlvs(k)
                  sig(i,nprlvs+1-k,j) = (grav*mtrphys/(tv(i,1)*z)+1.)**(-z/rdry)
               else
                  do while ( height(i,j,ii+1)<mtrlvs(k) .and. ii<nsgm1 )
                     ii = ii + 1
                  end do
                  !mtrphys = (mtrlvs(k)-height(i,j,ii))/mtrlvs(nprlvs)*(mtrlvs(nprlvs)-zs(i,j))
                  mtrphys = mtrlvs(k) - height(i,j,ii)
                  sig(i,nprlvs+1-k,j) = siglvs(nsglvs-ii+1) &
                      *exp(-2.*grav/rdry*mtrphys/(tv(i,ii)+tv(i,ii+1)))
               end if
            end do
         end do
      end do

      do j = 1,iy
!
!     mexup1, mexup2    mexdn1,mexdn2    min1,min2 set
!     the upper and lower limits on upward and downward
!     extrapolation and interpolation.
!     if there are no levels in a given class,the limits are
!     set to zero.
!
         do i = 1,ix
            mexup1(i,j) = 0
            mexup2(i,j) = 0
            do ii = 1,nprlvs
               if ( sig(i,ii,j) < siglvs(1) ) mexup2(i,j) = ii
            end do
            if ( mexup2(i,j) > 0 ) mexup1(i,j) = 1
            mexdn1(i,j) = 0
            mexdn2(i,j) = 0
            do ii = 1,nprlvs
               iii = nprlvs + 1 - ii
               if ( sig(i,iii,j) > siglvs(nsglvs) ) mexdn1(i,j) = iii
            end do
            if ( mexdn1(i,j) > 0 ) mexdn2(i,j) = nprlvs

            if ( mexup2(i,j) == 0 ) then
               !             no upward extrapolation
               if ( mexdn1(i,j) == 0) then
                  min1(i,j) = 1
                  min2(i,j) = nprlvs
               else
                  min1(i,j) = 1
                  min2(i,j) = mexdn1(i,j) - 1
                  if ( mexdn1(i,j) == 1 ) min1(i,j) = 0
               end if
            else 
!
!         upward extrapolation
!
               if ( mexdn1(i,j) == 0 ) then
                  min1(i,j) = mexup2(i,j) + 1
                  min2(i,j) = nprlvs
                  if ( mexup2(i,j) >= nprlvs ) min1(i,j) = 0
                  if ( mexup2(i,j) >= nprlvs ) min2(i,j) = 0
               else
                  min1(i,j) = mexup2(i,j) + 1
                  min2(i,j) = mexdn1(i,j) - 1
                  if ( mexdn1(i,j) == (mexup2(i,j)+1) ) min1(i,j) = 0
                  if ( mexdn1(i,j) == (mexup2(i,j)+1) ) min2(i,j) = 0
               end if
            end if
         end do

!        Calculate jaa, the index of the lowest sigma level above the
!        given height level.
         do kp = 1,nprlvs 
            jaa(:,kp,j) = search_fgt(siglvs(:nsgm1),sig(:,kp,j),nsgm1,ix)
         end do

      end do

      minlev = 1
      maxlev = nsglvs
      
   end subroutine mitop_setup

   
   subroutine titop_setup ( sigr, thetalvs, press, tg, maxlev, minlev )

      ! Setup for mitop. Calculate things that don't depend on the field
      ! being interpolated.
      use physparams
      use utils_m, only : assert, search_fgt
      real, dimension(:), intent(in)       :: sigr      ! Sigma levels
      real, dimension(:), intent(in)       :: thetalvs  ! Theta levels
      real, dimension(:,:,:), intent(in)   :: tg        ! Temperature
      real, dimension(:,:), intent(in)     :: press     ! Surface pressure
      real, dimension(size(tg,dim=1),size(tg,dim=3)) :: thetag
      real, parameter :: cmu1=0.5, clmdam=0.5

      integer, intent(inout) :: maxlev, minlev
      integer :: i, k, ii, ns, iii, kp, ks
      integer :: j, j1, jj 
      integer :: maxlev_g, minlev_g, ierr
      real :: v1, w1, v2, w2, h1, h2, h3, z
      real :: t_est, x

!----------------------------------------------------------------------
!
      if ( ix == 0 ) then
         ! Not initialised yet
         nsglvs = size(sigr)
         nprlvs = size(thetalvs)
         ix = size(press,dim=1)
         iy = size(press,dim=2)
         allocate ( jaa(ix,nprlvs,iy), sig(ix,nprlvs,iy) )
         allocate ( mexdn1(ix,iy), mexdn2(ix,iy), mexup1(ix,iy), &
                    mexup2(ix,iy), min1(ix,iy), min2(ix,iy) )
         allocate ( siglvs(nsglvs), h(nsglvs), x1(nsglvs), x2(nsglvs), &
                    a(nsglvs), c(nsglvs) )

         c(1) = cmu1
         a(nsglvs) = clmdam
!        Sigma levels in increasing order
         siglvs = sigr(nsglvs:1:-1)
         nsgm1 = nsglvs-1
         do i = 1,nsgm1
            h(i) = siglvs(i+1)-siglvs(i)
         end do
         do i = 2,nsgm1
            x1(i) = 0.5/(h(i)+h(i-1))
            x2(i) = h(i)/h(i-1)
            a(i) = h(i)*x1(i)
            c(i) = h(i-1)*x1(i)
         end do
         c(nsglvs) = 0.
      end if

!     Reverse data and theta levels.
      do j = 1,iy
         do k = 1,nsglvs  
            thetag(:,k) = tg(:,j,k)*(press(:,j)*sigr(k)/1.e3)**(-rdry/cp)
         end do  
         do i = 1,ix
            ii = 1         
            do k = 1,nprlvs
               do while ( thetag(i,ii+1)<thetalvs(k) .and. ii<nsgm1 )
                  ii = ii + 1
               end do
               ! estimate temperature based on interpolation of theta
               x = ( thetalvs(k) - thetag(i,ii) ) / ( thetag(i,ii+1) - thetag(i,ii) )
               x = min( max( x, 0. ), 1. )
               t_est = x*tg(i,j,ii+1) + (1.-x)*tg(i,j,ii)
               ! calculate sigma level
               ! press is in hPa so the reference pressure is 1000 hPa.
               sig(i,nprlvs+1-k,j) = (1.e3/press(i,j))*(t_est/thetalvs(k))**(cp/rdry)
            end do
         end do
      end do

      do j = 1,iy
!
!     mexup1, mexup2    mexdn1,mexdn2    min1,min2 set
!     the upper and lower limits on upward and downward
!     extrapolation and interpolation.
!     if there are no levels in a given class,the limits are
!     set to zero.
!
         do i = 1,ix
            mexup1(i,j) = 0
            mexup2(i,j) = 0
            do ii = 1,nprlvs
               if ( sig(i,ii,j) < siglvs(1) ) mexup2(i,j) = ii
            end do
            if ( mexup2(i,j) > 0 ) mexup1(i,j) = 1
            mexdn1(i,j) = 0
            mexdn2(i,j) = 0
            do ii = 1,nprlvs
               iii = nprlvs + 1 - ii
               if ( sig(i,iii,j) > siglvs(nsglvs) ) mexdn1(i,j) = iii
            end do
            if ( mexdn1(i,j) > 0 ) mexdn2(i,j) = nprlvs

            if ( mexup2(i,j) == 0 ) then
               !             no upward extrapolation
               if ( mexdn1(i,j) == 0) then
                  min1(i,j) = 1
                  min2(i,j) = nprlvs
               else
                  min1(i,j) = 1
                  min2(i,j) = mexdn1(i,j) - 1
                  if ( mexdn1(i,j) == 1 ) min1(i,j) = 0
               end if
            else 
!
!         upward extrapolation
!
               if ( mexdn1(i,j) == 0 ) then
                  min1(i,j) = mexup2(i,j) + 1
                  min2(i,j) = nprlvs
                  if ( mexup2(i,j) >= nprlvs ) min1(i,j) = 0
                  if ( mexup2(i,j) >= nprlvs ) min2(i,j) = 0
               else
                  min1(i,j) = mexup2(i,j) + 1
                  min2(i,j) = mexdn1(i,j) - 1
                  if ( mexdn1(i,j) == (mexup2(i,j)+1) ) min1(i,j) = 0
                  if ( mexdn1(i,j) == (mexup2(i,j)+1) ) min2(i,j) = 0
               end if
            end if
         end do

!        Calculate jaa, the index of the lowest sigma level above the
!        given height level.
         do kp = 1,nprlvs 
            jaa(:,kp,j) = search_fgt(siglvs(:nsgm1),sig(:,kp,j),nsgm1,ix)
         end do

      end do

      minlev = 1
      maxlev = nsglvs
      
   end subroutine titop_setup
  

   subroutine ditop_setup ( sigr, mtrlvs, zs )

      ! Setup for mitop. Calculate things that don't depend on the field
      ! being interpolated.
      use physparams
      use utils_m, only : assert, search_fgt
      include 'mpif.h'
      real, dimension(:), intent(in)       :: sigr      ! Sigma pr zstar levels
      real, dimension(:), intent(in)       :: mtrlvs    ! Height levels
      real, dimension(:,:), intent(in)     :: zs        ! bathymetry depth
      integer :: i, j, k, ii

!----------------------------------------------------------------------
!
      if ( oix == 0 ) then
         ! Not initialised yet
         onsglvs = size(sigr)
         onprlvs = size(mtrlvs)
         oix = size(zs,dim=1)
         oiy = size(zs,dim=2)
         allocate( osig(oix,onprlvs,oiy), osiglvs(onsglvs) )
         allocate( ojaa(oix,onprlvs,oiy), ow1(oix,onprlvs,oiy) )

!        Sigma or zstar levels in increasing order
         osiglvs(1:onsglvs) = sigr(1:onsglvs)
      end if  

!     Store data and heights.
      if ( all(osiglvs<=1.) ) then
         ! sigma levels 
         do k = 1,onprlvs
            osig(:,k,:) = mtrlvs(k)/max(zs,1.e-8)
         end do
      else
         ! zstar levels 
         do k = 1,onprlvs
            osig(:,k,:) = mtrlvs(k)
            do j = 1,oiy
              where ( mtrlvs(k)>zs(:,j) )
                osig(:,k,j) = 1.e8 ! below ocean floor
              end where
            end do  
         end do
      end if    
         
      do j = 1,oiy
         do k = 1,onprlvs
            do i = 1,oix
               ojaa(i,k,j) = -1. ! missing
               ow1(i,k,j) = 0.
               if ( osig(i,k,j)>9.e7 ) then
                  ojaa(i,k,j) = -1 ! missing
                  ow1(i,k,j) = 0.
               else if ( osig(i,k,j)<osiglvs(1) ) then
                  ojaa(i,k,j) = 1
                  ow1(i,k,j) = 1.
               else if ( osig(i,k,j)>zs(i,j) .or. osig(i,k,j)>osiglvs(onsglvs) ) then
                  ojaa(i,k,j) = -1 ! missing
                  ow1(i,k,j) = 0.
               else
                  do ii = 1,onsglvs-1
                     if ( osiglvs(ii+1)>zs(i,j) ) then
                        ojaa(i,k,j) = ii
                        ow1(i,k,j) = 1.
                        exit          
                     else if ( osig(i,k,j)>=osiglvs(ii) .and. osig(i,k,j)<=osiglvs(ii+1) ) then
                        ojaa(i,k,j) = ii
                        ow1(i,k,j) = (osig(i,k,j)-osiglvs(ii+1))/min(osiglvs(ii)-osiglvs(ii+1),-1.e-4)
                        ow1(i,k,j) = min( max( ow1(i,k,j), 0. ), 1. )
                        exit
                     end if   
                  end do
               end if
            end do
         end do 
      end do
               
   end subroutine ditop_setup

   
   subroutine ditop ( grids, gridp )


!     This routine interpolates a field "grids" on sigma surfaces
!     to "gridp" on height surfaces
!     This version uses linear interpolation due to missing values
!     for zstar depth coordinates

#ifdef usenc_mod
      use netcdf, only : NF90_FILL_FLOAT
#else
      use netcdf_m, only : NF90_FILL_FLOAT
#endif
      use utils_m, only : assert, search_fgt
      use physparams
      real, dimension(:,:,:), intent(in)   :: grids  ! Sigma field
      real, dimension(:,:,:), intent(out)  :: gridp  ! Height field

      real, dimension(size(grids,dim=1),size(grids,dim=3)) :: ogs
      real, dimension(size(gridp,dim=1),size(gridp,dim=3)) :: ogp
      integer :: i, j, ns, np, k
      real :: w1

!----------------------------------------------------------------------
!
      call assert ( oix /= 0, "Error, ditop called before ditop_setup")
!     Check all the array sizes match
      call assert ( onsglvs == size ( grids, dim=3 ), &
                  "Error, number of sigma levels doesn't match in ditop" )
      call assert ( onprlvs == size ( gridp, dim=3 ), &
                  "Error, number of height levels doesn't match in ditop" )

      do j = 1,oiy

         do ns = 1,onsglvs
            ogs(:,ns) = grids(:,j,ns)
         end do

         do np = 1,onprlvs
            do i = 1,oix
               k = ojaa(i,np,j)
               if ( k>=1 .and. k<onsglvs ) then
                  w1 = ow1(i,np,j)
                  if ( ogs(i,k)==NF90_FILL_FLOAT .or. ogs(i,k+1)==NF90_FILL_FLOAT ) then
                     ogp(i,np) = NF90_FILL_FLOAT
                  else   
                     ogp(i,np) = ogs(i,k)*w1 + ogs(i,k+1)*(1.-w1) ! linear interpolation
                  end if   
               else 
                  ogp(i,np) = NF90_FILL_FLOAT 
               end if    
            end do  
         end do 
      
         do np = 1,onprlvs
            gridp(:,j,np) = ogp(:,np)
         end do

      end do ! j=1,oiy  

   end subroutine ditop

   
end module sitop_m
