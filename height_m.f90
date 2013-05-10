module height_m
   implicit none
   real, private, allocatable, save, dimension(:,:) :: am
   real, private, allocatable, save, dimension(:) :: algf, alogg

   real, private, allocatable, save, dimension(:) :: algh, atm, atp, btm, btp
contains

subroutine initheight(nl,sig)

!  Calculate the A matrix needed to integrate hydrostatic equation
!  From spectral model inital.f.
!  Should really use a calculation compatible with the CC model's dynamics.

   use physparams

   integer, intent(in) :: nl
   real, intent(in), dimension(:) :: sig

   real, dimension(nl+1) :: sigh
   integer :: i, j, k

   allocate ( am(nl,nl), algf(nl), alogg(nl) )
   allocate ( algh(nl), atm(nl), atp(nl), btm(nl), btp(nl) )

!  Calculate sigh from sigma
   sigh(1) = 1.0
   do k=1,nl-1
      sigh(k+1) = 2*sig(k) - sigh(k)
   end do
   sigh(nl+1) = 0.

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
   return
end subroutine initheight

!-----------------------------------------------------------------
subroutine height ( tg, qg, zg, pg, sig, phistd, pstd )

!  Calculate geopotential height using virtual temperature
!  If the optional argument pstd is present it specifies the pressure levels
!  on which to calculate the height. Otherwise return sigma level heights.
   use physparams
   use utils_m, only : search_fgt
   real, parameter :: extrap_lapse = 6.5e-3
   real, parameter :: konst = extrap_lapse*rdry/grav
   real, intent(in), dimension(:,:,:) :: tg, qg
   real, intent(in), dimension(:,:)    :: zg, pg
   real, intent(in), dimension(:) :: sig
   real, intent(in), dimension(:), optional :: pstd
   real, intent(out), dimension(:,:,:) :: phistd

   real, dimension(size(tg,1),size(tg,3)) :: tv, phi
   real, dimension(size(sig)) :: isig
   integer :: j, k, lg, mg, kstd, kk
   integer, dimension(size(tg,1)) :: jlev
   real, dimension(size(tg,1)) :: siglev, slog
   real :: temp, qemp1, qemp2

   integer :: nx, ny, nl, nstd

   nx = size(tg,1)
   ny = size(tg,2)
   nl = size(sig)

   if ( present ( pstd ) ) then
      nstd = size(pstd)
   else
      nstd = nl
   end if

   if ( nstd /= size(phistd,3) ) then
      print*, "Error, Incorrect size of array phistd in routine height"
      stop
   end if

   do lg=1,ny

!     Calculate virtual temperature.
      tv = tg(:,lg,:) * (epsil+qg(:,lg,:))/(epsil*(1.+qg(:,lg,:)))

!     Calculate height on sigma levels
      do k=1,nl
         phi(:,k) = zg(:,lg) * grav
         do j=1,nl
            phi(:,k) = phi(:,k)+am(k,j)*tv(:,j)
         end do
      end do

!     Now calculate the height at the standard pressure levels.
      if ( present ( pstd ) ) then
         
         isig(1:nl) = sig(nl:1:-1)
         
         do kstd = 1,nstd
            siglev = pstd(kstd)/pg(:,lg)
            slog = alog(siglev)
            ! Here sig is in decreasing order
            
            jlev = nl+2 - search_fgt(isig, siglev, nl,nx)
            do mg=1,nx
               if ( jlev(mg) > 1 .and. jlev(mg) <= nl ) then
                  ! Interpolate the temperature in log sigma
                  qemp1 = ( slog(mg) - algf(jlev(mg)-1) )*tv(mg,jlev(mg))
                  qemp2 = ( algf(jlev(mg)) - slog(mg) )*tv(mg,jlev(mg)-1)
                  temp = (qemp1+qemp2)*alogg(jlev(mg))
               else if ( jlev(mg) == 1 ) then
                  ! Below bottom level
                  ! Extrapolate temperature using std lapse rate
                  temp = tv(mg,1) * (siglev(mg)/sig(1)) ** konst
               else
                  ! Above top level
                  jlev(mg) = nl
                  temp = tv(mg,nl)
               end if
!              Use average of interpolated and level temperature
               temp = 0.5 * ( temp + tv(mg,jlev(mg)) )
               phistd(mg,lg,kstd) = phi(mg,jlev(mg))-rdry*temp*  &
                                            (slog(mg) - algf(jlev(mg)))
            end do
         end do    ! Std level loop

      else ! Return sigma level heights
         
         phistd(:,lg,:) = phi(:,:)

      end if
   end do       ! ny loop

   phistd = phistd / grav

end subroutine height

end module height_m
