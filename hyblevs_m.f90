module hyblevs_m
#ifndef parnetcdf
   use netcdf_m
#else
   use pnetcdf_m
#endif
   implicit none
   public :: hyblevs, hybtop
   private :: gaussj
   real, parameter :: Pooeta = 1013.2
   real, save, dimension(:), allocatable :: sig, sigh, anf, bnf, anh, bnh, &
                                            dadnf, dbdnf

   ! For interpolation to pressure levels

   ! Controls whether vertical extrapolation below surface is done
   integer, parameter :: vextrap_default=0, vextrap_lin=1, vextrap_none=2, &
                         vextrap_missing=3, vextrap_t=4
   ! This calculation needs 64 bit precision.
   integer, private, parameter :: r8 = selected_real_kind(p=15, r=307)

contains

   subroutine hyblevs ( nl, oldhybrid )

!  Calculate coefficients for the model hybrid levels
   integer, intent(in) :: nl
   logical, intent(in), optional :: oldhybrid
   integer :: nlp
   integer, parameter :: neq=10, neq2=neq/2
   real(kind=r8), dimension(neq,1)   :: c
   real(kind=r8), dimension(neq,neq) :: uk2
   real(kind=r8) :: rk, alf
   integer :: i, j, k
   real(kind=r8) :: ptm1, ptm2, alpha, p500, p05, suma, sumb, ex, apb
   logical :: oldhyb
   real(kind=r8), dimension(nl)   :: sig_d, anf_d, bnf_d, dadnf_d, dbdnf_d
   real(kind=r8), dimension(nl+1) :: sigh_d, anh_d, bnh_d

   if ( present(oldhybrid) ) then
      oldhyb = oldhybrid
   else
      oldhyb = .false.
   end if

!      read(*,*) k

   nlp = nl + 1
   if ( allocated(sig) ) then
      print*, "Error, sig already allocated in hyblevs"
      stop
   end if
   allocate ( sig(nl), sigh(nlp), anf(nl), bnf(nl), anh(nlp), bnh(nlp), &
              dadnf(nl), dbdnf(nl) )
   alf = 1.0
   do k = 1,nl+1
      rk = (.5*nl-(k-1))/nl
      sigh_d(k) = 0.5 + 1.5*alf*rk - 2.*(3.*alf-2.)*rk**3
   enddo
   sigh_d(1) = 1.0
   sigh_d(nl+1) = 0.0

   sig_d = 0.5 * ( sigh_d(1:nl) + sigh_d(2:nlp) ) 

!  Create the elements of matrics (equations) and solve

!  Top two half levels to be constant pressure
   ptm2 = Pooeta*sigh_d(nlp-2)
   ptm1 = Pooeta*sigh_d(nlp-1)

!  Lowest two levels thickness at 1013.2mb = alpha*(lowest two levels
!  thickness at 500mb)
   alpha = (500.0-ptm2)/(Pooeta-ptm2)
   alpha = alpha*0.5    
   p500 = alpha*Pooeta-500.0
   p05 = (1.0-alpha)*Pooeta

   c   = 0.0
   uk2 = 0.0

!  condition 1
   uk2(1,1:neq2) = 1.0
!  condition 2
   uk2(2,1+neq2:neq) = 1.0
   c(2,1) = 1.0
!  condition 3
   do j = 1,neq2
      uk2(3,j+neq2) = sigh_d(nlp-2)**j
   end do
!  condition 4
   do j = 1,neq2
      uk2(4,j+neq2) = sigh_d(nlp-1)**j
   end do
!  condition 5
   do j = 1,neq2
      uk2(5,j) = sigh_d(nlp-2)**j
   end do
   c(5,1) = sigh_d(nlp-2)
!  condition 6
   do j = 1,neq2
      uk2(6,j)  = sigh_d(nlp-1)**j
   end do
   c(6,1) = sigh_d(nlp-1)
!  condition 7
   do j = 1,neq2
      uk2(7,j) = p05*(sigh_d(1)**j-sigh_d(2)**j)
      uk2(7,j+neq2) = -p500*(sigh_d(1)**j-sigh_d(2)**j)
   end do
!  condition 8
   do j = 1,neq2
      uk2(8,j) = p05*(sigh_d(2)**j-sigh_d(3)**j)
      uk2(8,j+neq2) = -p500*(sigh_d(2)**j-sigh_d(3)**j)
   end do
!  conditions 9 and 10 (ok for nl> = 7)
   if (nl < 7) then
      print *,'hybrid 9 & 10 not ok for nl<7'
      stop
   endif
   i = (4*nl)/9
   do j = 1,neq2
      uk2(9,j) = sigh_d(i)**j
      uk2(9,j+neq2) = sigh_d(i)**j
   end do
   c(9,1) = sigh_d(i)
   i = (6*nl)/9
   do j = 1,neq2
      uk2(10,j) = sigh_d(i)**j
      uk2(10,j+neq2) = sigh_d(i)**j
   end do
   c(10,1) = sigh_d(i)

!  FIND THE MATRIX SOLUTION: uk2(i,j) * xx(i,1) = c(i,1) ****
   call gaussj(uk2,neq,neq,c,1,1)
!  MATRIX c(i,1) WILL NOW BE REPLACED BY THE SOLUTION xx(i,1) 

!  A=anh(k) and B=bnh(k) at each half level 
   if ( oldhyb ) then

      do k = 1,nlp
         suma = 0.0
         sumb = 0.0
         do j = neq2,1,-1
            ex = 0.0
            if ( sigh_d(k) /= 0.0 ) ex = sigh_d(k)**j
            suma = suma + c(j,1)*ex
            sumb = sumb + c(j+neq2,1)*ex
         end do
         anh_d(k) = suma
         bnh_d(k) = sumb
         apb = suma+sumb  
      end do
!  Set A (half level 1)  =  B (half level 8,9)  =  0.0 (exactly) as 
!  defined in conditions 1,3 and 4 (see notes)
      anh_d(1) = 0.0
      bnh_d(nlp-1) = 0.0
      bnh_d(nlp-2) = 0.0

!  For consistency, and energy conservation, the values
!  at FULL levels must be computed from HALF level values
      do k = 1,nl
         anf_d(k) = 0.5*(anh_d(k)+anh_d(k+1))
         bnf_d(k) = 0.5*(bnh_d(k)+bnh_d(k+1))
         dadnf_d(k)=(anh_d(k)-anh_d(k+1))/(sigh_d(k)-sigh_d(k+1))
         dbdnf_d(k)=(bnh_d(k)-bnh_d(k+1))/(sigh_d(k)-sigh_d(k+1))
      end do

   else ! not oldhybrid Mk3++

!----
!---- More precise calculations for hybrid levels
!----  taking into account machine round off errors (Mk3++)
!----

!**** A=anh(k) and B=bnh(k) at each half level 
      do k=1,nlp
         suma=0.0
         sumb=0.0
!DIR$ NEXTSCALAR
!PDIR NOVECTOR
         do j=neq2,1,-1
            ex=0.0
            if(sigh_d(k).ne.0.0)ex=sigh_d(k)**j
            suma=suma+c(j,1)*ex
            sumb=sumb+c(j+neq2,1)*ex
         end do
         anh_d(k)=suma
         bnh_d(k)=sumb
      end do
!PDIR VECTOR

!**** Certain conditions should be met.
!**** These are :
!****  (1) A (half level 1) = B (half level nl & nl-1) = 0.0 (exactly)
!****      as defined in conditions 1,3 and 4 (see notes)
!****  (2) A (half level) + B (half level) = sigh_d (half level)
!****  (3) A (full level) + B (full level) = sig (full level)
!****  (4) dAdNf (full level) + dBdNf (full level) = 1.0
!****     with dAdNf(nl)=dAdNf(nl-1)=1.0 & dBdNf(nl)=dBdNf(nl-1)=0.0
!****
!**** However, due to the calculations involved, and round off,
!****  these criteria do not follow identically. 
!**** So adjustments will be undertaken to guarantee the above.
!**** These adjustments are minor (last decimal places), and
!****  will be applied first to both numbers equally, and if
!****  that is not sufficient, then the smallest number will
!****  be adjusted to guarantee the above conditions. See below.

!****  (1) A (half level 1) = B (half level nl & nl-1) = 0.0 (exactly)
      anh_d(1)=0.0
      bnh_d(nlp-1)=0.0
      bnh_d(nlp-2)=0.0
!****  (2) A (half level) + B (half level) = sigh_d (half level)
!****   Adjust (if needed) so that anh+bnh = sigh_d
      anh_d(nlp-1)=sigh_d(nlp-1)
      anh_d(nlp-2)=sigh_d(nlp-2)
      do k=2,nlp-3
         apb=(anh_d(k)+bnh_d(k))-sigh_d(k)
         anh_d(k)=anh_d(k)-apb/2
         bnh_d(k)=bnh_d(k)-apb/2
         apb=anh_d(k)+bnh_d(k)
         if(apb.ne.sigh_d(k))then
            if(anh_d(k).gt.bnh_d(k))then
               bnh_d(k)=sigh_d(k)-anh_d(k)
            else
               anh_d(k)=sigh_d(k)-bnh_d(k)
            endif
         endif
      end do
      
!****  (3) A (full level) + B (full level) = sig (full level)
!**** For consisentcy, and energy conservation, the values
!**** at FULL levels must be computed from HALF level values
      do k=1,nl
         anf_d(k)=0.5*(anh_d(k)+anh_d(k+1))
         bnf_d(k)=0.5*(bnh_d(k)+bnh_d(k+1))
         apb=(anf_d(k)+bnf_d(k))-sig_d(k)
         anf_d(k)=anf_d(k)-apb/2
         bnf_d(k)=bnf_d(k)-apb/2
         apb=anf_d(k)+bnf_d(k)
         if(apb.ne.sig_d(k))then
            if(anf_d(k).gt.bnf_d(k))then
               bnf_d(k)=sig_d(k)-anf_d(k)
            else
               anf_d(k)=sig_d(k)-bnf_d(k)
            endif
         endif
      end do

!****  (4) dAdNf (full level) + dBdNf (full level) = 1.0
!****     with dAdNf(nl)=dAdNf(nl-1)=1.0 & dBdNf(nl)=dBdNf(nl-1)=0.0
      do k=1,nl-2
         dadnf_d(k)=(anh_d(k)-anh_d(k+1))/(sigh_d(k)-sigh_d(k+1))
         dbdnf_d(k)=(bnh_d(k)-bnh_d(k+1))/(sigh_d(k)-sigh_d(k+1))
         apb=(dadnf_d(k)+dbdnf_d(k))-1.0
         dadnf_d(k)=dadnf_d(k)-apb/2
         dbdnf_d(k)=dbdnf_d(k)-apb/2
         apb=dadnf_d(k)+dbdnf_d(k)
         if(apb.ne.1.0)then
            if(dadnf_d(k).gt.dbdnf_d(k))then
               dbdnf_d(k)=1.0-dadnf_d(k)
            else
               dadnf_d(k)=1.0-dbdnf_d(k)
            endif
         endif
      end do
      dadnf_d(nl)=1.0
      dbdnf_d(nl)=0.0
      dadnf_d(nl-1)=1.0
      dbdnf_d(nl-1)=0.0

   endif ! (oldhybrid)

   anf_d = anf_d * Pooeta
   anh_d = anh_d * Pooeta
   dadnf_d = dadnf_d * Pooeta

   sig = real(sig_d)
   sigh = real(sigh_d)
   anf = real(anf_d)
   bnf = real(bnf_d)
   anh = real(anh_d)
   bnh = real(bnh_d)
   dadnf = real(dadnf_d)
   dbdnf = real(dbdnf_d)

end subroutine hyblevs

   subroutine gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL(kind=r8) a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL(kind=r8) big,dum,pivinv
      do 11 j=1,n
      ipiv(j)=0
11    continue
      do 22 i=1,n
      big=0.
      do 13 j=1,n
      if(ipiv(j).ne.1)then
      do 12 k=1,n
      if (ipiv(k).eq.0) then
      if (abs(a(j,k)).ge.big)then
      big=abs(a(j,k))
      irow=j
      icol=k
      endif
      else if (ipiv(k).gt.1) then
      print*, 'singular matrix in gaussj'
      stop
      endif
12    continue

      endif
13    continue
      ipiv(icol)=ipiv(icol)+1
      if (irow.ne.icol) then
      do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
14    continue
      do 15 l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
15    continue
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol).eq.0.) then
         print*, 'singular matrix in gaussj'
         stop
      end if
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do 16 l=1,n
      a(icol,l)=a(icol,l)*pivinv
16    continue
      do 17 l=1,m
      b(icol,l)=b(icol,l)*pivinv
17    continue
      do 21 ll=1,n
      if(ll.ne.icol)then
      dum=a(ll,icol)
      a(ll,icol)=0.
      do 18 l=1,n
      a(ll,l)=a(ll,l)-a(icol,l)*dum
18    continue
      do 19 l=1,m
      b(ll,l)=b(ll,l)-b(icol,l)*dum
19    continue
      endif
21    continue
22    continue
      do 24 l=n,1,-1
      if(indxr(l).ne.indxc(l))then
      do 23 k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
23    continue
      endif
24    continue
      return
   end subroutine gaussj

   subroutine hybtop ( grids, gridp, prelvs, press, vextrap )

!     This routine interpolates a field "grids" on the model hybrid
!     sigma-pressure surfaces to "gridp" on pressure surfaces

!     The hybrid level coefficients should have been set up by hyblevs.

!     The interpolation is performed using cubic splines
!     prelvs(nprlvs) are the pressure levels and siglvs(nsglvs) those
!     of sigma
!     "press" is the surface pressure and "ix" and "iy" specify the
!     size of the horizontal grid
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

      use netcdf, only : NCF90_FILL_FLOAT
      use utils_m, only : assert, search_fgt
      use physparams
      ! Lapse rate used in temperature extrapolation
      real, parameter :: extrap_lapse = 6.5e-3
      real, parameter :: konst = extrap_lapse*rdry/grav
      real, dimension(:,:,:), intent(in)   :: grids  ! Sigma field
      real, dimension(:,:,:), intent(out)  :: gridp  ! Pressure field
      real, dimension(:), intent(in)       :: prelvs ! Pressure levels
      real, dimension(:,:), intent(in)     :: press  ! Surface pressure
      integer, intent(in) :: vextrap       ! Controls vertical extrapolation

      real, dimension(size(grids,dim=1),size(grids,dim=3)) :: d, gd, gs
      integer, dimension(size(grids,dim=1)) :: mexdn1, mexdn2, mexup1, mexup2,&
                                               min1, min2
      integer, dimension(size(grids,dim=1),size(prelvs)) :: jaa
      real,    dimension(size(grids,dim=1),size(prelvs)) :: gp, sigout
      real, dimension(size(grids,dim=1),size(grids,dim=3)) :: sigin, h, x1, x2, a, c
      real, dimension(size(grids,dim=1)) :: presinv, x3
      real, parameter :: cmu1=0.5, clmdam=0.5
      integer :: nsglvs, nprlvs, ix, iy, i, k, kr, ii, nsgm1, ns, iii, kp, ks
      integer :: minmin1, maxmin1, maxmin2, j, j1, jj 
      real :: v1, w1, v2, w2, h1, h2, h3, z

!----------------------------------------------------------------------
!
      nprlvs = size(prelvs)
      ix = size(grids,dim=1)
      iy = size(grids,dim=2)
      nsglvs = size(grids,dim=3)
      if ( .not. allocated(sig) ) then
         print*, "Error, sig not allocated in hybtop. Call hyblevs first."
         stop
      end if

!     Check all the array sizes match
      call assert ( nsglvs == size ( grids, dim=3 ), &
                  "Error, number of sigma levels doesn't match in hybtop" )
      call assert ( nprlvs == size ( gridp, dim=3 ), &
                  "Error, number of pressure levels doesn't match in hybtop" )

      call assert ( ix == size(gridp,dim=1) .and. ix  == size(press,dim=1), &
                  "Error, nlon doesn't match in hybtop" )
      call assert ( iy == size(gridp,dim=2) .and. iy  == size(press,dim=2), &
                  "Error, nlat doesn't match in hybtop" )

      nsgm1 = nsglvs-1

      do j=1,iy
     
         presinv = 1.0/press(:,j)

         c(:,1)=cmu1
         a(:,nsglvs)=clmdam

         do k=1,nsglvs
            ! Effective sigma levels in increasing order
            kr = nsglvs+1-k
            sigin(:,k) = ( anf(kr) + press(:,j)*bnf(kr) ) *presinv
         end do

         do i=1,nsgm1
            h(:,i) = sigin(:,i+1)-sigin(:,i)
         end do
         do i=2,nsgm1
            x1(:,i) = 0.5/(h(:,i)+h(:,i-1))
            x2(:,i) = h(:,i)/h(:,i-1)
            a(:,i) = h(:,i)*x1(:,i)
            c(:,i) = h(:,i-1)*x1(:,i)
         end do
         c(:,nsglvs)=0.0

!        Reverse data and pressure levels.
         do k=1,nprlvs
            sigout(:,nprlvs+1-k) = prelvs(k)*presinv(:)
         end do
         do ns=1,nsglvs
            gs(:,nsglvs+1-ns) = grids(:,j,ns)
         end do
         d(:,1) = 1.5*(gs(:,2)-gs(:,1)) / h(:,1)
         d(:,nsglvs) = 1.5*(gs(:,nsglvs)-gs(:,nsglvs-1)) / h(:,nsglvs-1)
         do ns=2,nsgm1
            d(:,ns) = 3.0*x1(:,ns)*( (gs(:,ns)-gs(:,ns-1))*x2(:,ns) + &
                                   (gs(:,ns+1)-gs(:,ns))/x2(:,ns) )
         end do
!
!     calculate spline derivatives at orginal points
!
         gd(:,1)=-c(:,1)
         do ns=2,nsglvs
            x3 = 1.0 / (1.0+a(:,ns)*gd(:,ns-1))
            gd(:,ns) = -c(:,ns)*x3
            d(:,ns) = (d(:,ns)-a(:,ns)*d(:,ns-1)) * x3
         end do
         gd(:,nsglvs) = d(:,nsglvs)
         do ns=1,nsgm1
            k = nsglvs-ns
            gd(:,k) = gd(:,k)*gd(:,k+1)+d(:,k)
         end do
!
!     mexup1, mexup2    mexdn1,mexdn2    min1,min2 set
!     the upper and lower limits on upward and downward
!     extrapolation and interpolation.
!     if there are no levels in a given class,the limits are
!     set to zero.
!
         do i=1,ix
            mexup1(i)=0
            mexup2(i)=0
            do ii=1,nprlvs
               if(sigout(i,ii).lt.sigin(i,1)) mexup2(i)=ii
            end do
            if(mexup2(i).gt.0) mexup1(i)=1
            mexdn1(i)=0
            mexdn2(i)=0
            do ii=1,nprlvs
               iii=nprlvs+1-ii
               if(sigout(i,iii).gt.sigin(i,nsglvs)) mexdn1(i)=iii
            end do
            if(mexdn1(i).gt.0) mexdn2(i)=nprlvs

            if (mexup2(i).eq.0) then
               !             no upward extrapolation
               if(mexdn1(i).eq.0) then
                  min1(i)=1
                  min2(i)=nprlvs
               else
                  min1(i)=1
                  min2(i)=mexdn1(i)-1
                  if(mexdn1(i).eq.1) min1(i)=0
               end if
            else 
!
!         upward extrapolation
!
               if (mexdn1(i).eq.0) then
                  min1(i)=mexup2(i)+1
                  min2(i)=nprlvs
                  if(mexup2(i).ge.nprlvs) min1(i)=0
                  if(mexup2(i).ge.nprlvs) min2(i)=0
               else
                  min1(i)=mexup2(i)+1
                  min2(i)=mexdn1(i)-1
                  if(mexdn1(i).eq.(mexup2(i)+1)) min1(i)=0
                  if(mexdn1(i).eq.(mexup2(i)+1)) min2(i)=0
               end if
            end if
         end do
         
!     Find the extreme values of the interpolation. The assumption
!     used here is that most of the values will be interpolated and there
!     will be just a few odd extrapolations. Therefore it is worth 
!     vectorising the interpolation section, even if a few values are
!     calculated unnecessarily. These will be fixed up in the extrapolation
!     section
         minmin1 = minval(min1)
         maxmin1 = maxval(min1)
         maxmin2 = maxval(min2)
!
!     Do interpolation
!
         if(maxmin1.eq.0) go to 300

!        Calculate jaa, the index of the lowest sigma level above the
!        given pressure level.
         do kp=1,nprlvs 
            jaa(:,kp) = search_fgt(sigin(:,:nsgm1),sigout(:,kp),nsgm1,ix)
         end do

         do ii=max(1,minmin1),maxmin2
            do i=1,ix
!           Use max to ensure stay in array bounds. Any points where this
!           is used will be fixed in the extrapolation section.
               j1=max(jaa(i,ii)-1,1)
               jj=jaa(i,ii)
               v1=sigin(i,jaa(i,ii))-sigout(i,ii)
               w1=sigout(i,ii)-sigin(i,j1)
               v2=v1*v1
               w2=w1*w1
               h1=h(i,j1)
               h2=1.0/(h1*h1)
               h3=h2/h1
               z=(gd(i,j1)*v2*w1-gd(i,jaa(i,ii))*w2*v1)*h2
               gp(i,ii) = z + (gs(i,j1)*v2*(w1+w1+h1)+  &
                             gs(i,jaa(i,ii))*w2*(v1+v1+h1))*h3
            end do
         end do
!
!     linear extrapolation
!
!     here we use the options of extrapolation as explained above
!
!
300      continue
         select case ( vextrap )
         case ( vextrap_lin ) 
            do i=1,ix
               if(mexup1(i).ne.0) then
                  do ii=1,mexup2(i)
                     gp(i,ii) = gs(i,1)+gd(i,1)*(sigout(i,ii)-sigin(i,1))
                  end do
               end if
               if(mexdn1(i).ne.0) then
                  do ii=mexdn1(i),nprlvs
                     gp(i,ii) = gs(i,nsglvs)+gd(i,nsglvs)*(sigout(i,ii)-sigin(i,nsglvs))
                  end do
               end if
            end do
         case ( vextrap_none )
            do i=1,ix
               do ii=1,mexup2(i)
                  gp(i,ii) = gs(i,1)
               end do
               if(mexdn1(i).ne.0) then
                  do ii=mexdn1(i),nprlvs
                     gp(i,ii) = gs(i,nsglvs)
                  end do
               end if
            end do
         case ( vextrap_missing )
            do i=1,ix
               do ii=1,mexup2(i)
                  gp(i,ii) = NCF90_FILL_FLOAT
               end do
               if(mexdn1(i).ne.0) then
                  do ii=mexdn1(i),nprlvs
                     gp(i,ii) = NCF90_FILL_FLOAT
                  end do
               end if
            end do
         case ( vextrap_t )
            ! Extrapolate using the specified lapse rate
            ! temp = ts[0,0] * sig ** k
            do i=1,ix
               ! Linear extrapolation at the top.
               do ii=1,mexup2(i)
                  gp(i,ii) = gs(i,1)+gd(i,1)*(sigout(i,ii)-sigin(i,1))
               end do
               if(mexdn1(i).ne.0) then
                  do ii=mexdn1(i),nprlvs
                     gp(i,ii) = gs(i,nsglvs)*(sigout(i,ii)/sigin(i,nsglvs))**konst
                  end do
               end if
            end do
         case default
            print*, "Error in hybtop, unexpected value of vextrap", vextrap
            stop
         end select

         do ii=1,nprlvs
            gridp(:,j,nprlvs+1-ii)=gp(:,ii)
         end do

      end do

   end subroutine hybtop

end module hyblevs_m
   

