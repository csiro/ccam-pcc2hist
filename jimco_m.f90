module jimco_m

   use nfft_m
   use jim_utils
   use precis_m
   implicit none
   private

   public :: jimco
   private :: inoct, otoc1, indoct, ininmap, toct, cinvrt, otoc2, &
              ctoo, conv, getrot

!     conformal octagon - code of Jim Purser
!     JMcG usage:
!            1) calculate xx4, yy4, zz4 just for NE corner of regular octagon 1
!            2) in setxyz apply to all 14 panels
!            3) in setxyz convert panel values via schmidt
!            4) in setxyz save only xxy4, yy4 (after mapping via schmidta)

!  common / cstoct /
   complex(kind=rx), private, save :: ci, cir, ciq, ciap
   real(kind=rx), dimension(40), private, save :: a1, a2, b1, b2
   integer, private, save :: n1
   real(kind=rx), private, save :: a, b, cx, cy, ssr, sc, sci, zac, st, sti
   real(kind=rx), private, save, dimension(3,3) :: rotm

contains

subroutine jimco ( il, iquad, xx4, yy4, zz4 )
   integer, intent(in) :: il, iquad
   real(kind=rx), intent(out), dimension(:,:) :: xx4, yy4, zz4   ! (iquad,iquad)
   integer, parameter :: n = 128

   real(kind=rx), dimension(1200) :: work
   real(kind=rx), dimension(3) :: xe !, xep, dxedx, dxedy, gxedx, gxedy, hxedx, hxedy
   real(kind=rx), dimension(0:n) :: ra, qa
   real(kind=rx), dimension(n)   :: jumble
   integer :: i, j
   real(kind=rx) :: x, y

!  call inoct(n1grid,nagrid,n1,clon,clat,sc,ra,qa,n,jumble,work)
!             24      8    20  0.  90.  1.      128
   call inoct ( 6*il, 2*il, 20, 0.0_rx, 90.0_rx, 1.0_rx, ra, qa, 128, jumble, work )

!     print *,'il,jl,iquad,arrsize ',il,jl,iquad,
!    .     18*il*jl - iquad*iquad -1200
   do j=1,iquad
      y = (j-1)/real(iquad-1)
      do i=1,iquad
         x = (i-1)/real(iquad-1)
!        print *,'i,j,n1,iquad,x,y ',i,j,n1,iquad,x,y
         call otoc1 ( x, y, xe )
!        jlm chooses just a positive quadrant
!        print *,'xea,n1               ',xe,n1
         xx4(i,j) = xe(1)
!        print *,'xeb,n1               ',xe,n1
         yy4(i,j) = xe(2)
!        print *,'xec,n1               ',xe,n1
         zz4(i,j) = xe(3)
!        print *,'xed,n1               ',xe,n1
      enddo
   enddo

end subroutine jimco

!---------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1995
!                   subroutine inoct
!   initialize the constants of common/cstoct/ used to describe grid and to
!   perform transformations between map coordinates and the sphere or disk
!
! --> n1grid    number of grid-spaces from center to side of octagon
! --> nagrid    number of grid-spaces from center-line to nearest vertex
! --> ntay      number of taylor series coefficients retained (ntay<41)
! --> flon0     longitude of center of octagon-1
! --> flat0     latitude of center of octagon-1
! --> scen      scale-enhancement parameter = radius of image of octagon-1 in
!               a stereographic projection centered on (flat0,flon0) that maps
!               the concentric hemisphere to the unit disk
! --- ra,qa     work arrays of dimension n+1 (n defined below)
! --> n         number of the form 2**m (eg 128) and > 3*ntay defining the
!               period of fft used in the construction of taylor series
! --- jumble,work   work arrays of dimension n
!-----------------------------------------------------------------------------

subroutine inoct ( n1grid, nagrid, ntay, flon0, flat0, scen, ra, qa, n,  &
                   jumble, work )
   integer, intent(in) :: n1grid, nagrid, ntay, n
   real(kind=rx), intent(in) :: flon0, flat0, scen
   real(kind=rx), dimension(0:), intent(inout) :: ra, qa
   real(kind=rx), dimension(:), intent(inout) :: jumble, work
   integer :: m1g, mag

   n1 = ntay
   if ( 3*n1 > n ) then
      print*, "too small in inoct, choose a larger power of 2"
      stop
   end if
   m1g = n1grid ! grid intervals between center and side
   mag = nagrid ! grid intervals between center and first corner
   sc = scen
   sci = 1.0/sc
   zac = (sci-sc)/(sci+sc)    !   (1-c**2)/(1+c**2)  jlm
   a = real(mag,kind=rx)/m1g ! nice ratios: 1/3, 3/7, 7/17
!  print *,'calling indoct, n1= ',n1
   call indoct ( ra, qa, n, jumble, work )
!  print *,'calling ininmap'
   call ininmap ( flon0, flat0, rotm )

end subroutine inoct

subroutine indoct ( ra, qa, n, jumble, work )

   integer, intent(in) :: n
   real(kind=rx), intent(out), dimension(0:) :: ra, qa
   real(kind=rx), intent(out), dimension(:) :: jumble, work  ! (n)
   complex(kind=rx) :: cia, cirs2, z, w, cang, cdang1, cdang2
   integer :: nh, kit, i
   real(kind=rx) :: r2, orc, sw1, sw2, rs1, rs2, rs3, rs4, pio4, u, v, rs3p, rs4p

   st = 4.0/3.0
   sti = 3.0/4.0
   ci = cmplx(0.0, 1.0, kind=rx)
   cia = ci*a
   ciap = cia + 1.0
   cir = (-ci)**0.75_rx
   ciq = ci**0.25_rx
   nh = n/2
   r2 = sqrt(2.0_rx)
   orc = 0.5_rx            ! relaxation factor
   sw1 = sqrt(1.0+a*a)
   sw2 = min ( a*2, r2*(1.0-a) )  ! Really just min?????
   ssr = (sw1/sw2)**2
   rs1 = 0.95_rx*sw1
   rs2 = 0.95_rx*sw2
   rs3 = rs1**4
   rs4 = rs2**st

   a1(1:n1) = 0.0

!  First guess for first few taylor coefficients (based on a=1/3)
!  needed to begin the boot-strap procedure:
   a1(1) =  1.05060_rx
   a1(2) = -0.172629_rx
   a1(3) =  0.125716_rx
   a1(4) =  0.00015657_rx
   a1(5) = -0.0017866_rx
   a1(6) = -0.0031168_rx

   a2(1:n1) = 0.0
   
   a2(1) =  0.688411_rx
   a2(2) = -0.180924_rx
   a2(3) = -0.186593_rx
   a2(4) =  0.024811_rx
   a2(5) = -0.044888_rx
   a2(6) = -0.049190_rx
   b     = -0.147179_rx

   jumble(1) = 0.0
   pio4 = atan(1.0_rx)
   cdang1 = ci*pio4/nh
   cdang2 = -3.0*cdang1
   cirs2 = ci*rs2

   do kit = 1, 50
      do i = 0, nh
         cang = cdang1*i
         z = rs1*exp(cang)
         call toct ( z, w )
         w = w**4
         u = real(w)
         v = aimag(w)
         ra(i)   = u
         qa(i)   = -v
         ra(n-i) = u
         qa(n-i) = v
      enddo
      qa(0)  = 0.0
      qa(nh) = 0.0
      call cfft ( ra, qa, n, 1.0_rx, work, jumble )
      rs3p = 1.0
      do i = 1,n1
         rs3p = rs3p*rs3
         a1(i) = a1(i) + orc*(ra(i)/rs3p-a1(i))
      enddo

      do i = 0,nh
         cang = cdang2*i
         z = ciap-cirs2*exp(cang)
         call toct(z,w)
         w = ci*(w-1.0)/(w+1.0)
         u = real(w)
         v = aimag(w)
         ra(i) = u
         qa(i) = v
         ra(n-i) = u
         qa(n-i) = -v
      enddo
      qa(0) = 0.0
      qa(nh) = 0.0
      call cfft ( ra, qa, n, 1.0_rx, work, jumble )
      b = b + orc*(ra(0)-b)
      rs4p = 1.0
      do i=1,n1
         rs4p = rs4p*rs4
         a2(i) = a2(i) +orc*(ra(i)/rs4p-a2(i))
      enddo

   enddo
   w = b
   w = (w+ci)/(ci-w)
   cx = real(w)
   cy = aimag(w)
   call cinvrt(a1,b1,n1)
   call cinvrt(a2,b2,n1)

!     print'('' image of octagon-vertex on unit-circle:'')'
!     write(6,62)cx,cy
!     print'('' taylor coefficients, series a1,a2,b1,b2:'')'
!     write(6,63)b
!     do i=1,n1
!      write(6,64)i,a1(i),a2(i),b1(i),b2(i)
!     enddo
!62    format(2(1x,e12.6))
!63    format(4x,'0',13x,4(1x,e12.6))
!64    format(1x,i4,4(1x,e12.6))

end subroutine indoct


!--------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1995
!                   subroutine  otoc1
!   transform to earth-centered cartesians from octagon-1 (otoc1) or from
!   octagon-2 (otoc2).
!
! --> x,y   x and y map-coordinate of point in octagon-1 (otoc1) or
!           octagon-2 (otoc2)
! <-- xe    3-vector of earth-centered cartesian coordinates corresponding
!           to map location (x,y).
!--------------------------------------------------------------------------
subroutine otoc1 ( x, y, xe )
   real(kind=rx), intent(in) :: x, y
   real(kind=rx), intent(out), dimension(:) :: xe   ! (3)
   complex(kind=rx) :: z, w
   real(kind=rx) :: xa, ya, xx, s
   
   z = cmplx(x,y,kind=rx)
   call toct(z,w)
   w = w*sc
   xa = real(w)
   ya = aimag(w)
   xx = xa**2 + ya**2
   s = 2.0/(1.0+xx)
   xe(1) = s-1.0
   xe(2) = s*xa
   xe(3) = s*ya
!   call nmapt(rotm,xe,xe)
   xe = matmul ( transpose(rotm), xe )

end subroutine otoc1

subroutine otoc2 ( x, y, xe )
   real(kind=rx), intent(in) :: x, y
   real(kind=rx), intent(out), dimension(:) :: xe   ! (3)
   complex(kind=rx) :: z, w
   real(kind=rx) :: xa, ya, xx, s

   z = cmplx(-x,y,kind=rx)
   call toct(z,w)
   w = w*sci
   xa = real(w)
   ya = aimag(w)
   xx = xa**2+ya**2
   s = -2.0/(1.0+xx)
   xe(1) = 1.0 + s
   xe(2) = s*xa
   xe(3) = s*ya
!   call nmapt(rotm,xe,xe)
   xe = matmul ( transpose(rotm), xe )

end subroutine otoc2

!--------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1995
!                   subroutine  ctoo
!   transform from earth-centered cartesians to octagon-1 (kmap=1) or to
!   octagon-2 (kmap=-1)
!
! --> xe    earth-centered coordinates (3-vector) of a point
! <-- x,y   map-coordinates of this point
! <-- kmap  map indicator (=1 for octagon-1, =-1 for octagon-2)
!--------------------------------------------------------------------------
subroutine ctoo ( xe, x, y, kmap )
!  Note that Purser's comment of --> xe is wrong
   real(kind=rx), intent(inout), dimension(:) :: xe   ! (3)
   real(kind=rx), intent(out) :: x, y
   integer, intent(out) :: kmap
   logical :: kx, ky, kxy, ks
   complex(kind=rx) :: z, w
   real(kind=rx) :: xa, ya, za, zai, t, dd1, dd2
   
!   call nmap(rotm,xe,xe) 
!  This essentially uses same argument for in and out which is illegal?
   xe = matmul ( rotm, xe )
   xa = xe(2)
   ya = xe(3)
   za = xe(1)
   if ( za > zac ) then
      kmap = 1
      zai = sci/(1.0+za)
      x = xa*zai
      y = ya*zai
   else
      kmap = -1
      zai = sc/(1.0-za)
      x = xa*zai
      y = -ya*zai
   endif
   kx = x < 0.0
   if ( kx ) then
      x = -x
   end if
   ky = y < 0.0
   if ( ky ) then
      y = -y
   end if
   kxy = y > x
   if ( kxy ) then
      t = x
      x = y
      y = t
   endif
   dd1 = x**2 + y**2
   dd2 = (cx-x)**2 + (cy-y)**2
   w = cmplx(x,y,kind=rx)
   ks = dd1 < ssr*dd2
   if ( ks ) then
      w = w**4
      call tay(w,b1,n1,z)
      z = ciq*(-z*ci)**0.25_rx
   else
      w = ci*(w-1.0)/(w+1.0)
      w = w-b
      call tay(w,b2,n1,z)
      z = cir*(ci*z)**sti
      z = ciap - ci*z
   endif
   x = real(z)
   y = aimag(z)
   if ( kxy ) then
      t = x
      x = y
      y = t
   endif
   if ( ky ) then
      y = -y
   end if
   if ( kx ) then
      x = -x
   end if

end subroutine ctoo


!--------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1995
!                   subroutine  toct
!   transform from complex-z in the standard unit-octagon to complex-w in the
!   unit-circle
!----------------------------------------------------------------------------
subroutine toct ( z, w )
   complex(kind=rx), intent(in)  :: z
   complex(kind=rx), intent(out) :: w
   complex(kind=rx) :: zt
   logical :: kx, ky, kxy, kx1, kxy1, ks
   real(kind=rx) :: x, y, t, dd1, dd2
   
   x = real(z)
   y = aimag(z)
   kx = x < 0.0
   if ( kx ) then
      x = -x
   end if
   ky = y < 0.0
   if ( ky ) then
      y = -y
   end if
   kxy  =  y > x
   if ( kxy ) then
      t = x
      x = y
      y = t
   endif
   kx1 = x > 1.0
   kxy1 = x+y > 1.0+a
   if ( kx1 ) then
      x = 2.0 - x
   else if ( kxy1 ) then
      t = x
      x = 1.0 + a - y
      y = 1.0 + a - t
   endif
   dd1 = x**2 + y**2
   dd2 = (1.0-x)**2 + (a-y)**2
   
   zt = cmplx(x,y,kind=rx)
   ks = dd1 < ssr*dd2
   if ( ks ) then
      zt = zt**4
      call tay(zt,a1,n1,w)
      if ( w /= (0.0,0.0) ) then
         w = ciq*(-w*ci)**0.25_rx
      end if
   else
!     jlm: note that DEC alpha has trouble with (0.,0.)**4/3
      zt = a+ci*zt-ci
      if ( zt /= cmplx(0.0,0.0,kind=rx) ) then
         zt = zt**st   !  jlm fix for DEC computers
      end if
      call tay(zt,a2,n1,w)
      w = w+b
      w = (w+ci)/(ci-w)
   endif

   if ( kx1 .or. kxy1 ) then
      w = w/abs(w)**2
   end if

   x = real(w)
   y = aimag(w)
   if ( kxy ) then
      t = x
      x = y
      y = t
   endif
   if ( ky ) then
      y = -y
   end if
   if ( kx ) then
      x = -x
   end if
   w = cmplx(x,y,kind=rx)
end subroutine toct

!-----------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1995
!                   subroutine  ininmap
!  initialize the rotation matrix rot3 needed to transform standard
!  earth-centered cartesian components to the alternative cartesian frame
!  oriented so as to put geographical point (alat0,alon0) on the projection
!  axis.
!-----------------------------------------------------------------------
subroutine ininmap ( alon0, alat0, rot3 )
   real(kind=rx), intent(in) :: alon0, alat0
   real(kind=rx), intent(out), dimension(:,:) :: rot3  ! (3,3)

!   real, dimension(3) :: x, t, u
   real(kind=rx) :: pi, dr, blon0, blat0, clon0, clat0, slon0, slat0

!   print*, "entering ininmap, alat0 = ", alat0
   pi = 4*atan(1.0_rx)
   dr = pi/180.0_rx
   blon0 = dr*alon0
   blat0 = dr*alat0
   clon0 = cos(blon0)
   slon0 = sin(blon0)
   clat0 = cos(blat0)
   slat0 = sin(blat0)
!  Following as in J&J p.122 but for theta use -latitude (jlm)
   rot3(1,1) = clat0*clon0
   rot3(1,2) = clat0*slon0
   rot3(1,3) = slat0
   rot3(2,1) = -slon0
   rot3(2,2) = clon0
   rot3(2,3) = 0.0
   rot3(3,1) = -slat0*clon0
   rot3(3,2) = -slat0*slon0
   rot3(3,3) = clat0
!   print*, "within ininmap, alat0, blat0, clat0, slat0 = ", &
!                            alat0, blat0, clat0, slat0

end subroutine ininmap

!!$subroutine nmap ( rot3, x, t )
!!$   real, dimension(:,:), intent(in) :: rot3  ! (3,3)
!!$   real, dimension(:), intent(in) :: x ! (3)
!!$   real, dimension(:), intent(out) :: t ! (3)
!!$
!!$   t = matmul( rot3, x )
!!$
!!$end subroutine nmap
!!$
!!$subroutine nmapt ( rot3, x, t )
!!$   real, dimension(:,:), intent(in) :: rot3  ! (3,3)
!!$   real, dimension(:), intent(in) :: x ! (3)
!!$   real, dimension(:), intent(out) :: t ! (3)
!!$
!!$   t = matmul(transpose(rot3),x)
!!$
!!$end subroutine nmapt


!------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1994
!                   subroutine  cinvrt
!  compute the taylor series coefficients z for the functional-inverse of
!  the function whose taylor series coefficients are w, for the case where
!  the constant coefficients of both z and w are 0
!
!  --> w    taylor coefficients of original function (starting with linear)
!  <-- z    taylor coefficients of inverse function
!  --> m    number of taylor series coefficients computed
!  --- work workspace array consisting of at least 3*m elements
!--------------------------------------------------------------------------
subroutine cinvrt ( w, z, m )
   integer, intent(in) :: m
   real(kind=rx), intent(in), dimension(:) :: w
   real(kind=rx), intent(out), dimension(:) :: z
   real(kind=rx), dimension(m,3) :: work   ! (m,3)
   integer :: i, j
   real(kind=rx) :: zj
   
   do i = 1, m
      work(i,1) = w(i)
      work(i,2) = 0.0
   enddo
   work(1,2) = 1.0
   do j = 1,m
      zj = work(j,2)/work(j,1)
      z(j) = zj
      do i = j,m
         work(i,2) = work(i,2)-zj*work(i,1)
      enddo
      call conv ( work(:,1), w, work(:,3), m )
      do i = 1,m
         work(i,1) = work(i,3)
      enddo
   enddo
end subroutine cinvrt

!---------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1994
!                   subroutine  conv
!   convolve double-precision series a with b to form c, up to m terms
!   starting with element 1
!
! a,b   --> inputs (convolution factors)
! c     <-- output (convolution product)
! m     --> number of elements of a, b, c
!--------------------------------------------------------------------------
subroutine conv ( a, b, c, m )
   integer, intent(in) :: m
   real(kind=rx), intent(in),  dimension(:) :: a, b
   real(kind=rx), intent(out), dimension(:) :: c
   integer :: i, j, k

   c(1:m) = 0.0
   do i = 1,m
      do j = 1,m-i
         k = i+j
         c(k) = c(k) + a(i)*b(j)
      enddo
   enddo

end subroutine conv

subroutine getrot(rot1)
   real(kind=rx), intent(out), dimension(:,:) :: rot1  ! (3,3)
   rot1 = rotm
end subroutine getrot

end module jimco_m

