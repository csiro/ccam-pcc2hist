module setxyz_m
   use jimcc_m
   use jimco_m
   use xyzinfo_m
   use ind_m
   use precis_m
   use parm_m, only : pi, dtor, omega, rtod
   implicit none
   private

   public :: setxyz, setaxu
   private :: norm, cross3, vecpanel

!  emu and emv have to be contiguous for the out of bounds indexing to work.
!  This goes with arrays declared in map_m
   real(kind=rx), private, save, dimension(:), allocatable, target :: emuv

   real(kind=rx), private, save, dimension(:), allocatable, target :: abx, aby, abz

contains

subroutine setxyz ( il, jl, kl, npanels, ifull, iquad, diag, id, jd,        &
                    rlong0, rlat0, schmidt, schm13, ntang, rearth )

   use indices_m
   integer, intent(in) :: il, jl, kl, npanels, ifull, iquad
   integer, intent(in) :: diag, id, jd
   real(kind=rx), intent(in) :: rlong0, rlat0, schmidt, schm13
   integer, intent(in) :: ntang
!                           ! ntang=1 for tang. vectors by finite diffs
!                           ! ntang=2 for map factors by finite diffs too
   real(kind=rx), intent(in) :: rearth
!     schmidt included
!     *** to fix vecpanel & getting i & j from x,y,z; outfile; veg-files
!     this one bundles in jimcc.f at the end
!     sets up x, y, z on sphere and unit u,v vectors
!     note that x,y,z have been normalized by rearth, the radius of the
!     suffix 6 denotes hex (6)

   integer, parameter, dimension(0:5) :: npan6n = (/ 1, 103, 3, 105, 5, 101 /)
   integer, parameter, dimension(0:5) :: npan6e = (/ 102, 2, 104, 4, 100, 0 /)
   integer, parameter, dimension(0:5) :: npan6w = (/ 5, 105, 1, 101, 3, 103 /) 
   integer, parameter, dimension(0:5) :: npan6s = (/ 104, 0, 100, 2, 102, 4 /)
   integer, parameter, dimension(0:13):: npanno =                         &
      (/ 1, 2, 107, 4, 106, 6, 7, 109, 9, 112, 11, 12, 102, 101 /)  
   integer, parameter, dimension(0:13):: npaneo =                         &
      (/ 103, 3, 4, 105, 5, 110, 108, 8, 10, 11, 100, 113, 13, 0 /)  
   integer, parameter, dimension(0:13):: npanwo =                         &
      (/ 13, 113, 112, 1, 2, 4, 104, 102, 7, 107, 8, 9, 109, 12 /)
   integer, parameter, dimension(0:13):: npanso =                         &
      (/ 110, 0, 1, 100, 3, 103, 5, 6, 106, 8, 105, 10, 11, 111 /)
   integer, dimension(:), allocatable :: npann, npane, npanw, npans

   integer :: idjd, iq, n, ii, i, j, i0, j0, n0, in0, jn0, nn0, is0,  &
              js0, ns0, ie0, je0, ne0, iw0, jw0, nw0, inn0, jnn0, nnn0, iss0,&
              jss0, nss0, iee0, jee0, nee0, iww0, jww0, nww0, jj, jx, iy, iyy,&
              iqn, iq11, iq12, iq13, iq22, iq32 
   real(kind=rx) , dimension(3,3) :: rotpole 
   real(kind=rx) :: dsfact, xx, yy, alf, temp, xin, yin, zin, dot, &
         eps, dx2, dy2, sumwts, coslong, sinlong, coslat, sinlat, zz, demu11, &
         demv11, demu21, demv21, demu31, demv31 

   integer, dimension(ifull) :: i_nw, i_sw, i_es, i_ws 
   real(kind=rx), dimension(ifull) :: axx, ayy, azz, bxx, byy, bzz
   real(kind=rx), dimension(iquad,iquad) :: em4, ax4, ay4, az4, zz4
   real(kind=rx), dimension(ifull) :: cosa

   integer :: ijk
!  If abs(cos(lat)) < polelim than assume point is exactly at pole.
   real(kind=rx), parameter :: polelim = 10*epsilon(1.0)

!  Allocate all the public arrays
   allocate ( xx4(iquad,iquad), yy4(iquad,iquad) )
   allocate ( i_n(ifull),  i_s(ifull), i_w(ifull), i_e(ifull), i_nn(ifull),   &
              i_ss(ifull), i_ww(ifull), i_ee(ifull), i_ne(ifull),             &
              i_se(ifull), i_en(ifull), i_wn(ifull), i_wu(ifull),             &
              i_sv(ifull), i_wu2(ifull), i_sv2(ifull), i_eu2(ifull),          &
              i_nv2(ifull), i_ev2(ifull), i_nu2(ifull), i_eu(ifull),          &
              i_nv(ifull) )
   allocate ( lwws(0:npanels), lws(0:npanels), lwss(0:npanels),            &
              les(0:npanels), lees(0:npanels), less(0:npanels),            &
              lwwn(0:npanels), lwnn(0:npanels), leen(0:npanels),           &
              lenn(0:npanels), lsww(0:npanels), lsw(0:npanels),            &
              lssw(0:npanels), lsee(0:npanels), lsse(0:npanels),           &
              lnww(0:npanels), lnw(0:npanels), lnnw(0:npanels),            &
              lnee(0:npanels), lnne(0:npanels) )
   allocate ( em(ifull),                                                   &
              f(ifull), fu(ifull), fv(ifull),                              &
              dmdx(ifull), dmdy(ifull), dmdxv(ifull), dmdyu(ifull) )
   allocate ( x(ifull), y(ifull), z(ifull), rlat(ifull), rlong(ifull),     &
              wts(ifull) )

!  Pointer trick to mimic the effect of common.
   allocate ( emuv(-ifull+1:ifull) )
   emu => emuv(-ifull+1:ifull)
   emv => emuv
   allocate ( abx(-ifull+1:ifull) )
   ax => abx(-ifull+1:ifull)
   bx => abx
   allocate ( aby(-ifull+1:ifull) )
   ay => aby(-ifull+1:ifull)
   by => aby
   allocate ( abz(-ifull+1:ifull) )
   az => abz(-ifull+1:ifull)
   bz => abz

   ijk = il*jl*kl

!  When using the ifull notation: i_n, i_e, i_w and i_s give the
!  indices for the n, e, w, s neighbours respectively
!  a, b denote unit vectors in the direction of x, y (e & n) respectively

   idjd = id + il*(jd - 1) 
   do iq = 1, ifull 
      i_n(iq) = iq + il 
      i_s(iq) = iq - il 
      i_e(iq) = iq + 1 
      i_w(iq) = iq - 1 
   end do
 
   allocate ( npann(0:npanels), npane(0:npanels), npanw(0:npanels), &
              npans(0:npanels) )
   if (npanels == 5) then 
      npann = npan6n 
      npane = npan6e 
      npanw = npan6w 
      npans = npan6s 
   else
      npann = npanno 
      npane = npaneo 
      npanw = npanwo 
      npans = npanso 
   endif
 
   do n = 0, npanels 
!     print *,"ina il/2,n ",i_n(ind(il/2,il,n)),n
      if (npann(n) < 100) then 
         do ii = 1, il 
            i_n(ind(ii,il,n)) = ind(ii,1,npann(n)) 
         end do
      else 
         do ii = 1, il 
            i_n(ind(ii,il,n)) = ind(1,il + 1 - ii,npann(n)-100) 
         end do
      endif
!     print *,"inb il/2,n ",i_n(ind(il/2,il,n)),n
!     print *,"iea il/2,n ",i_e(ind(il,il/2,n)),n
      if (npane(n) < 100) then 
         do ii = 1, il 
            i_e(ind(il,ii,n)) = ind(1,ii,npane(n)) 
         end do
      else 
         do ii = 1, il 
            i_e(ind(il,ii,n)) = ind(il + 1 - ii,1,npane(n)-100) 
         end do
      endif
!     print *,"ieb il/2,n ",i_e(ind(il,il/2,n)),n
!     print *,"iwa il/2,n ",i_w(ind(1,il/2,n)),n
      if (npanw(n) < 100) then 
         do ii = 1, il 
            i_w(ind(1,ii,n)) = ind(il,ii,npanw(n)) 
         end do
      else 
         do ii = 1, il 
            i_w(ind(1,ii,n)) = ind(il + 1 - ii,il,npanw(n)-100) 
         end do
      endif
!     print *,"iwb il/2,n ",i_w(ind(1,il/2,n)),n
!     print *,"isa il/2,n ",is(ind(il/2,1,n)),n
      if (npans(n) < 100) then 
         do ii = 1, il 
            i_s(ind(ii,1,n)) = ind(ii,il,npans(n)) 
         end do
      else 
         do ii = 1, il 
            i_s(ind(ii,1,n)) = ind(il,il + 1 - ii,npans(n)-100) 
         end do
      endif
!     print *,"isb il/2,n ",is(ind(il/2,1,n)),n
   end do! n loop 
 
   i_nn = i_n(i_n) 
   i_ss = i_s(i_s) 
   i_ee = i_e(i_e) 
   i_ww = i_w(i_w) 
   i_ne = i_n(i_e) 
   i_se = i_s(i_e) 
   i_en = i_e(i_n) 
   i_wn = i_w(i_n) 
   i_wu = i_w 
   i_sv = i_s 

!  Use for unstaggered u,v in hordifg 
   i_wu2 = i_w      
   i_sv2 = i_s         
   i_eu2 = i_e         
   i_nv2 = i_n      

!  Use for staguv3
   i_eu = i_e
   i_nv = i_n

!  Following are extras not needed in model, just here for bdy values
   i_nw = i_n(i_w)               ! in temporary arrays 
   i_sw = i_s(i_w)               ! in temporary arrays 
   i_es = i_e(i_s)               ! in temporary arrays 
   i_ws = i_w(i_s)               ! in temporary arrays 

   do n = 0, npanels 
!     following treats unusual panel boundaries
!     print *,"ina il/2,n ",i_n(ind(il/2,il,n)),n
      if (npann(n) >= 100) then 
         do i = 1, il 
            iq = ind(i,il,n) 
            i_nn(iq) = i_e(i_n(iq)) 
            i_en(iq) = i_s(i_n(iq)) 
            i_wn(iq) = i_n(i_n(iq)) 
            i_nv2(iq) = i_n(iq) - ifull     ! converts 2D v array into u array 
            i_nv(iq) = i_n(iq) - ijk        ! converts 3D v array into u array 
         end do
      endif
!     print *,"inb il/2,n ",i_n(ind(il/2,il,n)),n
!     print *,"iea il/2,n ",i_e(ind(il,il/2,n)),n
      if (npane(n) >= 100) then 
         do j = 1, il 
            iq = ind(il,j,n) 
            i_ee(iq) = i_n(i_e(iq)) 
            i_ne(iq) = i_w(i_e(iq)) 
            i_se(iq) = i_e(i_e(iq)) 
            i_eu2(iq) = i_e(iq) + ifull     ! converts 2D u array into v array 
            i_eu(iq) = i_e(iq) + ijk        ! converts 3D u array into v array 
         end do
      endif
!     print *,"ieb il/2,n ",i_e(ind(il,il/2,n)),n
!     print *,"iwa il/2,n ",i_w(ind(1,il/2,n)),n
      if (npanw(n) >= 100) then 
         do j = 1, il 
            iq = ind(1,j,n) 
            i_ww(iq) = i_s(i_w(iq)) 
            i_nw(iq) = i_w(i_w(iq))         ! in temporary arrays 
            i_sw(iq) = i_e(i_w(iq))         ! in temporary arrays 
            i_wu2(iq) = i_w(iq) + ifull     ! converts 2D u array into v array 
            i_wu(iq) = i_w(iq) + ijk        ! converts 3D u array into v array 
         end do! j loop 
      endif
!     print *,"iwb il/2,n ",i_w(ind(1,il/2,n)),n
!     print *,"isa il/2,n ",i_s(ind(il/2,1,n)),n
      if (npans(n) >= 100) then
         do i = 1, il 
            iq = ind(i,1,n) 
            i_ss(iq) = i_w(i_s(iq)) 
            i_es(iq) = i_s(i_s(iq))         ! in temporary arrays 
            i_ws(iq) = i_n(i_s(iq))         ! in temporary arrays 
            i_sv2(iq) = i_s(iq) - ifull     ! converts 2D v array into u array 
            i_sv(iq) = i_s(iq) - ijk        ! converts 3D v array into u array 
         end do
      end if
!     print *,"isb il/2,n ",i_s(ind(il/2,1,n)),n
   end do! n loop 
 
!   print *,"lsw a  ",lsw
!   print *,"lnw a  ",lnw
!   print *,"lws a  ",lws
!   print *,"les a  ",les
!   print *,"leen a  ",leen
!   print *,"lenn a  ",lenn
!   print *,"lwwn a  ",lwwn
!   print *,"lwnn a  ",lwnn
!   print *,"lsee a  ",lsee
!   print *,"lsse a  ",lsse
!   print *,"lnee a  ",lnee
!   print *,"lnne a  ",lnne
!   print *,"lsww a  ",lsww
!   print *,"lssw a  ",lssw
!   print *,"lnww a  ",lnww
!   print *,"lnnw a  ",lnnw
!   print *,"lwws a  ",lwws
!   print *,"lwss a  ",lwss
!   print *,"lees a  ",lees
!   print *,"less a  ",less
   do n = 0, npanels 
      lsw(n) = i_sw(ind(1,1,n)) 
      lnw(n) = i_nw(ind(1,il,n)) 
      lws(n) = i_ws(ind(1,1,n)) 
      les(n) = i_es(ind(il,1,n)) 
      leen(n) = i_ee(i_n(ind(il,il,n))) 
      lenn(n) = i_en(i_n(ind(il,il,n))) 
      lwnn(n) = i_wn(i_n(ind(1,il,n))) 
      lsee(n) = i_se(i_e(ind(il,1,n))) 
      lnee(n) = i_ne(i_e(ind(il,il,n))) 
      lnne(n) = i_nn(i_e(ind(il,il,n))) 
      lsww(n) = i_sw(i_w(ind(1,1,n))) 
      lssw(n) = i_ss(i_w(ind(1,1,n))) 
      lnww(n) = i_nw(i_w(ind(1,il,n))) 
      lwws(n) = i_ww(i_s(ind(1,1,n))) 
      lwss(n) = i_ws(i_s(ind(1,1,n))) 
      less(n) = i_es(i_s(ind(il,1,n))) 
      lwwn(n) = i_ww(i_n(ind(1,il,n))) 
      lsse(n) = i_ss(i_e(ind(il,1,n))) 
      lnnw(n) = i_nn(i_w(ind(1,il,n))) 
      lees(n) = i_ee(i_s(ind(il,1,n))) 
      if ( npann(n) >= 100 ) then 
         leen(n) = i_ss(i_n(ind(il,il,n))) 
         lenn(n) = i_se(i_n(ind(il,il,n))) 
         lwnn(n) = i_ne(i_n(ind(1,il,n))) 
         lwwn(n) = i_nn(i_n(ind(1,il,n))) 
      endif! (npann(n).ge.100) 
      if ( npane(n) >= 100 ) then 
         lsee(n) = i_en(i_e(ind(il,1,n))) 
         lnee(n) = i_wn(i_e(ind(il,il,n))) 
         lnne(n) = i_ww(i_e(ind(il,il,n))) 
         lsse(n) = i_ee(i_e(ind(il,1,n))) 
      endif! (npane(n).ge.100) 
      if ( npanw(n) >= 100 ) then 
         lsww(n) = i_es(i_w(ind(1,1,n))) 
         lssw(n) = i_ee(i_w(ind(1,1,n))) 
         lnww(n) = i_ws(i_w(ind(1,il,n))) 
         lnnw(n) = i_ww(i_w(ind(1,il,n))) 
      endif! (npanw(n).ge.100) 
      if ( npans(n) >= 100 ) then
         lwws(n) = i_nn(i_s(ind(1,1,n))) 
         lwss(n) = i_nw(i_s(ind(1,1,n))) 
         less(n) = i_sw(i_s(ind(il,1,n))) 
         lees(n) = i_ss(i_s(ind(il,1,n))) 
      end if
   end do                                     ! n loop 
!   print *,"lsw b  ",lsw
!   print *,"lnw b  ",lnw
!   print *,"lws b  ",lws
!   print *,"les b  ",les
!   print *,"leen b  ",leen
!   print *,"lenn b  ",lenn
!   print *,"lwwn b  ",lwwn
!   print *,"lwnn b  ",lwnn
!   print *,"lsee b  ",lsee
!   print *,"lsse b  ",lsse
!   print *,"lnee b  ",lnee
!   print *,"lnne b  ",lnne
!   print *,"lsww b  ",lsww
!   print *,"lssw b  ",lssw
!   print *,"lnww b  ",lnww
!   print *,"lnnw b  ",lnnw
!   print *,"lwws b  ",lwws
!   print *,"lwss b  ",lwss
!   print *,"lees b  ",lees
!   print *,"less b  ",less
 
   if (diag == 3) then 
      do n = 0, npanels 
         do j = 1, il 
            do i = 1, il 
               iq = ind(i,j,n) 
               call indv (iq, i0, j0, n0) 
               call indv (i_n(iq), in0, jn0, nn0) 
               call indv (i_s(iq), is0, js0, ns0) 
               call indv (i_e(iq), ie0, je0, ne0) 
               call indv (i_w(iq), iw0, jw0, nw0) 
               call indv (i_nn(iq), inn0, jnn0, nnn0) 
               call indv (i_ss(iq), iss0, jss0, nss0) 
               call indv (i_ee(iq), iee0, jee0, nee0) 
               call indv (i_ww(iq), iww0, jww0, nww0) 
               write(unit=91,fmt="(9(i4,i2,i2))")                           &
                    i0, j0, n0, in0, jn0, nn0, is0, js0, ns0, ie0, je0,     &
                    ne0, iw0, jw0, nw0, inn0, jnn0, nnn0, iss0, jss0, nss0, &
                    iee0, jee0, nee0, iww0, jww0, nww0 
            end do! i loop 
         end do! j loop 
      end do! n loop 
   endif! (diag.eq.3) 
 
!----------------------------------------------------------------------------
!  calculate grid information using quadruple resolution grid
   if (npanels == 5) then 
      call jimcc ( iquad, xx4, yy4, em4, ax4, ay4, az4 )
      dsfact = 4*il/(2.0*pi)                   ! con-cube 
      ds = rearth/dsfact 
      if ( diag /= 0 ) then
         print *, "xx4 first & last ", xx4(1,1), xx4(iquad,iquad) 
         print *, "xx4 (5,5),(7,7),(9,9) ", xx4(5,5), xx4(7,7), xx4(9,9) 
         print *, "yy4 first & last ", yy4(1,1), yy4(iquad,iquad) 
         print *, "yy4 (5,5),(7,7),(9,9) ", yy4(5,5), yy4(7,7), yy4(9,9) 
         print *, "xx4, yy4 central", xx4(2*il+1,2*il+1), yy4(2*il+1,2*il+1) 
      end if
 
!     extend em4 to uppermost i and j rows
      em4(iquad,:4*il) = em4(1,:4*il) 
      em4(:4*il,iquad) = em4(:4*il,1) 

!     Move this out so it vectorises better.
!     The corner values are zero. They're not actually used but reset them
!     here to avoid a divide by zero.
      em4(1,1)         = 1.0
      em4(iquad,1)     = 1.0
      em4(1,iquad)     = 1.0
      em4(iquad,iquad) = 1.0
      em4 = pi/(2.0*em4) 
 
      do j = 1, il 
         do i = 1, il 
!           average Purser em is pi/2
            iq = ind(i,j,0)
!           Step of il*il loops over the faces
!!$            em(iq::il*il) = pi/(2.0*em4(4*i-1,4*j-1)) 
!!$!           em(i,j:5*il+j:il) = pi/(2.0*em4(4*i-1,4*j-1)) 
!!$            emu(iq::il*il) = pi/(2.0*em4(4*i+1,4*j-1)) 
!!$            emv(iq::il*il) = pi/(2.0*em4(4*i-1,4*j+1)) 
            em(iq::il*il) = em4(4*i-1,4*j-1)
            emu(iq::il*il) = em4(4*i+1,4*j-1)
            emv(iq::il*il) = em4(4*i-1,4*j+1)

            xx = xx4(4*i-1,4*j-1) 
            yy = yy4(4*i-1,4*j-1) 
            ax(iq) = ax4(4*i-1,4*j-1) 
            ay(iq) = ay4(4*i-1,4*j-1) 
            az(iq) = az4(4*i-1,4*j-1) 
 
!         set up x0, y0, z0 coords on cube -1 to 1
!         to save space have equivalenced x,x0  etc
            x(ind(i,j,0)) = 1.0
            y(ind(i,j,0)) = xx 
            z(ind(i,j,0)) = yy 
            x(ind(i,j,3)) = -1.0 
            z(ind(i,j,3)) = -xx 
            y(ind(i,j,3)) = -yy 
             
            x(ind(i,j,1)) = -yy 
            y(ind(i,j,1)) = xx 
            z(ind(i,j,1)) = 1.0 
            y(ind(i,j,4)) = -yy 
            x(ind(i,j,4)) = xx 
            z(ind(i,j,4)) = -1.0 
             
            x(ind(i,j,2)) = -yy 
            y(ind(i,j,2)) = 1.0 
            z(ind(i,j,2)) = -xx 
            z(ind(i,j,5)) = yy 
            y(ind(i,j,5)) = -1.0 
            x(ind(i,j,5)) = xx 
         end do! i loop 
      end do! j loop 
      if ( diag /= 0 ) then
         print *, "em (1,1) & (2,2) ", em(1), em(ind(2,2,0)) 
         print *, "ax6 (1,1,0) & (2,2,0) ", ax(ind(1,1,0)), ax(ind(2,2,0))
         print *, "ay6 (1,1,0) & (2,2,0) ", ay(ind(1,1,0)), ay(ind(2,2,0)) 
         print *, "az6 (1,1,0) & (2,2,0) ", az(ind(1,1,0)), az(ind(2,2,0))
      end if
!     x, y, z are coords on sphere  -1 to 1 
      call norm ( x, y, z ) 
   endif! (npanels.eq.5) 
 
   if (npanels == 13) then 
      if ( ntang /= 2 ) then 
         print*, " Error - octagon requires ntang=2 "
         stop
      end if
      call jimco ( il, iquad, xx4, yy4, zz4 )
      dsfact = (3.0 + 3.0*sqrt(2.0_rx))*il/(2.0*pi)  ! con-octagon 
      ds = rearth/dsfact 
!     Silicon Graphics has trouble with (1,1) so set here
      xx4(1,1) = 0.0 
      yy4(1,1) = 0.0 
      em4(1,1) = 1.0 
      if ( diag /= 0 ) then
         print *, "xyz (1,1) ", xx4(1,1), yy4(1,1), zz4(1,1) 
         print *, "xyz (iquad,1) ", xx4(iquad,1), yy4(iquad,1), zz4(iquad,1) 
         print *, "xyz (1,iquad) ", xx4(1,iquad), yy4(1,iquad), zz4(1,iquad) 
         print *, "xyz (iquad,iquad) ", xx4(iquad,iquad), yy4(iquad,iquad), &
                                        zz4(iquad,iquad) 
         print *, "xyz (iquad/2,iquad/2) ", xx4(iquad/2,iquad/2), &
            yy4(iquad/2,iquad/2), zz4(iquad/2,iquad/2) 
      end if

      do j=il/2+1,il
         jj=4*j -2*il -1
         do i=(il+1)/2,il        ! NE corner of panel 2
            ii=4*i -2*il -1
            x(ind(i,j,2)) = xx4(ii,jj)
            y(ind(i,j,2)) = yy4(ii,jj)
            z(ind(i,j,2)) = zz4(ii,jj)
!           define SE corner of panel 2
            x(ind(i,il+1-j,2)) = -x(ind(i,j,2))
            y(ind(i,il+1-j,2)) =  y(ind(i,j,2))
            z(ind(i,il+1-j,2)) =  z(ind(i,j,2))
         enddo  ! i loop
         do i=1,il               ! N half of panel 4
            ii=4*(i+il) -2*il -1
            x(ind(i,j,4)) = xx4(ii,jj)
            y(ind(i,j,4)) = yy4(ii,jj)
            z(ind(i,j,4)) = zz4(ii,jj)
!           define S half of panel 4
            x(ind(i,il+1-j,4)) = -x(ind(i,j,4))
            y(ind(i,il+1-j,4)) =  y(ind(i,j,4))
            z(ind(i,il+1-j,4)) =  z(ind(i,j,4))
         enddo! i loop
      enddo! j loop
      do jx=1,il
         jj=4*(jx+il) -2*il -1
         do iy=jx,il             ! SW half of panel 6
            iyy=il+1-iy
            ii=4*(iyy+il) -2*il -1
            x(ind(jx,iy,6)) = xx4(ii,jj)
            y(ind(jx,iy,6)) = yy4(ii,jj)
            z(ind(jx,iy,6)) = zz4(ii,jj)
!           define other half of panel 6
            x(ind(iy,jx,6)) =  x(ind(jx,iy,6))
            y(ind(iy,jx,6)) =  y(ind(jx,iy,6))
            z(ind(iy,jx,6)) = -z(ind(jx,iy,6))
         enddo  ! iyy loop
         do i=1,il/2
!           define W half of panel 2
            x(ind(i,jx,2)) =  x(ind(il+1-i,jx,2))
            y(ind(i,jx,2)) = -y(ind(il+1-i,jx,2))
            z(ind(i,jx,2)) =  z(ind(il+1-i,jx,2))
         enddo  ! i loop
      enddo! jx loop
!!$c                      nn=2
!!$c                      print *,'x for panel ',nn
!!$c                      do j=il,1,-1
!!$c                       print '(i3,30f6.3)',j,(x(ind(i,j,nn),i=1,il)
!!$c                      enddo
!!$c                      print *,'y for panel ',nn
!!$c                      do j=il,1,-1
!!$c                       print '(i3,30f6.3)',j,(y(ind(i,j,nn),i=1,il)
!!$c                      enddo
!!$c                      print *,'z for panel ',nn
!!$c                      do j=il,1,-1
!!$c                       print '(i3,30f6.3)',j,(z(ind(i,j,nn),i=1,il)
!!$c                      enddo
!       define remaining 11 panels
        do j=1,il
           do i=1,il
              x(ind(i,j,10)) =  y(ind(i,j,2))            ! C2
              y(ind(i,j,10)) =  x(ind(i,j,2))
              z(ind(i,j,10)) = -z(ind(i,j,2))
              x(ind(i,j,3))  =  y(ind(i,j,6))             ! SE
              y(ind(i,j,3))  = -x(ind(i,j,6))
              z(ind(i,j,3))  =  z(ind(i,j,6))
              x(ind(i,j,9))  = -y(ind(i,j,6))             ! NW
              y(ind(i,j,9))  =  x(ind(i,j,6))
              z(ind(i,j,9))  =  z(ind(i,j,6))
              x(ind(i,j,13)) = -x(ind(i,j,6))             ! SW
              y(ind(i,j,13)) = -y(ind(i,j,6))
              z(ind(i,j,13)) =  z(ind(i,j,6))
              x(ind(i,j,1))  =  y(ind(il+1-j,i,4))       ! SC1
              y(ind(i,j,1))  = -x(ind(il+1-j,i,4))
              z(ind(i,j,1))  =  z(ind(il+1-j,i,4))
              x(ind(i,j,7))  = -y(ind(i,j,4))            ! NC1
              y(ind(i,j,7))  =  x(ind(i,j,4))
              z(ind(i,j,7))  =  z(ind(i,j,4))
              x(ind(i,j,12)) = -x(ind(il+1-j,i,4))       ! WC1
              y(ind(i,j,12)) = -y(ind(il+1-j,i,4))
              z(ind(i,j,12)) =  z(ind(il+1-j,i,4))
              x(ind(i,j,5))  =  x(ind(il+1-i,j,4))       ! EC2
              y(ind(i,j,5))  =  y(ind(il+1-i,j,4))
              z(ind(i,j,5))  = -z(ind(il+1-i,j,4))
              x(ind(i,j,11)) =  x(ind(j,il+1-i,4))       ! WC2
              y(ind(i,j,11)) = -y(ind(j,il+1-i,4))
              z(ind(i,j,11)) = -z(ind(j,il+1-i,4))
              x(ind(i,j,8))  = -y(ind(il+1-i,j,4))     ! NC2
              y(ind(i,j,8))  =  x(ind(il+1-i,j,4))
              z(ind(i,j,8))  = -z(ind(il+1-i,j,4))
              x(ind(i,j,0))  =  y(ind(j,i,4))          ! SC2
              y(ind(i,j,0)) = -x(ind(j,i,4))
              z(ind(i,j,0)) = -z(ind(j,i,4))
           enddo  ! i loop
        enddo   ! j loop

!       for conformal octagon, save only xx4 and yy4
!           as absolute values and as if schmidt=.5
!           flip xx and yy (fix up later in jimco & above)
      alf = (1.0 - schm13**2)/(1.0 + schm13**2) 
      do jj = 1, iquad 
      do ii = 1, iquad 
         temp = abs(xx4(ii,jj))*schm13*(1.0 + alf)/(1.0 + alf*zz4(ii,jj)) 
         xx4(ii,jj) = abs(yy4(ii,jj))*schm13*(1.0 + alf) /             &
                         (1.0 + alf*zz4(ii,jj)) 
         yy4(ii,jj) = temp 
      end do! ii loop 
      end do! jj loop 
   endif! (npanels.eq.13) 
 
   if ( diag /= 0 ) then
      print *, "basic grid length ds =", ds 
   end if
   if (schmidt /= 1.0) then 
      alf = (1.0 - schmidt**2)/(1.0 + schmidt**2) 
      if ( diag /= 0 ) then
         print *, "doing schmidt with schmidt,alf: ", schmidt, alf 
      end if
      do iq = 1, ifull 
         xin = x(iq) 
         yin = y(iq) 
         zin = z(iq) 
         x(iq) = xin*schmidt*(1.0 + alf)/(1.0 + alf*zin) 
         y(iq) = yin*schmidt*(1.0 + alf)/(1.0 + alf*zin) 
         z(iq) = (alf + zin)/(1.0 + alf*zin) 
         if (ntang /= 2) then 
            em(iq) = em(iq)*schmidt*(1.0 + alf*zin)/(1.0 - alf) 
         end if
         do n = 0, npanels 
            iqn = ind((il + 1)/2,(il + 1)/2,n) 
            if (iq == iqn .and. diag /= 0 ) then
               print*, "After Schmidt at centre of face n:", n 
               write(unit=*,fmt="(a,3f7.3,2x,3f7.3)")  &
                    " xin,yin,zin,x,y,z ", xin, yin, zin, x(iq), y(iq), z(iq) 
            end if
         end do! n loop 
      end do! iq loop 
 
      if (ntang /= 2) then 
!        With schmidt must average em to get emu & emv
!        Note that explicit ranges required here because of the pointer trick.
         emu(1:ifull) = 0.5_rx * ( em + em(i_e) ) 
         emv(1:ifull) = 0.5_rx * ( em + em(i_n) ) 
      endif
   endif!  (schmidt.ne.1.0) 
 
!!$      if (diag == 2) call printp ("x   ", x) 
!!$      if (diag == 2) call printp ("y   ", y) 
!!$      if (diag == 2) call printp ("z   ", z) 
 
!  Set up vectors in direction of u and v
   if (ntang == 0) then 
!     define x-vectors on panels 1:5 from panel 0 
!     Pointer trick requires explicit dimensions here
      call vecpanel ( ax(1:ifull), ay(1:ifull), az(1:ifull) )
!     define y-vectors
      call cross3 ( bx(1:ifull), by(1:ifull), bz(1:ifull), x, y, z, &
                    ax(1:ifull), ay(1:ifull), az(1:ifull)) 
   else 
      do iq = 1, ifull 
!        first guess tang vectors by finite differences
         ax(iq) = x(i_e(iq)) - x(i_w(iq)) 
         ay(iq) = y(i_e(iq)) - y(i_w(iq)) 
         az(iq) = z(i_e(iq)) - z(i_w(iq)) 
         bx(iq) = x(i_n(iq)) - x(i_s(iq)) 
         by(iq) = y(i_n(iq)) - y(i_s(iq)) 
         bz(iq) = z(i_n(iq)) - z(i_s(iq)) 
         if ( diag /= 0 .and. ( iq == idjd .or. iq == i_n(idjd) ) ) then
            print *, "first guess values for ax,bx" 
            print *, "iq,ax,ay,az", iq, ax(iq), ay(iq), az(iq) 
            print *, "iq,bx,by,bz", iq, bx(iq), by(iq), bz(iq) 
            print *, "iq,x,y,z   ", iq, x(iq), y(iq), z(iq) 
            print *, "iq,x,y,z n ", iq, x(i_n(iq)), y(i_n(iq)), z(i_n(iq)) 
            print *, "iq,x,y,z e ", iq, x(i_e(iq)), y(i_e(iq)), z(i_e(iq)) 
            print *, "iq,x,y,z w ", iq, x(i_w(iq)), y(i_w(iq)), z(i_w(iq)) 
            print *, "iq,x,y,z s ", iq, x(i_s(iq)), y(i_s(iq)), z(i_s(iq)) 
         end if
      end do! iq loop 

!     Form axx and bxx tangential to the sphere
      call cross3 (axx, ayy, azz, bx(1:ifull), by(1:ifull), bz(1:ifull), &
                  x, y, z) 
      call cross3 (bxx, byy, bzz, x, y, z, ax(1:ifull), ay(1:ifull), &
                   az(1:ifull)) 
      call norm ( axx, ayy, azz ) 
      call norm ( bxx, byy, bzz ) 
      do iq = 1, ifull 
!        make sure they are perpendicular & normalize
         dot = axx(iq)*bxx(iq) + ayy(iq)*byy(iq) + azz(iq)*bzz(iq) 
         eps = -dot/(1.0 + sqrt(1.0 - dot*dot)) 
         ax(iq) = axx(iq) + eps*bxx(iq) 
         ay(iq) = ayy(iq) + eps*byy(iq) 
         az(iq) = azz(iq) + eps*bzz(iq) 
         bx(iq) = bxx(iq) + eps*axx(iq) 
         by(iq) = byy(iq) + eps*ayy(iq) 
         bz(iq) = bzz(iq) + eps*azz(iq) 
      end do
      call norm ( ax(1:ifull), ay(1:ifull), az(1:ifull) ) 
      call norm ( bx(1:ifull), by(1:ifull), bz(1:ifull) ) 

      if (ntang == 2) then 
         do iq = 1, ifull 
!          calculate inverse of emu & emv first
            dx2 = (x(i_e(iq))-x(iq))**2 + (y(i_e(iq))-y(iq))**2 +            &
                  (z(i_e(iq))-z(iq))**2 
            emu(iq) = sqrt(dx2)*(1.0 + dx2/24.0)*dsfact 
            dy2 = (x(i_n(iq))-x(iq))**2 + (y(i_n(iq))-y(iq))**2 +            &
                  (z(i_n(iq))-z(iq))**2 
            emv(iq) = sqrt(dy2)*(1.0 + dy2/24.0)*dsfact 
         end do! iq loop 

! based on inverse values of emu & emv
!        Note that the explicit ranges 1:ifull are required here because
!        the pointer trick has extended the arrays.
         em = 4.0 / ( emu(i_wu2) + emu(1:ifull) + emv(i_sv2) + emv(1:ifull) ) 
!        experimental option follows - only tiniest difference for il=20
!        em(iq) = 2.0 / sqrt( (emu(i_wu2)+emu(1:ifull))*(emv(isv2)+emv(1:ifull)) )
!        emu = 1.0 / emu
!        emv = 1.0 / emv
         emuv = 1.0 / emuv
      endif! (ntang.eq.2) 
   endif! (ntang.eq.0) 
   if ( diag /= 0 ) then
      do iq = il - 2, il 
         print *, "iq,em,emu,emv", iq, em(iq), emu(iq), emv(iq) 
      end do! iq loop 
      if (id<=il .and. jd<=jl) then 
         iq = id + il*(jd - 1) 
         print *, "values at idjd" 
         print *, "iq,ax,ay,az", iq, ax(iq), ay(iq), az(iq) 
         print *, "iq,bx,by,bz", iq, bx(iq), by(iq), bz(iq) 
         iq = i_n(id+il*(jd-1)) 
         print *, "values at i_n(idjd)" 
         print *, "iq,ax,ay,az", iq, ax(i_n(iq)), ay(i_n(iq)), az(i_n(iq)) 
         print *, "iq,bx,by,bz", iq, bx(i_n(iq)), by(i_n(iq)), bz(i_n(iq)) 
         print *, "values at i_e(idjd)" 
         print *, "iq,ax,ay,az", iq, ax(i_e(iq)), ay(i_e(iq)), az(i_e(iq)) 
         print *, "iq,bx,by,bz", iq, bx(i_e(iq)), by(i_e(iq)), bz(i_e(iq)) 
         print *, "values at i_w(idjd)" 
         print *, "iq,ax,ay,az", iq, ax(i_w(iq)), ay(i_w(iq)), az(i_w(iq)) 
         print *, "iq,bx,by,bz", iq, bx(i_w(iq)), by(i_w(iq)), bz(i_w(iq)) 
         print *, "values at i_s(idjd)" 
         print *, "iq,ax,ay,az", iq, ax(i_s(iq)), ay(i_s(iq)), az(i_s(iq)) 
         print *, "iq,bx,by,bz", iq, bx(i_s(iq)), by(i_s(iq)), bz(i_s(iq)) 
      endif
   end if
 
!  Calculate approx areas around each grid point
!  Just used for error diagnostics
!  Now use 1/(em**2) to cope with schmidt, rotated and ocatagon coordinates
   wts = 1.0/em**2 
   sumwts = sum(wts) 
!  cosa is dot product of unit vectors
!  *** only useful as diagnostic for gnew
   cosa = ax(1:ifull)*bx(1:ifull) + ay(1:ifull)*by(1:ifull) +  &
          az(1:ifull)*bz(1:ifull) 
!       call printp("dx  ",dx)
!       call printp("dy  ",dy)
!!$   if (diag == 2) then 
!!$      call printp ("cosa", cosa) 
!!$   endif
   if ( diag /= 0 ) then
      print *, "sumwts/ifull ", sumwts/ifull     ! ideally equals 4*pi ?? 
   end if
 
!  Previously calculates approx areas around each grid point using modulus of
!  cross-products; used for error diagnostics
!     do j=1,il
!      do i=1,il
!       iq=ind(i,j,0)
!       wts(iq)=.5*(crossmod(iq,4*i-3,4*j+1,4*i+1,4*j+1)  ! nw,ne  with
!    .             +crossmod(iq,4*i+1,4*j+1,4*i+1,4*j-3)  ! ne,se  with
!    .             +crossmod(iq,4*i+1,4*j-3,4*i-3,4*j-3)  ! se,sw  with
!    .             +crossmod(iq,4*i-3,4*j-3,4*i-3,4*j+1)) ! sw,nw  with
!       sumwts=sumwts+wts(iq)
!       do n=1,5
!        wts(iq+n*il*il)=wts(iq)
!       enddo  ! n loop
!      enddo   ! i loop
!     enddo    ! j loop
!     print *,"sumwts ",sumwts  ! ideally equals 4*pi/6, 2.09440
 
   coslong = cos(rlong0*dtor) 
   sinlong = sin(rlong0*dtor) 
   coslat = cos(rlat0*dtor) 
   sinlat = sin(rlat0*dtor) 
   rotpole(1,1) = coslong*sinlat 
   rotpole(1,2) = -sinlong 
   rotpole(1,3) = coslong*coslat 
   rotpole(2,1) = sinlong*sinlat 
   rotpole(2,2) = coslong 
   rotpole(2,3) = sinlong*coslat 
   rotpole(3,1) = -coslat 
   rotpole(3,2) = 0.0
   rotpole(3,3) = sinlat 
 
   if ( diag /= 0 ) then
      print*, "in setxyz rlong0,rlat0,schmidt ", rlong0, rlat0, schmidt 
   end if
   do iq = 1, ifull 
!      scale wts so sum over globe is 1.0
!      wts(iq)=wts(iq)/(6.*sumwts)  ! for old conf-cub defn
      wts(iq) = wts(iq)/sumwts 
!      also provide latitudes and longitudes (-pi to pi)
      if ( rlong0 == 0.0 .and. rlat0 == 90.0 ) then 
         xx = x(iq) 
         yy = y(iq) 
         zz = z(iq) 
      else 
         xx = rotpole(1,1)*x(iq) + rotpole(1,2)*y(iq) + rotpole(1,3)*z(iq) 
         yy = rotpole(2,1)*x(iq) + rotpole(2,2)*y(iq) + rotpole(2,3)*z(iq) 
         zz = rotpole(3,1)*x(iq) + rotpole(3,2)*y(iq) + rotpole(3,3)*z(iq) 
      endif
      rlat(iq) = asin(zz) 
      f(iq) = 2.0 *zz * omega
      if ( abs(yy) > polelim .or. abs(xx) > polelim ) then 
         rlong(iq) = atan2(yy,xx)             ! N.B. -pi to pi 
!                                             ! 0 to 2*pi  09-25-1997
         if ( rlong(iq) < 0.0 ) then
            rlong(iq) = rlong(iq) + 2.0*pi 
         end if
      else 
         rlong(iq) = 0.0                      ! a default value for NP/SP 
      endif
      if ( iq /= idjd ) then
         cycle  
      end if
      if ( diag /= 0 ) then
         print*, "iq,x,y,z,xx,yy,zz,rlong,rlat ", iq, x(iq), y(iq), z(iq),  &
                 xx, yy, zz, rlong(iq), rlat(iq) 
      end if
   end do! iq loop 
!!$   if (diag == 2) then 
!!$      cosa(:ifull) = 100.*wts(:ifull) 
!!$      call printp ("wts ", cosa) 
!!$      call printp ("lat ", rlat) 
!!$      call printp ("long", rlong) 
!!$   endif
   if ( diag /= 0 ) then
      print *, "At centre of the faces:" 
      do n = 0, npanels 
         iq = ind((il + 1)/2,(il + 1)/2,n) 
         write ( unit=*,                                                      &
                 fmt="(a,i3,i5,3f7.3,2f8.2,es12.5)") "n,iq,x,y,z,long,lat,f ",&
                     n, iq, x(iq), y(iq), z(iq), rlong(iq)*rtod,              &
                     rlat(iq)*rtod, f(iq) 
      end do
      print *, "On each panel map factors for (1,1),(1,2),(1,3),(2,2),(3,2)" 
      do n = 0, npanels 
         iq11 = ind(1,1,n) 
         iq12 = ind(1,2,n) 
         iq13 = ind(1,3,n) 
         iq22 = ind(2,2,n) 
         iq32 = ind(3,2,n) 
         write(unit=*, fmt="(i3,5f8.3)") &
           n, em(iq11), em(iq12), em(iq13), em(iq22), em(iq32) 
      end do
      print *, "On each panel demu  demv for (1,1),(1,2),(1,3)" 
      do n = 0, npanels 
         iq = ind(1,1,n) 
         demu11 = emu(i_wu2(iq)) - emu(iq) 
         demv11 = emv(i_sv2(iq)) - emv(iq) 
         iq = ind(2,1,n) 
         demu21 = emu(i_wu2(iq)) - emu(iq) 
         demv21 = emv(i_sv2(iq)) - emv(iq) 
         iq = ind(3,1,n) 
         demu31 = emu(i_wu2(iq)) - emu(iq) 
         demv31 = emv(i_sv2(iq)) - emv(iq) 
         write ( unit=*, fmt="(i3,6f8.3)" )                                 &
              n, demu11, demv11, demu21, demv21, demu31, demv31 
      end do
   end if

!  set up Coriolis
!  ok: [2,il-1;2,jl]  u point
   dmdx = ( em(i_e) - em ) / ds 
! ok: [2,il;2,jl-1]  v point
   dmdy = ( em(i_n) - em ) / ds 
   fu = ( f + f(i_e) ) * 0.5_rx
   fv = ( f + f(i_n) ) * 0.5_rx 
!  may need more thought  $$$
   dmdx = abs ( emu(i_wu2)-emu(1:ifull) ) + abs ( emv(i_sv2)-emv(1:ifull) ) 
   dmdyu = ( em(i_n(i_e)) - em(i_s(i_e)) + em(i_n) - em(i_s) ) / (4.0*ds) 
   dmdxv = ( em(i_e(i_n)) - em(i_w(i_n)) + em(i_e) - em(i_w) ) / (4.0*ds) 

   if ( diag /= 0 ) then
      print *, &
         "On each panel abs(demu..demv) for (1,1),(1,2),(1,3),(2,2),(3,2)" 
      do n = 0, npanels 
         iq11 = ind(1,1,n) 
         iq12 = ind(1,2,n) 
         iq13 = ind(1,3,n) 
         iq22 = ind(2,2,n) 
         iq32 = ind(3,2,n) 
         write(unit=*, fmt="(i3,5f8.3)" )                                    &
            n, dmdx(iq11), dmdx(iq12), dmdx(iq13), dmdx(iq22), dmdx(iq32) 
      end do
   end if
 
!   For calculating zonal and meridional wind components, use the
!   following information, where theta is the angle between the
!   (ax,ay,az) vector [along the xg axis] and the zonal-component-vector:
!   veczon = k x r, i.e. (-y,x,0)/sqrt(x**2 + y**2)
!   vecmer = r x veczon, i.e. (-xz,-yz,x**2 + y**2)/sqrt(x**2 + y**2)
!   costh is (veczon . a) = (-y*ax + x*ay)/sqrt(x**2 + y**2)
!   sinth is (vecmer . a) = [-xz*ax - yz*ay + (x**2 + y**2)*az]/sqrt
!    using (r . a)=0, sinth collapses to az/sqrt(x**2 + y**2)
!   With more effort than expected, it can be verified that
!   costh**2 + sinth**2 = 0, by again using r.a =0.
 
!!!contains

end subroutine setxyz

subroutine vecpanel(ax, ay, az) 
   use newmpar_m, only : il  ! This is required when it's not an internal fn.
   real(kind=rx), intent(inout), dimension(:) :: ax, ay, az
   integer :: s0, e0, s1, e1, s2, e2, s3, e3, s4, e4, s5, e5

!-----------------------------------------------

!  Define indices here for the start and end of each face
   s0 = ind(1,1,0)
   e0 = ind(il,il,0)
   s1 = ind(1,1,1)
   e1 = ind(il,il,1)
   s2 = ind(1,1,2)
   e2 = ind(il,il,2)
   s3 = ind(1,1,3)
   e3 = ind(il,il,3)
   s4 = ind(1,1,4)
   e4 = ind(il,il,4)
   s5 = ind(1,1,5)
   e5 = ind(il,il,5)

   ax(s1:e1) = -az(s0:e0) 
   ay(s1:e1) =  ay(s0:e0) 
   az(s1:e1) =  ax(s0:e0) 
   ax(s2:e2) = -az(s0:e0) 
   ay(s2:e2) =  ax(s0:e0) 
   az(s2:e2) = -ay(s0:e0) 
   ax(s3:e3) = -ax(s0:e0) 
   ay(s3:e3) = -az(s0:e0) 
   az(s3:e3) = -ay(s0:e0) 
   ax(s4:e4) =  ay(s0:e0) 
   ay(s4:e4) = -az(s0:e0) 
   az(s4:e4) = -ax(s0:e0) 
   ax(s5:e5) =  ay(s0:e0) 
   ay(s5:e5) = -ax(s0:e0) 
   az(s5:e5) =  az(s0:e0) 

end subroutine vecpanel 

!!! end subroutine setxyz

subroutine norm(a, b, c) 
   real(kind=rx), intent(inout), dimension(:) :: a 
   real(kind=rx), intent(inout), dimension(:)  :: b 
   real(kind=rx), intent(inout), dimension(:)  :: c 
   real(kind=rx), dimension(size(a)) :: den 
!-----------------------------------------------
   den = sqrt(a**2 + b**2 + c**2) 
   a = a/den 
   b = b/den 
   c = c/den 
end subroutine norm

!!$subroutine printp(name, s6) 
!!$   character , intent(in) :: name*4 
!!$   real, intent(in)  :: s6(il,il,0:5) 
!!$!-----------------------------------------------
!!$!   L o c a l   P a r a m e t e r s
!!$!-----------------------------------------------
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      integer :: j, i 
!!$      real, dimension(0:il + 1,3*il) :: s1f, s2f 
!!$!-----------------------------------------------
!!$ 
!!$!     s1 is Grenwich-NP section i.e.  0-1-3
!!$!     s2 is Oz-SP section i.e.  2-4-5
!!$      call strip2 (s6, s6, s1f, s2f) 
!!$      print *, name, "  013" 
!!$      do j = 3*il, 1, -1 
!!$         print 9, j, (s1f(i,j),i=0,il + 1) 
!!$      end do 
!!$    9 format(i3,1x,21f6.3) 
!!$      print * 
!!$      print *, name, "  245" 
!!$      do j = 3*il, 1, -1 
!!$         print 9, j, (s2f(i,j),i=0,il + 1) 
!!$      end do 
!!$      return  
!!$      end subroutine printp 
!!$
!!$
!!$ 
!!$      subroutine strip2(s, s6, s1, s2) 
!!$!-----------------------------------------------
!!$!   M o d u l e s 
!!$!-----------------------------------------------
!!$      use newmpar_h_m 
!!$      use indices_h_m, ONLY: in, is, iw, ie 
!!$!...Translated by PSUITE Trans90                  4.3ZA 13:56:07   3/05/98  
!!$!...Switches: -xb -ya             
!!$      implicit none
!!$!-----------------------------------------------
!!$!   G l o b a l   P a r a m e t e r s
!!$!-----------------------------------------------
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      real , intent(in) :: s(ifull) 
!!$      real , intent(in) :: s6(il,il,0:5) 
!!$      real , intent(out) :: s1(0:il + 1,il,3) 
!!$      real , intent(out) :: s2(0:il + 1,il,3) 
!!$!-----------------------------------------------
!!$!   L o c a l   P a r a m e t e r s
!!$!-----------------------------------------------
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      integer , dimension(il,il,0:5) :: in6, is6, iw6, ie6 
!!$      integer :: j, i 
!!$!-----------------------------------------------
!!$      equivalence (in, in6), (is, is6), (iw, iw6), (ie, ie6) 
!!$!     N.B.  s & s6 are equivalenced via the call
!!$!     dimension s1f(0:il+1,3*il),s2f(0:il+1,3*il)
!!$!     s1 is Grenwich-NP section i.e.  0-1-3
!!$!     s2 is Oz-SP section i.e.  2-4-5
!!$!     for gnewst, these are extended on the sides only (i=0 & il+1)
!!$      s1(1:il,:il,1) = s6(:il,:il,0) 
!!$      s1(1:il,:il,2) = s6(:il,:il,1) 
!!$      s1(1:il,:il,3) = transpose(s6(:il,il:1:(-1),3)) 
!!$      s2(1:il,:il,1) = transpose(s6(:il,il:1:(-1),2)) 
!!$      s2(1:il,:il,2) = s6(:il,:il,4) 
!!$      s2(1:il,:il,3) = s6(:il,:il,5) 
!!$!                                                ! i loop
!!$      s1(0,:il,1) = s(iw6(1,:il,0)) 
!!$!      print *,"j,iw6(1,j,0),s1(0,j,1) ",j,iw6(1,j,0),s1(0,j,1)
!!$      s1(0,:il,2) = s(iw6(1,:il,1)) 
!!$      s1(0,:il,3) = s(in6(:il,il,3)) 
!!$      s2(0,:il,1) = s(in6(:il,il,2)) 
!!$      s2(0,:il,2) = s(iw6(1,:il,4)) 
!!$      s2(0,:il,3) = s(iw6(1,:il,5)) 
!!$      s1(il+1,:il,1) = s(ie6(il,:il,0)) 
!!$      s1(il+1,:il,2) = s(ie6(il,:il,1)) 
!!$      s1(il+1,:il,3) = s(is6(:il,1,3)) 
!!$      s2(il+1,:il,1) = s(is6(:il,1,2)) 
!!$      s2(il+1,:il,2) = s(ie6(il,:il,4)) 
!!$      s2(il+1,:il,3) = s(ie6(il,:il,5)) 
!!$!                                                ! j loop
!!$      return  
!!$      end subroutine strip2 
!!$

 
subroutine cross3 ( c1, c2, c3, a1, a2, a3, b1, b2, b3 ) 
   real(kind=rx) , intent(out), dimension(:) :: c1, c2, c3
   real(kind=rx) , intent(in), dimension(:)  :: a1, a2, a3
   real(kind=rx) , intent(in), dimension(:)  :: b1, b2, b3

   c1 = a2*b3 - b2*a3
   c2 = a3*b1 - b3*a1
   c3 = a1*b2 - b1*a2

end subroutine cross3


!!$real function crossmod (iq, i4a, j4a, i4b, j4b) 
!!$!-----------------------------------------------
!!$!   M o d u l e s 
!!$!-----------------------------------------------
!!$      use newmpar_h_m 
!!$      use bigxy4_h_m, ONLY: xx4, yy4 
!!$      use xyzinfo_h_m, ONLY: x, y, z 
!!$!...Translated by PSUITE Trans90                  4.3ZA 13:56:07   3/05/98  
!!$!...Switches: -xb -ya             
!!$      implicit none
!!$!-----------------------------------------------
!!$!   G l o b a l   P a r a m e t e r s
!!$!-----------------------------------------------
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      integer , intent(in) :: iq 
!!$      integer , intent(in) :: i4a 
!!$      integer , intent(in) :: j4a 
!!$      integer , intent(in) :: i4b 
!!$      integer , intent(in) :: j4b 
!!$!-----------------------------------------------
!!$!   L o c a l   P a r a m e t e r s
!!$!-----------------------------------------------
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      real , dimension(il,il,0:5) :: x6, y6, z6, x06, y06, z06 
!!$      real :: y4a, z4a, x4a, y4b, z4b, x4b, den, vecax, vecay, vecaz, vecbx, &
!!$         vecby, vecbz, crossx, crossy, crossz 
!!$!-----------------------------------------------
!!$      y4a = xx4(i4a,j4a) 
!!$      z4a = yy4(i4a,j4a) 
!!$      x4a = 1.0 
!!$      y4b = xx4(i4b,j4b) 
!!$      z4b = yy4(i4b,j4b) 
!!$      x4b = 1.0 
!!$!     if(iq.eq.1.or.iq.eq.il*il)then
!!$!      print *,"iq,x,y,z ",iq,x(iq),y(iq),z(iq)
!!$!      print *,"i4a,j4a,y4a,y4a ",i4a,j4a,y4a,y4a
!!$!      print *,"i4b,j4b,z4b,z4b ",i4b,j4b,z4b,z4b
!!$!     endif
!!$      call norm (x4a, y4a, z4a, den)             ! converts xx4, yy4 to coords on sphere 
!!$      call norm (x4b, y4b, z4b, den)             ! converts xx4, yy4 to coords on sphere 
!!$!     if(iq.eq.1.or.iq.eq.il*il)then
!!$!      print *,"after norm"
!!$!      print *,"x4a,y4a,z4a ",x4a,y4a,z4a
!!$!      print *,"x4b,y4b,z4b ",x4b,y4b,z4b
!!$!     endif
!!$      vecax = x4a - x(iq) 
!!$      vecay = y4a - y(iq) 
!!$      vecaz = z4a - z(iq) 
!!$      vecbx = x4b - x(iq) 
!!$      vecby = y4b - y(iq) 
!!$      vecbz = z4b - z(iq) 
!!$      crossx = vecay*vecbz - vecby*vecaz 
!!$      crossy = vecaz*vecbx - vecbz*vecax 
!!$      crossz = vecax*vecby - vecbx*vecay 
!!$      crossmod = sqrt(crossx**2 + crossy**2 + crossz**2) 
!!$!     if(iq.eq.6.or.iq.eq.31)then
!!$!       print *,"iq,iqa,iqb,crossmod ",iq,iqa,iqb,crossmod
!!$!       print *,"veca ",vecax,vecay,vecaz
!!$!       print *,"vecb ",vecbx,vecby,vecbz
!!$!     endif
!!$      return  
!!$      end function crossmod 

!--------------------------------------------
   subroutine setaxu(ifull)
      use indices_m
!     Set up vectors for the A grid 
      integer, intent(in) :: ifull
!     Pointer trick to make these contiguous
      real(kind=rx), allocatable, dimension(:), target  :: wcuva, wcuvd
      real(kind=rx), dimension(:), pointer :: wcua, wcva, wcud, wcvd
!-----------------------------------------------

!     These are just working arrays
      allocate ( wcuva(-ifull+1:ifull), wcuvd(-ifull+1:ifull) )
      wcua => wcuva(-ifull+1:ifull)
      wcva => wcuva
      wcud => wcuvd(-ifull+1:ifull)
      wcvd => wcuvd
      allocate ( axu(ifull), bxv(ifull), ayu(ifull), byv(ifull), &
                 azu(ifull), bzv(ifull) )

!  ax, bx etc are overdimensioned via pointer trick so have to use 1:ifull
!  explicitly here

!  convert ax,bx to staggered positions
      wcua(1:ifull) = 0.5_rx*(ax(i_eu2)+ax(1:ifull))
      wcud(1:ifull) = ax(i_eu2) - ax(1:ifull)
      wcva(1:ifull) = 0.5_rx*(bx(i_nv2)+bx(1:ifull))
      wcvd(1:ifull) = bx(i_nv2) - bx(1:ifull)

!  store ax,bx at staggered positions
      axu = wcua(1:ifull) - (wcud(i_eu2)-wcud(i_wu2))/16.0
      bxv = wcva(1:ifull) - (wcvd(i_nv2)-wcvd(i_sv2))/16.0

!  convert ay,by to staggered positions
      wcua(1:ifull) = 0.5_rx*(ay(i_eu2)+ay(1:ifull))
      wcud(1:ifull) = ay(i_eu2) - ay(1:ifull)
      wcva(1:ifull) = 0.5_rx*(by(i_nv2)+by(1:ifull))
      wcvd(1:ifull) = by(i_nv2) - by(1:ifull)

!  store ay,by at staggered positions
      ayu = wcua(1:ifull) - (wcud(i_eu2)-wcud(i_wu2))/16.0
      byv = wcva(1:ifull) - (wcvd(i_nv2)-wcvd(i_sv2))/16.0

!  convert az,bz to staggered positions
      wcua(1:ifull) = 0.5_rx*(az(i_eu2)+az(1:ifull))
      wcud(1:ifull) = az(i_eu2) - az(1:ifull)
      wcva(1:ifull) = 0.5_rx*(bz(i_nv2)+bz(1:ifull))
      wcvd(1:ifull) = bz(i_nv2) - bz(1:ifull)

!  store az,bz at staggered positions
      azu = wcua(1:ifull) - (wcud(i_eu2)-wcud(i_wu2))/16.0
      bzv = wcva(1:ifull) - (wcvd(i_nv2)-wcvd(i_sv2))/16.0

      deallocate ( wcuva, wcuvd )
   end subroutine setaxu

end module setxyz_m
