! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module latltoij_m

contains
   subroutine latltoij ( rlongin, rlatin, xout, yout, nf,                   &
                         rlong0, rlat0, schmidt, schm13)
!     Given a pair of latitudes and longitudes (in degrees),
!     returns i and j values on the conformal-cubic grid as
!     xout ranging between .5 and   il +.5, and
!     yout ranging between .5 and 6*il +.5

      use newmpar_m
      use xyzinfo_m, only : xx4, yy4
      use precis_m, only : rx
      implicit none

      real, intent(in) :: rlongin, rlatin
      real, intent(out) :: xout, yout
      integer, intent(out) :: nf
      real(kind=rx), intent(in) :: rlong0, rlat0, schmidt, schm13
      integer, parameter :: diag=0

      real(kind=rx), save, dimension(3,3) :: rotpolei

!     following used by npanels=13
      integer, parameter, dimension(-1:1,-1:4) :: npanetab =                 &
        reshape ( (/ 9,12,13, 7,2,1, 6,4,3, 6,5,3, 8,10,0, 9,11,13 /),       &
                  (/ 3, 6 /) )
      real(kind=rx), parameter, dimension(-1:1,-1:4) :: acon =                        &
        reshape ( (/ 0,1,1, -1,0,0, -1,0,0, 0,0,-1, 1,1,0, 1,1,0 /),         &
                  (/ 3, 6 /) )
      real(kind=rx), parameter, dimension(-1:1,-1:4) :: bcon =                        &
        reshape ( (/ 1,0,0, 0,-1,-1, 0,-1,-1, -1,-1,0, 0,0,1, 0,0,1 /),      &
                  (/ 3, 6 /) )
      real(kind=rx), parameter, dimension(-1:1,-1:4) :: xadd =                        &
        reshape ( (/ -0.5, 0.5, -0.5,   -0.5, 0.5, 0.5,    -0.5,-0.5,-0.5,   &
                     -1.5,-1.5, 1.5,    1.5, 0.5, 3.5,      1.5, 0.5, 4.5 /),&
                  (/ 3, 6 /) )
      real(kind=rx), parameter, dimension(-1:1,-1:4) :: yadd =                        &
        reshape ( (/ 1.5, 1.5, 1.5,   0.5, 0.5, 1.5,    1.5, 0.5, 1.5,       &
                    -0.5, 0.5, 2.5,  -2.5,-2.5,-0.5,   -3.5,-3.5,-0.5 /),    &
                  (/ 3, 6 /) )

!        acon:   0    -1    -1       0     1     1     9  7  6   6  8  9
!                1     0     0       0     1     1    12  2  4   5 10 11
!                1     0     0      -1     0     0    13  1  3   3  0 13
!        bcon:   1     0     0      -1     0     0
!                0    -1    -1      -1     0     0
!                0    -1    -1       0     1     1
!        xadd: -0.5  -0.5  -0.5    -1.5   1.5   1.5      correct ones
!               0.5   0.5  -0.5    -1.5   0.5   0.5
!              -0.5   0.5  -0.5     1.5   3.5   4.5
!        yadd:  1.5   0.5   1.5    -0.5  -2.5  -3.5
!               1.5   0.5   0.5     0.5  -2.5  -3.5
!               1.5   1.5   1.5     2.5  -0.5  -0.5

      integer, save :: num = 0
      real(kind=rx), save :: alf, schmidtp, alfp, schmidtm, alfm
      integer, parameter :: nmaploop = 3,  numtst= 1063
      real(kind=rx), parameter :: pi = 3.1415926536 
      real(kind=rx) :: coslong, sinlong, coslat, sinlat
      real(kind=rx) :: xa, ya, za, x, y, z, x1, y1, z1, denxyz, xx, yy, zz
      real(kind=rx) :: xgrid, ygrid
      real(kind=rx) :: ri, rj, dxx, dyx, dxy, dyy, den, xsign, ysign, zsign
      real(kind=rx) :: xstr, ystr, xgr, ygr 
      integer :: i, j, ibox, jbox, loop, ig, jg, is, js
      real(kind=rx), save :: rmax
      real(kind=rx), dimension(3) :: vec


      if ( num == 0 ) then
         rmax = nearest(real(iquad),-1.0)
         alf = (1.0-schmidt**2)/(1.0+schmidt**2)
         if ( npanels == 13 ) then
            schmidtp = schm13/schmidt        ! for z3d above the "equator"
            alfp = (1.0-schmidtp**2)/(1.0+schmidtp**2)
            schmidtm = 1.0/(schm13*schmidt)   ! for z3d below the "equator"
            alfm = (1.0-schmidtm**2)/(1.0+schmidtm**2)
         end if
         if ( diag /= 0 ) then
            print *,"latltoij; rlong0,rlat0,schmidt,schm13,alf:",       &
                          rlong0,rlat0,schmidt,schm13,alf
            print *,"xx4(iquad,1) ",xx4(iquad,1)
         end if
         coslong = cos(rlong0*pi/180.)
         sinlong = sin(rlong0*pi/180.)
         coslat = cos(rlat0*pi/180.)
         sinlat = sin(rlat0*pi/180.)
         rotpolei(1,1) = coslong*sinlat
         rotpolei(2,1) = -sinlong
         rotpolei(3,1) = coslong*coslat
         rotpolei(1,2) = sinlong*sinlat
         rotpolei(2,2) = coslong
         rotpolei(3,2) = sinlong*coslat
         rotpolei(1,3) = -coslat
         rotpolei(2,3) = 0.
         rotpolei(3,3) = sinlat
      endif
      num = num+1
!     numtst = num
      if ( num == numtst .and. diag /= 0 ) then
         print *,"a rlongin,rlatin ",rlongin,rlatin
      end if
      xa = cos(rlongin*pi/180.)*cos(rlatin*pi/180.)
      ya = sin(rlongin*pi/180.)*cos(rlatin*pi/180.)
      za = sin(rlatin*pi/180.)
      if ( num == numtst .and. diag /= 0 ) then
         print *,"b xa,ya,za ",xa,ya,za
      end if
      x = rotpolei(1,1)*xa+rotpolei(1,2)*ya+rotpolei(1,3)*za
      y = rotpolei(2,1)*xa+rotpolei(2,2)*ya+rotpolei(2,3)*za
      z = rotpolei(3,1)*xa+rotpolei(3,2)*ya+rotpolei(3,3)*za
      if ( num == numtst .and. diag /= 0 ) then
         print *,"c x,y,z ",x,y,z
      end if

      if ( npanels== 5 ) then

!       If necessary, transform physical (x, y, z) to equivalent coordinates
!       on regular gnomonic panels

        if ( schmidt /= 1.0 ) then
           x1 = x
           y1 = y
           z1 = z
           x = x*(1.-alf)/(schmidt*(1.-alf*z))
           y = y*(1.-alf)/(schmidt*(1.-alf*z))
           z = (z-alf)/(1.-alf*z)
           if ( z1 > 0.82 .and. z1 .lt. 0.821 .and. diag /= 0 ) then
              print *,"latltoij: rlongin, rlatin ",rlongin, rlatin
              print *,"latltoij: xa,ya,za ",xa,ya,za
              print *,"latltoij: x1,y1,z1 ",x1,y1,z1
              print *,"latltoij: x,y,z ",x,y,z
           endif
        endif         ! (schmidt.ne.1.)

        vec = (/ x, z, y /)
        nf = sum(maxloc(abs(vec)))
        ! NF is 1, 2, 3 according to whether x, z, or y has the largest 
        ! magnitude. sum is just to reduce the dimension(1) result of maxloc
        ! to a scalar.
        ! If the value is negative add 3 to get the opposite face
        if ( vec(nf) < 0 ) then
           nf = nf + 3
        end if
        nf = nf - 1 ! Bring back to 0-5 range.
        denxyz = max( abs(x),abs(y),abs(z) )
        xx = x/denxyz
        yy = y/denxyz
        zz = z/denxyz
        select case ( nf ) 
        case ( 0 )
           xgrid = yy
           ygrid = zz
        case ( 3 )
           xgrid = -zz
           ygrid = -yy
        case ( 1 )
           xgrid =  yy
           ygrid = -xx
        case ( 4 )
           xgrid =  xx
           ygrid = -yy
        case ( 2 )
           xgrid = -zz
           ygrid = -xx
        case ( 5 )
           xgrid =  xx
           ygrid =  zz
        end select

!      Convert to grid point numbering
!      The xytoij routine follows

!      Use 4* resolution grid il --> 4*il
       xgrid = min(max(-.99999_rx,xgrid),.99999_rx)
       ygrid = min(max(-.99999_rx,ygrid),.99999_rx)
!      First guess for ri, rj and nearest ig,jg
       ri = 1.0 + (1.0+xgrid)*2*il
       rj = 1.0 + (1.0+ygrid)*2*il
       do loop = 1,nmaploop
          ig = nint(ri)
          jg = nint(rj)
          is = sign(1.0_rx,ri-ig)
          js = sign(1.0_rx,rj-jg)
!         predict new value for ri, rj
          dxx = xx4(ig+is,jg)-xx4(ig,jg)
          dyx = xx4(ig,jg+js)-xx4(ig,jg)
          dxy = yy4(ig+is,jg)-yy4(ig,jg)
          dyy = yy4(ig,jg+js)-yy4(ig,jg)
          den = dxx*dyy-dyx*dxy
          ri = ig+is*((xgrid-xx4(ig,jg))*dyy-(ygrid-yy4(ig,jg))*dyx)/den
          rj = jg+js*((ygrid-yy4(ig,jg))*dxx-(xgrid-xx4(ig,jg))*dxy)/den
          
          ri = min(ri,1.0_rx+1.99999_rx*real(2*il,rx))
          ri = max(ri,1.0_rx+0.00001_rx*real(2*il,rx))
          rj = min(rj,1.0_rx+1.99999_rx*real(2*il,rx))
          rj = max(rj,1.0_rx+0.00001_rx*real(2*il,rx))
       enddo! loop loop
       xout = 0.25*(ri+3.0) - 0.5  ! -.5 for stag; back to normal ri, rj defn
       yout = 0.25*(rj+3.0) - 0.5  ! -.5 for stag
!      Expect xout, yout (at this point) to range between .5 and il+.5

    elseif ( npanels == 13 ) then

!   First convert to equivalent of schmidt=.5 grid
       if ( z > alf ) then
          zsign = 1.0
          xstr = x*schmidtp*(1.0+alfp)/(1.0+alfp*z)
          ystr = y*schmidtp*(1.0+alfp)/(1.0+alfp*z)
       else
          zsign = -1.0
          xstr = x*schmidtm*(1.+alfm)/(1.+alfm*z)
          ystr = y*schmidtm*(1.+alfm)/(1.+alfm*z)
       endif!  (z.gt.alf)
!      could avoid above "if", by first doing 1/schmidt, then schmidt13
!      using abs(z). Extra operations may then increase errors slightly?
!      now remember departure quadrants
        xsign = sign(1.0_rx,xstr)
        ysign = sign(1.0_rx,ystr)
        if ( num == numtst .and. diag /= 0 ) then
           print *,"d xstr,ystr,xsign,ysign,zsign ",                      &
                    xstr,ystr,xsign,ysign,zsign
        end if

!       Use 4* resolution grid
!       N.B. for toij13, have xx4 and yy4 between 0 and .5 (after schmidtx )
        xgr = abs(xstr)
        ygr = abs(ystr)
!       first guess for ri, rj (1 to 6*il+1) and nearest i,j
        ri = 1.0+xgr*(iquad-1)/xx4(iquad,1)    ! divide by schm13 "equator" radius
        rj = 1.0+ygr*(iquad-1)/xx4(iquad,1)
        if ( num == numtst .and. diag /= 0 ) then
           print *,"e xgr,ygr,ri,rj",xgr,ygr,ri,rj
        end if
        do loop = 1,nmaploop
!        ri = max(1. , min(real(iquad)-.0001,ri) )  !  not needed
!        rj = max(1. , min(real(iquad)-.0001,rj) )  !  not needed
           ri = max(1.0_rx, min(rmax,ri))
           rj = max(1.0_rx, min(rmax,rj))
           i = nint(ri)
           j = nint(rj)
           is = sign(1.0_rx,ri-i)
           js = sign(1.0_rx,rj-j)
!          predict new value for ri, rj
           dxx = xx4(i+is,j)-xx4(i,j)
           dyx = xx4(i,j+js)-xx4(i,j)
           dxy = yy4(i+is,j)-yy4(i,j)
           dyy = yy4(i,j+js)-yy4(i,j)
           den = dxx*dyy-dyx*dxy
           ri = i+is*((xgr-xx4(i,j))*dyy-(ygr-yy4(i,j))*dyx)/den
           rj = j+js*((ygr-yy4(i,j))*dxx-(xgr-xx4(i,j))*dxy)/den
        enddo! loop loop
!       write ri,rj on a BIG grid (-1.5*il to 1.5*il, -1.5*il to 4.5*il)
!       where the Y variable is wrapping around the globe
!       and the "north pole" is now at (0.,0.)
        ri = 0.25*(ri-1.0)*xsign
        rj = 0.25*( (rj-1.0)*ysign*zsign + (1.0-zsign)*real(6*il) )
        if ( num == numtst .and. diag /= 0 ) then
           print *,"e2 bigy  ri,rj ",ri,rj
        end if
!       allocate to a box (-1:1, -1:4)
        ibox = max(-1,min(nint(ri/il),1))   ! allows for -1.5 or 1.5
        jbox = max(-1,min(nint(rj/il),4))   ! allows for -1.5 or 4.5
!       convert  xg, yg ( .5 to il+.5)
        if ( num == numtst .and. diag /= 0 ) then
           print *,"f ri,rj,ibox,jbox",ri,rj,ibox,jbox
        end if
        if( num == numtst .and. diag /= 0 ) then
           print *,"f1 nf,xadd,yadd,acon,bcon ", npanetab(ibox,jbox),       &
         xadd(ibox,jbox),yadd(ibox,jbox),acon(ibox,jbox),bcon(ibox,jbox)
        endif
        nf = npanetab(ibox,jbox)
        xout = 0.5 +xadd(ibox,jbox)*real(il) + acon(ibox,jbox)*ri -    &
               bcon(ibox,jbox)*rj
        yout = 0.5 +yadd(ibox,jbox)*real(il) + bcon(ibox,jbox)*ri +    &
               acon(ibox,jbox)*rj
        if ( num == numtst .and. diag /= 0 ) then
           print *,"f1 xadd,yadd,acon,bcon ",                               &
          xadd(ibox,jbox),yadd(ibox,jbox),acon(ibox,jbox),bcon(ibox,jbox)
           print *,"g nf,xout,yout ",nf,xout,yout
        end if

      endif  !  (npanels == 5) elseif(npanels == 13)

   end subroutine latltoij
end module latltoij_m
