module moistfuncs

   implicit none
!  Moisture related calculations, based on CSIRO model ESTABL.f

!  MKS table
!  TABLE OF ES fROM -150 C TO +70 C IN ONE-DEGREE INCREMENTS.
!  Temperature of last ES in each row is shown as a comment.
   real, private, dimension(0:220), parameter :: table = (/            &
     1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                & !-146C
     6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                             & !-141C
     36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                          & !-136C
     0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,  & !-131C
     0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,  & !-126C
     0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,   & !-121C
     0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,       & !-116C
     0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,         & !-111C
     0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,            & !-106C
     0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,             & !-101C
     .001403, .001719, .002101, .002561, .003117, .003784,             & !-95C
     .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658, & !-87C
     .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577, & !-78C
     .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032,   & !-69C
     .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080,    & !-60C
     1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476,    & !-51C
     3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,            & !-43C
     10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,    & !-34C
     27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85, & !-24C 
     77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67,      & !-16C
     171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78,  & !-8C
     353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78,   & !0C
     656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,   & !8C
     1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,   & !16C
     1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,   & !24C
     3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,   & !32C
     5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,   & !40C
     7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,         & !47C
     11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,    & !54C
     15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,    & !61C
     21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,    & !68C
     29845.0, 31169.0 /)                                                 !70C

   real, private, dimension(0:220), parameter :: tablei = (/             &
     & 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                & !-146C
     & 6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                             & !-141C
     & 36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                          & !-136C
     & 0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,  & !-131C
     & 0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,  & !-126C
     & 0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,   & !-121C
     & 0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,       & !-116C
     & 0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,         & !-111C
     & 0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,            & !-106C
     & 0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,             & !-101C
     & .001403, .001719, .002101, .002561, .003117, .003784,             & !-95C
     & .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658, & !-87C
     & .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577, & !-78C
     & .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032,   & !-69C
     & .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080,    & !-60C
     & 1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476,    & !-51C
     & 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,            & !-43C
     & 10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,    & !-34C
     & 27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85, & !-24C 
     & 77.09, 85.02, 93.70, 103.06, 113.40, 124.68, 136.98, 150.39,      & !-16C
     & 164.99, 180.88, 198.16, 216.94, 237.34, 259.47, 283.49, 309.51,   & !-8C
     & 337.71, 368.23, 401.25, 436.96, 475.54, 517.21, 562.19, 610.70,   & !0C
     & 656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,   & !8C
     & 1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,   & !16C
     & 1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,   & !24C
     & 3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,   & !32C
     & 5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,   & !40C
     & 7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,         & !47C
     & 11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,    & !54C
     & 15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,    & !61C
     & 21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,    & !68C
     & 29845.0, 31169.0/)                                                  !70C 

   real, private, dimension(-40:2), parameter :: esdiff = (/             &
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,       &
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,      &
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,      &
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,      &
     & 0.08, 0., 0. /)

   real, parameter, private :: epsil=0.622
   private :: tdiff, tdiffx
   public :: establ, estabi, qsat, qsati, relhum, esdiffx

contains

   elemental real function tdiff(t)
! T is temp in Kelvin, which should lie between 123.16 and 343.16;
! TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220
      real, intent(in) :: t
      tdiff = min(max( t-123.16, 0.), 219.)
   end function tdiff

   elemental real function tdiffx(t)
      real, intent(in) :: t
      tdiffx = min(max( t-273.1, -40.), 1.)
   end function tdiffx
   
   elemental real function establ(t)
      ! Saturation vapor pressure in Pa
      real, intent(in) :: t ! (K)
      establ = (1.-(tdiff(t)-int(tdiff(t))))*table(int(tdiff(t))) +    &
                  (tdiff(t)-int(tdiff(t)))*table(int(tdiff(t))+1)
   end function establ

   elemental real function estabi(t)
      real, intent(in) :: t ! (K)
      estabi = (1.-(tdiff(t)-aint(tdiff(t))))*tablei(int(tdiff(t))) +    &
                  (tdiff(t)-aint(tdiff(t)))*tablei(int(tdiff(t))+1)
   end function estabi

!  These functions all use mixing ratio rather than specific humidity
!  to match the CSIRO model
               
   elemental real function qsat(p,t)
      real, intent(in) :: p, t  ! p in Pa, t in K
!      qsat = epsil*establ(t)/p !Consistent with V4-5 to V4-7
      qsat = epsil*establ(t)/max(p-establ(t),0.1) ! MJT bug fix from JJK
   end function qsat

   elemental real function qsati(p,t)
      real, intent(in) :: p, t  ! p in Pa, t in K
      qsati = epsil*estabi(t)/max(p-estabi(t),0.1)
   end function qsati
   
   elemental real function relhum(p,q,t)
      real, intent(in) :: p, q, t ! p in Pa, q in kg/kg, t in K
      ! Result in range 0:1.
      relhum = q / qsat(p,t)
   end function relhum

   elemental real function tvirt(t,q)
      real, intent(in) :: t, q ! q in kg/kg, t in K
      tvirt = t * (epsil+q)/(epsil*(1.+q))
   end function tvirt

   elemental real function es_gg(t)
      ! Groff-Gratch formule, result in Pa
      real, intent(in) :: t ! K
      es_gg = 101324.6 * 10 ** ( -7.90298 * ( 373.16/t - 1 ) + &
                     5.02808 * log10 (373.16/t) - &
             1.3816e-7 * ( 10**(11.344*(1 - t/373.16)) -1 ) +  &
             8.1328e-3 * (10**(-3.49149*(373.16/t -1 )) -1 ))
   end function es_gg

   elemental real function esdiffx(t)
      real, intent(in) :: t
      esdiffx = (1.-(tdiffx(t)-aint(tdiffx(t))))*esdiff(int(tdiffx(t))) +    &
                  (tdiffx(t)-aint(tdiffx(t)))*esdiff(int(tdiffx(t))+1)
   end function esdiffx

   elemental function tdew(press, qmix, temp, aerr) RESULT(fn_val)
 
      ! Modified version of zeroin.f90 
      ! Now find zero of qsat(p,t) - q to calculate dew point temperature

      ! Code converted using TO_F90 by Alan Miller
      ! Date: 2003-07-14  Time: 12:32:54
      
      !-----------------------------------------------------------------------
      
      !         FINDING A ZERO OF THE FUNCTION F(X) IN THE INTERVAL (AX,BX)
      
      !                       ------------------------
      
      !  INPUT...
      
      !  F      FUNCTION SUBPROGRAM WHICH EVALUATES F(X) FOR ANY X IN THE
      !         CLOSED INTERVAL (AX,BX).  IT IS ASSUMED THAT F IS CONTINUOUS,
      !         AND THAT F(AX) AND F(BX) HAVE DIFFERENT SIGNS.
      !  AX     LEFT ENDPOINT OF THE INTERVAL
      !  BX     RIGHT ENDPOINT OF THE INTERVAL
      !  AERR   THE ABSOLUTE ERROR TOLERANCE TO BE SATISFIED
      !  RERR   THE RELATIVE ERROR TOLERANCE TO BE SATISFIED
      
      !  OUTPUT...
      
      !         ABCISSA APPROXIMATING A ZERO OF F IN THE INTERVAL (AX,BX)
      
      !-----------------------------------------------------------------------
      !  ZEROIN IS A SLIGHTLY MODIFIED TRANSLATION OF THE ALGOL PROCEDURE
      !  ZERO GIVEN BY RICHARD BRENT IN ALGORITHMS FOR MINIMIZATION WITHOUT
      !  DERIVATIVES, PRENTICE-HALL, INC. (1973).
      !-----------------------------------------------------------------------
      

      REAL, INTENT(IN)  :: press, qmix
      REAL, INTENT(IN), optional :: temp
      REAL, INTENT(IN), optional  :: aerr
      REAL :: rerr
      REAL              :: fn_val

      REAL :: a, b, c, d, e, eps, fa, fb, fc, tol, xm, p, q, r, s, atol, rtol

      !  COMPUTE EPS, THE RELATIVE MACHINE PRECISION

      eps = EPSILON(0.0)

      ! INITIALIZATION

      if ( present(temp) ) then
         ! Use this as starting guess
         a = temp + 5. ! Allow slight supersaturation
         b = max(150.,temp-50.)
      else
         ! Use full qsat range -150 to 70
         a = 123.
         b = 343.
      end if

      fa = qsat(press,a) - qmix
      fb = qsat(press,b) - qmix
      if ( present(aerr) ) then
         atol = 0.5 * aerr
      else
         atol = 1e-4
      end if
      rerr = 0.0
      rtol = MAX(0.5*rerr, 2.0*eps)

      ! BEGIN STEP

10    c = a
      fc = fa
      d = b - a
      e = d
20    IF (ABS(fc) < ABS(fb)) THEN
         a = b
         b = c
         c = a
         fa = fb
         fb = fc
         fc = fa
      END IF

      ! CONVERGENCE TEST

      tol = rtol * MAX(ABS(b),ABS(c)) + atol
      xm = 0.5 * (c-b)
      IF (ABS(xm) > tol) THEN
         IF (fb /= 0.0) THEN

            ! IS BISECTION NECESSARY

            IF (ABS(e) >= tol) THEN
               IF (ABS(fa) > ABS(fb)) THEN

                  ! IS QUADRATIC INTERPOLATION POSSIBLE

                  IF (a == c) THEN

                     ! LINEAR INTERPOLATION

                     s = fb / fc
                     p = (c-b) * s
                     q = 1.0 - s
                  ELSE

                     ! INVERSE QUADRATIC INTERPOLATION

                     q = fa / fc
                     r = fb / fc
                     s = fb / fa
                     p = s * ((c-b)*q*(q-r)-(b-a)*(r-1.0))
                     q = (q-1.0) * (r-1.0) * (s-1.0)
                  END IF

                  ! ADJUST SIGNS

                  IF (p > 0.0) q = -q
                  p = ABS(p)

                  ! IS INTERPOLATION ACCEPTABLE

                  IF (2.0*p < (3.0*xm*q-ABS(tol*q))) THEN
                     IF (p < ABS(0.5*e*q)) THEN
                        e = d
                        d = p / q
                        GO TO 30
                     END IF
                  END IF
               END IF
            END IF

            ! BISECTION

            d = xm
            e = d

            ! COMPLETE STEP

30          a = b
            fa = fb
            IF (ABS(d) > tol) b = b + d
            IF (ABS(d) <= tol) b = b + SIGN(tol,xm)
            fb = qsat(press,b) - qmix
            IF (fb*(fc/ABS(fc)) > 0.0) GO TO 10
            GO TO 20
         END IF
      END IF

      ! DONE

      fn_val = b
      RETURN
   END FUNCTION tdew

end module moistfuncs
