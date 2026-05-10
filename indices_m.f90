! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module indices_m

   use newmpar_m, only : il, ifull

   implicit none

   private
   public i_nn, i_ss, i_ww, i_ee
   public i_ne, i_se, i_en, i_wn
   public i_wu2, i_sv2
      
   public lwws, lwss, lees
   public less, lwwn, lwnn, leen
   public lenn, lsww, lssw
   public lsee, lsse, lnww
   public lnnw, lnee, lnne
   
   public i_n, i_s, i_e, i_w
   
   integer, parameter, dimension(0:5), private :: npann_g = (/ 1, 103, 3, 105, 5, 101 /)
   integer, parameter, dimension(0:5), private :: npans_g = (/ 104, 0, 100, 2, 102, 4 /)
   integer, parameter, dimension(0:5), private :: npane_g = (/ 102, 2, 104, 4, 100, 0 /)
   integer, parameter, dimension(0:5), private :: npanw_g = (/ 5, 105, 1, 101, 3, 103 /)
   
   integer, public, dimension(:), allocatable :: i_wu, i_sv, i_eu2, i_nv2
   integer, public, dimension(:), allocatable :: i_eu, i_nv
   
   !integer, public, dimension(:), allocatable :: i_nn, i_ss, i_ww, i_ee,     &
   !                                              i_ne, i_se, i_en, i_wn,     &
   !                                              i_wu, i_sv, i_wu2, i_sv2,   &
   !                                              i_eu2, i_nv2, i_ev2, i_nu2, &
   !                                              i_eu, i_nv 

   !integer, public, dimension(:), allocatable :: lwws, lws, lwss, les, lees,         &
   !                                              less, lwwn, lwnn, leen,             &
   !                                              lenn, lsww, lsw, lssw,              &
   !                                              lsee, lsse, lnww, lnw,              &
   !                                              lnnw, lnee, lnne
   
!#ifdef share_ifullg
!   integer, pointer, contiguous, dimension(:), public :: i_n, i_s, i_w, i_e
!   integer :: in_win, is_win, ie_win, iw_win
!#else
!   integer, public, dimension(:), allocatable :: i_n, i_s, i_w, i_e
!#endif
   
   contains
    
   function i_n(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (j==il) then
         if (npann_g(n)<100) then
            iqq=i+npann_g(n)*il*il
         else
            iqq=1+(il-i)*il+(npann_g(n)-100)*il*il
         end if
      else
         iqq=iq+il    
      end if
   end function i_n 
   
   function i_s(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (j==1) then
         if (npans_g(n)<100) then
            iqq=i+(il-1)*il+npans_g(n)*il*il
         else
            iqq=il+(il-i)*il+(npans_g(n)-100)*il*il
         end if
     else
         iqq=iq-il
      end if
   end function i_s

   function i_e(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (i==il) then
         if (npane_g(n)<100) then
            iqq=1+(j-1)*il+npane_g(n)*il*il
         else
            iqq=il+1-j+(npane_g(n)-100)*il*il
         end if
      else
         iqq=iq+1
      end if
   end function i_e

   function i_w(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (i==1) then
         if (npanw_g(n)<100) then
            iqq=il+(j-1)*il+npanw_g(n)*il*il
         else
            iqq=il+1-j+(il-1)*il+(npanw_g(n)-100)*il*il
         end if
      else
         iqq=iq-1
      end if
   end function i_w

   function i_nn(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npann_g(n)>=100.and.j==il) then
         iqq=i_e(i_n(iq))
      else
         iqq=i_n(i_n(iq))
      end if
   end function i_nn

   function i_ss(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npans_g(n)>=100.and.j==1) then
         iqq=i_w(i_s(iq))
      else
         iqq=i_s(i_s(iq))
      end if
   end function i_ss

   function i_ee(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npane_g(n)>=100.and.i==il) then
         iqq=i_n(i_e(iq))
      else
         iqq=i_e(i_e(iq))
      end if
   end function i_ee

   function i_ww(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npanw_g(n)>=100.and.i==1) then
         iqq=i_s(i_w(iq))
      else
         iqq=i_w(i_w(iq))
      end if
   end function i_ww

   function i_ne(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npane_g(n)>=100.and.i==il) then
         iqq=i_w(i_e(iq))
      else
         iqq=i_n(i_e(iq))
      end if
   end function i_ne

   function i_se(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npane_g(n)>=100.and.i==il) then
         iqq=i_e(i_e(iq))
      else
         iqq=i_s(i_e(iq))
      end if
   end function i_se

   function i_en(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npann_g(n)>=100.and.j==il) then
         iqq=i_s(i_n(iq))
      else
         iqq=i_e(i_n(iq))
      end if
   end function i_en

   function i_wn(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npann_g(n)>=100.and.j==il) then
         iqq=i_n(i_n(iq))
      else
         iqq=i_w(i_n(iq))
      end if
   end function i_wn
   
   function i_nw(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npanw_g(n)>=100.and.i==1) then
         iqq=i_w(i_w(iq))
      else
         iqq=i_n(i_w(iq))
      end if
   end function i_nw

   function i_sw(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npanw_g(n)>=100.and.i==1) then
         iqq=i_e(i_w(iq))
      else
         iqq=i_s(i_w(iq))
      end if
   end function i_sw

   function i_es(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npans_g(n)>=100.and.j==1) then
         iqq=i_s(i_s(iq))
      else
         iqq=i_e(i_s(iq))
      end if
   end function i_es

   function i_ws(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npans_g(n)>=100.and.j==1) then
         iqq=i_n(i_s(iq))
      else
         iqq=i_w(i_s(iq))
      end if
   end function i_ws

   function i_wu2(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npanw_g(n)>=100.and.i==1) then
         iqq=i_w(iq)+ifull
      else
         iqq=i_w(iq)
      end if
   end function i_wu2

   function i_sv2(iq) result(iqq)
      integer, intent(in) :: iq
      integer iqq
      integer n,i,j
      n=(iq-1)/(il*il)
      j=(iq-1-n*il*il)/il+1
      i=iq-(j-1)*il-n*il*il
      if (npans_g(n)>=100.and.j==1) then
         iqq=i_s(iq)-ifull
      else
         iqq=i_s(iq)
      end if
   end function i_sv2

   function leen(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npann_g(n)>=100) then
         iqq=i_ss(i_n(il+(il-1)*il+n*il*il))
      else
         iqq=i_ee(i_n(il+(il-1)*il+n*il*il))
      end if
   end function leen

   function lenn(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npann_g(n)>=100) then
         iqq=i_se(i_n(il+(il-1)*il+n*il*il))
      else
         iqq=i_en(i_n(il+(il-1)*il+n*il*il))
      end if
   end function lenn

   function lwnn(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npann_g(n)>=100) then
         iqq=i_ne(i_n(1+(il-1)*il+n*il*il))
      else
         iqq=i_wn(i_n(1+(il-1)*il+n*il*il))
      end if
   end function lwnn

   function lsee(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npane_g(n)>=100) then
         iqq=i_wn(i_e(il+n*il*il))
      else
         iqq=i_se(i_e(il+n*il*il))
      end if
   end function lsee

   function lnee(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npane_g(n)>=100) then
         iqq=i_wn(i_e(il+(il-1)*il+n*il*il))
      else
         iqq=i_ne(i_e(il+(il-1)*il+n*il*il))
      end if
   end function lnee

   function lnne(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npane_g(n)>=100) then
         iqq=i_ww(i_e(il+(il-1)*il+n*il*il))
      else
         iqq=i_nn(i_e(il+(il-1)*il+n*il*il))
      end if
   end function lnne

   function lsww(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npanw_g(n)>=100) then
         iqq=i_es(i_w(1+n*il*il))
      else
         iqq=i_sw(i_w(1+n*il*il))
      end if
   end function lsww

   function lssw(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npanw_g(n)>=100) then
         iqq=i_ee(i_w(1+n*il*il))
      else
         iqq=i_ss(i_w(1+n*il*il))
      end if
   end function lssw

   function lnww(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npanw_g(n)>=100) then
         iqq=i_ws(i_w(1+(il-1)*il+n*il*il))
      else
         iqq=i_nw(i_w(1+(il-1)*il+n*il*il))
      end if
   end function lnww

   function lwws(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npans_g(n)>=100) then
         iqq=i_nn(i_s(1+n*il*il))
      else
         iqq=i_ww(i_s(1+n*il*il))
      end if
   end function lwws

   function lwss(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npans_g(n)>=100) then
         iqq=i_nw(i_s(1+n*il*il))
      else
         iqq=i_ws(i_s(1+n*il*il))
      end if
   end function lwss

   function less(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npans_g(n)>=100) then
         iqq=i_sw(i_s(il+n*il*il))
      else
         iqq=i_es(i_s(il+n*il*il))
      end if
   end function less

   function lwwn(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npann_g(n)>=100) then
         iqq=i_nn(i_n(1+(il-1)*il+n*il*il))
      else
         iqq=i_ww(i_n(1+(il-1)*il+n*il*il))
      end if
   end function lwwn

   function lsse(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npane_g(n)>=100) then
         iqq=i_ee(i_e(il+n*il*il))
      else
         iqq=i_ss(i_e(il+n*il*il))
      end if
   end function lsse

   function lnnw(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npanw_g(n)>=100) then
         iqq=i_ww(i_w(1+(il-1)*il+n*il*il))
      else
         iqq=i_nn(i_w(1+(il-1)*il+n*il*il))
      end if
   end function lnnw

   function lees(n) result(iqq)
      integer, intent(in) :: n
      integer iqq
      if (npans_g(n)>=100) then
         iqq=i_ss(i_s(il+n*il*il))
      else
         iqq=i_ee(i_s(il+n*il*il))
      end if
   end function lees
   
end module indices_m
