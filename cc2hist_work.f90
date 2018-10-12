! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2018 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module work

   use mpidata_m
#ifdef usenc_mod
   use netcdf
#else
   use netcdf_m
#endif
   use ncutils_m, only : check_ncerr
   use gldata
   use precis_m, only : rx
   implicit none

   integer, parameter :: un_in = 10 ! Input unit number
   logical, private, save :: first_in = .true.
   
   integer, private, save :: resprocmode
   logical, private, save :: resprocformat

!  Passed between infile and main, previously in common infil
   integer :: ik, jk, kk, ok
   integer :: nlev, onlev  ! Number of levels in output
   integer :: ndt, ktau, nqg, nsd
   integer :: ksoil, kice
   integer :: ntrac, ilt ! Extra flag controlling tracers in file. Still needed?

!  Flag for presence of extra surface flux fields (epot_ave etc).
   logical :: extra_surfflux
!  Flag for presence of extra snow fields
   logical :: extra_snow
!  Flag for presence of extra 6, 12 and 18 hourly rain fields
   logical :: extra_rain, extra_rain_3hr
!  Record number for the netcdf input file
   integer, save :: nrec
   integer, save :: maxrec ! Length of record dimension.
   integer, save :: ncid  ! ID of the input file

   real, allocatable, dimension(:), public, save :: hlon, hlat
   integer, public, save :: nxhis, nyhis
   real, dimension(:,:), allocatable, private :: costh, sinth

   public :: initialise, fix_winds, final_init, check_cc2histfile
   public :: paraopen, paraclose, cc2hist_work_close
   private :: need3dfld, fill_cc
   interface fix_winds
      module procedure fix_winds2, fix_winds3
   end interface

   interface vread
      module procedure vread2, vread3
   end interface

   interface match
      module procedure matcha, matchc
   end interface
   
   interface ccmpi_scatter
      module procedure ccmpi_scatter_1d_r4_host, ccmpi_scatter_1d_r4_proc
      module procedure ccmpi_scatter_2d_r4_host, ccmpi_scatter_2d_r4_proc
   end interface   

   type input_var
      character(len=20) :: vname   ! Name used in input
      ! Initialise these so they don't end up with trailing nulls.
      character(len=100) :: long_name=""
      character(len=20) :: units=""
      real :: add_offset
      real :: scale_factor
      integer :: ndims ! Number of space dimensions
      logical :: fixed ! Does it have time dimension
      integer :: vid   ! netcdf vid
      logical :: daily ! Is it only valid once per day?
      logical :: vector ! Is it a vector component?
      logical :: xcmpnt ! Is it x-component of a vector?
      logical :: water  ! Is it a multi-level ocean variable
      integer :: othercmpnt ! Index of matching x or y component
   end type input_var

contains

   subroutine alloc_indata ( il, jl, kl, ol, ksoil, kice )

!     Allocates arrays for the variables that have to be stored rather than
!     being immediately processed.
      use history, only : needfld
      use s2p_m

      integer, intent(in) :: il, jl, kl, ol, ksoil, kice
      
      allocate ( psl(pil,pjl*pnpan*lproc),   zs(pil,pjl*pnpan*lproc) )
      allocate ( soilt(pil,pjl*pnpan*lproc), u(pil,pjl*pnpan*lproc,kl),  v(pil,pjl*pnpan*lproc,kl) )
      allocate ( t(pil,pjl*pnpan*lproc,kl) )
      allocate ( q(pil,pjl*pnpan*lproc,kl),  ql(pil,pjl*pnpan*lproc,kl), qf(pil,pjl*pnpan*lproc,kl) )
      allocate ( qs(pil,pjl*pnpan*lproc,kl), qg(pil,pjl*pnpan*lproc,kl), omega(pil,pjl*pnpan*lproc,kl) )
      allocate ( tgg(pil,pjl*pnpan*lproc,ksoil) )
      if ( needfld("zg") ) then
         if ( use_plevs .or. use_meters ) then
            allocate ( zstd(pil,pjl*pnpan*lproc,nplevs) )
         else
!           Sigma level height calculation requires all levels.
            allocate ( zstd(pil,pjl*pnpan*lproc,kl) )
         end if
      end if
      if ( needfld("press") .or. needfld("theta") .or. needfld("rh") .or. needfld("w") .or. needfld('td') ) then
        allocate( tmp3d(pil,pjl*pnpan*lproc,kl) )
      end if
      if ( use_meters ) then
         allocate( hstd(pil,pjl*pnpan*lproc,kl) )
      end if
      if ( ol > 0 ) then
         allocate( uo_tmp(pil,pjl*pnpan*lproc,ol), vo_tmp(pil,pjl*pnpan*lproc,ol) )
         allocate( thetao_tmp(pil,pjl*pnpan*lproc,ol), so_tmp(pil,pjl*pnpan*lproc,ol) )
         allocate( ocn_tmp(pil,pjl*pnpan*lproc,ol) )
      end if

   end subroutine alloc_indata

   subroutine getdate ( kdate, ktime, ieof) 
      use logging_m

!     Get record data from header
      integer, intent(out) :: kdate, ktime, ieof
      integer :: ik, jk, kk, ierr, vid
      integer :: iposa, iposb, yyyy, mm, dd, hh, mt
      character(len=80) :: datestring

      if ( nrec > maxrec ) then
         ieof = 1
         return
      end if
      call START_LOG(getdate_begin)      
      ierr = nf90_inq_varid (ncid, "time", vid )
      ierr = nf90_get_att(ncid, vid, "units", datestring )
      iposa = index(trim(datestring),'since')
      iposa = iposa + 5 ! skip 'since'
      iposb = index(trim(datestring(iposa:)),'-')
      iposb = iposa + iposb - 2 ! remove '-'
      read(datestring(iposa:iposb),FMT=*,iostat=ierr) yyyy
      if ( ierr/=0 ) then
        write(6,*) "ERROR reading time units.  Expecting year but found ",datestring(iposa:iposb)
        stop
      end if
      iposa = iposb + 2 ! skip '-'
      iposb = index(trim(datestring(iposa:)),'-')
      iposb = iposa + iposb - 2 ! remove '-'
      read(datestring(iposa:iposb),FMT=*,iostat=ierr) mm
      if ( ierr/=0 ) then
        write(6,*) "ERROR reading time units.  Expecting month but found ",datestring(iposa:iposb)
        stop
      end if
      iposa = iposb + 2 ! skip '-'
      iposb = index(trim(datestring(iposa:)),' ')
      iposb = iposa + iposb - 2 ! remove ' '
      read(datestring(iposa:iposb),FMT=*,iostat=ierr) dd
      if ( ierr/=0 ) then
        write(6,*) "ERROR reading time units.  Expecting day but found ",datestring(iposa:iposb)
        stop
      end if
      iposa = iposb + 2 ! skip ' '
      iposb = index(trim(datestring(iposa:)),':')
      iposb = iposa + iposb - 2 ! remove ':'
      read(datestring(iposa:iposb),FMT=*,iostat=ierr) hh
      if ( ierr/=0 ) then
        write(6,*) "ERROR reading time units.  Expecting hour but found ",datestring(iposa:iposb)
        stop
      end if
      iposa = iposb + 2 ! skip ':'
      iposb = index(trim(datestring(iposa:)),':')
      iposb = iposa + iposb - 2 ! remove ':'
      read(datestring(iposa:iposb),FMT=*,iostat=ierr) mt
      if ( ierr/=0 ) then
        write(6,*) "ERROR reading time units.  Expecting minutes but found ",datestring(iposa:iposb)
        stop
      end if
      kdate = yyyy*10000 + mm*100 + dd
      ktime = hh*100 + mt
      !! Get vid and then values for kdate, ktime, ktau
      !ierr = nf90_inq_varid (ncid, "kdate", vid )
      !call check_ncerr(ierr, "Error getting kdate id")
      !ierr = nf90_get_var ( ncid, vid, kdate, start=(/ nrec /) )
      !call check_ncerr(ierr, "Error getting kdate")
      !ierr = nf90_inq_varid (ncid, "ktime", vid )
      !call check_ncerr(ierr, "Error getting ktime id")
      !ierr = nf90_get_var ( ncid, vid, ktime, start=(/ nrec /) )
      !call check_ncerr(ierr, "Error getting ktime")
      ! Get ktau from time. Really should be renamed
      !ierr = nf90_inq_varid (ncid, "time", vid )
      !call check_ncerr(ierr, "Error getting time id")
      ierr = nf90_get_var ( ncid, vid, ktau, start=(/ nrec /) )
      call check_ncerr(ierr, "Error getting time")
      ieof = 0
            
      if ( myid == 0 ) then
         write(*,"(a,3i8)",advance="no") " kdate, ktime, ktau ", &
                                           kdate,  ktime, ktau
      end if
      call END_LOG(getdate_end)

   end subroutine getdate
   
   subroutine getstep(ktc)
   
      integer, intent(out) :: ktc
      integer :: vid, ierr, ktau0, ktau1
   
      ierr = nf90_inq_varid (ncid, "time", vid )
      call check_ncerr(ierr, "Error getting time id")
      
      ierr = nf90_get_var( ncid, vid, ktau0, start=(/ 1 /) )
      call check_ncerr(ierr, "Error getting time")
      ierr = nf90_get_var( ncid, vid, ktau1, start=(/ 2 /) )
      if ( ierr == nf90_noerr ) then
         ktc = max( ktau1 - ktau0, 1 )
      else
         print *,"WARN: Cannot locate second time-step in input file" 
         ktc = 1 
      end if    
   
   end subroutine getstep

   subroutine convert_date( mins, kdate, ktime, kt, basetime )
   
      integer, intent(out) :: mins
      integer, intent(in) :: kdate, ktime, kt
      integer :: jyear, jmonth, jday, jhour, jmin
      integer :: newdate, newtime, mtimer
      integer, dimension(12), parameter :: cdays=(/0,31,59,90,120,151,181,212,243,273,304,334/)
      character(len=*), intent(in) :: basetime
      
      jyear = kdate/10000
      newdate = kdate - jyear*10000
      jmonth = newdate/100
      newdate = newdate - jmonth*100
      jday = newdate
      jhour = ktime/100
      newtime = ktime - jhour*100
      jmin = newtime
         
      if ( basetime(1:13) == "minutes since" ) then
         mtimer = kt
      else if ( basetime(1:13) == "seconds since" ) then
         mtimer = kt/60 
      else
         write(6,*) "ERROR: Unknown time unit in input file"
         stop
      end if
      mins = (cdays(jmonth)+jday)*1440 + jhour*24 + jmin + mtimer

   end subroutine convert_date
   
   subroutine infile ( varlist, nvars, skip, mins )
      ! For netcdf input
      use history, only : savehist, needfld, cordex_compliant
      use physparams, only : grav, rdry, cp, pi
      use s2p_m
      use height_m
      use newmpar_m, only : ol
      use sitop_m
      use logging_m
      use parm_m, only : rlong0, rlat0
      use zenith_m
      type(input_var), dimension(:) :: varlist
      integer, intent(in) :: nvars, mins
      logical, intent(in)  :: skip
      integer :: k, ivar, ierr
      integer :: rad_day
      real, dimension(pil,pjl*pnpan*lproc) :: uten, udir, dtmp, ctmp, uastmp, vastmp
      real, dimension(pil,pjl*pnpan*lproc) :: mrso, mrfso
      real, dimension(pil,pjl*pnpan*lproc) :: sndw, egave, dpsdt
      real, dimension(pil,pjl*pnpan*lproc) :: tauxtmp, tauytmp
      real, dimension(pil,pjl*pnpan*lproc) :: rgn, rgd, sgn, sgd
      real, dimension(pil,pjl*pnpan*lproc) :: wind_norm
      real, dimension(1) :: rlong_a, rlat_a, cos_zen, frac
      real :: fjd, bpyear, r1, dlt, alp, slag, dhr
      character(len=10) :: name
      logical :: validvar
      real, parameter :: spval   = 999.
      real, parameter :: tfreeze = 271.38

      ! With a netcdf file there's no need to read in order to skip.
      if ( skip ) then
         nrec = nrec + 1
         return
      end if

      call START_LOG(infile_begin)
      if ( first_in ) then
         call alloc_indata ( ik, jk, kk, ol, ksoil, kice )
      end if
      
      mrso = 0.  ! total soil moisture
      mrfso = 0. ! total soil ice
      
      do ivar=1,nvars
         ! Just write the input with no further processing

         ! Files with fixed type only written on first pass
         if ( varlist(ivar)%fixed ) then
            if ( first_in ) then
               select case ( varlist(ivar)%vname )
               case ( "zs", "orog" )
                  ! This could also be done with the output_scale
                  call vread( "zht", zs )
                  zs = zs / grav
                  if ( needfld(varlist(ivar)%vname) ) then
                     call savehist ( varlist(ivar)%vname, zs )
                  end if   
               case ( "soilt" )
                  call vread( "soilt", soilt )
                  if ( needfld("soilt") ) then
                     call savehist("soilt", soilt)
                  end if   
                  if ( needfld("land_mask") ) then
                     where ( soilt > 0.5 )
                        dtmp = 1.
                     elsewhere
                        dtmp = 0.
                     end where
                     call savehist("land_mask", dtmp)
                  else if ( needfld("sftlf") ) then
                     where ( soilt > 0.5 )
                        dtmp = 100.
                     elsewhere
                        dtmp = 0.
                     end where
                     call savehist("sftlf", dtmp)
                  endif
               case ( "map" )
                  if ( needfld("map") .or. needfld("grid") ) then 
                     call vread( "map", dtmp )
                     if ( needfld("map") ) then
                        call savehist("map", dtmp )
                     end if   
                     if ( needfld("grid") ) then
                        dtmp = 90.*112./(real(pil_g)*dtmp)
                        call savehist("grid", dtmp)
                     end if   
                  end if
               case ( "ocndepth" )
                  call vread( "ocndepth", dtmp )
                  if ( needfld("ocndepth") ) then
                     call savehist( "ocndepth", dtmp )
                  end if   
                  if ( use_depth ) call ditop_setup( gosig, dlevs(1:onplevs), dtmp )
               case default
                  call readsave2 ( varlist(ivar)%vname )
               end select
            end if
            cycle ! Otherwise will match the following if test as well
         end if

        ! Only read these on a full day. ktau is really time in minutes
         if ( varlist(ivar)%daily .and. modulo(ktau,1440) /= 0 ) cycle

         if ( varlist(ivar)%ndims == 2 ) then
            select case ( varlist(ivar)%vname )
            case ( "clivi" )
               call readsave2 (varlist(ivar)%vname, input_name="iwp_ave")
            case ( "clh" )
               if ( needfld("clh") ) then 
                  call vread( "clh", dtmp )
                  if ( cordex_compliant ) then
                     where ( dtmp /= nf90_fill_float )
                        dtmp = dtmp*100.
                     end where  
                  end if
                  call savehist ( "clh", dtmp )
               end if   
            case ( "cll" )
               if ( needfld("cll") ) then 
                  call vread( "cll", dtmp )
                  if ( cordex_compliant ) then
                     where ( dtmp /= nf90_fill_float ) 
                        dtmp = dtmp*100.
                     end where
                  end if
                  call savehist ( "cll", dtmp )
               end if   
            case ( "clm" )
               if ( needfld("clm") ) then 
                  call vread( "clm", dtmp )
                  if ( cordex_compliant ) then
                     where ( dtmp /= nf90_fill_float ) 
                        dtmp = dtmp*100.
                     end where   
                  end if
                  call savehist ( "clm", dtmp )
               end if   
            case ( "clt" )
               if ( needfld("clt") ) then 
                  call vread( "cld", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp*100.
                  end where   
                  call savehist ( "clt", dtmp )
               end if   
            case ( "clwvi" )
               call readsave2 (varlist(ivar)%vname, input_name="lwp_ave")
            case ( "dpsdt" )
                if ( needfld("dpsdt") .or. needfld("w") ) then
                   call vread( "dpsdt", dpsdt )
                   if ( needfld("dpsdt") ) then
                      call savehist ( "dpsdt", dpsdt )
                   end if   
                end if   
            case ( "evspsblpot" )
               if ( needfld("evspsblpot") ) then 
                  call vread( "epot_ave", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/2.501e6 ! Latent heat of vaporisation (J kg^-1)
                  end where   
                  call savehist ( "evspsblpot", dtmp )
               end if   
            case ( "hfls" )
               if ( needfld("hfls") .or. needfld("evspsbl") ) then 
                  call vread( "eg_ave", egave )
                  if ( needfld("hfls") ) then
                     call savehist( "hfls", egave )
                  end if   
               end if   
            case ( "hfss" )
               call readsave2 (varlist(ivar)%vname, input_name="fg_ave")
            case ( "hurs" )
               call readsave2 (varlist(ivar)%vname, input_name="rhscrn")
            case ( "huss" )
               if ( needfld("huss") ) then 
                  call vread( "qgscrn", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/(dtmp+1.)
                  end where  
                  call savehist ( "huss", dtmp )
               end if   
            case ( "mrro" )
               if ( needfld("mrro") ) then 
                  call vread( "runoff", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/86400.
                  end where   
                  call savehist ( "mrro", dtmp )
               end if   
            case ( "mrros" )
               if ( needfld("mrros") ) then 
                  call vread( "mrros", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/86400.
                  end where   
                  call savehist ( "mrros", dtmp )
               end if   
            case ( "ocheight" )
               if ( needfld("ocheight") ) then 
                  call vread( "ocheight", dtmp )
                  where ( soilt > 0.5 )
                     dtmp = nf90_fill_float
                  end where   
                  call savehist ( "ocheight", dtmp )
               end if   
            case ( "pr" )
               if ( needfld("pr") ) then 
                  call vread( "rnd", dtmp )
                  dtmp = dtmp/86400.
                  call savehist ( "pr", dtmp )
               end if   
            case ( "prc" )
               if ( needfld("prc") ) then 
                  call vread( "rnc", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/86400.
                  end where   
                  call savehist ( "prc", dtmp )
               end if   
            case ( "prmax" )
               if ( needfld("prmax") ) then 
                  call vread( "maxrnd", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/86400.
                  end where   
                  call savehist ( "prmax", dtmp )
               end if   
            case ( "prsn" )
               if ( needfld("prsn") ) then 
                  call vread( "sno", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/86400.
                  end where   
                  call savehist ( "prsn", dtmp )
               end if   
            case ( "psl" )
               if ( needfld("psl") ) then 
                  call vread( "pmsl", dtmp )
                  if ( cordex_compliant ) then
                     where ( dtmp /= nf90_fill_float ) 
                        dtmp = dtmp*100.
                     end where   
                  end if
                  call savehist ( "psl", dtmp )
               end if   
            case ( "psl_ave" )
               if ( needfld("psl_ave") ) then 
                  call vread( "pmsl_ave", dtmp )
                  if ( cordex_compliant ) then
                     where ( dtmp /= nf90_fill_float ) 
                        dtmp = dtmp*100.
                     end where   
                  end if
                  call savehist ( "psl_ave", dtmp )
               end if   
            case ( "psf" )
               call vread( "psf", psl )
               where ( psl /= nf90_fill_float )
                  psl = 1.0e3 * exp(psl) ! hPa
               end where
               if ( needfld("ps") ) then
                  if ( cordex_compliant ) then
                     where ( dtmp /= nf90_fill_float ) 
                        dtmp = 100.*psl     ! Pa
                     end where   
                  else    
                     dtmp = psl          ! hPa
                  end if
                  call savehist ( "ps", dtmp )
               end if   
               ! This relies on surface pressure coming before the 3D variables
               if ( use_plevs ) call sitop_setup(sig, plevs(1:nplevs), psl, maxlev, minlev)
            case ( "rgdn_ave", "rlds" )
               if ( needfld("rgdn") .or. needfld("rlds") .or. needfld("rlus") ) then 
                  call vread( "rgdn_ave", rgd )
                  if ( needfld(varlist(ivar)%vname) ) then
                     call savehist( varlist(ivar)%vname, rgd )
                  end if   
               end if   
            case ( "rgn_ave" )
               if ( needfld("rgn_ave") .or. needfld("rlus") ) then 
                  call vread( "rgn_ave", rgn )
                  if ( needfld("rgn_ave") ) then
                     call savehist( "rgn_ave", rgn )
                  end if   
               end if   
            case ( "rlut" )
               call readsave2( varlist(ivar)%vname, input_name="rtu_ave" )
            case ( "rsdt" )
               call readsave2 (varlist(ivar)%vname, input_name="sint_ave")
            case ( "rsut" )
               call readsave2 (varlist(ivar)%vname, input_name="sot_ave") 
            case ( "sic" )
               if ( needfld("sic") ) then 
                  call vread2( "fracice", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp*100.
                  end where   
                  call savehist( "sic", dtmp )
               end if   
            case ( "sgdn_ave", "rsds" )
               if ( needfld("sgdn_ave") .or. needfld("rsds") .or. needfld("rsus") ) then 
                  call vread( "sgdn_ave", sgd )
                  if ( needfld(varlist(ivar)%vname) ) then
                     call savehist( varlist(ivar)%vname, sgd )
                  end if   
               end if   
            case ( "sgn_ave" )
               if ( needfld("sgn_ave") .or. needfld("rsus") ) then 
                  call vread( "sgn_ave", sgn )
                  if ( needfld("sgn_ave") ) then
                     call savehist( "sgn_ave", sgn )
                  end if   
               end if   
            case ( "snd" )
               if ( needfld("snd") .or. needfld("snc") .or. needfld("snw") ) then
                  call vread( "snd", sndw )
                  if ( needfld("snd") ) then
                     if ( cordex_compliant ) then
                        where ( dtmp /= nf90_fill_float )  
                           dtmp = sndw*10. ! change from equiv water to equiv snow
                        end where   
                     else
                        dtmp = sndw
                     end if
                     call savehist ( "snd", dtmp )
                  end if   
               end if   
            case ( "snm" )
               if ( needfld("snm") ) then 
                  call vread( "snm", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp/86400.
                  end where   
                  call savehist ( "snm", dtmp )
               end if   
            case ( "sund" )
               if ( needfld("sund") ) then 
                  call vread( "sunhours", dtmp )
                  where ( dtmp /= nf90_fill_float )
                     dtmp = dtmp*3600.
                  end where   
                  call savehist ( "sund", dtmp )
               end if   
            case ( "tas" )
               call readsave2 (varlist(ivar)%vname, input_name="tscrn")
            case ( "tasmax" )
               call readsave2 (varlist(ivar)%vname, input_name="tmaxscr")
            case ( "tasmin" )
               call readsave2 (varlist(ivar)%vname, input_name="tminscr")
            case ( "tauu", "taux" ) 
               call vread( "taux", tauxtmp ) 
            case ( "tauv", "tauy" )
               call vread( "tauy", tauytmp )  
            case ( "ts", "tsu" )
               if ( needfld("ts") .or. needfld("tsu") .or. needfld("tsea") ) then 
                  call vread( "tsu", dtmp )
                  ! Some (all?) initial conditions have tsu negative over ocean
                  where ( dtmp /= nf90_fill_float )
                     dtmp = abs(dtmp)
                  end where
                  if ( needfld(varlist(ivar)%vname) ) then
                     call savehist(varlist(ivar)%vname, dtmp)
                  end if   
                  if ( needfld("tsea") ) then
                     where ( dtmp /= nf90_fill_float )
                        dtmp = max(tfreeze, dtmp)
                     end where
                     ! Use soilt as a land-sea mask (integer but read as float)
                     where ( soilt > 0.5 )
                        ctmp = spval
                     elsewhere
                        ! Use the maxval to ignore ice points.
                        ctmp = dtmp
                     end where
                     call fill_cc(ctmp, spval)
                     call savehist ( "tsea", ctmp )
                  end if   
               end if
            case ( "u10", "sfcWind" )
               call vread( "u10", uten )
               if ( needfld(varlist(ivar)%vname) ) then
                  call savehist( varlist(ivar)%vname, uten )
               end if   
            case ( "uas" )
                call vread( "uas", uastmp ) ! only for high-frequency output
            case ( "vas" )
                call vread( "vas", vastmp ) ! only for high-frequency output
            case ( "zmla" )
               call readsave2 (varlist(ivar)%vname, input_name="pblh")
            case default
               if ( varlist(ivar)%vector) then
                  if (varlist(ivar)%xcmpnt) then
                     ! Process both x and y now, skip over y when we come to 
                     ! it later. Use taux and tauy as temporary arrays
                     ! name of other component
                     name = varlist(varlist(ivar)%othercmpnt)%vname
                     if ( needfld(varlist(ivar)%vname) .or. needfld(name) ) then
                        call vread( varlist(ivar)%vname, ctmp)
                        call vread( name, dtmp)
                        call fix_winds(ctmp, dtmp)
                        if ( needfld(varlist(ivar)%vname) ) then
                           call savehist ( varlist(ivar)%vname, ctmp )
                        end if
                        if ( needfld(name) ) then
                           call savehist ( name, dtmp )
                        end if   
                     end if
                  endif
               else if ( varlist(ivar)%vname(1:3)=='tgg' .or.     &
                         varlist(ivar)%vname(1:7)=='rooftgg' .or. &
                         varlist(ivar)%vname(1:7)=='waletgg' .or. &
                         varlist(ivar)%vname(1:7)=='walwtgg' .or. &
                         varlist(ivar)%vname(1:7)=='roadtgg' ) then
                  if ( needfld(varlist(ivar)%vname) ) then
                     ! Fix soil, ocean or urban temperature offset
                     call vread(varlist(ivar)%vname,dtmp)
                     where ( dtmp<100. .and. dtmp/=nf90_fill_float )
                        dtmp = dtmp + 290. ! reference temperature
                     end where
                     call savehist(varlist(ivar)%vname,dtmp)
                  end if   
               else
                  if ( needfld(varlist(ivar)%vname) ) then           
                     call vread( varlist(ivar)%vname, ctmp ) 
                     validvar = .true.
                     do k = 1,size(zse)
                        write(name,'(a,i1.1,a)') 'wb', k,'_ave'
                        if ( varlist(ivar)%vname == name ) then
                           mrso = mrso + ctmp*zse(k)*1000.
                        end if
                        write(name,'(a,i1.1,a)') 'wbice', k,'_ave'
                        if ( varlist(ivar)%vname == name ) then
                           mrfso = mrfso + ctmp*zse(k)*1000. 
                           validvar = .false.
                        end if    
                     end do
                     if ( validvar ) then
                        call savehist( varlist(ivar)%vname, ctmp )
                     end if
                  end if   
               end if
            end select
         else
            select case ( varlist(ivar)%vname )
            case ( "ta", "temp" )
               ! temp should be the first of the 3D fields
               if ( use_meters ) then
                  minlev = 1
                  maxlev = kk
                  ! assume that 2D zs is previously loaded
                  ! MJT notes - reading mixr skips ahead in the input file
                  ! possibly reorder temp, mixr, u and v in CCAM
                  call vread( "temp", t)
                  call vread( "mixr", q)
                  q = max( q, 1.e-20 )
                  ! psl will not be used in height
                  call height( t, q, zs, psl, sig, hstd )
                  do k = 1,size(hstd,dim=3)
                     hstd(:,:,k) = mlevs(nplevs)*(hstd(:,:,k) - zs)/(mlevs(nplevs)-zs)
                  end do
                  call mitop_setup( sig, mlevs(1:nplevs), hstd, zs, t, q, maxlev, minlev )
                  if ( needfld(varlist(ivar)%vname) ) then
                     call vsavehist ( varlist(ivar)%vname, t )
                  end if
               else if ( need3dfld(varlist(ivar)%vname) ) then
                  call vread( "temp", t)
                  if ( needfld(varlist(ivar)%vname) ) then
                     call vsavehist ( varlist(ivar)%vname, t )
                  end if   
               end if
             case ( "hus", "mixr" )
               if ( need3dfld(varlist(ivar)%vname)) then
                  if ( .not. use_meters ) then
                     call vread( "mixr", q )
                     q = max( q, 1.e-20 )
                  end if
                  if ( needfld(varlist(ivar)%vname) ) then
                     call vsavehist ( varlist(ivar)%vname, q )
                  end if   
               end if
            case ( "qlg" )
               if ( need3dfld("qlg")) then
                  call vread( "qlg", ql )
                  ql = max( ql, 0. )
                  if ( needfld("qlg") ) then
                     call vsavehist ( "qlg", ql )
                  end if   
               end if
            case ( "qfg" )
               if ( need3dfld("qfg")) then
                  call vread( "qfg", qf )
                  qf = max( qf, 0. )
                  if ( needfld("qfg") ) then
                     call vsavehist ( "qfg", qf )
                  end if   
               end if
            case ( "qsng" )
               if ( need3dfld("qsng")) then
                  call vread( "qsng", qs )
                  qs = max( qs, 0. )
                  if ( needfld("qsng") ) then
                     call vsavehist ( "qsng", qs )
                  end if   
               end if
            case ( "qgrg" )
               if ( need3dfld("qgrg")) then
                  call vread( "qgrg", qg )
                  qg = max( qg, 0. )
                  if ( needfld("qgrg") ) then
                     call vsavehist ( "qgrg", qg )
                  end if   
               end if
            case ( "so" )
               if ( need3dfld("so") ) then
                  call vread( "so", so_tmp )
                  do k = 1,size(so_tmp,3)
                     where ( soilt > 0.5 )
                        so_tmp(:,:,k) = nf90_fill_float
                     end where
                  end do
                  if ( needfld("so") ) then
                     call osavehist( "so", so_tmp )
                  end if   
               end if
            case ( "thetao" )
               if ( need3dfld("thetao") ) then
                  call vread( "thetao", thetao_tmp )
                  do k = 1,size(thetao_tmp,3)
                     where ( soilt > 0.5 )
                        thetao_tmp(:,:,k) = nf90_fill_float
                     end where
                  end do
                  if ( needfld("thetao") ) then
                     call osavehist( "thetao", thetao_tmp )
                  end if   
               end if    
            ! Should to u, v as above with vector flag, but this will do for now
            case ( "u", "ua" )
               if ( need3dfld(varlist(ivar)%vname) ) then
                  call vread( "u", u )
               end if
            case ( "v", "va" )
               if ( need3dfld(varlist(ivar)%vname) ) then
                  call vread( "v", v )
               end if
            case ( "uo" ) 
               if ( need3dfld("uo") ) then
                  call vread( "uo", uo_tmp ) 
                  do k = 1,size(uo_tmp,3)
                     where ( soilt > 0.5 )
                        uo_tmp(:,:,k) = nf90_fill_float
                     end where
                  end do   
               end if
            case ( "vo" ) 
               if ( need3dfld("vo") ) then
                  call vread( "vo", vo_tmp ) 
                  do k = 1,size(vo_tmp,3)
                     where ( soilt > 0.5 )
                        vo_tmp(:,:,k) = nf90_fill_float
                     end where
                  end do   
               end if
            case ( "wo" )
               if ( need3dfld("wo") ) then
                  call vread( "wo", ocn_tmp )
                  do k = 1,size(ocn_tmp,3)
                     where ( soilt > 0.5 )
                        ocn_tmp(:,:,k) = nf90_fill_float
                     end where
                  end do
                  if ( needfld("wo") ) then
                     call osavehist( "wo", ocn_tmp )
                  end if   
               end if   
            case ( "omega" )
               if ( need3dfld("omega") ) then
                  call vread( "omega", omega )
                  if ( needfld("omega") ) then
                     call vsavehist( "omega", omega )
                  end if   
               end if
            case ( "tgg" )
               if ( needfld("tgg") ) then 
                  ! Only in cf_compliant mode
                  do k = 1,ksoil
                     write(name,'(a,i1)') 'tgg', k
                     call vread(name,tgg(:,:,k))
                  end do
                  call savehist("tgg", tgg)
               end if   
            case ( "wetfrac" )
               if ( needfld("wetfrac") ) then 
                  ! Only in cf_compliant mode
                  do k = 1,ksoil
                     write(name,'(a,i1)') 'wetfrac', k
                     call vread(name,tgg(:,:,k))
                  end do
                  call savehist("wetfrac", tgg)
               end if   
            case default
               if ( varlist(ivar)%water ) then
                  write(6,*) "ERROR: Not expecting ocean scalar ",trim(varlist(ivar)%vname)
                  stop
               else
                  call readsave3 (varlist(ivar)%vname)
               end if   
            end select
         end if

      end do

      if ( needfld("evspsbl") ) then
         dtmp = egave/2.501e6 ! Latent heat of vaporisation (J kg^-1)
         call savehist( "evspsbl", dtmp )
      end if

      if ( needfld("mrfso") ) then
         ierr = nf90_inq_varid (ncid, "wbice1_ave", ivar )
         if ( ierr/=nf90_noerr ) then
            mrfso = NF90_FILL_FLOAT  
         end if
         call savehist( "mrfso", mrfso )
      end if
      
      if ( needfld("mrso") ) then
         call savehist( "mrso", mrso )
      end if   

      if ( needfld("rlus") ) then
         dtmp = rgn - rgd
         call savehist( "rlus", dtmp )
      end if

      if ( needfld("rsus") ) then
         dtmp = sgn - sgd
         call savehist( "rsus", dtmp )
      end if
      
      if ( needfld("snc") ) then
         where ( sndw>0. )
            dtmp = 100.
         elsewhere
            dtmp = 0.  
         end where
         call savehist( "snc", dtmp )
      end if
      
      if ( needfld("snw") ) then
         call savehist( "snw", sndw )
      end if   
      
      if ( needfld("tauu") .or. needfld("tauv") .or. &
           needfld("taux") .or. needfld("tauy") .or. &
           (needfld("uas").and.kk>1) .or.            &
           (needfld("vas").and.kk>1) .or.            &
           (needfld("d10").and.kk>1) ) then
         call fix_winds(tauxtmp, tauytmp)
         if ( needfld("tauu") ) then
            call savehist( "tauu", tauxtmp )
         end if
         if ( needfld("tauv") ) then
            call savehist( "tauv", tauytmp )
         end if
         if ( needfld("taux") ) then
            call savehist( "taux", tauxtmp )
         end if
         if ( needfld("tauy") ) then
            call savehist( "tauy", tauytmp )
         end if   
         wind_norm(:,:) = sqrt(tauxtmp*tauxtmp+tauytmp*tauytmp)
         if ( needfld("uas") .and. kk>1 ) then
            where ( wind_norm > 0. )
               dtmp = tauxtmp*uten/wind_norm
            elsewhere
               dtmp = 0. 
            end where    
            call savehist ( "uas", dtmp )    
         end if    
         if ( needfld("vas") .and. kk>1 ) then
            where ( wind_norm > 0. )
               dtmp = tauytmp*uten/wind_norm
            elsewhere
               dtmp = 0. 
            end where    
            call savehist ( "vas", dtmp )    
         end if
         if ( needfld("d10") .and. kk>1 ) then
            udir = atan2(-u(:,:,1),-v(:,:,1))*180./3.1415927
            where ( udir < 0. )
               udir = udir + 360.
            end where
            call savehist( "d10", udir )        
         end if    
      end if

      ! high-frequency output 
      if ( (needfld("uas").and.kk<=1) .or. &
           (needfld("vas").and.kk<=1) .or. &
           (needfld("u10").and.kk<=1) .or. &
           (needfld("d10").and.kk<=1) ) then
         call fix_winds( uastmp, vastmp )
         if ( needfld("uas") ) then
            call savehist( "uas", uastmp )
         end if
         if ( needfld("vas") ) then
            call savehist( "vas", vastmp )
         end if   
         if ( needfld("u10") ) then
            uten = sqrt( uastmp*uastmp + vastmp*vastmp )
            call savehist( "u10", uten )
         end if
         if ( needfld("d10") ) then
            udir = atan2(-uastmp,-vastmp)*180./3.1415927
            where ( udir < 0. )
               udir = udir + 360.
            end where
            call savehist( "d10", udir )
         end if
      end if

      if ( kk>1 ) then
         
         if ( needfld("qbot") ) then
            call savehist( "qbot", q(:,:,1))
         end if   
         
         if ( needfld("tbot") ) then
            call savehist( "tbot", t(:,:,1))
         end if   
          
         if ( needfld("press") ) then
            do k = 1,kk
               tmp3d(:,:,k) = psl*sig(k)
            end do
            call vsavehist( "press", tmp3d )
         end if

         if ( needfld("prw") ) then
            dtmp = 0.0
            do k = 1,kk
               dtmp = dtmp + dsig(k)*q(:,:,k)
            end do
            dtmp = 100.*psl/grav * dtmp
            call savehist ( "prw", dtmp )
         end if
         
         if ( needfld("pwc") ) then
            dtmp = 0.0
            do k = 1,kk
               dtmp = dtmp + dsig(k)*q(:,:,k)
            end do
            dtmp = 100.*psl/grav * dtmp
            call savehist ( "pwc", dtmp )
         end if
         
         if ( needfld("td") ) then
            call calc_td ( t, q, ql, qf, psl, sig, tmp3d )
            call vsavehist( "td", tmp3d ) 
         end if
         
         if ( needfld("rh") ) then
            call calc_rh ( t, q, ql, qf, psl, sig, tmp3d )
            call vsavehist ( "rh", tmp3d )
         end if
         
         if ( needfld("theta") ) then
            do k = 1,kk
               tmp3d(:,:,k) = t(:,:,k)*(psl*sig(k)/1.e3)**(-rdry/cp)
            end do
            call vsavehist ( "theta", tmp3d )
         end if
         
         if ( needfld("w") ) then
            do k = 1,kk
               tmp3d(:,:,k) = -(rdry/grav)*t(:,:,k)/(sig(k)*100.*psl) * &
                              ( omega(:,:,k) - sig(k)*dpsdt/864. ) ! Convert dpsdt to Pa/s
            end do
            call vsavehist ( "w", tmp3d )
         end if
         
      end if

      if ( needfld("zg") ) then
         if ( use_plevs ) then
            call height ( t, q, zs, psl, sig, zstd, plevs(1:nplevs) )
            call savehist ( "zg", zstd )
         else if ( use_meters ) then
            do k = 1,nplevs
               zstd(:,:,k) = mlevs(k)/mlevs(nplevs)*(mlevs(nplevs)-zs) + zs
            end do
            call savehist ( "zg", zstd )
         else
            call height ( t, q, zs, psl, sig, zstd )
            call savehist ( "zg", zstd(:,:,minlev:maxlev) )
         end if
      end if
      
      ! Wind vectors
      if ( kk > 1) then
         if ( needfld("u") .or. needfld("v")           .or. &
              needfld("ua") .or. needfld("va")         .or. &             
              needfld("vaveuq") .or. needfld("vavevq") .or. &
              needfld("vaveut") .or. needfld("vavevt") .or. &
              needfld("ubot")   .or. needfld("vbot") ) then
            call fix_winds(u, v)
            if ( cordex_compliant ) then
               if ( needfld("ua") ) then 
                  call vsavehist ( "ua", u )
               end if
               if ( needfld("va") ) then
                  call vsavehist ( "va", v )
               end if   
            else
               if ( needfld("u") ) then 
                  call vsavehist ( "u", u )
               end if
               if ( needfld("v") ) then
                  call vsavehist ( "v", v )
               end if   
            end if
            if ( needfld("ubot") ) then
               call savehist ( "ubot", u(:,:,1) )
            end if
            if ( needfld("vbot") ) then
               call savehist ( "vbot", v(:,:,1) )
            end if   
         end if
      end if       
      
      ! Note that these are just vertical averages, not vertical integrals
      ! Use the winds that have been rotatated to the true directions
      if ( kk>1 ) then
         if ( needfld("vaveuq") ) then
            dtmp = 0.0
            do k = 1,kk
               dtmp = dtmp + dsig(k)*u(:,:,k)*q(:,:,k)
            end do
            call savehist( "vaveuq", dtmp)
         end if
         if ( needfld("vavevq") ) then
            dtmp = 0.0
            do k = 1,kk
               dtmp = dtmp + dsig(k)*v(:,:,k)*q(:,:,k)
            end do
            call savehist( "vavevq", dtmp)
         end if
         if ( needfld("vaveut") ) then
            dtmp = 0.0
            do k = 1,kk
               dtmp = dtmp + dsig(k)*u(:,:,k)*t(:,:,k)
            end do
            call savehist( "vaveut", dtmp)
         end if
         if ( needfld("vavevt") ) then
            dtmp = 0.0
            do k = 1,kk
               dtmp = dtmp + dsig(k)*v(:,:,k)*t(:,:,k)
            end do
            call savehist( "vavevt", dtmp)
         end if
      end if
      
      ! ocean currents
      if ( ok>1 ) then
         if ( needfld("uos") .or. needfld("vos") .or. &
              needfld("uo") .or. needfld("vo") ) then
            call fix_winds(uo_tmp,vo_tmp)
            if ( needfld("uos") ) then
               call savehist( "uos", uo_tmp(:,:,1) )
            end if
            if ( needfld("vos") ) then
               call savehist( "vos", vo_tmp(:,:,1) )
            end if   
            if ( needfld("uo") ) then
               call osavehist( "uo", uo_tmp )
            end if
            if ( needfld("vo") ) then
               call osavehist( "vo", vo_tmp ) 
            end if
            if ( needfld("sos") ) then
               call savehist( "sos", so_tmp(:,:,1) )
            end if
            if ( needfld("tos") ) then
               call savehist( "tos", thetao_tmp(:,:,1) )
            end if   
         end if   
      end if
      
      ! only for TAPM output
      if ( needfld("cos_zen") ) then
         ! calculate zenith angle
         dhr = 1./3600.
         fjd = float(mod(mins, 525600))/1440. ! restrict to 365 day calendar
         ierr = nf90_get_att(ncid, nf90_global, "bpyear", bpyear )
         if ( ierr/=nf90_noerr ) bpyear = 0.
         call solargh(fjd,bpyear,r1,dlt,alp,slag)
         rlat_a(1) = rlat0
         rlong_a(1) = rlong0
         call zenith(fjd,r1,dlt,slag,rlat_a,rlong_a,dhr,1,cos_zen,frac)
         dtmp(:,:) = cos_zen(1)
         call savehist( "cos_zen", dtmp )
      end if
      
      first_in = .false.
      nrec = nrec + 1

      call END_LOG(infile_end)

   end subroutine infile

   function need3dfld(name) result ( needed )
      use history, only : needfld
      character(len=*), intent(in) :: name
      logical :: needed
      select case ( name )
      case ( "temp" )
         needed = needfld("temp") .or. needfld("zg") .or. needfld("rh") .or. &
                  needfld("tbot") .or. needfld("vaveut") .or.                &
                  needfld("vaveut") .or. needfld("theta") .or. needfld("w")
      case ( "ta" )
         needed = needfld("ta") .or. needfld("zg") .or. needfld("rh") .or.  &
                  needfld("tbot") .or. needfld("vaveut") .or.               &
                  needfld("vaveut") .or. needfld("theta") .or. needfld("w")
      case ( "mixr" )
         needed = needfld("mixr") .or. needfld("zg") .or. needfld("rh") .or. &
                  needfld("pwc") .or. needfld("qbot") .or.                   &
                  needfld("vaveuq") .or. needfld("vaveuq") .or.              &
                  needfld("prw")
      case ( "hus" )
         needed = needfld("hus") .or. needfld("zg") .or. needfld("rh") .or. &
                  needfld("pwc") .or. needfld("qbot") .or.                  &
                  needfld("vaveuq") .or. needfld("vaveuq") .or.             &
                  needfld("prw")
      case ( "u", "v" )
         needed = needfld("u") .or. needfld("v") .or. needfld("vaveuq") .or. &
                  needfld("vavevq") .or. needfld("vaveut") .or.              &
                  needfld("vavevt") .or. needfld("vbot") .or.                &
                  needfld("ubot")
      case ( "ua", "va" )
         needed = needfld("ua") .or. needfld("va") .or. needfld("vaveuq") .or. &
                  needfld("vavevq") .or. needfld("vaveut") .or.                &
                  needfld("vavevt") .or. needfld("vbot") .or.                  &
                  needfld("ubot")
     case ( "qlg" )
         needed = needfld("qlg") .or. needfld("rh")
      case ( "qfg" )
         needed = needfld("qfg") .or. needfld("rh")
      case ( "qsng" )
         needed = needfld("qsng")
      case ( "qgrg" )
         needed = needfld("qgrg")
      case ( "omega" )
         needed = needfld("omega") .or. needfld("w")
      case ( "uo", "vo" )
         needed = needfld("uo") .or. needfld("vo") .or. needfld("uos") .or. &
                  needfld("vos")
      case ( "so" ) 
         needed = needfld("so") .or. needfld("sos")
      case ( "thetao" )
         needed = needfld("thetao") .or. needfld("tos")
      case ( "wo" )
         needed = needfld("wo")
      case default
         print*, "Error - unsupported argument for need3dfld", name
         stop
      end select
   end function need3dfld

   subroutine vread2(name, var, required, vread_err)

      use logging_m
      ! Routine to read a variable from either a fortran binary or netcdf file.
      character(len=*), intent(in) :: name
      real, dimension(:,:), intent(out) :: var
      logical, intent(in), optional :: required
      logical req
      integer, intent(out), optional :: vread_err
      integer v_err

      call START_LOG(vread_begin)
      if ( present( required ) .and. present ( vread_err ) ) then
         req = required
      else
         req = .true.
      end if

      call paravar2a(name,var,nrec,req,v_err)
      
      if ( .not. req ) then
         vread_err = v_err
      else if ( present ( vread_err ) ) then
         ! If we get to this point, everything is ok.
         vread_err = 0
      end if
      call END_LOG(vread_end)

   end subroutine vread2
   
   subroutine vread3(name, var)
      use logging_m
      
      ! Routine to read a variable from either a fortran binary or netcdf file. 
      character(len=*), intent(in) :: name
      real, dimension(:,:,:), intent(out) :: var
      integer :: pkl

      call START_LOG(vread_begin)
      pkl = size(var,3)
      call paravar3a(name,var,nrec,pkl)
      call END_LOG(vread_end)

   end subroutine vread3
   
   subroutine readsave2 ( name, save_flag, input_name, required )
!     Read an array from the input file and call savehist
      use history, only : savehist, needfld
      use s2p_m
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: save_flag, required
      character(len=*), intent(in), optional :: input_name
      real, dimension(pil,pjl*pnpan*lproc) :: array
      character(len=50) :: nname
      integer :: ierr

      if ( present(save_flag) ) then
         if ( .not. save_flag ) then
!           If it's not to be saved just skip over the record.
!           This is used for general skipping of everything rather than
!           for specific variables.
            return
         end if
      end if

      if ( .not. needfld(name) ) return
      if ( present(input_name) ) then
         nname = input_name
      else
         nname = name
      end if
      call vread( nname, array, required=required, vread_err=ierr)
      if ( ierr == 0 ) then ! Only write it if it was present in the input file
         call savehist ( name, array )
      end if

   end subroutine readsave2

   subroutine readsave3 ( name, save_flag, input_name )
!     Read an array from the input file and call savehist
      use history, only : savehist, needfld
      use s2p_m
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: save_flag
      character(len=*), intent(in), optional :: input_name
      real, dimension(pil,pjl*pnpan*lproc,kk) :: array
      character(len=50) :: nname

      if ( present(save_flag) ) then
         if ( .not. save_flag ) then
!           If it's not to be saved just skip over the record.
            return
         end if
      end if
      if ( .not. needfld(name) ) return
      if ( present(input_name) ) then
         nname = input_name
      else
         nname = name
      end if
      call vread( nname, array)
      call vsavehist ( name, array )

   end subroutine readsave3
   
   subroutine initialise ( hres, minlon, maxlon, dlon, minlat, maxlat, dlat,  &
                           dx, dy, lx, ly,                                    &
                           kdate, ktime, ntracers, ksoil, kice, debug,        &
                           nqg )

#ifdef usenc_mod
      use netcdf
#else
      use netcdf_m
#endif
      use newmpar_m
      use history
      use xyzinfo_m
      use indices_m
      use latltoij_m
      use setxyz_m
      use interp_m
      use parm_m, only : rlong0, rlat0, schmidt ! Share with final_init
      use physparams, only : erad
      use vertutils_m, only : sig2ds
#ifdef usempi3
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
      use mpidata_m
      use shdata_m      
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif
#endif

      real, intent(inout)  :: hres
      real, intent(inout)  :: minlon, maxlon, dlon, minlat, maxlat, dlat
      real, intent(in)     :: dx, dy
      integer, intent(in)  :: lx, ly
      integer, intent(out) :: kdate, ktime, ntracers, ksoil, kice
      logical, intent(in)  :: debug
      integer, intent(out) :: nqg

      integer, parameter :: diag = 0, id = 1, jd = 1
      real(kind=rx), parameter :: schm13 = 0.1

      integer, parameter :: ntang=2, idiag=0
      integer :: i, j, k
      character(len=10) :: rundate*10
      integer :: m, meso, nx1, nps, mex, mup, nem, nx2, nmi,           &
                 npsav, nhor, nkuo, khdif, nx3, nx4, nvad,                  &
                 nx5, nrun, nrunx, khor, ksc, kountr, ndiur, nhort,         &
                 nhorps, nsoil, ms, ntsur, nrad, kuocb, nvmix, ntsea,       &
                 nonl, nextout, kwt
      real :: difknbd, rhkuo, du, tanl, timer, timeg, dslocal
      integer, dimension(:), allocatable :: int_header
      real, dimension(:), allocatable :: real_header
      integer :: hlen, ierr, vid, dimid, ieof
      integer :: ic, jc
      real :: grlong, grlat, grlongu, grlatu, grlongv, grlatv, tmp
      real :: rdt, lrlat0, lrlong0, lschmidt
      real :: hlon_tmp, hlat_tmp, hlon_dx, hlat_dy

#ifdef usempi3
      integer(kind=MPI_ADDRESS_KIND) :: ssize
      integer :: itest
#endif
!     Read the header here because doing the CC grid initialisation before
!     alloc_indata minimises the total memory requirements

      nrec = 1
      ! Get the total number of timesteps
      ierr = nf90_inq_dimid ( ncid, "time", dimid )
      call check_ncerr(ierr, "Error getting time dimension")
      ierr = nf90_inquire_dimension ( ncid, dimid, len=maxrec )
      call check_ncerr(ierr,"Error getting number of sets")
      ierr = nf90_get_att(ncid, nf90_global, "dt", rdt )
      if ( ierr==nf90_noerr ) then
         ! newer attribute method
         ndt = nint(rdt)
         ierr = nf90_get_att(ncid, nf90_global, "il_g", il )
         call check_ncerr(ierr, "Error getting jl_g attribute")
         ierr = nf90_get_att(ncid, nf90_global, "jl_g", jl )
         call check_ncerr(ierr, "Error getting jl_g attribute")
         ierr = nf90_get_att(ncid, nf90_global, "ms", ms )
         call check_ncerr(ierr, "Error getting ms attribute")
         ierr = nf90_get_att(ncid, nf90_global, "ntrac", ntrac )
         call check_ncerr(ierr, "Error getting ntrac attribute")
         ierr = nf90_inq_dimid(ncid, "lev", dimid )
         call check_ncerr(ierr, "Error getting lev dimension")
         ierr = nf90_inquire_dimension ( ncid, dimid, len=kl )
         call check_ncerr(ierr, "Error getting number of levels")
         ierr = nf90_inq_dimid(ncid, "olev", dimid )
         if ( ierr==nf90_noerr ) then
            ierr = nf90_inquire_dimension ( ncid, dimid, len=ol )
            call check_ncerr(ierr, "Error getting number of ocean levels")
         else
            ol = 0 
         end if    
      else
         ! older int_header method
         ! Get integer and real headers from attibutes. First check the 
         ! lengths of these.
         ierr = nf90_inquire_attribute(ncid, nf90_global, "int_header", len=hlen)
         call check_ncerr(ierr, "Error getting int_header length")
         if ( hlen < 43 ) then
            print*, "Error - insufficient header information: int_header too short"
            stop
         end if
         allocate (int_header(hlen))
         ierr = nf90_get_att(ncid, nf90_global, "int_header", int_header)
         call check_ncerr(ierr, "Error getting int_header")
         ! Only a few values are used
         il = int_header(1)
         jl = int_header(2)
         kl = int_header(3)
         ol = 0
         ndt = int_header(14)
         ms = int_header(34)
         ntrac = int_header(43)
      end if
      ik = il  ! These are only set once
      jk = jl
      kk = kl
      ok = ol
      nsd = 0
      nqg = 0
      ilt = il

      ierr = nf90_get_att(ncid, nf90_global, "rlat0", lrlat0 )
      if ( ierr==nf90_noerr ) then
         ! newer attribute method
         rlat0 = lrlat0
         ierr = nf90_get_att(ncid, nf90_global, "rlong0", lrlong0 )
         rlong0 = lrlong0
         call check_ncerr(ierr, "Error getting rlong0 attribute")
         ierr = nf90_get_att(ncid, nf90_global, "schmidt", lschmidt )
         schmidt = lschmidt
         call check_ncerr(ierr, "Error getting schmidt attribute")
      else
         ! older real_header method
         ierr = nf90_inquire_attribute(ncid, nf90_global, "real_header", len=hlen)
         call check_ncerr(ierr, "Error getting real_header length")
         if ( hlen < 8 ) then
            print*, "Error - insufficient header information: real_header too short"
            stop
         end if
         allocate (real_header(hlen))
         ierr = nf90_get_att(ncid, nf90_global, "real_header", real_header)
         call check_ncerr(ierr, "Error getting real_header")
         ! Only a few values are used
         rlong0 = real_header(5)
         rlat0 = real_header(6)
         schmidt = real_header(7)
         if(schmidt <= 0. .or. schmidt > 1.) then
            !  Some old model initial conditions
            rlong0 = real_header(6)
            rlat0  = real_header(7)
            schmidt = real_header(8)
         endif
      endif
      ! Need date for setting time origin. Call getdate with nrec=1
      nrec = 1
      call getdate(kdate, ktime, ieof)
      if ( myid == 0 ) then
         print *,"- Initialise"
      end if   
      if ( ieof /= 0 ) then
         print*, "Error in initialisation, empty netcdf file"
         stop
      end if

      if ( debug ) then 
         print*, "HEADER "
         write(*,'("kdate",i10," ktime",i6," ktau",i6)') kdate, ktime, ktau
         write(*,'("il",i4," jl",i4," kl",i4)') il, jl, kl
         write(*,'("nsd",i2," nqg",i3," ntrac",i3)') nsd, nqg, ntrac 
!          m,nsd,meso,nx1,nps  &
!         ,mex,mup,nem,nx2,nmi,ndt,npsav,rundate,nhor,nkuo,khdif,kwt         &
!         ,nx3,nx4,timer,timeg,dslocal,nvad,nqg                              &
!         ,nx5,nrun,nrunx,khor,ksc,kountr,ndiur,nhort,nhorps,nsoil           &
!         ,ms,ntsur,nrad,kuocb,nvmix,ntsea,nonl,nextout,ilt,ntrac         &
!         ,difknbd,rhkuo,du,tanl,
         write(*,'("rlong0",f8.2," rlat0",f8.2," schmidt",f6.3)') &
               rlong0, rlat0, schmidt
      end if
      !if ( nqg >= 8 ) then
         ksoil = ms
      !else
      !   ksoil = 2
      !end if
      !if ( nqg >= 11 ) then
         kice = ms
      !else
      !   kice = 0
      !end if


      if ( ilt > 1 ) then
         ntracers = ntrac
      else
         ntracers = 0
      end if

!     Have to read sig here because it's required for openhist.
      allocate ( sig(kl), dsig(kl), zsoil(ksoil), zse(ksoil) )

      if ( kl > 1 ) then
         ! Get sigma levels from level variable
         ierr = nf90_inq_varid (ncid, "lev", vid )
         call check_ncerr(ierr, "Error getting vid for lev")
         ierr = nf90_get_var ( ncid, vid, sig )
         call check_ncerr(ierr, "Error getting levels")
      else
         sig = 1.0
      end if
      call sig2ds(sig, dsig)
      if ( ol > 0 ) then
         allocate( gosig(ol) ) 
         ierr = nf90_inq_varid (ncid, "olev", vid )
         call check_ncerr(ierr, "Error getting vid for olev")
         ierr = nf90_get_var ( ncid, vid, gosig )
         call check_ncerr(ierr, "Error getting ocean levels")
      end if
      ! Note that some initial condition files don't have zsoil
      !if ( cf_compliant ) then
         ierr = nf90_inq_varid (ncid, "zsoil", vid )
         if ( ierr==nf90_noerr ) then
            ierr = nf90_get_var ( ncid, vid, zsoil)
            call check_ncerr(ierr, "Error getting zsoil")
         else
             ksoil = 0 ! missing value flag
             zsoil(:) = 0. 
         end if
      !end if
         
      ! calculate soil thickness   
      zse(1) = 2.*zsoil(1)
      zse(2) = 2.*(zsoil(2)-zse(1))
      do k = 3,ksoil
        zse(k) = 2.*(zsoil(k)-sum(zse(1:k-1)))  
      end do

!     Set all the resolution parameters
      npanels = jl/il - 1
      ifull=il*jl
      ij=il*jl
      ijk=il*jl*kl
      iquad=1+il*((8*npanels)/(npanels+4))

#ifdef usempi3
      if ( node_myid == 0 ) then
#else
      if ( myid == 0 ) then
#endif

          call setxyz ( il, jl, kl, npanels, ifull, iquad, idiag, id, jd,        &
                    rlong0, rlat0, schmidt, schm13, ntang, erad )
                    
!#ifdef usempi3
!      else
!          if ( node_myid == 0 ) then
!              ssize=ifull*10
!          else
!              ssize=0
!          endif
!          call allocshdata(i_ewns,ssize,(/ ifull, 10 /),i_ewns_win)
!          i_w => i_ewns(:,1)
!          i_ww => i_ewns(:,2)
!          i_e => i_ewns(:,3)
!          i_ee => i_ewns(:,4)
!          i_s => i_ewns(:,5)
!          i_ss => i_ewns(:,6)
!          i_n => i_ewns(:,7)
!          i_nn => i_ewns(:,8)
!          i_en => i_ewns(:,9)
!          i_wn => i_ewns(:,10)
!
!          if ( node_myid == 0 ) then
!              ssize=(npanels+1)*10
!          else
!              ssize=0
!          end if
!          call allocshdata(lewns,ssize,(/ npanels + 1, 10 /),lewns_win)
!          lwws(0:npanels) => lewns(:,1)
!          lws(0:npanels) => lewns(:,2)
!          lwss(0:npanels) => lewns(:,3)
!          lees(0:npanels) => lewns(:,4)
!          les(0:npanels) => lewns(:,5)
!          less(0:npanels) => lewns(:,6)
!          lwwn(0:npanels) => lewns(:,7)
!          lwnn(0:npanels) => lewns(:,8)
!          leen(0:npanels) => lewns(:,9)
!          lenn(0:npanels) => lewns(:,10)
!
!          itest = MPI_MODE_NOPRECEDE + MPI_MODE_NOSTORE
!          call MPI_Win_fence(itest,i_ewns_win,ierr)
!          call MPI_Win_fence(itest,lewns_win,ierr)
!!         the shared data is set from node_myid=0 between these fences
!          call MPI_Win_fence(MPI_MODE_NOSUCCEED,i_ewns_win,ierr)
!          call MPI_Win_fence(MPI_MODE_NOSUCCEED,lewns_win,ierr)
!#endif

      end if

      if ( int_default == int_none ) then
         nxhis = il
         nyhis = jl
      else if ( int_default == int_tapm ) then
         if ( lx==0 .or. ly==0 ) then
            print *,"Error, lx, ly, dx, dy need to be specified for TAPM output"
            stop
         end if
         nxhis = lx
         nyhis = ly
      else
!        If hres has been set use that, otherwise calculate a value
!        appropriate to the resolution. 360/(4*il) is the approx resolution
!        along the equator for the cube and 360/(6*il) for the octagon
         if ( hres == 0.0 ) then
            hres = ceiling ( 90.0*schmidt/il )
         end if
         if ( maxlon < minlon ) then
            print*, "Warning: specified maxlon < minlon - swapping"
            tmp = minlon
            minlon = maxlon
            maxlon = tmp
         end if
         ! If dlon set use that, otherwise use hres
         if ( dlon /= 0. ) then
            nxhis = nint ( (maxlon-minlon) / dlon )
         else
            nxhis = nint ( (maxlon-minlon) / hres )
         end if
!        If it's not periodic then need to add 1 to nxhis to get both
!        boundaries
         if ( maxlon - minlon /= 360. ) then
            nxhis = nxhis + 1
         end if
         if ( maxlat < minlat ) then
            print*, "Warning: specified maxlat < minlat - swapping"
            tmp = minlat
            minlat = maxlat
            maxlat = tmp
         end if
         if ( dlat /= 0. ) then
            nyhis = nint ( (maxlat-minlat) / dlat ) + 1
         else
            nyhis = nint ( (maxlat-minlat) / hres ) + 1
         end if
      end if


#ifdef usempi3
      if ( node_myid == 0 ) then
#else
      if ( myid == 0 ) then
#endif

!        To save memory de-allocate a number of arrays defined by setxyz
!        that aren't needed by cc2hist.
         deallocate ( f, fu, fv, dmdx, dmdy, dmdxv, dmdyu )
#ifdef usempi3
      end if

      if ( node_myid == 0 ) then
         ssize = nxhis*nyhis
      else
         ssize = 0
      end if
      call allocshdata(xyg,ssize*2,(/ nxhis, nyhis, 2 /),xyg_win)
      xg => xyg(:,:,1)
      yg => xyg(:,:,2)
      call allocshdata(nface,ssize,(/ nxhis, nyhis /),nface_win)

      itest = MPI_MODE_NOPRECEDE + MPI_MODE_NOSTORE
      call MPI_Win_fence(itest,xyg_win,ierr)
      call MPI_Win_fence(itest,nface_win,ierr)

      if ( node_myid == 0 ) then
#else
         allocate ( nface(nxhis,nyhis) )
         allocate ( xg(nxhis,nyhis), yg(nxhis,nyhis) )
#endif
         allocate ( hlon(nxhis), hlat(nyhis) )
         if ( int_default == int_none ) then
            hlat = (/ ( real(j), j=1,nyhis ) /)
            hlon = (/ ( real(i), i=1,nxhis ) /)
         else if ( int_default == int_tapm ) then
            jc = 0.5*real(ly+1) 
            ic = 0.5*real(lx+1)
            hlat = (/ ( (real(j)-jc)*dy, j=1,nyhis ) /) ! Y coordinate (not latitude)
            hlon = (/ ( (real(i)-ic)*dx, i=1,nxhis ) /) ! X coordinate (not longitude)
            hlat_dy = 1./(6.37e6*3.1415927/180.)
            hlon_dx = hlat_dy/cos(rlat0*3.1415927/180.)
            do j = 1,nyhis
               do i = 1,nxhis
                  hlon_tmp = rlong0 + hlon(i)*hlon_dx
                  hlat_tmp = rlat0 + hlat(j)*hlat_dy
                  if ( hlat_tmp>90. ) then
                     hlat_tmp = 180. - hlat_tmp
                     hlon_tmp = hlon_tmp + 180.
                  end if
                  if ( hlat_tmp<-90. ) then
                     hlat_tmp = -180. + hlat_tmp
                     hlon_tmp = hlon_tmp + 180.
                  end if
                  call latltoij ( hlon_tmp, hlat_tmp, xg(i,j), yg(i,j), nface(i,j),&
                                  rlong0, rlat0, schmidt, schm13 )
               enddo
            enddo 
         else
            !     Set lats and longs
            if ( nyhis == 1 ) then
               hlat(1) = minlat
            else  
               do j = 1,nyhis
                  hlat(j) = minlat + (j-1)*(maxlat-minlat)/(nyhis-1)
               end do
            end if   
            if ( maxlon - minlon == 360.0 ) then
               do i = 1,nxhis
                  hlon(i) = minlon + (i-1)*(maxlon-minlon)/nxhis
               end do
            else
               if ( nxhis == 1 ) then
                  hlon(1) = minlon 
               else    
                  do i = 1,nxhis
                     hlon(i) = minlon + (i-1)*(maxlon-minlon)/(nxhis-1)
                  end do
               end if   
            end if
            do j = 1,nyhis
               do i = 1,nxhis
                  call latltoij ( hlon(i), hlat(j), xg(i,j), yg(i,j), nface(i,j),&
                                  rlong0, rlat0, schmidt, schm13 )
               enddo
            enddo

         end if

!        These aren't needed now.
         deallocate ( xx4, yy4 )

      end if   

#ifdef usempi3
      call MPI_Win_fence(MPI_MODE_NOSUCCEED,xyg_win,ierr)
      call MPI_Win_fence(MPI_MODE_NOSUCCEED,nface_win,ierr)
#endif


   end subroutine initialise

   subroutine final_init(varlist, nvars)
      ! Some final initialisation that requires needfld, so can only be done
      ! after openhist has been called.
      use newmpar_m
      use history, only : needfld
      use xyzinfo_m
      use indices_m
      use parm_m, only : rlong0, rlat0
      use physparams, only : pi
      use logging_m      
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif
      type(input_var), dimension(:) :: varlist
      integer, intent(in) :: nvars
      real, dimension(:,:), allocatable :: costh_g, sinth_g
      real :: sinlong, coslong, sinlat, coslat
      real :: polenx, poleny, polenz, zonx, zony, zonz, den
      real :: theta_lc, rlongdeg, rlatdeg, ri, rj
      integer :: i, j, iq, ivar
      integer :: ierr, ip, n
      logical :: need_rotate
      
      call START_LOG(finalinit_begin)
#ifdef usempi3
      if ( node_myid == 0 ) then
#else
      if ( myid == 0 ) then
#endif
         deallocate ( em )
         deallocate ( i_wu, i_sv, i_eu, i_nv )
      end if

!     There may be extra fields that require rotation, not set in the standard
!     list. Check all names to see if any vector fields are set
!     Still require the test on names for non-netcdf inputs.
      need_rotate = .false.
      do ivar=1,nvars
         if ( varlist(ivar)%vector .and. needfld(varlist(ivar)%vname) ) then
            need_rotate = .true.
            exit
         end if
      end do

      if ( need_rotate .or. needfld("u") .or. needfld("v") .or.   &
           needfld("taux") .or. needfld("tauy") .or.              &
           needfld("u10max") .or. needfld("v10max") .or.          &
           needfld("uas") .or. needfld("vas") .or.                &
           needfld("d10") .or. needfld("ubot")  .or.              &
           needfld("vbot") ) then

         allocate ( costh(pil,pjl*pnpan*lproc), sinth(pil,pjl*pnpan*lproc) )
         
#ifdef usempi3
         if ( node_myid == 0 ) then
#else
         if ( myid == 0 ) then
#endif
         
           allocate ( costh_g(il,jl), sinth_g(il,jl) )

!     For calculating zonal and meridional wind components, use the
!     following information, where theta is the angle between the
!     (ax,ay,az) vector [along the xg axis] and the zonal-component-vector:
!     veczon = k x r, i.e. (-y,x,0)/sqrt(x**2 + y**2)
!     vecmer = r x veczon, i.e. (-xz,-yz,x**2 + y**2)/sqrt(x**2 + y**2)
!     costh is (veczon . a) = (-y*ax + x*ay)/sqrt(x**2 + y**2)
!     sinth is (vecmer . a) = [-xz*ax - yz*ay + (x**2 + y**2)*az]/sqrt
!      using (r . a)=0, sinth collapses to az/sqrt(x**2 + y**2)

!     For rotated coordinated version, see JMcG's notes
            coslong=cos(rlong0*pi/180.)
            sinlong=sin(rlong0*pi/180.)
            coslat=cos(rlat0*pi/180.)
            sinlat=sin(rlat0*pi/180.)
            polenx=-coslat
            poleny=0.
            polenz=sinlat
            do j=1,jl
               do i=1,il
                  iq = i + (j-1)*il
!                 Set up unit zonal vector components
                  zonx = poleny*z(iq)-polenz*y(iq)
                  zony = polenz*x(iq)-polenx*z(iq)
                  zonz = polenx*y(iq)-poleny*x(iq)
!                 Allow for poles by taking max
                  den = sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )
                  costh_g(i,j) =  (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
                  sinth_g(i,j) = -(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
               end do
            end do
            call ccmpi_scatter(costh,costh_g)
            call ccmpi_scatter(sinth,sinth_g)
            deallocate( costh_g, sinth_g )
     
         else
            call ccmpi_scatter(costh)
            call ccmpi_scatter(sinth)
         end if 

      end if
      
#ifdef usempi3
      if ( node_myid == 0 ) then
#else
      if ( myid == 0 ) then
#endif
!       x, y, z, ax, ay, az, bx, by, bz no longer needed.
!       ax etc are pointers to setxyz private arrays so the space can't be freed.
        deallocate ( x, y, z )
      end if
      call END_LOG(finalinit_end)


   end subroutine final_init
   
   subroutine fix_winds2 ( u, v )
!     Convert winds from grid directions to standard zonal and meridional.
      use history, only : int_default
      use interp_m, only : int_none
      use newmpar_m
      real, dimension(:,:), intent(inout) :: u, v
      real, dimension(pil,pjl*pnpan*lproc) :: uzon, vmer

      if ( int_default == int_none ) return      
      
      uzon = costh*u - sinth*v
      vmer = sinth*u + costh*v
!        Now save these back to the original arrays.
         u = uzon
         v = vmer
   end subroutine fix_winds2

   subroutine fix_winds3 ( u, v )
!     Convert winds from grid directions to standard zonal and meridional.
      use newmpar_m
      real, dimension(:,:,:), intent(inout) :: u, v
      integer :: k
      do k=1,size(u,3)
         call fix_winds2( u(:,:,k), v(:,:,k) ) 
      enddo
   end subroutine fix_winds3

   subroutine get_var_list(varlist, nvars)
      ! Get a list of the variables in the input file
      use history, only : addfld, int_default, cf_compliant, cordex_compliant
      use interp_m, only : int_nearest, int_none, int_tapm
      use newmpar_m, only : ol
      use physparams, only : grav, rdry
      use s2p_m, only: use_plevs, plevs, use_meters, mlevs
      type(input_var), dimension(:), pointer, contiguous :: varlist
      integer, intent(out) :: nvars
      integer :: ierr, ndimensions, nvariables, ndims, ivar, int_type, xtype
      integer :: londim, latdim, levdim, olevdim, procdim, timedim, vid, ihr, ind
      integer, dimension(nf90_max_var_dims) :: dimids
      logical :: procformat, ran_type
      character(len=10) :: substr
      character(len=20) :: vname
      character(len=100) :: long_name, tmpname, valid_att, std_name, cell_methods
      ! Perhaps should read these from the input?
      integer, parameter :: vmin=-32500, vmax=32500
      real :: xmin, xmax, aoff, sf, topsig, topheight
      real :: coord_height
      integer :: ivar_start

      ierr = nf90_inquire(ncid, ndimensions=ndimensions, nvariables=nvariables)
      call check_ncerr(ierr, "nf90_inquire error")
      ! This is slightly bigger than required because not all the variables
      ! in the netcdf file are processed.
      allocate ( varlist(nvariables) )

      ! Get the values of the dimension IDs
      ierr = nf90_inq_dimid(ncid, "longitude", londim)
      call check_ncerr(ierr,"Error getting lonid")
      ierr = nf90_inq_dimid(ncid, "latitude", latdim)
      call check_ncerr(ierr,"Error getting latid")
      ierr = nf90_inq_dimid(ncid, "lev", levdim)
      call check_ncerr(ierr,"Error getting levid")
      ierr = nf90_inq_dimid(ncid, "olev", olevdim)
      if ( ierr/=nf90_noerr ) olevdim = -1
      ierr = nf90_inq_dimid(ncid, "processor", procdim) ! only for procformat      
      procformat = (ierr==nf90_noerr)
      ierr = nf90_inq_dimid(ncid, "time", timedim)
      call check_ncerr(ierr,"Error getting timeid")

      nvars = 0
      ! For the sigma to pressure initialisation to work properly
      ! the surface pressure has to be read before any of the 3D variables
      ! Force it to be first by having an extra 0 iteration on the loop
      ! Some special purpose outputs don't have psf, so only treat cases where
      ! use_plevs=T in this way.
      !if ( use_plevs ) then
      !   ivar_start = 0
      !else
         ivar_start = 1
      !end if
      do ivar=ivar_start,nvariables
         !if ( ivar == 0 ) then
         !   ierr = nf90_inq_varid ( ncid, "psf", vid )
         !   call check_ncerr(ierr, "Error getting vid for psf")
         !   ierr = nf90_inquire_variable (ncid, vid, name=vname, ndims=ndims, dimids=dimids, xtype=xtype)
         !   call check_ncerr(ierr, "nf90_inquire_variable error")
         !else
            ierr = nf90_inquire_variable (ncid, ivar, name=vname, ndims=ndims, dimids=dimids, xtype=xtype)
            call check_ncerr(ierr, "nf90_inquire_variable error")
            !if ( use_plevs .and. vname == "psf" ) then
            !   ! This has already been handled above
            !   cycle
            !end if
         !end if
!         print*, ivar, vname, ndims, dimids(1:ndims)
         if ( ndims == 5 .and. procformat ) then   
            ! Should be lon, lat, lev, proc, time
            call check_ncerr(ierr,"Error getting processorid")
            if ( match ( dimids(1:ndims), (/ londim, latdim, levdim, procdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 3  ! Space only
            else if ( match ( dimids(1:ndims), (/ londim, latdim, olevdim, procdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 3
            else
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else if ( ndims == 4 .and. .not.procformat ) then
            ! Should be lon, lat, lev, time
            if ( match ( dimids(1:ndims), (/ londim, latdim, levdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 3  ! Space only
            else if ( match ( dimids(1:ndims), (/ londim, latdim, olevdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 3
            else
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else if ( ndims == 4 .and. procformat ) then
            ! Check for soil variables
            if ( cf_compliant .and. is_soil_var(vname) ) then
               cycle
            end if
            if ( match( dimids(1:ndims), (/ londim, latdim, procdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 2  ! Space only
            else
               ! 3D variables fixed in time aren't supported at the moment
               ! though there's no reason why they couldn't be
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else if ( ndims == 3 .and. .not.procformat ) then
            ! Check for soil variables
            if ( cf_compliant .and. is_soil_var(vname) ) then
               cycle
            end if
            if ( match( dimids(1:ndims), (/ londim, latdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 2  ! Space only
            else
               ! 3D variables fixed in time aren't supported at the moment
               ! though there's no reason why they couldn't be
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else if ( ndims == 3 .and. procformat ) then
            if ( match( dimids(1:ndims), (/ londim, latdim, procdim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .true.
               varlist(nvars)%ndims = 2  ! Space only
            else
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else if ( ndims == 2 .and. procformat ) then
            cycle
         else if ( ndims == 2 ) then
            if ( match( dimids(1:ndims), (/ londim, latdim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .true.
               varlist(nvars)%ndims = 2  ! Space only
            else
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else
            cycle
         end if

         varlist(nvars)%vname = vname
         ierr = nf90_inq_varid(ncid, vname, varlist(nvars)%vid)
         call check_ncerr(ierr, "Error getting vid for "//trim(vname))

         ierr = nf90_get_att(ncid, varlist(nvars)%vid, 'long_name', &
                             varlist(nvars)%long_name)
         call check_ncerr(ierr, "Error getting long_name for "//trim(vname))

         ierr = nf90_get_att(ncid, varlist(nvars)%vid, 'units',     &
                             varlist(nvars)%units)
         ! In initial condition files, some fields do not have units
         if ( vname /= "psf" .and. vname /= "wfg" .and. vname /= "wfb" .and. &
              vname /= "snd" .and. vname /= "fracice" ) then
            call check_ncerr(ierr, "Error getting units for "//trim(vname))
         end if

         if ( xtype == nf90_short ) then
            ierr = nf90_get_att(ncid, varlist(nvars)%vid, 'add_offset',   &
                 varlist(nvars)%add_offset)
            call check_ncerr(ierr, "Error getting add_offset for "//trim(vname))

            ierr = nf90_get_att(ncid, varlist(nvars)%vid, 'scale_factor',     &
                 varlist(nvars)%scale_factor)
            call check_ncerr(ierr, "Error getting scale_factor for "//trim(vname))
         else
            ! Should have some way of setting these. At the moment 16 bit
            ! output won't work with 32 bit input.
            varlist(nvars)%add_offset = 0.
            varlist(nvars)%scale_factor = 1.
         end if

      end do

      if ( cf_compliant ) then
         nvars = nvars + 1
         varlist(nvars)%vname = "tgg"
         varlist(nvars)%fixed = .false.
         varlist(nvars)%ndims = 4
         !nvars = nvars+1
         !varlist(nvars)%vname = "wb"
         !varlist(nvars)%fixed = .false.
         !varlist(nvars)%ndims = 4
         nvars = nvars+1
         varlist(nvars)%vname = "wetfrac"
         varlist(nvars)%fixed = .false.
         varlist(nvars)%ndims = 4
!        Don't need to set extra parameters because there's an explicit 
!        addfld call later.
!         call addfld('tgg','Soil temperature','K',100.,350.,ksoil,soil=.true.)
!         call addfld('wb','Soil moisture','frac',0.,1.,ksoil,soil=.true.)
      end if


      do ivar=1,nvars
         ! Check is variable is a component of a vector
         ! For now test both x-component and x-direction. Eventually
         ! only the former should be needed but older files use the
         ! latter for wind stress.
         if ( index(varlist(ivar)%long_name,"x-comp") /= 0 .or. &
              index(varlist(ivar)%long_name,"x-direction") /= 0) then
            varlist(ivar)%vector = .true.
            varlist(ivar)%xcmpnt = .true.
            ! Search for the matching y-component
            tmpname = varlist(ivar)%long_name
            ind = index(varlist(ivar)%long_name,"x-")
            tmpname(ind:ind) = "y"
            ! Some variables have name like u10 (should be 10m wind).
            ind = index(varlist(ivar)%long_name,"u10")
            if ( ind /= 0 ) tmpname(ind:ind) = "v"
            ! For now handle this specially.
            varlist(ivar)%othercmpnt = match_longname(varlist,tmpname)
            if ( varlist(ivar)%othercmpnt == 0 ) then
               print*, "Failed to find matching vector component for variable ", varlist(ivar)%vname
               stop
            end if
         else if ( index(varlist(ivar)%long_name,"y-comp") /= 0 .or. &
                   index(varlist(ivar)%long_name,"y-direction") /= 0) then
            varlist(ivar)%vector = .true.
            varlist(ivar)%xcmpnt = .false.
         else
            varlist(ivar)%vector = .false.
         end if

      end do

      ! Fix up the names of the vector fields. This comes after the checking
      ! to make sure that it doesn't muck up finding the matching component
      ! Some variables use x-compt rather than x-component which makes things
      ! more complicated.
      do ivar=1,nvars
         if ( varlist(ivar)%vector ) then
            if ( varlist(ivar)%xcmpnt ) then
               ind = index(varlist(ivar)%long_name,"x-comp")
               if ( ind /= 0 ) then
                  ! Search for space so that x-compt and x-component both work
                  ind = index(varlist(ivar)%long_name," ")
                  write( varlist(ivar)%long_name, "(a,a)" ) "Zonal",  &
                       varlist(ivar)%long_name(ind:len_trim(varlist(ivar)%long_name))
               end if
            else
               ind = index(varlist(ivar)%long_name,"y-comp")
               if ( ind /= 0 ) then
                  ind = index(varlist(ivar)%long_name," ")
                  write( varlist(ivar)%long_name, "(a,a)" ) "Meridional",  &
                       varlist(ivar)%long_name(ind:len_trim(varlist(ivar)%long_name))
               end if
            end if
         end if
      end do

      do ivar=1,nvars
         if ( varlist(ivar)%vname == "thetao" .or. varlist(ivar)%vname == "so" .or. &
              varlist(ivar)%vname == "uo" .or. varlist(ivar)%vname == "vo" .or.     &
              varlist(ivar)%vname == "wo" ) then    
            varlist(ivar)%water = .true. 
         else 
            varlist(ivar)%water = .false. 
         end if
      end do   
      
      do ivar=1,nvars
         xmin = varlist(ivar)%add_offset + varlist(ivar)%scale_factor*vmin
         xmax = varlist(ivar)%add_offset + varlist(ivar)%scale_factor*vmax
         ! As a check re-calc offset, scalef
         sf = (xmax - xmin)/(real(vmax)-real(vmin)) ! jlm fix for precision problems
         aoff = xmin - sf*vmin
!         print*, varlist(ivar)%vname, xmin, xmax, varlist(ivar)%add_offset, varlist(ivar)%scale_factor, aoff, sf

         ! Check if variable is only valid once per day. This is true for min
         ! and max variables and those with 3hr, 6hr .. 24hr in their long
         ! names
         varlist(ivar)%daily = .false.
         do ihr=3,24,3
            write(substr,"(i2,a)") ihr, "hr"
            if ( index(varlist(ivar)%long_name,trim(adjustl(substr))) /= 0 ) then
               varlist(ivar)%daily = .true.
               exit
            end if
         end do
         if ( varlist(ivar)%vname == "tmaxscr" .or.         &
              varlist(ivar)%vname == "tminscr" .or.         &
              varlist(ivar)%vname == "maxrnd" .or.          &
              varlist(ivar)%vname == "rhmaxscr" .or.        &
              varlist(ivar)%vname == "rhminscr" .or.        &
              varlist(ivar)%vname == "u10max" .or.          &
              varlist(ivar)%vname == "v10max" .or.          &
              varlist(ivar)%vname == "u1max" .or.           &
              varlist(ivar)%vname == "v1max" .or.           &
              varlist(ivar)%vname == "u2max" .or.           &
              varlist(ivar)%vname == "v2max" .or.           &
              varlist(ivar)%vname == "urbantasmax" .or.     &
              varlist(ivar)%vname == "urbantasmin" ) then
            varlist(ivar)%daily = .true.
         end if
         ! Also set daily if variable has a 'valid_time' attribute with the 
         ! value daily.
         valid_att = ""
         ierr = nf90_get_att(ncid, varlist(ivar)%vid, 'valid_time',valid_att)
         if ( ierr == 0 ) then
            varlist(ivar)%daily = valid_att == "daily"
         end if

         ! Is this really simpler than a string of if tests?
         if ( match ( varlist(ivar)%vname, &
               (/ "vegt                ", "soilt               ", "rsmin               ", "rs                  ", &
                  "zolnd               ", "sigmf               ", "wetfrac             ", "wetfrac?            ", &
                  "tgg?                ", "tgg??               ", "sal??               ", "roadtgg?            ", &
                  "rooftgg?            ", "waletgg?            ", "walwtgg?            ", "dmse_ave            ", &
                  "dmsso2_ave          ", "so2e_ave            ", "so2dd_ave           ", "so2wd_ave           ", &
                  "so2so4_ave          ", "so4e_ave            ", "so4dd_ave           ", "so4wd_ave           ", &
                  "bce_ave             ", "bcdd_ave            ", "bcwd_ave            ", "oce_ave             ", &
                  "ocdd_ave            ", "ocwd_ave            ", "duste_ave           ", "dustdd_ave          ", &
                  "dustwd_ave          ", "wb?_ave             ", "climate_biome       ", "climate_ivegt       ", &
                  "climate_min20       ", "climate_max20       ", "climate_alpha20     ", "climate_agdd5       ", &
                  "climate_gmd         ", "climate_dmoist_min20", "climate_dmoist_max20" /)) .and. int_type /= int_none ) then
            int_type = int_nearest
         else
            int_type = int_default
         end if
         if ( match ( varlist(ivar)%vname, &
               (/ "cape_ave  ", "cape_max  ", "cbas_ave  ", "cfrac     ", "cld       ", "clh       ", "cll       ", &
                  "clm       ", "convh_ave ", "ctop_ave  ", "lwp_ave   ", "maxrnd    ", "omega     ", "pblh      ", &
                  "pmsl      ", "psl       ", "qfg       ", "qgrg      ", "qgscrn    ", "qlg       ", "qrg       ", &
                  "qsng      ", "rfrac     ", "rhscrn    ", "rnc       ", "rnd       ", "rnd24     ", "sfcWindmax", &
                  "sunhours  ", "temp      ", "tscrn     ", "tsu       ", "u         ", "u10       ", "u10max    ", &
                  "uscrn     ", "v         ", "v10max    ", "zolnd     ", "zht       " /)) ) then
            ran_type = .true.
         else
            ran_type = .false.
         end if
         if ( varlist(ivar)%vname == "pmsl" ) then
            varlist(ivar)%vname = "psl"
            if ( cordex_compliant ) then
               varlist(ivar)%units = "Pa"
               xmin = 0.
               xmax = 120000.
            end if
         else if ( varlist(ivar)%vname == "psf" ) then
            cycle  ! Skip this one to avoid messages about it never being set
         else if ( varlist(ivar)%vname(1:5) == "wbice" ) then
            cycle  ! Skip this one to avoid messages about it never being set
         else if ( varlist(ivar)%vname == "zht" ) then
            varlist(ivar)%vname = "zs"
            varlist(ivar)%units = "m"
            varlist(ivar)%long_name = "Surface height"
            xmin = 0.
            xmax = 9000.
         end if
         if ( cordex_compliant ) then
            if ( varlist(ivar)%vname == "cld" ) then
               varlist(ivar)%vname = "clt"
               varlist(ivar)%units = "%"
               xmin = 0.
               xmax = 100.
            else if ( varlist(ivar)%vname == "clh" ) then   
               varlist(ivar)%units = "%"
               xmin = 0.
               xmax = 100.
            else if ( varlist(ivar)%vname == "cll" ) then
               varlist(ivar)%units = "%"
               xmin = 0.
               xmax = 100.
            else if ( varlist(ivar)%vname == "clm" ) then    
               varlist(ivar)%units = "%"
               xmin = 0.
               xmax = 100.
            else if ( varlist(ivar)%vname == "eg_ave" ) then
               varlist(ivar)%vname = "hfls"
            else if ( varlist(ivar)%vname == "epot_ave" ) then
               varlist(ivar)%vname = "evspsblpot"
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Potential Surface Evaporation"
               xmin = 0.
               xmax = 0.001
            else if ( varlist(ivar)%vname == "fg_ave" ) then
               varlist(ivar)%vname = "hfss"
            else if ( varlist(ivar)%vname == "fracice" ) then
               varlist(ivar)%vname = "sic"
               varlist(ivar)%units = "%"
               varlist(ivar)%long_name = "Sea Ice Area Fraction"
               xmin = 0.
               xmax = 100.
            else if ( varlist(ivar)%vname == "iwp_ave" ) then
               varlist(ivar)%vname = "clivi"
            else if ( varlist(ivar)%vname == "lwp_ave" ) then
               varlist(ivar)%vname = "clwvi"
            else if ( varlist(ivar)%vname == "mixr" ) then
               varlist(ivar)%vname = "hus"
            else if ( varlist(ivar)%vname == "mrros" ) then
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Surface runoff"
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "pblh" ) then
               varlist(ivar)%vname = "zmla"
            else if ( varlist(ivar)%vname == "pmsl_ave" ) then
               varlist(ivar)%vname = "psl_ave"
               varlist(ivar)%units = "Pa"
               xmin = 0.
               xmax = 120000.
            else if ( varlist(ivar)%vname == "qgscrn" ) then
               varlist(ivar)%vname = "huss"
               varlist(ivar)%units = "none"
               varlist(ivar)%long_name = "2m specific humidity"
               xmin = 0.
               xmax = 0.06
            else if ( varlist(ivar)%vname == "rgdn_ave" ) then
               varlist(ivar)%vname = "rlds"
            else if ( varlist(ivar)%vname == "rhscrn" ) then
               varlist(ivar)%vname = "hurs"
            else if ( varlist(ivar)%vname == "rnd" ) then
               varlist(ivar)%vname = "pr"
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Precipitation"
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "rnc" ) then
               varlist(ivar)%vname = "prc"
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Convective Precipitation"
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "maxrnd" ) then
               varlist(ivar)%vname = "prmax"
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Maximum precipitation"
               varlist(ivar)%daily = .true.
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "runoff" ) then
               varlist(ivar)%vname = "mrro"
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Total Runoff"
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "rtu_ave" ) then
               varlist(ivar)%vname = "rlut"
            else if ( varlist(ivar)%vname == "sgdn_ave" ) then
               varlist(ivar)%vname = "rsds"
            else if ( varlist(ivar)%vname == "sint_ave" ) then
               varlist(ivar)%vname = "rsdt"
            else if ( varlist(ivar)%vname == "snd" ) then
               varlist(ivar)%long_name = "Snow depth"
            else if ( varlist(ivar)%vname == "snm" ) then
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Snowmelt"
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "sno" ) then
               varlist(ivar)%vname = "prsn"
               varlist(ivar)%units = "kg/m2/s"
               varlist(ivar)%long_name = "Snowfall"
               xmin = 0.
               xmax = 0.013
            else if ( varlist(ivar)%vname == "sot_ave" ) then
               varlist(ivar)%vname = "rsut"
            else if ( varlist(ivar)%vname == "sunhours" ) then
               varlist(ivar)%vname = "sund"
               varlist(ivar)%units = "s"
               varlist(ivar)%long_name = "Sunshine Hours"
               xmin = 0.
               xmax = 86400.
            else if ( varlist(ivar)%vname == "taux" ) then
               varlist(ivar)%vname = "tauu"
            else if ( varlist(ivar)%vname == "tauy" ) then
               varlist(ivar)%vname = "tauv"
            else if ( varlist(ivar)%vname == "temp" ) then
               varlist(ivar)%vname = "ta"
            else if ( varlist(ivar)%vname == "tmaxscr" ) then
               varlist(ivar)%vname = "tasmax"
            else if ( varlist(ivar)%vname == "tminscr" ) then
               varlist(ivar)%vname = "tasmin"
            else if ( varlist(ivar)%vname == "tscrn" ) then
               varlist(ivar)%vname = "tas"
            else if ( varlist(ivar)%vname == "tsu" ) then
               varlist(ivar)%vname = "ts"
            else if ( varlist(ivar)%vname == "u10" ) then
               varlist(ivar)%vname = "sfcWind"
            else if ( varlist(ivar)%vname == "u" ) then
               varlist(ivar)%vname = "ua"
            else if ( varlist(ivar)%vname == "v" ) then
               varlist(ivar)%vname = "va"        
            else if ( varlist(ivar)%vname == "zs" ) then
               varlist(ivar)%vname = "orog"
            end if
         end if
         call cc_cfproperties(varlist(ivar), std_name, cell_methods)
         if ( varlist(ivar)%fixed ) then
            call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name,     &
                      varlist(ivar)%units, xmin, xmax, 1, ave_type="fixed", &
                      int_type=int_type, std_name=std_name,                 &
                      ran_type=ran_type )
         else if ( varlist(ivar)%ndims == 2 ) then
            if ( varlist(ivar)%vname(1:3) == "max" ) then
               ! Special check for maximum rainfall rate
               call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name, &
                      varlist(ivar)%units, xmin, xmax, 1, ave_type="max",  &
                      std_name=std_name, cell_methods=cell_methods,        &
                      ran_type=ran_type )
            else
               ! Check for screen and 10m variables
               coord_height = -huge(1.) ! Acts as a null value
               if ( cf_compliant ) then
                  if ( index(varlist(ivar)%long_name, "screen") /= 0 .or. &
                       index(varlist(ivar)%long_name, "Screen") /= 0 ) then
                     coord_height = 2.
                  end if
                  if ( index(varlist(ivar)%long_name, "10m") /= 0 ) then
                     coord_height = 10.
                  end if
               end if
               call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name,   &
                      varlist(ivar)%units, xmin, xmax, 1, std_name=std_name, &
                      coord_height=coord_height, cell_methods=cell_methods,  & 
                      int_type=int_type, ran_type=ran_type )
            end if
         else if ( varlist(ivar)%ndims == 3 ) then  
            if ( varlist(ivar)%water ) then
              call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name,        &
                        varlist(ivar)%units, xmin, xmax, onlev, multilev=.true., &
                        std_name=std_name, water=varlist(ivar)%water,            &
                        cell_methods=cell_methods, int_type=int_type,            &
                        ran_type=ran_type )
            else
              call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name,       &
                        varlist(ivar)%units, xmin, xmax, nlev, multilev=.true., &
                        std_name=std_name, cell_methods=cell_methods,           &
                        int_type=int_type, ran_type=ran_type )
            end if  
         end if
      end do

      ! Extra fields are handled explicitly
      if ( kk > 1 ) then
         call addfld ( "grid", "Grid resolution", "km", 0., 1000., 1, ave_type="fixed" )
         call addfld ( "tsea", "Sea surface temperature", "K", 150., 350., 1, &
                        std_name="sea_surface_temperature", ran_type=.true. )
         if ( cordex_compliant ) then
            call addfld ( "evspsbl", "Evaporation", "kg/m2/s", 0., 0.001, 1, std_name="water_evaporation_flux" )          
            call addfld ( "mrso", "Total soil moisture content", "kg/m2", 0., 100.0, 1, std_name="soil_moisture_content" )          
            call addfld ( "mrfso", "Soil frozen water content", "kg/m2", 0., 100.0, 1, std_name="soil_frozen_water_content" )   
            call addfld ( "prw", "Precipitable water column", "kg/m2", 0.0, 100.0, 1,  &
                           std_name="atmosphere_water_vapor_content", ran_type=.true. )
            call addfld ( "ps", "Surface pressure", "Pa", 0., 120000., 1, &
                           std_name="surface_air_pressure", ran_type=.true. )
            call addfld ( "rlus", "Upwelling Longwave radiation", "W/m2", -1000., 1000., 1 )
            call addfld ( "rsus", "Upwelling Shortwave radiation", "W/m2", -1000., 1000., 1 )
            if ( int_type /= int_none ) then
               call addfld ( "sftlf", "Land-sea mask", "%",  0.0, 100.0, 1, &
                              ave_type="fixed", int_type=int_nearest,       &
                              ran_type=.true. )
            else
               call addfld ( "sftlf", "Land-sea mask", "%",  0.0, 100.0, 1, &
                              ave_type="fixed", int_type=int_none,          &
                              ran_type=.true. )
            end if
            call addfld ( "snc",  "Snow area fraction", "%", 0., 6.5, 1 )
            call addfld ( "snw",  "Surface Snow Amount", "kg m-2", 0., 6.5, 1 )
         else
            call addfld ( "d10", "10m wind direction", "deg", 0.0, 360.0, 1, ran_type=.true. )          
            if ( int_type /= int_none ) then
               call addfld ( "land_mask", "Land-sea mask", "",  0.0, 1.0, 1, &
                              ave_type="fixed", int_type=int_nearest, ran_type=.true. )
            else
               call addfld ( "land_mask", "Land-sea mask", "",  0.0, 1.0, 1, &
                              ave_type="fixed", int_type=int_none, ran_type=.true. )
            end if
            call addfld ( "ps", "Surface pressure", "hPa", 0., 1200., 1, &
                           std_name="surface_air_pressure", ran_type=.true. )
            call addfld ( "pwc", "Precipitable water column", "kg/m2", 0.0, 100.0, 1, &
                          std_name="atmosphere_water_vapor_content", ran_type=.true. )           
         end if
         call addfld ( "uas", "x-component 10m wind", "m/s", -100.0, 100.0, 1, ran_type=.true. )
         call addfld ( "vas", "y-component 10m wind", "m/s", -100.0, 100.0, 1, ran_type=.true. )
         ! Packing is not going to work well in this case
         ! For height, estimate the height of the top level and use that
         ! for scaling
         if ( use_plevs ) then
            topsig = plevs(nlev) / 1000.
         else if ( use_meters ) then
            topsig = exp(-grav*mlevs(nlev)/(300.*rdry))
         else
            topsig = sig(kk)
         end if
         ! Assume isothermal 300K profile. This will lead to an overestimate
         ! of the height so is safe. Using moist adiabatic lapse rate can
         ! give an underestimate.
         topheight = -rdry*300./grav * log(topsig)
         ! Round to next highest 1000 m
         topheight = 1000*(floor(0.001*topheight)+1)
         call addfld ( "zg", "Geopotential height", "m", 0., topheight, nlev,    &
                        multilev=.true., std_name="geopotential_height",         &
                        ran_type=.true. )
         call addfld ( "press", "Air pressure", "hPa", 0., 1500., nlev,          &
                        multilev=.true., std_name="air_pressure",                &
                        ran_type=.true. )
         call addfld ( "rh", "Relative humidity", "%", 0., 110., nlev,           &
                        multilev=.true., std_name="relative_humidity",           &
                        ran_type=.true. )
         call addfld ( "theta", "Potential temperature", "K", 150., 1200., nlev, &
                        multilev=.true., std_name="potential_temperature" )
         call addfld ( "w", "Vertical velocity", "m/s", -1., 1., nlev,           &
                        multilev=.true., std_name="vertical_velocity" )
         call addfld ( "td", "Dew point temperature", "K", 100.0, 400.0, nlev,   &
                        multilev=.true., ran_type=.true. ) 
         
         ! If the output uses pressure levels save the lowest sigma level of
         ! the basic fields.
         call addfld ( "tbot", "Air temperature at lowest sigma level", "K", 100., 400., 1 )
         call addfld ( "ubot", "Eastward wind at lowest sigma level", "m/s", -100., 100., 1)
         call addfld ( "vbot", "Northward wind at lowest sigma level", "m/s", -100., 100., 1)
         call addfld ( "qbot", "Specific humidity at lowest sigma level", "kg/kg", 0., 0.1, 1, &
                       ran_type=.true. )

         ! Vertical averaged fluxes. Note that these are NOT vertical integrals
         call addfld ( "vaveuq", "Vertical average of zonal humidity flux", "m/s kg/kg", -0.5, 0.5, 1, std_name=std_name )
         call addfld ( "vavevq", "Vertical average of meridional humidity flux", "m/s kg/kg", -0.5, 0.5, 1, std_name=std_name )
         call addfld ( "vaveut", "Vertical average of zonal temperature flux", "m/s K", -1e4, 1e4, 1, std_name=std_name )
         call addfld ( "vavevt", "Vertical average of meridional temperature flux", "m/s K", -1e4, 2e4, 1, std_name=std_name )

         if ( cf_compliant ) then
            ! Define as an extra field for now
            call addfld('tgg','Soil temperature','K',100.,350.,ksoil,soil=.true.)
            ! Should have std_name = volume_fraction_of_water_in_soil, units=1
            ! Mentioned in Gregory email 2005-12-01. In official list?
            !call addfld('wb','Soil moisture','frac',0.,1.,ksoil,soil=.true.)
         end if
      
         if ( int_default == int_tapm ) then
            call addfld('cos_zen', 'Cosine of solar zenith angle', 'none', -1., 1., 1 )
         end if

      else
         ! high-frequency output
         call addfld ( "u10", "10m wind speed", "m/s", 0., 100.0, 1, ran_type=.true. ) 
         call addfld ( "d10", "10m wind direction", "deg", 0.0, 360.0, 1, ran_type=.true. )
      end if
      
      if ( ok > 1 ) then
         call addfld( "uos", "x-component surface current", "m/s", -100., 100., 1 )
         call addfld( "vos", "y-component surface current", "m/s", -100., 100., 1 )
         call addfld( "sos", "Surface ocean salinity", "PSU", 0., 100., 1 )
         call addfld( "tos", "Surface ocean temperature", "K", 150., 350., 1 )
      end if    

   end subroutine get_var_list

   logical function matcha(a1, a2)
      ! Function to check if a1 == a2 (which isn't legal directly)
      integer, dimension(:), intent(in) :: a1, a2
      integer :: i

      if ( size(a1) /= size(a2) ) then
         print*, "Error in match, lengths are different"
         stop
      end if
      matcha = .true.
      do i=1,size(a1)
         if ( a1(i) /= a2(i) ) matcha = .false.
      end do
   end function matcha

   logical function matchc(str, strarray)
      ! Function to check if str matches one of strarray
      ! Allows ? as single character wildcard
      character(len=*), intent(in) :: str
      character(len=*), dimension(:), intent(in) :: strarray
      integer :: i, j, len, start, jtest
      logical :: tmpmatch

      matchc = .false.
      do i=1,size(strarray)
         if ( len_trim(strarray(i)) /= len_trim(str) ) then
            cycle
         end if
         if (index(strarray(i), "?") > 0 ) then
            ! Treat this as single character wild card and do the match 
            ! excluding these
            len = len_trim(strarray(i))
            tmpmatch = .true.
            j = 0
            start = 1
            do 
               ! Look for next ?
               ! j is index into whole string
               jtest = index(strarray(i)(j+1:len),"?")
               j = j + jtest
               if ( jtest == 0 ) then
                  ! No more wildcards so match right to end
                  j = len+1
               end if
               ! Match from start to j
               tmpmatch = tmpmatch .and. strarray(i)(start:j-1) == str(start:j-1)
               start = j + 1
               if (start > len ) then
                  exit
               end if
            end do
            if (tmpmatch) then
               matchc = .true.
               exit
            end if
         else
            if ( strarray(i) == str ) then
               matchc = .true.
               exit
            end if
         end if
      end do
   end function matchc

   integer function match_longname(varlist,name)
      type(input_var), dimension(:) :: varlist
      character(len=*), intent(in) :: name
      integer :: ivar
      match_longname = 0
      do ivar=1,size(varlist)
         if ( name == varlist(ivar)%long_name ) then
            match_longname = ivar
            exit
         end if
      end do
   end function match_longname
   
   subroutine calc_rh ( t, q, ql, qf, psl, sig, rh )
      use moistfuncs
      real, dimension(pil,pjl*pnpan*lproc,kk), intent(in) :: t, q, ql, qf
      real, dimension(pil,pjl*pnpan*lproc), intent(in) :: psl
      real, dimension(kk), intent(in) :: sig
      real, dimension(pil,pjl*pnpan*lproc,kk), intent(out) :: rh
      real, dimension(pil,pjl*pnpan*lproc) :: p, tliq, fice, qsw, qc
      real, dimension(pil,pjl*pnpan*lproc) :: qsi, qsl, deles
      real, parameter :: tfrz  = 273.1
      real, parameter :: tice  = 233.15
      real, parameter :: cp    = 1004.64
      real, parameter :: hl    = 2.5104e6
      real, parameter :: hlf   = 3.36e5
      real, parameter :: epsil = 0.622
      integer :: k

      do k = 1,kk
         p = 100.*sig(k)*psl  ! p in Pa.
         qc = qf(:,:,k) + ql(:,:,k)
         where (qf(:,:,k)>1.E-12)
            fice = min(qf(:,:,k)/qc,1.)
         elsewhere
            fice = 0.
         end where
         tliq = t(:,:,k) - hl/cp*qc - hlf/cp*qf(:,:,k)
         qsi = qsati(p,tliq)
         deles = esdiffx(tliq)
         qsl = qsi + epsil*deles/p
         qsw = fice*qsi + (1.-fice)*qsl
         rh(:,:,k) = 100.*q(:,:,k)/qsw
         !rh(:,:,k) = 100.*relhum(p,q(:,:,k),t(:,:,k))
      end do
   end subroutine calc_rh

   subroutine calc_td ( t_in, q_in, ql_in, qf_in, psl_in, sig, td )
      real, dimension(pil,pjl*pnpan*lproc,kk), intent(out) :: td
      real, dimension(pil,pjl*pnpan*lproc,kk), intent(in) :: t_in, q_in, ql_in, qf_in
      real, dimension(pil,pjl*pnpan*lproc), intent(in) :: psl_in
      real, dimension(kk), intent(in) :: sig
      real, dimension(pil,pjl*pnpan*lproc) :: p, esl, esi, wsl, wsi
      real, dimension(pil,pjl*pnpan*lproc) :: ws, el, ei, cl, ci, tdl, tdi
      real, dimension(pil,pjl*pnpan*lproc) :: fice, t, q, ql, qf
      real, parameter :: tfrz = 273.1
      real, parameter :: tice = 233.16
      integer :: k
      
      do k = 1,kk
         p = min( max( 100.*psl_in*sig(k), 1.e3), 1.e7) ! p is Pa
         t = min( max( t_in(:,:,k), 100.), 400.)
         q = min( max( q_in(:,:,k), 1.e-10), 1.)
         ql = min( max( ql_in(:,:,k), 1.e-20), 1.)
         qf = min( max( qf_in(:,:,k), 1.e-20), 1.)
         where ( t>=tfrz )
            fice = 0.
         elsewhere ( t>=tice .and. qf>1.e-12 )
            fice = min(qf/(qf+ql), 1.)
         elsewhere( t>=tice )
            fice = 0.
         elsewhere
            fice = 1.
         end where
         esl = 611.2*exp(17.67*(t-273.15)/(t-29.65)) 
         esi = 611.2*exp(21.8745584*(t-273.15)/(t-7.66)) 
         wsl = (287.09/461.5)*esl/max(p-esl,0.1)
         wsi = (287.09/461.5)*esi/max(p-esi,0.1)
         ws = (1.-fice)*wsl + fice*wsi
         el = esl*q/ws
         ei = esi*q/ws
         cl = log(el/611.2)
         ci = log(ei/611.2)
         tdl = (17.67*273.15-29.65*cl)/(17.67-cl) 
         tdi = (21.8745584*273.15-7.66*ci)/(21.8745584-ci)
         ! MJT approximation
         where ( t_in(:,:,k)/=NF90_FILL_FLOAT .and. q_in(:,:,k)/=NF90_FILL_FLOAT .and.   &
                 qf_in(:,:,k)/=NF90_FILL_FLOAT .and. ql_in(:,:,k)/=NF90_FILL_FLOAT .and. &
                 psl_in(:,:)/=NF90_FILL_FLOAT )
            td(:,:,k) = (1.-fice)*tdl + fice*tdi
         elsewhere
            td(:,:,k) = NF90_FILL_FLOAT
         end where
      end do
      
   end subroutine calc_td
   
   subroutine check_cc2histfile()
      ! Check whether the input file was created by cc2hist. Without this
      ! check the error message that this gives rise to is rather obscure.
      integer :: ierr, attlen
      character(len=1000) :: source

      ! Look for the source attribute
      source = ""
      ierr = nf90_inquire_attribute ( ncid, nf90_global, "source", len=attlen)
      if ( ierr == nf90_noerr ) then
         if ( attlen < len(source) ) then
            ! This will be true if it came from cc2hist
            ierr = nf90_get_att ( ncid, nf90_global, "source", source )
            call check_ncerr(ierr, "Error getting source attribute")
            if ( index(source, "Processed by cc2hist") > 0 ) then
               print*, "Error - the input file is already processed by cc2hist"
               stop
            end if
         end if
      end if
   end subroutine check_cc2histfile

   subroutine cc_cfproperties(vinfo, stdname, cell_methods)
      use history, only : cf_compliant
      type(input_var), intent(in) :: vinfo
      character(len=80), intent(out) :: stdname
      character(len=80), intent(out) :: cell_methods
      character(len=20) :: vname

      ! Some fields like rnd03 can't really be made CF compliant. Their time
      ! bounds don't match those of the overall file properly.

      stdname = ""
      cell_methods = ""
      if ( .not. cf_compliant ) return
      ! Would some sort of external table be better
      ! Also return preferred variable name and units?
      vname = vinfo%vname
      select case (vname)
      case("ps")
         stdname = "surface_air_pressure"
      case("psl")
         stdname = "air_pressure_at_sea_level"
      case("u")
         stdname = "eastward_wind"
      case("v")
         stdname = "northward_wind"
      case("eg")
         stdname = "surface_upward_latent_heat_flux"
         cell_methods = "time: mean"
      case("fg")
         stdname = "surface_upward_sensible_heat_flux"
         cell_methods = "time: mean"
      case("zg")
         stdname = "geopotential_height"
      case ("zs")
         stdname = "surface_altitude"
      case ("alb")
         stdname = "surface_albedo"
      case ("cld")
         stdname = "cloud_area_fraction"
      case ("cfrac")
         stdname = "cloud_area_fraction_in_atmosphere_layer"
      case ("qfg")
         stdname = "mass_fraction_of_cloud_ice_in_air"
      case ("qlg")
         stdname = "mass_fraction_of_cloud_liquid_water_in_air"
      case ("rnc")
         stdname = "convective_precipitation_flux"
         cell_methods = "time: mean"
      case ("rnd")
         stdname = "precipitation_flux"
         cell_methods = "time: mean"
      case ("rnd24")
         stdname = "precipitation_flux"
         cell_methods = "time: mean"
      case ("snd")
         stdname = "SNOW ???????????????"
      case ("sno")
         stdname = "snowfall_flux"
      case ("tsu")
         stdname = "surface_temperature" 
      case ("u10")
         ! Perhaps poor choice in the model. This is wind speed rather than
         ! a component.
         stdname = "wind_speed"
      case ("evap")
         ! Over what period?
         stdname = "water_evaporation_amount"
      case ("mixr")
         stdname = "humidity_mixing_ratio"
      case ("pblh")
         stdname = "atmosphere_boundary_layer_thickness"
      case ("taux")
         stdname = "surface_downward_eastward_stress"
      case ("tauy")
         stdname = "surface_downward_northward_stress"
      case ("temp")
         stdname = "air_temperature"
      case ("omega")
         stdname = "lagrangian_tendency_of_air_pressure"
      case ("siced")
         stdname = "sea_ice_thickness"
      case ("sigmf")
         stdname = "vegetation_area_fraction"
      case ("tscrn")
         stdname = "air_temperature" 
      case ("tminscr")
         stdname = "air_temperature" 
         cell_methods = "time: minimum"
      case ("tmaxscr")
         stdname = "air_temperature" 
         cell_methods = "time: maximum"
      case ("tscr_ave")
         stdname = "air_temperature" 
         cell_methods = "time: mean"
      case ("uscrn")
         stdname = "wind_speed" 
      case ("zolnd")
         stdname = "surface_roughness_length" 
      case ("runoff")
         stdname = "surface_runoff_flux"
      case ("rgdn_ave")
         stdname = "surface_downwelling_longwave_flux_in_air"
      case ("sgdn_ave")
         stdname = "surface_downwelling_shortwave_flux_in_air"
      case ("sgn_ave")
         stdname = "surface_net_downward_shortwave_flux"
      case ("sint_ave")
         stdname = "toa_incoming_shortwave_flux"
      case ("cbas_ave")
         stdname = "air_pressure_at_cloud_base"
      case ("ctop_ave")
         stdname = "air_pressure_at_cloud_top"
      end select

      if ( vinfo%daily ) then
         cell_methods = cell_methods(1:len_trim(cell_methods)) // " over days"
      end if

   end subroutine cc_cfproperties

   logical function is_soil_var(vname) 
      character(len=*), intent(in) :: vname
      character(len=10) :: tmpname
      integer :: k

      is_soil_var = .false.
      do k = 1,ksoil
         write(tmpname,'(a,i1)') 'tgg', k
         if ( vname == tmpname ) then
            is_soil_var = .true.
            return
         end if
      end do
      !do k=1,ksoil
      !   write(tmpname,'(a,i1)') 'wb', k
      !   if ( vname == tmpname ) then
      !      is_soil_var = .true.
      !      return
      !   end if
      !end do
      do k = 1,ksoil
         write(tmpname,'(a,i1)') 'wetfrac', k
         if ( vname == tmpname ) then
            is_soil_var = .true.
            return
         end if
      end do
   end function is_soil_var

   ! From ccam infile.f
   subroutine fill_cc(b_io,value)
!     routine fills in interior of an array which has undefined points
      use logging_m
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif
      real, intent(inout) :: b_io(pil,pjl*pnpan*lproc)         ! input and output array
      real, intent(in)    :: value                             ! array value denoting undefined
      real, dimension(0,0,0) :: c_io
      integer :: ierr, lsize
      
      call START_LOG(fillcc_begin)
      if ( myid == 0 ) then
         call fill_cc0(b_io,value)
      else
         lsize = pil*pjl*pnpan*lproc
         call START_LOG(mpigather_begin)
         call MPI_Gather(b_io,lsize,MPI_REAL,c_io,lsize,MPI_REAL,0,comm_world,ierr)
         call END_LOG(mpigather_end)
         call START_LOG(mpiscatter_begin)
         call MPI_Scatter(c_io,lsize,MPI_REAL,b_io,lsize,MPI_REAL,0,comm_world,ierr)
         call END_LOG(mpiscatter_end)
      end if
      call END_LOG(fillcc_end)
      
   end subroutine fill_cc

   subroutine fill_cc0(b_io,value)
!     routine fills in interior of an array which has undefined points
      use newmpar_m
      use indices_m
      use logging_m      
#ifdef usempi_mod
      use mpi
#else
      include 'mpif.h'
#endif
      real, dimension(pil,pjl*pnpan*lproc), intent(inout) :: b_io ! input and output array
      real, intent(in)    :: value                                ! array value denoting undefined
      real, dimension(pil,pjl*pnpan,pnproc) :: c_io
      real, dimension(ifull) :: b, a
      integer, dimension(0:5) :: imin, imax, jmin, jmax
      integer iminb, imaxb, jminb, jmaxb
      integer :: nrem, iq, neighb, ierr, i, j
      integer :: ip, n, lsize
      real :: av, avx
      
      call START_LOG(fillcc0_begin)
      
      lsize = pil*pjl*pnpan*lproc
      call START_LOG(mpigather_begin)
      call MPI_Gather(b_io,lsize,MPI_REAL,c_io,lsize,MPI_REAL,0,comm_world,ierr)
      call END_LOG(mpigather_end)

      do ip = 0,pnproc-1   
         do n = 0,pnpan-1
            do j = 1,pjl
               do i = 1,pil
                  iq = i + ioff(ip,n) + (j+joff(ip,n)+n*pil_g-1)*il
                  a(iq) = c_io(i,j+n*pjl,ip+1)
               end do
            end do
         end do
      end do
      
      imin(:) = 1
      imax(:) = il
      jmin(:) = 1
      jmax(:) = il
      
      nrem = 1    ! Just for first iteration
!     nrem_gmin used to avoid infinite loops, e.g. for no sice
      do while ( nrem > 0 )
         nrem = 0
         b(:) = a(:)
         do n = 0,5
            iminb = il
            imaxb = 1
            jminb = il
            jmaxb = 1
            do j = jmin(n),jmax(n)
               do i = imin(n),imax(n)
                  iq = i + (j-1)*il + n*il*il
                  if ( a(iq) == value ) then
                     neighb = 0
                     av = 0.
                     if ( a(i_n(iq)) /= value ) then
                        neighb = neighb + 1
                        av = av + a(i_n(iq))
                     end if
                     if ( a(i_e(iq)) /= value ) then
                        neighb = neighb + 1
                        av = av + a(i_e(iq))
                     end if
                     if ( a(i_w(iq)) /= value ) then
                        neighb = neighb + 1
                        av = av + a(i_w(iq))
                     end if
                     if ( a(i_s(iq)) /= value ) then
                        neighb = neighb + 1
                        av = av + a(i_s(iq))
                     end if
                     if ( neighb > 0 ) then
                        b(iq) = av/neighb
                        avx = av
                     else
                        nrem = nrem + 1   ! current number of points without a neighbour
                        iminb = min( i, iminb )
                        imaxb = max( i, imaxb )
                        jminb = min( j, jminb )
                        jmaxb = max( j, jmaxb )
                     end if
                  end if
               end do
            end do
            imin(n) = iminb
            imax(n) = imaxb
            jmin(n) = jminb
            jmax(n) = jmaxb
         end do
         a(:) = b(:)
         if ( nrem == ifull ) then
            print*, "Error in fill_cc - no points defined"
            exit
         end if
      end do
 
      do ip = 0,pnproc-1   
         do n = 0,pnpan-1
            do j = 1,pjl
               do i = 1,pil
                 iq = i+ioff(ip,n) + (j+joff(ip,n)+n*pil_g-1)*il
                 c_io(i,j+n*pjl,ip+1) = a(iq)
               end do
            end do
         end do
      end do
      
      call START_LOG(mpiscatter_begin)
      call MPI_Scatter(c_io,lsize,MPI_REAL,b_io,lsize,MPI_REAL,0,comm_world,ierr)
      call END_LOG(mpiscatter_end)
      call END_LOG(fillcc0_end)
      
   end subroutine fill_cc0
   
   subroutine paraopen(ifile,nmode,ncid)
  
#ifdef usempi3
      use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
#endif
#ifdef usempi_mod
      use mpi
#endif
#ifdef usempi3
      use mpidata_m
      use shdata_m
#endif
      use logging_m
#ifndef usempi_mod
      include 'mpif.h'
#endif
  
      integer, intent(in) :: nmode
      integer, intent(out) :: ncid
      integer, dimension(5) :: jdum
      integer, dimension(54) :: int_header
      integer :: ier, ip, n, rip, ierr
      integer :: colour, vid
      integer, dimension(:), allocatable, save :: resprocmap_inv
      integer, dimension(:), allocatable, save :: procfileowner
      character(len=*), intent(in) :: ifile
      character(len=266) :: pfile, old_pfile
      character(len=8) :: sdecomp
      logical :: singlefile

#ifdef usempi3
      integer(kind=MPI_ADDRESS_KIND) :: ssize
      integer :: itest
#endif

      call START_LOG(paraopen_begin)

      nproc_orig = nproc
      
      if ( myid==0 ) then      
  
         ! parallel file input
         ip = 0
         singlefile = .false.
         write(pfile,"(a,'.',i6.6)") trim(ifile), ip
         ierr = nf90_open(pfile, nmode, ncid)
         if ( ierr /= nf90_noerr ) then
            old_pfile = pfile  
            write(pfile,"(a,'.',i4.4)") trim(ifile), ip
            ierr = nf90_open(pfile, nmode, ncid)
         end if  
         if ( ierr /= nf90_noerr ) then
            pfile = ifile
            singlefile = .true.
            ierr = nf90_open(pfile, nmode, ncid)
            call check_ncerr(ierr, "Error opening file "//trim(old_pfile)//" or "//trim(pfile))
            ierr = nf90_get_att(ncid, nf90_global, "nproc", pnproc)
            if ( ierr==nf90_noerr ) then
               write(6,*) "ERROR: Parallel file format found in ifile = ",trim(ifile)
               write(6,*) "       Try removing .000000 from ifile in the namelist."
               write(6,*) "       Alternatively, this could be a lat/lon output that"
               write(6,*) "       has already been processed by pcc2hist."
               stop
            end if
         end if    
         
         
         write(6,*) "Opening ",trim(pfile)
         if ( .not.singlefile ) then
            write(6,*) "Using parallel input files"
         end if   
      
         ! parallel metadata
         if ( .not.singlefile ) then
            ier = nf90_get_att(ncid, nf90_global, "nproc", pnproc)
            call check_ncerr(ier, "nproc")
            ier = nf90_get_att(ncid, nf90_global, "procmode", resprocmode)
            if ( ier==nf90_noerr ) then
               resprocformat = .true.  
            end if
         else
            pnproc = 1 
            resprocmode = 0
            resprocformat = .false. 
         end if    

         if ( resprocformat ) then
            allocate( resprocmap_inv(0:pnproc-1) )
            ierr = nf90_inq_varid (ncid, 'gprocessor', vid ) 
            call check_ncerr(ierr, "Error getting vid for gprocessor")
            ierr = nf90_get_var ( ncid, vid, resprocmap_inv, start=(/ 1 /), count=(/ pnproc /) )
            call check_ncerr(ierr, "Error getting var gprocessor")
         end if

!#ifndef usempi3
!         allocate( ioff(0:pnproc-1,0:5), joff(0:pnproc-1,0:5) )
!#endif
      
      end if
      
      call START_LOG(mpibcast_begin)
      if ( myid==0 ) then
        jdum(1) = pnproc
        if ( resprocformat ) then
          jdum(2) = 1 ! indicates true for resprocformat
          jdum(3) = resprocmode
        else
          jdum(2) = 0 ! indicates false for resprocformat
          jdum(3) = 0
        end if
      end if
      call MPI_Bcast(jdum(1:3), 3, MPI_INTEGER, 0, comm_world, ier)
      pnproc = jdum(1)
      resprocformat = (jdum(2)==1)
      if ( resprocformat ) then
         resprocmode = jdum(3)
      end if
      call END_LOG(mpibcast_end)

      
      if ( mod(pnproc,nproc) /= 0 ) then
          
         if ( myid == 0 ) then
            write(6,'(x,a,i0,a,i0,a)') "WARNING: Number of processors(",nproc,&
                                     ") is not a factor of the number of files(",pnproc,")"
         end if
         
         do n = nproc,1,-1
            if ( mod(pnproc,n) == 0 ) then
               nproc = n
               exit
            end if
         end do
         
         if ( myid == 0 ) then
            write(6,'(x,a,i0)') "WARNING: Using pcc2hist with the following number of processes: ",nproc
         end if

         if (myid < nproc ) then
            colour = 0
         else
            colour = 1
         end if
         call MPI_Comm_split(comm_world, colour, myid, comm_reduced, ierr)
         call MPI_Comm_size(comm_reduced, nproc_reduced, ierr)
         call MPI_Comm_rank(comm_reduced, myid_reduced, ierr)

         !redefine comm_world & myid
         comm_world = comm_reduced
         myid = myid_reduced
         nproc = nproc_reduced

         !exit color=1 ranks
         if ( colour == 1 ) then
            call MPI_Finalize(ierr)
            stop
         end if
         
      end if

#ifdef usempi3
      !redefine the per node communicator
      call MPI_Comm_split_type(comm_world, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, node_comm, ierr) ! Per node communictor
      call MPI_Comm_size(node_comm, node_nproc, ierr) ! Find number of processes on node
      call MPI_Comm_rank(node_comm, node_myid, ierr)  ! Find local processor id on node
#endif

      lproc = pnproc/nproc !number of files each mpi_proc will work on      
      allocate( ncid_in(0:lproc-1) )
      allocate( fown_in(0:lproc-1) )      
      fown_in(:) = .false.
       
!#ifdef usempi3
!      if ( node_myid == 0 ) then
!          ssize = pnproc*6*2
!      else
!          ssize = 0
!      end if
!      call allocshdata(ijoff,ssize,(/ pnproc, 6, 2 /),ijoff_win)
!      ioff(0:pnproc-1,0:5) => ijoff(:,:,1)
!      joff(0:pnproc-1,0:5) => ijoff(:,:,2)
!#endif
      
      allocate( ioff(0:pnproc-1,0:5), joff(0:pnproc-1,0:5) )

      if ( resprocformat ) then
         allocate( prid_in(0:lproc-1) )
         allocate( procfileowner(0:pnproc-1) )
         procfileowner(:) = -1
         if ( .not.allocated(resprocmap_inv) ) then
            allocate( resprocmap_inv(0:pnproc-1) )
         end if
         call MPI_Bcast(resprocmap_inv, pnproc, MPI_INTEGER, 0, comm_world, ier)
         do ip = 0,lproc-1
            rip = myid*lproc + ip
            rip = resprocmap_inv(rip)
            prid_in(ip) = mod(rip, resprocmode) + 1
         end do
      end if

      if ( myid /= 0 ) then
 
         if ( resprocformat ) then
            do ip = 0,lproc-1
               rip = myid*lproc + ip
               rip = resprocmap_inv(rip)
               rip = rip/resprocmode
               if ( procfileowner(rip) == -1 ) then
                  fown_in(ip) = .true.
                  procfileowner(rip) = ip
                  write(pfile,"(a,'.',i6.6)") trim(ifile), rip
                  ier = nf90_open ( pfile, nmode, ncid_in(ip) )
                  if ( ier /= nf90_noerr ) then
                     old_pfile = pfile
                     write(pfile,"(a,'.',i4.4)") trim(ifile), rip 
                     ier = nf90_open ( pfile, nmode, ncid_in(ip) )
                     if (ier /= nf90_noerr ) then
                        write(6,*) "ERROR: Cannot open ",trim(old_pfile)," or ",trim(pfile)
                        call check_ncerr(ier, "open")
                     end if                     
                  end if
               else
                  ncid_in(ip) = ncid_in(procfileowner(rip))
               end if
            end do
            deallocate(procfileowner)
            deallocate(resprocmap_inv)
         else
            do ip = 0,lproc-1
               rip = myid*lproc + ip
               fown_in(ip) = .true.
               write(pfile,"(a,'.',i6.6)") trim(ifile), rip
               ier = nf90_open ( pfile, nmode, ncid_in(ip) )
               if (ier /= nf90_noerr ) then
                  old_pfile = pfile
                  write(pfile,"(a,'.',i4.4)") trim(ifile), rip
                  ier = nf90_open ( pfile, nmode, ncid_in(ip) )
                  if (ier /= nf90_noerr ) then
                     write(6,*) "ERROR: Cannot open ",trim(old_pfile)," or ",trim(pfile)
                     call check_ncerr(ier, "open")
                  end if   
               end if
            end do
         end if
         ncid = ncid_in(0) 
      
      else
      
         ncid_in(0) = ncid
         fown_in(0) = .true.
         if ( resprocformat ) then
            procfileowner(0) = 0
            do ip = 1,lproc-1
               rip = resprocmap_inv(ip)
               rip = rip/resprocmode
               if ( procfileowner(rip) == -1 ) then
                  fown_in(ip) = .true.
                  procfileowner(rip) = ip
                  write(pfile,"(a,'.',i6.6)") trim(ifile), rip
                  ier = nf90_open ( pfile, nmode, ncid_in(ip) )
                  if ( ier /= nf90_noerr ) then
                     old_pfile = pfile
                     write(pfile,"(a,'.',i4.4)") trim(ifile), rip
                     ier = nf90_open ( pfile, nmode, ncid_in(ip) )
                     if ( ier /= nf90_noerr ) then
                        write(6,*) "ERROR: Cannot open ",trim(old_pfile)," or ",trim(pfile)
                        call check_ncerr(ier, "open")
                     end if   
                  end if
               else
                  ncid_in(ip) = ncid_in(procfileowner(rip)) 
               end if
            end do
            deallocate(procfileowner)
            deallocate(resprocmap_inv)
         else
            do ip = 1,lproc-1
               rip = ip
               fown_in(ip) = .true.
               write(pfile,"(a,'.',i6.6)") trim(ifile), rip
               ier = nf90_open ( pfile, nmode, ncid_in(ip) )
               if ( ier /= nf90_noerr ) then
                  old_pfile = pfile
                  write(pfile,"(a,'.',i4.4)") trim(ifile), rip
                  ier = nf90_open ( pfile, nmode, ncid_in(ip) )
                  if ( ier /= nf90_noerr ) then
                     write(6,*) "ERROR: Cannot open ",trim(old_pfile)," or ",trim(pfile)
                     call check_ncerr(ier, "open")
                  end if   
               end if
            end do
         end if
      
         ier = nf90_get_att(ncid, nf90_global, "il_g", pil_g )
         if ( ier==nf90_noerr ) then
            ierr = nf90_get_att(ncid, nf90_global, "jl_g", pjl_g )
            call check_ncerr(ierr, "Error getting jl_g attribute")
         else
            ! backwards compatibility option -------- 
            !  Get dimensions from int_header
            ier = nf90_get_att(ncid_in(0), nf90_global, "int_header", int_header)
            call check_ncerr(ier, "int_header")
            ! Only a few values are used
            pil_g = int_header(1)
            pjl_g = int_header(2)
         end if

         !  Calculate il, jl from these
         if ( .not.singlefile ) then
            sdecomp = ''
            ier = nf90_get_att(ncid_in(0), nf90_global, "decomp", sdecomp)
            call check_ncerr(ier, "decomp")
         else   
            sdecomp = 'face'
         end if    
!#ifdef usempi3
!      end if
!      call START_LOG(mpibcast_begin)
!      call MPI_Bcast(pil_g,1,MPI_INTEGER,0,comm_world,ierr)
!      call MPI_Bcast(pjl_g,1,MPI_INTEGER,0,comm_world,ierr)
!      call MPI_Bcast(sdecomp,8,MPI_CHAR,0,comm_world,ierr)
!      call END_LOG(mpibcast_end)
!      itest = MPI_MODE_NOPRECEDE + MPI_MODE_NOSTORE
!      call MPI_Win_fence(itest,ijoff_win,ierr)
!      if ( node_myid==0 ) then
!#endif
         select case(sdecomp)
            case ("uniform")
               do n = 0,5
                  call proc_setup_uniform(pil_g, pjl_g, pnproc, n, pil, pjl, ioff(:,n), &
	                                  joff(:,n), pnpan)
               end do
            case ("uniform1")
               do n = 0,5
                  call proc_setup_dix(pil_g, pjl_g, pnproc, n, pil, pjl, ioff(:,n), &
	                                  joff(:,n), pnpan)
               end do
            case ("face")
               call proc_setup(pil_g, pjl_g, pnproc, 0, pil, pjl, ioff(:,0), joff(:,0), pnpan)
               do n = 1,5
                  ioff(:,n) = ioff(:,0)
                  joff(:,n) = joff(:,0)
               end do
            case default
               print *,"ERROR: Unknown decomposition ",trim(sdecomp)
               stop
         end select
      
         jdum(1) = pil
         jdum(2) = pjl
         jdum(3) = pnpan
         jdum(4) = pil_g
         jdum(5) = pjl_g
      
      end if
!#ifdef usempi3
!      call MPI_Win_fence(MPI_MODE_NOSUCCEED,ijoff_win,ierr)
!#endif
      
      call START_LOG(mpibcast_begin)
      call MPI_Bcast(ioff(0:pnproc-1,0:5),pnproc*6,MPI_INTEGER,0,comm_world,ierr)
      call MPI_Bcast(joff(0:pnproc-1,0:5),pnproc*6,MPI_INTEGER,0,comm_world,ierr)
      call MPI_Bcast(jdum(1:5),5,MPI_INTEGER,0,comm_world,ier)
      call END_LOG(mpibcast_end)
      pil   = jdum(1)
      pjl   = jdum(2)
      pnpan = jdum(3)
      pil_g = jdum(4)
      pjl_g = jdum(5)
      
      call END_LOG(paraopen_end)
      
   end subroutine paraopen
   
   subroutine paraclose
      use logging_m
      integer ip, ierr
      
      call START_LOG(paraclose_begin)
      if ( allocated(prid_in) ) then
         deallocate(prid_in)
      end if
      do ip = 0,lproc-1
         if ( fown_in(ip) ) then
            ierr = nf90_close(ncid_in(ip))
         end if
      end do
      deallocate(fown_in)
      deallocate(ncid_in)
      call END_LOG(paraclose_end)
   
   end subroutine paraclose
   
   subroutine paravar2a(name,var,nrec,required,vread_err)
      use logging_m   
      integer, intent(in) :: nrec
      integer, intent(out) :: vread_err
      integer ip, n, vid, ierr, vartyp
      real, dimension(:,:), intent(out) :: var
      real, dimension(pil,pjl*pnpan) :: inarray2
      real addoff, sf
      logical, intent(in) :: required
      character(len=*), intent(in) :: name
   
      ! If the variable has the required flag set to false, and the 
      ! vread_err argument is present, then return an error flag rather
      ! then abort if the variable isn't found.
      if ( .not. required ) then
         vread_err = NF90_NOERR
         return
      end if
      
      call START_LOG(paravar2a_begin)
      
      do ip = 0,lproc-1
         
         ierr = nf90_inq_varid (ncid_in(ip), name, vid ) 
         call check_ncerr(ierr, "Error getting vid for "//name)
      
         if ( resprocformat ) then
            ierr = nf90_get_var ( ncid_in(ip), vid, inarray2(:,:), start=(/ 1, 1, prid_in(ip), nrec /), &
                                  count=(/ pil, pjl*pnpan, 1, 1 /) )
            call check_ncerr(ierr, "Error getting var "//name)
         else
            ierr = nf90_get_var ( ncid_in(ip), vid, inarray2(:,:), start=(/ 1, 1, nrec /), count=(/ pil, pjl*pnpan, 1 /) )
            call check_ncerr(ierr, "Error getting var "//name)
         end if

!        Check the type of the variable
         ierr = nf90_inquire_variable ( ncid_in(ip), vid, xtype=vartyp)
         if ( vartyp == NF90_SHORT ) then
            if ( all( inarray2 == -32501. ) ) then
               inarray2 = NF90_FILL_FLOAT
            else
               ierr = nf90_get_att ( ncid_in(ip), vid, "add_offset", addoff )
               call check_ncerr(ierr, "Error getting add_offset attribute")
               ierr = nf90_get_att ( ncid_in(ip), vid, "scale_factor", sf )
               call check_ncerr (ierr,"Error getting scale_factor attribute")
               inarray2 = addoff + inarray2*sf
            end if
         end if
   
         var(:,1+ip*pjl*pnpan:(ip+1)*pjl*pnpan) = inarray2(:,:)
         
      end do
      
      call END_LOG(paravar2a_end)
   
   end subroutine paravar2a

   subroutine paravar3a(name,var,nrec,pkl)
      use s2p_m, only : minlev, maxlev
      use logging_m
      integer, intent(in) :: nrec, pkl
      integer ip, n, vid, ierr, vartyp, k
      real, dimension(:,:,:), intent(out) :: var
      real, dimension(pil,pjl*pnpan,pkl) :: inarray3
      real addoff, sf
      character(len=*), intent(in) :: name

      call START_LOG(paravar3a_begin)

      do ip = 0,lproc-1

         ierr = nf90_inq_varid ( ncid_in(ip), name, vid )
         call check_ncerr(ierr, "Error getting vid for "//name)

         if ( resprocformat ) then
            ierr = nf90_get_var ( ncid_in(ip), vid, inarray3, start=(/ 1, 1, 1, prid_in(ip), nrec /), &
                                  count=(/ pil, pjl*pnpan, pkl, 1, 1 /) )
         else
            ierr = nf90_get_var ( ncid_in(ip), vid, inarray3, start=(/ 1, 1, 1, nrec /), &
                                  count=(/ pil, pjl*pnpan, pkl, 1 /) )
         end if
         call check_ncerr(ierr, "Error getting var "//name)
      
         ierr = nf90_inquire_variable ( ncid_in(ip), vid, xtype=vartyp )
         if ( vartyp == NF90_SHORT ) then
            if ( all( inarray3 == -32501. ) ) then
               inarray3 = NF90_FILL_FLOAT
            else
               ierr = nf90_get_att ( ncid_in(ip), vid, "add_offset", addoff )
               call check_ncerr(ierr, "Error getting add_offset attribute")
               ierr = nf90_get_att ( ncid_in(ip), vid, "scale_factor", sf )
               call check_ncerr (ierr,"Error getting scale_factor attribute")                
               inarray3 = addoff + inarray3*sf
            end if
         end if
      
         var(:,1+ip*pjl*pnpan:(ip+1)*pjl*pnpan,:) = inarray3(:,:,:)
         
      end do
      
      call END_LOG(paravar3a_end)
   
   end subroutine paravar3a

   subroutine proc_setup(il_g,jl_g,nproc,nin,il,jl,ioff,joff,npan)
      integer, intent(in) :: il_g, jl_g, nproc, nin
      integer, intent(out) :: il, jl, npan
      integer, dimension(0:) :: ioff, joff
      integer, parameter :: npanels=5
      integer :: i, j, n, nxproc, nyproc

      if ( nproc <= npanels+1 ) then
         if ( modulo(npanels+1,nproc) /= 0 ) then
            print*, "Error, number of processors must divide number of panels"
            stop
         end if
         npan = (npanels+1)/nproc
         il = il_g
         jl = il_g*npan
         nxproc = 1
         nyproc = nproc
         npan = 1 ! Global value
      else  ! nproc >= npanels+1
         if ( modulo (nproc, npanels+1) /= 0 ) then
            print*, "Error, number of processors must be a multiple of number of panels"
            stop
         end if
         npan = 1
         n = nproc / (npanels+1)
         !  n is the number of processors on each face
         !  Try to factor this into two values are close as possible.
         !  nxproc is the smaller of the 2.
         nxproc = nint(sqrt(real(n)))
         do nxproc = nint(sqrt(real(n))), 1, -1
            nyproc = n / nxproc
            if ( modulo(il_g,nxproc) == 0 .and. modulo(il_g,nyproc) == 0 .and. &
                 nxproc*nyproc == n ) exit
         end do
         nyproc = n / nxproc
         if ( nxproc*nyproc /= n ) then
            print*, "Error in splitting up faces"
            stop
         end if
         il = il_g/nxproc
         jl = il_g/nyproc

         nyproc = 6*nyproc ! Global value for offset calculation
      end if

      ! Offsets
      n = 0
      do j = 0, nyproc-1
         do i = 0, nxproc-1
            ioff(n) = i*il
            joff(n) = j*jl
            n = n + 1
         end do
      end do

   end subroutine proc_setup

   subroutine proc_setup_uniform(il_g,jl_g,nproc,nin,il,jl,ioff,joff,npan)
      integer, intent(in) :: il_g, jl_g, nproc, nin
      integer, intent(out) :: il, jl, npan
      integer, dimension(0:) :: ioff, joff
      integer :: i, j, n, nxproc, nyproc

      npan=6
      nxproc = nint(sqrt(real(nproc)))
      do nxproc = nint(sqrt(real(nproc))), 1, -1
         ! This will always exit eventually because it's trivially true 
         ! for nxproc=1
         nyproc = nproc / nxproc
         if ( modulo(nproc,nxproc) == 0 .and. &
              modulo(il_g,nxproc) == 0  .and. &
              modulo(il_g,nyproc) == 0 ) exit
      end do
      nyproc = nproc / nxproc
      if ( nxproc*nyproc /= nproc ) then
         print*, "Error in splitting up faces"
         stop
      end if
      il = il_g/nxproc
      jl = il_g/nyproc

      ! Offsets
      do n = 0,nproc-1
        select case(nin)
           case(0)
              joff(n) = (n/nxproc) * jl
              ioff(n) = modulo(n,nxproc)*il
              if (joff(n).ge.il_g/2) then
                 ioff(n)=il_g-ioff(n)-il
              end if
           case(1)
              joff(n) = (n/nxproc) * jl
              ioff(n) = modulo(n,nxproc)*il
           case(2)
              joff(n) = (n/nxproc) * jl
              ioff(n) = modulo(n,nxproc)*il
              if (ioff(n).ge.il_g/2) then
                 joff(n)=il_g-joff(n)-jl
              end if
           case(3)
              joff(n) = modulo(n,nyproc)*jl
              ioff(n) = (n/nyproc) * il
              if (ioff(n).ge.il_g/2) then
                 joff(n)=il_g-joff(n)-jl
              end if
           case(4)
              joff(n) = modulo(n,nyproc)*jl
              ioff(n) = (n/nyproc) * il
           case(5)
              joff(n) = modulo(n,nyproc)*jl
              ioff(n) = (n/nyproc) * il
              if (joff(n).ge.il_g/2) then
                 ioff(n)=il_g-ioff(n)-il
              end if        
        end select
      end do

   end subroutine proc_setup_uniform

   subroutine proc_setup_dix(il_g,jl_g,nproc,nin,il,jl,ioff,joff,npan)
      integer, intent(in) :: il_g, jl_g, nproc, nin
      integer, intent(out) :: il, jl, npan
      integer, dimension(0:) :: ioff, joff
      integer :: i, j, n, nxproc, nyproc

      npan=6
      nxproc = nint(sqrt(real(nproc)))
      do nxproc = nint(sqrt(real(nproc))), 1, -1
         ! This will always exit eventually because it's trivially true 
         ! for nxproc=1
         nyproc = nproc / nxproc
         if ( modulo(nproc,nxproc) == 0 .and. &
              modulo(il_g,nxproc) == 0  .and. &
              modulo(il_g,nyproc) == 0 ) exit
      end do
      nyproc = nproc / nxproc
      if ( nxproc*nyproc /= nproc ) then
         print*, "Error in splitting up faces"
         stop
      end if
      il = il_g/nxproc
      jl = il_g/nyproc

      ! Offsets
      do n = 0,nproc-1
         joff(n) = (n/nxproc) * jl
         ioff(n) = modulo(n,nxproc)*il
      end do

   end subroutine proc_setup_dix

   subroutine cc2hist_work_close
   
   use interp_m
#ifdef usempi3   
   use shdata_m 
#endif

#ifdef usempi3
   call freeshdata(xyg_win)
   call freeshdata(nface_win)
   !call freeshdata(ijoff_win)
   nullify(xyg)
   !nullify(ijoff)
   nullify(xg,yg)
   nullify(nface)
#else
   if ( allocated(xg) ) then
      deallocate(xg,yg)
      deallocate(nface)
   end if   
#endif
   
   end subroutine cc2hist_work_close

   subroutine ccmpi_scatter_1d_r4_host(data_l,data_g)
   
   use logging_m
#ifdef usempi_mod
   use mpi
#endif
   use newmpar_m, only : il, jl
   
   implicit none
   
#ifndef usempi_mod
   include 'mpif.h'
#endif
   
   integer ip, n, iq_a, iq_b, i, j, ierr
   real, dimension(pil*pjl*pnpan*lproc), intent(out) :: data_l
   real, dimension(il*jl), intent(in) :: data_g
   real, dimension(pil*pjl*pnpan,pnproc) :: data_s
   
   do ip = 0,pnproc-1   
      do n = 0,pnpan-1
         do j = 1,pjl
            do i = 1,pil
               iq_a = i + pil*(j-1) + n*pil*pjl
               iq_b = i + ioff(ip,n) + il*(j-1+joff(ip,n)) + n*il*il
               data_s(iq_a,ip+1) = data_g(iq_b)
            end do
         end do
      end do
   end do   

   call START_LOG(mpiscatter_begin)
   call MPI_Scatter(data_s,pil*pjl*pnpan*lproc,MPI_REAL,data_l,pil*pjl*pnpan*lproc,MPI_REAL,0,comm_world,ierr)
   call END_LOG(mpiscatter_end)
   
   end subroutine ccmpi_scatter_1d_r4_host
   
   subroutine ccmpi_scatter_1d_r4_proc(data_l)
   
   use logging_m
#ifdef usempi_mod
   use mpi
#endif
   use newmpar_m, only : il, jl
   
   implicit none
   
#ifndef usempi_mod
   include 'mpif.h'
#endif
   
   integer ierr
   real, dimension(pil*pjl*pnpan*lproc), intent(out) :: data_l
   real, dimension(0,0) :: data_s
   
   call START_LOG(mpiscatter_begin)
   call MPI_Scatter(data_s,pil*pjl*pnpan*lproc,MPI_REAL,data_l,pil*pjl*pnpan*lproc,MPI_REAL,0,comm_world,ierr)
   call END_LOG(mpiscatter_end)
   
   end subroutine ccmpi_scatter_1d_r4_proc
   
   subroutine ccmpi_scatter_2d_r4_host(data_l,data_g)
   
   use logging_m
#ifdef usempi_mod
   use mpi
#endif
   use newmpar_m, only : il, jl
   
   implicit none
   
#ifndef usempi_mod
   include 'mpif.h'
#endif
   
   integer :: ip, n, ierr
   real, dimension(pil,pjl*pnpan*lproc), intent(out) :: data_l
   real, dimension(il,jl), intent(in) :: data_g
   real, dimension(pil,pjl*pnpan,pnproc) :: data_s
   
   do ip = 0,pnproc-1   
      do n = 0,pnpan-1
         data_s(1:pil,1+n*pjl:(n+1)*pjl,ip+1) = &
              data_g(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g)
      end do
   end do

   call START_LOG(mpiscatter_begin)
   call MPI_Scatter(data_s,pil*pjl*pnpan*lproc,MPI_REAL,data_l,pil*pjl*pnpan*lproc,MPI_REAL,0,comm_world,ierr)
   call END_LOG(mpiscatter_end)

   
   end subroutine ccmpi_scatter_2d_r4_host         
         
   subroutine ccmpi_scatter_2d_r4_proc(data_l)

   use logging_m
#ifdef usempi_mod
   use mpi
#endif
   use newmpar_m, only : il, jl
   
   implicit none
   
#ifndef usempi_mod
   include 'mpif.h'
#endif
   
   integer ierr
   real, dimension(pil,pjl*pnpan*lproc), intent(out) :: data_l
   real, dimension(0,0,0) :: data_s
   
   call START_LOG(mpiscatter_begin)
   call MPI_Scatter(data_s,pil*pjl*pnpan*lproc,MPI_REAL,data_l,pil*pjl*pnpan*lproc,MPI_REAL,0,comm_world,ierr)
   call END_LOG(mpiscatter_end)
   
   end subroutine ccmpi_scatter_2d_r4_proc         
         
         
end module work
