module work

   use netcdf
   use ncutils_m, only : check_ncerr
   use gldata
   use precis_m, only : rx
   implicit none

   integer, parameter :: un_in = 10 ! Input unit number
   logical, private, save :: first_in = .true.

!  Passed between infile and main, previously in common infil
   integer :: ik, jk, kk
   integer :: nlev  ! Number of levels in output
   integer :: ndt, ktau, nqg, nsd
   integer :: ksoil, kice
   integer :: ntrac, ilt ! Extra flag controlling tracers in file. Still needed?
   real, private :: rlong0, rlat0, schmidt

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

!  Parallel input
   integer, dimension(:), save, allocatable, private :: ncid_in
   integer, dimension(:,:), save, allocatable, private :: ioff, joff
   integer, save :: pil, pjl, pnpan, pnproc, lproc

   real, allocatable, dimension(:), public, save :: hlon, hlat
   integer, public, save :: nxhis, nyhis
   real, dimension(:,:), allocatable, private :: costh, sinth
   real, dimension(:,:), allocatable, private :: costh_lc, sinth_lc

   public :: initialise, fix_winds, final_init, check_cc2histfile
   public :: paraopen, paraclose
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

   type input_var
      character(len=10) :: vname   ! Name used in input
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
      integer :: othercmpnt ! Index of matching x or y component
   end type input_var

contains

   subroutine alloc_indata ( il, jl, kl, ksoil, kice )

!     Allocates arrays for the variables that have to be stored rather than
!     being immediately processed.
      use history, only : needfld
      use s2p_m

      integer, intent(in) :: il, jl, kl, ksoil, kice
      
      ! MJT notes - possibly allocate for pil and pjl*pnpan for lproc files
      ! instead of using global arrays

      allocate ( psl(il,jl), pmsl(il,jl), zs(il,jl), tsu(il,jl), soilt(il,jl) )
      allocate ( u(il,jl,kl), v(il,jl,kl), t(il,jl,kl), q(il,jl,kl) )
      allocate ( pwc(il,jl), ql(il,jl,kl), qf(il,jl,kl) )
      allocate ( tgg(il,jl,ksoil), wetfrac(il,jl,ksoil), wbice(il,jl,kice) )
      allocate ( snowvar(il,jl,3), isflag(il,jl), tscrn3hr(il, jl, 8) )
      allocate ( taux(il,jl), tauy(il,jl) )
      if ( needfld("zg") ) then
         if ( use_plevs ) then
            allocate ( zstd(il,jl,nplevs) )
         else
!           Sigma level height calculation requires all levels.
            allocate ( zstd(il,jl,kl) )
         end if
      end if
      if ( needfld("rh") ) then
         allocate ( rh(il,jl,kl) )
      end if

   end subroutine alloc_indata

   subroutine getdate ( kdate, ktime, ieof) 

!     Get record data from header
      use mpidata_m
      integer, intent(out) :: kdate, ktime, ieof
      integer :: ik, jk, kk, ierr, vid

      if ( nrec > maxrec ) then
         ieof = 1
         return
      end if
      ! Get vid and then values for kdate, ktime, ktau
      ierr = nf90_inq_varid (ncid, "kdate", vid )
      call check_ncerr(ierr, "Error getting kdate id")
      ierr = nf90_get_var ( ncid, vid, kdate, start=(/ nrec /) )
      call check_ncerr(ierr, "Error getting kdate")
      ierr = nf90_inq_varid (ncid, "ktime", vid )
      call check_ncerr(ierr, "Error getting ktime id")
      ierr = nf90_get_var ( ncid, vid, ktime, start=(/ nrec /) )
      call check_ncerr(ierr, "Error getting ktime")
      ! Get ktau from time. Really should be renamed
      ierr = nf90_inq_varid (ncid, "time", vid )
      call check_ncerr(ierr, "Error getting time id")
      ierr = nf90_get_var ( ncid, vid, ktau, start=(/ nrec /) )
      call check_ncerr(ierr, "Error getting time")
      ieof = 0
            
      if ( myid == 0 ) then
         write(*,"(a,3i8)",advance="no") " kdate, ktime, ktau ", &
                                           kdate,  ktime, ktau
      end if

   end subroutine getdate

   subroutine infile ( varlist, nvars, skip )
      ! For netcdf input
      use mpidata_m
      use history, only : savehist, needfld
      use physparams, only : grav
      use s2p_m
      use height_m
      use sitop_m
      type(input_var), dimension(:) :: varlist
      integer, intent(in) :: nvars
      logical, intent(in)  :: skip
      integer :: k, ivar
      real, dimension(ik,jk) :: vave, lmask, tss_s, utmp, dtmp, mtmp
      real, dimension(ik,jk) :: wind_norm
      real, dimension(ik,jk,kk) :: ttmp
      character(len=10) :: name
      real, parameter :: spval=999.
      real, parameter :: tfreeze = 271.38

      ! With a netcdf file there's no need to read in order to skip.
      if ( skip ) then
         nrec = nrec + 1
         return
      end if

      if ( first_in ) then
         call alloc_indata ( ik, jk, kk, ksoil, kice )
      end if

      do ivar=1,nvars
         ! Just write the input with no further processing

         ! Files with fixed type only written on first pass
         if ( varlist(ivar)%fixed ) then
            if ( first_in ) then
               if ( varlist(ivar)%vname == "zs" ) then
                  ! This could also be done with the output_scale
                  call vread( "zht", zs )
                  zs = zs / grav
                  call savehist ( "zs", zs )
               else if ( varlist(ivar)%vname == "soilt" ) then
                  call vread( "soilt", soilt )
                  call savehist("soilt", soilt)
                  if ( needfld("land_mask") ) then
                     where ( soilt > 0.5 )
                        lmask = 1.
                     elsewhere
                        lmask = 0.
                     end where
                     call savehist("land_mask", lmask)
                  endif
               else
                  call readsave2 ( varlist(ivar)%vname )
               end if
            end if
            cycle ! Otherwise will match the following if test as well
         end if

        ! Only read these on a full day. ktau is really time in minutes
         if ( varlist(ivar)%daily .and. modulo(ktau,1440) /= 0 ) cycle

         if ( varlist(ivar)%ndims == 2 ) then
            select case ( varlist(ivar)%vname )
            case ( "psl" )
               call readsave2 (varlist(ivar)%vname, input_name="pmsl")
            case ( "tsu" )
               call vread( "tsu", tsu )
               ! Some (all?) initial conditions have tsu negative over ocean
               call savehist("tsu", abs(tsu))
               if ( needfld("tsea") ) then
                  ! Use soilt as a land-sea mask (integer but read as float)
                  where ( soilt > 0.5 )
                     tss_s = spval
                  elsewhere
                     ! Use the maxval to ignore ice points.
                     tss_s = max(tfreeze, abs(tsu))
                  end where
                  call fill_cc(tss_s, spval)
                  call savehist ( "tsea", tss_s )
               end if
             case ( "psf" )
               call vread( "psf", psl )
               psl = 1.0e3 * exp(psl)   ! hPa
               call savehist ( "ps", psl )
               ! This relies on surface pressure coming before the 3D variables
               if ( use_plevs ) call sitop_setup(sig, plevs(1:nplevs), psl)
             case ( "u10" )
               call vread( "u10", utmp )
               call savehist( "u10", utmp )
            case default
               if ( varlist(ivar)%vector) then
                  if (varlist(ivar)%xcmpnt) then
                     ! Process both x and y now, skip over y when we come to 
                     ! it later. Use taux and tauy as temporary arrays
                     ! name of other component
                     name = varlist(varlist(ivar)%othercmpnt)%vname
                     if ( needfld(varlist(ivar)%vname) .or. needfld(name) ) then
                        call vread( varlist(ivar)%vname, taux)
                        call vread( name, tauy)
                        call fix_winds(taux, tauy)
                        call savehist ( varlist(ivar)%vname, taux )
                        call savehist ( name, tauy )
                     end if
                  endif
               else
                  call readsave2 (varlist(ivar)%vname)
               end if
            end select
         else
            select case ( varlist(ivar)%vname )
            case ( "temp" )
               if ( need3dfld("temp")) then
                  call vread( "temp", t)
                  call vsavehist ( "temp", t )
               end if
            case ( "mixr" )
               if ( need3dfld("mixr")) then
                  call vread( "mixr", q )
                  call vsavehist ( "mixr", q )
               end if
            case ( "qlg" )
               if ( need3dfld("qlg")) then
                  call vread( "qlg", ql )
                  call vsavehist ( "qlg", ql )
               end if
            case ( "qfg" )
               if ( need3dfld("qfg")) then
                  call vread( "qfg", qf )
                  call vsavehist ( "qfg", qf )
               end if
            ! Should to u, v as above with vector flag, but this will do for now
            case ( "u" )
               if ( need3dfld("u")) then
                  call vread( "u", u )
               end if
            case ( "v" )
               if ( need3dfld("v")) then
                  call vread( "v", v )
               end if
            case ( "tgg" )
               ! Only in cf_compliant mode
               do k=1,ksoil
                  write(name,'(a,i1)') 'tgg', k
                  call vread(name,tgg(:,:,k))
               end do
               call savehist("tgg", tgg)
            !case ( "wb" )
            !   ! Only in cf_compliant mode
            !   do k=1,ksoil
            !      write(name,'(a,i1)') 'wb', k
            !      call vread(name,wb(:,:,k))
            !   end do
            !   call savehist("wb", wb)
            case ( "wetfrac" )
               ! Only in cf_compliant mode
               do k=1,ksoil
                  write(name,'(a,i1)') 'wetfrac', k
                  call vread(name,wetfrac(:,:,k))
               end do
               call savehist("wetfrac", wetfrac)
            case default
               call readsave3 (varlist(ivar)%vname)
            end select
         end if

      end do

      call savehist( "tbot", t(:,:,1))
      call savehist( "qbot", q(:,:,1))

      if ( needfld("u") .or. needfld("v")           .or. &
           needfld("vaveuq") .or. needfld("vavevq") .or. &
           needfld("vaveut") .or. needfld("vavevt") .or. &
           needfld("ubot")   .or. needfld("vbot")   .or. &
           needfld("uas")    .or. needfld("vas")    .or. &
           needfld("d10") ) then
         call fix_winds(u, v)
         call vsavehist ( "u", u )
         call vsavehist ( "v", v )
         call savehist ( "ubot", u(:,:,1) )
         call savehist ( "vbot", v(:,:,1) )
          if ( needfld("d10") ) then
            dtmp=atan2(-u(:,:,1),-v(:,:,1))*180./3.1415927
            where (dtmp.lt.0.)
              dtmp=dtmp+360.
            end where
            call savehist( "d10", dtmp )
          end if

          if ( needfld("uas") .or. needfld("vas") ) then
             wind_norm(:,:) = sqrt(u(:,:,1)*u(:,:,1)+v(:,:,1)*v(:,:,1))
          end if
          if ( needfld("uas") ) then

            where ( wind_norm > 0.0 )
               dtmp=u(:,:,1)*utmp/wind_norm
            elsewhere
               dtmp=0.0
            end where

            call savehist ( "uas", dtmp )
          end if

          if ( needfld("vas") ) then
            where ( wind_norm > 0.0 )
               dtmp=v(:,:,1)*utmp/wind_norm
            elsewhere
               dtmp=0.0
            end where
            call savehist ( "vas", dtmp )
          end if
      end if

      if ( needfld("pwc") ) then
         pwc = 0.0
         do k=1,kk
            pwc = pwc + dsig(k)*q(:,:,k)
         end do
         pwc = 100.*psl/grav * pwc
         call savehist ( "pwc", pwc )
      end if

      if ( needfld("zg") ) then
         if ( use_plevs ) then
            call height ( t, q, zs, psl, sig, zstd, plevs(1:nplevs) )
            call savehist ( "zg", zstd )
         else
            call height ( t, q, zs, psl, sig, zstd )
            call savehist ( "zg", zstd(:,:,minlev:maxlev) )
         end if
      end if

      if ( needfld("rh") ) then
         call calc_rh ( t, q, ql, qf, psl, sig, rh )
         call vsavehist ( "rh", rh )
      end if

      ! Note that these are just vertical averages, not vertical integrals
      ! Use the winds that have been rotatated to the true directions
      if ( needfld("vaveuq") ) then
         vave = 0.0
         do k=1,kk
            vave = vave + dsig(k)*u(:,:,k)*q(:,:,k)
         end do
         call savehist( "vaveuq", vave)
      end if
      if ( needfld("vavevq") ) then
         vave = 0.0
         do k=1,kk
            vave = vave + dsig(k)*v(:,:,k)*q(:,:,k)
         end do
         call savehist( "vavevq", vave)
      end if
      if ( needfld("vaveut") ) then
         vave = 0.0
         do k=1,kk
            vave = vave + dsig(k)*u(:,:,k)*t(:,:,k)
         end do
         call savehist( "vaveut", vave)
      end if
      if ( needfld("vavevt") ) then
         vave = 0.0
         do k=1,kk
            vave = vave + dsig(k)*v(:,:,k)*t(:,:,k)
         end do
         call savehist( "vavevt", vave)
      end if

      first_in = .false.
      nrec = nrec + 1

   end subroutine infile

   function need3dfld(name) result ( needed )
      use history, only : needfld
      character(len=*), intent(in) :: name
      logical :: needed
      select case ( name )
      case ( "temp" )
         needed = needfld("temp") .or. needfld("zg") .or. needfld("rh") .or. &
                  needfld("tbot") .or. needfld("vaveut") .or. needfld("vaveut")
      case ( "mixr" )
         needed = needfld("mixr") .or. needfld("zg") .or. needfld("rh") .or. &
                  needfld("pwc") .or. needfld("qbot") .or. &
                  needfld("vaveuq") .or. needfld("vaveuq")
      case ( "u", "v" )
         needed = needfld("u") .or. needfld("v") .or. needfld("vaveuq") .or. &
                  needfld("vavevq") .or. needfld("vaveut") .or.              &
                  needfld("vavevt") .or. needfld("vbot") .or.                &
                  needfld("ubot") .or. needfld("d10") .or.                   &
                  needfld("uas") .or. needfld("vas")
      case ( "qlg" )
         needed = needfld("qlg") .or. needfld("rh")
      case ( "qfg" )
         needed = needfld("qfg") .or. needfld("rh")
      case default
         print*, "Error - unsupported argument for need3dfld", name
         stop
      end select
   end function need3dfld

   subroutine inshallow(skip)
      use s2p_m 
      use history
      logical, intent(in) :: skip
      integer :: kdate, ktime, ieof

      read(un_in,iostat=ieof) ik, jk, kk, kdate, ktime, ktau, ndt,  &
                           rlong0, rlat0, schmidt

      if ( ieof /= 0 ) then
         ! This should never happen, EOF should be caught by getdate.
         print*, "Error, unexpected EOF in inshallow"
         stop
      end if

      if ( first_in ) then
         allocate ( u(ik,jk,kk), v(ik,jk,kk) )
      end if

      call readsave3 ( "zg", .not. skip )
      read(un_in) u
      read(un_in) v
      if ( .not. skip ) then
         call fix_winds(u, v)
         call vsavehist ( "u", u )
         call vsavehist ( "v", v )
      end if

      first_in = .false.

   end subroutine inshallow

   subroutine incom(skip)
      logical, intent(in) :: skip
      integer :: kdate, ktime, ieof
      integer :: ik, jk, kk

      read(un_in,iostat=ieof) ik, jk, kk, kdate, ktime, ktau, ndt,  &
                           rlong0, rlat0, schmidt
      if ( ieof /= 0 ) then
         ! This should never happen, EOF should be caught by getdate.
         print*, "Error, unexpected EOF in incom"
         stop
      end if

      read(un_in) sig

      call readsave2 ( "egm", .not. skip )
      call readsave2 ( "fgm", .not. skip ) 
      call readsave3 ( "dqc", .not. skip  )
      call readsave3 ( "dql", .not. skip  )
      call readsave3 ( "dqv", .not. skip  )
      call readsave3 ( "dtc", .not. skip  )
      call readsave3 ( "dtl", .not. skip  )
      call readsave3 ( "dtv", .not. skip  )
      call readsave3 ( "dtlw", .not. skip  )
      call readsave3 ( "dtsw", .not. skip  )
      call readsave3 ( "duv", .not. skip  )
      call readsave3 ( "dvv", .not. skip  )
      call readsave3 ( "omg", .not. skip  )

   end subroutine incom

   subroutine inrad(skip)
      logical, intent(in) :: skip
      integer :: kdate, ktime, ieof, ierr
      real, dimension(ik,jk) :: var

      read(un_in,iostat=ieof) ik, jk, kk, kdate, ktime, ktau, ndt,  &
                           rlong0, rlat0, schmidt
      if ( ieof /= 0 ) then
         ! This should never happen, EOF should be caught by getdate.
         print*, "Error, unexpected EOF in inrad"
         stop
      end if

      read(un_in) sig

      call readsave2 ( "sint", .not. skip  )
      call readsave2 ( "sot", .not. skip  ) 
      call readsave2 ( "soc", .not. skip  )
      call readsave2 ( "sgn", .not. skip  )
      call readsave2 ( "rtu", .not. skip  )
      call readsave2 ( "rtc", .not. skip  )
      call readsave2 ( "rgn", .not. skip  )
      call readsave2 ( "cld", .not. skip  )
      call readsave2 ( "cll", .not. skip  )
      call readsave2 ( "clm", .not. skip  )
      call readsave2 ( "clh", .not. skip  )
      call readsave2 ( "tsh", .not. skip  )
      call readsave2 ( "tsl", .not. skip  )
      ! Some files had qscr_ave as an extra field. Try reading this. If it
      ! fails then assume it's not present
      read ( unit=un_in, iostat=ierr ) var
      backspace ( unit=un_in )
      if ( ierr == 0 ) then
         call readsave2 ( "qscr_ave", .not. skip  )
      end if

   end subroutine inrad

   subroutine vread2(name, var, required, vread_err)

      ! Routine to read a variable from either a fortran binary or netcdf file.
      use mpidata_m
      character(len=*), intent(in) :: name
      real, dimension(:,:), intent(out) :: var
      logical, intent(in), optional :: required
      logical req
      integer, intent(out), optional :: vread_err
      integer v_err

      if ( present( required ) .and. present ( vread_err ) ) then
         req = required
      else
         req = .true.
      end if

      if ( myid == 0 ) then
         call paravar2a(name,var,nrec,req,v_err)
      else
         call paravar2b(name,nrec,req,v_err)
      end if
      
      if ( .not. req ) then
         vread_err = v_err
      else if ( present ( vread_err ) ) then
         ! If we get to this point, everything is ok.
         vread_err = 0
      end if

   end subroutine vread2
   
   subroutine vread3(name, var)

      ! Routine to read a variable from either a fortran binary or netcdf file. 
      use mpidata_m
      character(len=*), intent(in) :: name
      real, dimension(:,:,:), intent(out) :: var
      integer :: pkl

      pkl = size(var,3)
      if ( myid == 0 ) then
         call paravar3a(name,var,nrec,pkl)
      else
         call paravar3b(name,nrec,pkl)
      end if

   end subroutine vread3
   
   subroutine readsave2 ( name, save_flag, input_name, required )
!     Read an array from the input file and call savehist
      use history, only : savehist, needfld
      use s2p_m
      character(len=*), intent(in) :: name
      logical, intent(in), optional :: save_flag, required
      character(len=*), intent(in), optional :: input_name
      real, dimension(ik,jk) :: array
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
      real, dimension(ik,jk,kk) :: array
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
                           kdate, ktime, ntracers, ksoil, kice, debug,        &
                           nqg )

      use mpidata_m
      use netcdf
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

      real, intent(inout)  :: hres
      real, intent(inout)  :: minlon, maxlon, dlon, minlat, maxlat, dlat
      integer, intent(out) :: kdate, ktime, ntracers, ksoil, kice
      logical, intent(in)  :: debug
      integer, intent(out) :: nqg

      integer, parameter :: diag = 0, id = 1, jd = 1
      real(kind=rx), parameter :: schm13 = 0.1

      integer, parameter :: ntang=2, idiag=0
      integer :: i, j
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
      real :: grlong, grlat, grlongu, grlatu, grlongv, grlatv, tmp

!     Read the header here because doing the CC grid initialisation before
!     alloc_indata minimises the total memory requirements

      nrec = 1
      ! Get the total number of timesteps
      ierr = nf90_inq_dimid ( ncid, "time", dimid )
      call check_ncerr(ierr, "Error getting time dimension")
      ierr = nf90_inquire_dimension ( ncid, dimid, len=maxrec )
      call check_ncerr(ierr,"Error getting number of sets")
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
      ik = il  ! These are only set once
      jk = jl
      kk = kl
!     kwt = int_header(19)  ! No longer used
      nsd = int_header(5)
      ndt = int_header(14)
      ms = int_header(34)
      nqg = int_header(23)
      ilt = int_header(42)
      ntrac = int_header(43)

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
      ! Need date for setting time origin. Call getdate with nrec=1
      nrec = 1
      call getdate(kdate, ktime, ieof)
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
      if ( nqg >= 8 ) then
         ksoil = ms
      else
         ksoil = 2
      end if
      if ( nqg >= 11 ) then
         kice = ms
      else
         kice = 0
      end if


      if ( ilt > 1 ) then
         ntracers = ntrac
      else
         ntracers = 0
      end if

!     Have to read sig here because it's required for openhist.
      allocate ( sig(kl), dsig(kl), zsoil(ksoil) )

      if ( kl > 1 ) then
            ! Get sigma levels from level variable
            ierr = nf90_inq_varid (ncid, "lev", vid )
            call check_ncerr(ierr, "Error getting vid for lev")
            ierr = nf90_get_var ( ncid, vid, sig)
            call check_ncerr(ierr, "Error getting levels")
      else
         sig = 1.0
      end if
      call sig2ds(sig, dsig)
      ! Note that some initial condition files don't have zsoil
      if ( cf_compliant ) then
         ierr = nf90_inq_varid (ncid, "zsoil", vid )
         call check_ncerr(ierr, "Error getting vid for zsoil")
         ierr = nf90_get_var ( ncid, vid, zsoil)
         call check_ncerr(ierr, "Error getting zsoil")
      end if

!     Set all the resolution parameters
      npanels = jl/il - 1
      ifull=il*jl
      ij=il*jl
      ijk=il*jl*kl
      iquad=1+il*((8*npanels)/(npanels+4))

      call setxyz ( il, jl, kl, npanels, ifull, iquad, idiag, id, jd,        &
                rlong0, rlat0, schmidt, schm13, ntang, erad )

      if ( int_default == int_none ) then
         nxhis = il
         nyhis = jl
      else
!        If hres has been set use that, otherwise calculate a value
!        appropriate to the resolution. 360/(4*il) is the approx resolution
!        along the equator for the cube and 360/(6*il) for the octagon
         if ( hres == 0.0 ) then
            if ( npanels == 5 ) then
               hres = ceiling ( 90.0/il )
            else
               hres = ceiling ( 60.0/il )
            end if
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

!     To save memory de-allocate a number of arrays defined by setxyz
!     that aren't needed by cc2hist.
      deallocate ( f, fu, fv, dmdx, dmdy, dmdxv, dmdyu )

      if ( myid == 0 ) then

         allocate ( nface(nxhis,nyhis) )
         allocate ( xg(nxhis,nyhis), yg(nxhis,nyhis) )
         allocate ( hlon(nxhis), hlat(nyhis) )
         if ( int_default == int_none ) then
            hlat = (/ ( real(j), j=1,nyhis ) /)
            hlon = (/ ( real(i), i=1,nxhis ) /)
         else
            !     Set lats and longs
            do j=1,nyhis
               hlat(j) = minlat + (j-1)*(maxlat-minlat)/(nyhis-1)
            end do
            if ( maxlon - minlon == 360.0 ) then
               do i=1,nxhis
                  hlon(i) = minlon + (i-1)*(maxlon-minlon)/nxhis
               end do
            else
               do i=1,nxhis
                  hlon(i) = minlon + (i-1)*(maxlon-minlon)/(nxhis-1)
               end do
            end if
            do j=1,nyhis
               do i=1,nxhis
                  call latltoij ( hlon(i), hlat(j), xg(i,j), yg(i,j), nface(i,j),&
                                  rlong0, rlat0, schmidt, schm13 )
               enddo
            enddo

         end if

      end if   

!     These aren't needed now.
      deallocate ( xx4, yy4 )

   end subroutine initialise

   subroutine final_init(varlist, nvars)
      ! Some final initialisation that requires needfld, so can only be done
      ! after openhist has been called.
      use mpidata_m
      use newmpar_m
      use history, only : needfld
      use xyzinfo_m
      use indices_m
      use parm_m, only : rlong0, rlat0
      use physparams, only : pi
      type(input_var), dimension(:) :: varlist
      integer, intent(in) :: nvars
      real :: sinlong, coslong, sinlat, coslat
      real :: polenx, poleny, polenz, zonx, zony, zonz, den
      real :: theta_lc, rlongdeg, rlatdeg, ri, rj
      integer :: i, j, iq, ivar
      logical :: need_rotate

      deallocate ( em )
      deallocate ( i_wu, i_sv, i_eu, i_nv )

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
           needfld("taux") .or. needfld("tauy") .or.          &
           needfld("u10max") .or. needfld("v10max") .or.      &
           needfld("uas") .or.  needfld("vas") .or.           &
           needfld("d10") ) then

         allocate ( costh(il,jl), sinth(il,jl) )

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
!              Set up unit zonal vector components
               zonx = poleny*z(iq)-polenz*y(iq)
               zony = polenz*x(iq)-polenx*z(iq)
               zonz = polenx*y(iq)-poleny*x(iq)
!              Allow for poles by taking max
               den = sqrt( max(zonx**2 + zony**2 + zonz**2,1.e-7) )
               costh(i,j) =  (zonx*ax(iq)+zony*ay(iq)+zonz*az(iq))/den
               sinth(i,j) = -(zonx*bx(iq)+zony*by(iq)+zonz*bz(iq))/den
            end do
         end do

      end if
!     x, y, z, ax, ay, az, bx, by, bz no longer needed.
!     ax etc are pointers to setxyz private arrays so the space can't be freed.
      deallocate ( x, y, z ) 

   end subroutine final_init
   
   subroutine fix_winds2 ( u, v )
!     Convert winds from grid directions to standard zonal and meridional.
      use staguv_m
      use newmpar_m
      use xyzinfo_m
      real, dimension(:,:), intent(inout) :: u, v
      real, dimension(il,jl) :: uzon, vmer

      uzon = costh*u - sinth*v
      vmer = sinth*u + costh*v
!        Now save these back to the original arrays.
         u = uzon
         v = vmer
   end subroutine fix_winds2

   subroutine fix_winds3 ( u, v )
!     Convert winds from grid directions to standard zonal and meridional.
      use staguv_m
      use newmpar_m
      use xyzinfo_m
      real, dimension(:,:,:), intent(inout) :: u, v
      real, dimension(il,jl) :: uzon, vmer
      integer :: k

      do k=1,kl
         uzon(:,:) = costh*u(:,:,k) - sinth*v(:,:,k)
         vmer(:,:) = sinth*u(:,:,k) + costh*v(:,:,k)
!           Now save these back to the original arrays.
            u(:,:,k) = uzon
            v(:,:,k) = vmer
      enddo
   end subroutine fix_winds3

   subroutine get_var_list(varlist, nvars)
      ! Get a list of the variables in the input file
      use history, only : addfld, int_default, cf_compliant
      use interp_m, only : int_nearest, int_none
      use physparams, only : grav, rdry
      use s2p_m, only: use_plevs, plevs
      type(input_var), dimension(:), pointer :: varlist
      integer, intent(out) :: nvars
      integer :: ierr, ndimensions, nvariables, ndims, ivar, int_type, xtype
      integer :: londim, latdim, levdim, timedim, vid, ihr, ind
      integer, dimension(nf90_max_var_dims) :: dimids
      character(len=10) :: vname, substr
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
      ierr = nf90_inq_dimid(ncid, "time", timedim)
      call check_ncerr(ierr,"Error getting timeid")

      nvars = 0
      ! For the sigma to pressure initialisation to work properly
      ! the surface pressure has to be read before any of the 3D variables
      ! Force it to be first by having an extra 0 iteration on the loop
      ! Some special purpose outputs don't have psf, so only treat cases where
      ! use_plevs=T in this way.
      if ( use_plevs ) then
         ivar_start = 0
      else
         ivar_start = 1
      end if
      do ivar=ivar_start,nvariables
         if ( ivar == 0 ) then
            ierr = nf90_inq_varid ( ncid, "psf", vid )
            call check_ncerr(ierr, "Error getting vid for psf")
            ierr = nf90_inquire_variable (ncid, vid, name=vname, ndims=ndims, dimids=dimids, xtype=xtype)
            call check_ncerr(ierr, "nf90_inquire_variable error")
         else
            ierr = nf90_inquire_variable (ncid, ivar, name=vname, ndims=ndims, dimids=dimids, xtype=xtype)
            call check_ncerr(ierr, "nf90_inquire_variable error")
            if ( use_plevs .and. vname == "psf" ) then
               ! This has already been handled above
               cycle
            end if
         end if
!         print*, ivar, vname, ndims, dimids(1:ndims)
         if ( ndims == 4 ) then
            ! Should be lon, lat, lev, time
            if ( match ( dimids(1:ndims), (/ londim, latdim, levdim, timedim /) ) ) then
               nvars = nvars + 1
               varlist(nvars)%fixed = .false.
               varlist(nvars)%ndims = 3  ! Space only
            else
               print*, "Error, unexpected dimensions in input variable", vname
               stop
            end if
         else if ( ndims == 3 ) then

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
         nvars = nvars+1
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
         xmin = varlist(ivar)%add_offset + varlist(ivar)%scale_factor*vmin
         xmax = varlist(ivar)%add_offset + varlist(ivar)%scale_factor*vmax
         ! As a check re-calc offset, scalef
         sf = (xmax - xmin) /(real(vmax)-real(vmin)) ! jlm fix for precision problems
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
              varlist(ivar)%vname == "maxrnd") then
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
               (/ "vegt    ", "soilt   ", "rsmin   ", "rs      ", "zolnd   ", &
                  "sigmf   ", "wetfrac ", "wetfrac?", "tgg?    ", "tgg??   ", &
                  "sal??   ", "roadtgg?", "rooftgg?", "waletgg?", "walwtgg?"  &
                  /)) .and. int_type /= int_none ) then
            int_type = int_nearest
         else
            int_type = int_default
         end if
         if ( varlist(ivar)%vname == "pmsl" ) then
            varlist(ivar)%vname = "psl"
         else if ( varlist(ivar)%vname == "zht" ) then
            varlist(ivar)%vname = "zs"
            varlist(ivar)%units = "m"
            varlist(ivar)%long_name = "Surface height"
            xmin = 0.
            xmax = 9000.
         else if ( varlist(ivar)%vname == "psf" ) then
            cycle  ! Skip this one to avoid messages about it never being set
         end if
         call cc_cfproperties(varlist(ivar), std_name, cell_methods)
         if ( varlist(ivar)%fixed ) then
            call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name, &
                      varlist(ivar)%units, xmin, xmax, 1, ave_type="fixed", int_type=int_type, std_name=std_name )
         else if ( varlist(ivar)%ndims == 2 ) then
            if ( varlist(ivar)%vname(1:3) == "max" ) then
               ! Special check for maximum rainfall rate
               call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name, &
                      varlist(ivar)%units, xmin, xmax, 1, ave_type="max", std_name=std_name, cell_methods=cell_methods )
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
                      int_type=int_type )
            end if
         else if ( varlist(ivar)%ndims == 3 ) then
            call addfld ( varlist(ivar)%vname, varlist(ivar)%long_name, &
                      varlist(ivar)%units, xmin, xmax, nlev, multilev=.true., &
                      std_name=std_name, cell_methods=cell_methods, int_type=int_type )
         end if
      end do

      ! Extra fields are handled explicitly
      call addfld ( "ps", "Surface pressure", "hPa", 0., 1200., 1, &
                     std_name="surface_air_pressure" )
      call addfld ( "pwc", "Precipitable water column", "kg/m^2", 0.0, 100.0, 1, std_name="atmosphere_water_vapor_content")
      call addfld ( "tsea", "Sea surface temperature", "K", 150., 350., 1, &
                     std_name="sea_surface_temperature" )
      call addfld ( "land_mask", "Land-sea mask", "",  0.0, 1.0, 1, &
                     ave_type="fixed", int_type=int_nearest )
      call addfld ( "d10", "10m wind direction", "deg", 0.0, 360.0, 1 )
      call addfld ( "uas", "x-component 10m wind", "m/s", -100.0, 100.0, 1 )
      call addfld ( "vas", "y-component 10m wind", "m/s", -100.0, 100.0, 1 )
      ! Packing is not going to work well in this case
      ! For height, estimate the height of the top level and use that
      ! for scaling
      if ( use_plevs ) then
         topsig = plevs(nlev) / 1000.
      else
         topsig = sig(kk)
      end if
      ! Assume isothermal 300K profile. This will lead to an overestimate
      ! of the height so is safe. Using moist adiabatic lapse rate can
      ! give an underestimate.
      topheight = -rdry*300./grav * log(topsig)
      ! Round to next highest 1000 m
      topheight = 1000*(floor(0.001*topheight)+1)
      call addfld ( "zg", "Geopotential height", "m", 0., topheight, nlev, &
                     multilev=.true., std_name="geopotential_height" )

      ! If the output uses pressure levels save the lowest sigma level of
      ! the basic fields.
      call addfld ( "tbot", "Air temperature at lowest sigma level", "K", 100., 400., 1 )
      call addfld ( "ubot", "Eastward wind at lowest sigma level", "m/s", -100., 100., 1)
      call addfld ( "vbot", "Northward wind at lowest sigma level", "m/s", -100., 100., 1)
      call addfld ( "qbot", "Specific humidity at lowest sigma level", "kg/kg", 0., 0.1, 1)

      ! Vertical averaged fluxes. Note that these are NOT vertical integrals
      call addfld ( "vaveuq", "Vertical average of zonal humidity flux", "m/s kg/kg", -0.5, 0.5, 1, std_name=std_name )
      call addfld ( "vavevq", "Vertical average of meridional humidity flux", "m/s kg/kg", -0.5, 0.5, 1, std_name=std_name )
      call addfld ( "vaveut", "Vertical average of zonal temperature flux", "m/s K", -1e4, 1e4, 1, std_name=std_name )
      call addfld ( "vavevt", "Vertical average of meridional temperature flux", "m/s K", -1e4, 2e4, 1, std_name=std_name )
      call addfld ( "rh", "Relative humidity", "%", 0., 110., nlev,  multilev=.true., std_name="relative_humidity" )

      if ( cf_compliant ) then
         ! Define as an extra field for now
         call addfld('tgg','Soil temperature','K',100.,350.,ksoil,soil=.true.)
         ! Should have std_name = volume_fraction_of_water_in_soil, units=1
         ! Mentioned in Gregory email 2005-12-01. In official list?
         !call addfld('wb','Soil moisture','frac',0.,1.,ksoil,soil=.true.)
         call addfld('wetfrac','Soil moisture','frac',0.,1.,ksoil,soil=.true.)
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
      integer :: i, j, len, start
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
               j = j + index(strarray(i)(j+1:len),"?")
               if ( j == 0 ) then
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
      real, dimension(ik,jk,kk), intent(in) :: t, q, ql, qf
      real, dimension(ik,jk), intent(in) :: psl
      real, dimension(kk), intent(in) :: sig
      real, dimension(ik,jk,kk), intent(out) :: rh
      real, dimension(ik,jk) :: p, tliq, fice, qsw, qc
      real, dimension(ik,jk) :: qsi, qsl, deles
      real, parameter :: tfrz  = 273.1
      real, parameter :: tice  = 233.15
      real, parameter :: cp    = 1004.64
      real, parameter :: hl    = 2.5104e6
      real, parameter :: hlf   = 3.36e5
      real, parameter :: epsil = 0.622
      integer :: k

      do k=1,kk
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
      do k=1,ksoil
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
      do k=1,ksoil
         write(tmpname,'(a,i1)') 'wetfrac', k
         if ( vname == tmpname ) then
            is_soil_var = .true.
            return
         end if
      end do
   end function is_soil_var

   ! MJT notes - this could be modified for parallel files

   ! From ccam infile.f
   subroutine fill_cc(a_io,value)
!     routine fills in interior of an array which has undefined points
      use newmpar_m
      use indices_m
      real, intent(inout) :: a_io(il,jl)         ! input and output array
      real, intent(in)    :: value            ! array value denoting undefined
      real b(ifull), a(ifull)
      integer :: num, nrem, iq, neighb, ierr, i, j
      real :: av, avx
      
      ! Really just a reshape
      do j=1,jl
         do i=1,il
            iq = i + (j-1)*il
            a(iq) = a_io(i,j)
         end do
      end do
      num=0
      nrem = 1    ! Just for first iteration
!     nrem_gmin used to avoid infinite loops, e.g. for no sice
      do while ( nrem > 0 )
         nrem=0
         num=num+1
         do iq=1,ifull
            b(iq)=a(iq)
            if(a(iq)==value)then
               neighb=0
               av=0.
               if(a(i_n(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(i_n(iq))
               endif
               if(a(i_e(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(i_e(iq))
               endif
               if(a(i_w(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(i_w(iq))
               endif
               if(a(i_s(iq)).ne.value)then
                  neighb=neighb+1
                  av=av+a(i_s(iq))
               endif
               if(neighb>0)then
                  b(iq)=av/neighb
                  avx=av
               else
                  nrem=nrem+1   ! current number of points without a neighbour
               endif
            endif
         end do
         a(:)=b(:)
!         if(nrem>0)then
!            print*, "FILL", num, nrem, avx
!         endif                  ! (nrem>0)
         if ( nrem == ifull ) then
            print*, "Error in fill_cc - no points defined"
            return
         end if
      end do
      do j=1,jl
         do i=1,il
            iq = i + (j-1)*il
            a_io(i,j) = a(iq)
         end do
      end do
   end subroutine fill_cc
   
   subroutine paraopen(ifile,nmode,ncid)
  
      use mpidata_m
  
      include 'mpif.h'
  
      integer, intent(in) :: nmode
      integer, intent(out) :: ncid
      integer ier, ip, n, rip, ierr
      integer pil_g, pjl_g
      integer, dimension(3) :: jdum
      integer, dimension(54) :: int_header
      character(len=*), intent(in) :: ifile
      character(len=266) :: pfile
      character(len=7) :: sdecomp

      if (myid==0) then      
  
      ! parallel file input
      ip = 0
      write(pfile,"(a,'.',i6.6)") trim(ifile), ip
      ierr = nf90_open(pfile, nmode, ncid)
      call check_ncerr(ierr, "Error opening file")
      
      write(6,*) "Using parallel input files"
      
      ! parallel metadata
      ier = nf90_get_att(ncid, nf90_global, "nproc", pnproc)
      call check_ncerr(ier, "nproc")
      
      allocate(ioff(0:pnproc-1,0:5),joff(0:pnproc-1,0:5))
      
      end if
      
      call MPI_Bcast(pnproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      lproc=pnproc/nproc !number of files each mpi_proc will work on      
      allocate(ncid_in(0:lproc-1))
      
      if ( myid /= 0 ) then
      
         do ip = 0,lproc-1
            rip = myid*lproc + ip
            write(pfile,"(a,'.',i6.6)") trim(ifile), rip
            ier = nf90_open ( pfile, nmode, ncid_in(ip) )
            if (ier /= 0 ) then
               write(6,*) "ERROR: Cannot open ",trim(pfile)
               call check_ncerr(ier, "open")
            end if
         end do
         ncid = ncid_in(0) 
      
      else
      
         ncid_in(0) = ncid
         do ip = 1,lproc-1
            rip = ip
            write(pfile,"(a,'.',i6.6)") trim(ifile), rip
            ier = nf90_open ( pfile, nmode, ncid_in(ip) )
            if (ier /= 0 ) then
               write(6,*) "ERROR: Cannot open ",trim(pfile)
               call check_ncerr(ier, "open")
            end if
         end do
      
         !  Get dimensions from int_header
         ier = nf90_get_att(ncid_in(0), nf90_global, "int_header", int_header)
         call check_ncerr(ier, "int_header")
         ! Only a few values are used
         pil_g = int_header(1)
         pjl_g = int_header(2)

         !  Calculate il, jl from these
         sdecomp = ''
         ier = nf90_get_att(ncid_in(0), nf90_global, "decomp", sdecomp)
         call check_ncerr(ier, "decomp")
         select case(sdecomp)
            case ("uniform")
               do n = 0,5
                  call proc_setup_uniform(pil_g, pjl_g, pnproc, n, pil, pjl, ioff(:,n), &
	                                  joff(:,n), pnpan)
               end do
            case ("face")
               call proc_setup(pil_g, pjl_g, pnproc, 0, pil, pjl, ioff(:,0), joff(:,0), pnpan)
               do n = 1,5
                  ioff(:,n) = ioff(:,0)
                  joff(:,n) = joff(:,0)
               end do
            case default
               write(6,*) "ERROR: Unknown decomposition ",trim(sdecomp)
               stop
         end select
      
         jdum(1) = pil
         jdum(2) = pjl
         jdum(3) = pnpan
      
      end if
      
      call MPI_Bcast(jdum(1:3),3,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      pil   = jdum(1)
      pjl   = jdum(2)
      pnpan = jdum(3)
      
   end subroutine paraopen
   
   subroutine paraclose
      integer ip, ierr
      
      do ip=0,lproc-1
         ierr = nf90_close(ncid_in(ip))
      end do
      deallocate(ncid_in)
   
   end subroutine paraclose
   
   subroutine paravar2a(name,var,nrec,required,vread_err)
      use mpidata_m
      include 'mpif.h'
      integer, intent(in) :: nrec
      integer, intent(out) :: vread_err
      integer ip, n, pil_g, pjl_g, vid, ierr, vartyp
      real, dimension(:,:), intent(out) :: var
      real, dimension(pil,pjl*pnpan,0:lproc-1) :: inarray2
      real, dimension(pil,pjl*pnpan,pnproc) :: inarray2gather
      real addoff, sf
      logical, intent(in) :: required
      character(len=*), intent(in) :: name

      pil_g = size(var,1)
      pjl_g = size(var,2)
   
      ierr = nf90_inq_varid (ncid_in(0), name, vid )
      ! If the variable has the required flag set to false, and the 
      ! vread_err argument is present, then return an error flag rather
      ! then abort if the variable isn't found.
      if ( .not. required .and. ierr /= NF90_NOERR ) then
         vread_err = ierr
         return
      end if
      call check_ncerr(ierr, "Error getting vid for "//name)
      
      do ip = 0,lproc-1
         ierr = nf90_get_var ( ncid_in(ip), vid, inarray2(:,:,ip), start=(/ 1, 1, nrec /), count=(/ pil, pjl*pnpan, 1 /) )
         call check_ncerr(ierr, "Error getting var "//name)
      end do

!     Check the type of the variable
      ierr = nf90_inquire_variable ( ncid_in(0), vid, xtype=vartyp)
      if ( vartyp == NF90_SHORT ) then
         ierr = nf90_get_att ( ncid_in(0), vid, "add_offset", addoff )
         call check_ncerr(ierr, "Error getting add_offset attribute")
         ierr = nf90_get_att ( ncid_in(0), vid, "scale_factor", sf )
         call check_ncerr (ierr,"Error getting scale_factor attribute")
         if ( all( inarray2 == -32501. ) ) then
            inarray2 = NF90_FILL_FLOAT
         else
            inarray2 = addoff + inarray2*sf
         end if
      end if
      
      call MPI_GATHER(inarray2,pil*pjl*pnpan*lproc,MPI_REAL,inarray2gather,pil*pjl*pnpan*lproc,MPI_REAL,0,MPI_COMM_WORLD,ierr)
         
      do ip = 0,pnproc-1   
         do n = 0,pnpan-1
            var(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g) = &
               inarray2gather(1:pil,1+n*pjl:(n+1)*pjl,ip+1)
         end do
      end do
   
   end subroutine paravar2a

   subroutine paravar2b(name,nrec,required,vread_err)
      use mpidata_m
      include 'mpif.h'
      integer, intent(in) :: nrec
      integer, intent(out) :: vread_err
      integer ip, ierr, vid, vartyp
      real, dimension(pil,pjl*pnpan,0:lproc-1) :: inarray2
      real, dimension(0,0,0) :: inarray2gather
      real addoff, sf
      logical, intent(in) :: required
      character(len=*), intent(in) :: name
  
      ierr = nf90_inq_varid (ncid_in(0), name, vid )
      ! If the variable has the required flag set to false, and the 
      ! vread_err argument is present, then return an error flag rather
      ! then abort if the variable isn't found.
      if ( .not. required .and. ierr /= NF90_NOERR ) then
         vread_err = ierr
         return
      end if
      call check_ncerr(ierr, "Error getting vid for "//name)
      
      do ip = 0,lproc-1
         ierr = nf90_get_var ( ncid_in(ip), vid, inarray2(:,:,ip), start=(/ 1, 1, nrec /), count=(/ pil, pjl*pnpan, 1 /) )
         call check_ncerr(ierr, "Error getting var "//name)
      end do

!     Check the type of the variable
      ierr = nf90_inquire_variable ( ncid_in(0), vid, xtype=vartyp)
      if ( vartyp == NF90_SHORT ) then
         ierr = nf90_get_att ( ncid_in(0), vid, "add_offset", addoff )
         call check_ncerr(ierr, "Error getting add_offset attribute")
         ierr = nf90_get_att ( ncid_in(0), vid, "scale_factor", sf )
         call check_ncerr (ierr,"Error getting scale_factor attribute")
         if ( all( inarray2 == -32501. ) ) then
            inarray2 = NF90_FILL_FLOAT
         else
            inarray2 = addoff + inarray2*sf
         end if
      end if
      
      call MPI_GATHER(inarray2,pil*pjl*pnpan*lproc,MPI_REAL,inarray2gather,pil*pjl*pnpan*lproc,MPI_REAL,0,MPI_COMM_WORLD,ierr)
   
   end subroutine paravar2b
   
   subroutine paravar3a(name,var,nrec,pkl)
      use mpidata_m
      include 'mpif.h'
      integer, intent(in) :: nrec, pkl
      integer ip, n, pil_g, pjl_g, vid, ierr, vartyp
      real, dimension(:,:,:), intent(out) :: var
      real, dimension(pil,pjl*pnpan,pkl,0:lproc-1) :: inarray3
      real, dimension(pil,pjl*pnpan,pkl,pnproc) :: inarray3gather
      real addoff, sf
      character(len=*), intent(in) :: name

      pil_g = size(var,1)
      pjl_g = size(var,2)

      ierr = nf90_inq_varid (ncid_in(0), name, vid )
      call check_ncerr(ierr, "Error getting vid for "//name)
      
      do ip = 0,lproc-1
         ierr = nf90_get_var ( ncid_in(ip), vid, inarray3(:,:,:,ip), start=(/ 1, 1, 1, nrec /), count=(/ pil, pjl*pnpan, pkl, 1 /) )
         call check_ncerr(ierr, "Error getting var "//name)
      end do
      
      ierr = nf90_inquire_variable ( ncid_in(0), vid, xtype=vartyp)
      if ( vartyp == NF90_SHORT ) then
         ierr = nf90_get_att ( ncid_in(0), vid, "add_offset", addoff )
         call check_ncerr(ierr, "Error getting add_offset attribute")
         ierr = nf90_get_att ( ncid_in(0), vid, "scale_factor", sf )
         call check_ncerr (ierr,"Error getting scale_factor attribute")
         if ( all( inarray3 == -32501. ) ) then
           inarray3 = NF90_FILL_FLOAT
         else
           inarray3 = addoff + inarray3*sf
         end if
      end if
      
      call MPI_GATHER(inarray3,pil*pjl*pnpan*pkl*lproc,MPI_REAL,inarray3gather,pil*pjl*pkl*pnpan*lproc,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      do ip = 0,pnproc-1      
         do n = 0,pnpan-1
            var(1+ioff(ip,n):pil+ioff(ip,n),1+joff(ip,n)+n*pil_g:pjl+joff(ip,n)+n*pil_g,:) = &
               inarray3gather(1:pil,1+n*pjl:(n+1)*pjl,:,ip+1)
         end do
      end do
   
   end subroutine paravar3a

   subroutine paravar3b(name,nrec,pkl)
      use mpidata_m
      include 'mpif.h'
      integer, intent(in) :: nrec, pkl
      integer ip, ierr, vid, vartyp
      real, dimension(pil,pjl*pnpan,pkl,0:lproc-1) :: inarray3
      real, dimension(0,0,0,0) :: inarray3gather
      real addoff, sf
      character(len=*), intent(in) :: name

      ierr = nf90_inq_varid (ncid_in(0), name, vid )
      call check_ncerr(ierr, "Error getting vid for "//name)
      
      do ip = 0,lproc-1
         ierr = nf90_get_var ( ncid_in(ip), vid, inarray3(:,:,:,ip), start=(/ 1, 1, 1, nrec /), count=(/ pil, pjl*pnpan, pkl, 1 /) )
         call check_ncerr(ierr, "Error getting var "//name)
      end do

      ierr = nf90_inquire_variable ( ncid_in(0), vid, xtype=vartyp)
      if ( vartyp == NF90_SHORT ) then
         ierr = nf90_get_att ( ncid_in(0), vid, "add_offset", addoff )
         call check_ncerr(ierr, "Error getting add_offset attribute")
         ierr = nf90_get_att ( ncid_in(0), vid, "scale_factor", sf )
         call check_ncerr (ierr,"Error getting scale_factor attribute")
         if ( all( inarray3 == -32501. ) ) then
           inarray3 = NF90_FILL_FLOAT
         else
           inarray3 = addoff + inarray3*sf
         end if
      end if
      
      call MPI_GATHER(inarray3,pil*pjl*pnpan*pkl*lproc,MPI_REAL,inarray3gather,pil*pjl*pkl*pnpan*lproc,MPI_REAL,0,MPI_COMM_WORLD,ierr)

   end subroutine paravar3b

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
         nxproc=1
         nyproc=nproc
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
      integer, parameter :: npanels=5
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
      do n=0,nproc-1
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

end module work
