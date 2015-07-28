! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This example writes some surface pressure and temperatures. It is
! intended to illustrate the use of the netCDF fortran 90 API. The
! companion program sfc_pres_temp_rd.f shows how to read the netCDF
! data file created by this program.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: sfc_pres_temp_wr.f90,v 1.9 2007/01/24 19:32:10 russ Exp $

! Trivially modified by David.Benn@csiro.au to run under MPI and
! with conditional netcdf module to suit test purposes, July 2015.
! COPYRIGHT file taken from:
!  http://www.unidata.ucar.edu/software/netcdf/docs/copyright.html

program sfc_pres_temp_wr
#ifdef PARNETCDF
  use pnetcdf_m
#else
  use netcdf_m
#endif

  use mpi

  implicit none

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "sfc_pres_temp.nc"
  integer :: ncid

  ! We are writing 2D data, a 6 x 12 lat-lon grid. We will need two
  ! netCDF dimensions.
  integer, parameter :: NDIMS = 2
  integer, parameter :: NLATS = 6, NLONS = 12, SUBSET_NLONS = 3
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  integer :: lat_dimid, lon_dimid

  ! In addition to the latitude and longitude dimensions, we will also
  ! create latitude and longitude netCDF variables which will hold the
  ! actual latitudes and longitudes. Since they hold data about the
  ! coordinate system, the netCDF term for these is: "coordinate
  ! variables."
  real :: lats(NLATS), lons(NLONS)
  integer :: lat_varid, lon_varid
  real, parameter :: START_LAT = 25.0, START_LON = -125.0

  ! We will write surface temperature and pressure fields.
  character (len = *), parameter :: PRES_NAME="pressure"
  character (len = *), parameter :: TEMP_NAME="temperature"
  integer :: pres_varid, temp_varid
  integer :: dimids(NDIMS)
  integer :: start(1)
  integer :: count(1)

  ! It's good practice for each variable to carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PRES_UNITS = "hPa"
  character (len = *), parameter :: TEMP_UNITS = "celsius"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  ! We will create some pressure and temperature data to write out.
  real :: pres_out(NLONS, NLATS), temp_out(NLONS, NLATS)
  real, parameter :: SAMPLE_PRESSURE = 900.0
  real, parameter :: SAMPLE_TEMP = 9.0

  ! Loop indices
  integer :: lat, lon

  ! MPI variables
  integer :: ierror, my_rank, num_procs

  ! start up MPI
  call MPI_Init(ierror)

  ! find out process rank
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)

  ! find out number of processes
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierror)

  if (my_rank .eq. 0) then

      ! Create pretend data. If this wasn't an example program, we would
      ! have some real data to write, for example, model output.
      do lat = 1, NLATS
         lats(lat) = START_LAT + (lat - 1) * 5.0
      end do
      do lon = 1, NLONS
         lons(lon) = START_LON + (lon - 1) * 5.0
      end do
      do lon = 1, NLONS
         do lat = 1, NLATS
            pres_out(lon, lat) = SAMPLE_PRESSURE + (lon - 1) * NLATS + (lat - 1)
            temp_out(lon, lat) = SAMPLE_TEMP + .25 * ((lon - 1) * NLATS + (lat - 1))
         end do
      end do

      ! Create the file.
      call check( ncf90_create(FILE_NAME, ncf90_clobber, ncid) )
      print *,">> created"

      ! Define the dimensions.
      call check( ncf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
      call check( ncf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
      print *,">> created dims"

      ! Define the coordinate variables. They will hold the coordinate
      ! information, that is, the latitudes and longitudes. A varid is
      ! returned for each.
      call check( ncf90_def_var(ncid, LAT_NAME, NCF90_REAL, lat_dimid, lat_varid) )
      call check( ncf90_def_var(ncid, LON_NAME, NCF90_REAL, lon_dimid, lon_varid) )
      print *,">> created vars 1"

      ! Assign units attributes to coordinate var data. This attaches a
      ! text attribute to each of the coordinate variables, containing the
      ! units.
      call check( ncf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
      call check( ncf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
      print *,">> added attrs 1"

      ! Define the netCDF variables. The dimids array is used to pass the
      ! dimids of the dimensions of the netCDF variables.
      dimids = (/ lon_dimid, lat_dimid /)
      call check( ncf90_def_var(ncid, PRES_NAME, NCF90_REAL, dimids, pres_varid) )
      call check( ncf90_def_var(ncid, TEMP_NAME, NCF90_REAL, dimids, temp_varid) )
      print *,">> created vars 2"

      ! Assign units attributes to the pressure and temperature netCDF
      ! variables.
      call check( ncf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
      call check( ncf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) )
      print *,">> added attrs 2"

      ! End define mode.
      call check( ncf90_enddef(ncid) )
      print *,">> ended define mode"

      ! Write the coordinate variable data. This will put the latitudes
      ! and longitudes of our data grid into the netCDF file.
      start = (/ 1 /)
      count = (/ SUBSET_NLONS /)
      call check( ncf90_put_var(ncid, lat_varid, lats, start) )
      call check( ncf90_put_var(ncid, lon_varid, lons, start, count) )
      print *,">> added vars 1"

      ! Write the pretend data. This will write our surface pressure and
      ! surface temperature data. The arrays of data are the same size as
      ! the netCDF variables we have defined.
      call check( ncf90_put_var(ncid, pres_varid, pres_out) )
      call check( ncf90_put_var(ncid, temp_varid, temp_out) )
      print *,">> added vars 2"

      ! Close the file.
      call check( ncf90_close(ncid) )
      print *,">> closed file"

      ! If we got this far, everything worked as expected. Yipee!
      print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"

  end if

  ! shut down MPI
  call MPI_Finalize(ierror)

contains
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= ncf90_noerr) then
      print *, trim(ncf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program sfc_pres_temp_wr
