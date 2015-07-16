module pnetcdf_m

!   This module implements parallel netcdf functions, exposing
!   an interface that allows replacement with a different implementation.
!
!   Notes:
!   - Each function returns a status value.
!   - We can either have one module with N different compile-time
!     conditional blocks and less variation in Makefile (libs, incs)
!     or we can have N modules with more variation in Makefile (libs,
!     incs, dependencies). If N is small, the first may be better. On
!     the other hand, the second allows us to vary the module
!     implementation wildly if necessary. There is a fair amount of
!     minor (e.g. naming) variation, so the second option seems best.
!
!   References:
!   - http://trac.mcs.anl.gov/projects/parallel-netcdf
!   - http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

    use mpi
    use pnetcdf

    implicit none

    private

!   NetCDF functions
    public :: ncf90_close, ncf90_copy_att, ncf90_create, &
              ncf90_def_dim, ncf90_def_var, ncf90_enddef, &
              ncf90_get_att, ncf90_inq_attname, ncf90_inq_dimid, &
              ncf90_inq_varid, ncf90_get_var, &
              ncf90_inquire, ncf90_inquire_attribute, &
              ncf90_inquire_dimension, ncf90_inquire_variable, &
              ncf90_open, &
              ncf90_open, ncf90_put_att, ncf90_put_var, &
              ncf90_set_fill, ncf90_strerror, ncf90_sync

!   NetCDF parameters
    public :: NCF90_64BIT_OFFSET, NCF90_CHAR, NCF90_CLOBBER, &
              NCF90_DOUBLE, NCF90_ENOTATT, NCF90_FILL_FLOAT, &
              NCF90_FILL_SHORT, NCF90_FLOAT, NCF90_GLOBAL, &
              NCF90_INT, NCF90_INT2, NCF90_NOERR, NCF90_NOFILL, &
              NCF90_NOWRITE, NCF90_REAL, NCF90_SHORT, NCF90_UNLIMITED, &
              NCF90_WRITE

    integer NCF90_64BIT_OFFSET
    integer NCF90_CHAR
    integer NCF90_CLOBBER
    integer NCF90_DOUBLE
    integer NCF90_ENOTATT
    integer NCF90_FILL_FLOAT
    integer NCF90_FILL_SHORT
    integer NCF90_GLOBAL
    integer NCF90_INT
    integer NCF90_INT2
    integer NCF90_MAX_NAME
    integer NCF90_MAX_VAR_DIMS
    integer NCF90_NOERR
    integer NCF90_NOFILL
    integer NCF90_NOWRITE
    integer NCF90_REAL
    integer NCF90_FLOAT
    integer NCF90_SHORT
    integer NCF90_UNLIMITED
    integer NCF90_WRITE
    integer NCF90_NETCDF4

    parameter (NCF90_64BIT_OFFSET = NC_64BIT_OFFSET)
    parameter (NCF90_CHAR = NC_CHAR)
    parameter (NCF90_CLOBBER = NC_CLOBBER)
    parameter (NCF90_DOUBLE = NC_DOUBLE)
    parameter (NCF90_ENOTATT = NC_ENOTATT)
    parameter (NCF90_FILL_FLOAT = NC_FILL_FLOAT)
    parameter (NCF90_FILL_SHORT = NC_FILL_SHORT)
    parameter (NCF90_GLOBAL = NC_GLOBAL)
    parameter (NCF90_INT = NC_INT)
    parameter (NCF90_MAX_NAME = NC_MAX_NAME)
    parameter (NCF90_MAX_VAR_DIMS = NC_MAX_VAR_DIMS)
    parameter (NCF90_INT2 = NC_INT2)
    parameter (NCF90_NOERR = NC_NOERR)
    parameter (NCF90_NOFILL = NC_NOFILL)
    parameter (NCF90_NOWRITE = NC_NOWRITE)
    parameter (NCF90_REAL = NC_REAL)
    parameter (NCF90_FLOAT = NC_FLOAT)
    parameter (NCF90_SHORT = NC_SHORT)
    parameter (NCF90_UNLIMITED = NC_UNLIMITED)
    parameter (NCF90_WRITE = NC_WRITE)
    parameter (NCF90_NETCDF4 = NC_NETCDF4)

contains

    function ncf90_close(ncid)

        ! Closes an open netCDF dataset. If the dataset is in define
        ! mode, NCMPI_ENDDEF will be called before closing.

        integer, intent( in) :: ncid
        integer              :: ncf90_close

        ncf90_close = ncmpi_close(ncid)

    end function ncf90_close

    function ncf90_copy_att(ncid_in, varid_in, name, ncid_out, varid_out)
        integer,             intent( in) :: ncid_in,  varid_in
        character (len = *), intent( in) :: name
        integer,             intent( in) :: ncid_out, varid_out
        integer                          :: ncf90_copy_att

        ! Copies an attribute from one open netCDF dataset to another.

        ncf90_copy_att = ncmpi_copy_att(ncid_in, varid_in, name, &
                                        ncid_out, varid_out)
    end function ncf90_copy_att

    function ncf90_create(path, cmode, ncid, initialsize, bufrsize, &
                          cache_size, cache_nelems, cache_preemption, &
                          comm, info)

         ! Creates a new netCDF dataset, returning a netCDF ID that can
         ! subsequently be used to refer to the netCDF dataset in other
         ! netCDF function calls.

         implicit none
         character (len = *), intent(in) :: path
         integer, intent(in) :: cmode
         integer, intent(out) :: ncid
         integer :: ncf90_create

         ncf90_create = ncmpi_create(MPI_COMM_WORLD, path, cmode, &
                                     MPI_INFO_NULL , ncid)
    end function ncf90_create

    function ncf90_def_dim(ncid, name, len, dimid)

        ! Adds a new dimension to an open netCDF dataset in define mode
        ! and returns a dimension ID.

        integer,             intent( in) :: ncid
        character (len = *), intent( in) :: name
!        integer,             intent( in) :: len
        MPI_Offset,             intent( in) :: len
        integer,             intent(out) :: dimid
        integer                          :: ncf90_def_dim

        ncf90_def_dim = ncmpi_def_dim(ncid, name, len, dimid)

    end function ncf90_def_dim

    function ncf90_def_var(ncid, name, xtype, dimids, varid)

        ! Adds a new variable to an open netCDF dataset in define mode
        ! and returns a variable ID.

        integer, intent(in) :: ncid
        character (len = *), intent(in) :: name
        integer, intent( in) :: xtype
        integer, dimension(:), intent(in) :: dimids
        integer, intent(out) :: varid
        integer :: ncf90_def_var

        ncf90_def_var = ncmpi_def_var(ncid, name, xtype, dimids, varid)

    end function ncf90_def_var

    function ncf90_enddef(ncid)

        ! Takes an open netCDF dataset out of define mode.

        integer,           intent( in) :: ncid
        integer                        :: ncf90_enddef

        ncf90_enddef = ncmpi_enddef(ncid)

    end function ncf90_enddef

    function ncf90_get_att(ncid, varid, name, values)

        ! Gets the value(s) of a netCDF attribute, given
        ! its variable ID and name.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        ! any valid type, scalar or array of rank 1, &
                            intent(out) :: values
        integer                         :: ncf90_get_att

        ncf90_get_att = ncmpi_get_att(ncid, varid, name, values)

    end function ncf90_get_att

    function ncf90_get_var(ncid, varid, nc_type, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.
        !
        ! Note: this implementation currently only caters for
        !       integer and real scalars and arrays.

        integer,                         intent( in)    :: ncid, varid
        integer,                         intent( in)    :: nc_type
        ! any valid type, scalar or array of any rank, &
                                         intent(out)    :: values
!        integer, dimension(:), optional, intent( in) :: start, count
        MPI_Offset, dimension(:), optional, intent( in) :: start, count
        integer                                         :: ncf90_get_var

        ! TODO:
        ! - index is zero based according to pnetcdf library
        !   documentation; assume 1-based for Fortran?
        ! - check NCF90_FLOAT vs NCF90_REAL; same?

        if (present(count)) then
            ! Arrays
            select case (nc_type)
                case (NCF90_INT)
                    ncf90_get_var = &
                        ncmpi_get_vara_int(ncid, varid, &
                                           start, count, values)
                case (NCF90_REAL)
                    ncf90_get_var = &
                        ncmpi_get_vara_real(ncid, varid, &
                                            start, count, values)
            end select
        else
            ! Scalars
            select case (nc_type)
                case (NCF90_INT)
                    ncf90_get_var = &
                        ncmpi_get_var1_int(ncid, varid, &
                                           start, values)
                case (NCF90_REAL)
                    ncf90_get_var = &
                        ncmpi_get_var1_real(ncid, varid, &
                                            start,values)
            end select
        end if

    end function ncf90_get_var

     function ncf90_inq_attname(ncid, varid, attnum, name)

        ! Gets the name of an attribute, given its variable ID and
        ! number.

        integer,             intent( in) :: ncid, varid, attnum
        character (len = *), intent(out) :: name
        integer                          :: ncf90_inq_attname

        ncf90_inq_attname = ncmpi_inq_attname(ncid, varid, attnum, name)

    end function ncf90_inq_attname

    function ncf90_inquire_attribute(ncid, varid, name, xtype, len)

        ! Returns information about a netCDF attribute given the
        ! variable ID and attribute name.

        integer,             intent( in)           :: ncid, varid
        character (len = *), intent( in)           :: name
        integer,             intent(out), optional :: xtype
        MPI_Offset,          intent(out), optional :: len
        integer                                    :: ncf90_inquire_attribute

        ncf90_inquire_attribute = &
            ncmpi_inq_att(ncid, varid, name, xtype, len)

    end function ncf90_inquire_attribute

    function ncf90_inq_dimid(ncid, name, dimid)

        ! Returns (as an argument) the ID of a netCDF dimension,
        ! given the name of the dimension.

        integer,             intent( in) :: ncid
        character (len = *), intent( in) :: name
        integer,             intent(out) :: dimid
        integer                          :: ncf90_inq_dimid

        ncf90_inq_dimid = ncmpi_inq_dimid(ncid, name, dimid)

    end function ncf90_inq_dimid

    function ncf90_inq_varid(ncid, name, varid)

        ! Get the ID of a variable from the name.

        integer, intent(in) :: ncid
        character (len = *), intent( in) :: name
        integer, intent(out) :: varid
        integer :: ncf90_inq_varid

        ncf90_inq_varid = ncmpi_inq_varid(ncid, name, varid)

    end function ncf90_inq_varid

    function ncf90_inquire(ncid, nDimensions, nVariables, nAttributes)

        ! Returns information about an open netCDF dataset, given its
        ! netCDF ID.

        integer,           intent( in) :: ncid
        integer, optional, intent(out) :: nDimensions, nVariables, &
                                          nAttributes
        integer                        :: ncf90_inquire

        ncf90_inquire = ncmpi_inq(ncid, nDimensions, nVariables, &
                                  nAttributes)

    end function ncf90_inquire

    function ncf90_inquire_dimension(ncid, dimid, len)

        ! Returns information about a netCDF dimension given the
        ! dimension ID and attribute name.

        integer,                       intent( in) :: ncid, dimid
        MPI_Offset,          optional, intent(out) :: len
        integer                                    :: ncf90_inquire_dimension

        ncf90_inquire_dimension = ncmpi_inq_dim(ncid, dimid, len)

    end function ncf90_inquire_dimension

    function ncf90_inquire_variable(ncid, varid, name, xtype, ndims, &
                                    dimids, nAtts)

        ! Returns information about a netCDF variable given its ID.

        integer, intent(in) :: ncid, varid
        character (len = *), optional, intent(out) :: name
        integer, optional, intent(out) :: xtype, ndims
        integer, dimension(:), optional, intent(out) :: dimids
        integer, optional, intent(out) :: nAtts
        integer :: ncf90_inquire_variable

        ncf90_inquire_variable = &
            ncmpi_inq_var(ncid, varid, name, xtype, ndims, &
                          dimids, nAtts)

    end function ncf90_inquire_variable

    function ncf90_open(path, mode, ncid)

        ! Opens an existing netCDF dataset for access.

        implicit none
        character (len = *), intent(in) :: path
        integer, intent(in) :: mode
        integer, intent(out) :: ncid
!        integer, optional, intent(inout) :: bufrsize
!        integer, optional, intent(in) :: cache_size, cache_nelems
!        real, optional, intent(in) :: cache_preemption
!        integer, optional, intent(in) :: comm, info
        integer :: ncf90_open

        ncf90_open = ncmpi_open(MPI_COMM_WORLD, path, mode, &
                                MPI_INFO_NULL, ncid)

    end function ncf90_open

    function ncf90_put_att(ncid, varid, name, nc_type value)

        ! Adds or changes a variable attribute or global
        ! attribute of an open netCDF dataset.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        integer,            intent( in) :: nc_type
        character(len = *), intent( in) :: value
        integer                         :: ncf90_put_att

        ncf90_put_att = ncmpi_put_att(ncid, varid, name, &
                                      nc_type, 1, value)

    end function ncf90_put_att

    function ncf90_put_var(ncid, varid, nc_type, values, start, count)

        ! Puts one or more data values into the variable of an open
        ! netCDF dataset that is in data mode.

        ! Note: this implementation currently only caters for
        !       integer and real scalars and arrays.

        integer,                         intent( in) :: ncid, varid
        integer,                         intent( in) :: nc_type
        ! any valid type, scalar or array of any rank, &
                                         intent( in) :: values
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_put_var

        ! TODO:
        ! - index is zero based according to pnetcdf library
        !   documentation; assume 1-based for Fortran?
        ! - check NCF90_FLOAT vs NCF90_REAL; same?

        if (present(count)) then
            ! Arrays
            select case (nc_type)
                case (NCF90_INT)
                    ncf90_put_var = &
                        ncmpi_put_vara_int(ncid, varid, &
                                           start, count, values)
                case (NCF90_REAL)
                    ncf90_put_var = &
                        ncmpi_put_vara_real(ncid, varid, &
                                            start, count, values)
            end select
        else
            ! Scalars
            select case (nc_type)
                case (NCF90_INT)
                    ncf90_put_var = &
                        ncmpi_put_var1_int(ncid, varid, &
                                           start, values)
                case (NCF90_REAL)
                    ncf90_put_var = &
                        ncmpi_put_var1_real(ncid, varid, &
                                            start,values)
            end select
        end if

    end function ncf90_put_var

    function ncf90_set_fill(ncid, fillmode, old_mode)

        ! Sets the fill mode for a netCDF dataset open for
        ! writing and returns the current fill mode in a return
        ! parameter.

        integer, intent( in) :: ncid, fillmode
        integer, intent(out) :: old_mode
        integer              :: ncf90_set_fill

        ncf90_set_fill = ncmpi_set_fill(ncid, fillmode, old_mode)

    end function ncf90_set_fill

    function ncf90_strerror(ncerr)

        ! Returns a static reference to an error message
        ! string corresponding to an integer netCDF error status.

        integer, intent( in) :: ncerr
        character(len = 80)  :: ncf90_strerror

        ncf90_strerror = ncmpi_strerror(ncerr)

    end function ncf90_strerror

    function ncf90_sync(ncid)

        ! Offers a way to synchronize the disk copy of a
        ! netCDF dataset with in-memory buffers.

        integer, intent( in) :: ncid
        integer              :: ncf90_sync

        ncf90_sync = ncmpi_sync(ncid)

    end function ncf90_sync

end module netcdf_m
