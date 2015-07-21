module pnetcdf_m

!   This module implements parallel netcdf functions, exposing
!   an interface that allows replacement with a different implementation.
!
!   Notes:
!   - The prefix NCF90_ has deliberately been used for all functions
!     and parameters to distinguish between the standard netCDF
!     library and this implementation.
!   - Each function returns a status value.
!   - We can either have one module with N different compile-time
!     conditional blocks and less variation in Makefile (libs, incs)
!     or we can have N modules with more variation in Makefile (libs,
!     incs, dependencies). If N is small, the first may be better. On
!     the other hand, the second allows us to vary the module
!     implementation wildly if necessary. There is a fair amount of
!     minor (e.g. naming) variation, so the second option seems best.
!   - Only those functions (including overloadings) used in the code
!     are supported.
!
!   References:
!   - http://trac.mcs.anl.gov/projects/parallel-netcdf
!   - http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90.html

    use mpi
    use pnetcdf

    implicit none

    private

!   Generic function interfaces

    interface ncf90_def_var
        module procedure ncf90_def_var_with_dimids
        module procedure ncf90_def_var_one_dimid
        module procedure ncf90_def_var_no_dimids
    end interface ncf90_def_var

    interface ncf90_get_att
        module procedure ncf90_get_att_character
        module procedure ncf90_get_att_integer
        module procedure ncf90_get_att_integer_array1D
        module procedure ncf90_get_att_real
        module procedure ncf90_get_att_real_array1D
    end interface ncf90_get_att

    interface ncf90_put_att
        module procedure ncf90_put_att_character
        module procedure ncf90_put_att_integer
        module procedure ncf90_put_att_real
    end interface ncf90_put_att

    interface ncf90_get_var
        module procedure ncf90_get_var_integer
        module procedure ncf90_get_var_integer_array1D
        module procedure ncf90_get_var_integer_array2D
        module procedure ncf90_get_var_integer_array3D
        module procedure ncf90_get_var_real
        module procedure ncf90_get_var_real_array1D
        module procedure ncf90_get_var_real_array2D
        module procedure ncf90_get_var_real_array3D
    end interface ncf90_get_var

    interface ncf90_put_var
        module procedure ncf90_put_var_integer
        module procedure ncf90_put_var_integer_array1D
        module procedure ncf90_put_var_integer_array2D
        module procedure ncf90_put_var_real
        module procedure ncf90_put_var_real_array1D
        module procedure ncf90_put_var_real_array2D
    end interface ncf90_put_var

!   NetCDF function visbility
    public :: ncf90_close, ncf90_copy_att, ncf90_create, &
              ncf90_def_dim, ncf90_def_var, ncf90_enddef, &
              ncf90_get_att, ncf90_inq_attname, ncf90_inq_dimid, &
              ncf90_inq_varid, ncf90_get_var, &
              ncf90_inquire, ncf90_inquire_attribute, &
              ncf90_inquire_dimension, ncf90_inquire_variable, &
              ncf90_open, ncf90_put_att, ncf90_put_var, &
              ncf90_set_fill, ncf90_strerror, ncf90_sync

!   NetCDF parameters
    public :: NCF90_64BIT_OFFSET, NCF90_CHAR, NCF90_CLOBBER, &
              NCF90_DOUBLE, NCF90_ENOTATT, NCF90_FILL_FLOAT, &
              NCF90_FILL_SHORT, NCF90_FLOAT, NCF90_GLOBAL, &
              NCF90_INT, NCF90_INT2, &
              NCF90_MAX_NAME, NCF90_MAX_VAR_DIMS, &
              NCF90_NOERR, NCF90_NOFILL, &
              NCF90_NOWRITE, NCF90_WRITE, NCF90_NETCDF4, &
              NCF90_REAL, NCF90_SHORT, NCF90_UNLIMITED

    integer NCF90_64BIT_OFFSET
    integer NCF90_CHAR
    integer NCF90_CLOBBER
    integer NCF90_DOUBLE
    integer NCF90_ENOTATT
    real    NCF90_FILL_FLOAT
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

    parameter (NCF90_64BIT_OFFSET = NF_64BIT_OFFSET)
    parameter (NCF90_CHAR = NF_CHAR)
    parameter (NCF90_CLOBBER = NF_CLOBBER)
    parameter (NCF90_DOUBLE = NF_DOUBLE)
    parameter (NCF90_ENOTATT = NF_ENOTATT)
    parameter (NCF90_FILL_FLOAT = NF_FILL_FLOAT)
    parameter (NCF90_FILL_SHORT = NF_FILL_SHORT)
    parameter (NCF90_GLOBAL = NF_GLOBAL)
    parameter (NCF90_INT = NF_INT)
    parameter (NCF90_INT2 = NF_INT2)
    parameter (NCF90_MAX_NAME = NF_MAX_NAME)
    parameter (NCF90_MAX_VAR_DIMS = NF_MAX_VAR_DIMS)
    parameter (NCF90_NOERR = NF_NOERR)
    parameter (NCF90_NOFILL = NF_NOFILL)
    parameter (NCF90_NOWRITE = NF_NOWRITE)
    parameter (NCF90_REAL = NF_REAL)
    parameter (NCF90_FLOAT = NF_FLOAT)
    parameter (NCF90_SHORT = NF_SHORT)
    parameter (NCF90_UNLIMITED = NF_UNLIMITED)
    parameter (NCF90_WRITE = NF_WRITE)
    parameter (NCF90_NETCDF4 = NF_NETCDF4)

contains

    function ncf90_close(ncid)

        ! Closes an open netCDF dataset. If the dataset is in define
        ! mode, NF90_ENDDEF will be called before closing.

        integer, intent( in) :: ncid
        integer              :: ncf90_close

        ncf90_close = nfmpi_close(ncid)

    end function ncf90_close

    function ncf90_copy_att(ncid_in, varid_in, name, ncid_out, varid_out)
        integer,             intent( in) :: ncid_in,  varid_in
        character (len = *), intent( in) :: name
        integer,             intent( in) :: ncid_out, varid_out
        integer                          :: ncf90_copy_att

        ! Copies an attribute from one open netCDF dataset to another.

        ncf90_copy_att = nfmpi_copy_att(ncid_in, varid_in, name, &
                                        ncid_out, varid_out)
    end function ncf90_copy_att

    function ncf90_create(path, cmode, ncid)

         ! Creates a new netCDF dataset, returning a netCDF ID that can
         ! subsequently be used to refer to the netCDF dataset in other
         ! netCDF function calls.

         implicit none
         character (len = *), intent(in) :: path
         integer, intent(in) :: cmode
         integer, intent(out) :: ncid
!         integer, optional, intent(in) :: initialsize
!         integer, optional, intent(inout) :: bufrsize
!         integer, optional, intent(in) :: cache_size, cache_nelems
!         real, optional, intent(in) :: cache_preemption
!         integer, optional, intent(in) :: comm, info
         integer :: ncf90_create

         ncf90_create = nfmpi_create(MPI_COMM_WORLD, path, cmode, &
                                     MPI_INFO_NULL , ncid)

    end function ncf90_create

    function ncf90_def_dim(ncid, name, len, dimid)

        ! Adds a new dimension to an open netCDF dataset in define mode
        ! and returns a dimension ID.

        integer,             intent( in) :: ncid
        character (len = *), intent( in) :: name
        integer,             intent( in) :: len
        integer,             intent(out) :: dimid
        integer                          :: ncf90_def_dim

        ncf90_def_dim = nfmpi_def_dim(ncid, name, len, dimid)

    end function ncf90_def_dim

    function ncf90_def_var_with_dimids(ncid, name, xtype, dimids, varid)

        ! Adds a new variable to an open netCDF dataset in define mode
        ! and returns a variable ID.

        integer, intent(in) :: ncid
        character (len = *), intent(in) :: name
        integer, intent( in) :: xtype
        integer, dimension(:), intent(in) :: dimids
        integer, intent(out) :: varid
!        logical, optional, intent(in) :: contiguous
!        integer, optional, dimension(:), intent(in) :: chunksizes
!        integer, optional, intent(in) :: deflate_level
!        logical, optional, intent(in) :: shuffle, fletcher32
!        integer, optional, intent(in) :: endianness
!         integer, optional, intent(in) :: cache_size, cache_nelems, cache_preemption
        integer :: ncf90_def_var_with_dimids

        ncf90_def_var_with_dimids = &
            nfmpi_def_var(ncid, name, xtype, &
            size(dimids), dimids, varid)

    end function ncf90_def_var_with_dimids

    function ncf90_def_var_one_dimid(ncid, name, xtype, dimid, varid)

        ! Adds a new variable to an open netCDF dataset in define mode
        ! and returns a variable ID.

        integer, intent(in) :: ncid
        character (len = *), intent(in) :: name
        integer, intent( in) :: xtype
        integer,  intent(in) :: dimid
        integer, intent(out) :: varid
!        logical, optional, intent(in) :: contiguous
!        integer, optional, dimension(:), intent(in) :: chunksizes
!        integer, optional, intent(in) :: deflate_level
!        logical, optional, intent(in) :: shuffle, fletcher32
!        integer, optional, intent(in) :: endianness
!         integer, optional, intent(in) :: cache_size, cache_nelems, cache_preemption
        integer :: ncf90_def_var_one_dimid

        integer, dimension(:), intent(in) :: dimids

        dimids(1) = dimid

        ncf90_def_var_one_dimid = &
            nfmpi_def_var(ncid, name, xtype, 1, dimids, varid)

    end function ncf90_def_var_one_dimid

    function ncf90_def_var_no_dimids(ncid, name, xtype, varid)

        ! Adds a new variable, with no dimension IDs specified,
        ! to an open netCDF dataset in define mode and returns
        ! a variable ID.

        integer, intent(in) :: ncid
        character (len = *), intent(in) :: name
        integer, intent( in) :: xtype
        integer, intent(out) :: varid
!        logical, optional, intent(in) :: contiguous
!        integer, optional, dimension(:), intent(in) :: chunksizes
!        integer, optional, intent(in) :: deflate_level
!        logical, optional, intent(in) :: shuffle, fletcher32
!        integer, optional, intent(in) :: endianness
!         integer, optional, intent(in) :: cache_size, cache_nelems, cache_preemption
        integer :: ncf90_def_var_no_dimids

        integer, dimension(:), intent(in) :: dimids

        ncf90_def_var_no_dimids = &
            nfmpi_def_var(ncid, name, xtype, 0, dimids, varid)

    end function ncf90_def_var_no_dimids

    function ncf90_enddef(ncid)

        ! Takes an open netCDF dataset out of define mode.

        integer,           intent( in) :: ncid
!        integer, optional, intent( in) :: h_minfree, v_align, v_minfree, r_align
        integer                        :: ncf90_enddef

        ncf90_enddef = nfmpi_enddef(ncid)

    end function ncf90_enddef

    function ncf90_get_att_character(ncid, varid, name, value)

        ! Gets the value of a netCDF attribute, given
        ! its variable ID and name.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        ! any valid type, scalar or array of rank 1, &
        character(len = *), intent(out) :: value
        integer                         :: ncf90_get_att_character

        ncf90_get_att_character = &
            nfmpi_get_att(ncid, varid, name, value)

    end function ncf90_get_att_character

    function ncf90_get_att_integer(ncid, varid, name, value)

        ! Gets the value of a netCDF attribute, given
        ! its variable ID and name.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        ! any valid type, scalar or array of rank 1, &
        integer,            intent(out) :: value
        integer                         :: ncf90_get_att_integer

        ncf90_get_att_integer = &
            nfmpi_get_att(ncid, varid, name, value)

    end function ncf90_get_att_integer

    function ncf90_get_att_integer_array1D(ncid, varid, name, values)

        ! Gets the values of a netCDF attribute, given
        ! its variable ID and name.

        integer,               intent( in) :: ncid, varid
        character(len = *),    intent( in) :: name
        ! any valid type, scalar or array of rank 1, &
        integer, dimension(:), intent(out) :: values
        integer                            :: ncf90_get_att_integer_array1D

        ncf90_get_att_integer_array1D = &
            nfmpi_get_att(ncid, varid, name, values)

    end function ncf90_get_att_integer_array1D

    function ncf90_get_att_real(ncid, varid, name, values)

        ! Gets the value(s) of a netCDF attribute, given
        ! its variable ID and name.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        ! any valid type, scalar or array of rank 1, &
        real,               intent(out) :: values
        integer                         :: ncf90_get_att_real

        ncf90_get_att_real = &
            nfmpi_get_att(ncid, varid, name, values)

    end function ncf90_get_att_real

    function ncf90_get_att_real_array1D(ncid, varid, name, values)

        ! Gets the values of a netCDF attribute, given
        ! its variable ID and name.

        integer,               intent( in) :: ncid, varid
        character(len = *),    intent( in) :: name
        ! any valid type, scalar or array of rank 1, &
        real, dimension(:),    intent(out) :: values
        integer                            :: ncf90_get_att_real_array1D

        ncf90_get_att_real_array1D = &
            nfmpi_get_att(ncid, varid, name, values)

    end function ncf90_get_att_real_array1D

    function ncf90_get_var_integer(ncid, varid, value, start)

        ! Gets a single data value from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer,                         intent(out) :: value
        integer, dimension(:), optional, intent( in) :: start
        integer                                      :: ncf90_get_var_integer

        ! Note: this implementation ignores nc_type
        ncf90_get_var_integer = &
            nfmpi_get_var1_int(ncid, varid, start, values)

    end function ncf90_get_var_integer

    function ncf90_get_var_integer_array1D(ncid, varid, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer, dimension(:),           intent(out) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_get_var_integer_array1D

        ! Note: this implementation ignores nc_type
        ncf90_get_var_integer_array1D = &
            nfmpi_get_vara_int(ncid, varid, start, count, values)

    end function ncf90_get_var_integer_array1D

    function ncf90_get_var_integer_array2D(ncid, varid, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer, dimension(:,:),         intent(out) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_get_var_integer_array2D

        ! Note: this implementation ignores nc_type
        ncf90_get_var_integer_array2D = &
            nfmpi_get_vara_int(ncid, varid, start, count, values)

    end function ncf90_get_var_integer_array2D

    function ncf90_get_var_integer_array3D(ncid, varid, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer, dimension(:,:,:),       intent(out) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_get_var_integer_array3D

        ! Note: this implementation ignores nc_type
        ncf90_get_var_integer_array3D = &
            nfmpi_get_vara_int(ncid, varid, start, count, values)

    end function ncf90_get_var_integer_array3D

    function ncf90_get_var_real(ncid, varid, value, start)

        ! Gets a single data value from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real,                            intent(out) :: value
        integer, dimension(:), optional, intent( in) :: start
        integer                                      :: ncf90_get_var_real

        ! Note: this implementation ignores nc_type
        ncf90_get_var_real = &
            nfmpi_get_var1_real(ncid, varid, start, values)

    end function ncf90_get_var_real

    function ncf90_get_var_real_array1D(ncid, varid, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real, dimension(:),              intent(out) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_get_var_real_array1D

        ! Note: this implementation ignores nc_type
        ncf90_get_var_real_array1D = &
            nfmpi_get_vara_real(ncid, varid, start, count, values)

    end function ncf90_get_var_real_array1D

    function ncf90_get_var_real_array2D(ncid, varid, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real, dimension(:,:),            intent(out) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_get_var_real_array2D

        ! Note: this implementation ignores nc_type
        ncf90_get_var_real_array2D = &
            nfmpi_get_vara_real(ncid, varid, start, count, values)

    end function ncf90_get_var_real_array2D

    function ncf90_get_var_real_array3D(ncid, varid, values, start, count)

        ! Gets one or more data values from a netCDF variable of an
        ! open netCDF dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real, dimension(:,:,:),          intent(out) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_get_var_real_array3D

        ! Note: this implementation ignores nc_type
        ncf90_get_var_real_array3D = &
            nfmpi_get_vara_real(ncid, varid, start, count, values)

    end function ncf90_get_var_real_array3D

    function ncf90_inq_attname(ncid, varid, attnum, name)

        ! Gets the name of an attribute, given its variable ID and
        ! number.

        integer,             intent( in) :: ncid, varid, attnum
        character (len = *), intent(out) :: name
        integer                          :: ncf90_inq_attname

        ncf90_inq_attname = nfmpi_inq_attname(ncid, varid, attnum, name)

    end function ncf90_inq_attname

    function ncf90_inquire_attribute(ncid, varid, name, xtype, len)

        ! Returns information about a netCDF attribute given the
        ! variable ID and attribute name.
        integer,             intent( in)           :: ncid, varid
        character (len = *), intent( in)           :: name
!        integer,             intent(out), optional :: xtype, len, attnum
        integer,             intent(out), optional :: xtype, len
        integer                                    :: ncf90_inquire_attribute

        ncf90_inquire_attribute = &
            nfmpi_inq_att(ncid, varid, name, xtype, len)

    end function ncf90_inquire_attribute

    function ncf90_inq_dimid(ncid, name, dimid)

        ! Returns (as an argument) the ID of a netCDF dimension,
        ! given the name of the dimension.

        integer,             intent( in) :: ncid
        character (len = *), intent( in) :: name
        integer,             intent(out) :: dimid
        integer                          :: ncf90_inq_dimid

        ncf90_inq_dimid = nfmpi_inq_dimid(ncid, name, dimid)

    end function ncf90_inq_dimid

    function ncf90_inq_varid(ncid, name, varid)

        ! Get the ID of a variable from the name.

        integer, intent(in) :: ncid
        character (len = *), intent( in) :: name
        integer, intent(out) :: varid
        integer :: ncf90_inq_varid

        ncf90_inq_varid = nfmpi_inq_varid(ncid, name, varid)

    end function ncf90_inq_varid

    function ncf90_inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimId)

        ! Returns information about an open netCDF dataset, given its
        ! netCDF ID.

        integer,           intent( in) :: ncid
!        integer, optional, intent(out) :: nDimensions, nVariables, &
!                                          nAttributes, unlimitedDimId, &
!                                          formatNum
        integer, optional, intent(out) :: nDimensions, nVariables, &
                                          nAttributes, unlimitedDimId
        integer                        :: ncf90_inquire

        ncf90_inquire = nfmpi_inq(ncid, nDimensions, nVariables, &
                                  nAttributes, unlimitedDimId)

    end function ncf90_inquire

    function ncf90_inquire_dimension(ncid, dimid, len)

        ! Returns information about a netCDF dimension given the
        ! dimension ID and attribute name.

        integer,                       intent( in) :: ncid, dimid
!        character (len = *), optional, intent(out) :: name
        integer,             optional, intent(out) :: len
        integer                                    :: ncf90_inquire_dimension

        ncf90_inquire_dimension = nfmpi_inq_dim(ncid, dimid, len)

    end function ncf90_inquire_dimension

    function ncf90_inquire_variable(ncid, varid, name, xtype, ndims, &
                                    dimids, nAtts)

        ! Returns information about a netCDF variable given its ID.

        integer, intent(in) :: ncid, varid
        character (len = *), optional, intent(out) :: name
        integer, optional, intent(out) :: xtype, ndims
        integer, dimension(:), optional, intent(out) :: dimids
        integer, optional, intent(out) :: nAtts
!        logical, optional, intent(out) :: contiguous
!        integer, optional, dimension(:), intent(out) :: chunksizes
!        integer, optional, intent(out) :: deflate_level
!        logical, optional, intent(out) :: shuffle, fletcher32
!        integer, optional, intent(out) :: endianness
        integer :: ncf90_inquire_variable

        ncf90_inquire_variable = &
            nfmpi_inq_var(ncid, varid, name, xtype, ndims, &
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

        ncf90_open = nfmpi_open(MPI_COMM_WORLD, path, mode, &
                                MPI_INFO_NULL, ncid)

    end function ncf90_open

    function ncf90_put_att_character(ncid, varid, name, value)

        ! Adds or changes a variable attribute or global
        ! attribute of an open netCDF dataset.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        character(len = *), intent( in) :: value
        integer                         :: ncf90_put_att_character

        ncf90_put_att_character = &
            nfmpi_put_att(ncid, varid, name, NF_CHAR, 1, value)

    end function ncf90_put_att_character

    function ncf90_put_att_integer(ncid, varid, name, value)

        ! Adds or changes a variable attribute or global
        ! attribute of an open netCDF dataset.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        integer,            intent( in) :: value
        integer                         :: ncf90_put_att_integer

        ncf90_put_att_integer = &
            nfmpi_put_att(ncid, varid, name, NF_INT, 1, value)

    end function ncf90_put_att_integer

    function ncf90_put_att_real(ncid, varid, name, value)

        ! Adds or changes a variable attribute or global
        ! attribute of an open netCDF dataset.

        integer,            intent( in) :: ncid, varid
        character(len = *), intent( in) :: name
        real,               intent( in) :: value
        integer                         :: ncf90_put_att_real

        ncf90_put_att_real = &
            nfmpi_put_att(ncid, varid, name, NF_REAL, 1, value)

    end function ncf90_put_att_real

    function ncf90_put_var_integer(ncid, varid, value, start)

        ! Puts one data value into the variable of an open netCDF
        ! dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer,                         intent( in) :: value
        integer, dimension(:), optional, intent( in) :: start
        integer                                      :: ncf90_put_var_integer

        ncf90_put_var_integer = &
            nfmpi_put_var1_int(ncid, varid, start, value)

    end function ncf90_put_var_integer

    function ncf90_put_var_integer_array1D(ncid, varid, values, start, count)

        ! Puts one or more data values into the variable of an open netCDF
        ! dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer, dimension(:),           intent( in) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_put_var_integer_array1D

        ncf90_put_var_integer_array1D = &
            nfmpi_put_vara_int(ncid, varid, start, count, values)

    end function ncf90_put_var_integer_array1D

    function ncf90_put_var_integer_array2D(ncid, varid, values, start, count)

        ! Puts one or more data values into the variable of an open netCDF
        ! dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        integer, dimension(:,:),           intent( in) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_put_var_integer_array2D

        ncf90_put_var_integer_array2D = &
            nfmpi_put_vara_int(ncid, varid, start, count, values)

    end function ncf90_put_var_integer_array2D

    function ncf90_put_var_real(ncid, varid, value, start)

        ! Puts one data value into the variable of an open netCDF
        ! dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real,                            intent( in) :: value
        integer, dimension(:), optional, intent( in) :: start
        integer                                      :: ncf90_put_var_real

        ncf90_put_var_real = &
            nfmpi_put_var1_real(ncid, varid, start, value)

    end function ncf90_put_var_real

    function ncf90_put_var_real_array1D(ncid, varid, values, start, count)

        ! Puts one or more data values into the variable of an open netCDF
        ! dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real, dimension(:),              intent( in) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_put_var_real_array1D

        ncf90_put_var_real_array1D = &
            nfmpi_put_vara_real(ncid, varid, start, count, values)

    end function ncf90_put_var_real_array1D

    function ncf90_put_var_real_array2D(ncid, varid, values, start, count)

        ! Puts one or more data values into the variable of an open netCDF
        ! dataset that is in data mode.

        integer,                         intent( in) :: ncid, varid
        ! any valid type, scalar or array of any rank, &
        real, dimension(:,:),              intent( in) :: values
!        integer, dimension(:), optional, intent( in) :: start, count, stride, map
        integer, dimension(:), optional, intent( in) :: start, count
        integer                                      :: ncf90_put_var_real_array2D

        ncf90_put_var_real_array2D = &
            nfmpi_put_vara_real(ncid, varid, start, count, values)

    end function ncf90_put_var_real_array2D

    function ncf90_set_fill(ncid, fillmode, old_mode)

        ! Sets the fill mode for a netCDF dataset open for
        ! writing and returns the current fill mode in a return
        ! parameter.

        integer, intent( in) :: ncid, fillmode
        integer, intent(out) :: old_mode
        integer              :: ncf90_set_fill

        ncf90_set_fill = nfmpi_set_fill(ncid, fillmode, old_mode)

    end function ncf90_set_fill

    function ncf90_strerror(ncerr)

        ! Returns a static reference to an error message
        ! string corresponding to an integer netCDF error status.

        integer, intent( in) :: ncerr
        character(len = 80)  :: ncf90_strerror

        ncf90_strerror = nfmpi_strerror(ncerr)

    end function ncf90_strerror

    function ncf90_sync(ncid)

        ! Offers a way to synchronize the disk copy of a
        ! netCDF dataset with in-memory buffers.

        integer, intent( in) :: ncid
        integer              :: ncf90_sync

        ncf90_sync = nfmpi_sync(ncid)

    end function ncf90_sync

end module netcdf_m
