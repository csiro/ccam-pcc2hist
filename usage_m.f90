! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2019 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module usage_m
   implicit none
   contains
   subroutine usage()
   use mpidata_m
#ifdef usempi_mod
   use mpi
#else
   include 'mpif.h'
#endif
   integer ierr
   if (myid==0) then
      write(*,"(a)") &
"Usage: mpirun -np nproc pcc2hist [-h] [-r res] [-v] [input_file] [output_file]", &
"  pcc2hist -h for full list of options and more information."
   end if
   call mpi_barrier(comm_world,ierr)
   call mpi_abort(comm_world,-1,ierr)
   end subroutine usage

   subroutine help(version)
   use mpidata_m
#ifdef usempi_mod
   use mpi
#else
   include 'mpif.h'
#endif
   integer ierr
   character(len=*), intent(in) :: version
   if (myid==0) then
      write(*,"(a)") &
"pcc2hist is a program to convert from the CSIRO conformal cubic atmospheric", &
"model history file to a regular lat-lon netcdf file suitable",  &
"for processing by other programs.", &
"", &
"Usage: mpirun -np nproc pcc2hist [-h] [-r res] [-v] [input_file] [output_file]",&
"", &
"Command line options are", &
"", &
" -h for help (this message)", &
"", &
" -r res where res is the resolution in degrees of the output file", &
"   (Default resolution is approximately equal to model resolution).", &
"", &
" -v for version number", &
""

   print*, "There are also several options controlling the horizontal interpolation"
   print*, "and extrapolation to pressure levels below the surface"
   print*, " --interp=linear   Use bi-linear horizonal interpolation rather than standard bicubic"
   print*, " --interp=nearest  Use nearest value rather than interpolation"
   print*, " --interp=none     No interpolation.  Output on cubic grid"
   print*, " --interp=tapm     Interpolation to TAPM grid"
   print*, " --vextrap=linear  Linear extrapolation below surface"
   print*, " --vextrap=none    Use end values rather than extrapolation"
   print*, " --vextrap=missing Set values below surface as missing"
   print *," --cordex          Format output for CORDEX"
   print*, "Note that the missing option only works with linear or nearest "
   print*, "horizontal interpolation. Default is to extrapolate temperature"
   print*, "below surface with standard lapse rate and to use end values for"
   print*, "other fields."
   print*
   print*, "The program reads a control namelist file called cc.nml. The namelist"
   print*, "is called 'input' and the variables are"
   print*, ""
   print*, " NAME        TYPE     DEFAULT   FUNCTION"
   print*, ""
   print*, "       Selection of time steps"
   print*, " kta         int      0         First timestep to process"
   print*, " ktb         int      9999999   Last timestep to process"
   print*, " ktc         int      -1        Timestep increment (-1 means do all)"
   print*, ""
   print*, " sdate       int      -1        If these are set to positive integers they"
   print*, " stime       int      -1        specify the date/time range of the output"
   print*, " edate       int      -1        data. It's not necessary to set any of kta, "
   print*, " etime       int      -1        ktb or ktc in this case."
   print*, ""
   print*, " ndate       int      -1        Obsolete method for selecting a single"
   print*, " ntime       int      -1        time."
   print*, ""
   print*, "       Boundaries of output region "
   print*, " minlon      real       0 "
   print*, " maxlon      real     360"
   print*, " minlat      real     -90"
   print*, " maxlat      real      90"
   print*, ""
   print*, "       The boundary for the TAPM output region is"
   print*, " lx          int        >0"
   print*, " ly          int        >0"
   print*, " dx          real       >0"
   print*, " dy          real       >0"
   print*, ""
   print*, "       Output levels are the intersection of the following "
   print*, " minlev      int       1"
   print*, " maxlev      int       kl"
   print*, " minsig      real      0."
   print*, " maxsig      real      1."
   print*, ""
   print*, "By default the output variables are on the model sigma levels."
   print*, "It is also possible to interpolate to specified pressure levels by"
   print*, "setting the logical input variable use_plevs to T."
   print*, "The input array plevs specifies the levels to use."
   print*, "E.g."
   print*, "use_plevs = T"
   print*, "plevs = 1000, 850, 500, 200, 100"
   print*, "Note that this uses simple spline extrapolation above and below"
   print*, "the top and bottom model levels. This may give odd results so be"
   print*, "careful."
   print*, ""
   print*, "It is also possible to interpolate to height in meters using"
   print*, "use_meters = T"
   print*, "mlevs = 10, 50, 100, 1000"
   print*, ""
   print*, "For the ocean depths, interpolation can be specified according to"
   print*, "use_depth = T"
   print*, "dlevs = 5, 10, 50, 100, 500, 1000"
   print*, ""
   print*, "A second namelist section allows precise control over which"
   print*, "variables are saved. For example"
   print*, "&histnl"
   print*, 'hnames = "temp", "tsu", "psl"'
   print*, 'hfreq = 1, htype = "inst"'
   print*, "/"
   print*, "will save just temp, tsu and psl. Setting hnames=""all"" will save"
   print*, "all available variables. Particular variables can be excluded by"
   print*, "setting xnames. E.g."
   print*, 'hnames = "all", xnames="sdot"'
   print*, "will save everything but the vertical velocity."
   print*, "The hfreq and htype fields are required and should not be changed"
   print*, "for normal use."
   print*, ""
   print*, "Command line options can be replaced with the following namelist options"
   print*, "ifile                = input_file"
   print*, "ofile                = output_file"
   print*, "hres                 = res"
   print*, "int_default          = 0 (bicubic), 1 (nearest), 2 (bilinear), 5 (none), 9 (tapm)"
   print*, "vextrap              = 0 (default), 1 (linear), 2 (none), 3 (missing),"
   print*, "                       4 (lapse rate)"
   print*, "cf_compliant         = true or false"
   print*, "cordex_compliant     = true or false"
   print*, "save_ccam_parameters = true or false"
   print*, ""
   print*, "pcc2hist version ", trim(version)
 
   end if
   
   call mpi_barrier(comm_world,ierr)
   call mpi_abort(comm_world,-1,ierr)

   end subroutine help
end module usage_m
