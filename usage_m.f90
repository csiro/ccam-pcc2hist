module usage_m
   implicit none
   contains
   subroutine usage()
   use mpidata_m
   include 'mpif.h'
   integer ierr
   if (myid==0) then
      write(*,"(a)") &
"Usage: cc2hist [-r res] [-t type] input_file output_file", &
"  cc2hist -h for full list of options and more information."
   end if
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   call mpi_abort(MPI_COMM_WORLD,-1,ierr)
   end subroutine usage

   subroutine help(version)
   use mpidata_m
   include 'mpif.h'
   integer ierr
   character(len=*), intent(in) :: version
   if (myid==0) then
      write(*,"(a)") &
"cc2hist is a program to convert from the CSIRO conformal cubic atmospheric", &
"model history file to a regular lat-lon netcdf file suitable",  &
"for processing by other programs.", &
"", &
"Usage: cc2hist [-h] [-a] [-g grid] [-r res] [-t type] input_file output_file",&
"", &
"Command line options are", &
"", &
" -h for help (this message)", &
"", &
" -a to process entire file (all variables and times)", &
"", &
" -r res where res is the resolution in degrees of the output file", &
"   (Default resolution is approximately equal to model resolution).", &
""

   print*, "        Set type of the input file "
   print*, " -t h (default) Conformal model history file "
   print*, " -t r           Conformal model radstats file "
   print*, " -t c           Conformal model COMPARE diagnostics file "
   print*, " -t s           Shallow water model history file "
   print*, ""
   print*, " Note that the model history file may be a netcdf file. This is"
   print*, " recognised automatically."
   print*, " -g a (default) Model uses A grid"
   print*, " -g c           Model uses C grid"
   print*, ""
   print*, "There are also several options controlling the horizontal interpolation"
   print*, "and extrapolation to pressure levels below the surface"
   print*, " --interp=linear  Use bi-linear horizonal interpolation rather than standard bicubic"
   print*, " --interp=nearest  Use nearest value rather than interpolation"
   print*, " --vextrap=linear Linear extrapolation below surface"
   print*, " --vextrap=none  Use end values rather than extrapolation"
   print*, " --vextrap=missing  Set values below surface as missing"
   print*, "Note that the missing option only works with linear or nearest "
   print*, "horizontal interpolation. Default is to extrapolate temperature"
   print*, "below surface with standard lapse rate and to use end values for"
   print*, "other fields."
   print*
   print*, "If the -a option is not used then the program reads a control namelist"
   print*, "file from standard input. The namelist is called 'input' and the"
   print*, "variables are"
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
   print*, "for normal use. For more information on averaging options see the"
   print*, "header of the file sphere:~dix043/src/lib90/history.f90"
   print*, ""
   print*, "There is also an option for output on the DARLAM Lambert conformal"
   print*, "grid. To use this set darlam_grid=.true. in the input namelist."
   print*, "The program then expects a darlam grid namelist to follow immediately"
   print*, "(i.e. before histnl). This must set values for "
   print*, " il, jl, ds, du, tanl, rnml, stl1, stl2"
   print*, "Note that there are no defaults for these"
   print*, ""
   print*, "Complaints and suggestions to martin.dix@csiro.au"
   print*, ""
   print*, "cc2hist version ", trim(version)
 
   end if
   
   call mpi_barrier(MPI_COMM_WORLD,ierr)
   call mpi_abort(MPI_COMM_WORLD,-1,ierr)

   end subroutine help
end module usage_m
