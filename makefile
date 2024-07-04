FC = mpif90

# Intel compilter options
ifneq ($(CUSTOM),yes)
FHOST = -xHost
ifeq ($(BROADWELL),yes)
FHOST = -xCORE-AVX2 -align array32byte
endif
ifeq ($(SKYLAKE),yes)
FHOST = -xSKYLAKE-AVX512 -align array64byte
endif
ifeq ($(CASCADELAKE),yes)
FHOST = -xCASCADELAKE -align array64byte
endif
ifeq ($(SAPPHIRERAPIDS),yes)
FHOST = -xSAPPHIRERAPIDS -align array64byte
endif
ifeq ($(ZEN3),yes)
FHOST = -axCORE-AVX2 -align array32byte
endif
ifeq ($(MAGNUS),yes)
FC = ftn
FHOST = -xHost
endif
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dshare_ifullg
endif
FFLAGS = -O3 $(FHOST) -ftz -fp-model precise -traceback $(MPIFLAG)
INC = -I $(NETCDF_ROOT)/include
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf
PPFLAG90 = -fpp
DEBUGFLAG = -check all -debug all -fpe0
endif

# Gfortran compiler options
ifeq ($(GFORTRAN),yes)
MPIFC = gfortran
MPIF77 = gfortran
ifeq ($(NOMPI3),yes)
MPIFLAG =
else
MPIFLAG = -Dshare_ifullg
endif
FFLAGS = -O2 -fallow-argument-mismatch -march=native -fbacktrace $(MPIFLAG)
PPFLAG90 = -x f95-cpp-input
DEBUGFLAG = -g -Wall -Wextra -fbounds-check
endif

# Gfortran compiler options
ifeq ($(SETONIX),yes)
MPIFC = ftn
MPIF77 = ftn
FC = ftn
FFLAGS = -O2 -mtune=native -march=native -fbacktrace -fallow-argument-mismatch -Dsafe
INC =
LIBS = -lnetcdf
PPFLAG90 = -x f95-cpp-input
DEBUGFLAG = -g -Wall -Wextra -fbounds-check
endif

# Cray compiler options
ifeq ($(CRAY),yes)
FC = ftn
FFLAGS = -h noomp -Dusenc_mod
PPFLAG90 = -eZ
DEBUGFLAG =
endif

ifeq ($(MAUI),yes)
FC = ftn
FHOST = -xSKYLAKE-AVX512
MPIFLAG = 
FFLAGS = -O3 $(FHOST) -ftz -fp-model precise -traceback $(MPIFLAG)
INC = -I $(NETCDF_ROOT)/include
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
PPFLAG90 = -fpp
DEBUGFLAG = -check all -debug all -fpe0
endif


# Options for building with VAMPIRTrace
ifeq ($(VT),yes)
FC = vtfort -vt:fc mpif90 -vt:inst manual
FFLAGS += -Dvampir -DVTRACE
endif

# Testing - I/O and fpmodel
ifeq ($(TEST),yes)
FFLAGS += $(DEBUGFLAG)
endif

OBJ = pcc2hist.o cc2hist_work.o gldata_m.o height_m.o indices_m.o ind_m.o \
interp_m.o jimcc_m.o jimco_m.o jim_utils.o latltoij_m.o newmpar_m.o nfft_m.o \
parm_m.o precis_m.o s2p_m.o setxyz_m.o sitop_m.o staguv_m.o usage_m.o \
xyzinfo_m.o utilities.o history.o getopt_m.o utils_m.o ncutils_m.o \
kinds_m.o physparams.o vertutils_m.o moistfuncs.o hyblevs_m.o checkver_m.o \
mpidata_m.o shdata_m.o logging_m.o netcdf_m.o zenith.o

pcc2hist: $(OBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(OBJ) $(LIBS)


.SUFFIXES:.f90

.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $(LIBS) $(PPFLAG90) $<

# Remove mod rule from Modula 2
%.o : %.mod

clean:
	rm *.o *.mod pcc2hist tmpver

# Version string. Dummy dependency to force this to be checked every time
revision.h: FORCE
	echo "character(len=*), parameter :: cc2hist_revision= &" > tmpver
	echo "'GIT-`git log | head -3 | tail -1`" "`git log | head -1`'" >> tmpver
	# If string contains exported don't overwrite
	# Only update revision.h if it's different to avoid unncessary
	# recompilation
	cmp tmpver revision.h || mv tmpver revision.h

FORCE:

# Module dependencies
pcc2hist.o: parm_m.o checkver_m.o interp_m.o s2p_m.o height_m.o usage_m.o cc2hist_work.o newmpar_m.o getopt_m.o history.o revision.h mpidata_m.o logging_m.o
cc2hist_work.o: moistfuncs.o staguv_m.o vertutils_m.o parm_m.o interp_m.o setxyz_m.o latltoij_m.o indices_m.o xyzinfo_m.o newmpar_m.o sitop_m.o height_m.o physparams.o s2p_m.o history.o gldata_m.o ncutils_m.o precis_m.o mpidata_m.o shdata_m.o logging_m.o netcdf_m.o zenith.o
height_m.o: netcdf_m.o utils_m.o physparams.o s2p_m.o
ind_m.o: newmpar_m.o 
interp_m.o: ind_m.o indices_m.o newmpar_m.o precis_m.o logging_m.o netcdf_m.o
jimcc_m.o: parm_m.o precis_m.o 
jimco_m.o: precis_m.o jim_utils.o nfft_m.o 
jim_utils.o: precis_m.o 
latltoij_m.o: xyzinfo_m.o newmpar_m.o precis_m.o
ncutils_m.o: netcdf_m.o
nfft_m.o: precis_m.o 
parm_m.o: precis_m.o 
physparams.o : precis_m.o
s2p_m.o: gldata_m.o history.o usage_m.o sitop_m.o logging_m.o
setxyz_m.o: newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o
sitop_m.o: physparams.o utils_m.o netcdf_m.o
staguv_m.o: indices_m.o newmpar_m.o 
xyzinfo_m.o: precis_m.o 
utilities.o: precis_m.o 
history.o: interp_m.o logging_m.o mpidata_m.o ncutils_m.o netcdf_m.o newmpar_m.o parm_m.o utils_m.o
utils_m.o: kinds_m.o 
vertutils_m.o: hyblevs_m.o physparams.o 
usage_m.o : mpidata_m.o
shdata_m.o: mpidata_m.o 
logging_m.o: mpidata_m.o
