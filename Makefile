FC = mpif90
FFLAGS = -O -xHost -fpp -ipo -ftz
INC = -I $(NETCDF_ROOT)/include
LIBS = -L $(NETCDF_ROOT)/lib -L $(HDF5_HOME)/lib -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl


SRC = pcc2hist.f90 cc2hist_work.f90 gldata_m.f90 height_m.f90 indices_m.f90 \
ind_m.f90 interp_m.f90 jimcc_m.f90 jimco_m.f90 jim_utils.f90 latltoij_m.f90 \
newmpar_m.f90 nfft_m.f90 parm_m.f90 precis_m.f90 s2p_m.f90 setxyz_m.f90 \
sitop_m.f90 staguv_m.f90 usage_m.f90 xyzinfo_m.f90 utilities.f90 \
mpidata_m.f90

SRCLIB = history.f90 getopt_m.f90 utils_m.f90 ncutils_m.f90 kinds_m.f90 \
physparams.f90 vertutils_m.f90 moistfuncs.f90 hyblevs_m.f90 checkver_m.f90

OBJ = pcc2hist.o cc2hist_work.o gldata_m.o height_m.o indices_m.o ind_m.o \
interp_m.o jimcc_m.o jimco_m.o jim_utils.o latltoij_m.o newmpar_m.o nfft_m.o \
parm_m.o precis_m.o s2p_m.o setxyz_m.o sitop_m.o staguv_m.o usage_m.o \
xyzinfo_m.o utilities.o history.o getopt_m.o utils_m.o ncutils_m.o \
kinds_m.o physparams.o vertutils_m.o moistfuncs.o hyblevs_m.o checkver_m.o \
mpidata_m.o

pcc2hist: $(OBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(OBJ) $(LIBS)


.SUFFIXES:.f90

.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $(LIBS) $<

# Remove mod rule from Modula 2
%.o : %.mod

clean:
	rm *.o *.mod pcc2hist tmpver

# Version string. Dummy dependency to force this to be checked every time
revision.h: FORCE
	echo "   character(len=*), parameter :: cc2hist_revision='SVN-r`svnversion .`'" > tmpver
	# If string contains exported don't overwrite
	# Only update revision.h if it's different to avoid unncessary
	# recompilation
	grep exported tmpver || cmp tmpver revision.h || mv tmpver revision.h

FORCE:

# Module dependencies
pcc2hist.o: parm_m.o checkver_m.o interp_m.o s2p_m.o height_m.o usage_m.o cc2hist_work.o newmpar_m.o getopt_m.o history.o revision.h mpidata_m.o
cc2hist_work.o: moistfuncs.o staguv_m.o vertutils_m.o parm_m.o interp_m.o setxyz_m.o latltoij_m.o indices_m.o xyzinfo_m.o newmpar_m.o sitop_m.o height_m.o physparams.o s2p_m.o history.o gldata_m.o ncutils_m.o precis_m.o mpidata_m.o
height_m.o: utils_m.o physparams.o 
ind_m.o: newmpar_m.o 
interp_m.o: ind_m.o indices_m.o newmpar_m.o
jimcc_m.o: parm_m.o precis_m.o 
jimco_m.o: precis_m.o jim_utils.o nfft_m.o 
jim_utils.o: precis_m.o 
latltoij_m.o: xyzinfo_m.o newmpar_m.o precis_m.o
nfft_m.o: precis_m.o 
parm_m.o: precis_m.o 
physparams.o : precis_m.o
s2p_m.o: gldata_m.o history.o usage_m.o sitop_m.o 
setxyz_m.o: newmpar_m.o indices_m.o parm_m.o precis_m.o ind_m.o xyzinfo_m.o jimco_m.o jimcc_m.o 
sitop_m.o: physparams.o utils_m.o 
staguv_m.o: indices_m.o newmpar_m.o 
xyzinfo_m.o: precis_m.o 
utilities.o: precis_m.o 
history.o: utils_m.o ncutils_m.o mpidata_m.o
utils_m.o: kinds_m.o 
vertutils_m.o: hyblevs_m.o physparams.o 
hyblevs_m.o: physparams.o utils_m.o 
usage_m.o : mpidata_m.o

ccmerge.o: getopt_m.o ncutils_m.o
