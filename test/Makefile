FC = mpif90
FFLAGS = -O -xHost -fpp -ftz
INC = -I $(NETCDF_ROOT)/include -I ..
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff

# TODO: conditionally determine whether to:
# - use netcdf_m.o vs pnetcdf_m.o below (e.g. $(NETCDF_MOD))
# - define parnetcdf macro
# - -lnetcdf[f] vs -lpnetcdf (e.g. LIBS line above +/- INC)

OBJ = ../netcdf_m.o

TESTOBJ = sfc_pres_temp_wr.o sfc_pres_temp_rd.o

TESTPRG = ./sfc_pres_temp_wr ./sfc_pres_temp_rd

run: $(TESTPRG) $(TESTOBJ)
	mpirun -np 1 ./sfc_pres_temp_wr
	mpirun -np 1 ./sfc_pres_temp_rd

sfc_pres_temp_wr: $(OBJ) $(TESTOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(OBJ) $@.o $(LIBS)

sfc_pres_temp_rd: $(OBJ) $(TESTOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(OBJ) $@.o $(LIBS)

.SUFFIXES:.f90

.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $(LIBS) $<

.PHONY: clean
clean:
	rm -f $(TESTPRG)
	rm -f $(TESTOBJ)
	rm -f *.nc
	rm -f $(TESTOBJ)
	rm -f *.nc
