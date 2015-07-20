FC = mpif90
FFLAGS = -O -xHost -fpp -ftz
INC = -I $(NETCDF_ROOT)/include -I ..
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff

# TODO: conditionally determine whether to:
# - use netcdf_m.o vs pnetcdf_m.o below (e.g. $(NETCDF_MOD))
# - define parnetcdf macro
# - -lnetcdf[f] vs -lpnetcdf (e.g. LIBS line above +/- INC)

OBJ = $(shell find .. -name "*.o" | grep -v pcc2hist | grep -v "test/")

TESTOBJ = sfc_pres_temp_wr.o

TESTPRG = sfc_pres_temp_wr

all:
	run

run: $(TESTPRG)
	mpirun -np 1 ./sfc_pres_temp_wr

sfc_pres_temp_wr: $(OBJ) $(TESTOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(OBJ) $(TESTOBJ) $(LIBS)
	
.SUFFIXES:.f90

.f90.o:
	$(FC) -c $(FFLAGS) $(INC) $(LIBS) $<

.PHONY: clean
clean:
	rm -r $(TESTOBJ)
	rm -f *.nc
