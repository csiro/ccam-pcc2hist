FC = mpif90
FFLAGS = -O -xHost -fpp -ftz

# Use make PARNETCDF=1 for pnetcdf

ifndef PARNETCDF
INC = -I $(NETCDF_ROOT)/include -I ..
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
NETCDF_MOD = netcdf_m.o
else
INC = -I $(HOME)/pnetcdf/include -I ..
LIBS = -L $(HOME)/pnetcdf/lib -lpnetcdf
NETCDF_MOD = pnetcdf_m.o
endif

OBJ = ../$(NETCDF_MOD)

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
