# For pnetcdf build:
# - module load pnetcdf
# - make PARNETCDF=1 clean run

FC = mpif90

ifdef DEBUG 
FFLAGS = -g -xHost -fpp -ftz
else
FFLAGS = -O -xHost -fpp -ftz
endif

ifdef PARNETCDF
FFLAGS += -DPARNETCDF=1
endif

ifdef PARNETCDF
INC = -I $(PNETCDF_HOME)/include -I ..
LIBS = -L $(PNETCDF_HOME)/lib -lpnetcdf
NETCDF_MOD = pnetcdf_m.o
else
INC = -I $(NETCDF_ROOT)/include -I ..
LIBS = -L $(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
NETCDF_MOD = netcdf_m.o
endif

OBJ = ../$(NETCDF_MOD)

TESTOBJ = sfc_pres_temp_wr.o sfc_pres_temp_rd.o

TESTPRG = ./sfc_pres_temp_wr ./sfc_pres_temp_rd

run: $(TESTPRG) $(TESTOBJ)
	mpirun -np 1 ./sfc_pres_temp_wr
	mpirun -np 1 ./sfc_pres_temp_rd

test: sfc_pres_temp_wr sfc_pres_temp_rd # $(OBJ)

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
	rm -f *.o
	rm -f *.mod
	rm -f *.nc
	rm -f *~
