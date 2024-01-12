# PCC2Hist (Post-processing CCAM output files)

PCC2Hist is used to convert output files from the Conformal Cubic Atmospheric Model (CCAM) into lat/lon datasets.  PCC2Hist will also calculate additional
information from the CCAM output like relative humidity or dew point
temperature.

## Website

For documentation, see our website at

[https://confluence.csiro.au/display/CCAM/CCAM]

## Dependencies

PCC2hist requires the NetCDF C library and the Message Passing Interface (MPI)
library.

## Building PCC2Hist

To build PCC2Hist with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
