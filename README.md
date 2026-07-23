# pcc2hist (converting conformal cubic history file to lat/lon output)

Pcc2hist is used to post-process the Conformal Cubic Atmospheric Model (CCAM).  It converts the history output on the conformal cubic grid to a lat/lon grid typically used for subsequent workflows and analysis.

## Website

For documentation, see our website at

[https://research.csiro.au/ccam/]

## Dependencies

pcc2hist requires the NetCDF C library and the MPI library.

## Building CDFvidar

To build pcc2hist with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
