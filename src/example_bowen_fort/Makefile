# Build the TAM-derived planetary radiation model (J. M. Lora; uploaded by user).
#
# Uses the conda Fortran toolchain + netCDF installed in ../../.fortenv
# (see scripts/build_fortran.sh).  Override FC/NFCONFIG to use a system toolchain.

FC       ?= ../../.fortenv/bin/gfortran

# netCDF is not needed: the only `use netcdf` is dead (no nf90_ calls are
# reached in the prescribed-data configuration), so it is commented out and we
# link without netCDF/HDF5 (which also avoids a conda/system glibc mismatch).
# -ffree-line-length-none: the source has lines > 132 cols (free form)
FFLAGS  ?= -O2 -ffree-line-length-none -fallow-argument-mismatch
LDFLAGS ?=

OBJS = planetary_radiation.o reconstructed_stubs.o \
       planetary_radiation_driver.o run_planetary_radiation.o
BIN  = run_planetary_radiation

all: $(BIN)

$(BIN): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(LDFLAGS)

planetary_radiation.o: planetary_radiation.F90
	$(FC) $(FFLAGS) -c $<

reconstructed_stubs.o: reconstructed_stubs.F90
	$(FC) $(FFLAGS) -c $<

planetary_radiation_driver.o: planetary_radiation_driver.F90 planetary_radiation.o reconstructed_stubs.o
	$(FC) $(FFLAGS) -c $<

run_planetary_radiation.o: run_planetary_radiation.F90 planetary_radiation_driver.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(BIN)

.PHONY: all clean
