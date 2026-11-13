# Programs of "A Newton/Picard semi-implicit algorithm for finite-strain plasticity" (revised version).
# Tested with gfortran 15.3.              make            builds the eleven programs
#                                         make check      a quick check against reference values (check/)
#                                         ./run_paper.sh  the runs behind the figures and tables (outputs in out/)
# another compiler: make FC=... FFLAGS=...
FC        = gfortran
FFLAGS    = -O2 -fopenmp -fimplicit-none -ffree-line-length-none
PROGRAMS  = isomap maymap basemap limitscan countmap pctests wkmap symcheck raydiag dencheck nnsens
MODULES   = npicard_guo_kernel.o npicard_plasticity.o

all: $(PROGRAMS)

npicard_guo_kernel.o: npicard_guo_kernel.f90
	$(FC) $(FFLAGS) -c $<

npicard_plasticity.o: npicard_plasticity.f90 npicard_guo_kernel.o
	$(FC) $(FFLAGS) -c $<

$(PROGRAMS): %: %.f90 $(MODULES)
	$(FC) $(FFLAGS) -o $@ $< $(MODULES)

check: all
	./run_paper.sh check

clean:
	rm -f $(PROGRAMS) *.o *.mod

.PHONY: all check clean
