#$preamble

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
# FCFLAGS = -Ofast -Wl,-stack_size,0x10000000,-stack_addr,0xc0000000
FCFLAGS = -Ofast -fimplicit-none
#FCFLAGS = -fdefault-real-8 -O0 -g  -fbounds-check -Wall  -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -ftrapv -fimplicit-none   -ldislin
FCFLAG2 = -Ofast -fimplicit-none
#-Wall -Wextra -Wno-tabs
#-fopenmp
# flags forall (e.g. look for system .mod files, required in gfortran)
# FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
LDFLAGS = -ldislin

# List of executables to be built within the package
PROGRAMS = example_main

# "make" builds all
all: $(PROGRAMS)

#$intro
# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================

parameters_module.o: precis_mod.o
mod_tree.o: precis_mod.o parameters_module.o
init_module.o: precis_mod.o
mod_strat1.o: precis_mod.o parameters_module.o mod_tree.o
mod_sim1.o: precis_mod.o parameters_module.o init_module.o mod_tree.o mod_strat1.o
example_main.o: precis_mod.o mod_tree.o parameters_module.o init_module.o mod_strat1.o mod_sim1.o 
example_main: precis_mod.o mod_tree.o parameters_module.o init_module.o mod_strat1.o mod_sim1.o

#$conclusion
# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAG2) -o $@ $^ 

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
