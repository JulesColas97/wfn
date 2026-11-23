# Usage:
#   make                                  # gfortran, linux (default)
#   make COMP=intel MACHINE=marenostrum   # ifx on Mare Nostrum
#   make COMP=intel MACHINE=orion
#   make MACHINE=mac
#   make DEBUG=1                          # debug flags
#   make HDF5_INC=/my/inc HDF5_LIB=/my/lib   # manual override for any other machine
#   make clean

# COMP    : gfortran | intel
# MACHINE : linux | mac | orion | marenostrum | custom
# DEBUG   : 0 | 1
COMP    ?= gfortran
MACHINE ?= linux
DEBUG   ?= 0

EXE = 2D_EWAPE
SRC = main.f90

# ---------------------------------------------------------------
# HDF5 paths per machine (can be overridden from the command line)
# ---------------------------------------------------------------
ifeq ($(MACHINE),linux)
	HDF5_INC ?= /usr/include/hdf5/serial
	HDF5_LIB ?= /usr/lib/x86_64-linux-gnu/hdf5/serial
else ifeq ($(MACHINE),mac)
	HDF5_INC ?= /usr/local/include
	HDF5_LIB ?= /usr/local/lib
else ifeq ($(MACHINE),orion)
	HDF5_INC ?= $(EBROOTHDF5)/include
	HDF5_LIB ?= $(EBROOTHDF5)/lib
else ifeq ($(MACHINE),marenostrum)
	HDF5_INC ?= $(HDF5_DIR)/include
	HDF5_LIB ?= $(HDF5_DIR)/lib
endif

# ---------------------------------------------------------------
# Compiler, flags and linear algebra per compiler
# ---------------------------------------------------------------
ifeq ($(COMP),intel)
  FC      = ifx
  FCFLAGS = -r8 -O3 -mcmodel=medium -shared-intel -traceback -g -ipo -fpe0
  LAPACK  = -qmkl=sequential          # MKL replaces lapack/blas
  ifeq ($(DEBUG),1)
    FCFLAGS = -r8 -O0 -g -check all -fpe0 -init=snan,arrays -fp-model strict -traceback
  endif
else
  FC      = gfortran
  FCFLAGS = -fdefault-real-8 -O3 -march=native -fcheck=all -fbacktrace
  LAPACK  = -llapack -lblas
  ifeq ($(DEBUG),1)
    FCFLAGS = -fdefault-real-8 -O0 -g -fcheck=all -fbacktrace \
              -ffpe-trap=invalid,zero,overflow -finit-real=snan
  endif
endif

FCFLAGS += -I$(HDF5_INC)
LDFLAGS  = -L$(HDF5_LIB) -Wl,-rpath,$(HDF5_LIB)
LDLIBS   = -lhdf5_fortran -lhdf5 $(LAPACK)

# ---------------------------------------------------------------
# Modules, order matters (each one only depends on the previous ones)
# ---------------------------------------------------------------
MODS = kind_mod constants utils_module interp1D_mod interp_plane_mod \
       bspline_sub_module interp_spline_mod input_params atm_abs_module \
       mesh_module flow_module write_outputs iim_utils splitstep_utils \
       trisolve_mod splitstep_pe_wavy

OBJS = $(addsuffix .o,$(MODS)) main.o

.NOTPARALLEL:        # modules must be built in order
.PHONY: compil clean

compil: $(EXE)

$(EXE): $(OBJS)
		$(FC) $(FCFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

clean:
		rm -f *.o *.mod $(EXE)
