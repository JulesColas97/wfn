# make : compilatioi
# make omp : compilation avec -fopenmp
# possibilite de decommenter lignes pour ajouter d'autres options de compilation a FCFLAGS

# default compilation with gfortran 
# to compile with ifort run make intel
FC = gfortran
# FC = ifort

# on MAC 
# INC = /usr/local/include/
# LIB = /usr/local/lib/

# on Linux 
INC = /usr/include/hdf5/serial
LIB = usr/lib/x86_64-linux-gnu/hdf5/serial

# to compile on Orion uncomment these lines:
# IFC = ifx
# INC = ${EBROOTHDF5}/include
# LIB = ${EBROOTHDF5}/lib



# SHELL = /bin/bash
chemin_modules = './'

# HDF5 = /softs/eb/common/software/HDF5/1.8.20-intel-2018a/
# FC = gfortran
# IFC = ifort
#
# LIB = usr/lib/x86_64-linux-gnu/hdf5/serial
# INCN = /usr/include/hdf5/serial
#
# LIBI =$(HDF5)/lib
# INCI = $(HDF5)/include
# #INC = /usr/lib/x86_64-linux-gnu/hdf5/serial
# # INC = /usr/include/hdf5/serial/
# INC = /usr/local/include/
# LIB = usr/local/lib/
# INC = /usr/lib64/gfortran/modules

# Warnings flags
#FCFLAGS = -Waliasing
#FCFLAGS += -Wampersand
#FCFLAGS += -Wconversion
#FCFLAGS += -Wsurprising
#FCFLAGS += -Wc-binding-type
#FCFLAGS += -Wintrinsics-std
#FCFLAGS += -Wtabs -Wall
#FCFLAGS += -Wintrinsic-shadow
#FCFLAGS += -Wline-truncation
#FCFLAGS += -Wtarget-lifetime
#FCFLAGS += -Wreal-q-constant
#FCFLAGS += std=f2003
#FCFLAGS += -fimplicit-none

# Gfort Compiler Flage s
FCFLAGS += -fcheck=all -fbacktrace
FCFLAGS += -O3 -fdefault-real-8 #-finit-real=nan


# intel compiler Flags
# FCFLAGS += -traceback -g -r8 -O3  -mcmodel=medium -shared-intel
# FCFLAGS += -fpe0 -traceback -g -r8
#FCFLAGS += -traceback -g -r8 -O3 
# FCFLAGS += -O3
# FCFLAGS +=  -ipo
#
#
# Debug flag 
#FCFLAGS = -O0 -g -check all -fpe0 -init=snan,arrays -fp-model strict
#FCFLAGS += -check uninit

EXE = 2D_WAPE
SRC = main.f90
OUT_DIR = out
#all: compil simu

#omp : FCFLAGS += -fopenmp
omp : compil

compil :
	#Compilation des modules, l'ordre est important
	$(FC) -c $(FCFLAGS) $(chemin_modules)/types.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/bspline_kinds_module.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/bspline_sub_module.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/interpo.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/2D_interpolation.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/2D_spline_interpolation.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/1D_interpolation.f90
	$(FC) -c $(FCFLAGS) $(chemin_modules)/real_interpoland.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/lapack.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/utils.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/geometry_utils.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/types_splines.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/splines.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/read_write.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/init_mesh.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/init_flow.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/impedance_model.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/simple_wape.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/wape.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/pwape.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(chemin_modules)/gtpe.f90
	$(FC) -c $(FCFLAGS) -I$(INC) $(SRC)

	#Compilation du programme
	# -I donne un repertoire ou chercher les .mod
#	$(FC) $(SRC)  -I$(chemin_dossier_src_comp) $(chemin_dossier_src_comp)/*.o -o  $(EXE) $(FCFLAGS)
	$(FC) $(FCFLAGS) -L/$(LIB) ./*.o -o $(EXE) -lhdf5 -lhdf5_fortran

intel :

	$(IFC) -c $(IFCFLAGS) $(chemin_modules)/types.f90
	$(IFC) -c $(IFCFLAGS) $(chemin_modules)/bspline_kinds_module.f90
	$(IFC) -c $(IFCFLAGS) $(chemin_modules)/bspline_sub_module.f90
	$(IFC) -c $(IFCFLAGS) $(chemin_modules)/interpo.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/lapack.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/utils.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/types_splines.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/splines.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/read_write.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/init.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/wape.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(chemin_modules)/gtpe.f90
	$(IFC)  -c $(IFCFLAGS) -I$(INCN) $(SRC)
	$(IFC) $(IFCFLAGS) -L$(LIB) ./*.o -o $(EXE) -lhdf5 -lhdf5_fortran

clean :

	rm *.o
	rm *.mod
	rm 2D_WAPE
