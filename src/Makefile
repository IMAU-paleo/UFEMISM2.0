# Makefile for compiling the UFEMISM model
# For "make" syntax search on "make tutorial", or see:
# http://www.gnu.org/software/make/manual/make.html

# ================================
#       Compilation settings
# ================================

# Include settings from (uncommented) file below:

#include Makefile_include_local.txt
include Makefile_include_mac.txt
#include Makefile_include_gemini.txt
#include Makefile_include_snellius.txt

# The OBJ_PATH and the MOD_PATH are defined as:
OBJ_PATH = object-files
MOD_PATH = module-files

# ==========================
#       Building rules
# ==========================

# Pattern rules for creating object files from *.f90 (fortran 90) source files:
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

$(OBJ_PATH)/%.o : %.F90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



# ========================================
#       Building the UFEMISM_program
# ========================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# List of object files
ifeq ($(DO_SELEN),yes)
UFEMISM_program_objects = \
        basic/precisions.o \
        basic/parameters.o \
        basic/mpi_parallelisation.o \
        basic/utilities.o \
        climate/climate_main.o \
        UFEMISM_main_model.o \
        UFEMISM_program.o
else
UFEMISM_program_objects = \
        basic/precisions.o \
        basic/parameters.o \
        basic/mpi_parallelisation.o \
        basic/utilities.o \
        climate/climate_main.o \
        UFEMISM_main_model.o \
        UFEMISM_program.o
endif

# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)



# ==================================
#       Cleaning up everything
# ==================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod *.lst \
	object-files/basal_conditions/*.o \
	object-files/basal_mass_balance/*.o \
	object-files/basic/*.o \
	object-files/climate/*.o \
	object-files/glacial_isostatic_adjustment/*.o \
	object-files/ice/*.o \
	object-files/mesh/*.o \
	object-files/netcdf/*.o \
	object-files/ocean/*.o \
	object-files/surface_mass_balance/*.o \
	object-files/utilities/*.o



# ===============================
#       Building everything
# ===============================

.PHONY : all

all : $(all_programs)