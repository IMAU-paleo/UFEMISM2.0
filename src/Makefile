# Makefile for compiling the UFEMISM model
# For "make" syntax search on "make tutorial", or see:
# http://www.gnu.org/software/make/manual/make.html

# ================================
#       Compilation settings
# ================================

# Include settings from (uncommented) file below:

#include Makefile_include_local.txt
#include Makefile_include_mac.txt
include Makefile_include_nix.txt
#include Makefile_include_gemini.txt
#include Makefile_include_snellius.txt

# The OBJ_PATH and the MOD_PATH are defined as:
OBJ_PATH = object-files
MOD_PATH = module-files

# Compile with SELEN: anything else than yes (no quotes) will skip all
# SELEN-related stuff both in the Makefile and in the model code
DO_SELEN = no
ifeq ($(DO_SELEN),yes)
  F90FLAGS += -DDO_SELEN
endif

# ==========================
#       Building rules
# ==========================

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90 Makefile
	$(F90) $(F90FLAGS) -J${MOD_PATH} -fbackslash -c -o $@ $<



# ========================================
#       Building the UFEMISM_program
# ========================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# List of object files
ifeq ($(DO_SELEN),yes)

  UFEMISM_program_objects = \
    \
    basic/precisions.o \
    basic/parameters.o \
    basic/mpi_basic.o \
    basic/reallocate_mod.o \
    basic/control_resources_and_error_messaging.o \
    basic/model_configuration.o \
    basic/mpi_distributed_memory.o \
    basic/math_utilities.o \
    basic/analytical_solutions.o \
    basic/CSR_sparse_matrix_utilities.o \
    basic/petsc_basic.o \
    basic/grid_basic.o \
    basic/grid_lonlat_basic.o \
    \
    mesh/mesh_types.o \
    mesh/mesh_memory.o \
    mesh/mesh_utilities.o \
    mesh/mesh_Delaunay.o \
    mesh/mesh_edges.o \
    mesh/mesh_refinement.o \
    mesh/mesh_parallel_creation.o \
    mesh/mesh_zeta.o \
    mesh/mesh_operators.o \
    mesh/mesh_remapping.o \
    mesh/mesh_secondary.o \
    mesh/mesh_creation.o \
    \
    netcdf/netcdf_basic.o \
    netcdf/netcdf_input.o \
    netcdf/netcdf_output.o \
    netcdf/netcdf_debug.o \
    netcdf/netcdf_resource_tracking.o \
    \
    ice/ice_model_types.o \
    ice/ice_velocity_SIA.o \
    ice/ice_model_main.o \
    \
    basic/reference_geometries.o \
    \
    validation/unit_tests_mpi.o \
    validation/unit_tests_petsc.o \
    validation/unit_tests_mesh.o \
    validation/unit_tests_netcdf.o \
    validation/unit_tests_ice.o \
    validation/main_validation.o \
    \
    UFEMISM_main_model.o \
    UFEMISM_program.o

else

  UFEMISM_program_objects = \
    \
    basic/precisions.o \
    basic/parameters.o \
    basic/mpi_basic.o \
    basic/reallocate_mod.o \
    basic/control_resources_and_error_messaging.o \
    basic/model_configuration.o \
    basic/mpi_distributed_memory.o \
    basic/math_utilities.o \
    basic/analytical_solutions.o \
    basic/CSR_sparse_matrix_utilities.o \
    basic/petsc_basic.o \
    basic/grid_basic.o \
    basic/grid_lonlat_basic.o \
    \
    mesh/mesh_types.o \
    mesh/mesh_memory.o \
    mesh/mesh_utilities.o \
    mesh/mesh_Delaunay.o \
    mesh/mesh_edges.o \
    mesh/mesh_refinement.o \
    mesh/mesh_parallel_creation.o \
    mesh/mesh_zeta.o \
    mesh/mesh_operators.o \
    mesh/mesh_remapping.o \
    mesh/mesh_secondary.o \
    mesh/mesh_creation.o \
    \
    netcdf/netcdf_basic.o \
    netcdf/netcdf_input.o \
    netcdf/netcdf_output.o \
    netcdf/netcdf_debug.o \
    netcdf/netcdf_resource_tracking.o \
    \
    ice/ice_model_types.o \
    ice/ice_velocity_SIA.o \
    ice/ice_model_main.o \
    \
    basic/reference_geometries.o \
    \
    validation/unit_tests_mpi.o \
    validation/unit_tests_petsc.o \
    validation/unit_tests_mesh.o \
    validation/unit_tests_netcdf.o \
    validation/unit_tests_ice.o \
    validation/main_validation.o \
    \
    UFEMISM_main_model.o \
    UFEMISM_program.o

endif

# This rule is used to build an executable from the object files:
# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)
	$(F90) $(F90FLAGS)    -o $@ $^ $(LDFLAGS)



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
	object-files/utilities/*.o \
	object-files/validation/*.o



# ===============================
#       Building everything
# ===============================

.PHONY : all

all : $(all_programs)



# =============================================================
#      Make dependency graph, better re-compilation + multi
#      threaded compilation (-j flag)
#      Does nothing when makedepf90 is not present
# =============================================================

ifeq ($(shell command -v makedepf90 2>&1 > /dev/null; echo $$?), 0)
include .depend
.depend: Makefile *.f90
	  makedepf90 -b $(OBJ_PATH) *.f90 */*.f90 > .depend;
else
.NOTPARALLEL:
.depend:
	$(warning makedepf90 disabled, enabling .NOTPARALLEL: for the whole makefile )
	touch .depend;
endif
