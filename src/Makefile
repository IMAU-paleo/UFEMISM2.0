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
#include Makefile_include_book.txt

# The OBJ_PATH and the MOD_PATH are defined as:
OBJ_PATH = object-files
MOD_PATH = module-files

# Compile with SELEN: anything else than yes (no quotes) will skip all
# SELEN-related stuff both in the Makefile and in the model code
DO_SELEN = no
ifeq ($(DO_SELEN),yes)
  F90DEFINITIONS += -DDO_SELEN
endif

# ========================================
#       Building the UFEMISM_program
# ========================================

# Adding this program to the list of all programs:
all_programs += UFEMISM_program

# ===============================
#       Building everything
# ===============================

.PHONY : all

all : $(all_programs)


# ==========================
#       Building rules
# ==========================

$(OBJ_PATH) $(MOD_PATH):
	mkdir -p $@

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : */%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -fbackslash -c -o $@ $<

$(OBJ_PATH)/%.o : %.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -fbackslash -c -o $@ $<


# List of object files

UFEMISM_source = \
  \
  basic/precisions.f90 \
  basic/parameters.f90 \
  basic/mpi_basic.f90 \
  basic/reallocate_mod.f90 \
  basic/control_resources_and_error_messaging.f90 \
  basic/model_configuration.f90 \
  basic/mpi_distributed_memory.f90 \
  basic/math_utilities.f90 \
  basic/analytical_solutions.f90 \
  basic/CSR_sparse_matrix_utilities.f90 \
  basic/petsc_basic.f90 \
  basic/grid_basic.f90 \
  basic/grid_lonlat_basic.f90 \
  \
  types/mesh_types.f90 \
  types/scalar_types.f90 \
  types/ice_model_types.f90 \
  types/climate_model_types.f90 \
  types/ocean_model_types.f90 \
  types/SMB_model_types.f90 \
  types/BMB_model_types.f90 \
  types/reference_geometry_types.f90 \
  types/region_types.f90 \
  \
  mesh/mesh_memory.f90 \
  mesh/mesh_utilities.f90 \
  mesh/mesh_Delaunay.f90 \
  mesh/mesh_edges.f90 \
  mesh/mesh_refinement.f90 \
  mesh/mesh_parallel_creation.f90 \
  mesh/mesh_zeta.f90 \
  mesh/mesh_operators.f90 \
  mesh/mesh_remapping.f90 \
  mesh/mesh_secondary.f90 \
  mesh/mesh_creation.f90 \
  \
  netcdf/netcdf_basic.f90 \
  netcdf/netcdf_input.f90 \
  netcdf/netcdf_output.f90 \
  netcdf/netcdf_debug.f90 \
  netcdf/netcdf_resource_tracking.f90 \
  \
  reference_states/reference_geometries.f90 \
  \
  basal_conditions/basal_hydrology.f90 \
  basal_conditions/bed_roughness.f90 \
  basal_conditions/sliding_laws.f90 \
  basal_conditions/basal_conditions_main.f90 \
  \
  ice/ice_model_memory.f90 \
  ice/ice_model_scalars.f90 \
  ice/ice_model_utilities.f90 \
  ice/ice_flow_laws.f90 \
  ice/ice_velocity_SIA.f90 \
  ice/ice_velocity_SSA.f90 \
  ice/ice_velocity_DIVA.f90 \
  ice/ice_velocity_BPA.f90 \
  ice/ice_velocity_main.f90 \
  ice/ice_thickness.f90 \
  ice/ice_model_main.f90 \
  \
  thermodynamics/thermodynamics_utilities.f90 \
  thermodynamics/thermodynamics_3D_heat_equation.f90 \
  thermodynamics/thermodynamics_main.f90 \
  \
  climate/climate_main.f90 \
  \
  ocean/ocean_main.f90 \
  \
  surface_mass_balance/SMB_main.f90 \
  \
  basal_mass_balance/BMB_main.f90 \
  \
  main/main_regional_output.f90 \
  main/UFEMISM_main_model.f90 \
  \
  validation/unit_tests_mpi.f90 \
  validation/unit_tests_petsc.f90 \
  validation/unit_tests_mesh.f90 \
  validation/unit_tests_netcdf.f90 \
  validation/unit_tests_ice.f90 \
  validation/main_validation.f90 \
  \
  main/UFEMISM_program.f90

UFEMISM_source_files_=$(notdir $(UFEMISM_source))
UFEMISM_program_objects=$(UFEMISM_source_files_:%.f90=%.o)

# This rule is used to build an executable from the object files:
# The executable (the target) depends on the object files:
UFEMISM_program : $(UFEMISM_program_objects:%=$(OBJ_PATH)/%)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -o $@ $^ $(LDFLAGS)

# ==================================
#       Cleaning up everything
# ==================================

.PHONY : clean

clean :
	rm -frv $(all_programs) $(OBJ_PATH) $(MOD_PATH) *.mod *.lst



# =============================================================
#      Make dependency graph, better re-compilation + multi
#      threaded compilation (-j flag)
#      Does nothing when makedepf90 is not present
# =============================================================

ifeq ($(shell command -v makedepf90 2>&1 > /dev/null; echo $$?), 0)
include .depend
.depend: Makefile *.f90
	makedepf90 $(F90DEFINITIONS) -b $(OBJ_PATH)/ *.f90 */*.f90 > .depend;
else
.NOTPARALLEL:
.depend:
	$(warning makedepf90 disabled, enabling .NOTPARALLEL: for the whole makefile )
	touch .depend;
endif
