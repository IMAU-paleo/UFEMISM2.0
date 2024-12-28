# Makefile for compiling the UFEMISM model
# For "make" syntax search on "make tutorial", or see:
# http://www.gnu.org/software/make/manual/make.html

# ================================
#       Compilation settings
# ================================

# Include settings from (uncommented) file below:

include Makefile_include_local.txt
#include Makefile_include_ecmwf.txt
#include Makefile_include_nix.txt
#include Makefile_include_snellius.txt

# The OBJ_PATH and the MOD_PATH are defined as:
OBJ_PATH = object-files
MOD_PATH = module-files

# Compile with assertions: anything else than yes (no quotes) will
# result in the model being compiled without assertions
DO_ASSERTIONS = yes
ifeq ($(DO_ASSERTIONS),yes)
  F90DEFINITIONS += -DDO_ASSERTIONS
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
$(OBJ_PATH)/%.o : %.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<


# Enable source code to be organised into sub(sub,sub,sub,...)directories
$(OBJ_PATH)/%.o : */%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<

$(OBJ_PATH)/%.o : */*/%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<

$(OBJ_PATH)/%.o : */*/*/%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<

$(OBJ_PATH)/%.o : */*/*/*/%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<

$(OBJ_PATH)/%.o : */*/*/*/*/%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<

$(OBJ_PATH)/%.o : */*/*/*/*/*/%.f90 Makefile | $(OBJ_PATH) $(MOD_PATH)
	$(F90) $(F90FLAGS) $(F90DEFINITIONS) -J${MOD_PATH} -c -o $@ $<


# List of object files

UFEMISM_source = \
  \
  basic/precisions.f90 \
  \
  types/CSR_sparse_matrix_type.f90 \
  types/grid_types.f90 \
  types/mesh_types.f90 \
  types/remapping_types.f90 \
  types/reference_geometry_types.f90 \
  types/scalar_types.f90 \
  types/ice_model_types.f90 \
  types/climate_model_types.f90 \
  types/ocean_model_types.f90 \
  types/SMB_model_types.f90 \
  types/laddie_model_types.f90 \
  types/BMB_model_types.f90 \
  types/LMB_model_types.f90 \
  types/AMB_model_types.f90 \
  types/GIA_model_types.f90 \
  types/basal_inversion_types.f90 \
  types/tracer_tracking_model_types.f90 \
  types/region_types.f90 \
  \
  basic/parameters.f90 \
  basic/mpi_basic.f90 \
  basic/reallocate_mod.f90 \
  basic/control_resources_and_error_messaging.f90 \
  basic/model_configuration.f90 \
  \
  validation/basic_tests/tests_logical.f90 \
  validation/basic_tests/tests_int.f90 \
  validation/basic_tests/tests_dp.f90 \
  validation/basic_tests/tests_CSR.f90 \
  validation/basic_tests/tests_grid.f90 \
  validation/basic_tests/tests_grid_lonlat.f90 \
  validation/basic_tests/tests_mesh.f90 \
  validation/basic_tests/tests_main.f90 \
  validation/basic_tests/assertions_basic.f90 \
  \
  basic/mpi_distributed_memory.f90 \
  \
  basic/math_utilities/flip_mod.f90 \
  basic/math_utilities/permute_mod.f90 \
  basic/math_utilities/projections.f90 \
  basic/math_utilities/matrix_algebra.f90 \
  basic/math_utilities/plane_geometry.f90 \
  basic/math_utilities/ice_geometry_basics.f90 \
  basic/math_utilities/line_integrals.f90 \
  basic/math_utilities/sorting.f90 \
  basic/math_utilities/erf_mod.f90 \
  basic/math_utilities/interpolation.f90 \
  basic/math_utilities/tridiagonal_solver.f90 \
  basic/math_utilities/shape_functions.f90 \
  \
  basic/analytical_solutions.f90 \
  basic/CSR_sparse_matrix_utilities.f90 \
  basic/petsc_basic.f90 \
  \
  basic/grid/mpi_distributed_memory_grid.f90 \
  basic/grid/grid_basic.f90 \
  basic/grid_lonlat_basic.f90 \
  \
  mesh/mesh_memory.f90 \
  mesh/mesh_utilities.f90 \
  mesh/mesh_dummy_meshes.f90 \
  mesh/mesh_ROI_polygons.f90 \
  mesh/Delaunay/check_Delaunay_criterion.f90 \
  mesh/Delaunay/flip_triangles.f90 \
  mesh/Delaunay/move_vertices.f90 \
  mesh/Delaunay/split_border_edges.f90 \
  mesh/Delaunay/split_edges.f90 \
  mesh/Delaunay/split_triangles.f90 \
  mesh/mesh_edges.f90 \
  mesh/Voronoi_mesh/mesh_Voronoi.f90 \
  mesh/mesh_refinement_basic.f90 \
  mesh/mesh_refinement_basic_ROI.f90 \
  mesh/mesh_refinement_fun.f90 \
  mesh/mesh_contiguous_domains.f90 \
  mesh/mesh_Lloyds_algorithm.f90 \
  mesh/mesh_parallel_creation.f90 \
  mesh/mesh_zeta.f90 \
  mesh/discretisation/mesh_translation_tables.f90 \
  mesh/discretisation/mesh_disc_calc_matrix_operators_2D.f90 \
  mesh/discretisation/mesh_disc_calc_matrix_operators_3D.f90 \
  mesh/discretisation/mesh_disc_apply_operators.f90 \
  mesh/remapping/line_tracing_basic.f90 \
  mesh/remapping/line_tracing_grid.f90 \
  mesh/remapping/line_tracing_triangles.f90 \
  mesh/remapping/line_tracing_Voronoi.f90 \
  mesh/remapping/create_maps_grid_mesh.f90 \
  mesh/remapping/create_maps_gridlonlat_mesh.f90 \
  mesh/remapping/create_maps_mesh_mesh.f90 \
  mesh/remapping/apply_maps.f90 \
  mesh/remapping/remapping_main.f90 \
  mesh/mesh_data_smoothing.f90 \
  mesh/mesh_secondary.f90 \
  mesh/creation/reduce_ice_geometry.f90 \
  mesh/creation/mesh_creation_refine_in_ROIs.f90 \
  mesh/creation/mesh_creation_from_reduced_geometry.f90 \
  mesh/creation/mesh_creation_main.f90 \
  \
  netcdf/netcdf_basic.f90 \
  netcdf/netcdf_debug.f90 \
  netcdf/netcdf_input.f90 \
  netcdf/netcdf_output.f90 \
  netcdf/netcdf_resource_tracking.f90 \
  \
  reference_states/reference_geometries.f90 \
  \
  basal_conditions/geothermal_heat_flux.f90 \
  basal_conditions/basal_hydrology.f90 \
  basal_conditions/basal_inversion_H_dHdt_flowline.f90 \
  basal_conditions/basal_inversion_main.f90 \
  basal_conditions/bed_roughness.f90 \
  basal_conditions/sliding_laws.f90 \
  \
  surface_mass_balance/SMB_idealised.f90 \
  surface_mass_balance/SMB_prescribed.f90 \
  surface_mass_balance/SMB_main.f90 \
  \
  ocean/ocean_utilities.f90 \
  ocean/ocean_idealised.f90 \
  ocean/ocean_realistic.f90 \
  ocean/ocean_main.f90 \
  \
  laddie/laddie_utilities.f90 \
  laddie/laddie_main.f90 \
  \
  basal_mass_balance/BMB_idealised.f90 \
  basal_mass_balance/BMB_prescribed.f90 \
  basal_mass_balance/BMB_parameterised.f90 \
  basal_mass_balance/BMB_laddie.f90 \
  basal_mass_balance/BMB_main.f90 \
  \
  lateral_mass_balance/LMB_main.f90 \
  \
  artificial_mass_balance/AMB_main.f90 \
  \
  ice/ice_model_memory.f90 \
  ice/ice_model_utilities.f90 \
  ice/ice_flow_laws.f90 \
  ice/ice_velocity_SIA.f90 \
  ice/ice_velocity_SSA.f90 \
  ice/ice_velocity_DIVA.f90 \
  ice/ice_velocity_BPA.f90 \
  ice/ice_velocity_hybrid_DIVA_BPA.f90 \
  ice/ice_velocity_main.f90 \
  ice/ice_thickness.f90 \
  ice/ice_model_scalars.f90 \
  ice/ice_model_main.f90 \
  \
  thermodynamics/thermodynamics_utilities.f90 \
  thermodynamics/thermodynamics_3D_heat_equation.f90 \
  thermodynamics/thermodynamics_main.f90 \
  \
  climate/climate_idealised.f90 \
  climate/climate_realistic.f90 \
  climate/climate_main.f90 \
  \
  glacial_isostatic_adjustment/GIA_main.f90 \
  \
  tracer_tracking/tracer_tracking_particles.f90 \
  tracer_tracking/tracer_tracking_main.f90 \
  \
  main/main_regional_output.f90 \
  main/UFEMISM_main_model.f90 \
  \
  validation/unit_tests/ut_basic.f90 \
  validation/unit_tests/ut_mpi.f90 \
  validation/unit_tests/ut_petsc.f90 \
  validation/unit_tests/ut_math_utilities.f90 \
  validation/unit_tests/ut_mesh_Delaunay.f90 \
  validation/unit_tests/ut_mesh_refinement.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_grid_start.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_grid_a.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_grid_b.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_grid_cx.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_grid_cy.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_grid.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_tri_start.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_tri_ti.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_tri_vi.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_tri_ei.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_tri.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_Vor_start.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_Vor_vi.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_Vor_ti.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_Vor_ei.f90 \
  validation/unit_tests/ut_mesh_remapping_trace_line_Vor.f90 \
  validation/unit_tests/ut_mesh_remapping_mesh_to_mesh.f90 \
  validation/unit_tests/ut_mesh_remapping.f90 \
  validation/unit_tests/ut_mesh.f90 \
  validation/unit_tests/ut_tracer_tracking.f90 \
  validation/unit_tests/unit_tests.f90 \
  \
  validation/component_tests/ct_create_test_meshes.f90 \
  validation/component_tests/ct_discretisation_mapping_derivatives.f90 \
  validation/component_tests/ct_discretisation.f90 \
  validation/component_tests/ct_remapping_basic.f90 \
  validation/component_tests/ct_remapping_grid_to_mesh.f90 \
  validation/component_tests/ct_remapping_mesh_to_grid.f90 \
  validation/component_tests/ct_remapping_mesh_to_mesh.f90 \
  validation/component_tests/ct_remapping.f90 \
  validation/component_tests/component_tests.f90 \
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
.depend: Makefile */*.f90
	makedepf90 $(F90DEFINITIONS) -b $(OBJ_PATH)/ */*.f90 > .depend;
else
.NOTPARALLEL:
.depend:
	$(warning makedepf90 disabled, enabling .NOTPARALLEL: for the whole makefile )
	touch .depend;
endif
