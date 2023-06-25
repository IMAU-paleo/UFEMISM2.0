MODULE region_types

  ! The model region data type

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE grid_basic                                             , ONLY: type_grid
  USE mesh_types                                             , ONLY: type_mesh
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_types                                          , ONLY: type_climate_model
  USE SMB_types                                              , ONLY: type_SMB_model
  USE BMB_types                                              , ONLY: type_BMB_model
  USE scalar_types                                           , ONLY: type_regional_scalars

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_model_region

    ! Metadata
    CHARACTER(LEN=3)                        :: name                        ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name                   ! North America, Eurasia, Greenland, Antarctica

    ! The current time of this particular region.
    REAL(dp)                                :: time

    ! The mesh that all model components define their data on
    TYPE(type_mesh)                         :: mesh

    ! Reference geometries
    TYPE(type_reference_geometry)           :: refgeo_init
    TYPE(type_reference_geometry)           :: refgeo_PD
    TYPE(type_reference_geometry)           :: refgeo_GIAeq

    ! The ice dynamics model
    TYPE(type_ice_model)                    :: ice

    ! The climate model
    TYPE(type_climate_model)                :: climate

    ! The surface mass balance model
    TYPE(type_SMB_model)                    :: SMB

    ! The basal mass balance model
    TYPE(type_BMB_model)                    :: BMB

    ! Scalar data
    TYPE(type_regional_scalars)             :: scalars                     ! Scalar data (e.g. total area, volume, mass balance)

    ! Output
    TYPE(type_grid)                         :: output_grid                 ! The square grid used for gridded output files
    CHARACTER(LEN=256)                      :: output_filename_mesh        ! Name of NetCDF output file (mesh version)
    CHARACTER(LEN=256)                      :: output_filename_grid        ! Name of NetCDF output file (grid version)
    REAL(dp)                                :: output_t_next               ! Time when we should next write to output

  END TYPE type_model_region

CONTAINS

END MODULE region_types