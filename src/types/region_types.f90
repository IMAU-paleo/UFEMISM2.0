MODULE region_types

  ! The model region data type

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mesh_types                                             , ONLY: type_mesh
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  USE ice_model_types                                        , ONLY: type_ice_model
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

    ! The surface mass balance model
    TYPE(type_SMB_model)                    :: SMB

    ! The basal mass balance model
    TYPE(type_BMB_model)                    :: BMB

    ! Scalar data
    TYPE(type_regional_scalars)             :: scalars                     ! Scalar data (e.g. total area, volume, mass balance)

    ! Output
    CHARACTER(LEN=256)                      :: output_filename             ! Name of NetCDF output file
    REAL(dp)                                :: output_t_next               ! Time when we should next write to output

  END TYPE type_model_region

CONTAINS

END MODULE region_types