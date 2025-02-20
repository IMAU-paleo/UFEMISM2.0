MODULE SMB_model_types

  ! The different data types used in the SMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp

  IMPLICIT NONE

! ===== Types =====
! =================

  TYPE type_SMB_model
    ! The SMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: AlbedoSurf              ! Surface albedo underneath the snow layer (water, rock or ice)
    REAL(dp), DIMENSION(:  ),     ALLOCATABLE :: MeltPreviousYear        ! [m.w.e.] total melt in the previous year
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: FirnDepth               ! [m] depth of the firn layer
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: Rainfall                ! Monthly rainfall (m)
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: Snowfall                ! Monthly snowfall (m)
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: AddedFirn               ! Monthly added firn (m)
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: Melt                    ! Monthly melt (m)
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: Refreezing              ! Monthly refreezing (m)
    REAL(dp), DIMENSION(:  ),     ALLOCATABLE :: Refreezing_year         ! Yearly  refreezing (m)
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: Runoff                  ! Monthly runoff (m)
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: Albedo                  ! Monthly albedo
    REAL(dp), DIMENSION(:  ),     ALLOCATABLE :: Albedo_year             ! Yearly albedo
    REAL(dp), DIMENSION(:,:),     ALLOCATABLE :: SMB                     ! [m] Monthly SMB
    REAL(dp), DIMENSION(:  ),     ALLOCATABLE :: SMB_year                ! Yearly  SMB (m)

    ! Tuning parameters for the IMAU-ITM SMB model (different for each region, set from config)
    REAL(dp), ALLOCATABLE  :: C_abl_constant
    REAL(dp), ALLOCATABLE  :: C_abl_Ts
    REAL(dp), ALLOCATABLE  :: C_abl_Q
    REAL(dp), ALLOCATABLE  :: C_refr
    INTEGER :: wC_abl_constant, wC_abl_Ts, wC_abl_Q, wC_refr
    
    ! Sub-models
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: SMB_correction              ! [m.i.e./yr] Surface mass balance

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_SMB_model

CONTAINS

END MODULE SMB_model_types