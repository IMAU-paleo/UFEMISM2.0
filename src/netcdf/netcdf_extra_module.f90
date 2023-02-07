MODULE netcdf_extra_module

! This module contains routines for handling NetCDF files
! that don't conform to the typical x/y-grid, lon/lat-grid, or mesh,
! 2-D, 2-D monthly, or 3-D data fields.
!
! For example, the Laskar insolation solution.
!
! ===== Preamble =====
! ====================

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE petscksp
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
                                             allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                             allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                             allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                             allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                             allocate_shared_int_4D,   allocate_shared_dp_4D, &
                                             allocate_shared_bool_0D,  allocate_shared_bool_1D, &
                                             reallocate_shared_int_0D, reallocate_shared_dp_0D, &
                                             reallocate_shared_int_1D, reallocate_shared_dp_1D, &
                                             reallocate_shared_int_2D, reallocate_shared_dp_2D, &
                                             reallocate_shared_int_3D, reallocate_shared_dp_3D, &
                                             deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_forcing_data, type_SELEN_global, type_ocean_snapshot_global, type_highres_ocean_data
  USE data_types_netcdf_module,        ONLY: type_netcdf_resource_tracker
  USE netcdf_basic_module,             ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_z_ocean, &
                                             field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                                             field_name_options_dim_nV, field_name_options_dim_nTri, field_name_options_dim_nC_mem, &
                                             field_name_options_dim_nAc, field_name_options_dim_two, field_name_options_dim_three, &
                                             field_name_options_dim_six, field_name_options_V, field_name_options_Tri, field_name_options_nC, &
                                             field_name_options_C, field_name_options_niTri, field_name_options_iTri, &
                                             field_name_options_edge_index, field_name_options_Tricc, field_name_options_TriC, &
                                             field_name_options_Tri_edge_index, field_name_options_VAc, field_name_options_Aci, &
                                             field_name_options_iAci, field_name_options_A, field_name_options_R, &
                                             field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                                             field_name_options_SL, field_name_options_Ti, &
                                             open_existing_netcdf_file_for_reading, close_netcdf_file, &
                                             inquire_dim_multiple_options, inquire_var_multiple_options, &
                                             read_var_int_0D, read_var_int_1D, read_var_int_2D, read_var_int_3D, read_var_int_4D, &
                                             read_var_dp_0D , read_var_dp_1D , read_var_dp_2D , read_var_dp_3D , read_var_dp_4D

  IMPLICIT NONE

CONTAINS

! == Insolation solution (e.g. Laskar 2004)
! =========================================

  SUBROUTINE inquire_insolation_data_file( forcing)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data),             INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_insolation_data_file'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( forcing%netcdf_ins%filename, ncid)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( forcing%netcdf_ins%filename, ncid, field_name_options_time , &
      forcing%netcdf_ins%id_dim_time , dim_length = forcing%ins_nyears)
    CALL inquire_dim_multiple_options( forcing%netcdf_ins%filename, ncid, field_name_options_month, &
      forcing%netcdf_ins%id_dim_month)
    CALL inquire_dim_multiple_options( forcing%netcdf_ins%filename, ncid, field_name_options_lat  , &
      forcing%netcdf_ins%id_dim_lat  , dim_length = forcing%ins_nlat)

    ! Inquire variables
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, ncid, field_name_options_time , &
      forcing%netcdf_ins%id_var_time )
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, ncid, field_name_options_month, &
      forcing%netcdf_ins%id_var_month)
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, ncid, field_name_options_lat  , &
      forcing%netcdf_ins%id_var_lat  )
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, ncid, 'Q_TOA', &
      forcing%netcdf_ins%id_var_Q_TOA)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_insolation_data_file

  SUBROUTINE read_insolation_data_file_timeframes( forcing, ti0, ti1, ins_Q_TOA0, ins_Q_TOA1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),             INTENT(INOUT) :: forcing
    INTEGER,                             INTENT(IN)    :: ti0, ti1
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ins_Q_TOA0, ins_Q_TOA1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_insolation_data_file_timeframes'
    INTEGER                                            :: ncid
    INTEGER                                            :: mi, li
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  Q0_with_time,  Q1_with_time
    INTEGER                                            :: wQ0_with_time, wQ1_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Temporary memory to store the data read from the netCDF file
    CALL allocate_shared_dp_3D( 1, 12, forcing%ins_nlat, Q0_with_time, wQ0_with_time)
    CALL allocate_shared_dp_3D( 1, 12, forcing%ins_nlat, Q1_with_time, wQ1_with_time)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( forcing%netcdf_ins%filename, ncid)

    ! Read data
    CALL read_var_dp_3D( forcing%netcdf_ins%filename, ncid, forcing%netcdf_ins%id_var_Q_TOA, Q0_with_time, &
      start = (/ ti0, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /) )
    CALL read_var_dp_3D( forcing%netcdf_ins%filename, ncid, forcing%netcdf_ins%id_var_Q_TOA, Q1_with_time, &
      start = (/ ti1, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /) )

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Remove the time dimension
    IF (par%master) THEN
      DO mi = 1, 12
      DO li = 1, forcing%ins_nlat
        ins_Q_TOA0( li,mi) = Q0_with_time( 1,mi,li)
        ins_Q_TOA1( li,mi) = Q1_with_time( 1,mi,li)
      END DO
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Clean up temporary memory
    CALL deallocate_shared( wQ0_with_time)
    CALL deallocate_shared( wQ1_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_insolation_data_file_timeframes

  SUBROUTINE read_insolation_data_file_time_lat( forcing)
    ! Read only the time and latitude variables from the insolation file,
    !so that later on we know which time frames to read and how to map them.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),             INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_insolation_data_file_time_lat'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_reading( forcing%netcdf_ins%filename, ncid)

    ! Read the data
    CALL read_var_dp_1D( forcing%netcdf_ins%filename, ncid, forcing%netcdf_ins%id_var_time, forcing%ins_time)
    CALL read_var_dp_1D( forcing%netcdf_ins%filename, ncid, forcing%netcdf_ins%id_var_lat , forcing%ins_lat )

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_insolation_data_file_time_lat

! == SELEN
! ========

  ! Global topography input for SELEN
  SUBROUTINE inquire_SELEN_global_topo_file( SELEN)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_SELEN_global_topo_file'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Open the netcdf file
!    INQUIRE(EXIST=file_exists, FILE = TRIM( SELEN%netcdf_topo%filename))
!    IF (.NOT. file_exists) THEN
!      CALL crash('file "' // TRIM( SELEN%netcdf_topo%filename) // '" does not exist!')
!    ELSE
!      CALL open_netcdf_file( SELEN%netcdf_topo%filename, SELEN%netcdf_topo%ncid)
!    END IF
!
!    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
!    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_vi,    SELEN%mesh%nV,     SELEN%netcdf_topo%id_dim_vi)
!    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_ti,    SELEN%mesh%nTri,   SELEN%netcdf_topo%id_dim_ti)
!    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_ci,    SELEN%mesh%nC_mem, SELEN%netcdf_topo%id_dim_ci)
!    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_three, int_dummy,         SELEN%netcdf_topo%id_dim_three)
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
!    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_V,     (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_three /),  SELEN%netcdf_topo%id_var_V       )
!    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_Tri,   (/ SELEN%netcdf_topo%id_dim_ti, SELEN%netcdf_topo%id_dim_three /),  SELEN%netcdf_topo%id_var_Tri     )
!    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_nC,    (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_nC      )
!    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_C,     (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_ci    /),  SELEN%netcdf_topo%id_var_C       )
!    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_niTri, (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_niTri   )
!    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_iTri,  (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_ci    /),  SELEN%netcdf_topo%id_var_iTri    )
!    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_lat,   (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_lat     )
!    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_lon,   (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_lon     )
!    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_Hb,    (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_Hb      )
!    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_ianc,  (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_ianc    )
!
!    ! Close the netcdf file
!    CALL close_netcdf_file(SELEN%netcdf_topo%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_SELEN_global_topo_file

  SUBROUTINE read_SELEN_global_topo_file( SELEN)
    ! Read the init netcdf file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_SELEN_global_topo_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Open the netcdf file
!    CALL open_netcdf_file(SELEN%netcdf_topo%filename, SELEN%netcdf_topo%ncid)
!
!    ! Read the data
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_V,     SELEN%mesh%V,     start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_Tri,   SELEN%mesh%Tri,   start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_nC,    SELEN%mesh%nC,    start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_C,     SELEN%mesh%C,     start = (/ 1, 1 /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_niTri, SELEN%mesh%niTri, start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_iTri,  SELEN%mesh%iTri,  start = (/ 1, 1 /) ))
!
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_lat,   SELEN%mesh%lat,   start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_lon,   SELEN%mesh%lon,   start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_Hb,    SELEN%topo_ref,   start = (/ 1    /) ))
!    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_ianc,  SELEN%mesh%ianc,  start = (/ 1    /) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file(SELEN%netcdf_topo%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_SELEN_global_topo_file

  ! SELEN output file
  SUBROUTINE create_SELEN_output_file( SELEN)
    ! Create a new NetCDF output file for SELEN (on the irregular global mesh)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_SELEN_output_file'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, three, time, ki

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Set time frame index to 1
!    SELEN%output%ti = 1
!
!    ! Set output filename
!    SELEN%output%filename = TRIM(C%output_dir) // 'SELEN_output.nc'
!
!    ! Create a new output file if none exists and, to prevent loss of data,
!    ! stop with an error message if one already exists (not when differences are considered):
!    INQUIRE(EXIST=file_exists, FILE = TRIM(SELEN%output%filename))
!    IF(file_exists) THEN
!      CALL crash('file "' // TRIM( SELEN%output%filename) // '" already exists!')
!    END IF
!
!    ! Create netCDF file
!    CALL handle_error(nf90_create(SELEN%output%filename,IOR(nf90_clobber,nf90_share),SELEN%output%ncid))
!
!    ! Mesh data
!    ! =========
!
!    ! Define dimensions
!    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_vi,           SELEN%mesh%nV,           SELEN%output%id_dim_vi          ) ! Vertex indices
!    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ti,           SELEN%mesh%nTri,         SELEN%output%id_dim_ti          ) ! Triangle indices
!    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ci,           SELEN%mesh%nC_mem,       SELEN%output%id_dim_ci          ) ! Connection indices
!    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_three,        3,                       SELEN%output%id_dim_three       ) ! 3 (each vertex has three coordinates, each triangle has three vertices)
!
!    ! Placeholders for the dimension ID's, for shorter code
!    vi        = SELEN%output%id_dim_vi
!    ti        = SELEN%output%id_dim_ti
!    ci        = SELEN%output%id_dim_ci
!    three     = SELEN%output%id_dim_three
!
!    ! Define variables
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_V,                [vi,  three], SELEN%output%id_var_V,                long_name='Vertex coordinates', units='m')
!    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_Tri,              [ti,  three], SELEN%output%id_var_Tri,              long_name='Vertex indices')
!    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_nC,               [vi        ], SELEN%output%id_var_nC,               long_name='Number of connected vertices')
!    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_C,                [vi,  ci   ], SELEN%output%id_var_C,                long_name='Indices of connected vertices')
!    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_niTri,            [vi        ], SELEN%output%id_var_niTri,            long_name='Number of inverse triangles')
!    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_iTri,             [vi,  ci   ], SELEN%output%id_var_iTri,             long_name='Indices of inverse triangles')
!
!    ! Model output
!    ! ============
!
!    ! Define dimensions
!    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_time,  nf90_unlimited,         SELEN%output%id_dim_time ) ! Time frames
!    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ki,    C%SELEN_irreg_time_n+1, SELEN%output%id_dim_ki   ) ! Window frames
!
!    ! Placeholders for the dimension ID's, for shorter code
!    time  = SELEN%output%id_dim_time
!    ki    = SELEN%output%id_dim_ki
!
!    ! Define dimension variables
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_time,  [time  ], SELEN%output%id_var_time,  long_name='Time', units='years')
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_ki,    [ki    ], SELEN%output%id_var_ki,    long_name='Window frames', units='years')
!
!    ! Define model data variables
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_lat,              [vi             ], SELEN%output%id_var_lat,              long_name='Latitude', units='degrees north')
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_lon,              [vi             ], SELEN%output%id_var_lon,              long_name='Longtitude', units='degrees east')
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_Hi,               [vi,        time], SELEN%output%id_var_Hi,               long_name='Surface load', units='mie')
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_Hi_rel,           [vi,        time], SELEN%output%id_var_Hi_rel,           long_name='Relative surface load', units='mie')
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_U,                [vi,        time], SELEN%output%id_var_U,                long_name='Land surface change', units='m')
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_N,                [vi,        time], SELEN%output%id_var_N,                long_name='Sea surface change', units='m')
!    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_ocean_function,   [vi,        time], SELEN%output%id_var_ocean_function,   long_name='Ocean function (1 = ocean)')
!
!    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_load_history,     [vi,   ki,  time], SELEN%output%id_var_load_history,     long_name='Load history', units='mie')
!
!    ! Leave definition mode:
!    CALL handle_error(nf90_enddef( SELEN%output%ncid))
!
!    ! Write mesh data
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_V,               SELEN%mesh%V             ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Tri,             SELEN%mesh%Tri           ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_nC,              SELEN%mesh%nC            ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_C,               SELEN%mesh%C             ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_niTri,           SELEN%mesh%niTri         ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_iTri,            SELEN%mesh%iTri          ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_lat,             SELEN%mesh%lat           ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_lon,             SELEN%mesh%lon           ))
!
!    ! Window frames
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_ki, (/0._dp, C%SELEN_irreg_time_window/)))
!
!    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
!    CALL handle_error(nf90_sync( SELEN%output%ncid))
!
!    ! Close the file
!    CALL close_netcdf_file(SELEN%output%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_SELEN_output_file

  SUBROUTINE write_to_SELEN_output_file( SELEN, time)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_SELEN_output_file'
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN
    REAL(dp),                       INTENT(IN)    :: time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! DENK DROM
    CALL crash('fixme!')

!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Open the file for writing
!    CALL open_netcdf_file( SELEN%output%filename, SELEN%output%ncid)
!
!    ! Time
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_time, time, start = (/ SELEN%output%ti /)))
!
!    ! Model data
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Hi,             SELEN%Hi_glob,                        start = (/ 1,    SELEN%output%ti/) ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Hi_rel,         SELEN%Hi_rel_glob,                    start = (/ 1,    SELEN%output%ti/) ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_U,              SELEN%U_glob,                         start = (/ 1,    SELEN%output%ti/) ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_N,              SELEN%N_glob,                         start = (/ 1,    SELEN%output%ti/) ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_ocean_function, SELEN%of_glob,                        start = (/ 1,    SELEN%output%ti/) ))
!    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_load_history,   SELEN%ice_loading_history_irreg_glob, start = (/ 1, 1, SELEN%output%ti/) ))
!
!    ! Close the file
!    CALL close_netcdf_file(SELEN%output%ncid)
!
!    ! Increase time frame counter
!    SELEN%output%ti = SELEN%output%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_SELEN_output_file

! == Create and write to resource tracking file
! =============================================

  SUBROUTINE write_to_resource_tracking_file( netcdf, time, tcomp_tot)
    ! Write to the resource tracking output file

    USE configuration_module, ONLY: resource_tracker, mem_use_tot_max

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf
    REAL(dp),                           INTENT(IN)    :: time
    REAL(dp),                           INTENT(IN)    :: tcomp_tot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'write_to_resource_tracking_file'
    INTEGER                                           :: i,n
    INTEGER,  DIMENSION(1024)                         :: path_int_enc

!    IF (.NOT. par%master) RETURN
!
!    ! Open the file for writing
!    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)
!
!    ! Time
!    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, time, start = (/netcdf%ti/)))
!
!    ! Actual variables
!    ! ================
!
!    ! Total model resource use
!    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_tcomp_tot, tcomp_tot      , start = (/ netcdf%ti /) ))
!    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_mem_tot  , mem_use_tot_max, start = (/ netcdf%ti /) ))
!
!    ! Per-subroutine resource use
!
!    n = SIZE( resource_tracker)
!
!    DO i = 1, n
!
!      ! Subroutine name
!      CALL encode_subroutine_path_as_integer( resource_tracker( i)%routine_path, path_int_enc)
!      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_names( i), path_int_enc ))
!
!      ! Computation time
!      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_tcomp( i), resource_tracker( i)%tcomp      , start = (/ netcdf%ti /) ))
!
!      ! Memory use (defined as maximum over the preceding coupling interval)
!      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_mem(   i), resource_tracker( i)%mem_use_max, start = (/ netcdf%ti /) ))
!
!    END DO
!
!    ! Close the file
!    CALL close_netcdf_file( netcdf%ncid)
!
!    ! Increase time frame counter
!    netcdf%ti = netcdf%ti + 1

  END SUBROUTINE write_to_resource_tracking_file

  SUBROUTINE create_resource_tracking_file( netcdf)
    ! Create the resource tracking output file

    USE configuration_module, ONLY: resource_tracker

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'create_resource_tracking_file'
    LOGICAL                                           :: file_exists
    INTEGER                                           :: t,nl
    INTEGER                                           :: i,n
    CHARACTER(LEN=256)                                :: var_name, long_name

!    IF (.NOT. par%master) RETURN
!
!    ! Set time frame index to 1
!    netcdf%ti = 1
!
!    ! Create a new file if none exists and, to prevent loss of data,
!    ! stop with an error message if one already exists (not when differences are considered):
!    netcdf%filename = TRIM(C%output_dir) // '/resource_tracking.nc'
!    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
!    IF (file_exists) THEN
!      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
!    END IF
!
!    ! Create netCDF file
!    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
!    CALL handle_error( nf90_create( netcdf%filename, IOR( nf90_clobber, nf90_share), netcdf%ncid))
!
!    ! Define dimensions:
!    CALL create_dim( netcdf%ncid, netcdf%name_dim_time       , nf90_unlimited, netcdf%id_dim_time       )
!    CALL create_dim( netcdf%ncid, netcdf%name_dim_name_length, 1024          , netcdf%id_dim_name_length)
!
!    ! Placeholders for the dimension ID's, for shorter code
!    t  = netcdf%id_dim_time
!    nl = netcdf%id_dim_name_length
!
!    ! Define variables:
!    ! The order of the CALL statements for the different variables determines their
!    ! order of appearence in the netcdf file.
!
!    ! Dimension variables: time
!    CALL create_double_var( netcdf%ncid, netcdf%name_var_time , [t], netcdf%id_var_time, long_name='Time', units='years'   )
!
!    ! Actual variables
!    ! ================
!
!    ! Total model resource use
!    CALL create_double_var( netcdf%ncid, 'tcomp_tot', [t], netcdf%id_var_tcomp_tot, long_name='Computation time', units='s'    )
!    CALL create_double_var( netcdf%ncid, 'mem_tot'  , [t], netcdf%id_var_mem_tot  , long_name='Memory use'      , units='bytes')
!
!    ! Per-subroutine resource use
!
!    n = SIZE( resource_tracker)
!
!    ALLOCATE( netcdf%id_var_names( n))
!    ALLOCATE( netcdf%id_var_tcomp( n))
!    ALLOCATE( netcdf%id_var_mem(   n))
!
!    DO i = 1, n
!
!      ! Subroutine name
!      ! ===============
!
!      ! Generate variable name (name_00001, name_00002, etc.)
!      var_name(  1:256) = ' '
!      long_name( 1:256) = ' '
!      IF     (i < 10) THEN
!        WRITE( var_name ,'(A,I1)') 'name_0000', i
!      ELSEIF (i < 100) THEN
!        WRITE( var_name,'(A,I2)') 'name_000', i
!      ELSEIF (i < 1000) THEN
!        WRITE( var_name,'(A,I3)') 'name_00', i
!      ELSEIF (i < 10000) THEN
!        WRITE( var_name,'(A,I4)') 'name_0', i
!      ELSEIF (i < 100000) THEN
!        WRITE( var_name,'(A,I5)') 'name_', i
!      END IF
!
!      WRITE( long_name,'(A,I1)') 'Full name of subroutine #', i
!
!      ! Create the variable in the NetCDF file
!      CALL create_int_var( netcdf%ncid, var_name, [nl], netcdf%id_var_names( i),  long_name = long_name)
!
!      ! Computation time
!      ! ================
!
!      ! Generate variable name (tcomp_00001, tcomp_00002, etc.)
!      var_name(  1:256) = ' '
!      long_name( 1:256) = ' '
!      IF     (i < 10) THEN
!        WRITE( var_name ,'(A,I1)') 'tcomp_0000', i
!      ELSEIF (i < 100) THEN
!        WRITE( var_name,'(A,I2)') 'tcomp_000', i
!      ELSEIF (i < 1000) THEN
!        WRITE( var_name,'(A,I3)') 'tcomp_00', i
!      ELSEIF (i < 10000) THEN
!        WRITE( var_name,'(A,I4)') 'tcomp_0', i
!      ELSEIF (i < 100000) THEN
!        WRITE( var_name,'(A,I5)') 'tcomp_', i
!      END IF
!
!      WRITE( long_name,'(A,I5)') 'Computation time for subroutine #', i
!
!      ! Create the variable in the NetCDF file
!      CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_var_tcomp( i),  long_name = long_name, units = 's', missing_value = 0._dp)
!
!      ! Memory use
!      ! ==========
!
!      ! Generate variable name (mem_00001, mem_00002, etc.)
!      var_name(  1:256) = ' '
!      long_name( 1:256) = ' '
!      IF     (i < 10) THEN
!        WRITE( var_name ,'(A,I1)') 'mem_0000', i
!      ELSEIF (i < 100) THEN
!        WRITE( var_name,'(A,I2)') 'mem_000', i
!      ELSEIF (i < 1000) THEN
!        WRITE( var_name,'(A,I3)') 'mem_00', i
!      ELSEIF (i < 10000) THEN
!        WRITE( var_name,'(A,I4)') 'mem_0', i
!      ELSEIF (i < 100000) THEN
!        WRITE( var_name,'(A,I5)') 'mem_', i
!      END IF
!
!      WRITE( long_name,'(A,I5)') 'Memory use for subroutine #', i
!
!      ! Create the variable in the NetCDF file
!      CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_var_mem( i),  long_name = long_name, units = 'bytes', missing_value = 0._dp)
!
!    END DO
!
!    ! Leave definition mode:
!    CALL handle_error(nf90_enddef( netcdf%ncid))
!
!    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
!    CALL handle_error(nf90_sync( netcdf%ncid))
!
!    ! Close the file
!    CALL close_netcdf_file(netcdf%ncid)

  END SUBROUTINE create_resource_tracking_file

  SUBROUTINE encode_subroutine_path_as_integer( subroutine_path, path_int_enc)
    ! Encode the current subroutine path as an integer array so it can be saved as a NetCDF variable
    !
    ! Use the simplest possible encoding:
    !
    !  ' ' = -1 (empty character)
    !
    !    0 = 0
    !    1 = 1
    !    ...
    !    9 = 9
    !
    !    a = 10
    !    b = 11
    !    c = 12
    !    ...
    !    z = 36
    !
    !    A = 37
    !    B = 38
    !    C = 39
    !    ...
    !    Z = 62
    !
    !    _ = 63 (underscore)
    !    / = 64 (forward slash)
    !    ( = 65 (left  bracket)
    !    ) = 66 (right bracket)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=1024),                INTENT(IN)    :: subroutine_path
    INTEGER,  DIMENSION(1024),          INTENT(OUT)   :: path_int_enc

    ! Local variables:
    INTEGER                                           :: i

    path_int_enc = 0

    DO i = 1, 1024

      SELECT CASE ( subroutine_path( i:i))
      CASE( ' ')
        path_int_enc( i) = -1
      CASE( '0')
        path_int_enc( i) = 0
      CASE( '1')
        path_int_enc( i) = 1
      CASE( '2')
        path_int_enc( i) = 2
      CASE( '3')
        path_int_enc( i) = 3
      CASE( '4')
        path_int_enc( i) = 4
      CASE( '5')
        path_int_enc( i) = 5
      CASE( '6')
        path_int_enc( i) = 6
      CASE( '7')
        path_int_enc( i) = 7
      CASE( '8')
        path_int_enc( i) = 8
      CASE( '9')
        path_int_enc( i) = 9
      CASE( 'a')
        path_int_enc( i) = 11
      CASE( 'b')
        path_int_enc( i) = 12
      CASE( 'c')
        path_int_enc( i) = 13
      CASE( 'd')
        path_int_enc( i) = 14
      CASE( 'e')
        path_int_enc( i) = 15
      CASE( 'f')
        path_int_enc( i) = 16
      CASE( 'g')
        path_int_enc( i) = 17
      CASE( 'h')
        path_int_enc( i) = 18
      CASE( 'i')
        path_int_enc( i) = 19
      CASE( 'j')
        path_int_enc( i) = 20
      CASE( 'k')
        path_int_enc( i) = 21
      CASE( 'l')
        path_int_enc( i) = 22
      CASE( 'm')
        path_int_enc( i) = 23
      CASE( 'n')
        path_int_enc( i) = 24
      CASE( 'o')
        path_int_enc( i) = 25
      CASE( 'p')
        path_int_enc( i) = 26
      CASE( 'q')
        path_int_enc( i) = 27
      CASE( 'r')
        path_int_enc( i) = 28
      CASE( 's')
        path_int_enc( i) = 29
      CASE( 't')
        path_int_enc( i) = 30
      CASE( 'u')
        path_int_enc( i) = 31
      CASE( 'v')
        path_int_enc( i) = 32
      CASE( 'w')
        path_int_enc( i) = 33
      CASE( 'x')
        path_int_enc( i) = 34
      CASE( 'y')
        path_int_enc( i) = 35
      CASE( 'z')
        path_int_enc( i) = 36
      CASE( 'A')
        path_int_enc( i) = 37
      CASE( 'B')
        path_int_enc( i) = 38
      CASE( 'C')
        path_int_enc( i) = 39
      CASE( 'D')
        path_int_enc( i) = 40
      CASE( 'E')
        path_int_enc( i) = 41
      CASE( 'F')
        path_int_enc( i) = 42
      CASE( 'G')
        path_int_enc( i) = 43
      CASE( 'H')
        path_int_enc( i) = 44
      CASE( 'I')
        path_int_enc( i) = 45
      CASE( 'J')
        path_int_enc( i) = 46
      CASE( 'K')
        path_int_enc( i) = 47
      CASE( 'L')
        path_int_enc( i) = 48
      CASE( 'M')
        path_int_enc( i) = 49
      CASE( 'N')
        path_int_enc( i) = 50
      CASE( 'O')
        path_int_enc( i) = 51
      CASE( 'P')
        path_int_enc( i) = 52
      CASE( 'Q')
        path_int_enc( i) = 53
      CASE( 'R')
        path_int_enc( i) = 54
      CASE( 'S')
        path_int_enc( i) = 55
      CASE( 'T')
        path_int_enc( i) = 56
      CASE( 'U')
        path_int_enc( i) = 57
      CASE( 'V')
        path_int_enc( i) = 58
      CASE( 'W')
        path_int_enc( i) = 59
      CASE( 'X')
        path_int_enc( i) = 60
      CASE( 'Y')
        path_int_enc( i) = 61
      CASE( 'Z')
        path_int_enc( i) = 62
      CASE( '_')
        path_int_enc( i) = 63
      CASE( '/')
        path_int_enc( i) = 64
      CASE( '(')
        path_int_enc( i) = 65
      CASE( ')')
        path_int_enc( i) = 66
      CASE DEFAULT
        CALL crash('unknown character in routine_path "' // TRIM( subroutine_path) // '"!')
      END SELECT

    END DO

  END SUBROUTINE encode_subroutine_path_as_integer

END MODULE netcdf_extra_module