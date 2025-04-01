module ocean_realistic

  ! Realistic ocean models

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use ocean_model_types                                      , only: type_ocean_model
  use netcdf_io_main
  use ocean_extrapolation                                    , only: extrapolate_ocean_forcing

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine run_ocean_model_realistic( mesh, ice, ocean)
    ! Calculate the ocean
    !
    ! Use an realistic ocean scheme

    implicit none

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'run_ocean_model_realistic'

    ! Add routine to path
    call init_routine( routine_name)

    ! Run the chosen realistic ocean model
    select case (C%choice_ocean_model_realistic)
      case ('snapshot')

        ! Apply extrapolation method if required
        select case (C%choice_ocean_extrapolation_method)
          case('initialisation')
            ! nothing to do here 
          case default
            call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
        end select

      case default
        call crash('unknown choice_ocean_model_realistic "' // trim( C%choice_ocean_model_realistic) // '"')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_realistic

  subroutine initialise_ocean_model_realistic( mesh, ocean, region_name)
    ! Initialise the ocean model
    !
    ! Use an realistic ocean scheme

    implicit none

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ocean_model),                 intent(inout) :: ocean
    character(len=3),                       intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_realistic'
    character(len=256)                                    :: filename_ocean_snapshot

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '     Initialising realistic ocean model "' // &
      colour_string( trim( C%choice_ocean_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic ocean model
    select case (C%choice_ocean_model_realistic)
      case ('snapshot')
        ! Read single-time data from external file

        select case (region_name)
          case ('NAM')
            filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
          case ('EAS')
            filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
          case ('GRL')
            filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
          case ('ANT')
            filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
          case default
            call crash('unknown region_name "' // region_name // '"')
        end select

        ! Fill in  main variables
        call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean, mesh, ocean%T)
        call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean, mesh, ocean%S)

        ! Apply extrapolation method if required
        select case (C%choice_ocean_extrapolation_method)
          case('initialisation')
            call extrapolate_ocean_forcing( mesh, ocean%T)
            call extrapolate_ocean_forcing( mesh, ocean%S) 
          case default
            call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
        end select

      case default
        call crash('unknown choice_ocean_model_realistic "' // trim( C%choice_ocean_model_realistic) // '"')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_realistic

end module ocean_realistic
