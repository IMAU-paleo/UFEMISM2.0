module laddie_dummy_domain

  ! Some routines to set up a simple domain for testing

! ===== Preamble =====
! ====================

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use model_configuration, only: C
  use parameters, only: pi
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use ice_model_memory, only: allocate_ice_model
  use masks_mod, only: determine_masks
  use ice_shelf_base_slopes_onesided, only: calc_ice_shelf_base_slopes_onesided
  use ocean_model_types, only: type_ocean_model
  use ocean_main, only: initialise_ocean_vertical_grid
  use laddie_model_types, only: type_laddie_model, type_laddie_timestep
  use laddie_main, only: update_laddie_forcing
  use laddie_utilities, only: allocate_laddie_model, allocate_laddie_timestep
  use mesh_memory, only: allocate_mesh_primary
  use mesh_parallel_creation, only: broadcast_mesh
  use mesh_dummy_meshes, only: initialise_dummy_mesh_16
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use grid_basic, only: setup_square_grid
  use mesh_translation_tables, only: calc_field_to_vector_form_translation_tables

  implicit none

  private

  public :: create_dummy_domain_16

contains

! ===== Main routines =====
! =========================

  subroutine create_dummy_domain_16( mesh, ice, ocean, laddie)
    ! create the domain

    !   v1 ---------- v2 ---------- v3 ---------- v4
    !    |\            |\            |\            |
    !    |   \    t2   |   \    t4   |   \    t6   |
    !    |      \      |      \      |      \      |
    !    |   t1    \   |   t3    \   |   t5    \   |
    !    |            \|            \|            \|
    !   v5 ---------- v6 ---------- v7 ---------- v8
    !    |\            |\            |\            |
    !    |   \    t8   |   \    t10  |   \    t12  |
    !    |      \      |      \      |      \      |
    !    |   t7    \   |   t9    \   |   t11   \   |
    !    |            \|            \|            \|
    !   v9 ---------- v10 --------- v11 --------- v12
    !    |\            |\            |\            |
    !    |   \    t14  |   \    t16  |   \    t18  |
    !    |      \      |      \      |      \      |
    !    |   t13   \   |   t15   \   |   t17   \   |
    !    |            \|            \|            \|
    !   v13 --------- v14 --------- v15 --------- v16

    !   gr ---------- gr ---------- gr ---------- gr
    !    |\            |\            |\            |
    !    |   \         |   \         |   \         |
    !    |      \      |      \      |      \      |
    !    |         \   |         \   |         \   |
    !    |            \|            \|            \|
    !   gr ---------- fl ---------- fl ---------- oc
    !    |\            |\            |\            |
    !    |   \         |   \         |   \         |
    !    |      \      |      \      |      \      |
    !    |         \   |         \   |         \   |
    !    |            \|            \|            \|
    !   gr ---------- fl ---------- fl ---------- oc
    !    |\            |\            |\            |
    !    |   \         |   \         |   \         |
    !    |      \      |      \      |      \      |
    !    |         \   |         \   |         \   |
    !    |            \|            \|            \|
    !   gr ---------- fl ---------- oc ---------- oc

    ! In/output variables
    type(type_mesh),         intent(out  ) :: mesh
    type(type_ice_model),    intent(out  ) :: ice
    type(type_ocean_model),  intent(out  ) :: ocean
    type(type_laddie_model), intent(out  ) :: laddie

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_dummy_domain_16'
    real(dp)                       :: xmin, xmax, ymin, ymax
    character(len=1024)            :: name
    integer                        :: vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create the simple test mesh
    ! ===========================
    name = 'test_mesh'
    xmin = -50e3_dp
    xmax =  50e3_dp
    ymin = -50e3_dp
    ymax =  50e3_dp

    call allocate_mesh_primary( mesh, name, 16, 18)
    call initialise_dummy_mesh_16( mesh, xmin, xmax, ymin, ymax)
    call calc_all_secondary_mesh_data( mesh, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT)
    call calc_all_matrix_operators_mesh( mesh)

    ! Set up ice and bed geometry
    ! ===========================

    call allocate_ice_model( mesh, ice)

    ! Define Hi and Hb
    do vi = mesh%vi1, mesh%vi2
      if (vi == 1) then
        ice%Hi( vi) = 1000._dp
        ice%Hb( vi) = 0._dp
      elseif (vi == 2) then
        ice%Hi( vi) = 1000._dp
        ice%Hb( vi) = 0._dp
      elseif (vi == 3) then
        ice%Hi( vi) = 1000._dp
        ice%Hb( vi) = 0._dp
      elseif (vi == 4) then
        ice%Hi( vi) = 1000._dp
        ice%Hb( vi) = 0._dp
      elseif (vi == 5) then
        ice%Hi( vi) = 1000._dp
        ice%Hb( vi) = -800._dp
      elseif (vi == 6) then
        ice%Hi( vi) = 400._dp
        ice%Hb( vi) = -600._dp
      elseif (vi == 7) then
        ice%Hi( vi) = 100._dp
        ice%Hb( vi) = -400._dp
      elseif (vi == 8) then
        ice%Hi( vi) = 0._dp
        ice%Hb( vi) = -600._dp
      elseif (vi == 9) then
        ice%Hi( vi) = 400._dp
        ice%Hb( vi) = 0._dp
      elseif (vi == 10) then
        ice%Hi( vi) = 200._dp
        ice%Hb( vi) = -400._dp
      elseif (vi == 11) then
        ice%Hi( vi) = 100._dp
        ice%Hb( vi) = -400._dp
      elseif (vi == 12) then
        ice%Hi( vi) = 0._dp
        ice%Hb( vi) = -1000._dp
      elseif (vi == 13) then
        ice%Hi( vi) = 400._dp
        ice%Hb( vi) = 0._dp
      elseif (vi == 14) then
        ice%Hi( vi) = 100._dp
        ice%Hb( vi) = -600._dp
      elseif (vi == 15) then
        ice%Hi( vi) = 0._dp
        ice%Hb( vi) = -1000._dp
      elseif (vi == 16) then
        ice%Hi( vi) = 0._dp
        ice%Hb( vi) = 1000._dp
      end if
    end do

    ! Compute masks
    call determine_masks( mesh, ice)

    ! Extract dHib_dx_b and dHib_dy_b
    call calc_ice_shelf_base_slopes_onesided( mesh, ice)

    ! Define Ti
    do vi = mesh%vi1, mesh%vi2
      ice%Ti( vi, :) = 253._dp
    end do

    ! Set up ocean forcing
    ! ====================

    ! Vertical grid
    C%ocean_vertical_grid_max_depth = 5000._dp
    C%ocean_vertical_grid_dz = 100._dp

    if (allocated( C%z_ocean)) deallocate( C%z_ocean)
    call initialise_ocean_vertical_grid

    ! Ocean temperatures
    if (allocated( ocean%T)) deallocate( ocean%T)
    if (allocated( ocean%S)) deallocate( ocean%S)
    allocate( ocean%T( mesh%vi1:mesh%vi2,C%nz_ocean))
    allocate( ocean%S( mesh%vi1:mesh%vi2,C%nz_ocean))

    ! Define simple forcing
    do vi = mesh%vi1, mesh%vi2
      ocean%T( vi, :) = 0._dp
      ocean%S( vi, :) = 34._dp
    end do

    ! Port ice and ocean info to laddie
    ! =================================

    call allocate_laddie_model( mesh, laddie)

    call update_laddie_forcing( mesh, ice, ocean, laddie)


    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_dummy_domain_16

end module laddie_dummy_domain
