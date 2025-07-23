module mesh_creation_from_reduced_geometry

  ! Routines used to create a mesh.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use model_configuration, only: C
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform, refine_mesh_line, refine_mesh_polygon
  use mesh_refinement_fun, only: refine_CalvMIP_shelf_donut
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  use mesh_contiguous_domains, only: enforce_contiguous_process_domains
  use mesh_parallel_creation, only: broadcast_mesh
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use mesh_creation_refine_in_ROIs, only: refine_mesh_in_regions_of_interest

  implicit none

  private

  public :: create_mesh_from_reduced_geometry

contains

  subroutine create_mesh_from_reduced_geometry( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the reduced ice geometry

    ! In/output variables:
    character(len=3),           intent(in)        :: region_name
    character(len=256),         intent(in)        :: name
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in)        :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in)        :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_coastline
    real(dp),                   intent(in)        :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in)        :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(out)       :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'create_mesh_from_reduced_geometry'

    ! Add routine to path
    call init_routine( routine_name)

    ! Choose single-core or parallelised version
    if (mesh%do_singlecore_mesh_creation) then
      call create_mesh_from_reduced_geometry_singlecore( region_name, name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
        xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    else
      call crash('parallelised mesh creation not yet implemented')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_reduced_geometry

  subroutine create_mesh_from_reduced_geometry_singlecore( region_name, name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, &
    xmin, xmax, ymin, ymax, lambda_M, phi_M, beta_stereo, mesh)
    ! Create mesh from the ice geometry lines
    !
    ! Single-core version; all processes generate the same mesh independently

    ! In/output variables:
    character(len=3),           intent(in)        :: region_name
    character(len=256),         intent(in)        :: name
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in)        :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in)        :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in)        :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in)        :: p_line_coastline
    real(dp),                   intent(in)        :: xmin, xmax, ymin, ymax
    real(dp),                   intent(in)        :: lambda_M, phi_M, beta_stereo
    type(type_mesh),            intent(out)       :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'create_mesh_from_reduced_geometry_singlecore'
    real(dp)                                      :: res_max_uniform_applied
    integer                                       :: n1,nn,n2
    real(dp), dimension(:,:  ), allocatable       :: poly
    integer                                       :: i

    ! Add routine to path
    call init_routine( routine_name)

    ! Set mesh configuration
    mesh%resolution_tolerance = C%mesh_resolution_tolerance
    mesh%choice_zeta_grid     = C%choice_zeta_grid
    mesh%nz                   = C%nz
    mesh%zeta_irregular_log_R = C%zeta_irregular_log_R

    ! Single-core mesh generation: let the primary do this,
    ! and then broadcast its result to all the other processes.
    if (par%primary) then

      ! allocate mesh memory
      call allocate_mesh_primary( mesh, name, 1000, 2000)

      ! Initialise the dummy mesh
      call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

      ! == Refine to a uniform resolution; iteratively reduce this,
      ! == and smooth the mesh in between to get a nice, high-quality mesh
      ! ==================================================================

      res_max_uniform_applied = max( xmax-xmin, ymax-ymin)

      do while (.true.)

        ! Reduce the applied uniform resolution
        res_max_uniform_applied = max( res_max_uniform_applied / 2._dp, C%maximum_resolution_uniform)

        ! Refine the mesh to the applied uniform resolution
        call refine_mesh_uniform( mesh, res_max_uniform_applied, C%alpha_min)

        ! Stop refining once we've reached the desired resolution
        if (res_max_uniform_applied == C%maximum_resolution_uniform) exit

      end do ! do while (.true.)

      ! == Refine along the ice geometry lines (grounding line, calving front, etc.)
      ! ============================================================================

      ! Refine the mesh along the ice geometry lines
      call refine_mesh_line( mesh, p_line_grounding_line, C%maximum_resolution_grounding_line, C%grounding_line_width, C%alpha_min)
      call refine_mesh_line( mesh, p_line_calving_front , C%maximum_resolution_calving_front , C%calving_front_width , C%alpha_min)
      call refine_mesh_line( mesh, p_line_ice_front     , C%maximum_resolution_ice_front     , C%ice_front_width     , C%alpha_min)
      call refine_mesh_line( mesh, p_line_coastline     , C%maximum_resolution_coastline     , C%coastline_width     , C%alpha_min)

      ! == Refine along the ice geometry areas (sheet, shelf, etc.)
      ! ===========================================================

        ! Ice sheet
        ! =========

        n1 = 1
        n2 = 0

        do while (n2 < size( poly_mult_sheet,1))

          ! Copy a single polygon from poly_mult
          nn = nint( poly_mult_sheet( n1,1))
          n2 = n1 + nn
          allocate( poly( nn,2))
          poly = poly_mult_sheet( n1+1:n2,:)
          n1 = n2+1

          ! Refine mesh over this single polygon
          call refine_mesh_polygon( mesh, poly, C%maximum_resolution_grounded_ice, C%alpha_min)

          ! Clean up after yourself
          deallocate( poly)

        end do ! do while (n2 < size( poly_mult_sheet,1))

        ! Ice shelf
        ! =========

        n1 = 1
        n2 = 0

        do while (n2 < size( poly_mult_shelf,1))

          ! Copy a single polygon from poly_mult
          nn = nint( poly_mult_shelf( n1,1))
          n2 = n1 + nn
          allocate( poly( nn,2))
          poly = poly_mult_shelf( n1+1:n2,:)
          n1 = n2+1

          ! Refine mesh over this single polygon. Use the ice sheet polygon set as a
          ! "no-refinement" zone to avoid extreme cases where the ice shelf polygon
          ! encompases the ice sheet one (e.g. in the circular domain of CalvMIP)
          if (C%choice_refgeo_PD_ANT == 'idealised') then
            call refine_mesh_polygon( mesh, poly, C%maximum_resolution_floating_ice, C%alpha_min, poly_mult_sheet)
          else
            call refine_mesh_polygon( mesh, poly, C%maximum_resolution_floating_ice, C%alpha_min)
          end if

          ! Clean up after yourself
          deallocate( poly)

        end do ! do while (n2 < size( poly_mult_sheet,1))

      ! == Refine in regions of interest
      ! ================================

      call refine_mesh_in_regions_of_interest( region_name, poly_mult_sheet, poly_mult_shelf, &
        p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, mesh)

      ! == Special cases
      ! ================

      ! DENK DROM : Not very elegant; remove this later and generalise it
      if (C%do_ANT .AND. C%choice_refgeo_PD_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'calvmip_circular') then
        call refine_CalvMIP_shelf_donut( mesh, C%maximum_resolution_grounding_line*2._dp, 70000._dp)
      elseif (C%do_ANT .AND. C%choice_refgeo_PD_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'calvmip_Thule') then
        call refine_CalvMIP_shelf_donut( mesh, C%maximum_resolution_grounding_line*2._dp, 50000._dp)
      end if

      ! == Smooth the mesh
      ! ==================

      do i = 1, C%nit_Lloyds_algorithm
        call Lloyds_algorithm_single_iteration( mesh, C%alpha_min)
      end do

      ! == Enforce contiguous process domains
      ! =====================================

      call enforce_contiguous_process_domains( mesh)

    end if ! if (par%primary) then

    ! Broadcast the primary's mesh
    call broadcast_mesh( mesh)

    ! == Calculate secondary geometry data
    ! ====================================

    ! Calculate all secondary geometry data
    call calc_all_secondary_mesh_data( mesh, lambda_M, phi_M, beta_stereo)

    ! Calculate all matrix operators
    call calc_all_matrix_operators_mesh( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_from_reduced_geometry_singlecore

end module mesh_creation_from_reduced_geometry
