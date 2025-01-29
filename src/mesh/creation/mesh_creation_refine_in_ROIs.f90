module mesh_creation_refine_in_ROIs

  ! Routines used to create a mesh.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use model_configuration, only: C
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_refinement_basic, only: refine_mesh_polygon
  use mesh_ROI_polygons
  use mesh_refinement_basic_ROI, only: refine_mesh_polygon_ROI, refine_mesh_line_ROI

  implicit none

  private

  public :: refine_mesh_in_regions_of_interest

contains

  subroutine refine_mesh_in_regions_of_interest( region_name, poly_mult_sheet, poly_mult_shelf, &
    p_line_grounding_line, p_line_calving_front, p_line_ice_front, p_line_coastline, mesh)
    ! Refine the mesh in the specified regions of interest

    ! In/output variables:
    character(len=3),           intent(in   ) :: region_name
    real(dp), dimension(:,:  ), intent(in   ) :: poly_mult_sheet
    real(dp), dimension(:,:  ), intent(in   ) :: poly_mult_shelf
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_grounding_line
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_calving_front
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_ice_front
    real(dp), dimension(:,:  ), intent(in   ) :: p_line_coastline
    type(type_mesh),            intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'refine_mesh_in_regions_of_interest'
    integer                                       :: i
    character(len=256)                            :: all_names_ROI, name_ROI
    real(dp), dimension(:,:  ), allocatable       :: poly_ROI
    integer                                       :: n1,n2,nn, n3
    real(dp), dimension(:,:  ), allocatable       :: poly
    
    ! Add routine to path
    call init_routine( routine_name)

    ! if no regions of interest are specified, do nothing
    if (C%choice_regions_of_interest == '') then
      call finalise_routine( routine_name)
      return
    end if

    all_names_ROI = C%choice_regions_of_interest

    do while (.true.)

      ! == Parse list of input ROIs
      ! ===========================

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      if (i == 0) then
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      else
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      end if

      ! == Check validity of requested ROIs
      ! ===================================

      ! Check if current region is indeed defined in the model
      select case (name_ROI)
        case ('')
          ! No region requested: don't need to do anything
          exit
        case ('PineIsland','Thwaites','Amery','RiiserLarsen','SipleCoast', 'LarsenC','TransMounts','DotsonCrosson', & ! Antarctica
              'Narsarsuaq','Nuuk','Jakobshavn','NGIS','Qaanaaq', &                                                    ! Greenland
              'Patagonia', &                                                                                          ! Patagonia
              'Tijn_test_ISMIP_HOM_A','CalvMIP_quarter')                                                              ! Idealised
          ! List of known regions of interest: these pass the test
        case default
          ! Region not found
          call crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
      end select

      ! == Calculate ROIs
      ! =================

      ! Calculate the polygon describing the specified region of interest
      select case (region_name)
        case ('NAM')
          ! North america

          select case (name_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('EAS')
          ! Eurasia

          select case (name_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('GRL')
          ! Greenland

          select case (name_ROI)
            case ('Narsarsuaq')
              call calc_polygon_Narsarsuaq( poly_ROI)
            case ('Nuuk')
              call calc_polygon_Nuuk( poly_ROI)
            case ('Jakobshavn')
              call calc_polygon_Jakobshavn( poly_ROI)
            case ('NGIS')
              call calc_polygon_NGIS( poly_ROI)
            case ('Qaanaaq')
              call calc_polygon_Qaanaaq( poly_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('ANT')

          select case (name_ROI)
            case ('PineIsland')
              call calc_polygon_Pine_Island_Glacier( poly_ROI)
            case ('Thwaites')
              call calc_polygon_Thwaites_Glacier( poly_ROI)
            case ('Amery')
              call calc_polygon_Amery_ice_shelf( poly_ROI)
            case ('RiiserLarsen')
              call calc_polygon_Riiser_Larsen_ice_shelf( poly_ROI)
            case ('SipleCoast')
              call calc_polygon_Siple_Coast( poly_ROI)
            case ('LarsenC')
              call calc_polygon_Larsen_ice_shelf( poly_ROI)
            case ('TransMounts')
              call calc_polygon_Transantarctic_Mountains( poly_ROI)
            case ('DotsonCrosson')
              call calc_polygon_DotsonCrosson_ice_shelf( poly_ROI)
            case ('Patagonia')
              call calc_polygon_Patagonia( poly_ROI)
            case ('Tijn_test_ISMIP_HOM_A')
              call calc_polygon_Tijn_test_ISMIP_HOM_A( poly_ROI)
            case ('CalvMIP_quarter')
              call calc_polygon_CalvMIP_quarter( poly_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case default
          call crash('unknown region name "' // region_name // '"!')
      end select

      ! Refine the mesh in the specified region of interest
      ! ===================================================

      ! Uniform
      call refine_mesh_polygon( mesh, poly_ROI, C%ROI_maximum_resolution_uniform, C%alpha_min)

      ! Polygons: ice sheet, ice shelf

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
        call refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_grounded_ice, C%alpha_min, poly_ROI)

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
        ! encompases the ice sheet one (e.g. in the domains of CalvMIP)
        if (C%choice_refgeo_PD_ANT == 'idealised') then
          call refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_floating_ice, C%alpha_min, poly_ROI, poly_mult_sheet)
        else
          call refine_mesh_polygon_ROI( mesh, poly, C%ROI_maximum_resolution_floating_ice, C%alpha_min, poly_ROI)
        end if

        ! Clean up after yourself
        deallocate( poly)

      end do ! do while (n2 < size( poly_mult_sheet,1))

      ! Lines: grounding line, calving front, ice front, coastline
      call refine_mesh_line_ROI( mesh, p_line_grounding_line, C%ROI_maximum_resolution_grounding_line, C%ROI_grounding_line_width, C%alpha_min, poly_ROI)
      call refine_mesh_line_ROI( mesh, p_line_calving_front , C%ROI_maximum_resolution_calving_front , C%ROI_calving_front_width , C%alpha_min, poly_ROI)
      call refine_mesh_line_ROI( mesh, p_line_ice_front     , C%ROI_maximum_resolution_ice_front     , C%ROI_ice_front_width     , C%alpha_min, poly_ROI)
      call refine_mesh_line_ROI( mesh, p_line_coastline     , C%ROI_maximum_resolution_coastline     , C%ROI_coastline_width     , C%alpha_min, poly_ROI)

      ! Save poly_ROI in type mesh to use later to determine ROI mask; note this only works for the last one now
      mesh%poly_ROI = poly_ROI
      mesh%npoly_ROI = size(poly_ROI(:,1))
      print*, size(poly_ROI(:,1))

      deallocate( poly_ROI)

      ! if no names are left, we are finished
      if (all_names_ROI == '') exit

    end do ! do while (.true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_mesh_in_regions_of_interest

end module mesh_creation_refine_in_ROIs
