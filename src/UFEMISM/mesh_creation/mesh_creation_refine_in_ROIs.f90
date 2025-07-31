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
    integer                                       :: n1,n2,nn
    real(dp), dimension(:,:  ), allocatable       :: poly

    ! Add routine to path
    call init_routine( routine_name)

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
        case ('PineIsland','Thwaites','Amery','RiiserLarsen','SipleCoast', 'LarsenC', &
              'TransMounts','DotsonCrosson', 'Franka_WAIS', 'Dotson_channel', &                                       ! Antarctica
              'Narsarsuaq','Nuuk','Jakobshavn','NGIS','Qaanaaq', &                                                    ! Greenland
              'Patagonia', &                                                                                          ! Patagonia
              'CalvMIP_quarter')                                                              ! Idealised
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
            case ('CalvMIP_quarter')
              call calc_polygon_CalvMIP_quarter( poly_ROI)
            case ('Franka_WAIS')
              call calc_polygon_Franka_WAIS( poly_ROI)
            case ('Dotson_channel')
              call calc_polygon_Dotson_channel( poly_ROI)
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

      ! Clean up after yourself
      deallocate( poly_ROI)

      ! if no names are left, we are finished
      if (all_names_ROI == '') exit

    end do

    if (C%do_refine_TransAntMounts_glaciers) then
      call refine_mesh_over_TransAntarcticMountain_glaciers( mesh)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_mesh_in_regions_of_interest

  subroutine refine_mesh_over_TransAntarcticMountain_glaciers( mesh)
    !< Refine the mesh over the four big glaciers running from the East Antarctic plateau,
    !< through the Transantarctic Mountains, to the Ross ice shelf.

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'refine_mesh_over_TransAntarcticMountain_glaciers'
    real(dp), dimension(:,:), allocatable :: poly_ROI

    ! Add routine to path
    call init_routine( routine_name)

    ! Mulock glacier
    call calc_polygon_Mulock_glacier( poly_ROI)
    call refine_mesh_polygon( mesh, poly_ROI, C%max_res_TransAntMounts_glaciers, C%alpha_min)
    deallocate( poly_ROI)

    ! Byrd glacier
    call calc_polygon_Byrd_glacier( poly_ROI)
    call refine_mesh_polygon( mesh, poly_ROI, C%max_res_TransAntMounts_glaciers, C%alpha_min)
    deallocate( poly_ROI)

    ! Nimrod glacier
    call calc_polygon_Nimrod_glacier( poly_ROI)
    call refine_mesh_polygon( mesh, poly_ROI, C%max_res_TransAntMounts_glaciers, C%alpha_min)
    deallocate( poly_ROI)

    ! Beardmore glacier
    call calc_polygon_Beardmore_glacier( poly_ROI)
    call refine_mesh_polygon( mesh, poly_ROI, C%max_res_TransAntMounts_glaciers, C%alpha_min)
    deallocate( poly_ROI)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine refine_mesh_over_TransAntarcticMountain_glaciers

end module mesh_creation_refine_in_ROIs
