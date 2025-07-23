module remapping_transects

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use transect_types, only: type_transect
  use remapping_types, only: type_map
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, &
    finalise_matrix_CSR_dist
  use petsc_basic, only: mat_CSR2petsc
  use mesh_utilities, only: find_containing_triangle, find_containing_vertex
  use plane_geometry, only: triangle_area
  use apply_maps, only: Atlas
  use apply_maps_transects, only: apply_map_mesh_vertices_to_transect_2D, apply_map_mesh_vertices_to_transect_3D, &
    apply_map_mesh_triangles_to_transect_2D, apply_map_mesh_triangles_to_transect_3D

  implicit none

  private

  public :: create_map_from_mesh_vertices_to_transect_trilin, &
    create_map_from_mesh_vertices_to_transect_nearest_neighbour, &
    create_map_from_mesh_triangles_to_transect, &
    map_from_mesh_vertices_to_transect_2D, map_from_mesh_vertices_to_transect_3D, &
    map_from_mesh_triangles_to_transect_2D, map_from_mesh_triangles_to_transect_3D

contains

  ! From a mesh to a transect
  subroutine map_from_mesh_vertices_to_transect_2D( mesh, transect, d_mesh_partial, d_transect_partial, method)
    ! Map a 2-D data field from the vertices of a mesh to a transect

    ! In/output variables
    type(type_mesh),            intent(in   ) :: mesh
    type(type_transect),        intent(in   ) :: transect
    real(dp), dimension(:    ), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:    ), intent(  out) :: d_transect_partial
    character(len=*), optional, intent(in   ) :: method

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_from_mesh_vertices_to_transect_2D'
    integer                        :: mi, mi_valid
    logical                        :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == trim( mesh%name) // '_vertices' .and. &
          Atlas( mi)%name_dst == 'transect_' // trim( transect%name)) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          select case (method)
          case default
            call crash('invalid method "' // trim( method) // '"')
          case ('trilin')
            call create_map_from_mesh_vertices_to_transect_trilin( mesh, transect, Atlas( mi))
          case ('nearest_neighbour')
            call create_map_from_mesh_vertices_to_transect_nearest_neighbour( mesh, transect, Atlas( mi))
          end select
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_vertices_to_transect_2D( mesh, transect, &
      Atlas( mi), d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_vertices_to_transect_2D

  subroutine map_from_mesh_vertices_to_transect_3D( mesh, transect, d_mesh_partial, d_transect_partial, method)
    ! Map a 3-D data field from the vertices of a mesh to a transect

    ! In/output variables
    type(type_mesh),            intent(in   ) :: mesh
    type(type_transect),        intent(in   ) :: transect
    real(dp), dimension(:,:  ), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:,:  ), intent(  out) :: d_transect_partial
    character(len=*), optional, intent(in   ) :: method

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_from_mesh_vertices_to_transect_3D'
    integer                        :: mi, mi_valid
    logical                        :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == trim( mesh%name) // '_vertices' .and. &
          Atlas( mi)%name_dst == 'transect_' // trim( transect%name)) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_mesh_vertices_to_transect_trilin( mesh, transect, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_vertices_to_transect_3D( mesh, transect, &
      Atlas( mi), d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_vertices_to_transect_3D

  subroutine map_from_mesh_triangles_to_transect_2D( mesh, transect, d_mesh_partial, d_transect_partial, method)
    ! Map a 2-D data field from the triangles of a mesh to a transect

    ! In/output variables
    type(type_mesh),            intent(in   ) :: mesh
    type(type_transect),        intent(in   ) :: transect
    real(dp), dimension(:    ), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:    ), intent(  out) :: d_transect_partial
    character(len=*), optional, intent(in   ) :: method

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_from_mesh_triangles_to_transect_2D'
    integer                        :: mi, mi_valid
    logical                        :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == trim( mesh%name) // '_triangles' .and. &
          Atlas( mi)%name_dst == 'transect_' // trim( transect%name)) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_mesh_triangles_to_transect( mesh, transect, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_triangles_to_transect_2D( mesh, transect, &
      Atlas( mi), d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_triangles_to_transect_2D

  subroutine map_from_mesh_triangles_to_transect_3D( mesh, transect, d_mesh_partial, d_transect_partial, method)
    ! Map a 2-D data field from the triangles of a mesh to a transect

    ! In/output variables
    type(type_mesh),            intent(in   ) :: mesh
    type(type_transect),        intent(in   ) :: transect
    real(dp), dimension(:,:  ), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:,:  ), intent(  out) :: d_transect_partial
    character(len=*), optional, intent(in   ) :: method

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_from_mesh_triangles_to_transect_3D'
    integer                        :: mi, mi_valid
    logical                        :: found_map, found_empty_page

    ! Add routine to path
    call init_routine( routine_name)

    ! Browse the Atlas to see if an appropriate mapping object already exists.
    found_map = .false.
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == trim( mesh%name) // '_triangles' .and. &
          Atlas( mi)%name_dst == 'transect_' // trim( transect%name)) then
        ! if so specified, look for a mapping object with the correct method
        if (present( method)) then
          if (Atlas( mi)%method /= method) cycle
        end if
        found_map = .true.
        mi_valid  = mi
        exit
      end if
    end do

    ! if no appropriate mapping object could be found, create one.
    if (.not. found_map) then
      found_empty_page = .false.
      do mi = 1, size( Atlas,1)
        if (.not. Atlas( mi)%is_in_use) then
          found_empty_page = .true.
          call create_map_from_mesh_triangles_to_transect( mesh, transect, Atlas( mi))
          mi_valid = mi
          exit
        end if
      end do
      ! Safety
      if (.not. found_empty_page) call crash('No more room in Atlas - assign more memory!')
    end if

    ! Apply the appropriate mapping object
    call apply_map_mesh_triangles_to_transect_3D( mesh, transect, &
      Atlas( mi), d_mesh_partial, d_transect_partial)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_mesh_triangles_to_transect_3D

  subroutine create_map_from_mesh_vertices_to_transect_trilin( mesh, transect, map)

    ! In/output variables
    type(type_mesh),     intent(in   ) :: mesh
    type(type_transect), intent(in   ) :: transect
    type(type_map),      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_map_from_mesh_vertices_to_transect_trilin'
    integer                         :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp) :: M_CSR
    integer                         :: i
    real(dp), dimension(2)          :: p
    integer                         :: ti, via, vib, vic
    real(dp), dimension(2)          :: pa, pb, pc
    real(dp)                        :: Atri_abp, Atri_bcp, Atri_cap, Atri_abc, wa, wb, wc
    integer                         :: cola, colb, colc

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = trim( mesh%name) // '_vertices'
    map%name_dst  = 'transect_' // trim( transect%name)
    map%method    = 'trilin'

    ! Matrix size
    nrows           = transect%nV       ! to
    nrows_loc       = transect%nV_loc
    ncols           = mesh%nV           ! from
    ncols_loc       = mesh%nV_loc
    nnz_per_row_max = 3
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ti = 1
    do i = transect%vi1, transect%vi2

      p = transect%V( i,:)
      call find_containing_triangle( mesh, p, ti)

      ! Calculate the trilinear interpolation weights
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      pa  = mesh%V( via,:)
      pb  = mesh%V( vib,:)
      pc  = mesh%V( vic,:)

      Atri_abp = triangle_area( pa, pb, p)
      Atri_bcp = triangle_area( pb, pc, p)
      Atri_cap = triangle_area( pc, pa, p)
      Atri_abc = Atri_abp + Atri_bcp + Atri_cap

      wa = Atri_bcp / Atri_abc
      wb = Atri_cap / Atri_abc
      wc = Atri_abp / Atri_abc

      ! Matrix columns corresponding to these three vertices
      cola = mesh%vi2n( via)
      colb = mesh%vi2n( vib)
      colc = mesh%vi2n( vic)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, i, cola, wa)
      call add_entry_CSR_dist( M_CSR, i, colb, wb)
      call add_entry_CSR_dist( M_CSR, i, colc, wc)

    end do

    call finalise_matrix_CSR_dist( M_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_vertices_to_transect_trilin

  subroutine create_map_from_mesh_vertices_to_transect_nearest_neighbour( mesh, transect, map)

    ! In/output variables
    type(type_mesh),     intent(in   ) :: mesh
    type(type_transect), intent(in   ) :: transect
    type(type_map),      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_map_from_mesh_vertices_to_transect_nearest_neighbour'
    integer                         :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp) :: M_CSR
    integer                         :: i
    real(dp), dimension(2)          :: p
    integer                         :: vi, col

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = trim( mesh%name) // '_vertices'
    map%name_dst  = 'transect_' // trim( transect%name)
    map%method    = 'nearest_neighbour'

    ! Matrix size
    nrows           = transect%nV       ! to
    nrows_loc       = transect%nV_loc
    ncols           = mesh%nV           ! from
    ncols_loc       = mesh%nV_loc
    nnz_per_row_max = 3
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    vi = 1
    do i = transect%vi1, transect%vi2

      p = transect%V( i,:)
      call find_containing_vertex( mesh, p, vi)

      col = mesh%vi2n( vi)

      ! Add to the matrix
      call add_entry_CSR_dist( M_CSR, i, col, 1._dp)

    end do

    call finalise_matrix_CSR_dist( M_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_vertices_to_transect_nearest_neighbour

  subroutine create_map_from_mesh_triangles_to_transect( mesh, transect, map)

    ! In/output variables
    type(type_mesh),     intent(in   ) :: mesh
    type(type_transect), intent(in   ) :: transect
    type(type_map),      intent(inout) :: map

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_map_from_mesh_triangles_to_transect'
    integer                          :: ncols, nrows, nrows_loc, ncols_loc, nnz_per_row_max, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)  :: M_CSR
    integer                          :: i
    real(dp), dimension(2)           :: p
    integer                          :: vi, iti, ti, iti_nearest
    real(dp), dimension(mesh%nC_mem) :: ww
    real(dp)                         :: dist, dist_min

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = trim( mesh%name) // '_triangles'
    map%name_dst  = 'transect_' // trim( transect%name)
    map%method    = 'trilin'

    ! Matrix size
    nrows           = transect%nV       ! to
    nrows_loc       = transect%nV_loc
    ncols           = mesh%nTri         ! from
    ncols_loc       = mesh%nTri_loc
    nnz_per_row_max = mesh%nC_mem
    nnz_est_proc    = nnz_per_row_max * nrows_loc

    call allocate_matrix_CSR_dist( M_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! For all mesh_dst vertices, find the mesh_src vertex containing them
    vi = 1
    do i = transect%vi1, transect%vi2

      p = transect%V( i,:)
      call find_containing_vertex( mesh, p, vi)

      ww          = 0._dp
      dist_min    = huge( dist_min)
      iti_nearest = 0

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)

        dist = norm2( mesh%Tricc( ti,:) - p)
        ww( iti) = 1._dp / dist**2

        if (dist < dist_min) then
          dist_min    = dist
          iti_nearest = ti
        end if

      end do

      if (dist_min < mesh%tol_dist) then
        ! p lies on a triangle circumcentre
        ww = 0._dp
        ww( iti_nearest) = 1._dp
      else
        ww = ww / sum( ww)
      end if

      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        call add_entry_CSR_dist( M_CSR, i, ti, ww( iti))
      end do

    end do

    call finalise_matrix_CSR_dist( M_CSR)

    ! Convert matrices from Fortran to PETSc types
    call mat_CSR2petsc( M_CSR, map%M)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_mesh_triangles_to_transect

end module remapping_transects
