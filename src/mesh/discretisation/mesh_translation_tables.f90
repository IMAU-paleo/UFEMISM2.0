module mesh_translation_tables

  ! Calculate translation tables relating grid points to matrix rows.

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh

  implicit none

contains

  subroutine calc_field_to_vector_form_translation_tables( mesh)
    ! Calculate grid-cell-to-matrix-row translation tables

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'calc_field_to_vector_form_translation_tables'
    integer                       :: nz,vi,ti,ei,k,ks,uv,n

    ! Add routine to path
    call init_routine( routine_name)

    nz = mesh%nz

    ! Grid sizes
    mesh%nna     = mesh%nV
    mesh%nnauv   = mesh%nV            * 2
    mesh%nnak    = mesh%nV   *  nz
    mesh%nnakuv  = mesh%nV   *  nz    * 2
    mesh%nnaks   = mesh%nV   * (nz-1)
    mesh%nnaksuv = mesh%nV   * (nz-1) * 2

    mesh%nnb     = mesh%nTri
    mesh%nnbuv   = mesh%nTri          * 2
    mesh%nnbk    = mesh%nTri *  nz
    mesh%nnbkuv  = mesh%nTri *  nz    * 2
    mesh%nnbks   = mesh%nTri * (nz-1)
    mesh%nnbksuv = mesh%nTri * (nz-1) * 2

    mesh%nnc     = mesh%nE
    mesh%nncuv   = mesh%nE            * 2
    mesh%nnck    = mesh%nE   *  nz
    mesh%nnckuv  = mesh%nE   *  nz    * 2
    mesh%nncks   = mesh%nE   * (nz-1)
    mesh%nncksuv = mesh%nE   * (nz-1) * 2

    ! allocate shared memory
    allocate( mesh%n2vi(     mesh%nna       ), source = 0)
    allocate( mesh%n2viuv(   mesh%nnauv  , 2), source = 0)
    allocate( mesh%n2vik(    mesh%nnak   , 2), source = 0)
    allocate( mesh%n2vikuv(  mesh%nnakuv , 3), source = 0)
    allocate( mesh%n2viks(   mesh%nnaks  , 2), source = 0)
    allocate( mesh%n2viksuv( mesh%nnaksuv, 3), source = 0)

    allocate( mesh%n2ti(     mesh%nnb       ), source = 0)
    allocate( mesh%n2tiuv(   mesh%nnbuv  , 2), source = 0)
    allocate( mesh%n2tik(    mesh%nnbk   , 2), source = 0)
    allocate( mesh%n2tikuv(  mesh%nnbkuv , 3), source = 0)
    allocate( mesh%n2tiks(   mesh%nnbks  , 2), source = 0)
    allocate( mesh%n2tiksuv( mesh%nnbksuv, 3), source = 0)

    allocate( mesh%n2ei(     mesh%nnc       ), source = 0)
    allocate( mesh%n2eiuv(   mesh%nncuv  , 2), source = 0)
    allocate( mesh%n2eik(    mesh%nnck   , 2), source = 0)
    allocate( mesh%n2eikuv(  mesh%nnckuv , 3), source = 0)
    allocate( mesh%n2eiks(   mesh%nncks  , 2), source = 0)
    allocate( mesh%n2eiksuv( mesh%nncksuv, 3), source = 0)

    allocate( mesh%vi2n(     mesh%nV           ), source = 0)
    allocate( mesh%viuv2n(   mesh%nV        , 2), source = 0)
    allocate( mesh%vik2n(    mesh%nV  , nz     ), source = 0)
    allocate( mesh%vikuv2n(  mesh%nV  , nz  , 2), source = 0)
    allocate( mesh%viks2n(   mesh%nV  , nz-1   ), source = 0)
    allocate( mesh%viksuv2n( mesh%nV  , nz-1, 2), source = 0)

    allocate( mesh%ti2n(     mesh%nTri         ), source = 0)
    allocate( mesh%tiuv2n(   mesh%nTri      , 2), source = 0)
    allocate( mesh%tik2n(    mesh%nTri, nz     ), source = 0)
    allocate( mesh%tikuv2n(  mesh%nTri, nz  , 2), source = 0)
    allocate( mesh%tiks2n(   mesh%nTri, nz-1   ), source = 0)
    allocate( mesh%tiksuv2n( mesh%nTri, nz-1, 2), source = 0)

    allocate( mesh%ei2n(     mesh%nE           ), source = 0)
    allocate( mesh%eiuv2n(   mesh%nE        , 2), source = 0)
    allocate( mesh%eik2n(    mesh%nE  , nz     ), source = 0)
    allocate( mesh%eikuv2n(  mesh%nE  , nz  , 2), source = 0)
    allocate( mesh%eiks2n(   mesh%nE  , nz-1   ), source = 0)
    allocate( mesh%eiksuv2n( mesh%nE  , nz-1, 2), source = 0)

    ! == a-grid (vertices)

      ! == 2-D

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        n = n+1
        mesh%vi2n( vi) = n
        mesh%n2vi( n ) = vi
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do uv = 1, 2
          n = n+1
          mesh%viuv2n( vi,uv) = n
          mesh%n2viuv( n,1) = vi
          mesh%n2viuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        do k = 1, nz
          n = n+1
          mesh%vik2n( vi,k) = n
          mesh%n2vik( n,1) = vi
          mesh%n2vik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%vikuv2n( vi,k,uv) = n
            mesh%n2vikuv( n,1) = vi
            mesh%n2vikuv( n,2) = k
            mesh%n2vikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        do ks = 1, nz-1
          n = n+1
          mesh%viks2n( vi,ks) = n
          mesh%n2viks( n,1) = vi
          mesh%n2viks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%viksuv2n( vi,ks,uv) = n
            mesh%n2viksuv( n,1) = vi
            mesh%n2viksuv( n,2) = ks
            mesh%n2viksuv( n,3) = uv
          end do
        end do
      end do

    ! == b-grid (triangles)

      ! == 2-D

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        n = n+1
        mesh%ti2n( ti) = n
        mesh%n2ti( n ) = ti
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do uv = 1, 2
          n = n+1
          mesh%tiuv2n( ti,uv) = n
          mesh%n2tiuv( n,1) = ti
          mesh%n2tiuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        do k = 1, nz
          n = n+1
          mesh%tik2n( ti,k) = n
          mesh%n2tik( n,1) = ti
          mesh%n2tik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%tikuv2n( ti,k,uv) = n
            mesh%n2tikuv( n,1) = ti
            mesh%n2tikuv( n,2) = k
            mesh%n2tikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        do ks = 1, nz-1
          n = n+1
          mesh%tiks2n( ti,ks) = n
          mesh%n2tiks( n,1) = ti
          mesh%n2tiks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%tiksuv2n( ti,ks,uv) = n
            mesh%n2tiksuv( n,1) = ti
            mesh%n2tiksuv( n,2) = ks
            mesh%n2tiksuv( n,3) = uv
          end do
        end do
      end do

    ! == c-grid (edges)

      ! == 2-D

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        n = n+1
        mesh%ei2n( ei) = n
        mesh%n2ei( n ) = ei
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do uv = 1, 2
          n = n+1
          mesh%eiuv2n( ei,uv) = n
          mesh%n2eiuv( n,1) = ei
          mesh%n2eiuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        do k = 1, nz
          n = n+1
          mesh%eik2n( ei,k) = n
          mesh%n2eik( n,1) = ei
          mesh%n2eik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%eikuv2n( ei,k,uv) = n
            mesh%n2eikuv( n,1) = ei
            mesh%n2eikuv( n,2) = k
            mesh%n2eikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        do ks = 1, nz-1
          n = n+1
          mesh%eiks2n( ei,ks) = n
          mesh%n2eiks( n,1) = ei
          mesh%n2eiks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%eiksuv2n( ei,ks,uv) = n
            mesh%n2eiksuv( n,1) = ei
            mesh%n2eiksuv( n,2) = ks
            mesh%n2eiksuv( n,3) = uv
          end do
        end do
      end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_field_to_vector_form_translation_tables

end module mesh_translation_tables
