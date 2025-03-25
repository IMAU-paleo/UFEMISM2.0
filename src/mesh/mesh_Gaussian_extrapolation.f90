module mesh_Gaussian_extrapolation

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_MAX, &
    MPI_SUM, MPI_COMM_WORLD, MPI_WIN
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    gather_dist_shared_to_all

  implicit none

  private

  public :: extrapolate_Gaussian

contains

  subroutine extrapolate_Gaussian( mesh, mask_partial, d_partial, sigma, d_is_hybrid)
    !< Extrapolate the data field d into the area designated by the mask,
    !<using Gaussian extrapolation of sigma

    ! Note about the mask:
    !    2 = data provided
    !    1 = no data provided, fill allowed
    !    0 = no fill allowed
    !
    ! (so basically this routine extrapolates data from the area
    !  where mask == 2 into the area where mask == 1)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    integer,  dimension(:), intent(in   ) :: mask_partial
    real(dp), dimension(:), intent(inout) :: d_partial
    real(dp),               intent(in   ) :: sigma
    logical, optional,      intent(in   ) :: d_is_hybrid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'extrapolate_Gaussian'
    logical                        :: d_is_hybrid_

    ! Add routine to path
    call init_routine( routine_name)

    if (present( d_is_hybrid)) then
      d_is_hybrid_ = d_is_hybrid
    else
      d_is_hybrid_ = .false.
    end if

    if (d_is_hybrid_) then
      call extrapolate_Gaussian_hybrid( mesh, mask_partial, d_partial, sigma)
    else
      call extrapolate_Gaussian_dist( mesh, mask_partial, d_partial, sigma)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_Gaussian

  subroutine extrapolate_Gaussian_dist( mesh, mask_partial, d_partial, sigma)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_partial
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: d_partial
    real(dp),                               intent(in   ) :: sigma

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'extrapolate_Gaussian_dist'
    integer                                :: ierr
    integer,  dimension(mesh%vi1:mesh%vi2) :: mask_local
    integer,  dimension(mesh%nV)           :: mask_tot
    real(dp), dimension(mesh%nV)           :: d_tot
    integer                                :: it_floodfill
    integer                                :: vi
    logical,  dimension(mesh%vi1:mesh%vi2) :: do_fill_now
    integer                                :: n_do_fill_now
    LOGICAL                                :: has_filled_neighbour
    integer                                :: ci,vj
    integer,  dimension(mesh%nV)           :: map_neighbourhood_of_vi
    integer,  dimension(mesh%nV)           :: stack_front_around_vi
    integer                                :: stackN_front_around_vi
    integer,  dimension(mesh%nV)           :: stack_neighbourhood_of_vi
    integer                                :: stackN_neighbourhood_of_vi
    integer                                :: it_floodfill2
    integer                                :: cj,vk
    integer                                :: i
    real(dp)                               :: wj, w_sum, d_sum, d_av_of_neighbourhood

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    map_neighbourhood_of_vi    = 0
    stack_front_around_vi      = 0
    stackN_front_around_vi     = 0
    stack_neighbourhood_of_vi  = 0
    stackN_neighbourhood_of_vi = 0

    ! Copy mask to local, changeable array
    mask_local = mask_partial

    ! Gather complete fill mask and data to all processes
    call gather_to_all( mask_local, mask_tot)
    call gather_to_all(  d_partial , d_tot   )

    ! == Flood-fill iteration
    ! =======================

    it_floodfill = 0

    iterate_floodfill: do while (.true.)

      ! Safety
      it_floodfill = it_floodfill + 1
      if (it_floodfill > mesh%nV) call crash('main flood-fill iteration got stuck!')

      ! == Mark all vertices that should be filled now
      !    (defined as those that are allowed to be filled,
      !    and are next to at least one filled vertex).
      ! ===============================================

      n_do_fill_now = 0
      do_fill_now = .false.

      do vi = mesh%vi1, mesh%vi2
        if (mask_tot( vi) == 1) then
          ! Vertex vi is allowed to be filled

          has_filled_neighbour = .false.
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mask_tot( vj) == 2) then
              has_filled_neighbour = .true.
              exit
            end if
          end do

          if (has_filled_neighbour) then
            ! Vertex vi is allowed to be filled and has at least one filled neighbour,
            ! so it should be filled in this flood-fill iteration
            do_fill_now( vi) = .true.
            n_do_fill_now = n_do_fill_now + 1
          end if

        end if
      end do

      ! If no vertices can be filled anymore, end the flood-fill iteration
      call MPI_ALLREDUCE( MPI_IN_PLACE, n_do_fill_now, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (n_do_fill_now == 0) exit iterate_floodfill

      ! == Fill all vertices that can be filled now
      ! ===========================================

      do vi = mesh%vi1, mesh%vi2
        if (do_fill_now( vi)) then
          ! Fill vertex vi

          ! == Find all filled vertices within 3*sigma of vertex vi
          ! =======================================================

          ! Initialise the front with all filled neighbours of vi
          stackN_front_around_vi = 0
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mask_tot( vj) == 2) then
              stackN_front_around_vi = stackN_front_around_vi + 1
              stack_front_around_vi( stackN_front_around_vi) = vj
              map_neighbourhood_of_vi( vj) = 1
            end if
          end do

          ! Expand the front outward, flood-fill style
          it_floodfill2 = 0
          iterate_floodfill_around_vi: do while (stackN_front_around_vi > 0)

            ! Safety
            it_floodfill2 = it_floodfill2 + 1
            if (it_floodfill2 > mesh%nV) call crash('secondary flood-fill iteration got stuck!')

            ! Take the last vertex vj from the front stack
            vj = stack_front_around_vi( stackN_front_around_vi)
            stackN_front_around_vi = stackN_front_around_vi - 1

            ! Add vj to the neighbourhood list
            stackN_neighbourhood_of_vi = stackN_neighbourhood_of_vi + 1
            stack_neighbourhood_of_vi( stackN_neighbourhood_of_vi) = vj
            map_neighbourhood_of_vi( vj) = 2

            ! Add all remaining filled neighbours of vj to the front stack
            do cj = 1, mesh%nC( vj)
              vk = mesh%C( vj,cj)
              if (map_neighbourhood_of_vi( vk) == 0 .and. mask_tot( vk) == 2 .and. &
                  norm2( mesh%V( vi,:) - mesh%V( vk,:)) < (3._dp * sigma)) then
                ! Vertex vk is not yet marked as part of the neighbourhood of vi, not
                ! yet marked as part of the front stack, is filled, and lies within
                ! 3*sigma of vertex vi. Add it to the front stack.
                stackN_front_around_vi = stackN_front_around_vi + 1
                stack_front_around_vi( stackN_front_around_vi) = vk
                map_neighbourhood_of_vi( vk) = 1
              end if
            end do

          end do iterate_floodfill_around_vi

          ! Safety
          if (stackN_neighbourhood_of_vi == 0) call crash('couldnt find neighbourhood of vi!')

          ! == Extrapolate data to vi
          ! =========================

          ! Calculate Gaussian distance-weighted average of d over the filled neighbourhood of vi
          w_sum = 0._dp
          d_sum = 0._dp
          do i = 1, stackN_neighbourhood_of_vi
            vj = stack_neighbourhood_of_vi( i)
            wj = EXP( -0.5_dp * (NORM2( mesh%V( vj,:) - mesh%V( vi,:)) / sigma)**2)
            w_sum = w_sum + wj
            d_sum = d_sum + wj * d_tot( vj)
          end do
          d_av_of_neighbourhood = d_sum / w_sum

          ! Fill into data field
          d_partial( vi) = d_av_of_neighbourhood

          ! == Clean up lists for the neighbourhood of vi
          ! =============================================

          do i = 1, stackN_neighbourhood_of_vi
            vj = stack_neighbourhood_of_vi( i)
            map_neighbourhood_of_vi( vj) = 0
          end do
          stackN_neighbourhood_of_vi = 0
          stackN_front_around_vi     = 0

        end if
      end do

      ! Mark newly filled vertices in the mask, so they
      ! can contribute to the next extrapolation iteration
      ! ==================================================

      do vi = mesh%vi1, mesh%vi2
        if (do_fill_now( vi)) then
          mask_local( vi) = 2
        end if
      end do

      ! Exchange newly filled mask and data between the processes
      ! =========================================================

      call gather_to_all( mask_local, mask_tot)
      call gather_to_all(  d_partial , d_tot   )

    end do iterate_floodfill

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_Gaussian_dist

  subroutine extrapolate_Gaussian_hybrid( mesh, mask_partial, d_partial, sigma)

    ! In/output variables:
    type(type_mesh),                                  intent(in   ) :: mesh
    integer,  dimension(mesh%vi1_node:mesh%vi2_node), intent(in   ) :: mask_partial
    real(dp), dimension(mesh%vi1_node:mesh%vi2_node), intent(inout) :: d_partial
    real(dp),                                         intent(in   ) :: sigma

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'extrapolate_Gaussian_hybrid'
    integer                                :: ierr
    integer,  dimension(:), pointer        :: mask_local => null()
    integer,  dimension(:), pointer        :: mask_tot => null()
    real(dp), dimension(:), pointer        :: d_tot => null()
    type(MPI_WIN)                          :: wmask_local, wmask_tot, wd_tot
    integer                                :: it_floodfill
    integer                                :: vi
    logical,  dimension(:), pointer        :: do_fill_now => null()
    type(MPI_WIN)                          :: wdo_fill_now
    integer                                :: n_do_fill_now
    LOGICAL                                :: has_filled_neighbour
    integer                                :: ci,vj
    integer,  dimension(mesh%nV)           :: map_neighbourhood_of_vi
    integer,  dimension(mesh%nV)           :: stack_front_around_vi
    integer                                :: stackN_front_around_vi
    integer,  dimension(mesh%nV)           :: stack_neighbourhood_of_vi
    integer                                :: stackN_neighbourhood_of_vi
    integer                                :: it_floodfill2
    integer                                :: cj,vk
    integer                                :: i
    real(dp)                               :: wj, w_sum, d_sum, d_av_of_neighbourhood

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( mask_local , wmask_local , mesh%nV_node)
    call allocate_dist_shared( mask_tot   , wmask_tot   , mesh%nV)
    call allocate_dist_shared( d_tot      , wd_tot      , mesh%nV)
    call allocate_dist_shared( do_fill_now, wdo_fill_now, mesh%nV_node)
    mask_local ( mesh%vi1_node:mesh%vi2_node) => mask_local
    do_fill_now( mesh%vi1_node:mesh%vi2_node) => do_fill_now

    ! Initialise
    map_neighbourhood_of_vi    = 0
    stack_front_around_vi      = 0
    stackN_front_around_vi     = 0
    stack_neighbourhood_of_vi  = 0
    stackN_neighbourhood_of_vi = 0

    ! Copy mask to local, changeable array
    mask_local = mask_partial

    ! Gather complete fill mask and data to all processes
    call gather_dist_shared_to_all( mask_local, mask_tot)
    call gather_dist_shared_to_all(  d_partial , d_tot   )

    ! == Flood-fill iteration
    ! =======================

    it_floodfill = 0

    iterate_floodfill: do while (.true.)

      ! Safety
      it_floodfill = it_floodfill + 1
      if (it_floodfill > mesh%nV) call crash('main flood-fill iteration got stuck!')

      ! == Mark all vertices that should be filled now
      !    (defined as those that are allowed to be filled,
      !    and are next to at least one filled vertex).
      ! ===============================================

      n_do_fill_now = 0
      do_fill_now = .false.

      do vi = mesh%vi1, mesh%vi2
        if (mask_tot( vi) == 1) then
          ! Vertex vi is allowed to be filled

          has_filled_neighbour = .false.
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mask_tot( vj) == 2) then
              has_filled_neighbour = .true.
              exit
            end if
          end do

          if (has_filled_neighbour) then
            ! Vertex vi is allowed to be filled and has at least one filled neighbour,
            ! so it should be filled in this flood-fill iteration
            do_fill_now( vi) = .true.
            n_do_fill_now = n_do_fill_now + 1
          end if

        end if
      end do

      ! If no vertices can be filled anymore, end the flood-fill iteration
      call MPI_ALLREDUCE( MPI_IN_PLACE, n_do_fill_now, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (n_do_fill_now == 0) exit iterate_floodfill

      ! == Fill all vertices that can be filled now
      ! ===========================================

      do vi = mesh%vi1, mesh%vi2
        if (do_fill_now( vi)) then
          ! Fill vertex vi

          ! == Find all filled vertices within 3*sigma of vertex vi
          ! =======================================================

          ! Initialise the front with all filled neighbours of vi
          stackN_front_around_vi = 0
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mask_tot( vj) == 2) then
              stackN_front_around_vi = stackN_front_around_vi + 1
              stack_front_around_vi( stackN_front_around_vi) = vj
              map_neighbourhood_of_vi( vj) = 1
            end if
          end do

          ! Expand the front outward, flood-fill style
          it_floodfill2 = 0
          iterate_floodfill_around_vi: do while (stackN_front_around_vi > 0)

            ! Safety
            it_floodfill2 = it_floodfill2 + 1
            if (it_floodfill2 > mesh%nV) call crash('secondary flood-fill iteration got stuck!')

            ! Take the last vertex vj from the front stack
            vj = stack_front_around_vi( stackN_front_around_vi)
            stackN_front_around_vi = stackN_front_around_vi - 1

            ! Add vj to the neighbourhood list
            stackN_neighbourhood_of_vi = stackN_neighbourhood_of_vi + 1
            stack_neighbourhood_of_vi( stackN_neighbourhood_of_vi) = vj
            map_neighbourhood_of_vi( vj) = 2

            ! Add all remaining filled neighbours of vj to the front stack
            do cj = 1, mesh%nC( vj)
              vk = mesh%C( vj,cj)
              if (map_neighbourhood_of_vi( vk) == 0 .and. mask_tot( vk) == 2 .and. &
                  norm2( mesh%V( vi,:) - mesh%V( vk,:)) < (3._dp * sigma)) then
                ! Vertex vk is not yet marked as part of the neighbourhood of vi, not
                ! yet marked as part of the front stack, is filled, and lies within
                ! 3*sigma of vertex vi. Add it to the front stack.
                stackN_front_around_vi = stackN_front_around_vi + 1
                stack_front_around_vi( stackN_front_around_vi) = vk
                map_neighbourhood_of_vi( vk) = 1
              end if
            end do

          end do iterate_floodfill_around_vi

          ! Safety
          if (stackN_neighbourhood_of_vi == 0) call crash('couldnt find neighbourhood of vi!')

          ! == Extrapolate data to vi
          ! =========================

          ! Calculate Gaussian distance-weighted average of d over the filled neighbourhood of vi
          w_sum = 0._dp
          d_sum = 0._dp
          do i = 1, stackN_neighbourhood_of_vi
            vj = stack_neighbourhood_of_vi( i)
            wj = EXP( -0.5_dp * (NORM2( mesh%V( vj,:) - mesh%V( vi,:)) / sigma)**2)
            w_sum = w_sum + wj
            d_sum = d_sum + wj * d_tot( vj)
          end do
          d_av_of_neighbourhood = d_sum / w_sum

          ! Fill into data field
          d_partial( vi) = d_av_of_neighbourhood

          ! == Clean up lists for the neighbourhood of vi
          ! =============================================

          do i = 1, stackN_neighbourhood_of_vi
            vj = stack_neighbourhood_of_vi( i)
            map_neighbourhood_of_vi( vj) = 0
          end do
          stackN_neighbourhood_of_vi = 0
          stackN_front_around_vi     = 0

        end if
      end do

      ! Mark newly filled vertices in the mask, so they
      ! can contribute to the next extrapolation iteration
      ! ==================================================

      do vi = mesh%vi1, mesh%vi2
        if (do_fill_now( vi)) then
          mask_local( vi) = 2
        end if
      end do

      ! Exchange newly filled mask and data between the processes
      ! =========================================================

      call gather_dist_shared_to_all( mask_local, mask_tot)
      call gather_dist_shared_to_all(  d_partial , d_tot   )

    end do iterate_floodfill

    ! Clean up after yourself
    call deallocate_dist_shared( mask_local , wmask_local )
    call deallocate_dist_shared( mask_tot   , wmask_tot   )
    call deallocate_dist_shared( d_tot      , wd_tot      )
    call deallocate_dist_shared( do_fill_now, wdo_fill_now)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_Gaussian_hybrid

end module mesh_Gaussian_extrapolation