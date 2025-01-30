MODULE GIA_ELRA

  ! ELRA module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE GIA_model_types                                        , ONLY: type_GIA_model
  USE region_types                                           , ONLY: type_model_region
  USE grid_basic                                             , ONLY: setup_square_grid
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_master, distribute_gridded_data_from_master
  USE grid_types                                             , ONLY: type_grid
  USE reference_geometry_types                               , ONLY: type_reference_geometry
  use ice_geometry_basics, only: is_floating
  use remapping_main, only: map_from_mesh_to_xy_grid_2D, map_from_xy_grid_to_mesh_2D



  IMPLICIT NONE

CONTAINS


  ! The ELRA GIA model
  SUBROUTINE run_ELRA_model( region)
    ! Use the ELRA model to update bedrock elevation. Once every (dt_bedrock_ELRA) years,
    ! update deformation rates. In all other time steps, just incrementally add deformation.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ELRA_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! I will comment this bcs it will be controled by GIA_main
    !IF (region%do_ELRA) THEN
    CALL calculate_ELRA_bedrock_deformation_rate( region%mesh, region%GIA%grid, region%ice, region%GIA)
    !  region%t_last_ELRA = region%time
    !END IF
    
    ! Update bedrock with last calculated deformation rate
    ! this is done in GIA_main now, commented
    !DO vi = region%mesh%vi1, region%mesh%vi2
    !  region%ice%Hb(  vi) = region%ice%Hb( vi) + region%ice%dHb_dt( vi) * region%dt
    !  region%ice%dHb( vi) = region%ice%Hb( vi) - region%refgeo_GIAeq%Hb( vi) 
    !END DO
    !CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ELRA_model
  SUBROUTINE calculate_ELRA_bedrock_deformation_rate( mesh, grid, ice, GIA)
    ! Use the ELRA model to update bedrock deformation rates.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    ! maybe this one should be only IN and not INOUT? the changes in ice will be done in GIA_main
    ! CHECK if write the result of ELRA in GIA%Hb_next is right, if it is right, change ice to IN
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_GIA_model),                INTENT(INOUT) :: GIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_ELRA_bedrock_deformation_rate'
    INTEGER                                            :: vi,i,j,n,k,l,ii,jj
    REAL(dp)                                           :: Lr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp

    ! Calculate the absolute and relative surface loads on the mesh
    
    DO vi = mesh%vi1, mesh%vi2

! I am not really sure if this is needed, I think it is so I will keep it however, ice could be just IN not necessary INOUT
      ! Absolute surface load
      IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN
        GIA%surface_load_mesh( vi) = (ice%SL( vi) - ice%Hb( vi)) * grid%dx**2 * seawater_density
      ELSEIF (ice%Hi( vi) > 0._dp) THEN
        GIA%surface_load_mesh( vi) =  ice%Hi( vi) * grid%dx**2 * ice_density
      ELSE
        GIA%surface_load_mesh( vi) = 0._dp
      END IF


      ! Relative surface load
      !asegurarse que surface load PD mesh haya sido calculado
      GIA%relative_surface_load_mesh( vi) = GIA%surface_load_mesh( vi) - GIA%surface_load_PD_mesh( vi)

    END DO
    CALL sync

    ! Map relative surface load to the GIA grid
    
    !! ---
    ! changed function map_mesh2grid_2D() to map_from_mesh_to_xy_grid_2D
    !! ----
    !! THIS IS FINE, it goes from mesh to xy grid both partial, so the output is (n1:n2)
    CALL map_from_mesh_to_xy_grid_2D( mesh, grid, GIA%relative_surface_load_mesh, GIA%relative_surface_load_grid)
    
    !! HERE I WILL gather data to master goes from partial vec to total 2D
    !! CHECK IF THIS "LIBRARY"HAS TO BE ADDED at the begining
    call gather_gridded_data_to_master( grid, GIA%relative_surface_load_grid, GIA%relative_surface_load_grid_tot)
    
    
    ! GIA%flex_prof_rad added to GIA_model_types.f90
    n = GIA%flex_prof_rad
    
    ! Let the master do the actual work
    if (par%master) then
    
    do i = 1, grid%nx
    do j = 1, grid%ny
      GIA%dHb_eq_grid( i, j) = 0._dp
      do k = -n, n
      do l = -n, n
        ii = max( 1, min( grid%nx, i+k ))
        jj = max( 1, min( grid%ny, j+l ))
        GIA%dHb_eq_grid( i, j) = GIA%dHb_eq_grid( i, j) + &
          (0.5_dp * grav * Lr**2 /(pi * C%ELRA_lithosphere_flex_rigidity) * GIA%relative_surface_load_grid_tot( ii, jj) * GIA%flex_prof_grid( k+n+1,l+n+1))
      end do
      end do
    end do
    end do
          
    end if ! if (par%master) then
    
    ! Map the actual bedrock deformation to the grid
    ! ice%dHb_grid changed to GIA%dHb_grid_partial    
    call map_from_mesh_to_xy_grid_2D( mesh, grid, ice%dHb, GIA%dHb_grid_partial)
    
    ! gather data from all processors to master, from partial grid vec to total 2D grid 
    call gather_gridded_data_to_master( grid, GIA%dHb_grid_partial, GIA%dHb_grid_tot)    

	! Let the master do the actual work
    if (par%master) then

    ! Calculate the bedrock deformation rate from the difference between the current and the equilibrium deformation
    DO i = 1, grid%nx    
    DO j = 1, grid%ny
    ! ice%dHb_dt_grid changed to GIA%dHb_dt_grid
      GIA%dHb_dt_grid( i,j) = (GIA%dHb_eq_grid( i,j) - GIA%dHb_grid_tot( i,j)) / C%ELRA_bedrock_relaxation_time
    END DO
    END DO
    
    write(*,*) 'dHb_dt_grid, nx and ny:', GIA%dHb_dt_grid, grid%nx, grid%ny
!    GIA%relative_surface_load_grid_tot = GIA%dHb_dt_grid
    end if ! if (par%master) then
   
   ! distribute from 2D grid data on master to vector grid data on all processors
   call distribute_gridded_data_from_master( grid, GIA%dHb_dt_grid, GIA%dHb_dt_grid_partial)
!   call distribute_gridded_data_from_master( grid, GIA%relative_surface_load_grid_tot, GIA%relative_surface_load_grid)
!   print*, "before remapping xy grid to mesh GIA%dHb_dt_grid_partial", GIA%dHb_dt_grid_partial
   ! remap from partial grid vec data to mesh model
   ! add the changed dHb directly to GIA%dHb_next, this will be used later on GIA_main
   call map_from_xy_grid_to_mesh_2D( grid, mesh, GIA%dHb_dt_grid_partial, GIA%dHb_next)
!   call map_from_xy_grid_to_mesh_2D( grid, mesh, GIA%relative_surface_load_grid, GIA%relative_surface_load_mesh)
!   GIA%dHb_next=GIA%relative_surface_load_mesh

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_ELRA_bedrock_deformation_rate

  SUBROUTINE initialise_ELRA_model( mesh, grid, GIA, refgeo_GIAeq)
    ! Allocate and initialise the ELRA GIA model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_GIA_model),                INTENT(INOUT) :: GIA    
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ELRA_model'
    INTEGER                                            :: i,j,n,k,l
    REAL(dp)                                           :: Lr, r

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising ELRA GIA model...'

    ! Allocate memory
    ALLOCATE( GIA%surface_load_PD_mesh( mesh%vi1:mesh%vi2))
    ! relative_surface_load_mesh already allocated in GIA_main
    ! relative_surface_load_grid already allocated in GIA_main as     
    ! ALLOCATE( GIA%relative_surface_load_grid( GIA%grid%n1:GIA%grid%n2)) is a vector not 2D grid..
    ALLOCATE( GIA%relative_surface_load_grid_tot( grid%nx, grid%ny))
    ALLOCATE( GIA%surface_load_mesh( mesh%vi1:mesh%vi2))
    ALLOCATE( GIA%dHb_eq_grid( grid%nx, grid%ny))
    ALLOCATE( GIA%dHb_grid_partial( grid%n1:grid%n2))
    ALLOCATE( GIA%dHb_grid_tot( grid%nx, grid%ny))
    ALLOCATE( GIA%dHb_dt_grid( grid%nx, grid%ny))
    ALLOCATE( GIA%dHb_dt_grid_partial( grid%n1:grid%n2))    
    
!    ALLOCATE( GIA%dHb_grid( grid%nx, grid%ny))
!    ALLOCATE( GIA%dHb_dt_grid( grid%nx, grid%ny)) 
!    CALL allocate_shared_dp_1D( mesh%nV,          ice%surface_load_PD_mesh,  ice%wsurface_load_PD_mesh )
!    CALL allocate_shared_dp_1D( mesh%nV,          ice%surface_load_mesh,     ice%wsurface_load_mesh    )
!    CALL allocate_shared_dp_1D( mesh%nV,          ice%surface_load_rel_mesh, ice%wsurface_load_rel_mesh)
!    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%surface_load_rel_grid, ice%wsurface_load_rel_grid)
!    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%dHb_eq_grid,           ice%wdHb_eq_grid          )
!    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%dHb_grid,              ice%wdHb_grid             )
!    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%dHb_dt_grid,           ice%wdHb_dt_grid          )

    ! Fill in the 2D flexural profile (= Kelvin function), with which
    ! a surface load is convoluted to find the surface deformation
    ! ============================================================

    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp

    ! Calculate radius (in number of grid cells) of the flexural profile
    
    !! In the original here it comes an allocate of the flex_prof_rad but this is not done anymore in Ufe2
    !! I will omit it now but maybe I have to add something to make it work
    !CALL allocate_shared_int_0D( ice%flex_prof_rad, ice%wflex_prof_rad)
    
    IF (par%master) THEN
    GIA%flex_prof_rad = MIN( CEILING(grid%dx/2._dp), MAX(1, INT(6._dp * Lr / grid%dx) - 1))
    n = 2 * GIA%flex_prof_rad + 1
    END IF ! IF (par%master) THEN
    CALL sync
    
    ALLOCATE( GIA%flex_prof_grid( n, n))

! CHECK IF THIS IS RIGHT, DO NOT DELETE THE PREVIOUS    
!	ALLOCATE( GIA%flex_prof_grid( n, n))
!	ALLOCATE( GIA%relative_surface_load_ext_grid( grid%nx +n, grid%ny + n))
!    CALL allocate_shared_dp_2D( n,         n,         ice%flex_prof_grid,            ice%wflex_prof_grid           )
!    CALL allocate_shared_dp_2D( grid%nx+n, grid%ny+n, ice%surface_load_rel_ext_grid, ice%wsurface_load_rel_ext_grid)

    ! Calculate flexural profile
    IF (par%master) THEN
      DO i = -GIA%flex_prof_rad, GIA%flex_prof_rad
      DO j = -GIA%flex_prof_rad, GIA%flex_prof_rad
        l = i+GIA%flex_prof_rad+1
        k = j+GIA%flex_prof_rad+1
        r = grid%dx * SQRT( (REAL(i,dp))**2 + (REAL(j,dp))**2)
        GIA%flex_prof_grid( l,k) = kelvin(r / Lr)
      END DO
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Calculate the PD reference load
    ! ===============================

! HERE I ADDED GIA%grid
! should here ice be replaced by GIA?
! CHECK!
!!
    CALL initialise_ELRA_PD_reference_load( mesh, grid, GIA, refgeo_GIAeq)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ELRA_model
  SUBROUTINE initialise_ELRA_PD_reference_load( mesh, grid, GIA, refgeo_GIAeq)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    ! I changed ice for GIA here 
    TYPE(type_GIA_model),                INTENT(INOUT) :: GIA    
!    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ELRA_PD_reference_load'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (par%master) WRITE (0,*) '  Initialising ELRA PD reference load...'

    ! Calculate PD reference load on the mesh
    DO vi = mesh%vi1, mesh%vi2
!            print*, "refgeo Hi(vi)", refgeo_GIAeq%Hi(vi)
      IF (is_floating( refgeo_GIAeq%Hi( vi), refgeo_GIAeq%Hb( vi), 0._dp)) THEN
        GIA%surface_load_PD_mesh( vi) = -refgeo_GIAeq%Hb( vi) * grid%dx**2 * seawater_density
      ELSEIF (refgeo_GIAeq%Hi( vi) > 0._dp) THEN
        GIA%surface_load_PD_mesh( vi) =  refgeo_GIAeq%Hi( vi) * grid%dx**2 * ice_density
      ELSE
        GIA%surface_load_PD_mesh( vi) = 0._dp
      END IF
    END DO
    CALL sync
    
    IF (par%master) WRITE (0,*) '  Finalising ELRA PD reference load...'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ELRA_PD_reference_load
  SUBROUTINE remap_ELRA_model( mesh_old, mesh_new, GIA, refgeo_GIAeq, grid)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    ! changed ice for GIA type and INOUT
    TYPE(type_GIA_model),                INTENT(INOUT) :: GIA
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq
    ! this grid is GIA%grid on the main call
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ELRA_model'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)
! here I Will comment this and also the map type variable of the subroutine
! this was used to prevent compiler warnings for unused variables in previous version of Ufe
!    int_dummy = mesh_old%nV
!    int_dummy = map%int_dummy
	CALL reallocate_bounds( GIA%surface_load_PD_mesh, mesh_new%vi1, mesh_new%vi2)
	CALL reallocate_bounds( GIA%surface_load_mesh, mesh_new%vi1, mesh_new%vi2)
	! this one is also in GIA_main, is needed here?
	!CALL reallocate_bounds( GIA%relative_surface_load_mesh, mesh_new%vi1, mesh_new%vi2) commented is already in GIA_main
	
    !CALL reallocate_shared_dp_1D( mesh_new%nV, ice%surface_load_PD_mesh,  ice%wsurface_load_PD_mesh )
    !CALL reallocate_shared_dp_1D( mesh_new%nV, ice%surface_load_mesh,     ice%wsurface_load_mesh    )
    !CALL reallocate_shared_dp_1D( mesh_new%nV, ice%surface_load_rel_mesh, ice%wsurface_load_rel_mesh)

    ! Recalculate the PD reference load on the GIA grid
    ! I think the call here to ice should be GIA
    CALL initialise_ELRA_PD_reference_load( mesh_new, grid, GIA, refgeo_GIAeq)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ELRA_model

  ! The Kelvin function (just a wrapper for the Zhang & Jin "klvna" subroutine)
  FUNCTION kelvin(x) RESULT(kei)
    ! Evaluate the Kelvin function in a given point at a distance x from the load. Currently this is implemented
    ! using the klvna routine in mklvna.f (see that routine for details).

    IMPLICIT NONE

    ! Input variables:
    ! x is the distance between a considered load and the point in which we want to know the deformation, x is in units Lr
    REAL(dp), INTENT(IN) :: x

    ! Result variables:
    REAL(dp)             :: kei

    ! Local variables:
    REAL(dp)     :: xdp,berdp,beidp,gerdp,geidp,derdp,deidp,herdp,heidp

    xdp = x
    CALL klvna( xdp,berdp,beidp,gerdp,geidp,derdp,deidp,herdp,heidp)
    kei = geidp

  END FUNCTION kelvin
  SUBROUTINE klvna ( x, ber, bei, ger, gei, der, dei, her, hei )

!*****************************************************************************80
!
!! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real (dp) X, the argument.
!
!    Output, real (dp) BER, BEI, GER, GEI, DER, DEI, HER, HEI,
!    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
!

!  USE parameters_module, ONLY: pi

!  IMPLICIT NONE

  ! In/output variables:
  REAL(dp), INTENT(IN )   :: x
  REAL(dp), INTENT(OUT)   :: ber
  REAL(dp), INTENT(OUT)   :: bei
  REAL(dp), INTENT(OUT)   :: ger
  REAL(dp), INTENT(OUT)   :: gei
  REAL(dp), INTENT(OUT)   :: der
  REAL(dp), INTENT(OUT)   :: dei
  REAL(dp), INTENT(OUT)   :: her
  REAL(dp), INTENT(OUT)   :: hei

  ! Local variables:
  REAL(dp)   :: cn0
  REAL(dp)   :: cp0
  REAL(dp)   :: cs
  REAL(dp)   :: el
  REAL(dp)   :: eps
  REAL(dp)   :: fac
  REAL(dp)   :: gs

  INTEGER    :: k
  INTEGER    :: Km
  INTEGER    :: m

  REAL(dp)   :: pn0
  REAL(dp)   :: pn1
  REAL(dp)   :: pp0
  REAL(dp)   :: pp1
  REAL(dp)   :: qn0
  REAL(dp)   :: qn1
  REAL(dp)   :: qp0
  REAL(dp)   :: qp1
  REAL(dp)   :: r
  REAL(dp)   :: r0
  REAL(dp)   :: r1
  REAL(dp)   :: rc
  REAL(dp)   :: rs
  REAL(dp)   :: sn0
  REAL(dp)   :: sp0
  REAL(dp)   :: ss
  REAL(dp)   :: x2
  REAL(dp)   :: x4
  REAL(dp)   :: xc1
  REAL(dp)   :: xc2
  REAL(dp)   :: xd
  REAL(dp)   :: xe1
  REAL(dp)   :: xe2
  REAL(dp)   :: xt

  el = 0.5772156649015329_dp
  eps = 1.0D-15

  if ( x == 0.0_dp ) then
    ber = 1.0_dp
    bei = 0.0_dp
    ger = 1.0D+300
    gei = -0.25_dp * pi
    der = 0.0_dp
    dei = 0.0_dp
    her = -1.0D+300
    hei = 0.0_dp
    return
  end if

  x2 = 0.25_dp * x * x
  x4 = x2 * x2

  if ( abs ( x ) < 10.0_dp ) then

    ber = 1.0_dp
    r = 1.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) ** 2 * x4
      ber = ber + r
      if ( abs ( r ) < abs ( ber ) * eps ) then
        exit
      end if
    end do

    bei = x2
    r = x2
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      bei = bei + r
      if ( abs ( r ) < abs ( bei ) * eps ) then
        exit
      end if
    end do

    ger = - ( log ( x / 2.0_dp ) + el ) * ber + 0.25_dp * pi * bei
    r = 1.0_dp
    gs = 0.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m - 1.0_dp ) + 1.0_dp / ( 2.0_dp * m )
      ger = ger + r * gs
      if ( abs ( r * gs ) < abs ( ger ) * eps ) then
        exit
      end if
    end do

    gei = x2 - ( log ( x / 2.0_dp ) + el ) * bei - 0.25_dp * pi * ber
    r = x2
    gs = 1.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m ) + 1.0_dp / ( 2.0_dp * m + 1.0_dp )
      gei = gei + r * gs
      if ( abs ( r * gs ) < abs ( gei ) * eps ) then
        exit
      end if
    end do

    der = -0.25_dp * x * x2
    r = der
    do m = 1, 60
      r = -0.25_dp * r / m / ( m + 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      der = der + r
      if ( abs ( r ) < abs ( der ) * eps ) then
        exit
      end if
    end do

    dei = 0.5_dp * x
    r = dei
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) * x4
      dei = dei + r
      if ( abs ( r ) < abs ( dei ) * eps ) then
        exit
      end if
    end do

    r = -0.25_dp * x * x2
    gs = 1.5_dp
    her = 1.5_dp * r - ber / x &
      - ( log ( x / 2.0_dp ) + el ) * der + 0.25_dp * pi * dei
    do m = 1, 60
      r = -0.25_dp * r / m / ( m + 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2 * m + 1.0_dp ) + 1.0_dp &
        / ( 2 * m + 2.0_dp )
      her = her + r * gs
      if ( abs ( r * gs ) < abs ( her ) * eps ) then
        exit
      end if
    end do

    r = 0.5_dp * x
    gs = 1.0_dp
    hei = 0.5_dp * x - bei / x &
      - ( log ( x / 2.0_dp ) + el ) * dei - 0.25_dp * pi * der
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2 * m - 1.0_dp ) &
        / ( 2 * m + 1.0_dp ) * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m ) + 1.0_dp &
        / ( 2 * m + 1.0_dp )
      hei = hei + r * gs
      if ( abs ( r * gs ) < abs ( hei ) * eps ) then
        return
      end if
    end do

  else

    pp0 = 1.0_dp
    pn0 = 1.0_dp
    qp0 = 0.0_dp
    qn0 = 0.0_dp
    r0 = 1.0_dp

    if ( abs ( x ) < 40.0_dp ) then
      km = 18
    else
      km = 10
    end if

    fac = 1.0_dp
    do k = 1, km
      fac = -fac
      xt = 0.25_dp * k * pi - int ( 0.125_dp * k ) * 2.0_dp * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r0 = 0.125_dp * r0 * ( 2.0_dp * k - 1.0_dp ) ** 2 / k / x
      rc = r0 * cs
      rs = r0 * ss
      pp0 = pp0 + rc
      pn0 = pn0 + fac * rc
      qp0 = qp0 + rs
      qn0 = qn0 + fac * rs
    end do

    xd = x / sqrt (2.0_dp )
    xe1 = exp ( xd )
    xe2 = exp ( - xd )
    xc1 = 1.0_dp / sqrt ( 2.0_dp * pi * x )
    xc2 = sqrt ( 0.5_dp * pi / x )
    cp0 = cos ( xd + 0.125_dp * pi )
    cn0 = cos ( xd - 0.125_dp * pi )
    sp0 = sin ( xd + 0.125_dp * pi )
    sn0 = sin ( xd - 0.125_dp * pi )
    ger = xc2 * xe2 * (  pn0 * cp0 - qn0 * sp0 )
    gei = xc2 * xe2 * ( -pn0 * sp0 - qn0 * cp0 )
    ber = xc1 * xe1 * (  pp0 * cn0 + qp0 * sn0 ) - gei / pi
    bei = xc1 * xe1 * (  pp0 * sn0 - qp0 * cn0 ) + ger / pi
    pp1 = 1.0_dp
    pn1 = 1.0_dp
    qp1 = 0.0_dp
    qn1 = 0.0_dp
    r1 = 1.0_dp
    fac = 1.0_dp

    do k = 1, km
      fac = -fac
      xt = 0.25_dp * k * pi - int ( 0.125_dp * k ) * 2.0_dp * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r1 = 0.125_dp * r1 &
        * ( 4.0_dp - ( 2.0_dp * k - 1.0_dp ) ** 2 ) / k / x
      rc = r1 * cs
      rs = r1 * ss
      pp1 = pp1 + fac * rc
      pn1 = pn1 + rc
      qp1 = qp1 + fac * rs
      qn1 = qn1 + rs
    end do

    her = xc2 * xe2 * ( - pn1 * cn0 + qn1 * sn0 )
    hei = xc2 * xe2 * (   pn1 * sn0 + qn1 * cn0 )
    der = xc1 * xe1 * (   pp1 * cp0 + qp1 * sp0 ) - hei / pi
    dei = xc1 * xe1 * (   pp1 * sp0 - qp1 * cp0 ) + her / pi

  end if

  END SUBROUTINE klvna

END MODULE GIA_ELRA
