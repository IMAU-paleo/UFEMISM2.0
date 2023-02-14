MODULE math_utilities

  ! Some generally useful tools, and basic mathematical functions

! ===== Preamble =====
! ====================

  USE mpi
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, cerr, ierr, MPI_status, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, init_routine, finalise_routine
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate

  IMPLICIT NONE

  ! Interfaces to LAPACK, which are otherwise implicitly generated (taken from
  ! LAPACK source)
  !  *
  !  *  -- LAPACK routine (version 3.1) --
  !  *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  !  *     November 2006
  INTERFACE
    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
    END SUBROUTINE
    SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
    END SUBROUTINE
    SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
      INTEGER            INFO, LDB, N, NRHS
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
    END SUBROUTINE
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      INTEGER            INFO, LDA, LDB, N, NRHS
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
    END SUBROUTINE
  END INTERFACE

CONTAINS

! ===== Subroutinea ======
! ========================

! == Floatation criterion, surface elevation, and thickness above floatation

  PURE FUNCTION is_floating( Hi, Hb, SL) RESULT( isso)
    ! The floatation criterion

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: Hi     ! [m] Ice thickness
    REAL(dp)                                           , INTENT(IN)    :: Hb     ! [m] Bedrock elevation
    REAL(dp)                                           , INTENT(IN)    :: SL     ! [m] Water surface elevation

    ! Output variables:
    LOGICAL                                                            :: isso   ! Whether or not the ice will float

    isso = .FALSE.
    IF (Hi < (SL - Hb) * seawater_density/ice_density) isso = .TRUE.

  END FUNCTION is_floating

  PURE FUNCTION ice_surface_elevation( Hi, Hb, SL) RESULT( Hs)
    ! The ice surface elevation equation

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: Hi     ! [m] Ice thickness
    REAL(dp)                                           , INTENT(IN)    :: Hb     ! [m] Bedrock elevation
    REAL(dp)                                           , INTENT(IN)    :: SL     ! [m] Water surface elevation

    ! Output variables:
    REAL(dp)                                                           :: Hs     ! [m] Ice surface elevation

    Hs = Hi + MAX( SL - ice_density / seawater_density * Hi, Hb)

  END FUNCTION ice_surface_elevation

  PURE FUNCTION thickness_above_floatation( Hi, Hb, SL) RESULT( TAF)
    ! The thickness-above-floatation equation

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: Hi     ! [m] Ice thickness
    REAL(dp)                                           , INTENT(IN)    :: Hb     ! [m] Bedrock elevation
    REAL(dp)                                           , INTENT(IN)    :: SL     ! [m] Water surface elevation

    ! Output variables:
    REAL(dp)                                                           :: TAF    ! [m] Ice thickness above floatation

    TAF = Hi - MAX(0._dp, (SL - Hb) * (seawater_density / ice_density))

  END FUNCTION thickness_above_floatation

! == The error function

  PURE FUNCTION error_function( X) RESULT( ERR)
    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: X

    ! Output variables:
    REAL(dp)              :: ERR

    ! Local variables:
    REAL(dp)              :: EPS
    REAL(dp)              :: X2
    REAL(dp)              :: ER
    REAL(dp)              :: R
    REAL(dp)              :: C0
    INTEGER               :: k

    EPS = 1.0E-15_dp
    X2  = X * X
    IF(ABS(X) < 3.5_dp) THEN
     ER = 1.0_dp
     R  = 1.0_dp
     DO k = 1, 50
       R  = R * X2 / (REAL(k, dp) + 0.5_dp)
       ER = ER+R
       IF(ABS(R) < ABS(ER) * EPS) THEN
        C0  = 2.0_dp / SQRT(pi) * X * EXP(-X2)
        ERR = C0 * ER
        EXIT
       END IF
     END DO
    ELSE
     ER = 1.0_dp
     R  = 1.0_dp
     DO k = 1, 12
       R  = -R * (REAL(k, dp) - 0.5_dp) / X2
       ER = ER + R
       C0  = EXP(-X2) / (ABS(X) * SQRT(pi))
       ERR = 1.0_dp - C0 * ER
       IF(X < 0.0_dp) ERR = -ERR
     END DO
    ENDIF

    RETURN
  END FUNCTION error_function

! == The oblique stereographic projection

  PURE SUBROUTINE oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, beta_deg, x, y, k_P)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates to a rectangular coordinate system, with coordinates (x,y).
    !
    ! For more information about M, beta_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    ! For North and South Pole: lambda_M_deg = 0._dp, to generate the correct coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: lambda         ! [degrees east ] Longitude
    REAL(dp)                                           , INTENT(IN)    :: phi            ! [degrees north] Latitude

    ! Polar stereographic projection parameters
    REAL(dp)                                           , INTENT(IN)    :: lambda_M_deg   ! [degrees] Central parallel
    REAL(dp)                                           , INTENT(IN)    :: phi_M_deg      ! [degrees] Central meridian
    REAL(dp)                                           , INTENT(IN)    :: beta_deg       ! [degrees] Stereographic projection angle

    ! Output variables:
    REAL(dp)                                           , INTENT(OUT)   :: x              ! [m] x-coordinate
    REAL(dp)                                           , INTENT(OUT)   :: y              ! [m] y-coordinate
    REAL(dp)                                 , OPTIONAL, INTENT(OUT)   :: k_P            ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                                                           :: alpha_deg      ! [degrees]
    REAL(dp)                                                           :: phi_P          ! [radians]
    REAL(dp)                                                           :: lambda_P       ! [radians]
    REAL(dp)                                                           :: t_P_prime
    REAL(dp)                                                           :: lambda_M, phi_M, alpha

    ! Convert beta to alpha
    alpha_deg = 90._dp - beta_deg

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180._dp) * phi
    lambda_P = (pi / 180._dp) * lambda

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + COS(alpha)) / (1._dp + COS(phi_P) * COS(phi_M) * COS(lambda_P - lambda_M) + SIN(phi_P) * SIN(phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x =  earth_radius * (COS(phi_P) * SIN(lambda_P - lambda_M)) * t_P_prime
    y =  earth_radius * (SIN(phi_P) * COS(phi_M) - (COS(phi_P) * SIN(phi_M)) * COS(lambda_P - lambda_M)) * t_P_prime

!    ! See equation (21-4) on page 157 in Snyder (1987):
!    IF(PRESENT(k_P)) k_P = (1._dp + COS(alpha)) / (1._dp + SIN(phi_M) * SIN(phi_P) + COS(phi_M) * COS(phi_P) * COS(lambda_P - lambda_M))

  END SUBROUTINE oblique_sg_projection

  PURE SUBROUTINE inverse_oblique_sg_projection( x, y, lambda_M_deg, phi_M_deg, beta_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates to a longitude-latitude coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: x              ! [m] x-coordinate
    REAL(dp)                                           , INTENT(IN)    :: y              ! [m] y-coordinate

    ! Polar stereographic projection parameters
    REAL(dp)                                           , INTENT(IN)    :: lambda_M_deg   ! [degrees] Central parallel
    REAL(dp)                                           , INTENT(IN)    :: phi_M_deg      ! [degrees] Central meridian
    REAL(dp)                                           , INTENT(IN)    :: beta_deg       ! [degrees] Stereographic projection angle

    ! Output variables:
    REAL(dp)                                           , INTENT(OUT)   :: lambda_P       ! [degrees east ] Longitude
    REAL(dp)                                           , INTENT(OUT)   :: phi_P          ! [degrees north] Latitude

    ! Local variables:
    REAL(dp)                                                           :: alpha_deg      ! [degrees]
    REAL(dp)                                                           :: x_3D_P_prime   ! [m]
    REAL(dp)                                                           :: y_3D_P_prime   ! [m]
    REAL(dp)                                                           :: z_3D_P_prime   ! [m]
    REAL(dp)                                                           :: a
    REAL(dp)                                                           :: t_P
    REAL(dp)                                                           :: x_3D_P         ! [m]
    REAL(dp)                                                           :: y_3D_P         ! [m]
    REAL(dp)                                                           :: z_3D_P         ! [m]
    REAL(dp)                                                           :: lambda_M, phi_M, alpha

    ! Convert beta to alpha
    alpha_deg = 90._dp - beta_deg

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = earth_radius * COS(alpha) * COS(lambda_M) * COS(phi_M) - SIN(lambda_M) * x - COS(lambda_M) * SIN(phi_M) * y
    y_3D_P_prime = earth_radius * COS(alpha) * SIN(lambda_M) * COS(phi_M) + COS(lambda_M) * y - SIN(lambda_M) * SIN(phi_M) * y
    z_3D_P_prime = earth_radius * COS(alpha) *                 SIN(phi_M)                     +                 COS(phi_M) * y

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(lambda_M) * COS(phi_M) * x_3D_P_prime  +  SIN(lambda_M) * COS(phi_M) * y_3D_P_prime  +  SIN(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * earth_radius**2 + 2._dp * earth_radius * a) / (earth_radius**2 + 2._dp * earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  earth_radius * COS(lambda_M) * COS(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  earth_radius * SIN(lambda_M) * COS(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  earth_radius *                 SIN(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
      lambda_P = 180._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
      lambda_P =           (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 360._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
      lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
      lambda_P =   0._dp
    END IF

    ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
    IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
      phi_P = (180._dp / pi) * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
    ELSE IF(z_3D_P >  0._dp) THEN
      phi_P =   90._dp
    ELSE IF(z_3D_P <  0._dp) THEN
      phi_P =  -90._dp
    END IF

  END SUBROUTINE inverse_oblique_sg_projection

! == Line integrals used in conservative remapping

  PURE FUNCTION line_integral_xdy(   p, q, tol_dist) RESULT( I_pq)
    ! Calculate the line integral x dy from p to q

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(2)                             , INTENT(IN)    :: p, q
    REAL(dp)                                           , INTENT(IN)    :: tol_dist

    ! Output variables:
    REAL(dp)                                                           :: I_pq

    ! Local variables:
    REAL(dp)                                                           :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    IF (ABS( yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF

    dx = q( 1) - p( 1)
    dy = q( 2) - p( 2)

    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)

  END FUNCTION line_integral_xdy

  PURE FUNCTION line_integral_mxydx( p, q, tol_dist) RESULT( I_pq)
    ! Calculate the line integral -xy dx from p to q

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(2)                             , INTENT(IN)    :: p, q
    REAL(dp)                                           , INTENT(IN)    :: tol_dist

    ! Output variables:
    REAL(dp)                                                           :: I_pq

    ! Local variables:
    REAL(dp)                                                           :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    IF (ABS( xp-xq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF

    dx = q( 1) - p( 1)
    dy = q( 2) - p( 2)

    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))

  END FUNCTION line_integral_mxydx

  PURE FUNCTION line_integral_xydy(  p, q, tol_dist) RESULT( I_pq)
    ! Calculate the line integral xy dy from p to q

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(2)                             , INTENT(IN)    :: p, q
    REAL(dp)                                           , INTENT(IN)    :: tol_dist

    ! Output variables:
    REAL(dp)                                                           :: I_pq

    ! Local variables:
    REAL(dp)                                                           :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    IF (ABS( yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF

    dx = q( 1 ) - p( 1)
    dy = q( 2 ) - p( 2)

    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))

  END FUNCTION line_integral_xydy

! == Some wrappers for LAPACK matrix functionality

  FUNCTION tridiagonal_solve( ldiag, diag, udiag, b) RESULT(x)
    ! Solve the matrix equation Ax = b, where A is a tridiagonal matrix.
    ! Provide only the three diagonals of A as vectors.
    ! Uses the LAPACK function DGTSV to solve the equation.

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      )                       , INTENT(IN)    ::  diag    !       Diagonal of matrix A
    REAL(dp), DIMENSION(SIZE(diag)-1)                  , INTENT(IN)    :: ldiag    ! Lower diagonal of matrix A
    REAL(dp), DIMENSION(SIZE(diag)-1)                  , INTENT(IN)    :: udiag    ! Upper diagonal of matrix A
    REAL(dp), DIMENSION(SIZE(diag)  )                  , INTENT(IN)    :: b        ! Right-hand side of the equation

    ! Output variables:
    REAL(dp), DIMENSION(SIZE(diag)  )                                  :: x        ! Solution to the matrix equation

    ! Local variables:
    INTEGER                                                            :: info
    REAL(dp), DIMENSION(SIZE(diag)  )                                  :: diag_copy
    REAL(dp), DIMENSION(SIZE(diag)-1)                                  :: udiag_copy, ldiag_copy

    ! The LAPACK solver will overwrite the right-hand side b with the solution x. Therefore we
    ! first copy b to the solution vector x:
    x = b

    ! The LAPACK solver will change the elements in the matrix, therefore we copy them:
    diag_copy  =  diag
    udiag_copy = udiag
    ldiag_copy = ldiag

    CALL DGTSV( SIZE( diag), 1, ldiag_copy, diag_copy, udiag_copy, x, SIZE( diag), info)

    IF (info /= 0) THEN
      CALL crash('LAPACK solver DGTSV returned error message info = {int_01}', int_01 = info)
    END IF

  END FUNCTION tridiagonal_solve

! == Finding inverses of some small matrices

  PURE FUNCTION calc_matrix_inverse_2_by_2( A) RESULT( Ainv)
    ! Direct inversion of a 2-by-2 matrix

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(2,2    )                       , INTENT(IN)    :: A       ! The matrix A to be inverted

    ! Output variables:
    REAL(dp), DIMENSION(2,2    )                                       :: Ainv    ! Inverse of A

    ! Local variables:
    REAL(dp)                                                           :: detA

    ! Calculate the determinant of A
    detA = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Safety
    IF (ABS( detA) < TINY( detA)) THEN
      error stop 'calc_matrix_inverse_2_by_2 - matrix is singular to working precision!'
    END IF

    ! Calculate the inverse of A
    Ainv( 1,1) =  A( 2,2) / detA
    Ainv( 1,2) = -A( 1,2) / detA
    Ainv( 2,1) = -A( 2,1) / detA
    Ainv( 2,2) =  A( 1,1) / detA

  END FUNCTION calc_matrix_inverse_2_by_2

  PURE FUNCTION calc_matrix_inverse_3_by_3( A) RESULT( Ainv)
    ! Direct inversion of a 3-by-3 matrix
    !
    ! See: https://metric.ma.ic.ac.uk/metric_public/matrices/inverses/inverses2.html

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(3,3    )                       , INTENT(IN)    :: A       ! The matrix A to be inverted

    ! Output variables:
    REAL(dp), DIMENSION(3,3    )                                       :: Ainv    ! Inverse of A

    ! Local variables:
    REAL(dp)                                                           :: detA

    ! Calculate the minors of A
    Ainv( 1,1) = A( 2,2) * A( 3,3) - A( 2,3) * A( 3,2)
    Ainv( 1,2) = A( 2,1) * A( 3,3) - A( 2,3) * A( 3,1)
    Ainv( 1,3) = A( 2,1) * A( 3,2) - A( 2,2) * A( 3,1)
    Ainv( 2,1) = A( 1,2) * A( 3,3) - A( 1,3) * A( 3,2)
    Ainv( 2,2) = A( 1,1) * A( 3,3) - A( 1,3) * A( 3,1)
    Ainv( 2,3) = A( 1,1) * A( 3,2) - A( 1,2) * A( 3,1)
    Ainv( 3,1) = A( 1,2) * A( 2,3) - A( 1,3) * A( 2,2)
    Ainv( 3,2) = A( 1,1) * A( 2,3) - A( 1,3) * A( 2,1)
    Ainv( 3,3) = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Calculate the determinant of A
    detA = A( 1,1) * Ainv( 1,1) - A( 1,2) * Ainv( 1,2) + A( 1,3) * Ainv( 1,3)

    ! Safety
    IF (ABS( detA) < TINY( detA)) THEN
      error stop 'calc_matrix_inverse_3_by_3 - matrix is singular to working precision!'
    END IF

    ! Change matrix of minors to get the matrix of cofactors
    Ainv( 1,2) = -Ainv( 1,2)
    Ainv( 2,1) = -Ainv( 2,1)
    Ainv( 2,3) = -Ainv( 2,3)
    Ainv( 3,2) = -Ainv( 3,2)

    ! Transpose matrix of cofactors
    Ainv( 1,2) = Ainv( 1,2) + Ainv( 2,1)
    Ainv( 2,1) = Ainv( 1,2) - Ainv( 2,1)
    Ainv( 1,2) = Ainv( 1,2) - Ainv( 2,1)

    Ainv( 1,3) = Ainv( 1,3) + Ainv( 3,1)
    Ainv( 3,1) = Ainv( 1,3) - Ainv( 3,1)
    Ainv( 1,3) = Ainv( 1,3) - Ainv( 3,1)

    Ainv( 2,3) = Ainv( 2,3) + Ainv( 3,2)
    Ainv( 3,2) = Ainv( 2,3) - Ainv( 3,2)
    Ainv( 2,3) = Ainv( 2,3) - Ainv( 3,2)

    ! Divide by det(A)
    Ainv = Ainv / detA

  END FUNCTION calc_matrix_inverse_3_by_3

  PURE FUNCTION calc_matrix_inverse_5_by_5( A) RESULT( Ainv)
    ! Direct inversion of a 5-by-5 matrix
    !
    ! Source: https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt, accessed 2023-02-07

    USE iso_fortran_env, ONLY: real128

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(5,5    )                       , INTENT(IN)    :: A       ! The matrix A to be inverted

    ! Output variables:
    REAL(dp), DIMENSION(5,5    )                                       :: Ainv    ! Inverse of A

    ! Local variables:
    REAL(dp)                                                           :: detA    ! Determinant of A

    ! Local variables:
    REAL(real128) :: A11, A12, A13, A14, A15, A21, A22, A23, A24, &
         A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
         A51, A52, A53, A54, A55
    REAL(real128), DIMENSION(5,5) :: COFACTOR

    A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
    A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
    A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
    A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
    A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

    detA = REAL(                                                          &
       A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+       &
       A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-       &
       A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-       &
       A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+       &
       A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+       &
       A12*A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-       &
       A14*A23*A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-       &
       A12*A24*A33*A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-       &
       A15*A24*A33*A41*A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-       &
       A13*A25*A34*A41*A52-A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+       &
       A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-A15*A21*A34*A43*A52+       &
       A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*A24*A35*A43*A52-       &
       A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*A33*A44*A52-       &
       A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*A44*A52+       &
       A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*A52+       &
       A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+       &
       A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+       &
       A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-       &
       A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-       &
       A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+       &
       A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+       &
       A11*A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-       &
       A14*A22*A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-       &
       A11*A24*A32*A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-       &
       A15*A23*A32*A41*A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-       &
       A12*A25*A33*A41*A54-A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+       &
       A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-A15*A21*A33*A42*A54+       &
       A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*A23*A35*A42*A54-       &
       A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*A32*A43*A54-       &
       A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*A43*A54+       &
       A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*A54+       &
       A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+       &
       A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+       &
       A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-       &
       A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-       &
       A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+       &
       A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+       &
       A11*A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-       &
       A13*A22*A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-       &
       A11*A23*A32*A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55, dp)

    ! Safety
    IF (ABS( detA) < TINY( detA)) THEN
      error stop 'calc_matrix_inverse_3_by_3 - matrix is singular to working precision!'
    END IF

    COFACTOR(1,1) = A25*A34*A43*A52-A24*A35*A43*A52-A25*A33*A44*A52+      &
       A23*A35*A44*A52+A24*A33*A45*A52-A23*A34*A45*A52-A25*A34*A42*A53+   &
       A24*A35*A42*A53+A25*A32*A44*A53-A22*A35*A44*A53-A24*A32*A45*A53+   &
       A22*A34*A45*A53+A25*A33*A42*A54-A23*A35*A42*A54-A25*A32*A43*A54+   &
       A22*A35*A43*A54+A23*A32*A45*A54-A22*A33*A45*A54-A24*A33*A42*A55+   &
       A23*A34*A42*A55+A24*A32*A43*A55-A22*A34*A43*A55-A23*A32*A44*A55+   &
       A22*A33*A44*A55

    COFACTOR(2,1) = -A15*A34*A43*A52+A14*A35*A43*A52+A15*A33*A44*A52-     &
       A13*A35*A44*A52-A14*A33*A45*A52+A13*A34*A45*A52+A15*A34*A42*A53-   &
       A14*A35*A42*A53-A15*A32*A44*A53+A12*A35*A44*A53+A14*A32*A45*A53-   &
       A12*A34*A45*A53-A15*A33*A42*A54+A13*A35*A42*A54+A15*A32*A43*A54-   &
       A12*A35*A43*A54-A13*A32*A45*A54+A12*A33*A45*A54+A14*A33*A42*A55-   &
       A13*A34*A42*A55-A14*A32*A43*A55+A12*A34*A43*A55+A13*A32*A44*A55-   &
       A12*A33*A44*A55

    COFACTOR(3,1) = A15*A24*A43*A52-A14*A25*A43*A52-A15*A23*A44*A52+      &
       A13*A25*A44*A52+A14*A23*A45*A52-A13*A24*A45*A52-A15*A24*A42*A53+   &
       A14*A25*A42*A53+A15*A22*A44*A53-A12*A25*A44*A53-A14*A22*A45*A53+   &
       A12*A24*A45*A53+A15*A23*A42*A54-A13*A25*A42*A54-A15*A22*A43*A54+   &
       A12*A25*A43*A54+A13*A22*A45*A54-A12*A23*A45*A54-A14*A23*A42*A55+   &
       A13*A24*A42*A55+A14*A22*A43*A55-A12*A24*A43*A55-A13*A22*A44*A55+   &
       A12*A23*A44*A55

    COFACTOR(4,1) = -A15*A24*A33*A52+A14*A25*A33*A52+A15*A23*A34*A52-     &
       A13*A25*A34*A52-A14*A23*A35*A52+A13*A24*A35*A52+A15*A24*A32*A53-   &
       A14*A25*A32*A53-A15*A22*A34*A53+A12*A25*A34*A53+A14*A22*A35*A53-   &
       A12*A24*A35*A53-A15*A23*A32*A54+A13*A25*A32*A54+A15*A22*A33*A54-   &
       A12*A25*A33*A54-A13*A22*A35*A54+A12*A23*A35*A54+A14*A23*A32*A55-   &
       A13*A24*A32*A55-A14*A22*A33*A55+A12*A24*A33*A55+A13*A22*A34*A55-   &
       A12*A23*A34*A55

    COFACTOR(5,1) = A15*A24*A33*A42-A14*A25*A33*A42-A15*A23*A34*A42+      &
       A13*A25*A34*A42+A14*A23*A35*A42-A13*A24*A35*A42-A15*A24*A32*A43+   &
       A14*A25*A32*A43+A15*A22*A34*A43-A12*A25*A34*A43-A14*A22*A35*A43+   &
       A12*A24*A35*A43+A15*A23*A32*A44-A13*A25*A32*A44-A15*A22*A33*A44+   &
       A12*A25*A33*A44+A13*A22*A35*A44-A12*A23*A35*A44-A14*A23*A32*A45+   &
       A13*A24*A32*A45+A14*A22*A33*A45-A12*A24*A33*A45-A13*A22*A34*A45+   &
       A12*A23*A34*A45

    COFACTOR(1,2) = -A25*A34*A43*A51+A24*A35*A43*A51+A25*A33*A44*A51-     &
       A23*A35*A44*A51-A24*A33*A45*A51+A23*A34*A45*A51+A25*A34*A41*A53-   &
       A24*A35*A41*A53-A25*A31*A44*A53+A21*A35*A44*A53+A24*A31*A45*A53-   &
       A21*A34*A45*A53-A25*A33*A41*A54+A23*A35*A41*A54+A25*A31*A43*A54-   &
       A21*A35*A43*A54-A23*A31*A45*A54+A21*A33*A45*A54+A24*A33*A41*A55-   &
       A23*A34*A41*A55-A24*A31*A43*A55+A21*A34*A43*A55+A23*A31*A44*A55-   &
       A21*A33*A44*A55

    COFACTOR(2,2) = A15*A34*A43*A51-A14*A35*A43*A51-A15*A33*A44*A51+      &
       A13*A35*A44*A51+A14*A33*A45*A51-A13*A34*A45*A51-A15*A34*A41*A53+   &
       A14*A35*A41*A53+A15*A31*A44*A53-A11*A35*A44*A53-A14*A31*A45*A53+   &
       A11*A34*A45*A53+A15*A33*A41*A54-A13*A35*A41*A54-A15*A31*A43*A54+   &
       A11*A35*A43*A54+A13*A31*A45*A54-A11*A33*A45*A54-A14*A33*A41*A55+   &
       A13*A34*A41*A55+A14*A31*A43*A55-A11*A34*A43*A55-A13*A31*A44*A55+   &
       A11*A33*A44*A55

    COFACTOR(3,2) = -A15*A24*A43*A51+A14*A25*A43*A51+A15*A23*A44*A51-     &
       A13*A25*A44*A51-A14*A23*A45*A51+A13*A24*A45*A51+A15*A24*A41*A53-   &
       A14*A25*A41*A53-A15*A21*A44*A53+A11*A25*A44*A53+A14*A21*A45*A53-   &
       A11*A24*A45*A53-A15*A23*A41*A54+A13*A25*A41*A54+A15*A21*A43*A54-   &
       A11*A25*A43*A54-A13*A21*A45*A54+A11*A23*A45*A54+A14*A23*A41*A55-   &
       A13*A24*A41*A55-A14*A21*A43*A55+A11*A24*A43*A55+A13*A21*A44*A55-   &
       A11*A23*A44*A55

    COFACTOR(4,2) = A15*A24*A33*A51-A14*A25*A33*A51-A15*A23*A34*A51+      &
       A13*A25*A34*A51+A14*A23*A35*A51-A13*A24*A35*A51-A15*A24*A31*A53+   &
       A14*A25*A31*A53+A15*A21*A34*A53-A11*A25*A34*A53-A14*A21*A35*A53+   &
       A11*A24*A35*A53+A15*A23*A31*A54-A13*A25*A31*A54-A15*A21*A33*A54+   &
       A11*A25*A33*A54+A13*A21*A35*A54-A11*A23*A35*A54-A14*A23*A31*A55+   &
       A13*A24*A31*A55+A14*A21*A33*A55-A11*A24*A33*A55-A13*A21*A34*A55+   &
       A11*A23*A34*A55

    COFACTOR(5,2) = -A15*A24*A33*A41+A14*A25*A33*A41+A15*A23*A34*A41-     &
       A13*A25*A34*A41-A14*A23*A35*A41+A13*A24*A35*A41+A15*A24*A31*A43-   &
       A14*A25*A31*A43-A15*A21*A34*A43+A11*A25*A34*A43+A14*A21*A35*A43-   &
       A11*A24*A35*A43-A15*A23*A31*A44+A13*A25*A31*A44+A15*A21*A33*A44-   &
       A11*A25*A33*A44-A13*A21*A35*A44+A11*A23*A35*A44+A14*A23*A31*A45-   &
       A13*A24*A31*A45-A14*A21*A33*A45+A11*A24*A33*A45+A13*A21*A34*A45-   &
       A11*A23*A34*A45

    COFACTOR(1,3) = A25*A34*A42*A51-A24*A35*A42*A51-A25*A32*A44*A51+      &
       A22*A35*A44*A51+A24*A32*A45*A51-A22*A34*A45*A51-A25*A34*A41*A52+   &
       A24*A35*A41*A52+A25*A31*A44*A52-A21*A35*A44*A52-A24*A31*A45*A52+   &
       A21*A34*A45*A52+A25*A32*A41*A54-A22*A35*A41*A54-A25*A31*A42*A54+   &
       A21*A35*A42*A54+A22*A31*A45*A54-A21*A32*A45*A54-A24*A32*A41*A55+   &
       A22*A34*A41*A55+A24*A31*A42*A55-A21*A34*A42*A55-A22*A31*A44*A55+   &
       A21*A32*A44*A55

    COFACTOR(2,3) = -A15*A34*A42*A51+A14*A35*A42*A51+A15*A32*A44*A51-     &
       A12*A35*A44*A51-A14*A32*A45*A51+A12*A34*A45*A51+A15*A34*A41*A52-   &
       A14*A35*A41*A52-A15*A31*A44*A52+A11*A35*A44*A52+A14*A31*A45*A52-   &
       A11*A34*A45*A52-A15*A32*A41*A54+A12*A35*A41*A54+A15*A31*A42*A54-   &
       A11*A35*A42*A54-A12*A31*A45*A54+A11*A32*A45*A54+A14*A32*A41*A55-   &
       A12*A34*A41*A55-A14*A31*A42*A55+A11*A34*A42*A55+A12*A31*A44*A55-   &
       A11*A32*A44*A55

    COFACTOR(3,3) = A15*A24*A42*A51-A14*A25*A42*A51-A15*A22*A44*A51+      &
       A12*A25*A44*A51+A14*A22*A45*A51-A12*A24*A45*A51-A15*A24*A41*A52+   &
       A14*A25*A41*A52+A15*A21*A44*A52-A11*A25*A44*A52-A14*A21*A45*A52+   &
       A11*A24*A45*A52+A15*A22*A41*A54-A12*A25*A41*A54-A15*A21*A42*A54+   &
       A11*A25*A42*A54+A12*A21*A45*A54-A11*A22*A45*A54-A14*A22*A41*A55+   &
       A12*A24*A41*A55+A14*A21*A42*A55-A11*A24*A42*A55-A12*A21*A44*A55+   &
       A11*A22*A44*A55

    COFACTOR(4,3) = -A15*A24*A32*A51+A14*A25*A32*A51+A15*A22*A34*A51-     &
       A12*A25*A34*A51-A14*A22*A35*A51+A12*A24*A35*A51+A15*A24*A31*A52-   &
       A14*A25*A31*A52-A15*A21*A34*A52+A11*A25*A34*A52+A14*A21*A35*A52-   &
       A11*A24*A35*A52-A15*A22*A31*A54+A12*A25*A31*A54+A15*A21*A32*A54-   &
       A11*A25*A32*A54-A12*A21*A35*A54+A11*A22*A35*A54+A14*A22*A31*A55-   &
       A12*A24*A31*A55-A14*A21*A32*A55+A11*A24*A32*A55+A12*A21*A34*A55-   &
       A11*A22*A34*A55

    COFACTOR(5,3) = A15*A24*A32*A41-A14*A25*A32*A41-A15*A22*A34*A41+      &
       A12*A25*A34*A41+A14*A22*A35*A41-A12*A24*A35*A41-A15*A24*A31*A42+   &
       A14*A25*A31*A42+A15*A21*A34*A42-A11*A25*A34*A42-A14*A21*A35*A42+   &
       A11*A24*A35*A42+A15*A22*A31*A44-A12*A25*A31*A44-A15*A21*A32*A44+   &
       A11*A25*A32*A44+A12*A21*A35*A44-A11*A22*A35*A44-A14*A22*A31*A45+   &
       A12*A24*A31*A45+A14*A21*A32*A45-A11*A24*A32*A45-A12*A21*A34*A45+   &
       A11*A22*A34*A45

    COFACTOR(1,4) = -A25*A33*A42*A51+A23*A35*A42*A51+A25*A32*A43*A51-     &
       A22*A35*A43*A51-A23*A32*A45*A51+A22*A33*A45*A51+A25*A33*A41*A52-   &
       A23*A35*A41*A52-A25*A31*A43*A52+A21*A35*A43*A52+A23*A31*A45*A52-   &
       A21*A33*A45*A52-A25*A32*A41*A53+A22*A35*A41*A53+A25*A31*A42*A53-   &
       A21*A35*A42*A53-A22*A31*A45*A53+A21*A32*A45*A53+A23*A32*A41*A55-   &
       A22*A33*A41*A55-A23*A31*A42*A55+A21*A33*A42*A55+A22*A31*A43*A55-   &
       A21*A32*A43*A55

    COFACTOR(2,4) = A15*A33*A42*A51-A13*A35*A42*A51-A15*A32*A43*A51+      &
       A12*A35*A43*A51+A13*A32*A45*A51-A12*A33*A45*A51-A15*A33*A41*A52+   &
       A13*A35*A41*A52+A15*A31*A43*A52-A11*A35*A43*A52-A13*A31*A45*A52+   &
       A11*A33*A45*A52+A15*A32*A41*A53-A12*A35*A41*A53-A15*A31*A42*A53+   &
       A11*A35*A42*A53+A12*A31*A45*A53-A11*A32*A45*A53-A13*A32*A41*A55+   &
       A12*A33*A41*A55+A13*A31*A42*A55-A11*A33*A42*A55-A12*A31*A43*A55+   &
       A11*A32*A43*A55

    COFACTOR(3,4) = -A15*A23*A42*A51+A13*A25*A42*A51+A15*A22*A43*A51-     &
       A12*A25*A43*A51-A13*A22*A45*A51+A12*A23*A45*A51+A15*A23*A41*A52-   &
       A13*A25*A41*A52-A15*A21*A43*A52+A11*A25*A43*A52+A13*A21*A45*A52-   &
       A11*A23*A45*A52-A15*A22*A41*A53+A12*A25*A41*A53+A15*A21*A42*A53-   &
       A11*A25*A42*A53-A12*A21*A45*A53+A11*A22*A45*A53+A13*A22*A41*A55-   &
       A12*A23*A41*A55-A13*A21*A42*A55+A11*A23*A42*A55+A12*A21*A43*A55-   &
       A11*A22*A43*A55

    COFACTOR(4,4) = A15*A23*A32*A51-A13*A25*A32*A51-A15*A22*A33*A51+      &
       A12*A25*A33*A51+A13*A22*A35*A51-A12*A23*A35*A51-A15*A23*A31*A52+   &
       A13*A25*A31*A52+A15*A21*A33*A52-A11*A25*A33*A52-A13*A21*A35*A52+   &
       A11*A23*A35*A52+A15*A22*A31*A53-A12*A25*A31*A53-A15*A21*A32*A53+   &
       A11*A25*A32*A53+A12*A21*A35*A53-A11*A22*A35*A53-A13*A22*A31*A55+   &
       A12*A23*A31*A55+A13*A21*A32*A55-A11*A23*A32*A55-A12*A21*A33*A55+   &
       A11*A22*A33*A55

    COFACTOR(5,4) = -A15*A23*A32*A41+A13*A25*A32*A41+A15*A22*A33*A41-     &
       A12*A25*A33*A41-A13*A22*A35*A41+A12*A23*A35*A41+A15*A23*A31*A42-   &
       A13*A25*A31*A42-A15*A21*A33*A42+A11*A25*A33*A42+A13*A21*A35*A42-   &
       A11*A23*A35*A42-A15*A22*A31*A43+A12*A25*A31*A43+A15*A21*A32*A43-   &
       A11*A25*A32*A43-A12*A21*A35*A43+A11*A22*A35*A43+A13*A22*A31*A45-   &
       A12*A23*A31*A45-A13*A21*A32*A45+A11*A23*A32*A45+A12*A21*A33*A45-   &
       A11*A22*A33*A45

    COFACTOR(1,5) = A24*A33*A42*A51-A23*A34*A42*A51-A24*A32*A43*A51+      &
       A22*A34*A43*A51+A23*A32*A44*A51-A22*A33*A44*A51-A24*A33*A41*A52+   &
       A23*A34*A41*A52+A24*A31*A43*A52-A21*A34*A43*A52-A23*A31*A44*A52+   &
       A21*A33*A44*A52+A24*A32*A41*A53-A22*A34*A41*A53-A24*A31*A42*A53+   &
       A21*A34*A42*A53+A22*A31*A44*A53-A21*A32*A44*A53-A23*A32*A41*A54+   &
       A22*A33*A41*A54+A23*A31*A42*A54-A21*A33*A42*A54-A22*A31*A43*A54+   &
       A21*A32*A43*A54

    COFACTOR(2,5) = -A14*A33*A42*A51+A13*A34*A42*A51+A14*A32*A43*A51-     &
       A12*A34*A43*A51-A13*A32*A44*A51+A12*A33*A44*A51+A14*A33*A41*A52-   &
       A13*A34*A41*A52-A14*A31*A43*A52+A11*A34*A43*A52+A13*A31*A44*A52-   &
       A11*A33*A44*A52-A14*A32*A41*A53+A12*A34*A41*A53+A14*A31*A42*A53-   &
       A11*A34*A42*A53-A12*A31*A44*A53+A11*A32*A44*A53+A13*A32*A41*A54-   &
       A12*A33*A41*A54-A13*A31*A42*A54+A11*A33*A42*A54+A12*A31*A43*A54-   &
       A11*A32*A43*A54

    COFACTOR(3,5) = A14*A23*A42*A51-A13*A24*A42*A51-A14*A22*A43*A51+      &
       A12*A24*A43*A51+A13*A22*A44*A51-A12*A23*A44*A51-A14*A23*A41*A52+   &
       A13*A24*A41*A52+A14*A21*A43*A52-A11*A24*A43*A52-A13*A21*A44*A52+   &
       A11*A23*A44*A52+A14*A22*A41*A53-A12*A24*A41*A53-A14*A21*A42*A53+   &
       A11*A24*A42*A53+A12*A21*A44*A53-A11*A22*A44*A53-A13*A22*A41*A54+   &
       A12*A23*A41*A54+A13*A21*A42*A54-A11*A23*A42*A54-A12*A21*A43*A54+   &
       A11*A22*A43*A54

    COFACTOR(4,5) = -A14*A23*A32*A51+A13*A24*A32*A51+A14*A22*A33*A51-     &
       A12*A24*A33*A51-A13*A22*A34*A51+A12*A23*A34*A51+A14*A23*A31*A52-   &
       A13*A24*A31*A52-A14*A21*A33*A52+A11*A24*A33*A52+A13*A21*A34*A52-   &
       A11*A23*A34*A52-A14*A22*A31*A53+A12*A24*A31*A53+A14*A21*A32*A53-   &
       A11*A24*A32*A53-A12*A21*A34*A53+A11*A22*A34*A53+A13*A22*A31*A54-   &
       A12*A23*A31*A54-A13*A21*A32*A54+A11*A23*A32*A54+A12*A21*A33*A54-   &
       A11*A22*A33*A54

    COFACTOR(5,5) = A14*A23*A32*A41-A13*A24*A32*A41-A14*A22*A33*A41+      &
       A12*A24*A33*A41+A13*A22*A34*A41-A12*A23*A34*A41-A14*A23*A31*A42+   &
       A13*A24*A31*A42+A14*A21*A33*A42-A11*A24*A33*A42-A13*A21*A34*A42+   &
       A11*A23*A34*A42+A14*A22*A31*A43-A12*A24*A31*A43-A14*A21*A32*A43+   &
       A11*A24*A32*A43+A12*A21*A34*A43-A11*A22*A34*A43-A13*A22*A31*A44+   &
       A12*A23*A31*A44+A13*A21*A32*A44-A11*A23*A32*A44-A12*A21*A33*A44+   &
       A11*A22*A33*A44

    AINV = REAL( TRANSPOSE( COFACTOR),dp) / detA

  END FUNCTIOn calc_matrix_inverse_5_by_5

  SUBROUTINE calc_matrix_inverse_general( A, Ainv)
    ! Calculate the inverse Ainv of an n-by-n matrix A using LAPACK

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    )                       , INTENT(IN)    :: A       ! The matrix A to be inverted

    ! Output variables:
    REAL(dp), DIMENSION( SIZE( A,1), SIZE( A,1))       , INTENT(OUT)   :: Ainv    ! Inverse of A

    ! Local variables:
    REAL(dp), DIMENSION( SIZE( A,1))                                   :: work     ! work array for LAPACK
    INTEGER,  DIMENSION( SIZE( A,1))                                   :: ipiv     ! pivot indices
    INTEGER                                                            :: n,info

    n = size( A,1)

    ! Safety
    IF (SIZE (A,2) /= n) CALL crash('calc_matrix_inverse_general: matrix must be square!')

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
    CALL DGETRF( n, n, Ainv, n, ipiv, info)

    ! Safety
    IF (info /= 0) THEN
      CALL crash('calc_matrix_inverse_general: DGETRF error: matrix is singular to working precision!')
    END IF

    ! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
    CALL DGETRI( n, Ainv, n, ipiv, work, n, info)

    ! Safety
    IF (info /= 0) THEN
      CALL crash('calc_matrix_inverse_general: DGETRI error: matrix inversion failed!')
    END IF

  END SUBROUTINE calc_matrix_inverse_general

  PURE FUNCTION solve_Axb_2_by_2( A,b) RESULT( x)
    ! Direct solution of the 2-by-2 matrix equation Ax=b

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(2,2    )                       , INTENT(IN)    :: A       ! The matrix A
    REAL(dp), DIMENSION(2      )                       , INTENT(IN)    :: b       ! The right-hand side b

    ! Output variables:
    REAL(dp), DIMENSION(2      )                                       :: x       ! THe solution x

    ! Local variables:
    REAL(dp)                                                           :: detA
    REAL(dp), DIMENSION(2,2    )                                       :: Ainv

    ! Calculate the inverse of A
    Ainv = calc_matrix_inverse_2_by_2( A)

    ! Calculate x
    x( 1) = Ainv( 1,1) * b( 1) + Ainv( 1,2) * b( 2)
    x( 2) = Ainv( 2,1) * b( 1) + Ainv( 2,2) * b( 2)

  END FUNCTION solve_Axb_2_by_2

! == Some basic geometry

  PURE FUNCTION cross2( a,b) RESULT( z)
    ! Vector product z between 2-dimensional vectors a and b

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(2)                             , INTENT(IN)    :: a, b

    ! Output variables:
    REAL(dp)                                                           :: z

    z = (a( 1) * b( 2)) - (a( 2) * b( 1))

  END FUNCTION cross2

  PURE SUBROUTINE segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
    ! Find out if the line segments [pq] and [rs] intersect. If so, return
    ! the coordinates of the point of intersection

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(2)                             , INTENT(IN)    :: p, q, r, s
    REAL(dp)                                           , INTENT(IN)    :: tol_dist

    ! Output variables
    REAL(dp), DIMENSION(2)                             , INTENT(OUT)   :: llis
    LOGICAL                                            , INTENT(OUT)   :: do_cross

    ! Local variables:
    REAL(dp), DIMENSION(2,2)                                           :: A
    REAL(dp), DIMENSION(2)                                             :: x, b

    ! If pq and rs are colinear, define them as not intersecting
    IF ( ABS( cross2( [q(1)-p(1), q(2)-p(2)], [s(1)-r(1), s(2)-r(2)] )) < tol_dist) THEN
      llis = [0._dp, 0._dp]
      do_cross = .FALSE.
      RETURN
    END IF

    A( 1,:) = [ (p( 1) - q( 1)), (r( 1) - s( 1))]
    A( 2,:) = [ (p( 2) - q( 2)), (r( 2) - s( 2))]
    b = [ (r( 1) - q( 1)), (r( 2) - q( 2))]

    ! Solve for x
    x = solve_Axb_2_by_2( A,b)

    llis = [q( 1) + x( 1) * (p( 1) - q( 1)), q( 2) + x( 1) * (p( 2) - q( 2))]

    IF (x( 1) > 0._dp .AND. x( 1) < 1._dp .AND. x( 2) > 0._dp .AND. x( 2) < 1._dp) THEN
      do_cross = .TRUE.
    ELSE
      do_cross = .FALSE.
    END IF

  END SUBROUTINE segment_intersection

  PURE FUNCTION is_in_polygon( Vpoly, p) RESULT( isso)
    ! Use the ray-casting algorithm to check if the point p = [px,py]
    ! lies inside the polygon spanned by poly = [x1,y1; x2,y2; ...; xn,yn]

    IMPLICIT NONE

    ! Input variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    )                       , INTENT(IN)    :: Vpoly
    REAL(dp), DIMENSION(2)                             , INTENT(IN)    :: p

    ! Output variables:
    LOGICAL                                                            :: isso

    ! Local variables:
    REAL(dp), DIMENSION(2)                                             :: q,r,s,llis
    REAL(dp)                                                           :: xmin,xmax,ymin,ymax
    INTEGER                                                            :: n_intersects
    INTEGER                                                            :: vi,vj,n_vertices
    LOGICAL                                                            :: do_cross
    REAL(dp), PARAMETER                                                :: tol_dist = 1E-5_dp

    isso = .FALSE.

    xmin = MINVAL( Vpoly(:,1))
    xmax = MAXVAL( Vpoly(:,1))
    ymin = MINVAL( Vpoly(:,2))
    ymax = MAXVAL( Vpoly(:,2))

    ! Quick test
    IF (p(1) < xmin .OR. p(1) > xmax .OR. &
        p(2) < ymin .OR. p(2) > ymax) THEN
      isso = .FALSE.
      RETURN
    END IF

    ! Define the endpoint of the east-pointing ray
    q = [xmax + (xmax - xmin) / 10._dp, p(2)]

    ! Determine how often the ray intersects the polygon

    n_vertices   = SIZE( Vpoly,1)
    n_intersects = 0

    DO vi = 1, n_vertices

      ! Find vertices spanning a polygon line section
      IF (vi < n_vertices) THEN
        vj = vi + 1
      ELSE
        vj = 1
      END IF

      ! Define the line section
      r = Vpoly( vi,:)
      s = Vpoly( vj,:)

      ! Determine if the ray intersects the line section
      IF ((r(2) < p(2) .AND. s(2) < p(2)) .OR. (r(2) > p(2) .AND. s(2) > p(2))) THEN
        do_cross = .FALSE.
      ELSE
        CALL segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
      END IF

      IF (do_cross) n_intersects = n_intersects + 1

    END DO ! DO vi = 1, n_vertices

    ! If the number of intersections is odd, p lies inside the polygon
    IF (MOD( n_intersects,2) == 1) THEN
      isso = .TRUE.
    ELSE
      isso = .FALSE.
    END IF

  END FUNCTION is_in_polygon

  PURE FUNCTION lies_on_line_segment( pa, pb, pc, tol_dist) RESULT(isso)
    ! Test if the point pc lies within tol_dist of the line pa-pb

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: pa, pb, pc
    REAL(dp),                 INTENT(IN)          :: tol_dist
    LOGICAL                                       :: isso

    ! Local variables:
    REAL(dp), DIMENSION(2)                        :: d, e, d_norm, e_par, e_ort

    d = pb - pa
    e = pc - pa

    d_norm = d / NORM2( d)

    e_par = (e( 1) * d_norm( 1) + e( 2) * d_norm( 2)) * d_norm
    e_ort = e - e_par

    isso = .TRUE.
    IF (NORM2( e_ort) > tol_dist) THEN
      isso = .FALSE.
      RETURN
    END IF

    IF ((e( 1) * d( 1) + e( 2)* d ( 2)) > 0._dp) THEN
      IF (NORM2( e_par) > (NORM2( d))) THEN
        isso = .FALSE.
        RETURN
      END IF
    ELSE
      IF (NORM2( e_par) > 0._dp) THEN
        isso = .FALSE.
        RETURN
      END IF
    END IF

  END FUNCTION lies_on_line_segment

  PURE SUBROUTINE line_from_points( p, q, la, lb, lc)
    ! Find a,b,c such that the line ax + by = c passes through p and q

    IMPLICIT NONE

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q
    REAL(dp),                   INTENT(OUT)       :: la, lb, lc

    la = q( 2) - p( 2)
    lb = p( 1) - q( 1)
    lc = la * (p( 1))+ lb * (p( 2))

  END SUBROUTINE line_from_points

  PURE SUBROUTINE perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    ! Find a,b,c such that the line ax + by = c describes the perpendicular
    ! bisector to the line [pq]

    IMPLICIT NONE

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q
    REAL(dp),                   INTENT(IN)        :: la1, lb1
    REAL(dp),                   INTENT(OUT)       :: la2, lb2, lc2
    REAL(dp)                                      :: temp
    REAL(dp), DIMENSION(2)                        :: m

    m = (p+q)/2
    lc2 = -lb1*m(1) + la1*m(2)

    temp = la1
    la2 = -lb1
    lb2 = temp

  END SUBROUTINE perpendicular_bisector_from_line

  PURE SUBROUTINE line_line_intersection( la1, lb1, lc1, la2, lb2, lc2, llis)
    ! Find the intersection llis of the lines la1*x+lb1*y=lc1 and la2*x+lb2*y=lc2

    IMPLICIT NONE

    REAL(dp),                   INTENT(IN)        :: la1, lb1, lc1, la2, lb2, lc2
    REAL(dp), DIMENSION(2),     INTENT(OUT)       :: llis
    REAL(dp)                                      :: d

    d = la1*lb2 - la2*lb1
    IF (d == 0) THEN
      ! The lines are parallel.
      llis = [1E30, 1E30]
    ELSE
      llis = [(lb2*lc1 - lb1*lc2), (la1*lc2 - la2*lc1)]/d
    END IF

  END SUBROUTINE line_line_intersection

  PURE FUNCTION circumcenter( p, q, r) RESULT( cc)
    ! Find the circumcenter cc of the triangle pqr

    IMPLICIT NONE

    ! Some basic vector operations
    ! Find the circumcenter cc of the triangle pqr
    ! If pqr are colinear, returns [1e30,1e30]

    REAL(dp), DIMENSION(2),     INTENT(IN)        :: p, q, r
    REAL(dp), DIMENSION(2)                        :: cc
    REAL(dp)                                      :: la1,lb1,lc1,le1,lf1,lg1
    REAL(dp)                                      :: la2,lb2,lc2,le2,lf2,lg2

    cc = [0._dp, 0._dp]

    ! Line PQ is represented as ax + by = c, Line QR is represented as ex + fy = g
    CALL line_from_points( p, q, la1, lb1, lc1)
    CALL line_from_points( q, r, le1, lf1, lg1)

    ! Converting lines PQ and QR to perpendicular
    ! bisectors. After this, L = ax + by = c
    ! M = ex + fy = g
    CALL perpendicular_bisector_from_line( p, q, la1, lb1, la2, lb2, lc2)
    CALL perpendicular_bisector_from_line( q, r, le1, lf1, le2, lf2, lg2)

    ! The point of intersection of L and M gives
    ! the circumcenter
    CALL line_line_intersection( la2, lb2, lc2, le2, lf2, lg2, cc)

  END FUNCTION circumcenter

  PURE FUNCTION geometric_center( p, q, r) RESULT( gc)
    ! Calculate the geometric centre of triangle [pqr]

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2), INTENT(IN)            :: p, q, r
    REAL(dp), DIMENSION(2)                        :: gc

    gc = (p + q + r) / 3._dp

  END FUNCTION geometric_center

  PURE FUNCTION triangle_area( p, q, r) RESULT( A)
    ! Find the area of the triangle [p,q,r]

    IMPLICIT NONE

    REAL(dp), DIMENSION(2), INTENT(IN)  :: p, q, r
    REAL(dp)                            :: A

    A = ABS( cross2( [q(1)-p(1), q(2)-p(2)], [r(1)-p(1), r(2)-p(2)] )) / 2._dp

  END FUNCTION triangle_area

  PURE FUNCTION is_in_triangle( pa, pb, pc, p) RESULT(isso)
    ! Check if the point p lies inside the triangle abc, or within distance tol of its edges

    IMPLICIT NONE

    REAL(dp), DIMENSION(2),   INTENT(IN)          :: pa, pb, pc, p
    LOGICAL                                       :: isso
    REAL(dp)                                      :: as_x, as_y, s1, s2, s3
    REAL(dp), PARAMETER                           :: tol = 1E-8_dp

    as_x = p( 1) - pa( 1)
    as_y = p( 2) - pa( 2)

    s1 = ((pb( 1) - pa( 1)) * as_y - (pb( 2) - pa( 2)) * as_x)
    s2 = ((pc( 1) - pa( 1)) * as_y - (pc( 2) - pa( 2)) * as_x)
    s3 = ((pc( 1) - pb( 1)) * (p( 2) - pb( 2)) - (pc( 2) - pb( 2)) * (p( 1) - pb( 1)))

    isso = .FALSE.

    IF (s1 > -tol .AND. s2 < tol .AND. s3 > -tol) THEN
      isso = .TRUE.
      RETURN
    END IF

  END FUNCTION is_in_triangle

  PURE FUNCTION longest_triangle_leg( p, q, r) RESULT( d)
    ! Find the longest leg of the triangle [p,q,r]

    IMPLICIT NONE

    REAL(dp), DIMENSION(2), INTENT(IN)  :: p, q, r
    REAL(dp)                            :: d
    REAL(dp)                            :: d_pq, d_qr, d_rp

    d_pq = NORM2( p-q)
    d_qr = NORM2( q-r)
    d_rp = NORM2( r-p)
    d = MAX( MAX( d_pq, d_qr), d_rp)

  END FUNCTION longest_triangle_leg

  PURE FUNCTION smallest_triangle_angle( p, q, r) RESULT( alpha)
    ! Find the smallest internal angle of the triangle [p,q,r]

    IMPLICIT NONE

    REAL(dp), DIMENSION(2), INTENT(IN)  :: p, q, r
    REAL(dp)                            :: alpha
    REAL(dp), DIMENSION(2)              :: pq, qr, rp
    REAL(dp)                            :: ap, aq, ar

    ! Triangle legs
    pq = p-q
    qr = q-r
    rp = r-p

    ! Internal angles
    ap = ACOS(-(rp( 1) * pq( 1) + rp( 2) * pq( 2)) / (NORM2( rp) * NORM2( pq)))
    aq = ACOS(-(pq( 1) * qr( 1) + pq( 2) * qr( 2)) / (NORM2( pq) * NORM2( qr)))
    ar = ACOS(-(rp( 1) * qr( 1) + rp( 2) * qr( 2)) / (NORM2( rp) * NORM2( qr)))

    ! Smallest internal angle
    alpha = MIN( MIN( ap, aq), ar)

  END FUNCTION smallest_triangle_angle

  SUBROUTINE crop_line_to_domain( p, q, xmin, xmax, ymin, ymax, tol_dist, pp, qq, is_valid_line)
    ! Crop the line [pq] so that it lies within the specified domain;
    ! if [pq] doesn't pass through the domain at all, return is_valid_line = .FALSE.

    IMPLICIT NONE

    ! In/output variables
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p, q
    REAL(dp),                            INTENT(IN)    :: xmin, xmax, ymin, ymax, tol_dist
    REAL(dp), DIMENSION(2),              INTENT(OUT)   :: pp, qq
    LOGICAL,                             INTENT(OUT)   :: is_valid_line

    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: sw, se, nw, ne, llis
    INTEGER                                            :: edge_index_p, edge_index_q
    LOGICAL                                            :: do_cross

    pp = p
    qq = q
    is_valid_line = .TRUE.

    sw = [xmin,ymin]
    se = [xmax,ymin]
    nw = [xmin,ymax]
    ne = [xmax,ymax]

    ! Determine in which quadrants p and q lie
    ! (same as with edge_index; 1-8 clockwise starting north, 0 means inside

    IF     (pp(1) >= xmin .AND. pp(1) <= xmax .AND. pp(2) > ymax) THEN
      ! North
      edge_index_p = 1
    ELSEIF (pp(1) > xmax .AND. pp(2) > ymax) THEN
      ! Northeast
      edge_index_p = 2
    ELSEIF (pp(1) > xmax .AND. pp(2) >= ymin .AND. pp(2) <= ymax) THEN
      ! East
      edge_index_p = 3
    ELSEIF (pp(1) > xmax .AND. pp(2) < ymin) THEN
      ! Southeast
      edge_index_p = 4
    ELSEIF (pp(1) >= xmin .AND. pp(1) <= xmax .AND. pp(2) < ymin) THEN
      ! South
      edge_index_p = 5
    ELSEIF (pp(1) < xmin .AND. pp(2) < ymin) THEN
      ! Southwest
      edge_index_p = 6
    ELSEIF (pp(1) < xmin .AND. pp(2) >= ymin .AND. pp(2) <= ymax) THEN
      ! West
      edge_index_p = 7
    ELSEIF (pp(1) < xmin .AND. pp(2) > ymax) THEN
      ! Northwest
      edge_index_p = 8
    ELSE
      ! Inside the mesh domain
      edge_index_p = 0
    END IF

    IF     (qq(1) >= xmin .AND. qq(1) <= xmax .AND. qq(2) > ymax) THEN
      ! North
      edge_index_q = 1
    ELSEIF (qq(1) > xmax .AND. qq(2) > ymax) THEN
      ! Northeast
      edge_index_q = 2
    ELSEIF (qq(1) > xmax .AND. qq(2) >= ymin .AND. qq(2) <= ymax) THEN
      ! East
      edge_index_q = 3
    ELSEIF (qq(1) > xmax .AND. qq(2) < ymin) THEN
      ! Southeast
      edge_index_q = 4
    ELSEIF (qq(1) >= xmin .AND. qq(1) <= xmax .AND. qq(2) < ymin) THEN
      ! South
      edge_index_q = 5
    ELSEIF (qq(1) < xmin .AND. qq(2) < ymin) THEN
      ! Southwest
      edge_index_q = 6
    ELSEIF (qq(1) < xmin .AND. qq(2) >= ymin .AND. qq(2) <= ymax) THEN
      ! West
      edge_index_q = 7
    ELSEIF (qq(1) < xmin .AND. qq(2) > ymax) THEN
      ! Northwest
      edge_index_q = 8
    ELSE
      ! Inside the mesh domain
      edge_index_q = 0
    END IF

    IF (edge_index_p == 0 .AND. edge_index_q == 0) THEN
      ! Both p and q lie inside the mesh domain
      RETURN
    END IF

    IF (edge_index_p == 0 .AND. edge_index_q > 0) THEN
      ! p lies inside the mesh domain, q lies outside

      ! Check IF [pq] passes through any of the four corners
      IF     (lies_on_line_segment( pp, qq, sw, tol_dist)) THEN
        ! [pq] passes through the southwest corner of the mesh
        qq = sw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, se, tol_dist)) THEN
        ! [pq] passes through the southeast corner of the mesh
        qq = se
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, nw, tol_dist)) THEN
        ! [pq] passes through the northwest corner of the mesh
        qq = nw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, ne, tol_dist)) THEN
        ! [pq] passes through the northeast corner of the mesh
        qq = ne
        RETURN
      END IF

      ! Check IF [pq] crosses any of the four borders

      ! South
      CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the southern border
        qq = llis
        RETURN
      END IF

      ! West
      CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the western border
        qq = llis
        RETURN
      END IF

      ! North
      CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the northern border
        qq = llis
        RETURN
      END IF

      ! East
      CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the eastern border
        qq = llis
        RETURN
      END IF

      ! This point should be unreachable
      CALL crash('crop_line_to_mesh_domain - reached the unreachable point (p inside, q outside)!')

    END IF ! IF (edge_index_p == 0 .AND. edge_index_q > 0)

    IF (edge_index_q == 0 .AND. edge_index_p > 0) THEN
      ! q lies inside the mesh domain, p lies outside

      ! Check IF [pq] passes through any of the four corners
      IF     (lies_on_line_segment( pp, qq, sw, tol_dist)) THEN
        ! [pq] passes through the southwest corner of the mesh
        pp = sw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, se, tol_dist)) THEN
        ! [pq] passes through the southeast corner of the mesh
        pp = se
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, nw, tol_dist)) THEN
        ! [pq] passes through the northwest corner of the mesh
        pp = nw
        RETURN
      ELSEIF (lies_on_line_segment( pp, qq, ne, tol_dist)) THEN
        ! [pq] passes through the northeast corner of the mesh
        pp = ne
        RETURN
      END IF

      ! Check IF [pq] crosses any of the four borders

      ! South
      CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the southern border
        pp = llis
        RETURN
      END IF

      ! West
      CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the western border
        pp = llis
        RETURN
      END IF

      ! North
      CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the northern border
        pp = llis
        RETURN
      END IF

      ! East
      CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
      IF (do_cross) THEN
        ! [pq] crosses the eastern border
        pp = llis
        RETURN
      END IF

      ! This point should be unreachable
      CALL crash('crop_line_to_mesh_domain - reached the unreachable point (q inside, p outside)!')

    END IF ! IF (edge_index_q == 0 .AND. edge_index_p > 0)

    ! Both p and q lie outside the mesh domain

    IF     (pp(1) < xmin .AND. qq(1) < xmin) THEN
      ! Both p and q lie west of the western mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    ELSEIF (pp(1) > xmax .AND. qq(1) > xmax) THEN
      ! Both p and q lie east of the eastern mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    ELSEIF (pp(2) < ymin .AND. qq(2) < ymin) THEN
      ! Both p and q lie south of the southern mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    ELSEIF (pp(2) > ymax .AND. qq(2) > ymax) THEN
      ! Both p and q lie north of the northern mesh border; [pq] cannot pass through the mesh domain
      is_valid_line = .FALSE.
      RETURN
    END IF

    IF (edge_index_p == 1) THEN
      ! p lies in the northern quadrant

      IF (edge_index_q == 3) THEN
        ! q lies in the eastern quadrant; check IF [pq] cuts through the northeast corner

        IF (cross2( (ne - qq), (pp - qq)) > 0) THEN
          ! [pq] cuts through the northeast corner
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, ne, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect ne-se, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 7) THEN
        ! q lies in the western quadrant; check IF [pq] cuts through the northwest corner

        IF (cross2( (pp - qq), (ne - qq)) > 0) THEN
          ! [pq] cuts through the northeast corner
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, nw, sw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-sw, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    ELSEIF (edge_index_p == 3) THEN
      ! p lies in the eastern quadrant

      IF (edge_index_q == 1) THEN
        ! q lies in the northern quadrant; check IF [pq] cuts through the northeast corner

        IF (cross2( (ne - pp), (qq - pp)) > 0) THEN
          ! [pq] cuts through the northeast corner
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          qq = llis
          CALL segment_intersection( pp, qq, ne, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect ne-se, but it doesnt!')
          END IF
          pp = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 5) THEN
        ! q lies in the southern quadrant; check IF [pq] cuts through the southeast corner

        IF (cross2( (se - qq), (pp - qq)) > 0) THEN
          ! [pq] cuts through the southeast corner
          CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect se-ne, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    ELSEIF (edge_index_p == 5) THEN
      ! p lies in the southern quadrant

      IF (edge_index_q == 3) THEN
        ! q lies in the eastern quadrant; check IF [pq] cuts through the southeast corner

        IF (cross2( (se - pp), (qq - pp)) > 0) THEN
          ! [pq] cuts through the southwest corner
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, se, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect se-ne, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 7) THEN
        ! q lies in the western quadrant; check IF [pq] cuts through the southwest corner

        IF (cross2( (qq - pp), (sw - pp)) > 0) THEN
          ! [pq] cuts through the southwest corner
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-nw, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    ELSEIF (edge_index_p == 7) THEN
      ! p lies in the western quadrant

      IF (edge_index_q == 5) THEN
        ! q lies in the southern quadrant; check IF [pq] cuts through the southwest corner

        IF (cross2( (pp - qq), (sw - qq)) > 0) THEN
          ! [pq] cuts through the southwest corner
          CALL segment_intersection( pp, qq, sw, se, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-se, but it doesnt!')
          END IF
          qq = llis
          CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-nw, but it doesnt!')
          END IF
          pp = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSEIF (edge_index_q == 1) THEN
        ! q lies in the northern quadrant; check IF [pq] cuts through the northwest corner

        IF (cross2( (qq - pp), (nw - pp)) > 0) THEN
          ! [pq] cuts through the northwest corner
          CALL segment_intersection( pp, qq, sw, nw, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect sw-nw, but it doesnt!')
          END IF
          pp = llis
          CALL segment_intersection( pp, qq, nw, ne, llis, do_cross, tol_dist)
          IF (.NOT. do_cross) THEN
            CALL crash('crop_line_to_mesh_domain - [pq] should intersect nw-ne, but it doesnt!')
          END IF
          qq = llis
          RETURN
        ELSE
          ! [pq] does not pass through the mesh domain
          is_valid_line = .FALSE.
          RETURN
        END IF

      ELSE
        CALL crash('crop_line_to_mesh_domain - edge_index_p = {int_01}, edge_index_q = {int_02}', int_01 = edge_index_p, int_02 = edge_index_q)
      END IF

    END IF ! IF (edge_index_p == 1)

    ! Don't consider the other options...
    is_valid_line = .FALSE.

  END SUBROUTINE crop_line_to_domain

  PURE FUNCTION encroaches_upon( pa, pb, pc) RESULT( isso)
    ! Check if c encroaches upon segment ab

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2), INTENT(IN)            :: pa, pb, pc
    LOGICAL                                       :: isso

    isso = (NORM2( pc - (pa + pb) / 2._dp) < NORM2( pa - pb) / 2._dp)

  END FUNCTION encroaches_upon

! == Calculate contour lines for mesh generation from gridded/meshed data

  SUBROUTINE calc_grid_contour_as_line( x, y, d, f, line)
    ! Calculate a contour line at level f for data d on a square grid%
    ! Generate the contour line in UFEMISM line format (i.e. unordered
    ! individual line segments).

    IMPLICIT NONE

    ! In/output variables
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: x, y
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
    REAL(dp),                                INTENT(IN)    :: f
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE, INTENT(OUT)   :: line

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                          :: routine_name = 'mesh_add_smileyface'
    REAL(dp), PARAMETER                                    :: tol = 1E-5_dp
    INTEGER                                                :: nx, ny
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE                :: d_scaled
    REAL(dp)                                               :: d_scaled_min, d_scaled_max
    INTEGER                                                :: i,j
    INTEGER                                                :: n_max, n
    REAL(dp)                                               :: d_sw, d_nw, d_se, d_ne
    REAL(dp)                                               :: xw, xe, xs, xn, yw, ye, ys, yn
    LOGICAL                                                :: do_cross_w, do_cross_e, do_cross_s, do_cross_n

    ! Add routine to path
    CALL init_routine( routine_name)

    nx = SIZE( x,1)
    ny = SIZE( y,1)

    ! Safety
    IF (SIZE( d,1) /= nx .OR. SIZE( d,2) /= ny) CALL crash('d is not nx-by-ny!')

    ! Trivial case: if all values of d are greater or smaller than f,
    ! the contour line is empty
    IF (MINVAL( d) >= f .OR. MAXVAL( d) <= f) THEN
      ALLOCATE( line( 0,0))
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Scale d so that it runs from -1 to 1, with 0 representing f
    d_scaled = d - f
    d_scaled_min = MINVAL( d_scaled)
    d_scaled_max = MAXVAL( d_scaled)
    DO i = 1, nx
      DO j = 1, ny
        IF (d_scaled( i,j) > 0._dp) THEN
          d_scaled( i,j) = d_scaled( i,j) /           d_scaled_max
        ELSE
          d_scaled( i,j) = d_scaled( i,j) / (-1._dp * d_scaled_min)
        END IF
      END DO
    END DO

    n_max = 1000
    ALLOCATE( line( n_max, 4))

    n = 0
    DO i = 1, nx-1
      DO j = 1, ny-1

        ! Extend allocated memory IF needed
        IF (n > n_max - 10) THEN
          n_max = n + 1000
          CALL reallocate( line, n_max, 4)
        END IF

        ! The four corners of the b-grid cell
        d_sw = d_scaled( i  ,j  )
        d_nw = d_scaled( i  ,j+1)
        d_se = d_scaled( i+1,j  )
        d_ne = d_scaled( i+1,j+1)

        xw = x( i  )
        xe = x( i+1)
        ys = y( j  )
        yn = y( j+1)

        ! If all four corners are above/below the level, no line here
        IF ((d_sw >= 0._dp .AND. d_nw >= 0._dp .AND. d_se >= 0._dp .AND. d_ne >= 0._dp) .OR. &
            (d_sw <= 0._dp .AND. d_nw <= 0._dp .AND. d_se <= 0._dp .AND. d_ne <= 0._dp)) THEN
          CYCLE
        END IF

        ! Add tolerances to keep line lengths finite
        IF (d_sw >= 0._dp) THEN
          d_sw = MAX( d_sw, tol)
        ELSE
          d_sw = MIN( d_sw, tol)
        END IF
        IF (d_nw >= 0._dp) THEN
          d_nw = MAX( d_nw, tol)
        ELSE
          d_nw = MIN( d_nw, tol)
        END IF
        IF (d_se >= 0._dp) THEN
          d_se = MAX( d_se, tol)
        ELSE
          d_se = MIN( d_se, tol)
        END IF
        IF (d_ne >= 0._dp) THEN
          d_ne = MAX( d_ne, tol)
        ELSE
          d_ne = MIN( d_ne, tol)
        END IF

        ! Find boundary crossings

        do_cross_w = .FALSE.
        do_cross_e = .FALSE.
        do_cross_s = .FALSE.
        do_cross_n = .FALSE.

        IF (d_sw * d_nw < 0._dp) THEN
          ! The contour crosses the western boundary of this b-grid cell
          do_cross_w = .TRUE.
          yw = linint_points( ys, yn, d_sw, d_nw, 0._dp)
        END IF

        IF (d_se * d_ne < 0._dp) THEN
          ! The contour crosses the eastern boundary of this b-grid cell
          do_cross_e = .TRUE.
          ye = linint_points( ys, yn, d_se, d_ne, 0._dp)
        END IF

        IF (d_nw * d_ne < 0._dp) THEN
          ! The contour crosses the northern boundary of this b-grid cell
          do_cross_n = .TRUE.
          xn = linint_points( xw, xe, d_nw, d_ne, 0._dp)
        END IF

        IF (d_sw * d_se < 0._dp) THEN
          ! The contour crosses the southern boundary of this b-grid cell
          do_cross_s = .TRUE.
          xs = linint_points( xw, xe, d_sw, d_se, 0._dp)
        END IF

        ! Add line segments
        n = n + 1
        IF     (do_cross_w) THEN
          IF     (do_cross_e) THEN
            ! From west to east
            line( n,:) = [xw,yw,xe,ye]
          ELSEIF (do_cross_s) THEN
            ! From west to south
            line( n,:) = [xw,yw,xs,ys]
          ELSEIF (do_cross_n) THEN
            ! From west to north
            line( n,:) = [xw,yw,xn,yn]
          ELSE
            CALL crash('found only a crossing at the western boundary!')
          END IF
        ELSEIF (do_cross_e) THEN
          IF     (do_cross_s) THEN
            ! From east to south
            line( n,:) = [xe,ye,xs,ys]
          ELSEIF (do_cross_n) THEN
            ! From east to north
            line( n,:) = [xe,ye,xn,yn]
          ELSE
            CALL crash('found only a crossing at the eastern boundary!')
          END IF
        ELSEIF (do_cross_s) THEN
          IF     (do_cross_n) THEN
            ! From south to north
            line( n,:) = [xs,ys,xn,yn]
          ELSE
            CALL crash('found only a crossing at the southern boundary!')
          END IF
        ELSEIF (do_cross_n) THEN
            CALL crash('found only a crossing at the northern boundary!')
        ELSE
          CALL crash('whaa!')
        END IF

      END DO
    END DO

    ! Crop memory
    CALL reallocate( line, n, 4)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_grid_contour_as_line

  PURE FUNCTION linint_points( x1, x2, f1, f2, f0) RESULT( x0)
    ! Given a function f( x) and points x1, x2 such that f( x1) = f1, f( x2) = f2,
    ! interpolate f linearly to find the point x0 such that f( x0) = f0

    IMPLICIT NONE

    REAL(dp),               INTENT(IN)  :: x1, x2, f1, f2, f0
    REAL(dp)                            :: x0
    REAL(dp)                            :: lambda

    lambda = (f2 - f1) / (x2 - x1);
    x0 = x1 + (f0 - f1) / lambda;

  END FUNCTION linint_points

! == Debugging

  SUBROUTINE check_for_NaN_dp_0D(  d, d_name)
    ! Check if NaN values occur in the 0-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    IF     ( ISNAN( d)) THEN
      CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '"')
    ELSEIF (d /= d) THEN
      CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '"')
    ELSEIF (d > HUGE( d)) THEN
      CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '"')
    ELSEIF (d  < -HUGE( d)) THEN
      CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '"')
    END IF

  END SUBROUTINE check_for_NaN_dp_0D

  SUBROUTINE check_for_NaN_dp_1D(  d, d_name)
    ! Check if NaN values occur in the 1-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:      )                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)

      IF     ( ISNAN( d( i))) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      ELSEIF (d( i) /= d( i)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      ELSEIF (d( i)  < -HUGE( d( i))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      END IF

    END DO

  END SUBROUTINE check_for_NaN_dp_1D

  SUBROUTINE check_for_NaN_dp_2D(  d, d_name)
    ! Check if NaN values occur in the 2-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:    )                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i,j

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)
    DO j = 1, SIZE( d,2)

      IF     ( ISNAN( d( i,j))) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      ELSEIF (d( i,j) /= d( i,j)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      ELSEIF (d( i,j) > HUGE( d( i,j))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      ELSEIF (d( i,j)  < -HUGE( d( i,j))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      END IF

    END DO
    END DO

  END SUBROUTINE check_for_NaN_dp_2D

  SUBROUTINE check_for_NaN_dp_3D(  d, d_name)
    ! Check if NaN values occur in the 3-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:  )                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i,j,k

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)
    DO j = 1, SIZE( d,2)
    DO k = 1, SIZE( d,3)

      IF     ( ISNAN( d( i,j,k))) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      ELSEIF (d( i,j,k) /= d( i,j,k)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      ELSEIF (d( i,j,k) > HUGE( d( i,j,k))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      ELSEIF (d( i,j,k)  < -HUGE( d( i,j,k))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE check_for_NaN_dp_3D

  SUBROUTINE check_for_NaN_dp_4D(  d, d_name)
    ! Check if NaN values occur in the 4-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp), DIMENSION(:,:,:,:)                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i,j,k,r

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)
    DO j = 1, SIZE( d,2)
    DO k = 1, SIZE( d,3)
    DO r = 1, SIZE( d,4)

      IF     ( ISNAN( d( i,j,k,r))) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      ELSEIF (d( i,j,k,r) /= d( i,j,k,r)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      ELSEIF (d( i,j,k,r) > HUGE( d( i,j,k,r))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      ELSEIF (d( i,j,k,r)  < -HUGE( d( i,j,k,r))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      END IF

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE check_for_NaN_dp_4D

  SUBROUTINE check_for_NaN_int_0D(  d, d_name)
    ! Check if NaN values occur in the 0-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    REAL(dp)                                           , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    IF     (d /= d) THEN
      CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '"')
    ELSEIF (d > HUGE( d)) THEN
      CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '"')
    ELSEIF (d  < -HUGE( d)) THEN
      CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '"')
    END IF

  END SUBROUTINE check_for_NaN_int_0D

  SUBROUTINE check_for_NaN_int_1D(  d, d_name)
    ! Check if NaN values occur in the 1-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:      )                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)

      IF     (d( i) /= d( i)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      ELSEIF (d( i)  < -HUGE( d( i))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]', int_01 = i)
      END IF

    END DO

  END SUBROUTINE check_for_NaN_int_1D

  SUBROUTINE check_for_NaN_int_2D(  d, d_name)
    ! Check if NaN values occur in the 2-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:    )                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i,j

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)
    DO j = 1, SIZE( d,2)

      IF     (d( i,j) /= d( i,j)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      ELSEIF (d( i,j) > HUGE( d( i,j))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      ELSEIF (d( i,j)  < -HUGE( d( i,j))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]', int_01 = i, int_02 = j)
      END IF

    END DO
    END DO

  END SUBROUTINE check_for_NaN_int_2D

  SUBROUTINE check_for_NaN_int_3D(  d, d_name)
    ! Check if NaN values occur in the 3-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:  )                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i,j,k

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)
    DO j = 1, SIZE( d,2)
    DO k = 1, SIZE( d,3)

      IF     (d( i,j,k) /= d( i,j,k)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      ELSEIF (d( i,j,k) > HUGE( d( i,j,k))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      ELSEIF (d( i,j,k)  < -HUGE( d( i,j,k))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]', int_01 = i, int_02 = j, int_03 = k)
      END IF

    END DO
    END DO
    END DO

  END SUBROUTINE check_for_NaN_int_3D

  SUBROUTINE check_for_NaN_int_4D(  d, d_name)
    ! Check if NaN values occur in the 4-D dp data field d

    IMPLICIT NONE

    ! Inpit variables:
!   REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER , DIMENSION(:,:,:,:)                       , INTENT(IN)    :: d
    CHARACTER(LEN=*),                          OPTIONAL, INTENT(IN)    :: d_name

    ! Local variables:
    CHARACTER(LEN=256)                                                 :: d_name_loc
    INTEGER                                                            :: i,j,k,r

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = 1, SIZE( d,1)
    DO j = 1, SIZE( d,2)
    DO k = 1, SIZE( d,3)
    DO r = 1, SIZE( d,4)

      IF     (d( i,j,k,r) /= d( i,j,k,r)) THEN
        CALL crash( 'detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      ELSEIF (d( i,j,k,r) > HUGE( d( i,j,k,r))) THEN
        CALL crash( 'detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      ELSEIF (d( i,j,k,r)  < -HUGE( d( i,j,k,r))) THEN
        CALL crash( 'detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03},{int_04}]', int_01 = i, int_02 = j, int_03 = k, int_04 = r)
      END IF

    END DO
    END DO
    END DO
    END DO

  END SUBROUTINE check_for_NaN_int_4D

END MODULE math_utilities
