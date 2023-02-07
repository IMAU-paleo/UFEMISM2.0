MODULE utilities
  ! Some generally useful tools, and basic mathematical functions

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE parameters

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

! == Floatation criterion, surface elevation, and thickness above floatation

  FUNCTION is_floating( Hi, Hb, SL) RESULT( isso)
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

  FUNCTION ice_surface_elevation( Hi, Hb, SL) RESULT( Hs)
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

  FUNCTION thickness_above_floatation( Hi, Hb, SL) RESULT( TAF)
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

  SUBROUTINE error_function(X, ERR)
    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: X

    ! Output variables:
    REAL(dp), INTENT(OUT) :: ERR

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
  END SUBROUTINE error_function

! == The oblique stereographic projection

  SUBROUTINE oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, beta_deg, x, y, k_P)
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

    ! See equation (21-4) on page 157 in Snyder (1987):
    IF(PRESENT(k_P)) k_P = (1._dp + COS(alpha)) / (1._dp + SIN(phi_M) * SIN(phi_P) + COS(phi_M) * COS(phi_P) * COS(lambda_P - lambda_M))

  END SUBROUTINE oblique_sg_projection

  SUBROUTINE inverse_oblique_sg_projection( x, y, lambda_M_deg, phi_M_deg, beta_deg, lambda_P, phi_P)
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

!! == Line integrals used in conservative remapping
!
!  SUBROUTINE line_integral_xdy(   p, q, tol_dist, I_pq)
!    ! Calculate the line integral x dy from p to q
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
!    REAL(dp),                                INTENT(IN)    :: tol_dist
!    REAL(dp),                                INTENT(OUT)   :: I_pq
!
!    ! Local variables:
!    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
!
!    xp = p(1)
!    yp = p(2)
!    xq = q(1)
!    yq = q(2)
!
!    IF (ABS(yp-yq) < tol_dist) THEN
!      I_pq = 0._dp
!      RETURN
!    END IF
!
!    dx = q(1)-p(1)
!    dy = q(2)-p(2)
!
!    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)
!
!  END SUBROUTINE line_integral_xdy
!
!  SUBROUTINE line_integral_mxydx( p, q, tol_dist, I_pq)
!    ! Calculate the line integral -xy dx from p to q
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
!    REAL(dp),                                INTENT(IN)    :: tol_dist
!    REAL(dp),                                INTENT(OUT)   :: I_pq
!
!    ! Local variables:
!    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
!
!    xp = p(1)
!    yp = p(2)
!    xq = q(1)
!    yq = q(2)
!
!    IF (ABS(xp-xq) < tol_dist) THEN
!      I_pq = 0._dp
!      RETURN
!    END IF
!
!    dx = q(1)-p(1)
!    dy = q(2)-p(2)
!
!    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))
!
!  END SUBROUTINE line_integral_mxydx
!
!  SUBROUTINE line_integral_xydy(  p, q, tol_dist, I_pq)
!    ! Calculate the line integral xy dy from p to q
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
!    REAL(dp),                                INTENT(IN)    :: tol_dist
!    REAL(dp),                                INTENT(OUT)   :: I_pq
!
!    ! Local variables:
!    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
!
!    xp = p(1)
!    yp = p(2)
!    xq = q(1)
!    yq = q(2)
!
!    IF (ABS(yp-yq) < tol_dist) THEN
!      I_pq = 0._dp
!      RETURN
!    END IF
!
!    dx = q(1)-p(1)
!    dy = q(2)-p(2)
!
!    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))
!
!  END SUBROUTINE line_integral_xydy
!
!! == Some wrappers for LAPACK matrix functionality
!
!  FUNCTION tridiagonal_solve( ldiag, diag, udiag, rhs) RESULT(x)
!    ! Lapack tridiagnal solver (in double precision):
!    ! Matrix system solver for tridiagonal matrices.
!    ! Used e.g. in solving the ADI scheme.
!    ! ldiag = lower diagonal elements (j,j-1) of the matrix
!    ! diag  = diagonal elements (j,j) of the matrix
!    ! udiag = upper diagonal elements (j,j+1) of the matrix
!    ! rhs   = right hand side of the matrix equation in the ADI scheme
!
!    IMPLICIT NONE
!
!    ! Input variables:
!    REAL(dp), DIMENSION(:),            INTENT(IN) :: diag
!    REAL(dp), DIMENSION(SIZE(diag)-1), INTENT(IN) :: udiag, ldiag
!    REAL(dp), DIMENSION(SIZE(diag)),   INTENT(IN) :: rhs
!
!    ! Result variable:
!    REAL(dp), DIMENSION(SIZE(diag))               :: x
!
!    ! Local variables:
!    INTEGER                                       :: info
!    REAL(dp), DIMENSION(SIZE(diag))               :: diag_copy
!    REAL(dp), DIMENSION(SIZE(udiag))              :: udiag_copy, ldiag_copy
!
!    ! The LAPACK solver will overwrite the rhs with the solution x. Therefore we
!    ! first copy the rhs in the solution vector x:
!    x = rhs
!
!    ! The LAPACK solver will change the elements in the matrix, therefore we copy them:
!    diag_copy  =  diag
!    udiag_copy = udiag
!    ldiag_copy = ldiag
!
!    CALL DGTSV(SIZE(diag), 1, ldiag_copy, diag_copy, udiag_copy, x, SIZE(diag), info)
!
!    IF (info /= 0) THEN
!      WRITE(0,*) 'tridiagonal_solve - ERROR: LAPACK solver DGTSV returned error message info = ', info
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!  END FUNCTION tridiagonal_solve
!
!! == Finding inverses of some small matrices
!
!  SUBROUTINE calc_matrix_inverse_2_by_2( A, Ainv)
!    ! Direct inversion of a 2-by-2 matrix
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    REAL(dp), DIMENSION(2,2  ),          INTENT(IN)    :: A
!    REAL(dp), DIMENSION(2,2  ),          INTENT(INOUT) :: Ainv
!
!    ! Local variables
!    REAL(dp)                                           :: detA
!
!    ! Calculate the determinant of A
!    detA = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)
!
!    ! Safety
!    IF (ABS( detA) < TINY( detA)) THEN
!      !PRINT *, A(1,1), A(1,2)
!      !PRINT *, A(2,1), A(2,2)
!      WRITE(0,*) 'calc_matrix_inverse_2_by_2 - ERROR: matrix is numerically singular!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ! Find the inverse of A
!    Ainv( 1,1) =  A( 2,2) / detA
!    Ainv( 1,2) = -A( 1,2) / detA
!    Ainv( 2,1) = -A( 2,1) / detA
!    Ainv( 2,2) =  A( 1,1) / detA
!
!  END SUBROUTINE calc_matrix_inverse_2_by_2
!
!  SUBROUTINE calc_matrix_inverse_3_by_3( A, Ainv)
!    ! Direct inversion of a 3-by-3 matrix
!    !
!    ! See: https://metric.ma.ic.ac.uk/metric_public/matrices/inverses/inverses2.html
!
!    IMPLICIT NONE
!
!    ! In- and output variables:
!    REAL(dp), DIMENSION(3,3  ),          INTENT(IN)    :: A
!    REAL(dp), DIMENSION(3,3  ),          INTENT(INOUT) :: Ainv
!
!    ! Local variables
!    REAL(dp)                                           :: detA
!
!    ! Calculate the minors of A
!    Ainv( 1,1) = A( 2,2) * A( 3,3) - A( 2,3) * A( 3,2)
!    Ainv( 1,2) = A( 2,1) * A( 3,3) - A( 2,3) * A( 3,1)
!    Ainv( 1,3) = A( 2,1) * A( 3,2) - A( 2,2) * A( 3,1)
!    Ainv( 2,1) = A( 1,2) * A( 3,3) - A( 1,3) * A( 3,2)
!    Ainv( 2,2) = A( 1,1) * A( 3,3) - A( 1,3) * A( 3,1)
!    Ainv( 2,3) = A( 1,1) * A( 3,2) - A( 1,2) * A( 3,1)
!    Ainv( 3,1) = A( 1,2) * A( 2,3) - A( 1,3) * A( 2,2)
!    Ainv( 3,2) = A( 1,1) * A( 2,3) - A( 1,3) * A( 2,1)
!    Ainv( 3,3) = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)
!
!    ! Calculate the determinant of A
!    detA = A( 1,1) * Ainv( 1,1) - A( 1,2) * Ainv( 1,2) + A( 1,3) * Ainv( 1,3)
!
!    ! Safety
!    IF (ABS( detA) < TINY( detA)) THEN
!      ! PRINT *, A(1,1), A(1,2), A(1,3)
!      ! PRINT *, A(2,1), A(2,2), A(2,3)
!      ! PRINT *, A(3,1), A(3,2), A(3,3)
!      WRITE(0,*) 'calc_matrix_inverse_3_by_3 - ERROR: matrix is numerically singular!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ! Change matrix of minors to get the matrix of cofactors
!    Ainv( 1,2) = -Ainv( 1,2)
!    Ainv( 2,1) = -Ainv( 2,1)
!    Ainv( 2,3) = -Ainv( 2,3)
!    Ainv( 3,2) = -Ainv( 3,2)
!
!    ! Transpose matrix of cofactors
!    Ainv( 1,2) = Ainv( 1,2) + Ainv( 2,1)
!    Ainv( 2,1) = Ainv( 1,2) - Ainv( 2,1)
!    Ainv( 1,2) = Ainv( 1,2) - Ainv( 2,1)
!
!    Ainv( 1,3) = Ainv( 1,3) + Ainv( 3,1)
!    Ainv( 3,1) = Ainv( 1,3) - Ainv( 3,1)
!    Ainv( 1,3) = Ainv( 1,3) - Ainv( 3,1)
!
!    Ainv( 2,3) = Ainv( 2,3) + Ainv( 3,2)
!    Ainv( 3,2) = Ainv( 2,3) - Ainv( 3,2)
!    Ainv( 2,3) = Ainv( 2,3) - Ainv( 3,2)
!
!    ! Divide by det(A)
!    Ainv = Ainv / detA
!
!  END SUBROUTINE calc_matrix_inverse_3_by_3
!
!  SUBROUTINE calc_matrix_inverse_general( A, Ainv)
!    ! Calculate the inverse Ainv of an n-by-n matrix A using LAPACK
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: A
!    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: Ainv
!
!    ! Local variables:
!    REAL(dp), DIMENSION( SIZE( A,1))                   :: work     ! work array for LAPACK
!    INTEGER,  DIMENSION( SIZE( A,1))                   :: ipiv     ! pivot indices
!    INTEGER                                            :: n,info
!
!    n = size( A,1)
!
!    ! Store A in Ainv to prevent it from being overwritten by LAPACK
!    Ainv = A
!
!    ! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
!    CALL DGETRF( n, n, Ainv, n, ipiv, info)
!
!    ! Safety
!    IF (info /= 0) THEN
!      ! IF (n == 5) THEN
!        ! PRINT *, A(1,1), A(1,2), A(1,3), A(1,4), A(1,5)
!        ! PRINT *, A(2,1), A(2,2), A(2,3), A(2,4), A(2,5)
!        ! PRINT *, A(3,1), A(3,2), A(3,3), A(3,4), A(3,5)
!        ! PRINT *, A(4,1), A(4,2), A(4,3), A(4,4), A(4,5)
!        ! PRINT *, A(5,1), A(5,2), A(5,3), A(5,4), A(5,5)
!      ! END IF
!      WRITE(0,*) 'calc_matrix_inverse_general - DGETRF error: matrix is numerically singular!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
!    CALL DGETRI( n, Ainv, n, ipiv, work, n, info)
!
!    ! Safety
!    IF (info /= 0) THEN
!      WRITE(0,*) 'calc_matrix_inverse_general - DGETRI error: matrix inversion failed!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!  END SUBROUTINE calc_matrix_inverse_general
!
!! == Debugging
!
!  SUBROUTINE check_for_NaN_dp_1D(  d, d_name)
!    ! Check if NaN values occur in the 1-D dp data field d
!    ! NOTE: parallelised!
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d
!    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
!
!    ! Local variables:
!    INTEGER                                                :: nx,i,i1,i2
!    CHARACTER(LEN=256)                                     :: d_name_loc
!
!    ! Only do this when so specified in the config
!    IF (.NOT. C%do_check_for_NaN) RETURN
!
!    ! Field size
!    nx = SIZE(d,1)
!
!    ! Parallelisation range
!    CALL partition_list( nx, par%i, par%n, i1, i2)
!
!    ! Variable name and routine name
!    IF (PRESENT( d_name)) THEN
!      d_name_loc = TRIM(d_name)
!    ELSE
!      d_name_loc = '?'
!    END IF
!
!    ! Inspect data field
!    DO i = i1, i2
!
!      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
!      ! you use the property that a NaN is never equal to anything, including itself...
!
!      IF     (d( i) /= d( i)) THEN
!        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
!      ELSEIF (d( i) > HUGE( d( i))) THEN
!        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
!      ELSEIF (d( i)  < -HUGE( d( i))) THEN
!        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]'), int_01 = i)
!      END IF
!
!    END DO
!    CALL sync
!
!  END SUBROUTINE check_for_NaN_dp_1D
!
!  SUBROUTINE check_for_NaN_dp_2D(  d, d_name)
!    ! Check if NaN values occur in the 2-D dp data field d
!    ! NOTE: parallelised!
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
!    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
!
!    ! Local variables:
!    INTEGER                                                :: nx,ny,i,j,i1,i2
!    CHARACTER(LEN=256)                                     :: d_name_loc
!
!    ! Only do this when so specified in the config
!    IF (.NOT. C%do_check_for_NaN) RETURN
!
!    ! Field size
!    nx = SIZE(d,2)
!    ny = SIZE(d,1)
!
!    ! Parallelisation range
!    CALL partition_list( nx, par%i, par%n, i1, i2)
!
!    ! Variable name and routine name
!    IF (PRESENT( d_name)) THEN
!      d_name_loc = TRIM(d_name)
!    ELSE
!      d_name_loc = '?'
!    END IF
!
!    ! Inspect data field
!    DO i = i1, i2
!    DO j = 1, ny
!
!      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
!      ! you use the property that a NaN is never equal to anything, including itself...
!
!      IF     (d( j,i) /= d( j,i)) THEN
!        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
!      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
!        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
!      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
!        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
!      END IF
!
!    END DO
!    END DO
!    CALL sync
!
!  END SUBROUTINE check_for_NaN_dp_2D
!
!  SUBROUTINE check_for_NaN_dp_3D(  d, d_name)
!    ! Check if NaN values occur in the 3-D dp data field d
!    ! NOTE: parallelised!
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: d
!    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
!
!    ! Local variables:
!    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
!    CHARACTER(LEN=256)                                     :: d_name_loc
!
!    ! Only do this when so specified in the config
!    IF (.NOT. C%do_check_for_NaN) RETURN
!
!    ! Field size
!    nx = SIZE(d,3)
!    ny = SIZE(d,2)
!    nz = SIZE(d,1)
!
!    ! Parallelisation range
!    CALL partition_list( nx, par%i, par%n, i1, i2)
!
!    ! Variable name and routine name
!    IF (PRESENT( d_name)) THEN
!      d_name_loc = TRIM(d_name)
!    ELSE
!      d_name_loc = '?'
!    END IF
!
!    ! Inspect data field
!    DO i = i1, i2
!    DO j = 1, ny
!    DO k = 1, nz
!
!      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
!      ! you use the property that a NaN is never equal to anything, including itself...
!
!      IF     (d( k,j,i) /= d( k,j,i)) THEN
!        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
!      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
!        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
!      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
!        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
!      END IF
!
!    END DO
!    END DO
!    END DO
!    CALL sync
!
!  END SUBROUTINE check_for_NaN_dp_3D
!
!  SUBROUTINE check_for_NaN_int_1D( d, d_name)
!    ! Check if NaN values occur in the 1-D int data field d
!    ! NOTE: parallelised!
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: d
!    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
!
!    ! Local variables:
!    INTEGER                                                :: nx,i,i1,i2
!    CHARACTER(LEN=256)                                     :: d_name_loc
!
!    ! Only do this when so specified in the config
!    IF (.NOT. C%do_check_for_NaN) RETURN
!
!    ! Field size
!    nx = SIZE(d,1)
!
!    ! Parallelisation range
!    CALL partition_list( nx, par%i, par%n, i1, i2)
!
!    ! Variable name and routine name
!    IF (PRESENT( d_name)) THEN
!      d_name_loc = TRIM(d_name)
!    ELSE
!      d_name_loc = '?'
!    END IF
!
!    ! Inspect data field
!    DO i = i1, i2
!
!      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
!      ! you use the property that a NaN is never equal to anything, including itself...
!
!      IF     (d( i) /= d( i)) THEN
!        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
!      ELSEIF (d( i) > HUGE( d( i))) THEN
!        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
!      ELSEIF (d( i) < -HUGE( d( i))) THEN
!        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]'), int_01 = i)
!      END IF
!
!    END DO
!    CALL sync
!
!  END SUBROUTINE check_for_NaN_int_1D
!
!  SUBROUTINE check_for_NaN_int_2D( d, d_name)
!    ! Check if NaN values occur in the 2-D int data field d
!    ! NOTE: parallelised!
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: d
!    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
!
!    ! Local variables:
!    INTEGER                                                :: nx,ny,i,j,i1,i2
!    CHARACTER(LEN=256)                                     :: d_name_loc
!
!    ! Only do this when so specified in the config
!    IF (.NOT. C%do_check_for_NaN) RETURN
!
!    ! Field size
!    nx = SIZE(d,2)
!    ny = SIZE(d,1)
!
!    ! Parallelisation range
!    CALL partition_list( nx, par%i, par%n, i1, i2)
!
!    ! Variable name and routine name
!    IF (PRESENT( d_name)) THEN
!      d_name_loc = TRIM(d_name)
!    ELSE
!      d_name_loc = '?'
!    END IF
!
!    ! Inspect data field
!    DO i = i1, i2
!    DO j = 1, ny
!
!      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
!      ! you use the property that a NaN is never equal to anything, including itself...
!
!      IF     (d( j,i) /= d( j,i)) THEN
!        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
!      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
!        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
!      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
!        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
!      END IF
!
!    END DO
!    END DO
!    CALL sync
!
!  END SUBROUTINE check_for_NaN_int_2D
!
!  SUBROUTINE check_for_NaN_int_3D( d, d_name)
!    ! Check if NaN values occur in the 3-D int data field d
!    ! NOTE: parallelised!
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    INTEGER,  DIMENSION(:,:,:),              INTENT(IN)    :: d
!    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
!
!    ! Local variables:
!    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
!    CHARACTER(LEN=256)                                     :: d_name_loc
!
!    ! Only do this when so specified in the config
!    IF (.NOT. C%do_check_for_NaN) RETURN
!
!    ! Field size
!    nx = SIZE(d,3)
!    ny = SIZE(d,2)
!    nz = SIZE(d,1)
!
!    ! Parallelisation range
!    CALL partition_list( nx, par%i, par%n, i1, i2)
!
!    ! Variable name and routine name
!    IF (PRESENT( d_name)) THEN
!      d_name_loc = TRIM(d_name)
!    ELSE
!      d_name_loc = '?'
!    END IF
!
!    ! Inspect data field
!    DO i = i1, i2
!    DO j = 1, ny
!    DO k = 1, nz
!
!      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
!      ! you use the property that a NaN is never equal to anything, including itself...
!
!      IF     (d( k,j,i) /= d( k,j,i)) THEN
!        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
!      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
!        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
!      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
!        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
!      END IF
!
!    END DO
!    END DO
!    END DO
!    CALL sync
!
!  END SUBROUTINE check_for_NaN_int_3D
!
!  ! == Basic array operations
!  SUBROUTINE permute_2D_dp(  d, wd, map)
!    ! Permute a 2-D array
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
!    INTEGER,                             INTENT(INOUT) :: wd
!    INTEGER,  DIMENSION(2),              INTENT(IN)    :: map
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'permute_2D_dp'
!    INTEGER                                            :: i,j,n1,n2,i1,i2
!    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_temp
!    INTEGER                                            :: wd_temp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (map( 1) == 1 .AND. map( 2) == 2) THEN
!      ! Trivial
!      RETURN
!    ELSEIF (map( 1) == 2 .AND. map( 2) == 1) THEN
!      ! 2-D transpose, as expected
!    ELSE
!      CALL crash('invalid permutation!')
!    END IF
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Allocate temporary memory
!    CALL allocate_shared_dp_2D( n1, n2, d_temp, wd_temp)
!
!    ! Copy data to temporary memory
!    d_temp( i1:i2,:) = d( i1:i2,:)
!    CALL sync
!
!    ! Deallocate memory
!    CALL deallocate_shared( wd)
!
!    ! Reallocate transposed memory
!    CALL allocate_shared_dp_2D( n2, n1, d, wd)
!
!    ! Copy and transpose data from temporary memory
!    DO i = i1, i2
!    DO j = 1, n2
!      d( j,i) = d_temp( i,j)
!    END DO
!    END DO
!    CALL sync
!
!    ! Deallocate temporary memory
!    CALL deallocate_shared( wd_temp)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE permute_2D_dp
!
!  SUBROUTINE permute_2D_int( d, wd, map)
!    ! Permute a 2-D array
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
!    INTEGER,                             INTENT(INOUT) :: wd
!    INTEGER,  DIMENSION(2),              INTENT(IN)    :: map
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'permute_2D_int'
!    INTEGER                                            :: i,j,n1,n2,i1,i2
!    INTEGER,  DIMENSION(:,:  ), POINTER                ::  d_temp
!    INTEGER                                            :: wd_temp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (map( 1) == 1 .AND. map( 2) == 2) THEN
!      ! Trivial
!      RETURN
!    ELSEIF (map( 1) == 2 .AND. map( 2) == 1) THEN
!      ! 2-D transpose, as expected
!    ELSE
!      CALL crash('invalid permutation!')
!    END IF
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Allocate temporary memory
!    CALL allocate_shared_int_2D( n1, n2, d_temp, wd_temp)
!
!    ! Copy data to temporary memory
!    d_temp( i1:i2,:) = d( i1:i2,:)
!    CALL sync
!
!    ! Deallocate memory
!    CALL deallocate_shared( wd)
!
!    ! Reallocate transposed memory
!    CALL allocate_shared_int_2D( n2, n1, d, wd)
!
!    ! Copy and transpose data from temporary memory
!    DO i = i1, i2
!    DO j = 1, n2
!      d( j,i) = d_temp( i,j)
!    END DO
!    END DO
!    CALL sync
!
!    ! Deallocate temporary memory
!    CALL deallocate_shared( wd_temp)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE permute_2D_int
!
!  SUBROUTINE permute_3D_dp(  d, wd, map)
!    ! Permute a 3-D array
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: d
!    INTEGER,                             INTENT(INOUT) :: wd
!    INTEGER,  DIMENSION(3),              INTENT(IN)    :: map
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'permute_3D_dp'
!    INTEGER                                            :: i,j,k,n1,n2,n3,i1,i2
!    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_temp
!    INTEGER                                            :: wd_temp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    n3 = SIZE( d,3)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Allocate temporary memory
!    CALL allocate_shared_dp_3D( n1, n2, n3, d_temp, wd_temp)
!
!    ! Copy data to temporary memory
!    d_temp( i1:i2,:,:) = d( i1:i2,:,:)
!    CALL sync
!
!    ! Deallocate memory
!    CALL deallocate_shared( wd)
!
!    ! Different permutation options
!    IF (map( 1) == 1 .AND. map( 2) == 2 .AND. map( 3) == 3) THEN
!      ! [i,j,k] -> [i,j,k] (trivial...)
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_dp_3D( n1, n2, n3, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( i,j,k) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 1 .AND. map( 2) == 3 .AND. map( 3) == 2) THEN
!      ! [i,j,k] -> [i,k,j]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_dp_3D( n1, n3, n2, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( i,k,j) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 2 .AND. map( 2) == 1 .AND. map( 3) == 3) THEN
!      ! [i,j,k] -> [j,i,k]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_dp_3D( n2, n1, n3, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( j,i,k) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 2 .AND. map( 2) == 3 .AND. map( 3) == 1) THEN
!      ! [i,j,k] -> [j,k,i]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_dp_3D( n2, n3, n1, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( j,k,i) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 3 .AND. map( 2) == 1 .AND. map( 3) == 2) THEN
!      ! [i,j,k] -> [k,i,j]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_dp_3D( n3, n1, n2, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( k,i,j) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 3 .AND. map( 2) == 2 .AND. map( 3) == 1) THEN
!      ! [i,j,k] -> [k,j,i]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_dp_3D( n3, n2, n1, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( k,j,i) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSE
!      CALL crash('invalid permutation!')
!    END IF
!
!    ! Deallocate temporary memory
!    CALL deallocate_shared( wd_temp)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE permute_3D_dp
!
!  SUBROUTINE permute_3D_int( d, wd, map)
!    ! Permute a 3-D array
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: d
!    INTEGER,                             INTENT(INOUT) :: wd
!    INTEGER,  DIMENSION(3),              INTENT(IN)    :: map
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'permute_3D_int'
!    INTEGER                                            :: i,j,k,n1,n2,n3,i1,i2
!    INTEGER,  DIMENSION(:,:,:), POINTER                ::  d_temp
!    INTEGER                                            :: wd_temp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    n3 = SIZE( d,3)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Allocate temporary memory
!    CALL allocate_shared_int_3D( n1, n2, n3, d_temp, wd_temp)
!
!    ! Copy data to temporary memory
!    d_temp( i1:i2,:,:) = d( i1:i2,:,:)
!    CALL sync
!
!    ! Deallocate memory
!    CALL deallocate_shared( wd)
!
!    ! Different permutation options
!    IF (map( 1) == 1 .AND. map( 2) == 2 .AND. map( 3) == 3) THEN
!      ! [i,j,k] -> [i,j,k] (trivial...)
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_int_3D( n1, n2, n3, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( i,j,k) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 1 .AND. map( 2) == 3 .AND. map( 3) == 2) THEN
!      ! [i,j,k] -> [i,k,j]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_int_3D( n1, n3, n2, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( i,k,j) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 2 .AND. map( 2) == 1 .AND. map( 3) == 3) THEN
!      ! [i,j,k] -> [j,i,k]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_int_3D( n2, n1, n3, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( j,i,k) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 2 .AND. map( 2) == 3 .AND. map( 3) == 1) THEN
!      ! [i,j,k] -> [j,k,i]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_int_3D( n2, n3, n1, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( j,k,i) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 3 .AND. map( 2) == 1 .AND. map( 3) == 2) THEN
!      ! [i,j,k] -> [k,i,j]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_int_3D( n3, n1, n2, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( k,i,j) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSEIF (map( 1) == 3 .AND. map( 2) == 2 .AND. map( 3) == 1) THEN
!      ! [i,j,k] -> [k,j,i]
!
!      ! Reallocate permuted memory
!      CALL allocate_shared_int_3D( n3, n2, n1, d, wd)
!
!      ! Copy and permuted data from temporary memory
!      DO i = i1, i2
!      DO j = 1, n2
!      DO k = 1, n3
!        d( k,j,i) = d_temp( i,j,k)
!      END DO
!      END DO
!      END DO
!      CALL sync
!
!    ELSE
!      CALL crash('invalid permutation!')
!    END IF
!
!    ! Deallocate temporary memory
!    CALL deallocate_shared( wd_temp)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE permute_3D_int
!
!  SUBROUTINE flip_1D_dp( d)
!    ! Flip a 1-D array
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_1D_dp'
!    INTEGER                                            :: i,nx,iopp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    nx = SIZE( d,1)
!
!    ! Flip the data
!    CALL sync
!    IF (par%master) THEN
!      DO i = 1, nx
!        iopp = nx + 1 - i
!        IF (iopp <= i) EXIT         ! [a  ] [b  ]
!        d( i   ) = d( i) + d( iopp) ! [a+b] [b  ]
!        d( iopp) = d( i) - d( iopp) ! [a+b] [a  ]
!        d( i   ) = d( i) - d( iopp) ! [b  ] [a  ]
!      END DO
!    END IF ! IF (par%master) THEN
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE flip_1D_dp
!
!  SUBROUTINE flip_2D_x1_dp( d)
!    ! Flip a 2-D array along the first dimension
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_2D_x1_dp'
!    INTEGER                                            :: i,j,n1,n2,j1,j2,iopp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    CALL partition_list( n2, par%i, par%n, j1, j2)
!
!    ! Flip the data
!    DO j = j1,j2
!      DO i = 1, n1
!        iopp = n1 + 1 - i
!        IF (iopp <= i) EXIT               ! [a  ] [b  ]
!        d( i   ,j) = d( i,j) + d( iopp,j) ! [a+b] [b  ]
!        d( iopp,j) = d( i,j) - d( iopp,j) ! [a+b] [a  ]
!        d( i   ,j) = d( i,j) - d( iopp,j) ! [b  ] [a  ]
!      END DO
!    END DO
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE flip_2D_x1_dp
!
!  SUBROUTINE flip_2D_x2_dp( d)
!    ! Flip a 2-D array along the second dimension
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_2D_x2_dp'
!    INTEGER                                            :: i,j,n1,n2,i1,i2,jopp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Flip the data
!    DO i = i1,i2
!      DO j = 1, n2
!        jopp = n2 + 1 - j
!        IF (jopp <= j) EXIT               ! [a  ] [b  ]
!        d( i,j   ) = d( i,j) + d( i,jopp) ! [a+b] [b  ]
!        d( i,jopp) = d( i,j) - d( i,jopp) ! [a+b] [a  ]
!        d( i,j   ) = d( i,j) - d( i,jopp) ! [b  ] [a  ]
!      END DO
!    END DO
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE flip_2D_x2_dp
!
!  SUBROUTINE flip_3D_x1_dp( d)
!    ! Flip a 3-D array along the first dimension
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_3D_x1_dp'
!    INTEGER                                            :: i,j,n1,n2,j1,j2,iopp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    CALL partition_list( n2, par%i, par%n, j1, j2)
!
!    ! Flip the data
!    DO j = j1,j2
!      DO i = 1, n1
!        iopp = n1 + 1 - i
!        IF (iopp <= i) EXIT                     ! [a  ] [b  ]
!        d( i   ,j,:) = d( i,j,:) + d( iopp,j,:) ! [a+b] [b  ]
!        d( iopp,j,:) = d( i,j,:) - d( iopp,j,:) ! [a+b] [a  ]
!        d( i   ,j,:) = d( i,j,:) - d( iopp,j,:) ! [b  ] [a  ]
!      END DO
!    END DO
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE flip_3D_x1_dp
!
!  SUBROUTINE flip_3D_x2_dp( d)
!    ! Flip a 3-D array along the second dimension
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_3D_x2_dp'
!    INTEGER                                            :: i,j,n1,n2,i1,i2,jopp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Flip the data
!    DO i = i1,i2
!      DO j = 1, n2
!        jopp = n2 + 1 - j
!        IF (jopp <= j) EXIT                     ! [a  ] [b  ]
!        d( i,j   ,:) = d( i,j,:) + d( i,jopp,:) ! [a+b] [b  ]
!        d( i,jopp,:) = d( i,j,:) - d( i,jopp,:) ! [a+b] [a  ]
!        d( i,j   ,:) = d( i,j,:) - d( i,jopp,:) ! [b  ] [a  ]
!      END DO
!    END DO
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE flip_3D_x2_dp
!
!  SUBROUTINE flip_3D_x3_dp( d)
!    ! Flip a 3-D array along the third dimension
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'flip_3D_x3_dp'
!    INTEGER                                            :: i,j,k,n1,n2,n3,i1,i2,kopp
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    n1 = SIZE( d,1)
!    n2 = SIZE( d,2)
!    n3 = SIZE( d,2)
!    CALL partition_list( n1, par%i, par%n, i1, i2)
!
!    ! Flip the data
!    DO i = i1,i2
!      DO j = 1, n2
!        DO k = 1, n3
!          kopp = n3 + 1 - k
!          IF (kopp <= k) EXIT                     ! [a  ] [b  ]
!          d( i,j,k   ) = d( i,j,k) + d( i,j,kopp) ! [a+b] [b  ]
!          d( i,j,kopp) = d( i,j,k) - d( i,j,kopp) ! [a+b] [a  ]
!          d( i,j,k   ) = d( i,j,k) - d( i,j,kopp) ! [b  ] [a  ]
!        END DO
!      END DO
!    END DO
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE flip_3D_x3_dp
!
!  ! == Check if a point lies inside a polygon (used in basin definition)
!  FUNCTION is_in_polygon( Vpoly, p) RESULT( isso)
!    ! Use the ray-casting algorithm to check if the point p = [px,py]
!    ! lies inside the polygon spanned by poly = [x1,y1; x2,y2; ...; xn,yn]
!
!    IMPLICIT NONE
!
!    ! In- and output variables
!    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Vpoly
!    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
!    LOGICAL                                            :: isso
!
!    ! Local variables:
!    REAL(dp), DIMENSION(2)                             :: q,r,s,llis
!    REAL(dp)                                           :: xmin,xmax,ymin,ymax
!    INTEGER                                            :: n_intersects
!    INTEGER                                            :: vi,vj,n_vertices
!    LOGICAL                                            :: do_cross
!    REAL(dp), PARAMETER                                :: tol_dist = 1E-5_dp
!
!    isso = .FALSE.
!
!    xmin = MINVAL( Vpoly(:,1))
!    xmax = MAXVAL( Vpoly(:,1))
!    ymin = MINVAL( Vpoly(:,2))
!    ymax = MAXVAL( Vpoly(:,2))
!
!    ! Quick test
!    IF (p(1) < xmin .OR. p(1) > xmax .OR. &
!        p(2) < ymin .OR. p(2) > ymax) THEN
!      isso = .FALSE.
!      RETURN
!    END IF
!
!    ! Define the endpoint of the east-pointing ray
!    q = [xmax + (xmax - xmin) / 10._dp, p(2)]
!
!    ! Determine how often the ray intersects the polygon
!
!    n_vertices   = SIZE( Vpoly,1)
!    n_intersects = 0
!
!    DO vi = 1, n_vertices
!
!      ! Find vertices spanning a polygon line section
!      IF (vi < n_vertices) THEN
!        vj = vi + 1
!      ELSE
!        vj = 1
!      END IF
!
!      ! Define the line section
!      r = Vpoly( vi,:)
!      s = Vpoly( vj,:)
!
!      ! Determine if the ray intersects the line section
!      IF ((r(2) < p(2) .AND. s(2) < p(2)) .OR. (r(2) > p(2) .AND. s(2) > p(2))) THEN
!        do_cross = .FALSE.
!      ELSE
!        CALL segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
!      END IF
!
!      IF (do_cross) n_intersects = n_intersects + 1
!
!    END DO ! DO vi = 1, n_vertices
!
!    ! If the number of intersections is odd, p lies inside the polygon
!    IF (MOD( n_intersects,2) == 1) THEN
!      isso = .TRUE.
!    ELSE
!      isso = .FALSE.
!    END IF
!
!  END FUNCTION is_in_polygon
!
!  SUBROUTINE segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
!    ! Find out if the line segments [pq] and [rs] intersect. If so, return
!    ! the coordinates of the point of intersection
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p, q, r, s
!    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: llis
!    LOGICAL,                  INTENT(OUT)         :: do_cross
!    REAL(dp),                 INTENT(IN)          :: tol_dist
!
!    ! Local variables:
!    REAL(dp), DIMENSION(2,2)                      :: A
!    REAL(dp), DIMENSION(2)                        :: x, b
!    INTEGER,  DIMENSION(2)                        :: IPIV
!    INTEGER                                       :: info
!
!    ! If pq and rs are colinear, define them as not intersecting
!    IF ( ABS( cross2( [q(1)-p(1), q(2)-p(2)], [s(1)-r(1), s(2)-r(2)] )) < tol_dist) THEN
!      llis = [0._dp, 0._dp]
!      do_cross = .FALSE.
!      RETURN
!    END IF
!
!    A(1,:) = [(p(1)-q(1)), (r(1)-s(1))]
!    A(2,:) = [(p(2)-q(2)), (r(2)-s(2))]
!    b = [(r(1)-q(1)), (r(2)-q(2))]
!
!    ! The LAPACK solver will overwrite the right-hand side b with the solution x. Therefore we
!    ! first copy the rhs in the solution vector x:
!    x = b
!
!    ! Solve Ax = b using LAPACK routine that solves matrix equation Ax=b for x (in double precision)
!    CALL DGESV( 2, 1, A, 2, IPIV, x, 2, info)
!
!    llis = [q(1) + x(1) * (p(1)-q(1)), q(2) + x(1) * (p(2)-q(2))]
!
!    IF (x(1)>0._dp .AND. x(1)<1._dp .AND. x(2)>0._dp .AND. x(2)<1._dp) THEN
!      do_cross = .TRUE.
!    ELSE
!      do_cross = .FALSE.
!    END IF
!
!  END SUBROUTINE segment_intersection
!
!  FUNCTION cross2( a,b) RESULT(z)
!    ! Vector product z between 2-dimensional vectors a and b
!
!    IMPLICIT NONE
!
!    REAL(dp), DIMENSION(2),     INTENT(IN)        :: a, b
!    REAL(dp)                                      :: z
!
!    z = (a(1)*b(2)) - (a(2)*b(1))
!
!  END FUNCTION cross2

  ! == Extras

!  SUBROUTINE time_display( region, t_end, dt_ave, it)
!    ! Little time display for the screen
!
!    implicit none
!
!    ! Input/Ouput variables
!    type(type_model_region), intent(in)  :: region
!    real(dp),                intent(in)  :: t_end
!    real(dp),                intent(in)  :: dt_ave
!    integer,                 intent(in)  :: it
!
!    ! Local variables
!    character(len=9)                     :: r_time, r_step, r_adv, r_ave
!
!    if (region%time + region%dt < t_end) then
!      r_adv = "no"
!      write(r_time,"(F8.3)") min(region%time,t_end) / 1000._dp
!      write(r_step,"(F6.3)") max(region%dt,0.001_dp)
!      write(*,"(A)",advance=trim(r_adv)) repeat(c_backspace,999) // &
!              "   t = " // trim(r_time) // " kyr - dt = " // trim(r_step) // " yr"
!    else
!      r_adv = "yes"
!      write(r_time,"(F8.3)") min(region%time,t_end) / 1000._dp
!      write(r_step, "(F6.3)") dt_ave / real(it-1,dp)
!      write(*,"(A)",advance=trim(r_adv)) repeat(c_backspace,999) // &
!            "   t = " // trim(r_time) // " kyr - dt_ave = " // trim(r_step) // " yr"
!    end if
!    if (region%do_output) then
!      r_adv = "no"
!      write(*,"(A)",advance=trim(r_adv)) repeat(c_backspace,999)
!    end if
!
!  END SUBROUTINE time_display

END MODULE utilities
