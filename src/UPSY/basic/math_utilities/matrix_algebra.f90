module matrix_algebra

  ! Directly calculate determinants and inverses of small matrices

  use precisions, only: dp

  implicit none

contains

  pure function calc_determinant_2_by_2( A) result( detA)
    ! Determinant of a 2-by-2 matrix

    real(dp), dimension(2,2), intent(in) :: A
    real(dp)                             :: detA

    detA = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

  end function calc_determinant_2_by_2

  pure function calc_determinant_3_by_3( A) result( detA)
    ! Determinant of a 3-by-3 matrix

    real(dp), dimension(3,3), intent(in) :: A
    real(dp)                             :: detA
    real(dp), dimension(3,3)             :: Am

    ! Calculate the minors of A
    Am( 1,1) = A( 2,2) * A( 3,3) - A( 2,3) * A( 3,2)
    Am( 1,2) = A( 2,1) * A( 3,3) - A( 2,3) * A( 3,1)
    Am( 1,3) = A( 2,1) * A( 3,2) - A( 2,2) * A( 3,1)
    ! Am( 2,1) = A( 1,2) * A( 3,3) - A( 1,3) * A( 3,2)  ! Don't need these
    ! Am( 2,2) = A( 1,1) * A( 3,3) - A( 1,3) * A( 3,1)
    ! Am( 2,3) = A( 1,1) * A( 3,2) - A( 1,2) * A( 3,1)
    ! Am( 3,1) = A( 1,2) * A( 2,3) - A( 1,3) * A( 2,2)
    ! Am( 3,2) = A( 1,1) * A( 2,3) - A( 1,3) * A( 2,1)
    ! Am( 3,3) = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Calculate the determinant of A
    detA = A( 1,1) * Am( 1,1) - A( 1,2) * Am( 1,2) + A( 1,3) * Am( 1,3)

  end function calc_determinant_3_by_3

  pure function calc_determinant_5_by_5( A) result( detA)
    ! Determinant of a 5-by-5 matrix
    !
    ! Source: https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt, accessed 2023-02-07

    use iso_fortran_env, only: real64

    ! Input variables:
    real(dp), dimension(5,5), intent(in) :: A
    real(dp)                             :: detA

    ! Local variables:
    real(real64) :: A11, A12, A13, A14, A15, A21, A22, A23, A24, &
        A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
        A51, A52, A53, A54, A55

    A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
    A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
    A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
    A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
    A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

    detA = real(                                                         &
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

  end function calc_determinant_5_by_5

  pure function calc_matrix_inverse_2_by_2( A) result( Ainv)
    ! Direct inversion of a 2-by-2 matrix

    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(2,2)             :: Ainv
    real(dp)                             :: detA

    ! Calculate the determinant of A
    detA = calc_determinant_2_by_2( A)

    ! Safety
    if (abs( detA) < tiny( detA)) then
      error stop 'calc_matrix_inverse_2_by_2 - matrix is singular to working precision!'
    end if

    ! Calculate the inverse of A
    Ainv( 1,1) =  A( 2,2) / detA
    Ainv( 1,2) = -A( 1,2) / detA
    Ainv( 2,1) = -A( 2,1) / detA
    Ainv( 2,2) =  A( 1,1) / detA

  end function calc_matrix_inverse_2_by_2

  pure function calc_matrix_inverse_3_by_3( A) result( Ainv)
    ! Direct inversion of a 3-by-3 matrix
    !
    ! See: https://metric.ma.ic.ac.uk/metric_public/matrices/inverses/inverses2.html

    real(dp), dimension(3,3), intent(in) :: A
    real(dp), dimension(3,3)             :: Ainv
    real(dp)                             :: detA

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
    if (abs( detA) < tiny( detA)) then
      error stop 'calc_matrix_inverse_3_by_3 - matrix is singular to working precision!'
    end if

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

  end function calc_matrix_inverse_3_by_3

  pure function calc_matrix_inverse_5_by_5( A) result( Ainv)
    ! Direct inversion of a 5-by-5 matrix
    !
    ! Source: https://caps.gsfc.nasa.gov/simpson/software/m55inv_f90.txt, accessed 2023-02-07

    use iso_fortran_env, only: real64

    real(dp), dimension(5,5), intent(in) :: A
    real(dp), dimension(5,5)             :: Ainv
    real(dp)                             :: detA

    ! Local variables:
    real(real64) :: A11, A12, A13, A14, A15, A21, A22, A23, A24, &
        A25, A31, A32, A33, A34, A35, A41, A42, A43, A44, A45,   &
        A51, A52, A53, A54, A55
    real(real64), dimension(5,5) :: COFACTOR

    A11=A(1,1); A12=A(1,2); A13=A(1,3); A14=A(1,4); A15=A(1,5)
    A21=A(2,1); A22=A(2,2); A23=A(2,3); A24=A(2,4); A25=A(2,5)
    A31=A(3,1); A32=A(3,2); A33=A(3,3); A34=A(3,4); A35=A(3,5)
    A41=A(4,1); A42=A(4,2); A43=A(4,3); A44=A(4,4); A45=A(4,5)
    A51=A(5,1); A52=A(5,2); A53=A(5,3); A54=A(5,4); A55=A(5,5)

    detA = real(                                                          &
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
    if (abs( detA) < tiny( detA)) then
      error stop 'calc_matrix_inverse_3_by_3 - matrix is singular to working precision!'
    end if

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

    AINV = real( transpose( COFACTOR),dp) / detA

  end function calc_matrix_inverse_5_by_5

  pure function solve_Axb_2_by_2( A,b) result( x)
    ! Direct solution of the 2-by-2 matrix equation Ax=b

    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(2  ), intent(in) :: b
    real(dp), dimension(2  )             :: x
    real(dp), dimension(2,2)             :: Ainv

    ! Calculate the inverse of A
    Ainv = calc_matrix_inverse_2_by_2( A)

    ! Calculate x
    x( 1) = Ainv( 1,1) * b( 1) + Ainv( 1,2) * b( 2)
    x( 2) = Ainv( 2,1) * b( 1) + Ainv( 2,2) * b( 2)

  end function solve_Axb_2_by_2

end module matrix_algebra
