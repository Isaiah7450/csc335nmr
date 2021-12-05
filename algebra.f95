! Created by: Isaiah Hoffman
! Created on: December 4, 2021
! Functions for algebraic manipulations.
module linear_algebra_library
use numeric_type_library
implicit none
contains
  ! Multiplies a matrix by a vector.
  ! A : 2-D double complex array : An n * n matrix stored in row-major order.
  ! x : double complex array : A column vector with n entries.
  ! b : double complex array : The result is written here. It should be a
  !   column vector with n entries and preallocated.
  ! n : integer : The dimensions of the matrix.
  subroutine matrix_vector_multiply(A, x, b, n)
    implicit none
    complex(kind = 8), dimension(n, n), intent(in) :: A
    complex(kind = 8), dimension(n), intent(in) :: x
    complex(kind = 8), dimension(n), intent(out) :: b
    integer, intent(in) :: n
    integer :: i, j, k
    complex(kind = 8) :: temp_sum
    do i = 1, n
      temp_sum = dcmplx(0D0, 0D0)
      do j = 1, n
        temp_sum = temp_sum + A(i, j) * x(j)
      enddo
      b(i) = temp_sum
    enddo
  end subroutine matrix_vector_multiply

  ! Solves a linear system of equations using Gaussian elimination
  ! and partial pivoting along with back-substitution.
  ! A : 2-D double complex array : An n * (n + 1) augmented matrix
  !   [C b] where C is an n * n matrix specifying the coefficients
  !   of the system and b is the solution vector. The matrix should
  !   be in row-major order. The matrix may be modified as part of
  !   the solution process.
  ! x : double complex array : An n * 1 column vector where the solution
  !   will be written. It should be preallocated. If an error occurs,
  !   zeroes will be written to this vector.
  ! n : integer : The number of equations involved in the system.
  ! err : logical : A flag indicating whether or not an error occurred.
  subroutine solve_matrix_partial_pivoting(A, x, n, err)
    implicit none
    complex(kind = 8), dimension(n, n), intent(inout) :: A
    complex(kind = 8), dimension(n), intent(out) :: x
    integer, intent(in) :: n
    logical, intent(out) :: err

    integer :: i, j, k
    ! rp is for "row pointer"; we will simulate row changes.
    integer, dimension(n) :: rp

    err = .false.
    ! Adapted from algorithm 6.2 in the textbook (pg. 378)
    ! (This was something cool I saw online; it's a convenient way
    ! to initialize an array.)
    rp = (/ (i, i=1, n) /)
    print *, rp
  end subroutine solve_matrix_partial_pivoting
end module linear_algebra_library
