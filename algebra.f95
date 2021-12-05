! Created by: Isaiah Hoffman
! Created on: December 4, 2021
! Functions for algebraic manipulations.
module linear_algebra_library
  use numeric_type_library
  implicit none
  interface swap
    module procedure swap_dcomplex, swap_integer
  end interface
  private :: swap, row_op_scale_plus
contains
  ! Multiplies a matrix by a vector.
  ! A : 2-D double complex array : An n * n matrix stored in row-major order.
  ! x : double complex array : A column vector with n entries.
  ! b : double complex array : The result is written here. It should be a
  !   column vector with n entries and preallocated.
  ! n : integer : The dimensions of the matrix.
  subroutine matrix_vector_multiply(A, x, b, n)
    implicit none
    integer, intent(in) :: n
    complex(kind = 8), dimension(n, n), intent(in) :: A
    complex(kind = 8), dimension(n), intent(in) :: x
    complex(kind = 8), dimension(n), intent(out) :: b
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
    integer, intent(in) :: n
    complex(kind = 8), dimension(n, n + 1), intent(inout) :: A
    complex(kind = 8), dimension(n), intent(out) :: x
    logical, intent(out) :: err

    integer :: i, j, k, p
    real(kind = 8) :: max_magnitude
    ! rp is for "row pointer"; we will simulate row changes.
    integer, dimension(n) :: rp

    err = .false.
    ! Adapted from algorithm 6.2 in the textbook (pg. 378)
    ! (This was something cool I saw online; it's a convenient way
    ! to initialize an array.)
    rp = (/ (i, i=1, n) /)
    do i = 1, n - 1
      ! Partial pivoting
      max_magnitude = 0D0
      do j = i, n
        if (abs(A(rp(j), i)) > max_magnitude) then
          max_magnitude = abs(A(rp(j), i))
          p = j
        endif
      enddo
      if (rp(i) .ne. rp(p)) call swap(rp(i), rp(p))
      ! Elimination
      do j = i + 1, n
        call row_op_scale_plus(A, n, A(rp(j), i) / A(rp(i), i), &
          rp(i), rp(j))
      enddo
    enddo
    ! Check for unique solution.
    if (A(rp(n), n) .eq. 0) then
      err = .true.
      x = dcmplx(0D0, 0D0)
      return
    endif
    x = dcmplx(0D0, 0D0)
    ! Start back-substitution
  end subroutine solve_matrix_partial_pivoting

  ! Utility subroutines.
  subroutine swap_integer(a, b)
    integer, intent(inout) :: a, b
    integer :: temp
    temp = a
    a = b
    b = temp
  end subroutine swap_integer

  subroutine swap_dcomplex(a, b)
    complex(kind = 8), intent(inout) :: a, b
    complex(kind = 8) :: temp
    temp = a
    a = b
    b = temp
  end subroutine swap_dcomplex

  ! Performs the row operation A_j = lambda * A_i + A_j.
  ! A : double complex 2-D array : The augmented n * 1 matrix
  !   to perform the row operation to.
  ! n : integer : The number of rows in the matrix.
  ! lambda : double complex : The scaling factor.
  ! i : integer : The row to add.
  ! j : integer : The row to store the result.
  subroutine row_op_scale_plus(A, n, lambda, i, j)
    implicit none
    integer, intent(in) :: n
    complex(kind = 8), dimension(n, n + 1), intent(inout) :: A
    complex(kind = 8), intent(in) :: lambda
    integer, intent(in) :: i, j
    integer :: k
    do k = 1, n + 1
      A(j, k) = lambda * A(i, k) + A(j, k)
    enddo
  end subroutine row_op_scale_plus
end module linear_algebra_library
