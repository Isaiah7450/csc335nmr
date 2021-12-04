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
end module linear_algebra_library
