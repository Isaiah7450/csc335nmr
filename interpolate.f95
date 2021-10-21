! Created by: Isaiah Hoffman
! Created on: October 1, 2021
module interpolation_library
use numeric_type_library
implicit none
contains

! Uses Lagrange interpolation to estimate the value of a function
! at the provided x-value.
! x : double : The point to approximate the function value of using the provided
!   data set.
! points : PointList : The list of known x and y values to interpolate
!   between.
! tol : double : This parameter is ignored but kept for interface convenience.
! returns : double : The approximation for f(x) using Lagrange interpolation.
pure function lagrange_interpolation(x, points, tol) result(out)
  implicit none
  type(PointList), intent(in) :: points
  real(kind = 8), intent(in) :: x, tol
  real(kind = 8) :: out

  integer :: i, k, n
  real(kind = 8) :: outer, inner

  ! Wow, I spent >2 hours trying to figure why my code didn't work
  ! in certain cases only to realize I forgot this statement...
  out = 0D0
  ! Don't forget Fortran starts arrays at 1...
  n = points%length - 1
  do k = 1, n + 1
    outer = points%y(k); inner = 1D0
    do i = 1, n + 1
      if (i .eq. k) then
        cycle
      endif
      inner = (x - points%x(i)) / (points%x(k) - points%x(i)) * inner
    enddo
    outer = outer * inner
    out = out + outer
  enddo
end function lagrange_interpolation

! Uses Newton's Divided Differences method to estimate the value of a
! function at a given point.
! x : double : The point to approximate.
! points : PointList : The list of known x and y values to interpolate
!   between.
! returns : double : The approximated function value at the point x
!   using the divided differences algorithm.
pure function divided_differences(x, points) result(out)
  type(PointList), intent(in) :: points
  real(kind = 8), intent(in) :: x
  real(kind = 8) :: out

  integer :: i, j, n
  real(kind = 8) :: inner
  real(kind = 8), dimension(:), allocatable :: coefficients

  n = points%length
  allocate(coefficients(n))
  call divided_differences_coefficients(points, coefficients)
  out = 0
  do i = 1, n
    inner = coefficients(i)
    do j = 1, i - 1
      inner = inner * (x - points%x(j))
    enddo
    out = out + inner
  enddo
  deallocate(coefficients)
end function divided_differences

! Obtains the coefficients of the divided difference interpolating
! polynomial.
! points : The set of points to interpolate between.
! out : The coefficients are written to this parameter.
!   It should already have been allocated to be the same size
!   as the number of points.
pure subroutine divided_differences_coefficients(points, out)
  type(PointList), intent(in) :: points
  real(kind = 8), dimension(:), intent(out) :: out
  
  real(kind = 8), dimension(:,:), allocatable :: table

  integer :: i, j, n
  n = points%length
  allocate(table(n, n))
  ! Algorithm from textbook, page 124
  do i = 1, n
    table(i, 1) = points%y(i)
  enddo
  do i = 2, n
    do j = 2, i
      table(i, j) = (table(i, j - 1) - table(i - 1, j - 1))
      table(i, j) = table(i, j) / (points%x(i) - points%x(i - j + 1))
    enddo
  enddo
  ! Copy coefficients to output parameter.
  do i = 1, n
    out(i) = table(i, i)
  enddo
  ! Clean up.
  deallocate(table)
end subroutine divided_differences_coefficients

! Uses Hermite polynomial interpolation to estimate the value of a
! function at a given point.
! x : double : The point to approximate.
! points : PointList : The list of known x and y values to interpolate
!   between. The first derivatives at the given x points should also
!   be included.
! returns : double : The approximated function value at the point x
!   using Hermite interpolating polynomials.
function hermite_interpolation(x, points) result(out)
  type(PointList), intent(in) :: points
  real(kind = 8), intent(in) :: x
  real(kind = 8) :: out

  real(kind = 8), dimension(:), allocatable :: coefficients
  real(kind = 8) :: subtotal
  integer :: i, j, n

  n = points%length
  allocate(coefficients(2 * n))
  call hermite_coefficients(points, coefficients)
  out = 0
  do i = 1, 2 * n
    subtotal = coefficients(i)
    do j = 2, i
      subtotal = subtotal * (x - points%x((j - 2) / 2 + 1))
    enddo
    out = out + subtotal
  enddo
  deallocate(coefficients)
end function hermite_interpolation

! Obtains the coefficients of the divided difference interpolating
! polynomial.
! points : The set of points to interpolate between.
! out : The coefficients are written to this parameter.
!   It should already have been allocated to be the same size
!   as the number of points * 2 + 1.
pure subroutine hermite_coefficients(points, out)
  type(PointList), intent(in) :: points
  real(kind = 8), dimension(:), intent(out) :: out

  real(kind = 8), dimension(:,:), allocatable :: table
  real(kind = 8), dimension(:), allocatable :: z
  integer :: n, i, j

  ! Algorithm from textbook page 139.
  ! The extra +1s are due to Fortran being base-1 indexing
  ! versus base-0 indexing.
  n = points%length - 1
  allocate(table(2 * n + 2, 2 * n + 2))
  allocate(z(2 * n + 2))
  do i = 0, n
    z(2 * i + 1) = points%x(i + 1)
    z(2 * i + 2) = points%x(i + 1)
    table(2 * i + 1, 1) = points%y(i + 1)
    table(2 * i + 2, 1) = points%y(i + 1)
    table(2 * i + 2, 2) = points%y_prime(i + 1)
    if (i .ne. 0) then
      table(2 * i + 1, 2) = (table(2 * i + 1, 1) &
        - table(2 * i, 1)) &
        / (z(2 * i + 1) - z(2 * i))
    endif
  enddo
  do i = 2, 2 * n + 1
    do j = 2, i
      table(i + 1, j + 1) = (table(i + 1, j + 1 - 1) &
        - table(i - 1 + 1, j - 1 + 1)) &
        / (z(i + 1) - z(i - j + 1))
    enddo
  enddo
  ! Copy to output. Also, I am pretty sure the textbook
  ! should say Q_{0,0} to Q_{2*(n+1)+1,2*(n+1)+1} instead
  ! of Q_{0,0} to Q_{2*n+1,2*n+1}. The former worked
  ! properly whereas the latter did not.
  do i = 0, 2 * n + 1
    out(i + 1) = table(i + 1, i + 1)
  enddo
  ! Clean up.
  deallocate(z)
  deallocate(table)
end subroutine hermite_coefficients

! This subroutine obtains all of the natural cubic spline coefficients
! for the given dataset.
! points : The set of points to interpolate between. Note that it is
!   assumed the points are sorted by x-value ascending.
! out : The coefficients a, b, c, d as well as the value x_j for each
!   spline. It should be preallocated to the appropriate size (5
!   * length).
pure subroutine natural_cubic_spline_all_coefficients(points, out)
  type(PointList), intent(in) :: points
  real(kind = 8), dimension(:), intent(out) :: out

  integer :: i, j, n
  real(kind = 8), dimension(:), allocatable :: a, b, c, d, h, alpha, l, mu, z
  n = points%length - 1
  allocate(a(n + 1))
  allocate(b(n + 1))
  allocate(c(n + 1))
  allocate(d(n + 1))
  allocate(h(n))
  allocate(alpha(n + 1))
  allocate(l(n + 1))
  allocate(mu(n + 1))
  allocate(z(n + 1))
  ! Zero-initialize everything.
  a = 0; b = 0; c = 0; d = 0; h = 0; alpha = 0; l = 0; mu = 0; z = 0
  ! Algorithm adapted from algorithm 3.4 (pg. 147) in textbook.
  do i = 0, n
    a(i + 1) = points%y(i + 1)
  enddo
  do i = 0, n - 1
    h(i + 1) = points%x(i + 1 + 1) - points%x(i + 1)
  enddo
  do i = 1, n - 1
    alpha(i + 1) = (3D0 / h(i + 1)) * (a(i + 1 + 1) - a(i + 1)) &
      - (3D0 / h(i + 1 - 1)) * (a(i + 1) - a(i + 1 - 1))
  enddo
  ! Now we solve a tridiagonal linear system...
  l(1) = 1D0
  mu(1) = 0D0
  z(1) = 0D0
  do i = 1, n - 1
    l(i + 1) = 2D0 * (points%x(i + 1 + 1) - points%x(i + 1 - 1)) &
      - h(i + 1 - 1) * mu(i + 1 - 1)
    mu(i + 1) = h(i + 1) / l(i + 1)
    z(i + 1) = (alpha(i + 1) - h(i + 1 - 1) * z(i + 1 - 1)) / l(i + 1)
  enddo
  l(n + 1) = 1D0
  z(n + 1) = 0D0
  c(n + 1) = 0D0
  do j = n - 1, 0, -1
    c(j + 1) = z(j + 1) - mu(j + 1) * c(j + 1 + 1)
    b(j + 1) = (a(j + 1 + 1) - a(j + 1)) / h(j + 1) &
      - h(j + 1) * (c(j + 1 + 1) + 2D0 * c(j + 1)) / 3D0
    d(j + 1) = (c(j + 1 + 1) - c(j + 1)) / (3D0 * h(j + 1))
  enddo
  ! Copy coefficients to coefficient array.
  do j = 0, n
    out(j * 5 + 1) = a(j + 1)
    out(j * 5 + 2) = b(j + 1)
    out(j * 5 + 3) = c(j + 1)
    out(j * 5 + 4) = d(j + 1)
    out(j * 5 + 5) = points%x(j + 1)
  enddo
  ! Clean up.
  deallocate(z)
  deallocate(mu)
  deallocate(l)
  deallocate(alpha)
  deallocate(h)
  deallocate(d)
  deallocate(c)
  deallocate(b)
  deallocate(a)
end subroutine natural_cubic_spline_all_coefficients

! Obtains the coefficients for the natural cubic spline that includes
! the point to be approximated.
! x : The point to approximate using a cubic spline.
! points : The set of points to interpolate between. Note that it is
!   assumed the points are sorted by x-value ascending.
! out : The coefficients a, b, c, d as well as the value x_j are
!   written to this value. It should be already allocated to be of
!   size 5.
pure subroutine natural_cubic_spline_coefficients(x, points, out)
  type(PointList), intent(in) :: points
  real(kind = 8), intent(in) :: x
  ! Dimension is : for consistency across other interpolation functions
  real(kind = 8), dimension(:), intent(out) :: out

  real(kind = 8), dimension(:), allocatable :: coefficients
  integer :: j, n
  allocate(coefficients(5 * points%length))
  call natural_cubic_spline_all_coefficients(points, coefficients)
  n = points%length - 1
  ! Copy to output (after finding the right spline).
  ! It doesn't look like we can approximate values outside the points
  ! domain.
  do j = 2, n + 1
    if (x < points%x(j)) then
      out(1) = coefficients(5 * j - 9)
      out(2) = coefficients(5 * j - 8)
      out(3) = coefficients(5 * j - 7)
      out(4) = coefficients(5 * j - 6)
      out(5) = coefficients(5 * j - 5)
      exit
    endif
  enddo
  if (x > points%x(n + 1)) then
    out(1) = coefficients(5 * (n + 2) - 9)
    out(2) = coefficients(5 * (n + 2) - 8)
    out(3) = coefficients(5 * (n + 2) - 7)
    out(4) = coefficients(5 * (n + 2) - 6)
    out(5) = coefficients(5 * (n + 2) - 5)
  endif
  ! Clean up.
  deallocate(coefficients)
end subroutine natural_cubic_spline_coefficients

! Uses a natural cubic spline to estimate the value of a
! function at a given point.
! x : double : The point to approximate.
! points : PointList : The list of known x and y values to interpolate
!   between.
! tol : double : This parameter is ignored but kept for interface convenience.
! returns : double : The approximated function value at the point x
!   using the 
pure function natural_cubic_spline_interpolation(x, points, tol) result(out)
  type(PointList), intent(in) :: points
  real(kind = 8), intent(in) :: x, tol
  real(kind = 8) :: out

  real(kind = 8), dimension(:), allocatable :: coefficients

  allocate(coefficients(5))
  call natural_cubic_spline_coefficients(x, points, coefficients)
  out = coefficients(1)
  out = out + coefficients(2) * (x - coefficients(5))
  out = out + coefficients(3) * (x - coefficients(5)) ** 2
  out = out + coefficients(4) * (x - coefficients(5)) ** 3
end function natural_cubic_spline_interpolation

! This subroutine obtains all of the natural cubic spline coefficients
! for the given dataset.
! points : The set of points to interpolate between. Note that it is
!   assumed the points are sorted by x-value ascending. The values
!   y_prime(1) and y_prime(length) should be set as well.
! out : The coefficients a, b, c, d as well as the value x_j for each
!   spline. It should be preallocated to the appropriate size (5
!   * length).
subroutine clamped_cubic_spline_all_coefficients(points, out)
  type(PointList), intent(in) :: points
  real(kind = 8), dimension(:), intent(out) :: out

  integer :: i, j, n
  real(kind = 8) :: fpo, fpn
  real(kind = 8), dimension(:), allocatable :: a, b, c, d, h, alpha, l, mu, z
  n = points%length - 1
  allocate(a(n + 1))
  allocate(b(n + 1))
  allocate(c(n + 1))
  allocate(d(n + 1))
  allocate(h(n))
  allocate(alpha(n + 1))
  allocate(l(n + 1))
  allocate(mu(n + 1))
  allocate(z(n + 1))
  ! Zero-initialize everything.
  a = 0; b = 0; c = 0; d = 0; h = 0; alpha = 0; l = 0; mu = 0; z = 0
  ! Algorithm adapted from algorithm 3.5 (pg. 152) in textbook.
  fpo = points%y_prime(1)
  fpn = points%y_prime(points%length)
  do i = 0, n
    a(i + 1) = points%y(i + 1)
  enddo
  do i = 0, n - 1
    h(i + 1) = points%x(i + 1 + 1) - points%x(i + 1)
  enddo
  alpha(1) = 3D0 * (a(1 + 1) - a(0 + 1)) / h(0 + 1) - 3D0 * fpo
  alpha(n + 1) = 3D0 * fpn - 3D0 * (a(n + 1) - a(n + 1 - 1)) / h(n + 1 - 1)
  do i = 1, n - 1
    alpha(i + 1) = (3D0 / h(i + 1)) * (a(i + 1 + 1) - a(i + 1)) &
      - (3D0 / h(i + 1 - 1)) * (a(i + 1) - a(i + 1 - 1))
  enddo
  ! Now we solve a tridiagonal linear system...
  l(1) = 2D0 * h(0 + 1)
  mu(1) = 0.5D0
  z(1) = alpha(0 + 1) / l(0 + 1)
  do i = 1, n - 1
    l(i + 1) = 2D0 * (points%x(i + 1 + 1) - points%x(i + 1 - 1)) &
      - h(i + 1 - 1) * mu(i + 1 - 1)
    mu(i + 1) = h(i + 1) / l(i + 1)
    z(i + 1) = (alpha(i + 1) - h(i + 1 - 1) * z(i + 1 - 1)) / l(i + 1)
  enddo
  l(n + 1) = h(n + 1 - 1) * (2D0 - mu(n + 1 - 1))
  z(n + 1) = (alpha(n + 1) - h(n + 1 - 1) * z(n + 1 - 1)) / l(n + 1)
  c(n + 1) = z(n + 1)
  do j = n - 1, 0, -1
    c(j + 1) = z(j + 1) - mu(j + 1) * c(j + 1 + 1)
    b(j + 1) = (a(j + 1 + 1) - a(j + 1)) / h(j + 1) &
      - h(j + 1) * (c(j + 1 + 1) + 2D0 * c(j + 1)) / 3D0
    d(j + 1) = (c(j + 1 + 1) - c(j + 1)) / (3D0 * h(j + 1))
  enddo
  ! Copy coefficients to coefficient array.
  do j = 0, n
    out(j * 5 + 1) = a(j + 1)
    out(j * 5 + 2) = b(j + 1)
    out(j * 5 + 3) = c(j + 1)
    out(j * 5 + 4) = d(j + 1)
    out(j * 5 + 5) = points%x(j + 1)
  enddo
  ! Clean up.
  deallocate(z)
  deallocate(mu)
  deallocate(l)
  deallocate(alpha)
  deallocate(h)
  deallocate(d)
  deallocate(c)
  deallocate(b)
  deallocate(a)
end subroutine clamped_cubic_spline_all_coefficients

end module interpolation_library

