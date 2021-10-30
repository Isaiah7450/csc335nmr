! Created by: Isaiah Hoffman
! Created on: October 1, 2021
! This module stores various data types. This allows the types to be
! used in multiple modules while reducing dependencies and helping
! prevent cyclic dependencies.
module numeric_type_library
implicit none
  ! Constants
  ! Mathematical constants.
  real(kind = 8), parameter :: pi = 3.1415926535897932384626433832795D0

  ! Constants for filter type.
  integer, parameter :: No_Filter = 0
  integer, parameter :: Boxcar_Filter = 1
  integer, parameter :: SG_Filter = 2
  ! Constants for integration method.
  integer, parameter :: Newton_Cotes_Integration_Method = 0
  integer, parameter :: Romberg_Integration_Method = 1
  integer, parameter :: Adaptive_Integration_Method = 2
  integer, parameter :: Quadrature_Integration_Method = 3

  ! This type is used to store a list of points for interpolation.
  type PointList
    ! The x-coordinates of the points
    real(kind = 8), dimension(:), allocatable :: x
    ! The corresponding function values. A 1-1 correspondence between
    ! x and y should exist (i.e.: they should be the same size.)
    real(kind = 8), dimension(:), allocatable :: y
    ! The derivative information (if it is available) for each point.
    real(kind = 8), dimension(:), allocatable :: y_prime
    ! The size of the list. This should match the dimensions of x and y.
    integer :: length
  end type PointList
contains
  ! No functions or subroutines.
end module numeric_type_library
