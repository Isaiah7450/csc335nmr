! Created by: Isaiah Hoffman
! Created on: October 1, 2021
! This module stores various data types. This allows the types to be
! used in multiple modules while reducing dependencies and helping
! prevent cyclic dependencies.
module numeric_type_library
implicit none
  ! This type is used to store a list of points for interpolation.
  type PointList
    ! The x-coordinates of the points
    real(kind = 8), dimension(:), allocatable :: x
    ! The corresponding function values. A 1-1 correspondence between
    ! x and y should exist (i.e.: they should be the same size.)
    real(kind = 8), dimension(:), allocatable :: y
    ! The size of the list. This should match the dimensions of x and y.
    integer :: length
  end type PointList
contains
  ! No functions or subroutines.
end module numeric_type_library