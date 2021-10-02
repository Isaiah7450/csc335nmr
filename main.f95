! Created by: Isaiah Hoffman
! Created on: September 28, 2021
program main
use numeric_type_library
use root_finder_library
implicit none

interface
  ! Checks if the given parameters are valid. If not,
  ! prints an error message indicating such.
  ! tol : double : Tolerance. Should be positive.
  ! ft : FilterType : Filter type. Should be 0, 1, or 2.
  ! fsize : integer : Filter size. Should be odd; for SG, should be 5,
  !   11, or 17.
  ! passes : integer : Number of filter passes. Should be nonnegative.
  ! im : IntegrationMethod : Integration method. Should be 0, 1, 2, or 3.
  subroutine validate_parameters(tol, ft, fsize, passes, im)
    implicit none
    real(kind = 8), intent(in) :: tol
    integer, intent(in) :: ft, fsize, passes, im
  end subroutine validate_parameters
end interface


! We will get around to function-izing stuff later (hopefully).
! To start off, let's try reading program options from standard input.
character(len = 40) :: input_name
real(kind = 8) :: baseline_adjust, tolerance
integer :: filter_type, filter_size, filter_passes, integration_method
character(len = 40) :: output_name

read *, input_name, baseline_adjust, tolerance
read *, filter_type, filter_size, filter_passes
read *, integration_method, output_name
call validate_parameters(tolerance, filter_type, &
  filter_size, filter_passes, integration_method) 

end program main

subroutine validate_parameters(tol, ft, fsize, passes, im)
  use numeric_type_library
  implicit none
  real(kind = 8), intent(in) :: tol
  integer, intent(in) :: ft, fsize, passes, im
  logical :: err = .false.

  if (tol <= 0D0) then
    err = .true.
    print *, "Error: Tolerance should be positive."
  endif
  if (ft .ne. No_Filter .and. ft .ne. SG_Filter .and. &
    ft .ne. Boxcar_Filter) then
    err = .true.
    print *, "Error: Filter should be 0, 1, or 2."
  endif
  if (ft .ne. No_Filter) then
    if (fsize .eq. 0) then
      err = .true.
      print *, "Error: Filter size should be positive and odd."
    elseif (mod(fsize, 2) == 0) then
      err = .true.
      print *, "Error: Filter size should be odd."
    elseif (ft .eq. SG_Filter .and. fsize .ne. 5 .and. &
      fsize .ne. 11 .and. fsize .ne. 17) then
      err = .true.
      print *, "Error: Savitzsky-Golay filter size should be 5, 11, or 17."
    endif
    if (passes <= 0) then
      err = .true.
      print *, "Error: Number of filter passes should be positive."
    endif
  else if (passes .ne. 0) then
    err = .true.
    print *, "Error: Number of filter passes should be zero."
  endif
  if (im .ne. Newton_Cotes_Integration .and. im .ne. Romberg_Integration .and. &
    im .ne. Adaptive_Integration .and. im .ne. Quadrature_Integration) then
    err = .true.
    print *, "Error: Integration method should be 0, 1, 2, or 3."
  endif
  if (err) then
    print *, "Note: Aborting program..."
    call exit(1)
  endif
end subroutine validate_parameters
