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

! To start off, let's try reading program options from standard input.
character(len = 40) :: input_name
real(kind = 8) :: baseline_adjust, tolerance
integer :: filter_type, filter_size, filter_passes, integration_method
character(len = 40) :: output_name

integer :: io_status
integer, parameter :: input_unit = 17, output_unit = 18

read *, input_name, baseline_adjust, tolerance
read *, filter_type, filter_size, filter_passes
read *, integration_method, output_name
call validate_parameters(tolerance, filter_type, &
  filter_size, filter_passes, integration_method) 

! Open input file for reading.
open(unit = input_unit, status = "old", access = "direct", &
  form = "unformatted", recl = 1, file = input_name)
! Open output file for writing. This will be function-ized later.
! (First, check if it exists already and delete it if so.)
open(unit = output_unit, status = "old", access = "sequential", &
  form = "unformatted", recl = 1, iostat = io_status, file = output_name)
if (io_status == 0) close(output_unit, status="delete")
! Then open it for real.
open(unit = output_unit, status = "new", file = output_name)
! Print program options.
write(output_unit, *) "==> NMR Analysis <=="
write(output_unit, *) ""
write(output_unit, *) "Program Options"
write(output_unit, *) "Baseline Adjustment", ":", baseline_adjust
write(output_unit, *) "Tolerance", ":", tolerance
if (filter_type .eq. No_Filter) then
  write(output_unit, *) "No Filtering"
else
  if (filter_type .eq. Boxcar_Filter) then
    write(output_unit, *) "Boxcar Filtering"
    write(output_unit, *) "Boxcar Size (Cyclic)", ":", filter_size
    write(output_unit, *) "Boxcar Passes", ":", filter_passes
  else
    write(output_unit, *) "Savitzky-Golay Filtering"
    write(output_unit, *) "Filter Size", ":", filter_size
    write(output_unit, *) "Filter Passes", ":", filter_passes
  endif
endif
write(output_unit, *) ""
write(output_unit, *) "Integration Method"
if (integration_method .eq. Newton_Cotes_Integration) then
  write(output_unit, *) "Newton-Cotes Integration"
elseif (integration_method .eq. Romberg_Integration) then
  write(output_unit, *) "Romberg Integration"
elseif (integration_method .eq. Adaptive_Integration) then
  write(output_unit, *) "Adaptive Integration"
else
  write(output_unit, *) "Quadrature Integration"
endif
write(output_unit, *) ""
write(output_unit, *) "Plot File Data"
write(output_unit, *) "File", ":", input_name
! Close files.
close(output_unit)
close(input_unit)
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
