! Created by: Isaiah Hoffman
! Created on: September 28, 2021
program main
use numeric_type_library
use root_finder_library
use interpolation_library
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

  ! This subroutine reads in the point data from an input file
  ! and stores the data.
  ! input_name : string : The name of the input file.
  ! points : PointList : The data points information will be stored
  !   in this parameter. The subroutine will allocate the needed memory;
  !   the caller is responsible for handling clean up.
  subroutine read_input(input_name, points)
    use numeric_type_library
    implicit none
    character(len = 40), intent(in) :: input_name
    type(PointList), intent(out) :: points
  end subroutine read_input

  ! This subroutine sorts the data points by x-value
  ! in ascending order.
  ! points : The list of node points to sort. The updated
  !   list is then written back to this parameter.
  subroutine sort_points(points)
    use numeric_type_library
    implicit none
    type(PointList), intent(inout) :: points
  end subroutine sort_points
end interface

character(len = 40) :: input_name
real(kind = 8) :: baseline_adjust, tolerance
integer :: filter_type, filter_size, filter_passes, integration_method
character(len = 40) :: output_name

integer :: io_status
integer, parameter :: output_unit = 18

type(PointList) :: points

! To start off, let's try reading program options from standard input.
read *, input_name, baseline_adjust, tolerance
read *, filter_type, filter_size, filter_passes
read *, integration_method, output_name
call validate_parameters(tolerance, filter_type, &
  filter_size, filter_passes, integration_method) 
call read_input(input_name, points)
!print *, points%x
! Open output file for writing. This will be function-ized later.
! (First, check if it exists already and delete it if so.)
open(unit = output_unit, status = "old", access = "sequential", &
  form = "unformatted", recl = 1, iostat = io_status, file = output_name)
if (io_status == 0) close(output_unit, status="delete")
! Then open it for real.
open(unit = output_unit, status = "new", file = output_name)
! Print program options.
! @TODO: Place in subroutine.
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
if (integration_method .eq. Newton_Cotes_Integration_Method) then
  write(output_unit, *) "Newton-Cotes Integration"
elseif (integration_method .eq. Romberg_Integration_Method) then
  write(output_unit, *) "Romberg Integration"
elseif (integration_method .eq. Adaptive_Integration_Method) then
  write(output_unit, *) "Adaptive Integration"
else
  write(output_unit, *) "Quadrature Integration"
endif
write(output_unit, *) ""
write(output_unit, *) "Plot File Data"
write(output_unit, *) "File", ":", input_name
! Clean up resources.
close(output_unit)
deallocate(points%y_prime)
deallocate(points%y)
deallocate(points%x)
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
  if (im .ne. Newton_Cotes_Integration_Method &
    .and. im .ne. Romberg_Integration_Method &
    .and. im .ne. Adaptive_Integration_Method &
    .and. im .ne. Quadrature_Integration_Method) then
    err = .true.
    print *, "Error: Integration method should be 0, 1, 2, or 3."
  endif
  if (err) then
    print *, "Note: Aborting program..."
    call exit(1)
  endif
end subroutine validate_parameters

subroutine read_input(input_name, points)
  use numeric_type_library
  implicit none
  character(len = 40), intent(in) :: input_name
  type(PointList), intent(out) :: points

  integer, parameter :: input_unit = 17
  integer :: n = 0
  integer :: error
  real(kind = 8) :: throwaway, throwaway2
  character(len = 100) :: error_msg

  ! Open input file for reading.
  open(unit = input_unit, status = "old", access = "sequential", &
    form = "formatted", recl=1, file = input_name)
  ! The first time, we are just getting the number of points.
  error = 0
  do while (.true.)
    read(input_unit, *, iostat=error, iomsg=error_msg) throwaway, throwaway2
    if (error .eq. 0) then
      n = n + 1
    elseif (error .lt. 0) then
      ! Reached end of file.
      exit
    else
      ! Some other I/O error so just quit the program.
      print *, "Error: An error occurred while reading ", input_name
      print *, error_msg
      call exit(1)
    endif
  enddo
  ! Do allocations
  ! Note the +1 since the last line triggers EOF and doesn't increment n.
  allocate(points%x(n))
  allocate(points%y(n))
  allocate(points%y_prime(n))
  points%length = n
  close(input_unit)
  ! Reopen file; there's probably a way to reset the file pointer
  ! but I'm not concerned with looking that up right now.
  open(unit = input_unit, status = "old", access = "sequential", &
    form = "formatted", recl = 1, file = input_name)
  error = 0
  n = 1
  do while (error .ne. -1)
    read (input_unit, *, iostat=error) points%x(n), points%y(n)
    ! Derivative data isn't currently available, but maybe it
    ! would in the future.
    if (n .eq. points%length) then
      exit
    elseif (error .eq. 0) then
      n = n + 1
    elseif (error .ne. 0) then
      ! The first case should cover EOF, so this only happens
      ! if some other I/O error occurs.
      print *, "Error: An error occurred while reading ", input_name
      print *, error_msg
      call exit(1)
      exit
    endif
  enddo
  call sort_points(points)
end subroutine read_input

subroutine sort_points(points)
  use numeric_type_library
  implicit none
  type(PointList), intent(inout) :: points
  integer :: i, j
  real(kind = 8) :: temp_x, temp_y, temp_y_prime
  ! We will implement a simple bubble sort for now.
  do i = 1, points%length
    do j = 1, points%length
      if (points%x(i) < points%x(j)) then
        temp_x = points%x(i)
        temp_y = points%y(i)
        temp_y_prime = points%y_prime(i)
        points%x(i) = points%x(j)
        points%y(i) = points%y(j)
        points%y_prime(i) = points%y_prime(j)
        points%x(j) = temp_x
        points%y(j) = temp_y
        points%y_prime(j) = temp_y_prime
      endif
    enddo
  enddo
end subroutine sort_points

