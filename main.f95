! Created by: Isaiah Hoffman
! Created on: September 28, 2021
module utility_module
implicit none
  interface resize_list
    module procedure resize_list_auto, resize_list_man
  end interface resize_list
contains
  ! Resizes an array to have greater capacity. The new
  ! size is determined automatically.
  ! list : double array : The list to resize.
  subroutine resize_list_auto(list)
    real(kind = 8), dimension(:), allocatable, intent(inout) :: list
    real(kind = 8), parameter :: growth_factor = 2D0
    call resize_list(list, floor(size(list) * growth_factor))
  end subroutine resize_list_auto

  ! Resizes an array to have greater capacity. The new
  ! size is passed as a parameter.
  ! list : double array : The list to resize.
  ! new_size : integer : The new size for the array. It should
  !   be greater than the old size, but no checking is done.
  subroutine resize_list_man(list, new_size)
    real(kind = 8), dimension(:), allocatable, intent(inout) :: list
    integer, intent(in) :: new_size
    real(kind = 8), dimension(:), allocatable :: temp_list
    integer :: i, old_size
    old_size = size(list)
    allocate(temp_list(new_size))
    temp_list = 0D0
    do i = 1, min(old_size, new_size)
      temp_list(i) = list(i)
    enddo
    deallocate(list)
    allocate(list(new_size))
    list = 0D0
    do i = 1, min(old_size, new_size)
      list(i) = temp_list(i)
    enddo
    deallocate(temp_list)
  end subroutine resize_list_man
end module utility_module

module main_module
use numeric_type_library
use utility_module
use root_finder_library
use interpolation_library
use numerical_calculus_library
implicit none
contains
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

  ! This subroutine sorts the data points by x-value
  ! in ascending order.
  ! points : The list of node points to sort. The updated
  !   list is then written back to this parameter.
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

  ! Finds the TMS peak of the provided points list.
  ! points : PointList : The set of points to search.
  ! baseline : double : The location of the baseline.
  ! returns : double : The location of the TMS peak.
  function find_tms(points, baseline) result(out)
    type(PointList), intent(in) :: points
    real(kind = 8), intent(in) :: baseline
    real(kind = 8) :: out
    real(kind = 8) :: start, finish
    logical :: found_start
    integer :: i
    ! Search backwards for the starting and ending points.
    do i = points%length, 0, -1
      if (.not. found_start) then
        if (points%y(i) > baseline) then
          finish = points%x(i)
          found_start = .true.
        endif
      else
        if (points%y(i) < baseline) then
          start = points%x(i)
          exit
        endif
      endif
    enddo
    out = (start + finish) / 2D0
  end function find_tms

  ! This subroutine adjusts the x-values of the points
  ! based on the location of the TMS peak.
  ! points : PointList : The list of points to modify.
  ! tms : double : The x-coordinate of the TMS peak.
  subroutine adjust_tms(points, tms)
    type(PointList), intent(inout) :: points
    real(kind = 8), intent(in) :: tms
    integer :: i
    do i = 1, points%length
      points%x(i) = points%x(i) - tms
    enddo
  end subroutine adjust_tms

  ! This subroutine applies the provided filter to the input points.
  ! points : PointList : The list of node points.
  ! filter_type : integer : The type of filter to apply.
  ! filter_size : integer : The number of points to use in the filter.
  ! filter_passes : integer : The number of times to recursively apply
  !   the filter.
  subroutine apply_filter(points, filter_type, filter_size, filter_passes)
    use numeric_type_library
    implicit none
    type(PointList), intent(inout) :: points
    integer, intent(in) :: filter_type, filter_size, filter_passes
    if (filter_type .eq. Boxcar_Filter) then
      call apply_boxcar_filter(points, filter_size, filter_passes)
    elseif (filter_type .eq. SG_Filter) then
      call apply_sg_filter(points, filter_size, filter_passes)
    endif
  end subroutine apply_filter

  ! This subroutine applies a boxcar filter to the input points.
  ! points : PointList : The list of node points.
  ! filter_size : integer : The size of the filter to apply; should be odd.
  ! filter_passes : integer : The number of passes to apply.
  subroutine apply_boxcar_filter(points, filter_size, filter_passes)
    use numeric_type_library
    implicit none
    type(PointList), intent(inout) :: points
    integer, intent(in) :: filter_size, filter_passes
    integer :: i, j, n
    real(kind = 8), dimension(:), allocatable :: filtered_values
    allocate(filtered_values(points%length))
    do i = 1, filter_passes
      filtered_values = 0D0
      do j = 1, points%length
        do n = 1, filter_size
          ! Apply filter.
          filtered_values(j) = filtered_values(j) &
            + points%y(mod(j + n - 2 - filter_size / 2 &
            + points%length, points%length) + 1)
        enddo
        filtered_values(j) = filtered_values(j) / filter_size
      enddo
      ! Copy transformed points back.
      do n = 1, points%length
        points%y(n) = filtered_values(n)
      enddo
    enddo
    deallocate(filtered_values)
  end subroutine apply_boxcar_filter

  ! This subroutine applies a Savitzsky-Golay filter to the input points.
  ! points : PointList : The list of node points.
  ! filter_size : integer : The size of the filter; should be 5, 11, or 17.
  ! filter_passes : integer : The number of passes to apply.
  subroutine apply_sg_filter(points, filter_size, filter_passes)
    use numeric_type_library
    implicit none
    type(PointList), intent(inout) :: points
    integer, intent(in) :: filter_size, filter_passes
    integer :: i, j, iters
    real(kind = 8), dimension(:), allocatable :: filtered_values
    ! So np came from their sample code which I later learned
    ! is not the code we needed to implement, but I didn't want
    ! to have to massively rewrite this section, so I kept it as-is.
    real(kind = 8), dimension(filter_size) :: np
    allocate(filtered_values(points%length))
    do iters = 1, filter_passes
      do i = 1, points%length
        do j = 1, filter_size
          np(j) = points%y(mod(i + j - 2 - filter_size / 2 &
            + points%length, points%length) + 1)
        enddo
        ! 5, 11, 17
        if (filter_size .eq. 5) then
          ! (The weights come from table 1 in the paper.)
          filtered_values(i) = -3D0 * (np(1) + np(5)) + 12D0 * (np(2) + np(4)) &
            + 17D0 * np(3)
          filtered_values(i) = filtered_values(i) / 35D0
        elseif (filter_size .eq. 9) then
          filtered_values(i) = 59D0 * np(5) + 54D0 * (np(4) + np(6)) &
            + 39D0 * (np(3) + np(7)) + 14D0 * (np(2) + np(8)) &
            - 21D0 * (np(1) + np(9))
          ! (Just sum up the weights used.)
          filtered_values(i) = filtered_values(i) / 231D0
        elseif (filter_size .eq. 11) then
          filtered_values(i) = -36D0 * (np(1) + np(11)) &
            + 9D0 * (np(2) + np(10)) + 44D0 * (np(3) + np(9)) &
            + 69D0 * (np(4) + np(8)) + 84D0 * (np(5) + np(7)) &
            + 89D0 * (np(6))
          filtered_values(i) = filtered_values(i) / 429D0
        elseif (filter_size .eq. 17) then
          filtered_values(i) = -21D0 * (np(1) + np(17)) &
            - 6D0 * (np(2) + np(16)) + 7D0 * (np(3) + np(15)) &
            + 18D0 * (np(4) + np(14)) + 27D0 * (np(5) + np(13)) &
            + 34D0 * (np(6) + np(12)) + 39D0 * (np(7) + np(11)) &
            + 42D0 * (np(8) + np(10)) + 43D0 * (np(9))
          filtered_values(i) = filtered_values(i) / 323D0
        endif
      enddo
      do i = 1, points%length
        points%y(i) = filtered_values(i)
      enddo
    enddo
    deallocate(filtered_values)
  end subroutine apply_sg_filter

  ! This subroutine adjusts all the points so that
  ! only points above the baseline have a positive value.
  ! points : PointList : The list of points to modify.
  ! baseline : double : The location of the baseline.
  subroutine adjust_baseline(points, baseline)
    type(PointList), intent(inout) :: points
    real(kind = 8), intent(in) :: baseline
    integer :: i
    do i = 1, points%length
      points%y(i) = points%y(i) - baseline
    enddo
  end subroutine adjust_baseline

  ! Finds the start and end locations of the cubic spline's peaks.
  ! points : PointList : The list of points that are used to create
  !   the spline.
  ! peak_list : double array : The set of peak locations are written
  !   to this array in pairs. Odd indices correspond to starts and
  !   even indices correspond to finishes.
  ! tol : double : The tolerance for numerical algorithms.
  subroutine find_peaks(points, peak_list, tol)
    type(PointList), intent(in) :: points
    real(kind = 8), dimension(:), allocatable, intent(out) :: peak_list
    real(kind = 8), intent(in) :: tol
    real(kind = 8) :: start, finish, guess, root
    integer :: iterations, peak_index, i
    logical :: err
    integer, parameter :: initial_size = 26
    allocate(peak_list(initial_size))
    peak_index = 1
    i = 1
    do while (.true.)
      ! Search for start.
      do i = i, points%length
        if (points%y(i) >= 0D0) exit
      enddo
      if (i > points%length) exit
      ! Find start using bisection.
      if (i .ne. 1) then
        iterations = 0
        guess = (points%x(i) + points%x(i - 1)) / 2D0
        root = bisect(natural_cubic_spline_interpolation, &
          points, points%x(i - 1), points%x(i), guess, tol, err, iterations)
        start = root
      else
        start = points%x(1)
      endif
      ! Search for end.
      do i = i, points%length
        if (points%y(i) < 0D0) exit
      enddo
      ! Find actual finish using bisection.
      iterations = 0
      guess = (points%x(i) + points%x(i - 1)) / 2D0
      root = bisect(natural_cubic_spline_interpolation, &
        points, points%x(i - 1), points%x(i), guess, tol, err, iterations)
      if (err) then
        finish = points%x(points%length)
      else
        finish = root
      endif
      peak_list(peak_index) = start
      peak_list(peak_index + 1) = finish
      peak_index = peak_index + 2
      if (peak_index + 1 .ge. size(peak_list)) then
        call resize_list(peak_list)
      endif
    enddo
    call resize_list(peak_list, peak_index - 1)
  end subroutine find_peaks
end module main_module


program main
use numeric_type_library
use root_finder_library
use interpolation_library
use numerical_calculus_library
use main_module
implicit none

character(len = 40) :: input_name
real(kind = 8) :: baseline_adjust, tolerance
integer :: filter_type, filter_size, filter_passes, integration_method
character(len = 40) :: output_name
real(kind = 8) :: tms_location
real(kind = 8), dimension(:), allocatable :: peak_list

integer :: io_status
! @TODO: Delete when finished testing.
integer :: i, num_points
real(kind = 8) :: data_range

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
tms_location = find_tms(points, baseline_adjust)
call adjust_tms(points, tms_location)
call apply_filter(points, filter_type, filter_size, filter_passes)
call adjust_baseline(points, baseline_adjust)
call find_peaks(points, peak_list, tolerance)
! @TODO: Remove later: Testing code.
num_points = 10000
data_range = points%x(points%length) - points%x(1)
!print *, points%x
!do i = 1, num_points
!  print *, points%x(1) + i * data_range / num_points, &
!    natural_cubic_spline_interpolation(points%x(1) + i * data_range &
!      / num_points, points, tolerance)
!enddo
print *, peak_list
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
write(output_unit, *) "File", ": ", input_name
write(output_unit, *) "Plot shifted ", tms_location, " ppm for TMS calibration."
write(output_unit, *) ""
write(output_unit, *) ""
! If nothing else, this definitely should be put in some subroutine.
write(output_unit, *) "Peak | Begin | End | Location | Top | Area | Hydrogen"
do i = 1, size(peak_list) / 2
  write(output_unit, *) i, peak_list(i * 2 - 1), peak_list(i * 2), &
    (peak_list(i * 2) + peak_list(i * 2 - 1)) / 2D0, &
    ! Will write the correct values for these later.
    0D0, 0D0, 0
enddo
! Clean up resources.
close(output_unit)
deallocate(points%y_prime)
deallocate(points%y)
deallocate(points%x)
end program main
