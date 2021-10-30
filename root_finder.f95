! Created by: Isaiah Hoffman
! Created on: September 11, 2021
! This program helps one to find the roots of functions and equations.
module root_finder_library
use numeric_type_library
  implicit none

  private :: secant_slope
contains
  ! Uses the bisection algorithm to attempt to find a root.
  ! f : function : The function to find the root of. The provided interface
  !   definition allows one to use an interpolation function on a list of
  !   points. (For standard functions, the additional two arguments can be
  !   ignored.)
  ! points : PointList : The list of points to interpolate between. If
  !   an analytic function is used, this argument can be whatever.
  ! x_min : double : The minimum domain value.
  ! x_max : double : The maximum domain value.
  ! guess : double : The initial guess for the root.
  ! tolerance : double : The acceptable difference between the result
  !   and the true root.
  ! err : logical : This argument is set to true if an error occurs.
  ! iters : integer : This argument is set to the number of iterations
  !   performed. If it is initially set to a nonzero value, the algorithm
  !   will stop after that many iterations are performed regardless of the
  !   the accuracy.
  ! returns : double : The value of the root within the provided tolerance.
  !   If an error occurs, the error flag is set. If the error is not due
  !   to reaching the maximum iteration number, 0 is returned.
  function bisect(f, points, x_min, x_max, guess, tolerance, err, iters) result(out)
    implicit none
    ! It took me way too long to figure out how to specify the type
    ! of a function as an argument. I wonder if there is an easier
    ! way to specify this though as this is rather a lot.
    interface
      pure function f(x, points, tolerance) result(out)
        import PointList
        double precision, intent(in) :: x
        type(PointList), intent(in) :: points
        double precision, intent(in) :: tolerance
        double precision :: out
      end function f
    end interface
    type(PointList), intent(in) :: points
    double precision, intent(in) :: x_min, x_max, guess, tolerance
    logical, intent(out) :: err
    integer, intent(inout) :: iters
    double precision :: out

    double precision :: low, high, curr
    double precision :: f_low, f_high, f_curr
    integer :: max_iters
    ! Validate parameters. I might move some of these to the driver
    ! program.
    err = .false.
    !if (x_min > x_max) then
    !  print *, "Warning: Invalid domain specified."
    !  print *, "...Swapping values of x_min and x_max."
      ! We're using low here as a temporary variable.
    !  low = x_min
    !  x_min = x_max
    !  x_max = low
    !endif
    if (guess > x_max .or. guess < x_min) then
      print *, "Error: Bad guess."
      err = .true.
    endif
    if (tolerance <= 0D0) then
      print *, "Error: Tolerance must be nonzero."
      err = .true.
    endif
    if (iters <= 5 .and. iters .ne. 0) then
      print *, "Error: Max number of iterations too low."
      print *, "  For no limit, set this parameter to 0."
      err = .true.
    endif
    max_iters = iters
    iters = 0
    ! Check for a sign change.
    if ((f(x_min, points, tolerance) > 0D0 .and. &
      f(x_max, points, tolerance) > 0D0) .or. &
      (f(x_min, points, tolerance) < 0D0 .and. &
      f(x_max, points, tolerance) < 0D0)) then
      print *, "Error: No sign change on domain."
      err = .true.
    endif
    ! Leave if an error occurred.
    if (err) then
      out = 0
      return
    endif
    curr = guess
    low = x_min
    high = x_max
    do while ((high - low) * 2D0 > tolerance)
      f_low = f(low, points, tolerance)
      f_high = f(high, points, tolerance)
      f_curr = f(curr, points, tolerance)
      ! If we have by some miracle stumbled upon the root,
      ! then we should just exit.
      !if (f_curr .eq. 0D0) then
      !  exit
      !elseif (f_low .eq. 0D0) then
      !  exit
      !elseif (f_high .eq. 0D0) then
      !  exit
      !endif

      !print *, iters, low, curr, high, f_curr

      iters = iters + 1
      if (iters > max_iters .and. max_iters .ne. 0) then
        err = .true.
        exit
      endif
      ! Determine on which interval we have the sign change.
      if (f_low < 0D0 .and. f_curr > 0D0) then
        high = curr
      elseif (f_low > 0D0 .and. f_curr < 0D0) then
        high = curr
      elseif (f_high < 0D0 .and. f_curr > 0D0) then
        low = curr
      elseif (f_high > 0D0 .and. f_curr < 0D0) then
        low = curr
      else
        ! This case shouldn't be possible, but it is probably
        ! good to keep.
        print *, "Error: No sign change detected."
        err = .true.
        out = 0
        return
      endif
      ! Set new guess to be in the middle.
      curr = (high + low) / 2.0D0
    enddo
    out = curr
  end function bisect

  ! Uses the fixed point method to attempt to find a root.
  ! f : function : The real-valued function to find the root of.
  ! guess : double : The initial guess for the root.
  ! tolerance : double : The acceptable difference between successive
  !   iterations of the root.
  ! err : logical : This is set to true if an error occurs.
  ! iters : integer : The maximum number of iterations. If 0, the
  !   function will run indefinitely. Returns the number of iterations
  !   that were actually ran.
  ! returns : double : The value of the root to within the given tolerance.
  !   If an error occurs, the error flag is set. If the error is due to
  !   any reason besides reaching the maximum iteration number, 0 is returned.
  function fixed_point(f, guess, tolerance, err, iters) result(out)
    implicit none
    interface
      pure function f(x) result(out)
        double precision, intent(in) :: x
        double precision :: out
      end function f
    end interface
    double precision, intent(in) :: guess, tolerance
    logical, intent(out) :: err
    integer, intent(inout) :: iters
    double precision :: out
    
    double precision :: curr, prev, diff, prev_diff
    integer :: max_iters, threshold, max_threshold

    err = .false.
    prev = guess
    curr = f(guess) + guess
    diff = abs(curr - prev)
    max_iters = iters
    iters = 1
    threshold = 0
    max_threshold = 1
    print *, iters, prev, curr
    do while (diff > tolerance)
      prev_diff = diff
      prev = curr
      curr = f(curr) + curr
      diff = abs(curr - prev)
      iters = iters + 1
      print *, iters, prev, curr
      if (iters >= max_iters .and. iters .ne. 0) then
        err = .true.
        exit
      endif
      ! Probably not converging...
      if (diff > prev_diff) then
        threshold = threshold + 1
        if (threshold .ge. max_threshold) then
          print *, "Error: Function diverges."
          print *, "  Try using the inverse function or a different method."
          err = .true.
          out = 0
          return
        endif
      endif
    enddo
    out = curr
  end function fixed_point

  ! Uses Newton's method to find a root of the function.
  ! f : function : The real-valued function to find the root of.
  ! df : function : The derivative of f in analytic form.
  ! guess : double : An initial approximation for the root.
  ! tolerance : double : The minimum acceptable difference between
  !   iterations.
  ! err : logical : This is set to true if an error occurs.
  ! iters : integer : The maximum number of iterations to perform.
  !   If this is 0, then there is no limit. The function will
  !   write the actual number of iterations performed back to this
  !   parameter.
  ! returns : double : If successful, returns the calculated root.
  !   If unsuccessful due to reaching the maximum iterations, returns
  !   the current calculated value. Otherwise, returns 0.
  function newton_method(f, df, guess, tolerance, err, iters) result(out)
    implicit none
    interface
      double precision pure function f(x) result(out)
        double precision, intent(in) :: x
      end function
      double precision pure function df(x) result(out)
        double precision, intent(in) :: x
      end function
    end interface
    double precision, intent(in) :: guess, tolerance
    logical, intent(out) :: err
    integer, intent(inout) :: iters
    double precision :: out

    double precision :: x1, x0
    integer :: max_iters

    err = .false.
    max_iters = iters
    iters = 0
    x0 = guess
    x1 = x0 - f(x0) / df(x0)
    print *, iters, x0, f(x0), df(x0)
    iters = iters + 1
    do while (abs(x1 - x0) > tolerance)
      x0 = x1
      x1 = x0 - f(x0) / df(x0)
      print *, iters, x0, f(x0), df(x0)
      if (iters >= max_iters .and. max_iters .ne. 0) then
        err = .true.
        exit
      endif
      iters = iters + 1
    enddo
    out = x1
  end function newton_method

  function secant_method(f, guess, tolerance, err, iters) result(out)
    implicit none
    interface
      double precision pure function f(x) result(out)
        double precision, intent(in) :: x
      end function
    end interface
    double precision, intent(in) :: guess, tolerance
    logical, intent(out) :: err
    integer, intent(inout) :: iters
    double precision :: out

    double precision :: x0, x1
    integer :: max_iters

    err = .false.
    max_iters = iters
    iters = 0
    x0 = guess
    x1 = x0 - f(x0) / secant_slope(f, x0, 2D0 * tolerance)
    print *, iters, x0, f(x0), secant_slope(f, x0, 2D0 * tolerance)
    iters = iters + 1
    do while (abs(x1 - x0) > tolerance)
      x0 = x1
      x1 = x0 - f(x0) / secant_slope(f, x0, 2D0 * tolerance)
      print *, iters, x0, f(x0), secant_slope(f, x0, 2D0 * tolerance)
      if (iters >= max_iters .and. max_iters .ne. 0) then
        err = .true.
        exit
      endif
      iters = iters + 1
    enddo
    out = x1
  end function secant_method

  ! This function uses Steffenson's method to find a root.
  ! f : function : The function to find the root of. Please note that
  !   x is added to this function during the process of fixed point iteration.
  ! guess : double : An initial estimate for the root's location.
  ! tolerance : double : The desired accuracy.
  ! err : logical : This is set to true if an error occurs.
  ! iters : integer : The max number of iterations to perform (or 0 for no bound).
  !   The actual number of iterations performed is written to this parameter.
  ! returns : double : The root as obtained by Steffenson's method to the
  !   desired accuracy.
  function steffenson_method(f, guess, tolerance, err, iters) result(out)
    interface
      double precision pure function f(x) result(out)
        double precision, intent(in) :: x
      end function
    end interface
    double precision, intent(in) :: guess, tolerance
    logical, intent(out) :: err
    integer, intent(inout) :: iters
    double precision :: out

    double precision :: x0, x1, x2, x
    integer :: max_iters

    err = .false.
    max_iters = iters
    iters = 0

    ! Code copied and modified from example code given to us
    ! drlbs/converge/steff.f90
    x0 = guess
    x1 = f(x0) + x0
    do while ( iters .lt. max_iters) 
      x1 = f(x0) + x0
      x2 = f(x1) + x1
      x = x0 - (x1-x0)**2/(x2-2*x1+x0)
      if ( abs(x-x0) .le.  tolerance)  then 
        out = x
        exit
      endif
      iters=iters+1
      print *, iters, '  ', x0, '  ', x1 
      x0 = x
    enddo
  end function steffenson_method

  ! Uses Newton's method to find a root of the function.
  ! f : function : The real-valued function to find the root of.
  ! df : function : The derivative of f in analytic form.
  ! ddf : function : The second derivative of f in analytic form.
  ! guess : double : An initial approximation for the root.
  ! tolerance : double : The minimum acceptable difference between
  !   iterations.
  ! err : logical : This is set to true if an error occurs.
  ! iters : integer : The maximum number of iterations to perform.
  !   If this is 0, then there is no limit. The function will
  !   write the actual number of iterations performed back to this
  !   parameter.
  ! returns : double : If successful, returns the calculated root.
  !   If unsuccessful due to reaching the maximum iterations, returns
  !   the current calculated value. Otherwise, returns 0.
  function modified_newton_method(f, df, ddf, guess, tolerance, err, iters) result(out)
    implicit none
    interface
      double precision pure function f(x) result(out)
        double precision, intent(in) :: x
      end function
      double precision pure function df(x) result(out)
        double precision, intent(in) :: x
      end function
      double precision pure function ddf(x) result(out)
        double precision, intent(in) :: x
      end function
    end interface
    double precision, intent(in) :: guess, tolerance
    logical, intent(out) :: err
    integer, intent(inout) :: iters
    double precision :: out

    double precision :: x1, x0
    integer :: max_iters

    err = .false.
    max_iters = iters
    iters = 0
    x0 = guess
    x1 = x0 - (f(x0) * df(x0)) / (df(x0) * df(x0) - f(x0) * ddf(x0))
    print *, iters, x0, f(x0), df(x0), ddf(x0)
    iters = iters + 1
    do while (abs(x1 - x0) > tolerance)
      x0 = x1
      x1 = x0 - (f(x0) * df(x0)) / (df(x0) * df(x0) - f(x0) * ddf(x0))
      print *, iters, x0, f(x0), df(x0), ddf(x0)
      if (iters >= max_iters .and. max_iters .ne. 0) then
        err = .true.
        exit
      endif
      iters = iters + 1
    enddo
    out = x1
  end function modified_newton_method

  ! Converts from degrees to radians.
  pure function convert_to_radians(deg) result(out)
    double precision, intent(in) :: deg
    double precision :: out
    out = deg * pi / 180D0
  end function convert_to_radians

  ! Converts from radians to degrees.
  pure function convert_to_degrees(rad) result(out)
    double precision, intent(in) :: rad
    double precision :: out
    out = rad * 180D0 / pi
  end function

  ! Utility function that determines an approximation of
  ! the derivative of f at x using a secant line.
  ! f : function : The function to approximate the derivative of.
  ! x : double : The point on the function to approximate the value
  !   of the derivative at.
  ! interval : double : The width of the interval to use for the
  !   secant slope calculation. The interval is centered at x.
  ! returns : double : The approximation of the derivative of f at x
  !   using the secant slope.
  function secant_slope(f, x, interval) result(out)
    interface
      double precision pure function f(x) result(out)
        double precision, intent(in) :: x
      end function
    end interface
    double precision, intent(in) :: x, interval
    double precision :: out

    out = f(x + interval / 2D0) - f(x - interval / 2D0)
    out = out / interval
  end function secant_slope
end module root_finder_library

