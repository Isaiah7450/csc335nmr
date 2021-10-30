! Created by: Isaiah Hoffman
! Created on: October 17, 2021
module numerical_calculus_library
use numeric_type_library
implicit none
  ! Using generic interfaces here to make specifying the number of intervals
  ! optional. Being able to manually specify is good if high performance is
  ! desirable or if a homework problem specifically asks for a certain number to
  ! be used.
  interface composite_trapezoid_rule
    module procedure composite_trapezoid_rule_man, composite_trapezoid_rule_auto
  end interface composite_trapezoid_rule

  interface composite_simpson_rule
    module procedure composite_simpson_rule_man, composite_simpson_rule_auto
  end interface composite_simpson_rule

  interface romberg_integration_table
    module procedure romberg_integration_table_man
    module procedure romberg_integration_table_auto
  end interface

  interface romberg_integration
    module procedure romberg_integration_man, romberg_integration_auto
  end interface
contains
  ! Part I: Numerical Differentiation
  ! Uses the midpoint formula to find the derivative of a function
  ! at the given point.
  ! f : function : The function to differentiate. Interpolating
  !   functions are also valid.
  ! x0 : double : The point to evaluate the derivative at.
  ! points : PointList : The list of interpolation points (if relevant).
  ! tol : double : The tolerance to use when evaluating f (if needed).
  ! h : double : The step size.
  ! n : integer : The number of points to use. This should be 3 or 5.
  ! returns : double : The approximation for the derivative.
  pure function midpoint_differentiation(f, x0, points, tol, h, n) result(out)
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
    real(kind = 8), intent(in) :: tol, x0, h
    integer, intent(in) :: n
    real(kind = 8) :: out

    if (n .eq. 3) then
      out = f(x0 + h, points, tol) - f(x0 - h, points, tol)
      out = out / (2D0 * h)
    elseif (n .eq. 5) then
      out = f(x0 - 2D0 * h, points, tol) - 8D0 * f(x0 - h, points, tol)
      out = out + 8D0 * f(x0 + h, points, tol) - f(x0 + 2D0 * h, points, tol)
      out = out / (12D0 * h)
    else
      ! Throw an error perhaps? Can't print because pure.
      out = 0D0
    endif
  end function midpoint_differentiation

  ! Uses the endpoint formula to find the derivative of a function
  ! at the given point.
  ! f : function : The function to differentiate.
  ! x0 : double : The point to evaluate the derivative at.
  ! points : PointList : The list of points to interpolate between (if
  !   needed.)
  ! tol : double : The tolerance to use when evaluating f (if needed).
  ! h : double : The step size.
  ! n : integer : The number of points to use. This should be 3 or 5.
  ! start : logical : Set this to true if this is the left endpoint
  !   of the function.
  ! returns : double : The approximation for the derivative.
  pure function endpoint_differentiation(f, x0, points, tol, h, n, start) result(out)
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
    real(kind = 8), intent(in) :: tol, x0, h
    integer, intent(in) :: n
    logical, intent(in) :: start
    real(kind = 8) :: out

    real(kind = 8) :: h0
    if (start) then
      h0 = h
    else
      h0 = -h
    endif
    if (n .eq. 3) then
      out = -3D0 * f(x0, points, tol) + 4D0 * f(x0 + h0, points, tol)
      out = out - f(x0 + 2D0 * h0, points, tol)
      out = out / (2D0 * h0)
    elseif (n .eq. 5) then
      out = -25D0 * f(x0, points, tol) + 48D0 * f(x0 + h0, points, tol)
      out = out - 36D0 * f(x0 + 2D0 * h0, points, tol)
      out = out + 16D0 * f(x0 + 3D0 * h0, points, tol)
      out = out - 3D0 * f(x0 + 4D0 * h0, points, tol)
      out = out / (12D0 * h0)
    else
      ! Throw an error perhaps? Can't print because pure.
      out = 0D0
    endif
  end function endpoint_differentiation

  ! @NOTE: This might actually be better as a subroutine
  !   that returns the table of calculated values.
  ! Performs Richardson's extrapolation to find a derivative.
  ! f : function : The original function.
  ! x0 : double : The point to evaluate the derivative of.
  ! points : PointList : The set of points to interpolate between (if needed).
  ! tol : double : The desired accuracy. Currently, this parameter is only
  !   used by df and does not affect other operations.
  ! h : double : The step size to use.
  ! n : integer : The number of iterations to perform. In the future, this
  !   may become optional (i.e.: allowed to be zero) if tol is nonzero.
  ! returns : double : The derivative approximation.
  function richardson_extrapolation(f, x0, points, tol, h, n) result(out)
    interface
      pure function f(x, points, tolerance) result(out)
        import PointList
        double precision, intent(in) :: x
        type(PointList), intent(in) :: points
        double precision, intent(in) :: tolerance
        double precision :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: x0, tol, h
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8) :: out

    real(kind = 8), dimension(:,:), allocatable :: table
    integer :: i, j
    allocate(table(n, n))
    do j = 1, n
      ! Calculate the N_1 values needed.
      ! Using the forward difference formula.
      table(1, j) = f(x0 + h / (2D0 ** (j - 1D0)), points, tol) &
        - f(x0, points, tol)
      table(1, j) = table(1, j) * (2D0 ** (j - 1D0)) / h
      print *, 1, j, table(1, j)
    enddo
    do i = 2, n
      do j = 1, n - i + 1
        table(i, j) = table(i - 1, j + 1) &
          + (table(i - 1, j + 1) - table(i - 1, j)) &
          / (4D0 ** (i - 1) - 1)
        print *, i, j, table(i, j)
      enddo
    enddo
    out = table(n, 1)
    deallocate(table)
  end function richardson_extrapolation
  ! Part II: Numerical Integration
  ! Unless otherwise specified, it is assumed that $a<b$.
  ! Estimates the value of an integral using the composite trapezoid rule.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional list of points if needed by an
  !   interpolating function $f$.
  ! tol : double : The desired tolerance. This is passed to $f$ and is otherwise
  !   not used.
  ! b : double : The point to stop integrating at.
  ! n : integer : The number of subintervals to use.
  ! returns : double : The estimated value of the integral from a to b
  !   using the composite trapezoid rule.
  pure function composite_trapezoid_rule_man(f, a, points, tol, b, n) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8) :: out

    real(kind = 8) :: h
    integer :: i

    h = (b - a) / n
    out = f(a, points, tol) + f(b, points, tol)
    do i = 1, n - 1
      out = out + 2D0 * f(a + i * h, points, tol)
    enddo
    out = out * h / 2D0
  end function composite_trapezoid_rule_man
  ! Estimates the value of a definite integral using the composite trapezoid
  ! rule.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional list of points if needed by an
  !   interpolating function $f$.
  ! tol : double : The desired tolerance. This is passed to $f$ and is also used
  !   to determine how many intervals to divide the integral into.
  ! b : double : The point to stop integrating at.
  ! returns : double : The estimated value of the integral from a to b
  !   using the composite trapezoid rule.
  pure function composite_trapezoid_rule_auto(f, a, points, tol, b) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    real(kind = 8) :: out
    integer :: n
    real(kind = 8) :: h, previous, current
    ! The 1000 here is rather arbitarily.
    integer, parameter :: max_n = 1000
    ! Come up with an estimate for $n$ first. We don't have
    ! the second derivative, but that's okay; we'll address
    ! that later.
    do n = 1, max_n
      h = (b - a) / n
      if ((b - a) / 12D0 * h * h < tol) then
        exit
      endif
    enddo
    previous = 0D0
    current = composite_trapezoid_rule_man(f, a, points, tol, b, n)
    ! This is how we're addressing the fact we don't have a second
    ! derivative.
    do while ((previous - current) > tol)
      ! This honestly may be too much of a jump.
      n = n * 2
      if (n > max_n) exit
      previous = current
      current = composite_trapezoid_rule_man(f, a, points, tol, b, n)
    enddo
    out = current
  end function composite_trapezoid_rule_auto

  ! Estimate the value of an integral using the composite form of
  ! Simpson's Rule.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional list of points if needed by an
  !   interpolating function $f$.
  ! tol : double : The desired tolerance. This is passed to $f$ and is otherwise
  !   not used.
  ! b : double : The point to stop integrating at.
  ! n : integer : The number of subintervals to use.
  ! returns : double : The estimated value of the integral from a to b
  !   using the composite Simpson rule.
  pure function composite_simpson_rule_man(f, a, points, tol, b, n) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8) :: out

    real(kind = 8) :: h
    integer :: i
    h = (b - a) / n
    out = f(a, points, tol) + f(b, points, tol)
    do i = 1, n / 2 - 1
      out = out + 2D0 * f(a + (2D0 * i) * h, points, tol)
    enddo
    do i = 1, n / 2
      out = out + 4D0 * f(a + (2D0 * i - 1D0) * h, points, tol)
    enddo
    out = out * h / 3D0
  end function composite_simpson_rule_man

  ! Estimates the value of an integral using the composite form of Simpson's
  ! Rule.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional list of points if needed by an
  !   interpolating function $f$.
  ! tol : double : The desired tolerance. This is passed to $f$ and is also used
  !   to determine how many intervals to divide the integral into.
  ! b : double : The point to stop integrating at.
  ! returns : double : The estimated value of the integral from a to b
  !   using the composite Simpson's Rule.
  pure function composite_simpson_rule_auto(f, a, points, tol, b) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    real(kind = 8) :: out
    integer :: n
    real(kind = 8) :: h, previous, current
    ! The 1000 here is rather arbitarily.
    integer, parameter :: max_n = 1000
    ! Come up with an estimate for $n$ first. We don't have
    ! the fourth derivative, but that's okay; we'll address
    ! that later.
    do n = 2, max_n, 2
      h = (b - a) / n
      if ((b - a) / 180D0 * h * h * h * h < tol) then
        exit
      endif
    enddo
    previous = 0D0
    current = composite_simpson_rule_man(f, a, points, tol, b, n)
    do while ((previous - current) > tol)
      ! This honestly may be too much of a jump.
      n = n * 2
      if (n > max_n) exit
      previous = current
      current = composite_simpson_rule_man(f, a, points, tol, b, n)
    enddo
    out = current
  end function composite_simpson_rule_auto

  ! Performs Romberg integration and returns the intermediary results.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : The tolerance to use. This is only used when calling $f$.
  ! b : double : The point to integrate to.
  ! n : integer : The number of iterations of Romberg integration to perform.
  ! table : 2D array : The results are stored here. The table should be
  !   preallocated to be of size n * n, and the caller is responsible for
  !   deallocating.
  pure subroutine romberg_integration_table_man(f, a, points, tol, b, n, table)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8), dimension(:,:), intent(inout) :: table

    ! Algorithm adapted from Algorithm 4.2 in textbook (pg. 216)
    real(kind = 8) :: h
    integer :: i, j, k
    table = 0D0
    table(1,1) = composite_trapezoid_rule(f, a, points, tol, b, 1)
    h = (b - a)
    do i = 2, n
      ! Trapezoid rule
      do k = 1, 2 ** (i - 2)
        table(i, 1) = table(i, 1) + f(a + (k - 0.5D0) * h, points, tol)
      enddo
      table(i, 1) = (table(i, 1) * h + table(i - 1, 1)) / 2D0
      ! Extrapolation
      do j = 2, i
        table(i, j) = table(i, j - 1) &
          + (table(i, j - 1) - table(i - 1, j - 1)) &
          / (4D0 ** (j - 1) - 1D0)
      enddo
      h = h / 2D0
      ! No need to overwrite because we're storing the whole table.
    enddo
  end subroutine romberg_integration_table_man

  ! Performs Romberg integration and returns the intermediary results.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : The tolerance to use. This is used when calling $f$ and
  !   to determine when to stop computing rows.
  ! b : double : The point to integrate to.
  ! table : 2D array : The results are stored here. The procedure will
  !   handle its allocation; however, the caller is responsible for
  !   deallocating it.
  pure subroutine romberg_integration_table_auto(f, a, points, tol, b, table)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    real(kind = 8), dimension(:,:), allocatable, intent(inout) :: table

    ! Algorithm adapted from Algorithm 4.2 in textbook (pg. 216)
    integer :: n
    real(kind = 8) :: h
    integer :: i, j, k, r, s
    ! For reallocations.
    real(kind = 8), dimension(:,:), allocatable :: temp_table
    ! These control the initial table allocation and reallocation
    ! growth rate.
    integer, parameter :: Starting_Size = 20
    real(kind = 8), parameter :: Growth_Factor = 1.5D0
    n = Starting_Size
    allocate(table(n, n))
    table = 0D0
    table(1,1) = composite_trapezoid_rule(f, a, points, tol, b, 1)
    h = (b - a)
    i = 2
    ! We won't take the slightly fancier approach described in the textbook.
    do while (.true.)
      ! Trapezoid rule
      do k = 1, 2 ** (i - 2)
        table(i, 1) = table(i, 1) + f(a + (k - 0.5D0) * h, points, tol)
      enddo
      table(i, 1) = (table(i, 1) * h + table(i - 1, 1)) / 2D0
      ! Extrapolation
      do j = 2, i
        table(i, j) = table(i, j - 1) &
          + (table(i, j - 1) - table(i - 1, j - 1)) &
          / (4D0 ** (j - 1) - 1D0)
      enddo
      h = h / 2D0
      if (abs(table(i - 1, i - 1) - table(i, i)) < tol) then
        exit
      endif
      i = i + 1
      if (i > n) then
        ! Time to resize. This admittedly should be a subroutine...
        n = floor(n * Growth_Factor)
        allocate(temp_table(n, n))
        do r = 1, i - 1
          do s = 1, i - 1
            temp_table(r, s) = table(r, s)
          enddo
        enddo
        deallocate(table)
        allocate(table(n, n))
        table = 0D0
        do r = 1, i - 1
          do s = 1, i - 1
            table(r, s) = temp_table(r, s)
          enddo
        enddo
        deallocate(temp_table)
      endif
    enddo
    n = i
    ! Do some shuffling.
    allocate(temp_table(n, n))
    do i = 1, n
      do j = 1, n
        temp_table(i, j) = table(i, j)
      enddo
    enddo
    deallocate(table)
    allocate(table(n, n))
    table = 0D0
    do i = 1, n
      do j = 1, n
        table(i, j) = temp_table(i, j)
      enddo
    enddo
    deallocate(temp_table)
  end subroutine romberg_integration_table_auto

  ! Performs Romberg integration and returns the integral estimate.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : The tolerance to use. This is only used when calling $f$.
  ! b : double : The point to integrate to.
  ! n : integer : The number of iterations of Romberg integration to perform.
  ! returns : double : The estimate for the integral using Romberg integration.
  pure function romberg_integration_man(f, a, points, tol, b, n) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8) :: out

    real(kind = 8), dimension(:,:), allocatable :: table
    allocate(table(n, n))
    call romberg_integration_table(f, a, points, tol, b, n, table)
    out = table(n, n)
    deallocate(table)
  end function romberg_integration_man

  ! Performs Romberg integration and returns the integral estimate.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : This determines when to stop computing rows.
  ! b : double : The point to integrate to.
  ! returns : double : The estimate for the integral using Romberg integration.
  pure function romberg_integration_auto(f, a, points, tol, b) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    real(kind = 8) :: out

    real(kind = 8), dimension(:,:), allocatable :: table
    call romberg_integration_table(f, a, points, tol, b, table)
    out = table(size(table, 1), size(table, 1))
    deallocate(table)
  end function romberg_integration_auto

  ! Performs Adaptive Quadrature integration using Simpson's Rule and returns the
  ! integral estimate.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : The desired accuracy. The program will stop subdividing
  !   based on 10 * tol.
  ! b : double : The point to integrate to.
  ! returns : double : The estimate for the integral using Adaptive Quadrature.
  function adaptive_quadrature(f, a, points, tol, b) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    real(kind = 8) :: out
    integer :: i
    integer, parameter :: Max_Level = 1000
    real(kind = 8), dimension(Max_Level) :: S, aa, toltol, h, fa, fb, fc
    real(kind = 8), dimension(4) :: v
    integer, dimension(Max_Level) :: L
    real(kind = 8) :: S1, S2, fd, fe
    S = 0D0; aa = 0D0; toltol = 0D0; h = 0D0; fa = 0D0; fb = 0D0; fc = 0D0
    v = 0D0; L = 0; S1 = 0D0; S2 = 0D0; fd = 0D0; fe = 0D0
    ! Algorithm adapted from Algorithm 4.3 on page 224 in the textbook.
    out = 0D0
    i = 1
    toltol(i) = 10D0 * tol
    aa(i) = a;
    h(i) = (b - a) / 2D0
    fa(i) = f(a, points, tol)
    fb(i) = f(b, points, tol)
    fc(i) = f(a + h(i), points, tol)
    S(i) = h(i) * (fa(i) + 4D0 * fc(i) + fb(i)) / 3D0
    L(i) = 1
    do while (i .gt. 0)
      fd = f(aa(i) + h(i) / 2D0, points, tol)
      fe = f(aa(i) + 3D0 * h(i) / 2D0, points, tol)
      S1 = (h(i) / 6D0) * (fa(i) + 4D0 * fd + fc(i))
      S2 = (h(i) / 6D0) * (fc(i) + 4D0 * fe + fb(i))
      ! I feel like the variable saves here are unnecessary
      ! and just obfuscate the code.
      ! But I don't really have the time to try to figure out
      ! the logic involved in getting rid of them.
      v(2) = fa(i)
      v(3) = fc(i)
      v(4) = fb(i)
      i = i - 1
      if (abs(S1 + S2 - S(i + 1)) < toltol(i + 1)) then
        ! We're within our tolerance or have done too many levels already.
        out = out + S1 + S2
      elseif (L(i + 1) .ge. Max_Level) then
        out = 0D0; exit
      else
        ! Consider right subinterval.
        i = i + 1
        aa(i) = aa(i) + h(i)
        fa(i) = v(3)
        fc(i) = fe
        fb(i) = v(4)
        h(i) = h(i) / 2D0
        toltol(i) = toltol(i) / 2D0
        S(i) = S2
        L(i) = L(i) + 1
        ! Consider left subinterval.
        i = i + 1
        aa(i) = aa(i - 1) - h(i - 1) * 2D0
        fa(i) = v(2)
        fc(i) = fd
        fb(i) = v(3)
        h(i) = h(i - 1)
        toltol(i) = toltol(i - 1)
        S(i) = S1
        L(i) = L(i - 1)
      endif
    enddo
  end function adaptive_quadrature

  ! Performs Gaussian-Legendra quadrature integration to estimate an integral
  ! with finite limits.
  ! f : function : The function to integrate.
  ! a : double : The point to start integrating from.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : This parameter is passed to $f$ and otherwise ignored.
  ! b : double : The point to integrate to.
  ! n : integer : The number of points to use.
  ! returns : double : The estimate for the integral using Gaussian Quadrature.
  pure function gauss_legendre_quadrature(f, a, points, tol, b, n) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: a, tol, b
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8) :: out

    real(kind = 8), dimension(n) :: w
    real(kind = 8), dimension(n) :: z
    integer :: i
    ! Set up weights and abscissas
    select case (n)
    case (2)
       z = (/ 1D0 / sqrt(3D0), -1D0 / sqrt(3D0) /)
       w = (/ 1D0, 1D0 /)
    case (3)
       z = (/ 0D0, sqrt(3D0 / 5D0), -sqrt(3D0 / 5D0) /)
       w = (/ 8D0 / 9D0, 5D0 / 9D0, 5D0 / 9D0 /)
    case (4)
       z = (/ -8.611363115940526D-1, -3.399810435848563D-1, &
         3.399810435848563D-1, 8.611363115940526D-1 /)
       w = (/ 3.478548451374477D-1, 6.521451548625462D-1, &
         6.521451548625462D-1, 3.478548451374477D-1 /)
    case (5)
      z = (/ 0D0, (1D0 / 3D0) * sqrt(5D0 - 2D0 * sqrt(10D0 / 7D0)), &
        -(1D0 / 3D0) * sqrt(5D0 - 2D0 * sqrt(10D0 / 7D0)), &
        (1D0 / 3D0) * sqrt(5D0 + 2D0 * sqrt(10D0 / 7D0)), &
        -(1D0 / 3D0) * sqrt(5D0 + 2D0 * sqrt(10D0 / 7D0)) /)
      w = (/ 128D0 / 225D0, (322D0 + 13D0 * sqrt(70D0)) / 900D0, &
        (322D0 + 13D0 * sqrt(70D0)) / 900D0, &
        (322D0 - 13D0 * sqrt(70D0)) / 900D0, &
        (322D0 - 13D0 * sqrt(70D0)) / 900D0 /)
    case default
      ! Not supported.
      out = 0D0
      return
    end select
    out = 0D0
    do i = 1, n
      out = out + w(i) * f(((b - a) * z(i) + (a + b)) / 2D0, points, tol) &
        * (b - a) / 2D0
    enddo
  end function gauss_legendre_quadrature

  ! Performs Gaussian-Laguerre quadrature integration to estimate an integral
  ! that goes from 0 to infinity.
  ! f : function : The function to integrate.
  ! points : PointList : An optional set of points to interpolate between.
  ! tol : double : This parameter is passed to $f$ and otherwise ignored.
  ! n : integer : The number of points to use.
  ! returns : double : The estimate for the integral using Gaussian Quadrature.
  function gauss_laguerre_quadrature(f, points, tol, n) result(out)
    interface
      pure function f(x, points, tol) result(out)
        import PointList
        real(kind = 8), intent(in) :: x, tol
        type(PointList), intent(in) :: points
        real(kind = 8) :: out
      end function f
    end interface
    real(kind = 8), intent(in) :: tol
    type(PointList), intent(in) :: points
    integer, intent(in) :: n
    real(kind = 8) :: out

    real(kind = 8), dimension(n) :: w
    real(kind = 8), dimension(n) :: z
    integer :: i

    select case (n)
    case (3)
      z = (/ 4.15774556783479D-1, 2.29428036027904D0, 6.28994508293748D0 /)
      w = (/ 7.11093009929173D-1, 2.78517733569241D-1, 1.03892565015861D-2 /)
    case (7)
      z = (/ 1.93043676560362D-1, 1.02666489533919D0, 1.02666489533919D0, &
        4.90035308452648D0, 8.18215344456286D0, 1.27341802917978D+1, &
        1.93957278622625D+1 /)
      w = (/ 4.09318951701273D-1, 4.21831277861720D-1, 1.47126348657505D-1, &
        2.06335144687169D-2, 1.07401014328075D-3, 1.58654643485641D-5, &
        3.17031547899558D-8 /)
    case (15)
      z = (/ 9.33078120172818D-2, 4.92691740301884D-1, 1.21559541207095D0, &
        2.26994952620374D0, 3.66762272175144D0, 5.42533662741355D0, &
        7.56591622661307D0, 1.01202285680191D+1, 1.31302824821757D+1, &
        1.66544077083300D+1, 2.07764788994488D+1, 2.56238942267288D+1, &
        3.14075191697539D+1, 3.85306833064860D+1, 4.80260855726858D+1 /)
      w = (/ 2.18234885940085D-1, 3.42210177922883D-1, 2.63027577941712D-1, &
        1.26425818105934D-1, 4.02068649210009D-2, 8.56387780361183D-3, &
        1.21243614721425D-3, 1.11674392344252D-4, 6.45992676202291D-6, &
        2.22631690709627D-7, 4.22743038497937D-9, 3.92189726704109D-11, &
        1.45651526407312D-13, 1.48302705111330D-16, 1.60059490621113D-20 /)
    case default
      ! Not supported.
      out = 0D0
      return
    end select
    out = 0D0
    do i = 1, n
      out = out + exp(z(i)) * w(i) * f(z(i), points, tol)
    enddo
  end function gauss_laguerre_quadrature
end module numerical_calculus_library
