module endfTable_test
  use numPrecision
  use endfConstants
  use endfTable_class, only : endfTable, endf_bin_integral
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test ENDF integration over a bin
  !!
!@Test
  subroutine testEndfBinIntegral()
    real(defReal) :: x1, x0, y1, y0, val
    real(defReal), parameter :: TOL = 1.0E-5

    !! *** Test histogram
    x1 = 1.0_defReal
    x0 = 0.0_defReal
    y1 = 15.0_defReal
    y0 = 3.3_defReal
    ! @0.3
    val = 3.3_defReal * 0.3_defReal
#line 26 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 0.3_defReal, histogramInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! @1.0
    val = 3.3_defReal
#line 30 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 1.0_defReal, histogramInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    !! *** Test lin-lin Interpolation
    x1 = 2.0_defReal
    x0 = 1.0_defReal
    y1 = -1.0_defReal
    y0 = 1.0_defReal
    ! @1.5
    val = 0.25_defReal
#line 39 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 1.5_defReal, linLinInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! @2.0
    val = 0.0_defReal
#line 43 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 2.0_defReal, linLinInterpolation), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! *** Test log-lin Interpolation (log in x-axis)
    ! y = 2*log(x)
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = TWO * log(x1)
    y0 = 0.0_defReal
    ! @1.5
    val = 0.216395_defReal
#line 53 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 1.5_defReal, logLinInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! @5.0
    val = 8.09438_defReal
#line 57 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, logLinInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! *** Test lin-log Interpolation (lin in x-axis)
    ! y = 2*exp(x)
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = TWO * exp(x1)
    y0 = TWO * exp(x0)
    ! @2.5
    val = TWO * (exp(2.5_defReal) - exp(x0))
#line 67 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 2.5_defReal, linLogInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! @5.0
    val = TWO * (exp(5.0_defReal) - exp(x0))
#line 71 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, linLogInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"


    ! *** Test log-log Interpolation
    ! Slope /= -1
    ! y = x^7
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = x1**7
    y0 = x0**7
    ! @3.0
    val = (3.0_defReal**8 - x0**8) / 8.0_defReal
#line 83 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 3.0_defReal, logLogInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! @5.0
    val = (5.0_defReal**8 - x0**8) / 8.0_defReal
#line 87 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, logLogInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Slope == -1
    ! y = 1/x
    x1 = 5.0_defReal
    x0 = 1.0_defReal
    y1 = ONE / x1
    y0 = ONE / x0
    ! @2.5
    val = log(2.5_defReal/x0)
#line 97 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 2.5_defReal, logLogInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 97) )
  if (anyExceptions()) return
#line 98 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! @5.0
    val = log(5.0_defReal/x0)
#line 101 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, endf_bin_integral(x0, x1, y0, y1, 5.0_defReal, logLogInterpolation), TOL*val, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  end subroutine testEndfBinIntegral

  !!
  !! Test table with lin-lin interpolation
  !!
!@Test
  subroutine testLinLinTable()
    real(defReal),dimension(*),parameter :: x = [-1.0_defReal, 2.0_defReal, 2.5_defReal, 3.5_defReal]
    real(defReal),dimension(*),parameter :: y = [17.0_defReal, -2.0_defReal, 0.0_defReal, 1.5_defReal]
    type(endfTable)                      :: tab
    real(defReal),parameter :: TOL = 1.0E-6

    call tab % init(x, y)

    ! Check values at few points (including start and end of the table)
#line 117 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(17.0_defReal, tab % at(-1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 118 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(1.5_defReal, tab % at(3.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 119 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(-2.0_defReal, tab % at(2.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Check interpolated values
#line 122 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(32.0_defReal/3, tab % at(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 123 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(0.75_defReal, tab % at(3.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Test reload
    call tab % reloadY([ONE, ONE, ONE, ONE])
#line 127 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(ONE, tab % at(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 128 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(ONE, tab % at(3.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

  end subroutine testLinLinTable

  !!
  !! Test integration in lin-lin table
  !!
!@Test
  subroutine testLinLinIntegral()
    real(defReal), dimension(*), parameter :: x = [1.0_defReal, 2.0_defReal, 2.5_defReal, 3.0_defReal]
    real(defReal), dimension(*), parameter :: y = [0.0_defReal, 1.0_defReal, 0.5_defReal, 0.0_defReal]
    type(endfTable)                        :: tab
    real(defReal), dimension(4)            :: out, val, points
    real(defReal),parameter :: TOL = 1.0E-6

    call tab % init(x, y)

    ! Integral at grid points
    out = tab % integral(x)
    val = [0.0_defReal, 0.5_defReal, 0.875_defReal, 1.0_defReal]
#line 148 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, out, TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Integral of grid points
    points = [1.25_defReal, 1.75_defReal, 2.6_defReal, 2.9_defReal]
    val    = [3.125E-2_defReal, 0.28125_defReal, 0.92_defReal, 0.995_defReal]
    out = tab % integral(points)
#line 154 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, out, TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"


  end subroutine testLinLinIntegral


  !!
  !! Test table with a single interpolation region that is not lin-lin
  !!
!@Test
  subroutine testOneRegionInterpolation()
    real(defReal),dimension(*),parameter :: x = [-1.0_defReal, 2.0_defReal, 2.5_defReal, 3.5_defReal]
    real(defReal),dimension(*),parameter :: y = [17.0_defReal, -2.0_defReal, 0.0_defReal, 1.5_defReal]
    type(endfTable)                      :: tab
    real(defReal),parameter :: TOL = 1.0E-6

    ! Initialise with histogram interpolation
    call tab % init(x, y, [4], [histogramInterpolation])

    ! Check values at few points (including start and end of the table)
#line 174 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(17.0_defReal, tab % at(-1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 175 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(0.0_defReal, tab % at(3.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 176 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(-2.0_defReal, tab % at(2.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 176) )
  if (anyExceptions()) return
#line 177 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Check interpolated values
#line 179 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(17.0_defReal, tab % at(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 180 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(0.0_defReal, tab % at(3.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

  end subroutine testOneRegionInterpolation

  !!
  !! Test integration with single (log-log) intrepolation egion
  !!
  !! Polynomial functions
  !!  x^2 if x < 1
  !!  x^(-1) if x >= 1
  !!
!@Test
  subroutine testOneRegionIntegral()
    real(defReal), dimension(6), parameter :: x = [0.01_defReal, 0.5_defReal, 1.0_defReal,&
                                                   2.0_defReal, 4.0_defReal, 8.0_defReal]
    real(defReal), dimension(6), parameter :: y = [1.0E-4_defReal, 0.25_defReal, 1.0_defReal, &
                                                   0.5_defReal, 0.25_defReal, 0.125_defReal ]
    type(endfTable)                        :: tab
    real(defReal), dimension(6)            :: out, val, p
    real(defReal),parameter                :: TOL = 1.0E-6

    call tab % init(x, y, [6], [loglogInterpolation])

    ! Integral at grid points
    val(1) = ZERO
    val(2) = (x(2)**3 - x(1)**3) / 3.0_defReal
    val(3) = (ONE - x(1)**3) / 3.0_defReal
    val(4) = val(3) + log(x(4)/x(3))
    val(5) = val(3) + log(x(5)/x(3))
    val(6) = val(3) + log(x(6)/x(3))
    out = tab % integral(x)
#line 211 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, out, TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 211) )
  if (anyExceptions()) return
#line 212 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Integral in between
    p = [0.1_defReal, 0.75_defReal, 0.9_defReal, 1.1_defReal, 1.4_defReal, 6.0_defReal]
    val(1) = (p(1)**3 - x(1)**3) / 3.0_defReal
    val(2) = (p(2)**3 - x(1)**3) / 3.0_defReal
    val(3) = (p(3)**3 - x(1)**3) / 3.0_defReal
    val(4) = ONE/3.0_defReal + log(p(4)/x(3)) - x(1)**3 / 3.0_defReal
    val(5) = ONE/3.0_defReal + log(p(5)/x(3)) - x(1)**3 / 3.0_defReal
    val(6) = ONE/3.0_defReal + log(p(6)/x(3)) - x(1)**3 / 3.0_defReal
    out = tab % integral(p)
#line 222 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, out, TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

  end subroutine testOneRegionIntegral

  !!
  !! Test table with multiple interpolation regions
  !!
!@Test
  subroutine testMultipleInterRegions()
    real(defReal),dimension(*),parameter :: x = [-1.0_defReal, 2.0_defReal, 2.5_defReal, 2.5_defReal, 3.5_defReal]
    real(defReal),dimension(*),parameter :: y = [17.0_defReal, 2.0_defReal, 0.0_defReal, 1.0_defReal, 1.5_defReal]
    type(endfTable)                      :: tab
    real(defReal),parameter :: TOL = 1.0E-6

    ! Initialise with linLog, logLin and logLog interpolations
    call tab % init(x, y, [2,4,5], [loglinInterpolation, linLogInterpolation, logLogInterpolation])

    ! Check values at few points (including start and end of the table)
#line 240 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(17.0_defReal, tab % at(-1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 241 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(1.5_defReal, tab % at(3.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
#line 242 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(2.0_defReal, tab % at(2.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Check y-log x-lin interpolation
#line 245 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(8.329954_defReal ,tab % at(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Check y-lin x-log interpolation
#line 248 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(1.5627016_defReal, tab % at(2.1_defReal),TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 248) )
  if (anyExceptions()) return
#line 249 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Check y-log x-log interpolation
#line 251 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(1.346458993_defReal, tab % at(3.2_defReal), TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 251) )
  if (anyExceptions()) return
#line 252 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

  end subroutine testMultipleInterRegions

  !!
  !! Test Integration with multiple interpolation regions
  !!  3 regions hist, lin and
  !!   Histogram for x < 2.0
  !!   Lin-Lin   for x < 3.0  y = x - 2
  !!   Log-Log   for x > 3.0  y = 3/x
  !!
  !!
!@Test
  subroutine testMultiRegionIntegral()
    real(defReal), dimension(7), parameter :: x = [1.0_defReal, 2.0_defReal, 2.0_defReal, &
                                                   2.5_defReal, 3.0_defReal, 4.0_defReal, 7.0_defReal]
    real(defReal), dimension(7)            :: y
    type(endfTable)                        :: tab
    real(defReal), dimension(7)            :: out, val
    real(defReal), dimension(3)            :: p, out2, val2
    real(defReal),parameter                :: TOL = 1.0E-6

    ! Calculate y values
    y(1) = 2.0_defReal
    y(2) = -4.0_defReal  ! Irrelevant for hist-interp
    y(3:5) = x(3:5) - 2.0_defReal
    y(6:7) = 3.0_defReal / x(6:7)

    ! Create table
    call tab % init(x, y, [2, 5, 7], [histogramInterpolation, linLinInterpolation, loglogInterpolation])

    ! Integral at grid points
    val(1) = ZERO
    val(2) = TWO
    val(3) = TWO
    val(4) = TWO + 0.125_defReal
    val(5) = TWO + HALF
    val(6) = TWO + HALF + 3*log(x(6)/x(5))
    val(7) = TWO + HALF + 3*log(x(7)/x(5))
    out = tab % integral(x)
#line 291 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"
  call assertEqual(val, out, TOL, &
 & location=SourceLocation( &
 & 'endfTable_test.f90', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/endfTable_test.f90"

    ! Integral out grid points
    p = [1.5_defReal, 2.9_defReal, 3.2_defReal]
    val2(1) = ONE
    val2(2) = TWO + 0.405_defReal
    val2(3) = TWO + HALF + 3*log(p(3)/x(5))

  end subroutine testMultiRegionIntegral


end module endfTable_test

module WrapendfTable_test
   use pFUnit_mod
   use endfTable_test
   implicit none
   private

contains


end module WrapendfTable_test

function endfTable_test_suite() result(suite)
   use pFUnit_mod
   use endfTable_test
   use WrapendfTable_test
   type (TestSuite) :: suite

   suite = newTestSuite('endfTable_test_suite')

   call suite%addTest(newTestMethod('testEndfBinIntegral', testEndfBinIntegral))

   call suite%addTest(newTestMethod('testLinLinTable', testLinLinTable))

   call suite%addTest(newTestMethod('testLinLinIntegral', testLinLinIntegral))

   call suite%addTest(newTestMethod('testOneRegionInterpolation', testOneRegionInterpolation))

   call suite%addTest(newTestMethod('testOneRegionIntegral', testOneRegionIntegral))

   call suite%addTest(newTestMethod('testMultipleInterRegions', testMultipleInterRegions))

   call suite%addTest(newTestMethod('testMultiRegionIntegral', testMultiRegionIntegral))


end function endfTable_test_suite

