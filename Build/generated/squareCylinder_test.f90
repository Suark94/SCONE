module squareCylinder_test
  use numPrecision
  use universalVariables
  use dictionary_class,     only : dictionary
  use squareCylinder_class, only : squareCylinder
  use pfUnit_mod

  implicit none

  !!
  !! Test parameter wrapper around AN INTEGER (bit of boilerplate)
  !!
  !!
  !@testParameter(constructor=newParam)
  type, extends (AbstractTestParameter) :: dirParam
     integer(shortInt) :: dir
  contains
     procedure :: toString
  end type dirParam

  !!
  !! Cylinder Test case
  !!
  !@TestCase(constructor=newTestCase)
    type, extends(ParameterizedTestCase) :: test_squareCylinder
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      type(squareCylinder)            :: surf
    contains
      procedure :: tearDown
    end type test_squareCylinder

contains

  !!
  !! Test parameter constructor
  !!
  function newParam(i) result(param)
     integer(shortInt), intent(in) :: i
     type (dirParam) :: param

     param % dir = i

  end function newParam

  !!
  !! Print parameter to string for more verbose description
  !!
  function toString(this) result(string)
     class (dirParam), intent(in) :: this
     character(:), allocatable :: string

     select case(this % dir)
       case(X_AXIS)
         string = 'xSquareCylinder'
       case(Y_AXIS)
         string = 'ySquareCylinder'
       case(Z_AXIS)
         string = 'zSquareCylinder'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_squareCylinder test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !!            axis p1   p2
  !! Origin     2.0, 1.0, 2.0
  !! Halfwidths      2.0, 3.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_squareCylinder)  :: tst
    type(dictionary)           :: dict
    character(nameLen)         :: type
    real(defReal), dimension(3) :: origin, hw

    ! Select type of squareCylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xSquareCylinder'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'ySquareCylinder'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zSquareCylinder'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"

    end select

    ! Set origin & halfwidth
    origin = TWO
    origin(tst % plane(1)) = ONE
    origin(tst % plane(2)) = TWO

    hw = ZERO
    hw(tst % plane(1)) = TWO
    hw(tst % plane(2)) = 3.0_defReal

    ! Build surface
    call dict % init(4)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', origin)
    call dict % store('halfwidth', hw)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_squareCylinder), intent(inout) :: this

    call this % surf % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Proper tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test Misc functionality
  !!
  !! Directions must be given as integers for pFUnit parser to work
  !!
!@Test(cases = [1, 2, 3])
  subroutine testMisc(this)
    class(test_squareCylinder), intent(inout) :: this
    real(defReal), dimension(6) :: aabb, ref
    character(nameLen)          :: name
    real(defReal), parameter    :: TOL = 1.0E-6_defReal

    ! Test ID
#line 149 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(75, this % surf % id(), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Change ID
    call this % surf % setID(1)
#line 153 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(1, this % surf % id(), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xSquareCylinder'
      case(Y_AXIS)
        name = 'ySquareCylinder'
      case(Z_AXIS)
        name = 'zSquareCylinder'
    end select
#line 164 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(name, this % surf % myType(), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Bounding Box
    ref = [ZERO, ZERO, ZERO, 4.0_defReal, 4.0_defReal, 4.0_defReal]
    ref(this % plane) = [-1.0_defReal, -1.0_defReal]
    ref(this % plane+3) = [3.0_defReal, 5.0_defReal]
    ref(this % axis) = -INF
    ref(this % axis+3) = INF

    aabb = this % surf % boundingBox()
#line 174 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
!@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_squareCylinder), intent(inout) :: this
    integer(shortInt), dimension(6)  :: BC
    integer(shortInt)                :: ax, p1, p2
    real(defReal), dimension(3)      :: r, u, r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Set BC
    BC = VACUUM_BC
    BC(p1*2) = REFLECTIVE_BC
    BC(p1*2-1) = VACUUM_BC
    BC(p2*2) = PERIODIC_BC
    BC(p2*2-1) = PERIODIC_BC

    call this % surf % setBC(BC)

    ! Explicit BC

    ! Vacuum face
    r([ax, p1, p2]) = [ZERO, -1.0_defReal, 3.0_defReal]
    u([ax, p1, p2]) = [ZERO, -ONE, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
#line 211 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 211) )
  if (anyExceptions()) return
#line 212 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 212 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 212) )
  if (anyExceptions()) return
#line 213 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Reflection face
    r([ax, p1, p2]) = [ZERO, 3.0_defReal, 3.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    r_ref = r
    u_ref = u
    u_ref(p1) = -ONE
    call this % surf % explicitBC(r, u)
#line 221 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 222 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Periodic face
    r([ax, p1, p2]) = [ZERO, 2.0_defReal, 5.0_defReal]
    u([ax, p1, p2]) = [ZERO, ZERO, ONE]
    r_ref([ax, p1, p2]) = [ZERO, 2.0_defReal, -1.0_defReal]
    u_ref = u
    call this % surf % explicitBC(r, u)
#line 230 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 231 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 231) )
  if (anyExceptions()) return
#line 232 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Corner
    r([ax, p1, p2]) = [ZERO, 3.0_defReal, -1.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)
    r_ref([ax, p1, p2]) = [ZERO, 3.0_defReal, 5.0_defReal]
    u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)
    call this % surf % explicitBC(r, u)
#line 241 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 242 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Transform BC
    r([ax, p1, p2]) = [ZERO, 20.0_defReal, -13.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, -14.0_defReal, 5.0_defReal]
    u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)

    call this % surf % transformBC(r, u)
#line 254 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 254) )
  if (anyExceptions()) return
#line 255 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 255 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 255) )
  if (anyExceptions()) return
#line 256 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Try the corner
    r([ax, p1, p2]) = [ZERO, 3.0_defReal, -1.0_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, -ONE]
    u = u/norm2(u)

    r_ref([ax, p1, p2]) = [ZERO, 3.0_defReal, 5.0_defReal]
    u_ref([ax, p1, p2]) = [ZERO, -ONE, -ONE]
    u_ref = u_ref/norm2(u_ref)

    call this % surf % transformBC(r, u)
#line 267 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 267) )
  if (anyExceptions()) return
#line 268 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 268 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 268) )
  if (anyExceptions()) return
#line 269 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
!@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_squareCylinder), intent(inout) :: this
    integer(shortInt)                         :: ax, p1, p2
    real(defReal), dimension(3)               :: r, u, u2
    real(defReal)                             :: eps

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! ** Well inside the square cylinder
    r = ZERO
    u = ZERO
    u(p2) = ONE
    r([p1, p2]) = [2.0_defReal, 0.0_defReal]
#line 292 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 292) )
  if (anyExceptions()) return
#line 293 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Different Quadrant
    r([p1, p2]) = [-0.5_defReal, 2.0_defReal]
#line 296 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 296) )
  if (anyExceptions()) return
#line 297 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! ** Well outside the square cylinder
    r([p1, p2]) = [-1.5_defReal, 2.0_defReal]
#line 300 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 300) )
  if (anyExceptions()) return
#line 301 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Diffrent Quadrant
    r([p1, p2]) = [0.5_defReal, 5.2_defReal]
#line 304 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 304) )
  if (anyExceptions()) return
#line 305 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! ** Proximity of the surface
    r([p1, p2]) = [-1.0_defReal, 3.0_defReal]
    u = ZERO
    u(p1) = -ONE
#line 310 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 310) )
  if (anyExceptions()) return
#line 311 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 311 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, -u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 311) )
  if (anyExceptions()) return
#line 312 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Inside within Surface Tolerance
    eps = -HALF * SURF_TOL
#line 315 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 315) )
  if (anyExceptions()) return
#line 316 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! A bit above Surface Tolerance
    eps = -1.0001_defReal * SURF_TOL
#line 319 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 319) )
  if (anyExceptions()) return
#line 320 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Parallel to the surface
    u2 = ZERO
    u2(p2) = ONE
    eps = HALF * SURF_TOL

#line 326 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u2), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 326) )
  if (anyExceptions()) return
#line 327 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 327 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r - eps*u, u2), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 327) )
  if (anyExceptions()) return
#line 328 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! At the surface in diffrent quadrant



  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
!@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_squareCylinder), intent(inout) :: this
    integer(shortInt)                   :: ax, p1, p2
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    !** Outside the square cylinder
    r = ZERO
    r(p1) = -2.0_defReal
    r(p2) = ZERO

    ! Direct hit
    u([ax, p1, p2]) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    ref = SQRT2
#line 360 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 360) )
  if (anyExceptions()) return
#line 361 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Moving away
#line 363 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 363) )
  if (anyExceptions()) return
#line 364 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Oblique Hit
    u([ax, p1, p2]) = [ONE, ONE, ONE]
    u = u/norm2(u)
    ref = sqrt(3.0_defReal)
#line 369 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 369) )
  if (anyExceptions()) return
#line 370 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Miss
    u([ax, p1, p2]) = [ONE, ONE, -TWO]
    u = u/norm2(u)
#line 374 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 374) )
  if (anyExceptions()) return
#line 375 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Corner skim
    ! Use dirty values
    r([ax, p1, p2]) = [ZERO, -1.3_defReal, 0.3_defReal]
    u([ax, p1, p2])  = [ZERO, 0.3_defReal, -1.3_defReal]
    u = u/norm2(u)
#line 381 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 381) )
  if (anyExceptions()) return
#line 382 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Parallel
    r([ax, p1, p2]) = [ZERO, -1.3_defReal, 0.3_defReal]
    u([ax, p1, p2])  = [ZERO, ZERO, ONE]
    u = u/norm2(u)
#line 387 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 387) )
  if (anyExceptions()) return
#line 388 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! ** At the surface
    r([ax, p1, p2]) = [ZERO, -1.0_defReal, 0.5_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal
#line 393 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 393) )
  if (anyExceptions()) return
#line 394 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Outside within Surface Tolerance
    r([ax, p1, p2]) = [ZERO, -1.0_defReal - HALF * SURF_TOL, 0.5_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal + HALF * SURF_TOL
#line 399 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 399) )
  if (anyExceptions()) return
#line 400 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Inside within Surface Tolerance
    r([ax, p1, p2]) = [ZERO, -1.0_defReal + HALF * SURF_TOL, 0.5_defReal]
    u([ax, p1, p2]) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal - HALF * SURF_TOL
#line 405 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 405) )
  if (anyExceptions()) return
#line 406 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! ** Inside
    r([ax, p1, p2]) = [9.0_defReal, 1.15_defReal, 1.0_defReal]
    u([ax, p1, p2]) = [ZERO, ZERO, ONE]
    ref = 4.0_defReal
#line 411 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 411) )
  if (anyExceptions()) return
#line 412 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

    ! Other direction
    ref = 2.0_defReal
#line 415 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, -u), ref * TOL, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 415) )
  if (anyExceptions()) return
#line 416 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

  end subroutine testDistance

  !!
  !! Test Edge Cases
  !!
  !! Test unlikley cases to make sure that halfspace + distance
  !! procedures allow particle to escape
  !!
!@Test(cases=[1, 2, 3])
  subroutine testEdgeCases(this)
    class(test_squareCylinder), intent(inout) :: this
    real(defReal), dimension(3) :: r, u, u2
    integer(shortInt)           :: ax, p1, p2
    real(defReal)               :: eps, d
    logical(defBool)            :: hs

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! ** Corner
    ! * Particle is almost at the corner
    !   Either it is outside or can escape with a short movment in next step
    eps =  5.0_defReal * epsilon(eps)
    r([ax, p1, p2]) = [ZERO, -1.0_defReal+eps, -1.0_defReal+eps]
    u([ax, p1, p2]) = [ZERO, TWO, -ONE]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
#line 448 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 448) )
  if (anyExceptions()) return
#line 449 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 449 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 449) )
  if (anyExceptions()) return
#line 450 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
    end if

    ! Try in other direction
    u([ax, p1, p2]) = [ZERO, -ONE, TWO]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
#line 458 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 458) )
  if (anyExceptions()) return
#line 459 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 459 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 459) )
  if (anyExceptions()) return
#line 460 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
    end if

    ! Try asymertic corner
    ! Point is not exactly at the diagonal
    r([ax, p1, p2]) = [ZERO, -1.0_defReal+eps, -1.0_defReal+eps*TWO]
    u([ax, p1, p2]) = [ZERO, TWO, -ONE]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
#line 470 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 470) )
  if (anyExceptions()) return
#line 471 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 471 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 471) )
  if (anyExceptions()) return
#line 472 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
    end if

    ! Try in other direction
    r([ax, p1, p2]) = [ZERO, -1.0_defReal+eps*TWO, -1.0_defReal+eps]
    u([ax, p1, p2]) = [ZERO, -ONE, TWO]
    u = u/norm2(u)
    hs = this % surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = this % surf % distance(r, u)
#line 481 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 481) )
  if (anyExceptions()) return
#line 482 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
#line 482 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 482) )
  if (anyExceptions()) return
#line 483 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
    end if

  end subroutine testEdgeCases

  !!
  !! Test encountered problems
  !!
  !! Contains test related to bugs found at some point
  !! TODO: Move some of this tests to main test procedures
  !!
!@Test(cases=[1])
  subroutine test_problems(this)
    class(test_squareCylinder), intent(inout) :: this ! Ignore this
    type(squareCylinder) :: surf
    type(dictionary)     :: dict
    real(defReal), dimension(3) :: r, u

    call dict % init(5)
    call dict % store('type','zSquareCylinder')
    call dict % store('id', 7)
    call dict % store('origin', [ZERO, ZERO, ZERO])
    call dict % store('halfwidth', [8.00_defReal, 1.26_defReal, 0.0_defReal])
    call surf % init(dict)

    r = [-7.63_defReal, 1.26_defReal, ZERO]
    u = [ZERO, -ONE, ZERO]
#line 509 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"
  call assertFalse(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'squareCylinder_test.f90', &
 & 509) )
  if (anyExceptions()) return
#line 510 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/squareCylinder_test.f90"

  end subroutine test_problems



end module squareCylinder_test

module WrapsquareCylinder_test
   use pFUnit_mod
   use squareCylinder_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_squareCylinder) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use squareCylinder_test
        class (test_squareCylinder), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod, testParameter) result(aTest)
#ifdef INTEL_13
      use pfunit_mod, only: testCase
#endif
      type (WrapUserTestCase) :: aTest
#ifdef INTEL_13
      target :: aTest
      class (WrapUserTestCase), pointer :: p
#endif
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      type (dirParam), intent(in) :: testParameter
      aTest%test_squareCylinder = newTestCase(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WrapsquareCylinder_test

function squareCylinder_test_suite() result(suite)
   use pFUnit_mod
   use squareCylinder_test
   use WrapsquareCylinder_test
   type (TestSuite) :: suite

   type (dirParam), allocatable :: testParameters(:)
   type (dirParam) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('squareCylinder_test_suite')

   cases = [1, 2, 3]
   testParameters = [(newParam(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testMisc', testMisc, testParameter))
   end do

   cases = [1,2,3]
   testParameters = [(newParam(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testBC', testBC, testParameter))
   end do

   cases = [1,2,3]
   testParameters = [(newParam(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testHalfspace', testHalfspace, testParameter))
   end do

   cases = [1, 2, 3]
   testParameters = [(newParam(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testDistance', testDistance, testParameter))
   end do

   cases = [1, 2, 3]
   testParameters = [(newParam(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testEdgeCases', testEdgeCases, testParameter))
   end do

   cases = [1]
   testParameters = [(newParam(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('test_problems', test_problems, testParameter))
   end do


end function squareCylinder_test_suite

