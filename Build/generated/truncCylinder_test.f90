module truncCylinder_test
  use numPrecision
  use universalVariables
  use dictionary_class,     only : dictionary
  use truncCylinder_class,  only : truncCylinder
  use pfUnit_mod

  implicit none

  !!
  !! Test parameter wrapper around AN INTEGER (bit of boilerplate)
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
    type, extends(ParameterizedTestCase) :: test_truncCylinder
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      type(truncCylinder)            :: surf
    contains
      procedure :: tearDown
    end type test_truncCylinder

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
         string = 'xTruncCylinder'
       case(Y_AXIS)
         string = 'yTruncCylinder'
       case(Z_AXIS)
         string = 'zTruncCylinder'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_truncCylinder test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !!            axis p1   p2
  !! Origin     2.0, 1.0, 2.0
  !! Halfwidth  1.5
  !! Radius     2.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in)  :: dir
    type(test_truncCylinder)    :: tst
    type(dictionary)            :: dict
    character(nameLen)          :: type
    real(defReal), dimension(3) :: origin
    real(defReal)               :: hw, r

    ! Select type of truncCylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xTruncCylinder'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yTruncCylinder'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zTruncCylinder'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"

    end select

    ! Set origin & halfwidth
    origin = TWO
    origin(tst % plane(1)) = ONE
    origin(tst % plane(2)) = TWO

    hw = 1.5_defReal
    r = 2.0_defReal

    ! Build surface
    call dict % init(5)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', origin)
    call dict % store('halfwidth', hw)
    call dict % store('radius', r)
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_truncCylinder), intent(inout) :: this

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
    class(test_truncCylinder), intent(inout) :: this
    real(defReal), dimension(6) :: aabb, ref
    character(nameLen)          :: name
    integer(shortInt)                :: ax, p1, p2
    real(defReal), parameter    :: TOL = 1.0E-6_defReal

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Test ID
#line 156 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(75, this % surf % id(), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Change ID
    call this % surf % setID(1)
#line 160 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(1, this % surf % id(), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xTruncCylinder'
      case(Y_AXIS)
        name = 'yTruncCylinder'
      case(Z_AXIS)
        name = 'zTruncCylinder'
    end select
#line 171 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(name, this % surf % myType(), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Bounding Box
    ref([ax, p1, p2]) = [0.5_defReal, -1.0_defReal, 0.0_defReal]
    ref([ax, p1, p2]+3) = [3.5_defReal, 3.0_defReal, 4.0_defReal]
    aabb = this % surf % boundingBox()
#line 177 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
!@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_truncCylinder), intent(inout) :: this
    integer(shortInt), dimension(6)  :: BC
    integer(shortInt)                :: ax, p1, p2
    integer(shortInt), dimension(3)  :: pe
    real(defReal), dimension(3)      :: r, u, r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)
    pe = [ax, p1, p2] ! Common permutation

    ! *** CASE 1 VACUUM + REFLECTIVE
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC])

    ! ** Explicit BC
    ! Vacuum face
    r(pe) = [0.5_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [-ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    call this % surf % explicitBC(r, u)
#line 209 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 210 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 210) )
  if (anyExceptions()) return
#line 211 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Reflective face
    r(pe) = [3.5_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    u_ref(ax) = -u(ax)
    call this % surf % explicitBC(r, u)
#line 219 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 219) )
  if (anyExceptions()) return
#line 220 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 220 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 220) )
  if (anyExceptions()) return
#line 221 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** Transform BC
    r(pe) = [30.0_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref(pe) = [-23.0_defReal, 0.0_defReal, 2.3_defReal]
    u_ref(pe) = [-ONE, ZERO, ZERO]
    call this % surf % transformBC(r, u)
#line 228 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 228) )
  if (anyExceptions()) return
#line 229 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 229 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 229) )
  if (anyExceptions()) return
#line 230 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! **** CASE 2 PERIODIC
    call this % surf % setBC([PERIODIC_BC, PERIODIC_BC])

    ! Explicit BC
    ! Periodic face
    r(pe) = [3.5_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref(pe) = [0.5_defReal, 0.0_defReal, 2.3_defReal]
    u_ref = u
    call this % surf % explicitBC(r, u)
#line 241 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 242 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** Transform BC
    r(pe) = [12.0_defReal, 0.0_defReal, 2.3_defReal]
    u(pe) = [ONE, ZERO, ZERO]
    r_ref(pe) = [3.0_defReal, 0.0_defReal, 2.3_defReal]
    u_ref = u     
    call this % surf % transformBC(r, u)
#line 250 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 250) )
  if (anyExceptions()) return
#line 251 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 251 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 251) )
  if (anyExceptions()) return
#line 252 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
!@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_truncCylinder), intent(inout) :: this
    integer(shortInt)                         :: ax, p1, p2
    integer(shortInt), dimension(3)           :: pe
    real(defReal), dimension(3)               :: r, u, u2
    real(defReal)                             :: eps

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)
    pe = [ax, p1, p2] ! Common permutation

    ! ** Outside
    r(pe) = [ 5.0_defReal, 3.0_defReal, 4.0_defReal ]
    u(pe) = [ZERO, ONE, ZERO]
#line 275 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 275) )
  if (anyExceptions()) return
#line 276 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Diffrent position
    r(pe) = [1.3_defReal, 1.0_defReal, 5.0_defReal]
#line 279 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 279) )
  if (anyExceptions()) return
#line 280 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! One more
    r(pe) = [1.3_defReal, -2.0_defReal, 3.3_defReal]
#line 283 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 283) )
  if (anyExceptions()) return
#line 284 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** Inside the surface
    ! Close to plane
    r(pe) = [3.3_defReal, 1.1_defReal, 2.2_defReal]
#line 288 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 288) )
  if (anyExceptions()) return
#line 289 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Close to cylinder
    r(pe) = [2.1_defReal, -0.7_defReal, 2.1_defReal]
#line 292 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 292) )
  if (anyExceptions()) return
#line 293 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** Close to surface
    ! Cylinder
    r(pe) = [3.3_defReal, -1.0_defReal, 2.0_defReal]
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
#line 299 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 300 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, -u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 300) )
  if (anyExceptions()) return
#line 301 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Outside within surface tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
#line 307 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 307) )
  if (anyExceptions()) return
#line 308 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Outside outside surface tol
    eps = 1.001_defReal * SURF_TOL
    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
#line 314 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 314) )
  if (anyExceptions()) return
#line 315 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Parallel within tolerance
    eps = HALF * SURF_TOL
    u(pe) = [ONE, ZERO, ZERO]

    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
#line 321 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 321) )
  if (anyExceptions()) return
#line 322 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    r(pe) = [3.3_defReal, -1.0_defReal+eps, 2.0_defReal]
#line 324 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 324) )
  if (anyExceptions()) return
#line 325 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Plane
    r(pe) = [3.5_defReal, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ZERO, ONE]
    u = u/norm2(u)
#line 330 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 330) )
  if (anyExceptions()) return
#line 331 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 331 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, -u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 331) )
  if (anyExceptions()) return
#line 332 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Outside within surface tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ONE, ZERO]
    u = u/norm2(u)
#line 338 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 338) )
  if (anyExceptions()) return
#line 339 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Outside outside surface tol
    eps = 1.001_defReal * SURF_TOL
    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ONE, ZERO]
    u = u/norm2(u)
#line 345 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 345) )
  if (anyExceptions()) return
#line 346 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Parallel within tolerance
    eps = HALF * SURF_TOL
    u(pe) = [ZERO, -ONE, ONE]
    u = u/norm2(u)

    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
#line 353 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 353) )
  if (anyExceptions()) return
#line 354 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    r(pe) = [3.5_defReal-eps, 1.3_defReal, 1.8_defReal]
#line 356 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 356) )
  if (anyExceptions()) return
#line 357 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
!@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_truncCylinder), intent(inout) :: this
    integer(shortInt)                   :: ax, p1, p2
    integer(shortInt), dimension(3)     :: pe
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: ref, eps
    real(defReal), parameter :: TOL = 1.0E-7

    ! Get axis and diffrent planar directions
    ax = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)
    pe = [ax, p1, p2] ! Common permutation

    ! ** Outside
    ! *Plane hits
    r(pe) = [5.0_defReal, 1.0_defReal, 2.0_defReal ]

    ! Direct hit
    u(pe) = [-ONE, ZERO, ZERO]
    ref = 1.5_defReal
#line 385 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 385) )
  if (anyExceptions()) return
#line 386 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Oblique hit
    u(pe) = [-ONE, -ONE, ZERO]
    u = u/norm2(u)
    ref = 1.5_defReal * SQRT2
#line 391 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 391) )
  if (anyExceptions()) return
#line 392 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Miss
    u(pe) = [-ONE, -TWO, -TWO]
    u = u/norm2(u)
#line 396 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 396) )
  if (anyExceptions()) return
#line 397 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Parallel miss
    u(pe) = [ZERO, -TWO, -TWO]
    u = u/norm2(u)
#line 401 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 401) )
  if (anyExceptions()) return
#line 402 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Wrong direction
    u(pe) = [ONE, -ONE, -ONE]
    u = u/norm2(u)
#line 406 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 406) )
  if (anyExceptions()) return
#line 407 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! * Cylinder hits
    r(pe) = [1.5_defReal, -2.0_defReal, 2.0_defReal ]

    ! Direct hit
    u(pe) = [ZERO, ONE, ZERO]
    ref = 1.0_defReal
#line 414 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 414) )
  if (anyExceptions()) return
#line 415 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Oblique hit 1
    u(pe) = [ONE, ONE, ZERO]
    u = u/norm2(u)
    ref = 1.0_defReal * SQRT2
#line 420 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 420) )
  if (anyExceptions()) return
#line 421 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Oblique hit 2
    u(pe) = [ZERO, 1.5_defReal, sqrt(7.0_defReal)* HALF]
    ref = norm2(u)
    u = u/norm2(u)
#line 426 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 426) )
  if (anyExceptions()) return
#line 427 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Miss
    u(pe) = [ONE, ONE, TWO]
    u = u/norm2(u)
#line 431 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 431) )
  if (anyExceptions()) return
#line 432 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Parallel miss
    u(pe) = [ONE, ZERO, ZERO]
#line 435 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 435) )
  if (anyExceptions()) return
#line 436 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Wrong direction
    u(pe) = [ONE, -ONE, ZERO]
#line 439 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 439) )
  if (anyExceptions()) return
#line 440 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** At the surface - Plane
    ! Outside within surf tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.5_defReal+eps, 1.3_defReal, 1.8_defReal]
    u(pe) = [-ONE, ZERO, ZERO]
    ref = 3.0_defReal
#line 447 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 447) )
  if (anyExceptions()) return
#line 448 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Inside going outside
    r(pe) = [3.5_defReal-eps, 1.3_defReal, 1.8_defReal]
#line 451 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 451) )
  if (anyExceptions()) return
#line 452 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** At the surface - Cylinder
    ! Outside within surface tolerance
    eps = HALF * SURF_TOL
    r(pe) = [3.3_defReal, -1.0_defReal-eps, 2.0_defReal]
    u(pe) = [ZERO, ONE, ZERO]
    ref = 4.0_defReal
#line 459 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 459) )
  if (anyExceptions()) return
#line 460 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Inside going outside
    r(pe) = [3.3_defReal, -1.0_defReal+eps, 2.0_defReal]
#line 463 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 463) )
  if (anyExceptions()) return
#line 464 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! ** Inside
    r(pe) = [3.0_defReal, 0.0_defReal, 2.0_defReal]

    ! Oblique hits
    u(pe) = [-ONE, -ONE, ZERO]
    u = u/norm2(u)
    ref = SQRT2
#line 472 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 472) )
  if (anyExceptions()) return
#line 473 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ref = HALF * SQRT2
#line 475 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, -u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 475) )
  if (anyExceptions()) return
#line 476 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Parallel Planes
    u(pe) = [ZERO, -ONE, ZERO]
    ref = 1.0_defReal
#line 480 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 480) )
  if (anyExceptions()) return
#line 481 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

    ! Parallel Cylinder
    u(pe) = [-ONE, ZERO, ZERO]
    ref = 2.5_defReal
#line 485 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), ref * TOL, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 485) )
  if (anyExceptions()) return
#line 486 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"

  end subroutine testDistance

   !!
   !! Test Edge Cases
   !!
   !! Test unlikley cases to make sure that halfspace + distance
   !! procedures allow particle to escape
   !!
 !@Test(cases=[1, 2, 3])
   subroutine testEdgeCases(this)
     class(test_truncCylinder), intent(inout) :: this
     real(defReal), dimension(3)    :: r, u, u2
    integer(shortInt), dimension(3) :: pe
     integer(shortInt)              :: ax, p1, p2
     real(defReal)                  :: eps, d
     logical(defBool)               :: hs

     ! Get axis and diffrent planar directions
     ax = this % axis
     p1 = this % plane(1)
     p2 = this % plane(2)
     pe = [ax, p1, p2] ! Common permutation

     ! ** Try to escape the surface beeing in a corner
     ! * Particle is almost at the corner
     !   Either it is outside or can escape with a short movment in next step
     eps =  5.0_defReal * epsilon(eps)
     r([ax, p1, p2]) = [3.5_defReal-eps, -1.0_defReal+eps, 2.0_defReal]
     u([ax, p1, p2]) = [-TWO, -ONE, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
#line 520 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 520) )
  if (anyExceptions()) return
#line 521 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 521 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 521) )
  if (anyExceptions()) return
#line 522 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
     end if

     ! Try other direction
     u([ax, p1, p2]) = [ONE, TWO, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
#line 530 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 530) )
  if (anyExceptions()) return
#line 531 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 531 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 531) )
  if (anyExceptions()) return
#line 532 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
     end if

     ! * Asymetric corner
     eps =  5.0_defReal * epsilon(eps)
     r([ax, p1, p2]) = [3.5_defReal-eps, -1.0_defReal+TWO*eps, 2.0_defReal]
     u([ax, p1, p2]) = [-TWO, -ONE, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
#line 542 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 542) )
  if (anyExceptions()) return
#line 543 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 543 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 543) )
  if (anyExceptions()) return
#line 544 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
     end if

     ! Try other direction
     r([ax, p1, p2]) = [3.5_defReal-TWO*eps, -1.0_defReal+eps, 2.0_defReal]
     u([ax, p1, p2]) = [ONE, TWO, ZERO]
     u = u/norm2(u)
     hs = this % surf % halfspace(r, u)
     if (.not.hs) then ! Perform small movment
       d = this % surf % distance(r, u)
#line 553 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 553) )
  if (anyExceptions()) return
#line 554 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
#line 554 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'truncCylinder_test.f90', &
 & 554) )
  if (anyExceptions()) return
#line 555 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/truncCylinder_test.f90"
     end if

   end subroutine testEdgeCases

end module truncCylinder_test

module WraptruncCylinder_test
   use pFUnit_mod
   use truncCylinder_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_truncCylinder) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use truncCylinder_test
        class (test_truncCylinder), intent(inout) :: this
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
      aTest%test_truncCylinder = newTestCase(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WraptruncCylinder_test

function truncCylinder_test_suite() result(suite)
   use pFUnit_mod
   use truncCylinder_test
   use WraptruncCylinder_test
   type (TestSuite) :: suite

   type (dirParam), allocatable :: testParameters(:)
   type (dirParam) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('truncCylinder_test_suite')

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


end function truncCylinder_test_suite

