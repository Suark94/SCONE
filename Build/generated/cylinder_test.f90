module cylinder_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use cylinder_class,    only : cylinder
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
    type, extends(ParameterizedTestCase) :: test_cylinder
      integer(shortInt)               :: axis
      integer(shortInt), dimension(2) :: plane
      type(cylinder)                  :: surf
    contains
      procedure :: tearDown
    end type test_cylinder

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
         string = 'xCylinder'
       case(Y_AXIS)
         string = 'yCylinder'
       case(Z_AXIS)
         string = 'zCylinder'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_cylinder test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !! Origin 2.0, 2.0, 2.0
  !! Radius 2.0
  !! ID 75
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_cylinder)        :: tst
    type(dictionary)      :: dict
    character(nameLen)    :: type

    ! Select type of cylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        tst % plane = [Y_AXIS, Z_AXIS]
        type = 'xCylinder'

      case(Y_AXIS)
        tst % axis = Y_AXIS
        tst % plane = [X_AXIS, Z_AXIS]
        type = 'yCylinder'

      case(Z_AXIS)
        tst % axis = Z_AXIS
        tst % plane = [X_AXIS, Y_AXIS]
        type = 'zCylinder'

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"

    end select

    ! Build surface
    call dict % init(4)
    call dict % store('id', 75)
    call dict % store('type', type)
    call dict % store('origin', [TWO, TWO, TWO])
    call dict % store('radius', TWO)
    call tst % surf % init(dict)
  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_cylinder), intent(inout) :: this

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
    class(test_cylinder), intent(inout) :: this
    real(defReal), dimension(6) :: aabb, ref
    character(nameLen)          :: name
    real(defReal), parameter    :: TOL = 1.0E-6_defReal

    ! Test ID
#line 137 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(75, this % surf % id(), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Change ID
    call this % surf % setID(1)
#line 141 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(1, this % surf % id(), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xCylinder'
      case(Y_AXIS)
        name = 'yCylinder'
      case(Z_AXIS)
        name = 'zCylinder'
    end select
#line 152 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(name, this % surf % myType(), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Bounding Box
    ref = [ZERO, ZERO, ZERO, 4.0_defReal, 4.0_defReal, 4.0_defReal]
    ref(this % axis) = -INF
    ref( 3 + this % axis) = INF
    aabb = this % surf % boundingBox()
#line 159 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
!@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_cylinder), intent(inout) :: this
    real(defReal), dimension(3) :: r, u, r_pre, u_pre

    ! Set Boundary Contidions
    ! Should ignore extra entries
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BC
    r = [TWO, TWO, TWO]
    u = ZERO
    ! Moving out at the surface in one planar direction
    r(this % plane(1)) = 4.0_defReal
    u(this % plane(1)) = ONE

    r_pre = r
    u_pre = u

    ! Explicit
    call this % surf % explicitBC(r, u)
#line 187 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 187) )
  if (anyExceptions()) return
#line 188 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
#line 188 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Transform
    call this % surf % transformBC(r, u)
#line 192 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
#line 193 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
!@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_cylinder), intent(inout) :: this
    integer(shortInt)                   :: a, p1, p2
    real(defReal), dimension(3)    :: r, u
    real(defReal)                  :: eps

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! Choose point at the surface moving in
    ! DIrection is at 45deg to the plane
    r = TWO
    r(p1) = ZERO
    u(p2) = ZERO
    u(p1) = ONE
    u(a)  = ONE
    u = u / norm2(u)

    ! At the surface
#line 222 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertFalse(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Out within SURF_TOL
    eps = -SURF_TOL
#line 226 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertFalse(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Out outside SURF_TOL
    eps = -SQRT2 * 1.00001_defReal * SURF_TOL
#line 230 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Well Outside
    eps = -TWO
#line 234 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 234) )
  if (anyExceptions()) return
#line 235 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Well withn
    eps = TWO
#line 238 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertFalse(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 238) )
  if (anyExceptions()) return
#line 239 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Tangent particle should be outside
    u(p2) =  ONE
    u(p1) = ZERO
    u(a)  = ONE
    u = u /norm2(u)
#line 245 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
!@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_cylinder), intent(inout) :: this
    integer(shortInt)                   :: a, p1, p2
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! Set axis and plane axis indices
    a = this % axis
    p1 = this % plane(1)
    p2 = this % plane(2)

    ! **Outside the cylinder
    r = TWO
    r(p1) = -ONE

    ! Perpendicular (in plane) impact
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = SQRT2
#line 275 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 275) )
  if (anyExceptions()) return
#line 276 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Oblique impact
    ! at 30deg
    u(a)  = ONE
    u(p1) = sqrt(3.0_defReal) / TWO
    u(p2) = HALF
    u = u/norm2(u)
    ref = 1.275200556_defReal * SQRT2
#line 284 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 284) )
  if (anyExceptions()) return
#line 285 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Almost Parallel
    u(a) = 1.0E+20_defReal
    u(p1) = 1.0_defReal
    u(p2) = ZERO
    u = u/norm2(u)
#line 291 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! **Exactly at the surface
    r = TWO
    r(p1) = ZERO

    ! Particle going inside
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 4.0_defReal * SQRT2
#line 303 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 303) )
  if (anyExceptions()) return
#line 304 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Particle going outside
#line 306 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 306) )
  if (anyExceptions()) return
#line 307 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! Tangent particle
    u(a) = ONE
    u(p1) = ZERO
    u(p2) = ONE
    u = u/norm2(u)
#line 313 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(INF, this % surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 313) )
  if (anyExceptions()) return
#line 314 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! **Outside within surface tolerance
    r = TWO
    r(p1) = ZERO - HALF * SURF_TOL
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = (4.0_defReal + HALF * SURF_TOL) * SQRT2
#line 323 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 323) )
  if (anyExceptions()) return
#line 324 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! ** Inside the surface
    r = TWO
    r(p1) = ONE

    ! +ve direction
    u(a) = ONE
    u(p1) = ONE
    u(p2) = ZERO
    u = u/norm2(u)
    ref = 3.0_defReal * SQRT2
#line 335 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 335) )
  if (anyExceptions()) return
#line 336 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

    ! -ve direction
    ref = 1.0_defReal * SQRT2
#line 339 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"
  call assertEqual(ref, this % surf % distance(r, -u), TOL * ref, &
 & location=SourceLocation( &
 & 'cylinder_test.f90', &
 & 339) )
  if (anyExceptions()) return
#line 340 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/cylinder_test.f90"

  end subroutine testDistance

end module cylinder_test

module Wrapcylinder_test
   use pFUnit_mod
   use cylinder_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_cylinder) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use cylinder_test
        class (test_cylinder), intent(inout) :: this
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
      aTest%test_cylinder = newTestCase(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module Wrapcylinder_test

function cylinder_test_suite() result(suite)
   use pFUnit_mod
   use cylinder_test
   use Wrapcylinder_test
   type (TestSuite) :: suite

   type (dirParam), allocatable :: testParameters(:)
   type (dirParam) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('cylinder_test_suite')

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


end function cylinder_test_suite

