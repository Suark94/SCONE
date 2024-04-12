module aPlane_test
  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use aPlane_class,      only : aPlane
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
    type, extends(ParameterizedTestCase) :: test_aPlane
      integer(shortInt)             :: axis
      type(aPlane)                  :: surf
    contains
      procedure :: tearDown
    end type test_aPlane

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
         string = 'xPlane'
       case(Y_AXIS)
         string = 'yPlane'
       case(Z_AXIS)
         string = 'zPlane'
       case default
         string ="Unknown"
      end select
  end function toString

  !!
  !! Build new test_aPlane test case
  !! Given integer direction X_AXIS, Y_AXIS or Z_AXIS
  !!
  !! ID 75
  !! x0/y0/z0 4.3;
  !!
  function newTestCase(dir) result(tst)
    type(dirParam), intent(in) :: dir
    type(test_aPlane)        :: tst
    type(dictionary)      :: dict
    character(nameLen)    :: type

    ! Start dictionary
    ! Build surface
    call dict % init(4)
    call dict % store('id', 75)

    ! Select type of cylinder and axis
    select case(dir % dir)
      case(X_AXIS)
        tst % axis = X_AXIS
        call dict % store('type','xPlane')
        call dict % store('x0',4.3_defReal)

      case(Y_AXIS)
        tst % axis = Y_AXIS
        call dict % store('type','yPlane')
        call dict % store('y0',4.3_defReal)

      case(Z_AXIS)
        tst % axis = Z_AXIS
        call dict % store('type','zPlane')
        call dict % store('z0',4.3_defReal)

      case default
        print *, "Should not happen. Wrong direction in testcase constructor"

    end select

    ! Build surface
    call tst % surf % init(dict)

  end function newTestCase

  !!
  !! Deconstruct the test case
  !!
  subroutine tearDown(this)
    class(test_aPlane), intent(inout) :: this

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
    class(test_aPlane), intent(inout) :: this
    real(defReal), dimension(6) :: aabb, ref
    character(nameLen)          :: name
    real(defReal), parameter    :: TOL = 1.0E-6_defReal

     ! Test ID
#line 136 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(75, this % surf % id(), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Change ID
    call this % surf % setID(1)
#line 140 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(1, this % surf % id(), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Name
    select case(this % axis)
      case(X_AXIS)
        name = 'xPlane'
      case(Y_AXIS)
        name = 'yPlane'
      case(Z_AXIS)
        name = 'zPlane'
    end select
#line 151 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(name, this % surf % myType(), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Bounding Box
    ref = [-INF, -INF, -INF, INF, INF, INF]
    ref(this % axis) = 4.3_defReal
    ref( 3 + this % axis) = 4.3_defReal
    aabb = this % surf % boundingBox()
#line 158 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
!@Test(cases=[1,2,3])
  subroutine testBC(this)
    class(test_aPlane), intent(inout) :: this
    real(defReal), dimension(3) :: r, u, r_pre, u_pre

    ! Set Boundary Contidions
    ! Should ignore extra entries
    call this % surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BC
    r = [TWO, TWO, TWO]
    u = ZERO

    ! Put atthe surface
    r(this % axis) = 4.3_defReal
    u(this % axis) = ONE

    r_pre = r
    u_pre = u

    ! Explicit
    call this % surf % explicitBC(r, u)
#line 187 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 187) )
  if (anyExceptions()) return
#line 188 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
#line 188 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Transform
    call this % surf % transformBC(r, u)
#line 192 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
#line 193 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

  end subroutine testBC

  !!
  !! Test Halfspaces membership
  !!
!@Test(cases=[1,2,3])
  subroutine testHalfspace(this)
    class(test_aPlane), intent(inout) :: this
    integer(shortInt)                 :: a, p1, p2
    real(defReal), dimension(3)       :: r, u, u2
    real(defReal)                     :: eps

    r = ZERO
    u = ZERO

    ! At the surface
    r(this % axis) = 4.3_defReal
    u(this % axis) = ONE
#line 213 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertTrue(this % surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 213) )
  if (anyExceptions()) return
#line 214 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! In within SURF_TOL
    eps = -HALF * SURF_TOL
#line 217 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! In outside SURF_TOL
    eps = -1.0001_defReal * SURF_TOL
#line 221 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertFalse(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Well Outside
    eps = TWO
#line 225 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 225) )
  if (anyExceptions()) return
#line 226 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Well withn
    eps = -TWO
#line 229 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertFalse(this % surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 229) )
  if (anyExceptions()) return
#line 230 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! Tangent particle should use position
    eps = HALF * SURF_TOL
    u2 = cshift(u, 1) ! Point to an orthogonal direction

#line 235 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertTrue(this % surf % halfspace(r + eps*u, u2), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
#line 236 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertFalse(this % surf % halfspace(r - eps*u, u2), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 236) )
  if (anyExceptions()) return
#line 237 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

  end subroutine testHalfspace

  !!
  !! Test distance calculation
  !!
!@Test(cases=[1, 2, 3])
  subroutine testDistance(this)
    class(test_aPlane), intent(inout) :: this
    integer(shortInt)                   :: a, p1, p2
    real(defReal), dimension(3)         :: r, u, u2
    real(defReal)                       :: ref
    real(defReal), parameter :: SQRT3 = sqrt(3.0_defReal)
    real(defReal), parameter :: TOL = 1.0E-7

    ! ** Inside
    r = TWO
    u = ZERO
    r(this % axis) = ZERO
    u = ONE   ! 45deg to each axis
    u = u/norm2(u)

    ref = 4.3_defReal * SQRT3
#line 260 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(ref, this % surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 260) )
  if (anyExceptions()) return
#line 261 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! ** Exactly at surface
    r(this % axis) = 4.3_defReal
#line 264 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 264) )
  if (anyExceptions()) return
#line 265 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! ** Within Surface tolerance
    r(this % axis) = 4.3_defReal - HALF * SURF_TOL
#line 268 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 268) )
  if (anyExceptions()) return
#line 269 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! ** Well outside
    ! +ve direction
    r(this % axis) = 5.0_defReal
#line 273 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(INF, this % surf % distance(r, u), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 273) )
  if (anyExceptions()) return
#line 274 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! -ve direction
    ref = 0.7_defReal * SQRT3
#line 277 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(ref, this % surf % distance(r, -u), TOL * ref, &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 277) )
  if (anyExceptions()) return
#line 278 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

    ! ** Parallel to plane
    r = TWO
    r(this % axis) = ZERO
    u2 = ZERO
    u2(this % axis) = ONE
    u2 = cshift(u2, 1) ! Point to an orthogonal direction
#line 285 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"
  call assertEqual(INF, this % surf % distance(r, u2), &
 & location=SourceLocation( &
 & 'aPlane_test.f90', &
 & 285) )
  if (anyExceptions()) return
#line 286 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/aPlane_test.f90"

  end subroutine testDistance

end module aPlane_test

module WrapaPlane_test
   use pFUnit_mod
   use aPlane_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_aPlane) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use aPlane_test
        class (test_aPlane), intent(inout) :: this
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
      aTest%test_aPlane = newTestCase(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WrapaPlane_test

function aPlane_test_suite() result(suite)
   use pFUnit_mod
   use aPlane_test
   use WrapaPlane_test
   type (TestSuite) :: suite

   type (dirParam), allocatable :: testParameters(:)
   type (dirParam) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('aPlane_test_suite')

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


end function aPlane_test_suite

