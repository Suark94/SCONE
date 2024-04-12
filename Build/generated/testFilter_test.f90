module testFilter_test

  use numPrecision
  use testFilter_class, only : testFilter
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_testFilter
    private
    type(testFilter) :: filter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testFilter


contains

  !!
  !! Sets up test_testFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testFilter), intent(inout) :: this
    type(dictionary)                      :: tempDict

    call tempDict % init(2)
    call tempDict % store('minIdx', 4)
    call tempDict % store('maxIdx', 6)

    call this % filter % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_testFilter object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testFilter), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
!@Test
  subroutine testFiltering(this)
    class(test_testFilter), intent(inout) :: this
    type(particleState)                   :: state
    logical(defBool)                      :: filterRes

    state % matIdx = 1
    filterRes = this % filter % isPass(state)
#line 61 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'testFilter_test.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"

    state % matIdx = 4
    filterRes = this % filter % isPass(state)
#line 65 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'testFilter_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"

    state % matIdx = 6
    filterRes = this % filter % isPass(state)
#line 69 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'testFilter_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"

    state % matIdx = 7
    filterRes = this % filter % isPass(state)
#line 73 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'testFilter_test.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/testFilter_test.f90"

  end subroutine testFiltering

end module testFilter_test

module WraptestFilter_test
   use pFUnit_mod
   use testFilter_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_testFilter) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use testFilter_test
        class (test_testFilter), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod) result(aTest)
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
      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
   end function makeCustomTest

end module WraptestFilter_test

function testFilter_test_suite() result(suite)
   use pFUnit_mod
   use testFilter_test
   use WraptestFilter_test
   type (TestSuite) :: suite

   suite = newTestSuite('testFilter_test_suite')

   call suite%addTest(makeCustomTest('testFiltering', testFiltering))


end function testFilter_test_suite

