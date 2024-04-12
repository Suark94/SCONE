module testResponse_test

  use numPrecision
  use testResponse_class,    only : testResponse
  use particle_class,        only : particle
  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_testResponse
    private
    type(testResponse) :: response
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testResponse


contains

  !!
  !! Sets up test_testResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testResponse), intent(inout) :: this
    type(dictionary)                        :: tempDict

    call tempDict % init(1)
    call tempDict % store('value', 1.3_defReal)
    call this % response % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_testResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testResponse), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
!@Test
  subroutine testResponseing(this)
    class(test_testResponse), intent(inout) :: this
    type(particle)                          :: p
    class(nuclearDatabase),pointer          :: xsData

#line 58 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/testResponse_test.f90"
  call assertEqual(1.3_defReal, this % response % get(p, xsData), 1.0E-9_defReal, &
 & location=SourceLocation( &
 & 'testResponse_test.f90', &
 & 58) )
  if (anyExceptions()) return
#line 59 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/testResponse_test.f90"

  end subroutine testResponseing

end module testResponse_test

module WraptestResponse_test
   use pFUnit_mod
   use testResponse_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_testResponse) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use testResponse_test
        class (test_testResponse), intent(inout) :: this
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

end module WraptestResponse_test

function testResponse_test_suite() result(suite)
   use pFUnit_mod
   use testResponse_test
   use WraptestResponse_test
   type (TestSuite) :: suite

   suite = newTestSuite('testResponse_test_suite')

   call suite%addTest(makeCustomTest('testResponseing', testResponseing))


end function testResponse_test_suite

