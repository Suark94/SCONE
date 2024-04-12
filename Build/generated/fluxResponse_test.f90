module fluxResponse_test

  use numPrecision
  use fluxResponse_class,    only : fluxResponse
  use particle_class,        only : particle
  use dictionary_class,      only : dictionary
  use nuclearDatabase_inter, only : nuclearDatabase
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_fluxResponse
    private
    type(fluxResponse) :: response
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_fluxResponse


contains

  !!
  !! Sets up test_fluxResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_fluxResponse), intent(inout) :: this
    type(dictionary)                      :: tempDict

  end subroutine setUp

  !!
  !! Kills test_fluxResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_fluxResponse), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
!@Test
  subroutine fluxResponseing(this)
    class(test_fluxResponse), intent(inout) :: this
    type(particle)                          :: p
    class(nuclearDatabase),pointer          :: xsData

#line 54 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/fluxResponse_test.f90"
  call assertEqual(ONE, this % response % get(p, xsData), 1.0E-9_defReal, &
 & location=SourceLocation( &
 & 'fluxResponse_test.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/fluxResponse_test.f90"

  end subroutine fluxResponseing

end module fluxResponse_test

module WrapfluxResponse_test
   use pFUnit_mod
   use fluxResponse_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_fluxResponse) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use fluxResponse_test
        class (test_fluxResponse), intent(inout) :: this
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

end module WrapfluxResponse_test

function fluxResponse_test_suite() result(suite)
   use pFUnit_mod
   use fluxResponse_test
   use WrapfluxResponse_test
   type (TestSuite) :: suite

   suite = newTestSuite('fluxResponse_test_suite')

   call suite%addTest(makeCustomTest('fluxResponseing', fluxResponseing))


end function fluxResponse_test_suite

