module macroResponse_test

  use numPrecision
  use endfConstants
  use macroResponse_class,            only : macroResponse
  use particle_class,                 only : particle, P_NEUTRON
  use dictionary_class,               only : dictionary
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_macroResponse
    private
    type(macroResponse)        :: response_total
    type(macroResponse)        :: response_capture
    type(macroResponse)        :: response_fission
    type(macroResponse)        :: response_nuFission
    type(macroResponse)        :: response_absorbtion
    type(testNeutronDatabase)  :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_macroResponse


contains

  !!
  !! Sets up test_macroResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_macroResponse), intent(inout) :: this
    type(dictionary)                         :: tempDict

    ! Allocate and initialise test nuclearData

    ! Cross-sections:         Total        eScatering   IeScatter  Capture     Fission       nuFission
    call this % xsData % build(6.0_defReal, 3.0_defReal, ZERO,     2.0_defReal, 1.0_defReal, 1.5_defReal)

    ! Set up responses
    ! Total
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroTotal)
    call this % response_total % init(tempDict)
    call tempDict % kill()

    ! Capture
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroCapture)
    call this % response_capture % init(tempDict)
    call tempDict % kill()

    ! Fission
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroFission)
    call this % response_fission % init(tempDict)
    call tempDict % kill()

    ! nuFission
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroNuFission)
    call this % response_nuFission % init(tempDict)
    call tempDict % kill()

    ! Absorbtion
    call tempDict % init(2)
    call tempDict % store('type','macroResponse')
    call tempDict % store('MT', macroAbsorbtion)
    call this % response_absorbtion % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_macroResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_macroResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
!@Test
  subroutine testGettingResponse(this)
    class(test_macroResponse), intent(inout) :: this
    type(particle)                           :: p
    real(defReal), parameter :: TOL = 1.0E-9

    p % type = P_NEUTRON

    ! Test response values
#line 107 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
  call assertEqual(6.0_defReal, this % response_total % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'macroResponse_test.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
#line 108 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
  call assertEqual(2.0_defReal, this % response_capture % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'macroResponse_test.f90', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
#line 109 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
  call assertEqual(1.0_defReal, this % response_fission % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'macroResponse_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
#line 110 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
  call assertEqual(1.5_defReal, this % response_nuFission % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'macroResponse_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
#line 111 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"
  call assertEqual(3.0_defReal, this % response_absorbtion % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'macroResponse_test.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/macroResponse_test.f90"

  end subroutine testGettingResponse

end module macroResponse_test

module WrapmacroResponse_test
   use pFUnit_mod
   use macroResponse_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_macroResponse) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use macroResponse_test
        class (test_macroResponse), intent(inout) :: this
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

end module WrapmacroResponse_test

function macroResponse_test_suite() result(suite)
   use pFUnit_mod
   use macroResponse_test
   use WrapmacroResponse_test
   type (TestSuite) :: suite

   suite = newTestSuite('macroResponse_test_suite')

   call suite%addTest(makeCustomTest('testGettingResponse', testGettingResponse))


end function macroResponse_test_suite

