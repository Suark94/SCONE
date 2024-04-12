module microResponse_test

  use numPrecision
  use endfConstants
  use microResponse_class,            only : microResponse
  use particle_class,                 only : particle, P_NEUTRON
  use dictionary_class,               only : dictionary
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use materialMenu_mod,               only : init
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_microResponse
    private
    type(microResponse)        :: response_total
    type(microResponse)        :: response_eScatter
    type(microResponse)        :: response_capture
    type(microResponse)        :: response_fission
    type(microResponse)        :: response_absorbtion
    type(testNeutronDatabase)  :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_microResponse


contains

  !!
  !! Sets up test_microResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_microResponse), intent(inout) :: this
    type(dictionary)                         :: tempDict, dictMat1, dictMat2, dictMat3

    ! Allocate and initialise test nuclearData

    ! Cross-sections:         Total        eScatering   IeScatter  Capture     Fission       nuFission
    call this % xsData % build(6.0_defReal, 3.0_defReal, ZERO,     2.0_defReal, 1.0_defReal, 1.5_defReal)

    ! Set dictionaries to initialise material
    call dictMat1 % init(1)
    call dictMat2 % init(2)
    call dictMat3 % init(1)

    call dictMat3 % store('54135.03', 2.0_defReal)

    call dictMat2 % store('temp', 300.0_defReal)
    call dictMat2 % store('composition', dictMat3)

    call dictMat1 % store('Xenon', dictMat2)

    ! Initialise material
    call init(dictMat1)

    ! Set up responses
    ! Total
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_TOTAL)
    call tempDict % store('material', 'Xenon')
    call this % response_total % init(tempDict)
    call tempDict % kill()

    ! Capture
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_GAMMA)
    call tempDict % store('material', 'Xenon')
    call this % response_capture % init(tempDict)
    call tempDict % kill()

    ! Fission
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_FISSION)
    call tempDict % store('material', 'Xenon')
    call this % response_fission % init(tempDict)
    call tempDict % kill()

    ! nuFission
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_N_ELASTIC)
    call tempDict % store('material', 'Xenon')
    call this % response_eScatter % init(tempDict)
    call tempDict % kill()

    ! Absorbtion
    call tempDict % init(3)
    call tempDict % store('type', 'microResponse')
    call tempDict % store('MT', N_ABSORPTION)
    call tempDict % store('material', 'Xenon')
    call this % response_absorbtion % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_microResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_microResponse), intent(inout) :: this

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
    class(test_microResponse), intent(inout) :: this
    type(particle)                           :: p
    real(defReal), parameter :: TOL = 1.0E-9

    p % type = P_NEUTRON

    ! Test response values
#line 128 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
  call assertEqual(3.0_defReal, this % response_total % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'microResponse_test.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
#line 129 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
  call assertEqual(1.0_defReal, this % response_capture % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'microResponse_test.f90', &
 & 129) )
  if (anyExceptions()) return
#line 130 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
#line 130 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
  call assertEqual(0.5_defReal, this % response_fission % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'microResponse_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
#line 131 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
  call assertEqual(1.5_defReal, this % response_eScatter % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'microResponse_test.f90', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
#line 132 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"
  call assertEqual(1.5_defReal, this % response_absorbtion % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'microResponse_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/Tallies/TallyResponses/Tests/microResponse_test.f90"

  end subroutine testGettingResponse

end module microResponse_test

module WrapmicroResponse_test
   use pFUnit_mod
   use microResponse_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_microResponse) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use microResponse_test
        class (test_microResponse), intent(inout) :: this
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

end module WrapmicroResponse_test

function microResponse_test_suite() result(suite)
   use pFUnit_mod
   use microResponse_test
   use WrapmicroResponse_test
   type (TestSuite) :: suite

   suite = newTestSuite('microResponse_test_suite')

   call suite%addTest(makeCustomTest('testGettingResponse', testGettingResponse))


end function microResponse_test_suite

