module weightWindowsField_iTest
  use numPrecision
  use pFUnit_mod
  use particle_class,           only : particle
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use weightWindowsField_class, only : weightWindowsField

  implicit none


!@testCase
  type, extends(TestCase) :: test_weightWindows
    private
    type(weightWindowsField) :: wwField
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_weightWindows


  !!
  !! Weight Windows Definition
  !!
  character(*), parameter :: DICT_DEF = &
  & "file ./IntegrationTestFiles/testWW ;"

contains

  !!
  !! Sets up test_weightWindows object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_weightWindows), intent(inout) :: this
    type(dictionary)                         :: dict

    call charToDict(dict, DICT_DEF)

    call this % wwField % init(dict)

  end subroutine setUp

  !!
  !! Kills test_weightWindows object
  !!
  subroutine tearDown(this)
    class(test_weightWindows), intent(inout) :: this

    call this % wwField % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving the weight window values
  !!
!@Test
  subroutine testGetValue(this)
    class(test_weightWindows), intent(inout) :: this
    type(particle)                           :: p
    integer(shortInt), dimension(3)          :: bins, EXPECTED_BINS

    p % isMG = .false.
    p % coords % lvl(1) % r = [0.5, 7.0, 0.0]
    p % E = 10.0

    bins = this % wwField % at(p)
    EXPECTED_BINS = [0.4, 1.5, 0.8]

#line 73 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/weightWindowsField_iTest.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'weightWindowsField_iTest.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/weightWindowsField_iTest.f90"

    p % isMG = .false.
    p % coords % lvl(1) % r = [-0.5, 7.0, 0.0]
    p % E = 10.0

    bins = this % wwField % at(p)
    EXPECTED_BINS = ZERO

#line 82 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/weightWindowsField_iTest.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'weightWindowsField_iTest.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/weightWindowsField_iTest.f90"

  end subroutine testGetValue


end module weightWindowsField_iTest

module WrapweightWindowsField_iTest
   use pFUnit_mod
   use weightWindowsField_iTest
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_weightWindows) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use weightWindowsField_iTest
        class (test_weightWindows), intent(inout) :: this
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

end module WrapweightWindowsField_iTest

function weightWindowsField_iTest_suite() result(suite)
   use pFUnit_mod
   use weightWindowsField_iTest
   use WrapweightWindowsField_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('weightWindowsField_iTest_suite')

   call suite%addTest(makeCustomTest('testGetValue', testGetValue))


end function weightWindowsField_iTest_suite

