module testMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use testMap_class,           only : testMap

  implicit none


!@testCase
  type, extends(TestCase) :: test_testMap
    private
    type(testMap) :: map

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_testMap

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_testMap), intent(inout) :: this
    type(dictionary) :: tempDict

    call tempDict % init(1)
    call tempDict % store('maxIdx',4)
    call this % map % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_testMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that maps performs as expected
  !!
!@Test
  subroutine testMapping(this)
    class(test_testMap), intent(inout) :: this
    type(particleState) :: state

    state % matIdx = 0
#line 59 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(0,this % map % map(state),'Invalid idx case:', &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 59) )
  if (anyExceptions()) return
#line 60 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

    state % matIdx = 1
#line 62 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(1,this % map % map(state),'Normal idx case:', &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

    state % matIdx = 2
#line 65 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(2,this % map % map(state),'Normal idx case:', &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

    state % matIdx = 4
#line 68 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(4,this % map % map(state),'Normal idx case:', &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

    state % matIdx = 5
#line 71 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(0,this % map % map(state),'Invalid idx case:', &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

    state % matIdx = -1
#line 74 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(0,this % map % map(state),'Invalid idx case:', &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

  end subroutine testMapping

  !!
  !! Test getting number of bins and dimensions
  !!
!@Test
  subroutine testNumBin(this)
    class(test_testMap), intent(inout) :: this

    ! Test number of bins
#line 86 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(4, this % map % bins(0), &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
#line 87 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(4, this % map % bins(1), &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
#line 88 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(0, this % map % bins(2), &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
#line 89 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(0, this % map % bins(-2), &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

    ! Test dimension
#line 92 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"
  call assertEqual(1, this % map % dimensions(), &
 & location=SourceLocation( &
 & 'testMap_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/testMap_test.f90"

  end subroutine testNumBin

end module testMap_test

module WraptestMap_test
   use pFUnit_mod
   use testMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_testMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use testMap_test
        class (test_testMap), intent(inout) :: this
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

end module WraptestMap_test

function testMap_test_suite() result(suite)
   use pFUnit_mod
   use testMap_test
   use WraptestMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('testMap_test_suite')

   call suite%addTest(makeCustomTest('testMapping', testMapping))

   call suite%addTest(makeCustomTest('testNumBin', testNumBin))


end function testMap_test_suite

