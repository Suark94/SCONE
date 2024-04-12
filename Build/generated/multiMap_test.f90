module multiMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use dictParser_func,         only : charToDict
  use multiMap_class,          only : multiMap
  use outputFile_class,        only : outputFile
  implicit none


!@testCase
  type, extends(TestCase) :: test_multiMap
    type(multiMap) :: map

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_multiMap

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_multiMap), intent(inout) :: this
    type(dictionary)                    :: tempDict
    character(*),parameter :: def = "                                       &
      type multiMap;                                                        &
      maps (map1 map2 map3);                                                &
      map1 {type spaceMap; axis x; grid unstruct; bins (0.0 1.0 2.0); }     &
      map2 {type spaceMap; axis y; grid unstruct; bins (0.0 2.0 4.0 6.0); } &
      map3 {type spaceMap; axis z; grid unstruct; bins (0.0 3.0 6.0); }     "


    call charToDict(tempDict, def)
    call this % map % init(tempDict)

    call tempDict % kill()
  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_multiMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that map correctly returns number of bins and dimensions
  !!
!@Test
  subroutine testBinAndDimension(this)
    class(test_multiMap), intent(inout) :: this

    ! Test approperiate dimension
#line 65 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(3, this % map % dimensions(), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

    ! Test Bin number
#line 68 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(12, this % map % bins(0), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
#line 69 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(2, this % map % bins(1), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
#line 70 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(3, this % map % bins(2), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
#line 71 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(2, this % map % bins(3), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

#line 73 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(0, this % map % bins(-1), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
#line 74 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(0, this % map % bins(4), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

  end subroutine testBinAndDimension

  !!
  !! Test mapping
  !!
!@Test
  subroutine testMapping(this)
    class(test_multiMap), intent(inout) :: this
    type(particleState)                 :: state

    ! Map to bin (1 1 1) -> idx == 1
    state % r = [0.1_defReal, 0.1_defReal, 0.1_defReal]
#line 88 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(1, this % map % map(state), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

    ! Map to bin (2 3 2) -> idx == 12
    state % r = [1.1_defReal, 5.1_defReal, 5.1_defReal]
#line 92 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(12, this % map % map(state), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

    ! Map to bin (2, 2, 1) -> idx == 4
    state % r = [1.1_defReal, 3.1_defReal, 2.1_defReal]
#line 96 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(4, this % map % map(state), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

    ! Map outside division -> idx == 0
    state % r = [-1.1_defReal, 5.1_defReal, 5.1_defReal]
#line 100 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(0, this % map % map(state), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 100) )
  if (anyExceptions()) return
#line 101 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

    state % r = [1.1_defReal, 50.1_defReal, 5.1_defReal]
#line 103 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(0, this % map % map(state), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 103) )
  if (anyExceptions()) return
#line 104 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

    state % r = [1.1_defReal, 5.1_defReal, -5.1_defReal]
#line 106 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertEqual(0, this % map % map(state), &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"

  end subroutine testMapping




  !!
  !! Test correctness of calls when printing
  !!
!@Test
  subroutine testPrint(this)
    class(test_multiMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map % print(out)
#line 124 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
  call assertTrue(out % isValid(),'Incorrect printing sequence: ', &
 & location=SourceLocation( &
 & 'multiMap_test.f90', &
 & 124) )
  if (anyExceptions()) return
#line 125 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/multiMap_test.f90"
    call out % reset() 

  end subroutine testPrint


end module multiMap_test

module WrapmultiMap_test
   use pFUnit_mod
   use multiMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_multiMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use multiMap_test
        class (test_multiMap), intent(inout) :: this
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

end module WrapmultiMap_test

function multiMap_test_suite() result(suite)
   use pFUnit_mod
   use multiMap_test
   use WrapmultiMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('multiMap_test_suite')

   call suite%addTest(makeCustomTest('testBinAndDimension', testBinAndDimension))

   call suite%addTest(makeCustomTest('testMapping', testMapping))

   call suite%addTest(makeCustomTest('testPrint', testPrint))


end function multiMap_test_suite

