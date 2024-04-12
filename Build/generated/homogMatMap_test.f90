module homogMatMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use dictParser_func,         only : charToDict
  use outputFile_class,        only : outputFile

  ! May not be ideal but there is a dependance on Global materialMenu
  use materialMenu_mod,        only : mm_init => init, mm_kill => kill

  use homogMatMap_class,       only : homogMatMap

  implicit none


!@testCase
  type, extends(TestCase) :: test_homogMatMap
    private
    type(homogMatMap),allocatable :: map_noUndef
    type(homogMatMap),allocatable :: map_Undef
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_homogMatMap


  !!
  !! Test parameters
  !!
  character(*), dimension(*), parameter :: BIN_NAMES   = ['bin1','bin2','bin3']
  character(*), dimension(*), parameter :: MAT_IN_BIN1 = ['mat3','mat4']
  character(*), dimension(*), parameter :: MAT_IN_BIN2 = ['mat1']
  character(*), dimension(*), parameter :: MAT_IN_BIN3 = ['mat5']

  !!
  !! Material Definitions
  !!
  character(*), parameter :: DICT_DEF = &
  " mat1 { temp 17; composition {} } &
  & mat2 { temp 17; composition {} } &
  & mat3 { temp 17; composition {} } &
  & mat4 { temp 17; composition {} } &
  & mat5 { temp 17; composition {} } "



contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_homogMatMap), intent(inout) :: this
    type(dictionary)                  :: dict
    type(dictionary)                  :: mapDict1

    ! Build nuclear data
    call charToDict(dict, DICT_DEF)
    call mm_init(dict)

    ! Initialise dictionaries
    call mapDict1 % init(5)

    ! Build material map definition
    call mapDict1 % store('bins', BIN_NAMES)
    call mapDict1 % store(BIN_NAMES(1), MAT_IN_BIN1)
    call mapDict1 % store(BIN_NAMES(2), MAT_IN_BIN2)
    call mapDict1 % store(BIN_NAMES(3), MAT_IN_BIN3)
    allocate(this % map_noUndef, source = homogMatMap(mapDict1))

    call mapDict1 % store('undefBin','true')
    allocate(this % map_undef, source = homogMatMap(mapDict1))

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_homogMatMap), intent(inout) :: this

    call mm_kill()
    call this % map_noUndef % kill()
    call this % map_Undef % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Mapping test without undefined bin
  !!
!@Test
  subroutine testMappingNoUndefined(this)
    class(test_homogMatMap), intent(inout)   :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(6)           :: bins
    integer(shortInt),dimension(6),parameter :: EXPECTED_BINS = [2, 0, 1, 1, 3, 0]

    do i = 1,6
      state % matIdx = i
      bins(i) = this % map_noUndef % map(state)
    end do

#line 109 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"

  end subroutine testMappingNoUndefined


  !!
  !! Mapping test with undefined bin
  !!
!@Test
  subroutine testMappingUndefined(this)
    class(test_homogMatMap), intent(inout)   :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(6)           :: bins
    integer(shortInt),dimension(6),parameter :: EXPECTED_BINS = [2, 4, 1, 1, 3, 4]

    do i = 1,6
      state % matIdx = i
      bins(i) = this % map_undef % map(state)
    end do

#line 130 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"

  end subroutine testMappingUndefined


  !!
  !! Test number of bins inquiry
  !!
!@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_homogMatMap), intent(inout) :: this

#line 142 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(3, this % map_noUndef % bins(1), 'homogMatMap without undefined bin', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
#line 143 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(4, this % map_undef % bins(1),   'homogMatMap with undefined bin', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 143) )
  if (anyExceptions()) return
#line 144 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
#line 144 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(3, this % map_noUndef % bins(0), 'Number of all bins', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 144) )
  if (anyExceptions()) return
#line 145 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
#line 145 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(0, this % map_noUndef % bins(2), 'higher dimension', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
#line 146 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertEqual(0, this % map_noUndef % bins(-2),'invalid dimension', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"

  end subroutine testNumberOfBinsInquiry

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
!@Test
  subroutine testPrint(this)
    class(test_homogMatMap), intent(inout) :: this
    type(outputFile)                       :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_noUndef % print(out)
#line 162 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertTrue(out % isValid(),'For map with no undefined material bin: ', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
    call out % reset()

    call this % map_undef % print(out)
#line 166 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
  call assertTrue(out % isValid(),'For map with undefined material bin: ', &
 & location=SourceLocation( &
 & 'homogMatMap_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/homogMatMap_test.f90"
    call out % reset()

  end subroutine testPrint



end module homogMatMap_test

module WraphomogMatMap_test
   use pFUnit_mod
   use homogMatMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_homogMatMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use homogMatMap_test
        class (test_homogMatMap), intent(inout) :: this
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

end module WraphomogMatMap_test

function homogMatMap_test_suite() result(suite)
   use pFUnit_mod
   use homogMatMap_test
   use WraphomogMatMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('homogMatMap_test_suite')

   call suite%addTest(makeCustomTest('testMappingNoUndefined', testMappingNoUndefined))

   call suite%addTest(makeCustomTest('testMappingUndefined', testMappingUndefined))

   call suite%addTest(makeCustomTest('testNumberOfBinsInquiry', testNumberOfBinsInquiry))

   call suite%addTest(makeCustomTest('testPrint', testPrint))


end function homogMatMap_test_suite

