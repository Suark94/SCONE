module materialMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use dictParser_func,         only : charToDict
  use outputFile_class,        only : outputFile

  ! May not be ideal but there is a dependance on Global materialMenu
  use materialMenu_mod,        only : mm_init => init, mm_kill => kill

  use materialMap_class,       only : materialMap

  implicit none


!@testCase
  type, extends(TestCase) :: test_materialMap
    private
    type(materialMap),allocatable :: map_noUndef
    type(materialMap),allocatable :: map_Undef
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_materialMap


  !!
  !! Test parameters
  !!
  character(*),dimension(*),parameter :: MAT_NAMES=['mat1','mat2','mat3','mat4','mat5']
  character(*),dimension(*),parameter :: MAT_IN_MAP =['mat2','mat3','mat5']

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
    class(test_materialMap), intent(inout) :: this
    type(dictionary)                  :: dict
    type(dictionary)                  :: mapDict1



!    !*** Provisional -> allow registry without handles
!    call tempDict2 % store('myMat','datalessMaterials')
!
!    ! Store empty- handles dictionary
!    call dict % store('handles', tempDict2)
!
!    ! Create materials dictionary of empty dictionaries
!    do i=1,size(MAT_NAMES)
!      call tempDict1 % store(MAT_NAMES(i), tempDict3)
!
!    end do
!
!    ! Store materials dictionary in nuclearData dict
!    call dict % store('materials',tempDict1)


    ! Build nuclear data
    call charToDict(dict, DICT_DEF)
    call mm_init(dict)

    ! Initialise dictionaries
    call mapDict1 % init(2)

    ! Build material map definition
    call mapDict1 % store('materials', MAT_IN_MAP)
    allocate(this % map_noUndef, source = materialMap(mapDict1))

    call mapDict1 % store('undefBin','true')
    allocate(this % map_undef, source = materialMap(mapDict1))

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_materialMap), intent(inout) :: this

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
    class(test_materialMap), intent(inout)   :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [0, 1, 2, 0, 3]

    do i=1,5
      state % matIdx = i
      bins(i) = this % map_noUndef % map(state)
    end do

#line 122 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"

  end subroutine testMappingNoUndefined


  !!
  !! Mapping test with undefined bin
  !!
!@Test
  subroutine testMappingUndefined(this)
    class(test_materialMap), intent(inout)   :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [4, 1, 2, 4, 3]

    do i=1,5
      state % matIdx = i
      bins(i) = this % map_undef % map(state)
    end do

#line 143 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 143) )
  if (anyExceptions()) return
#line 144 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"

  end subroutine testMappingUndefined


  !!
  !! Test number of bins inquiry
  !!
!@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_materialMap), intent(inout) :: this

#line 155 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(3, this % map_noUndef % bins(1),'materialMap without undefined bin', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
#line 156 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(4, this % map_undef % bins(1),'materialMap with undefined bin', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
#line 157 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(3, this % map_noUndef % bins(0), 'Number of all bins', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
#line 158 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(0, this % map_noUndef % bins(2), 'higher dimension', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
#line 159 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertEqual(0, this % map_noUndef % bins(-2),'invalid dimension', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"

  end subroutine testNumberOfBinsInquiry

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
!@Test
  subroutine testPrint(this)
    class(test_materialMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_noUndef % print(out)
#line 175 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertTrue(out % isValid(),'For map with no undefined material bin: ', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
    call out % reset()

    call this % map_undef % print(out)
#line 179 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
  call assertTrue(out % isValid(),'For map with undefined material bin: ', &
 & location=SourceLocation( &
 & 'materialMap_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/materialMap_test.f90"
    call out % reset()

  end subroutine testPrint



end module materialMap_test

module WrapmaterialMap_test
   use pFUnit_mod
   use materialMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_materialMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use materialMap_test
        class (test_materialMap), intent(inout) :: this
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

end module WrapmaterialMap_test

function materialMap_test_suite() result(suite)
   use pFUnit_mod
   use materialMap_test
   use WrapmaterialMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('materialMap_test_suite')

   call suite%addTest(makeCustomTest('testMappingNoUndefined', testMappingNoUndefined))

   call suite%addTest(makeCustomTest('testMappingUndefined', testMappingUndefined))

   call suite%addTest(makeCustomTest('testNumberOfBinsInquiry', testNumberOfBinsInquiry))

   call suite%addTest(makeCustomTest('testPrint', testPrint))


end function materialMap_test_suite

