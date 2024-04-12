module uniFissSitesField_test
  use numPrecision
  use pFUnit_mod
  use particle_class,           only : particle, particleState
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use geometry_inter,           only : geometry
  use RNG_class,                only : RNG
  use uniFissSitesField_class,  only : uniFissSitesField

  implicit none


!@testCase
  type, extends(TestCase) :: test_uniFissSitesField
    private
    type(uniFissSitesField) :: ufsField
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_uniFissSitesField

  !!
  !! Map definition
  !!
  character(*), parameter :: DICT_DEF = &
  " type spaceMap;  axis z;  grid unstruct; &
    bins (0.0 20.0 40.0 60.0 80.0); "

contains

  !!
  !! Sets up test_weightWindows object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_uniFissSitesField), intent(inout) :: this
    type(dictionary)                             :: dict, dictMap
    class(geometry), pointer                     :: geom
    class(RNG), pointer                          :: rand
    integer(shortInt)                            :: type

    call charToDict(dictMap, DICT_DEF)

    ! Initialise dictionaries
    call dict % init(3)

    ! Build material map definition
    call dict % store('type', 'uniFissSitesField')
    call dict % store('uniformMap', 1)
    call dict % store('map', dictMap)

    call this % ufsField % init(dict)
    call this % ufsField % estimateVol(geom, rand, type)

  end subroutine setUp

  !!
  !! Kills test_weightWindows object
  !!
  subroutine tearDown(this)
    class(test_uniFissSitesField), intent(inout) :: this

    call this % ufsField % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving the ufs values
  !!
!@Test
  subroutine testGetValue(this)
    class(test_uniFissSitesField), intent(inout) :: this
    type(particle)                               :: p
    type(particleState)                          :: state
    integer(shortInt), dimension(3)              :: bins, EXPECTED_BINS

    ! Test case in the map
    p % coords % lvl(1) % r = [0.5, 7.0, 50.0]

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [0.25, 0.25, 0.0]
#line 86 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniFissSitesField_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'uniFissSitesField_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniFissSitesField_test.f90"

    ! Test case outside the map
    p % coords % lvl(1) % r = [0.5, 7.0, 100.0]

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [1.0, 1.0, 1.0]
#line 93 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniFissSitesField_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'uniFissSitesField_test.f90', &
 & 93) )
  if (anyExceptions()) return
#line 94 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniFissSitesField_test.f90"

    ! Modify the map by storing fission sites
    state % r   = [0.5, 7.0, 12.0]
    state % wgt = 0.2

    call this % ufsField % storeFS(state)

    state % r   = [0.5, 7.0, 23.2]
    state % wgt = 0.8
    call this % ufsField % storeFS(state)

    call this % ufsField % updateMap()

    ! Test case in the updated map
    p % coords % lvl(1) % r = [0.5, 7.0, 18.1]

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [0.25, 0.06666666667, 0.0]
#line 112 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniFissSitesField_test.f90"
  call assertEqual(EXPECTED_BINS,bins, &
 & location=SourceLocation( &
 & 'uniFissSitesField_test.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniFissSitesField_test.f90"

  end subroutine testGetValue


end module uniFissSitesField_test

module WrapuniFissSitesField_test
   use pFUnit_mod
   use uniFissSitesField_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_uniFissSitesField) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use uniFissSitesField_test
        class (test_uniFissSitesField), intent(inout) :: this
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

end module WrapuniFissSitesField_test

function uniFissSitesField_test_suite() result(suite)
   use pFUnit_mod
   use uniFissSitesField_test
   use WrapuniFissSitesField_test
   type (TestSuite) :: suite

   suite = newTestSuite('uniFissSitesField_test_suite')

   call suite%addTest(makeCustomTest('testGetValue', testGetValue))


end function uniFissSitesField_test_suite

