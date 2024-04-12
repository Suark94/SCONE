module spaceMap_test
  use numPrecision
  use universalVariables
  use pFUnit_mod
  use particle_class,       only : particleState
  use dictionary_class,     only : dictionary
  use outputFile_class,     only : outputFile
  use spaceMap_class,       only : spaceMap

  implicit none


!@testParameter
  type, extends(AbstractTestParameter) :: dirPar
    integer(shortInt) :: dir
  contains
    procedure :: toString
  end type dirPar


!@testCase(constructor=newTest)
  type, extends(ParameterizedTestCase) :: test_spaceMap
    private
    integer(shortInt) :: dir
    type(spaceMap) :: map_struct
    type(spaceMap) :: map_unstruct

  end type test_spaceMap

contains

  !!
  !! Returns an array of test parameters
  !!
  function getParameters() result(params)
    type(dirPar), dimension(3) :: params

    params(1) % dir = X_AXIS
    params(2) % dir = Y_AXIS
    params(3) % dir = Z_AXIS

  end function getParameters

  !!
  !! Returns only direction x
  !!
  function XdirParameter() result(params)
    type(dirPar), dimension(1) :: params

    params(1) % dir = X_AXIS

  end function XdirParameter

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(dirPar), intent(in) :: this
    character(:), allocatable :: string
    character(1)              :: str

    select case(this % dir)
      case(X_AXIS)
        str ='x'
      case(Y_AXIS)
        str='y'
      case(Z_AXIS)
        str='z'
      case default
        str='?'
    end select
    string = str

  end function toString


  !!
  !! Construct test case
  !!
  function newTest(testParam) result(tst)
    type(dirPar), intent(in) :: testParam
    type(test_spaceMap)      :: tst
    type(dictionary)         :: tempDict
    real(defReal),dimension(*),parameter :: BIN_DIV = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, &
                                                        2.0,   4.0,  6.0,  8.0,  10.0]

    ! Load direction
    tst % dir = testParam % dir

    ! Create structured grid
    call tempDict % init(5)
    call tempDict % store('axis',testParam % toString())
    call tempDict % store('grid','lin')
    call tempDict % store('min', -10.0_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 20)

    call tst % map_struct % init(tempDict)
    call tempDict % kill()

    ! Create unstructured grid
    call tempDict % init(3)
    call tempDict % store('axis',testParam % toString())
    call tempDict % store('grid','unstruct')
    call tempDict % store('bins', BIN_DIV)

    call tst % map_unstruct % init(tempDict)
    call tempDict % kill()

  end function newTest

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test structured grid
  !!
!@Test(testParameters={getParameters()})
  subroutine testStructuredGrid(this)
    class(test_spaceMap), intent(inout)      :: this
    real(defReal),dimension(2),parameter     :: POS = [0.5_defReal, -10.1_defReal]
    integer(shortInt),dimension(2),parameter :: RES = [11, 0]
    integer(shortInt),dimension(2)           :: idx
    type(particleState),dimension(2)         :: states

    states % r(this % dir) = POS
    idx = this % map_struct % map(states)

#line 130 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(RES, idx, &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"

  end subroutine testStructuredGrid

  !!
  !! Test unstructured grid
  !!
!@Test(testParameters={getParameters()})
  subroutine testUnstructuredGrid(this)
    class(test_spaceMap), intent(inout)      :: this
    real(defReal),dimension(2),parameter     :: POS = [0.5_defReal, -10.1_defReal]
    integer(shortInt),dimension(2),parameter :: RES = [6, 0]
    integer(shortInt),dimension(2)           :: idx
    type(particleState),dimension(2)         :: states

    states % r(this % dir) = POS
    idx = this % map_unstruct % map(states)

#line 148 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(RES, idx, &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"

  end subroutine testUnstructuredGrid

 !!
 !! Test bin output
 !!
!@Test(testParameters ={XdirParameter()})
  subroutine testBins(this)
    class(test_spaceMap), intent(inout) :: this

    ! Structured grid
#line 160 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(20, this % map_struct % bins(1),'Normal use', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
#line 161 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(20, this % map_struct % bins(0),'All binbs', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
#line 162 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(0, this % map_struct % bins(-2),'Invalid dimension', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"

    ! Unstructured grid
#line 165 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(10, this % map_unstruct % bins(1),'Normal use', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
#line 166 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(10, this % map_unstruct % bins(0),'All bins', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
#line 167 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertEqual(0, this % map_unstruct % bins(-2),'Invalid dimension', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 167) )
  if (anyExceptions()) return
#line 168 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"

  end subroutine testBins

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
!@Test(testParameters ={XdirParameter()})
  subroutine testPrint(this)
    class(test_spaceMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_struct % print(out)
#line 183 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertTrue(out % isValid(),'For map with structured grid: ', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
    call out % reset()

    call this % map_unstruct % print(out)
#line 187 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
  call assertTrue(out % isValid(),'For map with unstructured grid: ', &
 & location=SourceLocation( &
 & 'spaceMap_test.f90', &
 & 187) )
  if (anyExceptions()) return
#line 188 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/spaceMap_test.f90"
    call out % reset() 

  end subroutine testPrint


end module spaceMap_test

module WrapspaceMap_test
   use pFUnit_mod
   use spaceMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_spaceMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use spaceMap_test
        class (test_spaceMap), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod, testParameter) result(aTest)
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
      type (dirPar), intent(in) :: testParameter
      aTest%test_spaceMap = newTest(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WrapspaceMap_test

function spaceMap_test_suite() result(suite)
   use pFUnit_mod
   use spaceMap_test
   use WrapspaceMap_test
   type (TestSuite) :: suite

   type (dirPar), allocatable :: testParameters(:)
   type (dirPar) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('spaceMap_test_suite')

   testParameters = getParameters()

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testStructuredGrid', testStructuredGrid, testParameter))
   end do

   testParameters = getParameters()

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testUnstructuredGrid', testUnstructuredGrid, testParameter))
   end do

   testParameters = XdirParameter()

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testBins', testBins, testParameter))
   end do

   testParameters = XdirParameter()

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testPrint', testPrint, testParameter))
   end do


end function spaceMap_test_suite

