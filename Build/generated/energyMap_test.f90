module energyMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use outputFile_class,        only : outputFile

  use energyMap_class,       only : energyMap

  implicit none


!@testCase
  type, extends(TestCase) :: test_energyMap
    private
    type(energyMap) :: map_lin
    type(energyMap) :: map_log
    type(energyMap) :: map_predef
    type(energyMap) :: map_unstruct
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyMap

  real(defReal),dimension(*), parameter :: UNSTRUCT_GRID = [ 0.00000000001_defReal, &
               0.00000003_defReal, 0.000000058_defReal, 0.00000014_defReal, 0.00000028_defReal, &
               0.00000035_defReal, 0.000000625_defReal, 0.000000972_defReal, 0.00000102_defReal,&
               0.000001097_defReal, 0.00000115_defReal, 0.000001855_defReal, 0.000004_defReal,&
               0.000009877_defReal, 0.000015968_defReal, 0.000148728_defReal, 0.00553_defReal,&
               0.009118_defReal, 0.111_defReal, 0.5_defReal, 0.821_defReal, 1.353_defReal, &
               2.231_defReal, 3.679_defReal, 6.0655_defReal, 10.0_defReal]



contains

  !!
  !! Sets up test_energyMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_energyMap), intent(inout) :: this
    type(dictionary)                     :: tempDict

    ! Build map lin
    call tempDict % init(4)
    call tempDict % store('grid','lin')
    call tempDict % store('min', 0.01_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 20)

    call this % map_lin % init(tempDict)
    call tempDict % kill()

    ! Build map log
    call tempDict % init(4)
    call tempDict % store('grid','log')
    call tempDict % store('min', 1.0E-7_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 20)

    call this % map_log % init(tempDict)
    call tempDict % kill()

    ! Build map predef
    call tempDict % init(2)
    call tempDict % store('grid','predef')
    call tempDict % store('name', 'casmo23')

    call this % map_predef % init(tempDict)
    call tempDict % kill()

    ! Build map log
    call tempDict % init(2)
    call tempDict % store('grid','unstruct')
    call tempDict % store('bins', UNSTRUCT_GRID)

    call this % map_unstruct % init(tempDict)
    call tempDict % kill()


  end subroutine setUp

  !!
  !! Kills test_energyMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_energyMap), intent(inout) :: this

    call this % map_lin % kill()
    call this % map_log % kill()
    call this % map_predef % kill()
    call this % map_unstruct % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test Linear grid
  !!
!@Test
  subroutine testLinearGrid(this)
    class(test_energyMap), intent(inout) :: this
    real(defReal),dimension(6),parameter :: E = [7.5774_defReal, 9.3652_defReal, 3.9223_defReal, &
                                                 6.5548_defReal, 1.7119_defReal, 20.0_defReal]
    integer(shortInt),dimension(6),parameter :: RES_IDX = [16, 19, 8, 14, 4, 0]
    integer(shortInt),dimension(6)           :: idx
    type(particleState),dimension(6)         :: states

    states % E = E
    idx = this % map_lin % map(states)
#line 114 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

  end subroutine testLinearGrid

  !!
  !! Test Log grid
  !!
!@Test
  subroutine testLogGrid(this)
    class(test_energyMap), intent(inout) :: this
    real(defReal),dimension(6),parameter :: E = [0.0445008907555061_defReal,   &
                                                 1.79747463687278e-07_defReal, &
                                                 1.64204055725811e-05_defReal, &
                                                 2.34083673923110e-07_defReal, &
                                                 5.98486350302033e-07_defReal, &
                                                 20.00000000000000000_defReal]
    integer(shortInt),dimension(6),parameter :: RES_IDX = [15, 1, 6, 1, 2, 0]
    integer(shortInt),dimension(6)           :: idx
    type(particleState),dimension(6)         :: states

    states % E = E
    idx = this % map_log % map(states)
#line 136 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

  end subroutine testLogGrid

  !!
  !! Test Predefined grid
  !!
!@Test
  subroutine testPredefGrid(this)
    class(test_energyMap), intent(inout) :: this
    real(defReal),dimension(6),parameter :: E = [0.0445008907555061_defReal,   &
                                                 1.79747463687278e-07_defReal, &
                                                 1.64204055725811e-05_defReal, &
                                                 2.34083673923110e-07_defReal, &
                                                 5.98486350302033e-07_defReal, &
                                                 20.00000000000000000_defReal]
    integer(shortInt),dimension(6),parameter :: RES_IDX = [16, 4, 13, 4, 6, 0]
    integer(shortInt),dimension(6)           :: idx
    type(particleState),dimension(6)         :: states

    states % E = E
    idx = this % map_predef % map(states)
#line 158 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

  end subroutine testPredefGrid

  !!
  !! Test Unstruct grid
  !!
!@Test
  subroutine testUnstructGrid(this)
    class(test_energyMap), intent(inout) :: this
    real(defReal),dimension(6),parameter :: E = [0.0761191517392624_defReal,   &
                                                 0.00217742635754091_defReal,  &
                                                 6.38548311340975e-08_defReal, &
                                                 2.52734532533842_defReal,     &
                                                 2.59031729968032e-11_defReal, &
                                                 20.00000000000000000_defReal]
    integer(shortInt),dimension(6),parameter :: RES_IDX = [18, 16, 3, 23, 1, 0]
    integer(shortInt),dimension(6)           :: idx
    type(particleState),dimension(6)         :: states

    states % E = E
    idx = this % map_unstruct % map(states)
#line 180 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

  end subroutine testUnstructGrid

  !!
  !! Test MG particle behaviour
  !!
!@Test
  subroutine testMGParticle(this)
    class(test_energyMap), intent(inout) :: this
    type(particleState)                  :: state
    integer(shortInt)                    :: idx

    state % isMG = .true.

    ! Linear energyMap
    idx = this % map_lin % map(state)
#line 197 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0, idx,'Linear energy Map', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 197) )
  if (anyExceptions()) return
#line 198 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

    ! Log energyMap
    idx = this % map_log % map(state)
#line 201 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0, idx,'Log energy Map', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

    ! Predef energyMap
    idx = this % map_predef % map(state)
#line 205 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0, idx,'Predef energy Map', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 205) )
  if (anyExceptions()) return
#line 206 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

    ! Unstructured energyMap
    idx = this % map_unstruct % map(state)
#line 209 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0, idx,'Unstructured energy Map', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

  end subroutine testMGParticle

  !!
  !! Test bin number retrival
  !!
!@Test
  subroutine testBinNumber(this)
    class(test_energyMap), intent(inout) :: this

    ! Linear energyMap
#line 221 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(20, this % map_lin % bins(1),'1st Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 222 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(20, this % map_lin % bins(0),'All bins', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 223 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0,  this % map_lin % bins(-3),'Invalid Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 223) )
  if (anyExceptions()) return
#line 224 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

    ! Log energyMap
#line 226 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(20, this % map_log % bins(1),'1st Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 227 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(20, this % map_log % bins(0),'All bins', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 227) )
  if (anyExceptions()) return
#line 228 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 228 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0,  this % map_log % bins(-3),'Invalid Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 228) )
  if (anyExceptions()) return
#line 229 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

    ! Predef energyMap
#line 231 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(23, this % map_predef % bins(1),'1st Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 231) )
  if (anyExceptions()) return
#line 232 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 232 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(23, this % map_predef % bins(0),'All bins', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 232) )
  if (anyExceptions()) return
#line 233 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 233 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0,  this % map_predef % bins(-3),'Invalid Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 233) )
  if (anyExceptions()) return
#line 234 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

    ! Unstructured energyMap
#line 236 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(25, this % map_unstruct % bins(1),'1st Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 236) )
  if (anyExceptions()) return
#line 237 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 237 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(25, this % map_unstruct % bins(0),'All bins', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
#line 238 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertEqual(0,  this % map_unstruct % bins(-3),'Invalid Dimension', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 238) )
  if (anyExceptions()) return
#line 239 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
!@Test
  subroutine testPrint(this)
    class(test_energyMap), intent(inout) :: this
    type(outputFile)                     :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_lin % print(out)
#line 254 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertTrue(out % isValid(),'Linear map case', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 254) )
  if (anyExceptions()) return
#line 255 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
    call out % reset()

    call this % map_log % print(out)
#line 258 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertTrue(out % isValid(),'Logarithmic map case', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 258) )
  if (anyExceptions()) return
#line 259 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
    call out % reset()

    call this % map_predef % print(out)
#line 262 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertTrue(out % isValid(),'Predefined map case', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 262) )
  if (anyExceptions()) return
#line 263 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
    call out % reset()

    call this % map_unstruct % print(out)
#line 266 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
  call assertTrue(out % isValid(),'Unstructured map case', &
 & location=SourceLocation( &
 & 'energyMap_test.f90', &
 & 266) )
  if (anyExceptions()) return
#line 267 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/energyMap_test.f90"
    call out % reset()

  end subroutine testPrint


end module energyMap_test

module WrapenergyMap_test
   use pFUnit_mod
   use energyMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_energyMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use energyMap_test
        class (test_energyMap), intent(inout) :: this
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

end module WrapenergyMap_test

function energyMap_test_suite() result(suite)
   use pFUnit_mod
   use energyMap_test
   use WrapenergyMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('energyMap_test_suite')

   call suite%addTest(makeCustomTest('testLinearGrid', testLinearGrid))

   call suite%addTest(makeCustomTest('testLogGrid', testLogGrid))

   call suite%addTest(makeCustomTest('testPredefGrid', testPredefGrid))

   call suite%addTest(makeCustomTest('testUnstructGrid', testUnstructGrid))

   call suite%addTest(makeCustomTest('testMGParticle', testMGParticle))

   call suite%addTest(makeCustomTest('testBinNumber', testBinNumber))

   call suite%addTest(makeCustomTest('testPrint', testPrint))


end function energyMap_test_suite

