module sphericalMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use outputFile_class,        only : outputFile

  use sphericalMap_class,      only : sphericalMap

  implicit none


!@testCase
  type, extends(TestCase) :: test_sphericalMap
    private
    type(sphericalMap) :: map_lin_from_zero
    type(sphericalMap) :: map_lin_from_min
    type(sphericalMap) :: map_unstruct
    type(sphericalMap) :: map_equivol

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_sphericalMap

contains

  !!
  !! Sets up test_sphericalMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_sphericalMap), intent(inout) :: this
    type(dictionary)                     :: tempDict

    ! Build map with default origin & minimum radius
    call tempDict % init(4)
    call tempDict % store('grid','lin')
    call tempDict % store('Rmax', 10.0_defReal)
    call tempDict % store('N', 20)

    call this % map_lin_from_zero % init(tempDict)
    call tempDict % kill()

    ! Build map with diffrent origin & minimum radius
    call tempDict % init(5)
    call tempDict % store('origin', [ONE, ONE, ONE])
    call tempDict % store('grid', 'lin')
    call tempDict % store('Rmin', 5.0_defReal)
    call tempDict % store('Rmax', 10.0_defReal)
    call tempDict % store('N', 5)

    call this % map_lin_from_min % init(tempDict)
    call tempDict % kill()

    ! Build map with unstruct bins
    call tempDict % init(2)
    call tempDict % store('grid', 'unstruct')
    call tempDict % store('bins', [4.0_defReal, 5.8_defReal, 7.3_defReal, 10.5_defReal, 15.0_defReal])

    call this % map_unstruct % init(tempDict)
    call tempDict % kill()

    ! Build map with equivolume bins
    call tempDict % init(4)
    call tempDict % store('grid', 'equivolume')
    call tempDict % store('Rmin', 2.0_defReal)
    call tempDict % store('Rmax', 20.0_defReal)
    call tempDict % store('N', 8)

    call this % map_equivol % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_sphericalMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_sphericalMap), intent(inout) :: this

    call this % map_lin_from_zero % kill()
    call this % map_lin_from_min % kill()
    call this % map_unstruct % kill()
    call this % map_equivol % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test default-initialised grid
  !!
!@Test
  subroutine testFromOrigin(this)
    class(test_sphericalMap), intent(inout) :: this
    real(defReal),dimension(4),parameter     :: r = [0.4_defReal, 3.58_defReal, 8.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [1, 8, 18, 0]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    idx = this % map_lin_from_zero % map(states)
#line 111 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"

  end subroutine testFromOrigin

  !!
  !! Test grid with shifted origin & minimum radius
  !!
!@Test
  subroutine testFromMin(this)
    class(test_sphericalMap), intent(inout) :: this
    real(defReal),dimension(4),parameter     :: r = [1.5_defReal, 5.5_defReal, 8.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [0, 1, 4, 0]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    ! Shift the origin
    states(:) % r(1) = states(:) % r(1) + ONE
    states(:) % r(2) = states(:) % r(2) + ONE
    states(:) % r(3) = states(:) % r(3) + ONE

    idx = this % map_lin_from_min % map(states)
#line 139 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"

  end subroutine testFromMin

  !!
  !! Test grid with unstruct bins
  !!
!@Test
  subroutine testUnstruct(this)
    class(test_sphericalMap), intent(inout) :: this
    real(defReal),dimension(4),parameter     :: r = [4.5_defReal, 15.5_defReal, 8.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [1, 0, 3, 4]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    idx = this % map_unstruct % map(states)
#line 162 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"

  end subroutine testUnstruct

  !!
  !! Test grid with equivolume bins
  !!
!@Test
  subroutine testEquiVol(this)
    class(test_sphericalMap), intent(inout) :: this
    real(defReal),dimension(4),parameter     :: r = [1.5_defReal, 5.5_defReal, 18.9_defReal, 11.0_defReal]
    real(defReal), dimension(4),parameter    :: phi = [1.4_defReal, 3.98_defReal, 0.5_defReal, PI/2]
    real(defReal), dimension(4), parameter   :: tht = [ZERO, PI/2, PI/4, -PI/2]
    integer(shortInt),dimension(4),parameter :: RES_IDX = [0, 1, 7, 2]
    integer(shortInt),dimension(4)           :: idx
    type(particleState),dimension(4)         :: states

    ! Initialise states
    states(:) % r(1) = r * cos(phi) * sin(tht)
    states(:) % r(2) = r * sin(phi) * sin(tht)
    states(:) % r(3) = r * cos(tht)

    idx = this % map_equivol % map(states)
#line 185 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(RES_IDX, idx, &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"

  end subroutine testEquiVol

  !!
  !! Test bin number retrival
  !!
!@Test
  subroutine testBinNumber(this)
    class(test_sphericalMap), intent(inout) :: this

#line 196 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(20, this % map_lin_from_zero % bins(1),'1st Dimension', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 196) )
  if (anyExceptions()) return
#line 197 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
#line 197 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(4,  this % map_unstruct % bins(1),'1st Dimension', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 197) )
  if (anyExceptions()) return
#line 198 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
#line 198 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(20, this % map_lin_from_zero % bins(0),'All bins', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 198) )
  if (anyExceptions()) return
#line 199 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
#line 199 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(0,  this % map_lin_from_zero % bins(-3),'Invalid Dimension', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 199) )
  if (anyExceptions()) return
#line 200 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"

    ! Get dimensionality
#line 202 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertEqual(1, this % map_lin_from_zero % dimensions(), &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"

  end subroutine testBinNumber

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
!@Test
  subroutine testPrint(this)
    class(test_sphericalMap), intent(inout) :: this
    type(outputFile)                        :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_lin_from_zero % print(out)
#line 218 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertTrue(out % isValid(),'Default-initialised map case', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
    call out % reset()

    call this % map_lin_from_min % print(out)
#line 222 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertTrue(out % isValid(),'Map with minimum R', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
    call out % reset()

    call this % map_unstruct % print(out)
#line 226 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertTrue(out % isValid(),'Map with unstruct bins', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
    call out % reset()

    call this % map_equivol % print(out)
#line 230 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
  call assertTrue(out % isValid(),'Map with equivolume bins', &
 & location=SourceLocation( &
 & 'sphericalMap_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/Tallies/TallyMaps/Tests/sphericalMap_test.f90"
    call out % reset()

  end subroutine testPrint


end module sphericalMap_test

module WrapsphericalMap_test
   use pFUnit_mod
   use sphericalMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_sphericalMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use sphericalMap_test
        class (test_sphericalMap), intent(inout) :: this
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

end module WrapsphericalMap_test

function sphericalMap_test_suite() result(suite)
   use pFUnit_mod
   use sphericalMap_test
   use WrapsphericalMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('sphericalMap_test_suite')

   call suite%addTest(makeCustomTest('testFromOrigin', testFromOrigin))

   call suite%addTest(makeCustomTest('testFromMin', testFromMin))

   call suite%addTest(makeCustomTest('testUnstruct', testUnstruct))

   call suite%addTest(makeCustomTest('testEquiVol', testEquiVol))

   call suite%addTest(makeCustomTest('testBinNumber', testBinNumber))

   call suite%addTest(makeCustomTest('testPrint', testPrint))


end function sphericalMap_test_suite

