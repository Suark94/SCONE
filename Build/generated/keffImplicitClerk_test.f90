module keffImplicitClerk_test

  use numPrecision
  use tallyCodes
  use endfConstants
  use tallyResult_class,       only : tallyResult
  use keffAnalogClerk_class,   only : keffResult
  use keffImplicitClerk_class, only : keffImplicitClerk
  use particle_class,          only : particle
  use particleDungeon_class,   only : particleDungeon
  use dictionary_class,        only : dictionary
  use scoreMemory_class,       only : scoreMemory
  use outputFile_class,        only : outputFile
  use pFUnit_mod

  use testNeutronDatabase_class, only : testNeutronDatabase

  implicit none

!@testCase
  type, extends(TestCase) :: test_keffImplicitClerk
    private
    type(keffImplicitClerk)   :: clerk
    type(testNeutronDatabase) :: nucData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_keffImplicitClerk

contains

  !!
  !! Sets up test_keffImplicitClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(dictionary)                             :: dict
    character(nameLen)                           :: name

    call dict % init(2)
    call this % clerk % init(dict, name)
    call dict % kill()

    call this % nucData % build(ONE, captureXS = 2.0_defReal, fissionXS = ONE, nuFissionXS = 3.0_defReal)

  end subroutine setUp

  !!
  !! Kills test_keffImplicitClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_keffImplicitClerk), intent(inout) :: this

    call this % nucData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test for 1 cycle batches
  !!
!@Test
  subroutine test1CycleBatch(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(particle)                               :: p
    type(particleDungeon)                        :: pit
    type(scoreMemory)                            :: mem
    class(tallyResult),allocatable               :: res
    real(defReal), parameter                     :: TOL = 1.0E-9_defReal

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call pit % init(4)


    ! Configure particle
    p % fate = leak_FATE

    !*** Start cycle 1
    ! Score implicit reaction rates
    p % w = 0.7_defReal
    call this % clerk % reportInColl(p, this % nucData, mem)

    ! Score analog production
    p % preCollision % wgt = 0.1_defReal
    call this % clerk % reportOutColl(p, N_2N, 0.5_defReal, this % nucData, mem)

    ! Score leakage
    p % w = 0.3_defReal
    call this % clerk % reportHist(p, this % nucData, mem)

    ! End cycle
    call pit % detain(p)
    call this % clerk % reportCycleEnd(pit, mem)
    call pit % release(p)
    call mem % closeCycle(ONE)

    !*** Start cycle 2
    ! Score implicit reaction rates
    p % w = 0.6_defReal
    call this % clerk % reportInColl(p, this % nucData, mem)

    ! Score analog production
    p % preCollision % wgt = 0.1_defReal
    call this % clerk % reportOutColl(p, N_2N, 0.5_defReal, this % nucData, mem)

    ! Score leakage
    p % w = 0.3_defReal
    call this % clerk % reportHist(p, this % nucData, mem)

    ! End cycle
    call pit % detain(p)
    call this % clerk % reportCycleEnd(pit, mem)
    call pit % release(p)
    call mem % closeCycle(ONE)

    ! Verify result
    call this % clerk % getResult(res, mem)
    select type(res)
      type is(keffResult)
#line 125 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"
  call assertEqual(0.906521739130435_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result:', &
 & location=SourceLocation( &
 & 'keffImplicitClerk_test.f90', &
 & 125) )
  if (anyExceptions()) return
#line 126 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"
#line 126 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"
  call assertEqual(0.006521739130435_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result:', &
 & location=SourceLocation( &
 & 'keffImplicitClerk_test.f90', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"

      class default
#line 129 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"
  call assertTrue(.false.,'Result is not a keffResult', &
 & location=SourceLocation( &
 & 'keffImplicitClerk_test.f90', &
 & 129) )
  if (anyExceptions()) return
#line 130 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"

    end select
  end subroutine test1CycleBatch

  !!
  !! Test getSize() and print
  !!
!@Test
  subroutine testMisc(this)
    class(test_keffImplicitClerk), intent(inout) :: this
    type(scoreMemory)                            :: mem
    type(outputFile)                             :: out

    ! Configure memory
    call mem % init(10_longInt, 1)
    call this % clerk % setMemAddress(1_longInt)
    call out % init('dummyPrinter', fatalErrors = .false.)

    ! Test getting size
#line 149 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"
  call assertEqual(5, this % clerk % getSize(),'Test getSize():', &
 & location=SourceLocation( &
 & 'keffImplicitClerk_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"

    ! Test correctness of output calls
    call this % clerk % print(out, mem)
#line 153 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"
  call assertTrue(out % isValid(), 'Test print():', &
 & location=SourceLocation( &
 & 'keffImplicitClerk_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffImplicitClerk_test.f90"

  end subroutine testMisc

end module keffImplicitClerk_test

module WrapkeffImplicitClerk_test
   use pFUnit_mod
   use keffImplicitClerk_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_keffImplicitClerk) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use keffImplicitClerk_test
        class (test_keffImplicitClerk), intent(inout) :: this
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

end module WrapkeffImplicitClerk_test

function keffImplicitClerk_test_suite() result(suite)
   use pFUnit_mod
   use keffImplicitClerk_test
   use WrapkeffImplicitClerk_test
   type (TestSuite) :: suite

   suite = newTestSuite('keffImplicitClerk_test_suite')

   call suite%addTest(makeCustomTest('test1CycleBatch', test1CycleBatch))

   call suite%addTest(makeCustomTest('testMisc', testMisc))


end function keffImplicitClerk_test_suite

