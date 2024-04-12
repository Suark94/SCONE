module keffAnalogClerk_test

  use numPrecision
  use tallyResult_class,     only : tallyResult
  use keffAnalogClerk_class, only : keffAnalogClerk, keffResult
  use particle_class,        only : particle
  use particleDungeon_class, only : particleDungeon
  use dictionary_class,      only : dictionary
  use scoreMemory_class,     only : scoreMemory
  use outputFile_class,      only : outputFile
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_keffAnalogClerk
    private
    type(keffAnalogClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_keffAnalogClerk

contains

  !!
  !! Sets up test_keffAnalogClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(dictionary)                           :: dict
    character(nameLen)                         :: name

    call dict % init(2)
    call this % clerk % init(dict, name)

    call dict % kill()

  end subroutine setUp

  !!
  !! Kills test_keffAnalogClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_keffAnalogClerk), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test keffAnalogClerk for a 1 Cycle Batch case
  !!
!@Test
  subroutine test1CycleBatch(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(scoreMemory)                          :: mem
    type(particleDungeon)                      :: pit
    type(particle)                             :: p
    real(defReal)                              :: k, STD
    class(tallyResult), allocatable            :: res
    real(defReal),parameter :: TOL = 1.0E-9

    ! Initialise objects
    call pit % init(4)
    call mem % init(2_longInt,1)
    call this  % clerk % setMemAddress(1_longInt)

    ! Start cycle 1
    p % w = 1000.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 1
    call pit % release(p)
    p % w = 1200.0_defReal
    call pit % detain(p)
    pit % k_eff = ONE

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 2
    call pit % release(p)
    p % w = 1000.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 2
    call pit % release(p)
    p % w = 900.0_defReal
    call pit % detain(p)
    pit % k_eff = 1.2_defReal

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Validate results

    ! Directly from memory
    call mem % getResult(k, STD, 1_longInt)
#line 104 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(1.1400_defReal, k, TOL, '1 Cycle Batch, keff from memory:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
#line 105 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(0.0600_defReal, STD, TOL, '1 Cycle Batch, keff STD from memory:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

    ! From result
    call this % clerk % getResult(res, mem)

    select type(res)
      type is(keffResult)
#line 112 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(1.1400_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
#line 113 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(0.0600_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

      class default
#line 116 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertTrue(.false.,'Result is not a keffResult', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 116) )
  if (anyExceptions()) return
#line 117 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

    end select

  end subroutine test1CycleBatch


  !!
  !! Test keffAnalogClerk for a 2 Cycle Batch case
  !!
!@Test
  subroutine test2CycleBatch(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(scoreMemory)                          :: mem
    type(particleDungeon)                      :: pit
    type(particle)                             :: p
    real(defReal)                              :: k, STD
    class(tallyResult), allocatable            :: res
    real(defReal),parameter :: TOL = 1.0E-9

    ! Initialise objects
    call pit % init(4)
    call mem % init(2_longInt,1, batchSize = 2 )
    call this  % clerk % setMemAddress(1_longInt)

    ! Start cycle 1
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 1
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    pit % k_eff = ONE

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 2
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 2
    call pit % release(p)
    p % w = 700.0_defReal
    call pit % detain(p)
    pit % k_eff = ONE

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 3
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 3
    call pit % release(p)
    p % w = 400.0_defReal
    call pit % detain(p)
    pit % k_eff = 1.2_defReal

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)

    ! Start cycle 4
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    call this % clerk % reportCycleStart(pit, mem)

    ! End cycle 4
    call pit % release(p)
    p % w = 500.0_defReal
    call pit % detain(p)
    pit % k_eff = 1.2_defReal

    call this % clerk % reportCycleEnd(pit,mem)
    call mem % closeCycle(0.8_defReal)


    ! Validate results

    ! Directly from memory
    call mem % getResult(k, STD, 1_longInt)
#line 205 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(1.1400_defReal, k, TOL, '1 Cycle Batch, keff from memory:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 205) )
  if (anyExceptions()) return
#line 206 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
#line 206 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(0.0600_defReal, STD, TOL, '1 Cycle Batch, keff STD from memory:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 206) )
  if (anyExceptions()) return
#line 207 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

    ! From result
    call this % clerk % getResult(res, mem)

    select type(res)
      type is(keffResult)
#line 213 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(1.1400_defReal, res % keff(1), TOL, '1 Cycle Batch, keff from result:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 213) )
  if (anyExceptions()) return
#line 214 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
#line 214 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(0.0600_defReal, res % keff(2), TOL, '1 Cycle Batch, keff STD from result:', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 214) )
  if (anyExceptions()) return
#line 215 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

      class default
#line 217 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertTrue(.false.,'Result is not a keffResult', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

    end select

  end subroutine test2CycleBatch

  !!
  !! Test functions: print & getSize
  !!
!@Test
  subroutine testMisc(this)
    class(test_keffAnalogClerk), intent(inout) :: this
    type(outputFile)                           :: out
    type(scoreMemory)                          :: mem

    ! Initialise objects
    call mem % init(2_longInt,1)
    call this  % clerk % setMemAddress(1_longInt)

    ! Test getting size
#line 237 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertEqual(1, this % clerk % getSize(),'Test getSize() :', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

    ! Test output printing correctness
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print(out, mem)
#line 242 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"
  call assertTrue(out % isValid(), 'Test print():', &
 & location=SourceLocation( &
 & 'keffAnalogClerk_test.f90', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/keffAnalogClerk_test.f90"

  end subroutine testMisc


end module keffAnalogClerk_test

module WrapkeffAnalogClerk_test
   use pFUnit_mod
   use keffAnalogClerk_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_keffAnalogClerk) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use keffAnalogClerk_test
        class (test_keffAnalogClerk), intent(inout) :: this
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

end module WrapkeffAnalogClerk_test

function keffAnalogClerk_test_suite() result(suite)
   use pFUnit_mod
   use keffAnalogClerk_test
   use WrapkeffAnalogClerk_test
   type (TestSuite) :: suite

   suite = newTestSuite('keffAnalogClerk_test_suite')

   call suite%addTest(makeCustomTest('test1CycleBatch', test1CycleBatch))

   call suite%addTest(makeCustomTest('test2CycleBatch', test2CycleBatch))

   call suite%addTest(makeCustomTest('testMisc', testMisc))


end function keffAnalogClerk_test_suite

