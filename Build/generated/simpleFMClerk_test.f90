module simpleFMClerk_test

  use numPrecision
  use tallyResult_class,              only : tallyResult
  use simpleFMClerk_class,            only : simpleFMClerk, FMResult
  use particle_class,                 only : particle, particleState, P_NEUTRON
  use particleDungeon_class,          only : particleDungeon
  use dictionary_class,               only : dictionary
  use scoreMemory_class,              only : scoreMemory
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use outputFile_class,               only : outputFile
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_simpleFMClerk
    private
    type(simpleFMClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_simpleFMClerk

contains

  !!
  !! Sets up test_simpleFMClerk object we can use in a number of tests
  !!
  !! Simple 3x3 fission matrix divided with test map
  !!
  subroutine setUp(this)
    class(test_simpleFMClerk), intent(inout) :: this
    type(dictionary)                         :: dict
    type(dictionary)                         :: mapDict
    character(nameLen)                       :: name

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',3)


    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type','simpleFMClerk')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict,name)


    call mapDict % kill()
    call dict % kill()
  end subroutine setUp

  !!
  !! Kills test_simpleFMClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_simpleFMClerk), intent(inout) :: this

    call this % clerk % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


  !!
  !! Test correctness in a simple use case
  !!
!@Test
  subroutine testSimpleUseCase(this)
    class(test_simpleFMClerk), intent(inout) :: this
    type(scoreMemory)                        :: mem
    type(particle)                           :: p
    type(particleState)                      :: phase
    type(particleDungeon)                    :: pop
    type(testNeutronDatabase)                :: xsData
    real(defReal)                            :: val
    class(tallyResult), allocatable          :: res
    real(defReal), parameter :: TOL = 1.0E-7


    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

    ! Crate dungeon of original events
    ! One particle born in matIdx 1 and other in 2
    call pop % init(3)

    phase % wgt = ONE
    phase % matIdx = 2
    call pop % detain(phase)

    phase % wgt = ONE
    phase % matIdx = 1
    call pop % detain(phase)

    call this % clerk % reportCycleStart(pop, mem)

    ! Score some events
    p % type = P_NEUTRON

    call p % setMatIdx(2)
    p % w = 0.7
    p % preHistory % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    call p % setMatIdx(1)
    p % w = 1.1
    p % preHistory % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    call p % setMatIdx(1)
    p % w = 1.0
    p % preHistory % matIdx = 1
    call this % clerk % reportInColl(p, xsData, mem)

    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results

    ! Fission matrix
    ! 1 -> 1 Transition
    call mem % getResult(val, 1_longInt)
#line 136 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(1.818181818181_defReal ,val, TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

    ! 1 -> 2 Transition
    call mem % getResult(val, 2_longInt)
#line 140 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(ZERO, val, TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

    ! 1 -> 3 Transition
    call mem % getResult(val, 3_longInt)
#line 144 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(ZERO, val, TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 144) )
  if (anyExceptions()) return
#line 145 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

    ! 2 -> 1 Transition
    call mem % getResult(val, 4_longInt)
#line 148 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(2.0_defReal, val, TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

    ! 2 -> 2 Transition
    call mem % getResult(val, 5_longInt)
#line 152 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(1.27272727272727_defReal, val, TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

    ! Verify run-time result
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (FMresult)
#line 159 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(3, res % N, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! 1 -> 1 Transition
#line 162 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(1.818181818181_defReal ,res % FM(1,1,1) , TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! 1 -> 2 Transition
#line 165 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(ZERO, res  % FM(2,1,1), TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! 1 -> 3 Transition
#line 168 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(ZERO, res % FM(3,1,1), TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 168) )
  if (anyExceptions()) return
#line 169 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! 2 -> 1 Transition
#line 171 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(2.0_defReal, res % FM(1,2,1), TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! 2 -> 2 Transition
#line 174 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(1.27272727272727_defReal, res % FM(2,2,1), TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! Clean all entries
        res % FM = ZERO

      class default
#line 180 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(1,2, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
    end select

    ! Get result again -> verify correcness of reallocation logic by code coverage
    call this % clerk % getResult(res, mem)
    select type(res)
      class is (FMresult)
        ! 1 -> 1 Transition
#line 188 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(1.818181818181_defReal ,res % FM(1,1,1) , TOL, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

        ! Change size of matrix
        res % N = 2
        deallocate(res % FM)
        allocate(res % FM(2,2,1))

    end select
    ! Get result yet again. This time with wrong size
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (FMresult)
#line 201 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual(3, res % N, &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
#line 202 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertEqual([3,3,2], shape(res % FM), &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
    end select

    ! Clean
    call xsData % kill()
    call pop % kill()

  end subroutine testSimpleUseCase



  !!
  !! Test correctness of the printing calls
  !!
!@Test
  subroutine testPrintingCorrectness(this)
    class(test_simpleFMClerk), intent(inout) :: this
    type(outputFile)                         :: outF
    type(scoreMemory)                        :: mem

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print (outF, mem)

#line 230 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"
  call assertTrue(outF % isValid(), &
 & location=SourceLocation( &
 & 'simpleFMClerk_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/simpleFMClerk_test.f90"

  end subroutine testPrintingCorrectness



end module simpleFMClerk_test

module WrapsimpleFMClerk_test
   use pFUnit_mod
   use simpleFMClerk_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_simpleFMClerk) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use simpleFMClerk_test
        class (test_simpleFMClerk), intent(inout) :: this
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

end module WrapsimpleFMClerk_test

function simpleFMClerk_test_suite() result(suite)
   use pFUnit_mod
   use simpleFMClerk_test
   use WrapsimpleFMClerk_test
   type (TestSuite) :: suite

   suite = newTestSuite('simpleFMClerk_test_suite')

   call suite%addTest(makeCustomTest('testSimpleUseCase', testSimpleUseCase))

   call suite%addTest(makeCustomTest('testPrintingCorrectness', testPrintingCorrectness))


end function simpleFMClerk_test_suite

