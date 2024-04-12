module collisionProbabilityClerk_test

  use numPrecision
  use tallyResult_class,               only : tallyResult
  use collisionProbabilityClerk_class, only : collisionProbabilityClerk, CPMResult
  use particle_class,                  only : particle, particleState, P_NEUTRON
  use particleDungeon_class,           only : particleDungeon
  use dictionary_class,                only : dictionary
  use scoreMemory_class,               only : scoreMemory
  use testNeutronDatabase_class,       only : testNeutronDatabase
  use outputFile_class,                only : outputFile
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_collisionProbabilityClerk
    private
    type(collisionProbabilityClerk) :: clerk
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_collisionProbabilityClerk

contains

  !!
  !! Sets up test_collisionProbabilityClerk object we can use in a number of tests
  !!
  !! Simple 2x2 collision probability score (with outside region) divided with test map
  !!
  subroutine setUp(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this
    type(dictionary)                                     :: dict
    type(dictionary)                                     :: mapDict
    character(nameLen)                                   :: name

    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',2)


    ! Build intput dictionary
    call dict % init(2)
    call dict % store('type','collisionProbabilityClerk')
    call dict % store('map', mapDict)

    name = 'testClerk'
    call this % clerk % init(dict,name)


    call mapDict % kill()
    call dict % kill()
  end subroutine setUp

  !!
  !! Kills test_collisionProbabilityClerk object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_collisionProbabilityClerk), intent(inout) :: this

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
    class(test_collisionProbabilityClerk), intent(inout) :: this
    type(scoreMemory)                                    :: mem
    type(particle)                                       :: p
    type(particleState)                                  :: phase
    type(particleDungeon)                                :: pop
    type(testNeutronDatabase)                            :: xsData
    real(defReal)                                        :: val
    class(tallyResult), allocatable                      :: res
    real(defReal), parameter :: TOL = 1.0E-7


    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1, batchSize = 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Create test transport Nuclear Data
    call xsData % build(1.1_defReal, fissionXS = 1.1_defReal, nuFissionXS = 2.0_defReal)

    ! Create one particle that can be made to collide
    ! repeatedly in several materials

    ! Score some events
    p % type = P_NEUTRON

    ! Particle starts in material 2 and collides in material 2
    call p % setMatIdx(2)
    p % w = 0.7
    p % preCollision % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 1 and collides in material 2
    call p % setMatIdx(2)
    p % w = 1.1
    p % preCollision % matIdx = 1
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 1 and collides in material 1
    call p % setMatIdx(1)
    p % w = 1.0
    p % preCollision % matIdx = 1
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 2 and collides in material 1
    call p % setMatIdx(1)
    p % w = 1.4
    p % preCollision % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in material 2 and collides in another, unknown material
    call p % setMatIdx(7)
    p % w = 1.0
    p % preCollision % matIdx = 2
    call this % clerk % reportInColl(p, xsData, mem)

    ! Particle starts in an unknown material and collides in material 1
    call p % setMatIdx(1)
    p % w = 0.9
    p % preCollision % matIdx = 88
    call this % clerk % reportInColl(p, xsData, mem)
    call this % clerk % reportCycleEnd(pop, mem)

    ! Close cycle
    call mem % closeCycle(ONE)

    ! Verify results

    ! Collision probability matrix

    ! outside -> outside Transition
    call mem % getResult(val, 1_longInt)
#line 146 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ZERO ,val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! outside -> 1 Transition
    call mem % getResult(val, 2_longInt)
#line 150 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ONE ,val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 150) )
  if (anyExceptions()) return
#line 151 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! outside -> 2 Transition
    call mem % getResult(val, 3_longInt)
#line 154 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ZERO ,val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! 1 -> outside Transition
    call mem % getResult(val, 4_longInt)
#line 158 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ZERO ,val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! 1 -> 1 Transition
    call mem % getResult(val, 5_longInt)
#line 162 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.47619047619_defReal ,val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! 1 -> 2 Transition
    call mem % getResult(val, 6_longInt)
#line 166 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.52380952381_defReal, val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! 2 -> outside Transition
    call mem % getResult(val, 7_longInt)
#line 170 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.32258064516_defReal, val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! 2 -> 1 Transition
    call mem % getResult(val, 8_longInt)
#line 174 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.45161290322_defReal ,val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! 2 -> 2 Transition
    call mem % getResult(val, 9_longInt)
#line 178 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.22580645161_defReal, val, TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

    ! Verify run-time result
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (CPMresult)
#line 185 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(3, res % N, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! outside -> outside Transition
#line 188 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ZERO, res  % CPM(1,1,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! outside -> 1 Transition
#line 191 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ONE, res  % CPM(2,1,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! outside -> 2 Transition
#line 194 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ZERO, res  % CPM(3,1,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 194) )
  if (anyExceptions()) return
#line 195 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! 1 -> outside Transition
#line 197 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(ZERO, res  % CPM(1,2,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 197) )
  if (anyExceptions()) return
#line 198 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! 1 -> 1 Transition
#line 200 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.47619047619_defReal, res  % CPM(2,2,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 200) )
  if (anyExceptions()) return
#line 201 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! 1 -> 2 Transition
#line 203 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.52380952381_defReal, res  % CPM(3,2,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 203) )
  if (anyExceptions()) return
#line 204 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! 2 -> outside Transition
#line 206 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.32258064516, res  % CPM(1,3,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 206) )
  if (anyExceptions()) return
#line 207 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! 2 -> 1 Transition
#line 209 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.45161290322, res  % CPM(2,3,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! 2 -> 2 Transition
#line 212 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.22580645161, res  % CPM(3,3,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 212) )
  if (anyExceptions()) return
#line 213 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! Clean all entries
        res % CPM = ZERO

      class default
#line 218 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(1,2, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
    end select

    ! Get result again -> verify correcness of reallocation logic by code coverage
    call this % clerk % getResult(res, mem)
    select type(res)
      class is (CPMresult)
        ! 1 -> 1 Transition
#line 226 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(0.47619047619_defReal, res  % CPM(2,2,1), TOL, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

        ! Change size of matrix
        res % N = 2
        deallocate(res % CPM)
        allocate(res % CPM(2,2,1))

    end select

    ! Get result yet again to ensure the size was not incorrectly modified
    call this % clerk % getResult(res, mem)

    select type(res)
      class is (CPMresult)
#line 240 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual(3, res % N, &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
#line 241 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertEqual([3,3,2], shape(res % CPM), &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
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
    class(test_collisionProbabilityClerk), intent(inout) :: this
    type(outputFile)                         :: outF
    type(scoreMemory)                        :: mem

    ! Create score memory
    call mem % init(int(this % clerk % getSize(), longInt) , 1)
    call this % clerk % setMemAddress(1_longInt)

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk % print (outF, mem)

#line 269 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"
  call assertTrue(outF % isValid(), &
 & location=SourceLocation( &
 & 'collisionProbabilityClerk_test.f90', &
 & 269) )
  if (anyExceptions()) return
#line 270 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/collisionProbabilityClerk_test.f90"

  end subroutine testPrintingCorrectness

end module collisionProbabilityClerk_test

module WrapcollisionProbabilityClerk_test
   use pFUnit_mod
   use collisionProbabilityClerk_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_collisionProbabilityClerk) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use collisionProbabilityClerk_test
        class (test_collisionProbabilityClerk), intent(inout) :: this
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

end module WrapcollisionProbabilityClerk_test

function collisionProbabilityClerk_test_suite() result(suite)
   use pFUnit_mod
   use collisionProbabilityClerk_test
   use WrapcollisionProbabilityClerk_test
   type (TestSuite) :: suite

   suite = newTestSuite('collisionProbabilityClerk_test_suite')

   call suite%addTest(makeCustomTest('testSimpleUseCase', testSimpleUseCase))

   call suite%addTest(makeCustomTest('testPrintingCorrectness', testPrintingCorrectness))


end function collisionProbabilityClerk_test_suite

