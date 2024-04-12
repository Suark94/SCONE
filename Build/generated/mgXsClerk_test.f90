module mgXsClerk_test

  use numPrecision
  use endfConstants
  use genericProcedures,         only : numToChar
  use mgXsClerk_class,           only : mgXsClerk
  use particle_class,            only : particle
  use particleDungeon_class,     only : particleDungeon
  use dictionary_class,          only : dictionary
  use scoreMemory_class,         only : scoreMemory
  use testNeutronDatabase_class, only : testNeutronDatabase
  use outputFile_class,          only : outputFile
  use pFUnit_mod

  implicit none

  !@testCase
    type, extends(TestCase) :: test_mgXsClerk
      private
      type(mgXsClerk)  :: clerk_test1
      type(mgXsClerk)  :: clerk_test2
      type(testNeutronDatabase) :: nucData
    contains
      procedure :: setUp
      procedure :: tearDown
  end type test_mgXsClerk

contains

  !!
  !! Sets up test_mgXsClerk object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_mgXsClerk), intent(inout) :: this
    type(dictionary)                     :: tempDict, energyDict, spaceDict
    character(nameLen)                   :: temp

    ! Build energy map
    call energyDict % init(3)
    call energyDict % store('type', 'energyMap')
    call energyDict % store('grid', 'unstruct')
    call energyDict % store('bins', [1.0E-03_defReal, ONE, 10.0_defReal])

    ! Build material map
    call spaceDict % init(2)
    call spaceDict % store('type', 'testMap')
    call spaceDict % store('maxIdx', 2)

    ! Define first clerk, with high order scattering and spatial map
    call tempDict % init(2)
    call tempDict % store('energyMap', energyDict)
    call tempDict % store('spaceMap', spaceDict)

    temp = 'MGxs1'
    call this % clerk_test1 % init(tempDict, temp)
    call spaceDict % kill()
    call energyDict % kill()
    call tempDict % kill()

    ! Build energy map
    call energyDict % init(3)
    call energyDict % store('type', 'energyMap')
    call energyDict % store('grid', 'unstruct')
    call energyDict % store('bins', [1.0E-11_defReal, 0.6_defReal, 1.2_defReal, 20.0_defReal])

    ! Define second clerk, without spatial map and high order scattering
    call tempDict % init(2)
    call tempDict % store('energyMap', energyDict)
    call tempDict % store('PN', 0)

    temp = 'MGxs2'
    call this % clerk_test2 % init(tempDict, temp)
    call energyDict % kill()
    call tempDict % kill()

    ! Build test neutronDatabase
    call this % nucData % build(ONE, captureXS = 2.0_defReal, &
                                fissionXS = 1.5_defReal, nuFissionXS = 3.0_defReal)

  end subroutine setUp

  !!
  !! Kills test case object
  !!
  subroutine tearDown(this)
    class(test_mgXsClerk), intent(inout) :: this

    call this % clerk_test1 % kill()
    call this % clerk_test2 % kill()
    call this % nucData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Scoring test for clerk 1
  !!
!@Test
  subroutine testScoring_clerk1(this)
    class(test_mgXsClerk), intent(inout) :: this
    character(:),allocatable             :: case
    type(scoreMemory)                    :: mem
    type(particle)                       :: p
    type(particleDungeon)                :: pit
    type(outputFile)                     :: out
    real(defReal), dimension(:,:), allocatable :: fiss, capt, transFL, transOS, &
                                                  nu, chi, P0, P1, P2, P3, P4,  &
                                                  P5, P6, P7, prod
    real(defReal), parameter :: TOL = 1.0E-9

    ! Configure memory
    call mem % init(1000_longInt, 1)
    call this % clerk_test1 % setMemAddress(1_longInt)

    ! Configure particle dungeon
    call pit % init(3)

    p % isMG = .false.
    p % w    = 0.5_defReal
    p % E    = 0.3_defReal
    call p % setMatIdx(2)
    call pit % detain(p)

    p % E = 3.0_defReal
    call pit % detain(p)

    ! Scoring
    call this % clerk_test1 % reportInColl(p, this % nucData, mem)
    call this % clerk_test1 % reportCycleEnd(pit, mem)

    p % preCollision % wgt = 0.2_defReal
    p % preCollision % E   = 3.0_defReal
    p % preCollision % matIdx = 2
    p % preCollision % isMG   = .false.
    p % E = 0.1_defReal

    call this % clerk_test1 % reportOutColl(p, N_2N, 0.75_defReal, this % nucData, mem)
    call mem % closeCycle(ONE)

    ! Process and get results
    call this % clerk_test1 % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)
    call this % clerk_test1 % processPN(mem, P2, P3, P4, P5, P6, P7)

    ! Verify results of scoring
#line 148 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, TWO, ZERO], capt(1,:), TOL, 'Capture XS' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 149 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, 1.5_defReal, ZERO], fiss(1,:), TOL, 'Fission XS' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 150 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, TWO, ZERO], nu(1,:), TOL, 'NuFission XS' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 150) )
  if (anyExceptions()) return
#line 151 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 151 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, HALF, HALF], chi(1,:), TOL, 'Chi' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 152 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, 4.0_defReal, ZERO], transOS(1,:), TOL, 'Transport XS O.S.' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 153 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, 5.5_defReal, ZERO], transFL(1,:), TOL, 'Transport XS F.L.' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 154 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, TWO, ZERO, ZERO], P0(1,:), TOL, 'P0' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 155 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, 1.5_defReal, ZERO, ZERO], P1(1,:), TOL, 'P1' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 156 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ONE, ONE, ONE, ONE, ONE, TWO, ONE, ONE], prod(1,:), TOL, 'prod' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

#line 158 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, 0.6875_defReal, ZERO, ZERO], P2(1,:), TOL, 'P2' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 159 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.140625_defReal, ZERO, ZERO], P3(1,:), TOL, 'P3' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 160 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.7001953125_defReal, ZERO, ZERO], P4(1,:), TOL, 'P4' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 161 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.8327636719_defReal, ZERO, ZERO], P5(1,:), TOL, 'P5' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 162 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.5615539551_defReal, ZERO, ZERO], P6(1,:), TOL, 'P6' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 163 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, -0.0683670044_defReal, ZERO, ZERO], P7(1,:), TOL, 'P7' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

    ! Test getting size
#line 166 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual(100, this % clerk_test1 % getSize(),'Test getSize():', &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

    ! Test correctness of output calls
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk_test1 % print(out, mem)
#line 171 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertTrue(out % isValid(), 'Test print():', &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

  end subroutine testScoring_clerk1

  !!
  !! Scoring test for clerk 2
  !!
!@Test
  subroutine testScoring_clerk2(this)
    class(test_mgXsClerk), intent(inout) :: this
    character(:),allocatable             :: case
    type(scoreMemory)                    :: mem
    type(particle)                       :: p
    type(particleDungeon)                :: pit
    type(outputFile)                     :: out
    real(defReal), dimension(:,:), allocatable :: fiss, capt, transFL, transOS, &
                                                  nu, chi, P0, P1, prod
    real(defReal), parameter :: TOL = 1.0E-9

    ! Configure memory
    call mem % init(1000_longInt, 1)
    call this % clerk_test2 % setMemAddress(1_longInt)

    ! Configure particle dungeon
    call pit % init(3)

    p % isMG = .false.
    p % w    = 0.5_defReal
    p % E    = 3.0_defReal
    call pit % detain(p)

    p % E = 0.3_defReal
    call pit % detain(p)

    ! Scoring
    call this % clerk_test2 % reportInColl(p, this % nucData, mem)
    call this % clerk_test2 % reportCycleEnd(pit, mem)

    p % preCollision % wgt = 0.2_defReal
    p % preCollision % E   = 0.3_defReal
    p % preCollision % isMG   = .false.
    p % E = 1.1_defReal

    call this % clerk_test2 % reportOutColl(p, N_2N, 0.75_defReal, this % nucData, mem)
    call mem % closeCycle(ONE)

    ! Process and get results
    call this % clerk_test2 % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)

    ! Verify results of scoring
#line 221 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, TWO], capt(1,:), TOL, 'Capture XS' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 222 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, 1.5_defReal], fiss(1,:), TOL, 'Fission XS' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 223 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, TWO], nu(1,:), TOL, 'NuFission XS' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 223) )
  if (anyExceptions()) return
#line 224 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 224 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([HALF, ZERO, HALF], chi(1,:), TOL, 'Chi' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 224) )
  if (anyExceptions()) return
#line 225 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 225 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, 4.0_defReal], transOS(1,:), TOL, 'Transport XS O.S.' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 225) )
  if (anyExceptions()) return
#line 226 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 226 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, 5.5_defReal], transFL(1,:), TOL, 'Transport XS F.L.' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 227 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, TWO, ZERO], P0(1,:), TOL, 'P0' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 227) )
  if (anyExceptions()) return
#line 228 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 228 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, 1.5_defReal, ZERO], P1(1,:), TOL, 'P1' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 228) )
  if (anyExceptions()) return
#line 229 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
#line 229 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual([ONE, ONE, ONE, ONE, ONE, ONE, ONE, TWO, ONE], prod(1,:), TOL, 'prod' , &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 229) )
  if (anyExceptions()) return
#line 230 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

    ! Test getting size
#line 232 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertEqual(48, this % clerk_test2 % getSize(),'Test getSize():', &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 232) )
  if (anyExceptions()) return
#line 233 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

    ! Test correctness of output calls
    call out % init('dummyPrinter', fatalErrors = .false.)
    call this % clerk_test2 % print(out, mem)
#line 237 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"
  call assertTrue(out % isValid(), 'Test print():', &
 & location=SourceLocation( &
 & 'mgXsClerk_test.f90', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/home/mhgk4/SCONE/Tallies/TallyClerks/Tests/mgXsClerk_test.f90"

  end subroutine testScoring_clerk2

end module mgXsClerk_test

module WrapmgXsClerk_test
   use pFUnit_mod
   use mgXsClerk_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_mgXsClerk) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use mgXsClerk_test
        class (test_mgXsClerk), intent(inout) :: this
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

end module WrapmgXsClerk_test

function mgXsClerk_test_suite() result(suite)
   use pFUnit_mod
   use mgXsClerk_test
   use WrapmgXsClerk_test
   type (TestSuite) :: suite

   suite = newTestSuite('mgXsClerk_test_suite')

   call suite%addTest(makeCustomTest('testScoring_clerk1', testScoring_clerk1))

   call suite%addTest(makeCustomTest('testScoring_clerk2', testScoring_clerk2))


end function mgXsClerk_test_suite

