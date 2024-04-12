module scoreMemory_test
  use numPrecision
  use genericProcedures, only : numToChar
  use scoreMemory_class, only : scoreMemory
  use pFUnit_mod

  implicit none

!@testParameter(constructor = new_testNumber)
  type, extends(AbstractTestParameter) :: testNumber
    integer(shortInt) :: i
  contains
    procedure :: toString
  end type testNumber

!@testCase(constructor=newTest)
  type, extends(ParameterizedTestCase) :: test_scoreMemory
    private
    integer(longInt)                            :: Ncycles
    integer(shortInt)                           :: batchSize
    real(defReal),dimension(:), allocatable     :: scores
    integer(shortInt), dimension(:),allocatable :: scoresInt

  end type test_scoreMemory


contains

  !!
  !! Build new test parameter form integer
  !!
  function new_testNumber(i) result (tstNum)
    integer(shortInt) :: i
    type(testNumber)  :: tstNum

    tstNum % i = i

  end function new_testNumber

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(testNumber), intent(in) :: this
    character(:), allocatable :: string
    character(nameLen)        :: str

    write (str,*) this % i
    string = str

  end function toString

  !!
  !! Construct test case
  !!
  !!
  !!
  function newTest(testParam) result(tst)
    type(testNumber), intent(in)     :: testParam
    type(test_scoreMemory)           :: tst
    real(defReal),dimension(200)     :: random
    integer(shortInt)                :: seed, i
    integer(shortInt),parameter      :: A = 2469   ! Multiplier of LC PRNG
    integer(shortInt),parameter      :: M = 65521  ! Modulus of PRNG

    ! Load batchSize
    tst % batchSize = testParam % i
    tst % Ncycles   = 10 * tst % batchSize

    ! Generate a vector of 20 pseudo-random numbers in <0;1>
    ! Generator is not sophisticated but robust
    seed = 9294
    do i=1,200
      seed = mod(A * seed , M)
      random(i)    = seed / real(M,defReal)
    end do

    ! Generate some scores and calculate their sum and sum of squares
    tst  % scores    = TWO + sin(PI * random - PI/2)
    tst % scoresInt = int(random * 100, shortInt)

  end function newTest

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test acoring for a case with batchSize == 1
  !! Look at the end of the file to find MATLAB SCRIPT used to generate reference values
  !!
!@Test(cases=[1])
  subroutine testScoring(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem
    integer(shortInt)                      :: i, j
    real(defReal)                          :: res1, res2, STD
    real(defReal), parameter :: TOL = 1.0E-9

    ! Initialise score memory
    call mem % init(7_longInt, 1, batchSize = this % batchSize)

    ! Test getting batchSize
#line 104 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(this % batchSize, mem % getBatchSize(),'Test getBatchSize() :', &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Score in
    do i=1,10
      ! Score
      do j=20*(i-1)+1,20 * i
        call mem % score(this % scores(j), 1_longInt)
        call mem % score(this % scoresInt(j), 2_longInt)
        call mem % score(int(this % scoresInt(j),longInt),3_longInt)
        call mem % accumulate(this % scores(j), 4_longInt)
        call mem % accumulate(this % scoresInt(j), 5_longInt)
        call mem % accumulate(int(this % scoresInt(j),longInt),6_longInt)

      end do
      ! Close a single bin with diffrent normalisation
      call mem % closeBin(1.2_defReal, 3_longInt)

      ! Close Cycle
      call mem % closeCycle(0.7_defReal)

    end do

    ! Get results from bin 1
    call mem % getResult(res1, 1_longInt)
    call mem % getResult(res2, STD, 1_longInt)

#line 130 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(26.401471259728442_defReal, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 131 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(26.401471259728442_defReal, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 132 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(0.645969443981583_defReal, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from bin 2
    call mem % getResult(res1, 2_longInt)
    call mem % getResult(res2, STD, 2_longInt)

#line 138 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(623.0_defReal, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 139 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(623.0_defReal, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 140 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(27.982494527829360_defReal, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from bin 3
    call mem % getResult(res1, 3_longInt)
    call mem % getResult(res2, STD, 3_longInt)

#line 146 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(1068.0_defReal, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 147 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(1068.0_defReal, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 147) )
  if (anyExceptions()) return
#line 148 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 148 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(47.969990619136050_defReal, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from bin 4
    call mem % getResult(res1, 4_longInt, 200)
    call mem % getResult(res2, STD, 4_longInt, 200)

#line 154 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(1.885819375694888_defReal, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 155 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(1.885819375694888_defReal, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 156 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(0.049102082638055_defReal, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from bin 5
    call mem % getResult(res1, 5_longInt, 200)
    call mem % getResult(res2, STD, 5_longInt, 200)

#line 162 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(44.500000000000000_defReal, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 163 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(44.500000000000000_defReal, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 164 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(2.015580019267494_defReal, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from bin 6
    call mem % getResult(res1, 6_longInt, 200)
    call mem % getResult(res2, STD, 6_longInt, 200)

#line 170 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(44.500000000000000_defReal, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 171 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(44.500000000000000_defReal, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 172 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(2.015580019267494_defReal, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 172) )
  if (anyExceptions()) return
#line 173 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from an empty bin 7
    call mem % getResult(res1, 7_longInt)
    call mem % getResult(res2, STD, 7_longInt)

#line 178 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 179 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 180 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Get results from invalid bins
    call mem % getResult(res1, -7_longInt)
    call mem % getResult(res2, STD, -7_longInt)

#line 186 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 186) )
  if (anyExceptions()) return
#line 187 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 187 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 187) )
  if (anyExceptions()) return
#line 188 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 188 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    call mem % getResult(res1, 8_longInt)
    call mem % getResult(res2, STD, 8_longInt)

#line 193 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, res1, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 194 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, res2, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 194) )
  if (anyExceptions()) return
#line 195 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 195 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, STD, TOL, &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 195) )
  if (anyExceptions()) return
#line 196 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    ! Free memor y
    call mem % kill()

  end subroutine testScoring

  !!
  !! Test lastCycle
  !! Ignors test parametrisation
  !!
!@Test(cases=[1])
  subroutine testLastCycle(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem
    integer(shortInt)                      :: i

    call mem % init(1_longInt, 1, batchSize = 8)

    ! Test getting batchSize
#line 215 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(8, mem % getBatchSize(),'Test getBatchSize() :', &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 215) )
  if (anyExceptions()) return
#line 216 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

    do i=1,16
      if(i == 8 .or. i == 16) then
#line 219 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertTrue( mem % lastCycle(), 'In cycle num: '//numToChar(i), &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 219) )
  if (anyExceptions()) return
#line 220 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
      else
#line 221 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertFalse( mem % lastCycle(), 'In cycle num: '//numToChar(i), &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
      end if
      call mem % closeCycle(ONE)
    end do

    call mem % kill()

  end subroutine testLastCycle

  !!
  !! Test get score
  !! Ignore test parametrisation
  !!
!@Test(cases=[1])
  subroutine testGetScore(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem
    real(defReal),parameter :: TOL = 1.0E-9

    call mem % init(1_longInt, 1)

    call mem % score(ONE,1_longInt)
    call mem % score(ONE,1_longInt)
    call mem % score(ONE,1_longInt)

#line 246 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(3*ONE, mem % getScore(1_longInt), TOL, 'Test getScore, valid bin:', &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 246) )
  if (anyExceptions()) return
#line 247 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 247 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, mem % getScore(0_longInt), TOL, 'Test getScore, not +ve bin:', &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 247) )
  if (anyExceptions()) return
#line 248 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
#line 248 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"
  call assertEqual(ZERO, mem % getScore(2_longInt), TOL, 'Test getScore, too large bin:', &
 & location=SourceLocation( &
 & 'scoreMemory_test.f90', &
 & 248) )
  if (anyExceptions()) return
#line 249 "/home/mhgk4/SCONE/Tallies/Tests/scoreMemory_test.f90"

  end subroutine testGetScore

  !!
  !! Test killing uninitialised scoreMemory
  !!
!@Test(cases=[1])
  subroutine testKillUnalloc(this)
    class(test_scoreMemory), intent(inout) :: this
    type(scoreMemory)                      :: mem

    call mem % kill()

  end subroutine testKillUnalloc

end module scoreMemory_test
!! MATLAB SCRIPT USED TO GENERATE REFERENCE VALUES
!clear
!rand = zeros(20,1);
!seed = 9294;
!
!%LCG Params
!A = 2469;
!M = 65521;
!
!for i=1:1:200
!  seed = mod(A * seed, M);
!  rand(i) = seed/M;
!end
!
!% Calculate scores vector
!scores = 2.0 + sin(pi() .* rand - pi()/2);
!scoresInt = floor(100.*rand);
!
!% Accumulate results
!resAcc = mean(scores)
!stdAcc = sqrt(var(scores)./200)
!
!resAccInt = mean(scoresInt)
!stdAccInt = sqrt(var(scoresInt)./200)
!
!% Reshape scores
!scores = reshape(scores,[20,10]);
!scores = sum(scores,1)* 0.7;
!res = mean(scores)
!std = sqrt(var(scores)./10)
!
!% Reshape scores
!scoresInt = reshape(scoresInt,[20,10]);
!scoresInt = sum(scoresInt,1)* 0.7;
!resInt = mean(scoresInt)
!stdInt = sqrt(var(scoresInt)./10)

module WrapscoreMemory_test
   use pFUnit_mod
   use scoreMemory_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_scoreMemory) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use scoreMemory_test
        class (test_scoreMemory), intent(inout) :: this
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
      type (testNumber), intent(in) :: testParameter
      aTest%test_scoreMemory = newTest(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WrapscoreMemory_test

function scoreMemory_test_suite() result(suite)
   use pFUnit_mod
   use scoreMemory_test
   use WrapscoreMemory_test
   type (TestSuite) :: suite

   type (testNumber), allocatable :: testParameters(:)
   type (testNumber) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('scoreMemory_test_suite')

   cases = [1]
   testParameters = [(new_testNumber(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testScoring', testScoring, testParameter))
   end do

   cases = [1]
   testParameters = [(new_testNumber(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testLastCycle', testLastCycle, testParameter))
   end do

   cases = [1]
   testParameters = [(new_testNumber(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testGetScore', testGetScore, testParameter))
   end do

   cases = [1]
   testParameters = [(new_testNumber(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testKillUnalloc', testKillUnalloc, testParameter))
   end do


end function scoreMemory_test_suite

