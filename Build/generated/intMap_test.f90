module intMap_test
  use numPrecision
  use intMap_class, only : intMap
  use pFUnit_mod

  implicit none


!@testCase
  type, extends(TestCase) :: test_intMap
    private
    type(intMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_intMap

  ! Parameters
  integer(shortInt), parameter :: VAL1 = 9
  integer(shortInt), parameter :: VAL2 = 1
  integer(shortInt), parameter :: VAL3 = -1
  integer(shortInt), parameter :: VAL4 =  huge(VAL1)
  integer(shortInt), parameter :: VAL5 = -huge(VAL1)
  integer(shortInt), parameter :: VAL6 =  133316666 ! Number of deamons of hell
  integer(shortInt), parameter :: VAL7 = -133316666

  integer(shortInt), parameter :: KEY1 = 0
  integer(shortInt), parameter :: KEY2 = 1
  integer(shortInt), parameter :: KEY3 = -1
  integer(shortInt), parameter :: KEY4 =  huge(VAL1)
  integer(shortInt), parameter :: KEY5 = -huge(VAL1)
  integer(shortInt), parameter :: KEY6 =  133316666 ! Number of deamons of hell
  integer(shortInt), parameter :: KEY7 = -133316666




contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp
    integer(shortInt)                 :: N, i

    ! Load entries
     call this % map % add(KEY1, VAL1)
     call this % map % add(KEY2, VAL2)
     call this % map % add(KEY3, VAL3)
     call this % map % add(KEY4, VAL4)
     call this % map % add(KEY5, VAL5)
     call this % map % add(KEY6, VAL6)
     call this % map % add(KEY7, VAL7)

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_intMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving values from the map.
  !! All values are present
  !!
!@Test
  subroutine testRetrieval(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp

    temp = this % map % get(KEY1)
#line 83 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL1, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY2)
#line 86 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL2, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY3)
#line 89 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL3, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY4)
#line 92 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL4, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY5)
#line 95 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL5, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY6)
#line 98 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY7)
#line 101 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL7, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testRetrieval

  !!
  !! Test Retrivial with deletions
  !!
!@Test
  subroutine testRetrievalWithDel(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp

    ! Delate Elements
    call this % map % del(KEY1)
    call this % map % del(KEY2)
    call this % map % del(KEY3)
    call this % map % del(KEY4)
    call this % map % del(KEY5)

    ! Obtain some elements
    temp = this % map % get(KEY6)
#line 122 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % get(KEY7)
#line 125 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL7, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 125) )
  if (anyExceptions()) return
#line 126 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testRetrievalWithDel

  !!
  !! Test getting entery with default for absent keys
  !!
!@Test
  subroutine testGetOrDefault(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp
    integer(shortInt)                 :: default = 8

    ! Google says its Chinese Lucky number. We do need luck
    default = 8

    ! Retrieve present entries. Do all cases to spot possible
    ! small differences between get and getOrDefault implementations
    temp = this % map % getOrDefault(KEY1, default)
#line 144 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL1, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 144) )
  if (anyExceptions()) return
#line 145 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % getOrDefault(KEY2, default)
#line 147 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL2, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 147) )
  if (anyExceptions()) return
#line 148 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % getOrDefault(KEY3, default)
#line 150 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL3, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 150) )
  if (anyExceptions()) return
#line 151 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % getOrDefault(KEY4, default)
#line 153 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL4, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % getOrDefault(KEY5, default)
#line 156 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL5, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % getOrDefault(KEY6, default)
#line 159 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    temp = this % map % getOrDefault(KEY7, default)
#line 162 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL7, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    ! Get default. (Yes I love Tolkien - MAK)
#line 165 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(default, this % map % getOrDefault(1973,default), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testGetOrDefault

  !!
  !! Test deletions
  !!
!@Test
  subroutine testDel(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp

    ! Before deletion
    temp = this % map % get(KEY6)
#line 179 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    call this % map % del(KEY6)

    ! After deletion
    temp = this % map % getOrDefault(KEY6, VAL1)
#line 185 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(VAL1, temp, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testDel

  !!
  !! Test Looping over map
  !!
!@Test
  subroutine testLooping(this)
    class(test_intMap), intent(inout)        :: this
    integer(shortInt),dimension(6),parameter :: VALS = [VAL1, VAL2, VAL3, VAL4, VAL5, VAL7]
    integer(shortInt),dimension(6),parameter :: KEYS = [KEY1, KEY2, KEY3, KEY4, KEY5, KEY7]
    integer(shortInt),dimension(6)           :: KEYS_PAST
    integer(shortInt)                        :: counter, it, tVal, tKey

    ! Initialise parameters
    KEYS_PAST = 7
    counter = 0

    ! Delete Entry 6 from map
    call this % map % del(KEY6)

    ! Loop over remaining elements
    it = this % map % begin()
    do while( it /= this % map % end() )
      tVal = this % map % atVal(it)
      tKey = this % map % atKey(it)

#line 213 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertTrue(any(VALS == tVal),"Wrong Value", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 213) )
  if (anyExceptions()) return
#line 214 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
#line 214 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertTrue(any(KEYS == tKey),"Wrong Key", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 214) )
  if (anyExceptions()) return
#line 215 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

      ! Change value and read it
      call this % map % atSet(7, it)
#line 218 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(7, this % map % atVal(it), "Ups. Wrong Value after change", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

      ! Make shure that KEY is not getting repeated
#line 221 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertFalse(any(KEYS_PAST(1:counter) == tKey),"Repeated KEY", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
      counter = counter + 1
      KEYS_PAST(counter) = tKey

      it = this % map % next(it)
    end do

#line 228 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(6, counter, &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 228) )
  if (anyExceptions()) return
#line 229 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testLooping

  !!
  !! Test Looping Edge Cases
  !!
!@Test
  subroutine testLoopingEdgeCases(this)
    class(test_intMap), intent(inout) :: this
    type(intMap)                      :: locMap
    integer(shortInt)                 :: it

    ! Loop over uninitialised map
    it = locMap % begin()
    do while(it /= locMap % end())
#line 244 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertTrue(.false.,"Should not enter the loop", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 244) )
  if (anyExceptions()) return
#line 245 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
      it = locMap % next(it)
    end do

    ! Loop over empty map
    call locMap % init(8)
    it = locMap % begin()
    do while(it /= locMap % end())
#line 252 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertTrue(.false.,"Should not enter the loop", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 252) )
  if (anyExceptions()) return
#line 253 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
      it = locMap % next(it)
    end do

    ! Loop over map with deleted elements
    call locMap % add(1,3)
    call locMap % add(2,3)
    call locMap % add(7,3)

    call locMap % del(1)
    call locMap % del(2)
    call locMap % del(7)

    it = locMap % begin()
    do while(it /= locMap % end())
#line 267 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertTrue(.false.,"Should not enter the loop", &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 267) )
  if (anyExceptions()) return
#line 268 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
      it = locMap % next(it)
    end do

  end subroutine testLoopingEdgeCases

  !!
  !! Test putting new entry under exoisting key
  !!
!@Test
  subroutine testOverwriting(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp

    temp = 7

    ! Store over existing keyword
    call this % map % add(KEY4, temp)

    ! Verify correctness
#line 287 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(temp, this % map % get(KEY4), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 287) )
  if (anyExceptions()) return
#line 288 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testOverwriting

  !!
  !! Test getting length
  !!
!@Test
  subroutine testGetLength(this)
    class(test_intMap), intent(inout) :: this

#line 298 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(7, this % map % length(), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 298) )
  if (anyExceptions()) return
#line 299 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    ! Deleate some elements
    call this % map % del(KEY1)
    call this % map % del(KEY6)

#line 304 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(5, this % map % length(), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 304) )
  if (anyExceptions()) return
#line 305 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testGetLength

  !!
  !! Test getOrDefaoult from uninitialised map
  !!
!@Test
  subroutine testGetOrDefaultUninitialised(this)
    class(test_intMap), intent(inout) :: this
    type(intMap)                      :: locMap

#line 316 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(7, locMap % getOrDefault(3, 7), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 316) )
  if (anyExceptions()) return
#line 317 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

  end subroutine testGetOrDefaultUninitialised

  !!
  !! Test getOrDefault from initialised but empty map
  !!
!@Test
  subroutine testGetOrDefaultEmpty(this)
    class(test_intMap), intent(inout) :: this
    type(intMap)                      :: locMap

    ! Pure empty map
    call locMap % init(2)
#line 330 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(7, locMap % getOrDefault(3, 7), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 330) )
  if (anyExceptions()) return
#line 331 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"

    ! Map with deleted element
    call locMap % add(3, 7)
    call locMap % del(3)
    call locMap % init(2)
#line 336 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"
  call assertEqual(2, locMap % getOrDefault(3, 2), &
 & location=SourceLocation( &
 & 'intMap_test.f90', &
 & 336) )
  if (anyExceptions()) return
#line 337 "/home/mhgk4/SCONE/DataStructures/Tests/intMap_test.f90"


  end subroutine testGetOrDefaultEmpty


end module intMap_test

module WrapintMap_test
   use pFUnit_mod
   use intMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_intMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use intMap_test
        class (test_intMap), intent(inout) :: this
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

end module WrapintMap_test

function intMap_test_suite() result(suite)
   use pFUnit_mod
   use intMap_test
   use WrapintMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('intMap_test_suite')

   call suite%addTest(makeCustomTest('testRetrieval', testRetrieval))

   call suite%addTest(makeCustomTest('testRetrievalWithDel', testRetrievalWithDel))

   call suite%addTest(makeCustomTest('testGetOrDefault', testGetOrDefault))

   call suite%addTest(makeCustomTest('testDel', testDel))

   call suite%addTest(makeCustomTest('testLooping', testLooping))

   call suite%addTest(makeCustomTest('testLoopingEdgeCases', testLoopingEdgeCases))

   call suite%addTest(makeCustomTest('testOverwriting', testOverwriting))

   call suite%addTest(makeCustomTest('testGetLength', testGetLength))

   call suite%addTest(makeCustomTest('testGetOrDefaultUninitialised', testGetOrDefaultUninitialised))

   call suite%addTest(makeCustomTest('testGetOrDefaultEmpty', testGetOrDefaultEmpty))


end function intMap_test_suite

