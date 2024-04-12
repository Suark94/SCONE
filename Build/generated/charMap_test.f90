module charMap_test
  use numPrecision
  use charMap_class, only : charMap
  use pFUnit_mod

  implicit none


!@testCase
  type, extends(TestCase) :: test_charMap
    private
    type(charMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_charMap

  ! Parameters
  integer(shortInt), parameter :: VAL1 = 0
  integer(shortInt), parameter :: VAL2 = 1
  integer(shortInt), parameter :: VAL3 = -1
  integer(shortInt), parameter :: VAL4 =  huge(VAL1)
  integer(shortInt), parameter :: VAL5 = -huge(VAL1)
  integer(shortInt), parameter :: VAL6 =  133316666 ! Number of deamons of hell
  integer(shortInt), parameter :: VAL7 = -133316666

  character(nameLen), parameter :: KEY1 = 'Since you are reluctant '
  character(nameLen), parameter :: KEY2 = 'to provide us with the'
  character(nameLen), parameter :: KEY3 = 'location of the rebel base,'
  character(nameLen), parameter :: KEY4 = 'I have chosen to test this'
  character(nameLen), parameter :: KEY5 = 'station destructive power'
  character(nameLen), parameter :: KEY6 = 'on your home planet of'
  character(nameLen), parameter :: KEY7 = 'Alderaan'

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_charMap), intent(inout) :: this
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
    class(test_charMap), intent(inout) :: this

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
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    temp = this % map % get(KEY1)
#line 80 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL1, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY2)
#line 83 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL2, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY3)
#line 86 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL3, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY4)
#line 89 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL4, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY5)
#line 92 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL5, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY6)
#line 95 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY7)
#line 98 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL7, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testRetrieval

  !!
  !! Test Retrivial with deletions
  !!
!@Test
  subroutine testRetrievalWithDel(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    ! Delate Elements
    call this % map % del(KEY1)
    call this % map % del(KEY2)
    call this % map % del(KEY3)
    call this % map % del(KEY4)
    call this % map % del(KEY5)

    ! Obtain some elements
    temp = this % map % get(KEY6)
#line 119 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % get(KEY7)
#line 122 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL7, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testRetrievalWithDel

  !!
  !! Test getting entery with default for absent keys
  !!
!@Test
  subroutine testGetOrDefault(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp
    integer(shortInt)                  :: default = 8
    character(nameLen)                 :: tempChar

    ! Google says its Chinese Lucky number. We do need luck
    default = 8

    ! Retrieve present entries. Do all cases to spot possible
    ! small differences between get and getOrDefault implementations
    temp = this % map % getOrDefault(KEY1, default)
#line 142 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL1, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % getOrDefault(KEY2, default)
#line 145 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL2, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % getOrDefault(KEY3, default)
#line 148 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL3, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % getOrDefault(KEY4, default)
#line 151 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL4, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % getOrDefault(KEY5, default)
#line 154 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL5, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % getOrDefault(KEY6, default)
#line 157 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    temp = this % map % getOrDefault(KEY7, default)
#line 160 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL7, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    ! Get default.
    tempChar = 'this key is not present'
#line 164 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(default, this % map % getOrDefault(tempChar, default), &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testGetOrDefault

  !!
  !! Test deletions
  !!
!@Test
  subroutine testDel(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    ! Before deletion
    temp = this % map % get(KEY6)
#line 178 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL6, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    call this % map % del(KEY6)

    ! After deletion
    temp = this % map % getOrDefault(KEY6, VAL1)
#line 184 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(VAL1, temp, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testDel

  !!
  !! Test Looping over map
  !!
!@Test
  subroutine testLooping(this)
    class(test_charMap), intent(inout)        :: this
    integer(shortInt),dimension(6),parameter  :: VALS = [VAL1, VAL2, VAL3, VAL4, VAL5, VAL7]
    character(nameLen),dimension(6),parameter :: KEYS = [KEY1, KEY2, KEY3, KEY4, KEY5, KEY7]
    character(nameLen),dimension(6)           :: KEYS_PAST
    integer(shortInt)                         :: counter, it, tVal
    character(nameLen)                        :: tKey

    ! Initialise parameters
    KEYS_PAST = "This is not a Key. It's picture of a key!"
    counter = 0

    ! Delete Entry 6 from map
    call this % map % del(KEY6)

    ! Loop over remaining elements
    it = this % map % begin()
    do while( it /= this % map % end() )
      tVal = this % map % atVal(it)
      tKey = this % map % atKey(it)

#line 213 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertTrue(any(VALS == tVal),"Wrong Value", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 213) )
  if (anyExceptions()) return
#line 214 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
#line 214 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertTrue(any(KEYS == tKey),"Wrong Key", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 214) )
  if (anyExceptions()) return
#line 215 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

      ! Change value and read it
      call this % map % atSet(7, it)
#line 218 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(7, this % map % atVal(it), "Ups. Wrong Value after change", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

      ! Make shure that KEY is not getting repeated
#line 221 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertFalse(any(KEYS_PAST(1:counter) == tKey),"Repeated KEY", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
      counter = counter + 1
      KEYS_PAST(counter) = tKey

      it = this % map % next(it)
    end do

#line 228 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(6, counter, &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 228) )
  if (anyExceptions()) return
#line 229 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testLooping

  !!
  !! Test Looping Edge Cases
  !!
!@Test
  subroutine testLoopingEdgeCases(this)
    class(test_charMap), intent(inout) :: this
    type(charMap)                      :: locMap
    integer(shortInt)                  :: it

    ! Loop over uninitialised map
    it = locMap % begin()
    do while(it /= locMap % end())
#line 244 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertTrue(.false.,"Should not enter the loop", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 244) )
  if (anyExceptions()) return
#line 245 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
      it = locMap % next(it)
    end do

    ! Loop over empty map
    call locMap % init(8)
    it = locMap % begin()
    do while(it /= locMap % end())
#line 252 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertTrue(.false.,"Should not enter the loop", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 252) )
  if (anyExceptions()) return
#line 253 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
      it = locMap % next(it)
    end do

    ! Loop over map with deleted elements
    call locMap % add(KEY1,3)
    call locMap % add(KEY2,3)
    call locMap % add(KEY3,3)

    call locMap % del(KEY1)
    call locMap % del(KEY2)
    call locMap % del(KEY3)

    it = locMap % begin()
    !print *, it, locMap % atKey(it), locMap % map(:) % status
    do while(it /= locMap % end())
#line 268 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertTrue(.false.,"Should not enter the loop", &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 268) )
  if (anyExceptions()) return
#line 269 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
      it = locMap % next(it)
    end do

  end subroutine testLoopingEdgeCases

  !!
  !! Test putting new entry under exoisting key
  !!
!@Test
  subroutine testOverwriting(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    temp = 7

    ! Store over existing keyword
    call this % map % add(KEY4, temp)

    ! Verify correctness
#line 288 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(temp, this % map % get(KEY4), &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 288) )
  if (anyExceptions()) return
#line 289 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testOverwriting

  !!
  !! Test getting length
  !!
!@Test
  subroutine testGetLength(this)
    class(test_charMap), intent(inout) :: this

#line 299 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(7, this % map % length(), &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

    ! Deleate some elements
    call this % map % del(KEY1)
    call this % map % del(KEY6)

#line 305 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(5, this % map % length(), &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 305) )
  if (anyExceptions()) return
#line 306 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"


  end subroutine testGetLength

  !!
  !! Test getOrDefaoult from uninitialised map
  !!
!@Test
  subroutine testGetOrDefaultUninitialised(this)
    class(test_charMap), intent(inout) :: this
    character(nameLen)                 :: key
    type(charMap)                      :: locMap


    key = 'Key of admiral Dodanna'
#line 321 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"
  call assertEqual(7, locMap % getOrDefault(key, 7), &
 & location=SourceLocation( &
 & 'charMap_test.f90', &
 & 321) )
  if (anyExceptions()) return
#line 322 "/home/mhgk4/SCONE/DataStructures/Tests/charMap_test.f90"

  end subroutine testGetOrDefaultUninitialised

end module charMap_test

module WrapcharMap_test
   use pFUnit_mod
   use charMap_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_charMap) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use charMap_test
        class (test_charMap), intent(inout) :: this
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

end module WrapcharMap_test

function charMap_test_suite() result(suite)
   use pFUnit_mod
   use charMap_test
   use WrapcharMap_test
   type (TestSuite) :: suite

   suite = newTestSuite('charMap_test_suite')

   call suite%addTest(makeCustomTest('testRetrieval', testRetrieval))

   call suite%addTest(makeCustomTest('testRetrievalWithDel', testRetrievalWithDel))

   call suite%addTest(makeCustomTest('testGetOrDefault', testGetOrDefault))

   call suite%addTest(makeCustomTest('testDel', testDel))

   call suite%addTest(makeCustomTest('testLooping', testLooping))

   call suite%addTest(makeCustomTest('testLoopingEdgeCases', testLoopingEdgeCases))

   call suite%addTest(makeCustomTest('testOverwriting', testOverwriting))

   call suite%addTest(makeCustomTest('testGetLength', testGetLength))

   call suite%addTest(makeCustomTest('testGetOrDefaultUninitialised', testGetOrDefaultUninitialised))


end function charMap_test_suite

