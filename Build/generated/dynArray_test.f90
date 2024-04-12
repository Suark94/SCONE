module dynArray_test
  use numPrecision
  use dynArray_class, only : dynIntArray
  use pFUnit_mod

  implicit none

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour on unallocated dynArray
  !!
!@Test
  subroutine testUnallocatedInt()
    type(dynIntArray) :: array

#line 21 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(0, array % getSize(),'Size of unallocated array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 21) )
  if (anyExceptions()) return
#line 22 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 22 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(0, array % capacity(), 'Capacity of unallocated array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 22) )
  if (anyExceptions()) return
#line 23 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 23 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertTrue(array % isEmpty(),'isEmpty on unallocated array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 23) )
  if (anyExceptions()) return
#line 24 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    ! Make shure no SEG ERR happens
    call array % shrink()

  end subroutine testUnallocatedInt

  !!
  !! Test correct behaviour on empty dynArray
  !!
!@Test
  subroutine testEmptyInt()
    type(dynIntArray) :: array

    call array % add([1,2,3,4,5])
    call array % empty()

#line 40 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(0, array % getSize(),'Size of empty array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 41 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertTrue(array % isEmpty(),'isEmpty on empty array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    ! Make shure no SEG ERR happens
    call array % shrink()
#line 45 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(0, array % capacity(),'Capacity of shrunk empty array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

  end subroutine testEmptyInt

  !!
  !! Test resizing
  !!
!@Test
  subroutine testResizeInt()
    type(dynIntArray) :: array

    call array % resize(2)
#line 57 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertLessThanOrEqual(2, array % capacity(),'Resize to 2 from 0', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    call array % resize(5)
#line 60 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertLessThanOrEqual(5, array % capacity(),'Resize to 5 from 2', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    call array % resize(20)
#line 63 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertLessThanOrEqual(20, array % capacity(),'Resize to 20 from 5 ', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    call array % kill()
#line 66 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(0, array % capacity(),'Capacity of a killed array', &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

  end subroutine testResizeInt

  !!
  !! Test usage
  !!
!@Test
  subroutine testUsage()
    type(dynIntArray) :: array
    integer(shortInt) :: i

    ! Build from vector by assignment
    array = [1,2,0,4,2]

    ! Test getting elements by index
#line 82 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(1, array % get(1), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 83 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(2, array % get(2), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 84 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(0, array % get(3), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 85 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(4, array % get(4), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 86 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(2, array % get(5), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    call array % empty()

    ! Build by elements
    do i=1,7
      call array % add(2 * i)
    end do

    ! Test isEmpty and getSize
#line 96 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertFalse(array % isEmpty(), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 97 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(7, array % getSize(), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 97) )
  if (anyExceptions()) return
#line 98 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

    ! Test popping
    do i=7,1,-1
#line 101 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(2*i, array % pop(), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
    end do

   ! Test building by emelent and vector
   call array % add(1)
   call array % add([8,3])

   ! Test expose
   associate (ar => array % expose())
#line 110 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(1, ar(1), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 111 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(8, ar(2), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
#line 112 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(3, ar(3), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
   end associate

   ! Shrink non-empty array
   call array % shrink()
#line 117 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"
  call assertEqual(3, array % capacity(), &
 & location=SourceLocation( &
 & 'dynArray_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/DataStructures/Tests/dynArray_test.f90"

  end subroutine testUsage

  !!
  !! Test kill of unallocated array
  !!
!@Test
  subroutine testKillUnalloc()
    type(dynIntArray) :: array

    call array % kill()

  end subroutine testKillUnalloc

end module dynArray_test

module WrapdynArray_test
   use pFUnit_mod
   use dynArray_test
   implicit none
   private

contains


end module WrapdynArray_test

function dynArray_test_suite() result(suite)
   use pFUnit_mod
   use dynArray_test
   use WrapdynArray_test
   type (TestSuite) :: suite

   suite = newTestSuite('dynArray_test_suite')

   call suite%addTest(newTestMethod('testUnallocatedInt', testUnallocatedInt))

   call suite%addTest(newTestMethod('testEmptyInt', testEmptyInt))

   call suite%addTest(newTestMethod('testResizeInt', testResizeInt))

   call suite%addTest(newTestMethod('testUsage', testUsage))

   call suite%addTest(newTestMethod('testKillUnalloc', testKillUnalloc))


end function dynArray_test_suite

