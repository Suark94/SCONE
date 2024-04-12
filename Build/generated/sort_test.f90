module sort_test
  use numPrecision
  use genericProcedures, only : swap, quickSort
  use pfUnit_mod

  implicit none


contains

  !!
  !! Test swaping
  !!
!@Test
  subroutine testSwaping()
    real(defReal)     :: r1, r2, d1, d2
    integer(shortInt) :: i1, i2

    ! Swap diffrent integers
    i1 = -7
    i2 =  4

    call swap(i1, i2)

#line 25 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(4_shortInt, i1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 25) )
  if (anyExceptions()) return
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(-7_shortInt, i2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
    
    ! Swap same integers (If one tries XOR swap this can fail)
    i1 = 4
    i2 = 4
    call swap(i1, i2)
#line 32 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(4_shortInt, i1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 32) )
  if (anyExceptions()) return
#line 33 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 33 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(4_shortInt, i2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Swap reals
    r1 = -7.3_defReal
    r2 = 2.3_defReal

    call swap(r1, r2)

#line 41 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(2.3_defReal, r1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 42 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(-7.3_defReal, r2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 42) )
  if (anyExceptions()) return
#line 43 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Swap same reals
    r1 = 2.3_defReal
    r2 = 2.3_defReal

    call swap(r1, r2)

#line 50 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(2.3_defReal, r1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 51 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(2.3_defReal, r2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 51) )
  if (anyExceptions()) return
#line 52 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Swap pair of reals
    r1 = 2.3_defReal
    r2 = 0.3_defReal
    d1 = 1.7_defReal
    d2 = 0.4_defReal

    call swap(r1, r2, d1, d2)

#line 61 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(1.7_defReal, r1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 62 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(0.4_defReal, r2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 63 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(2.3_defReal, d1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 64 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(0.3_defReal, d2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 64) )
  if (anyExceptions()) return
#line 65 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Swap same pair of reals
    r1 = 2.3_defReal
    r2 = 0.3_defReal
    d1 = 2.3_defReal
    d2 = 0.3_defReal

    call swap(r1, r2, d1, d2)

#line 74 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(2.3_defReal, r1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 75 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(0.3_defReal, r2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 75) )
  if (anyExceptions()) return
#line 76 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 76 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(2.3_defReal, d1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 77 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(0.3_defReal, d2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

  end subroutine testSwaping

  !!
  !! Test sorting integer array
  !!
!@Test
  subroutine testQuickSortInt()
    integer(shortInt), dimension(4) :: I = [2, 1, 4, 3]
    integer(shortInt), dimension(:), allocatable :: IA

    ! Ordinary sort
    call quickSort(I)
#line 91 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1, 2, 3, 4], I, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 91) )
  if (anyExceptions()) return
#line 92 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Sorted array
    I = [1, 2, 3, 4]
    call quickSort(I)
#line 96 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1, 2, 3, 4], I, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Uniform Array
    I = [1, 1, 1, 1]
    call quickSort(I)
#line 101 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1, 1, 1, 1], I, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! 1 Element array
    IA = [1]
    call quickSort(IA)
#line 106 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1],IA, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Empty array
    deallocate(IA)
    allocate(IA(0))

    call quickSort(IA)
#line 113 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(0, size(IA), &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

  end subroutine testQuickSortInt

  !!
  !! Test sorting real array
  !!
!@Test
  subroutine testQuickSortReal()
    real(defReal), dimension(4) :: R = [2.1_defReal, 1.4_defReal, 4.9_defReal, 3.4_defReal]
    real(defReal), dimension(:), allocatable :: RA

    ! Ordinary sort
    call quickSort(R)
#line 127 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Sorted array
    R = [1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal]
    call quickSort(R)
#line 132 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Uniform Array
    R = [1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal]
    call quickSort(R)
#line 137 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal], R, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! 1 Element array
    RA = [1.0_defReal]
    call quickSort(RA)
#line 142 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.0_defReal],RA, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Empty array
    deallocate(RA)
    allocate(RA(0))

    call quickSort(RA)
#line 149 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual(0, size(RA), &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

  end subroutine testQuickSortReal

  !!
  !! Test sorting two arrays of reals by keys in first array
  !!
!@Test
  subroutine testQuickSortRealReal()
    real(defReal), dimension(4) :: R1 = [2.1_defReal, 1.4_defReal, 4.9_defReal, 3.4_defReal]
    real(defReal), dimension(4) :: R2 = [1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal]
    real(defReal), dimension(3) :: R = ZERO

    ! Ordinary sort
    call quickSort(R1, R2)
#line 164 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 165 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([2.0_defReal, 1.0_defReal, 4.0_defReal, 3.0_defReal], R2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Sorted array
    R1 = [1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal]
    R2 = [1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal]

    call quickSort(R1, R2)
#line 172 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.4_defReal, 2.1_defReal, 3.4_defReal, 4.9_defReal], R1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 172) )
  if (anyExceptions()) return
#line 173 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 173 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal], R2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 173) )
  if (anyExceptions()) return
#line 174 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

    ! Uniform array
    R1 = [1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal]
    R2 = [1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal]

    call quickSort(R1, R2)
#line 180 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.0_defReal, 1.0_defReal, 1.0_defReal, 1.0_defReal], R1, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
#line 181 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"
  call assertEqual([1.0_defReal, 2.0_defReal, 3.0_defReal, 4.0_defReal], R2, &
 & location=SourceLocation( &
 & 'sort_test.f90', &
 & 181) )
  if (anyExceptions()) return
#line 182 "/home/mhgk4/SCONE/SharedModules/Tests/sort_test.f90"

  end subroutine testQuickSortRealReal


end module sort_test

module Wrapsort_test
   use pFUnit_mod
   use sort_test
   implicit none
   private

contains


end module Wrapsort_test

function sort_test_suite() result(suite)
   use pFUnit_mod
   use sort_test
   use Wrapsort_test
   type (TestSuite) :: suite

   suite = newTestSuite('sort_test_suite')

   call suite%addTest(newTestMethod('testSwaping', testSwaping))

   call suite%addTest(newTestMethod('testQuickSortInt', testQuickSortInt))

   call suite%addTest(newTestMethod('testQuickSortReal', testQuickSortReal))

   call suite%addTest(newTestMethod('testQuickSortRealReal', testQuickSortRealReal))


end function sort_test_suite

