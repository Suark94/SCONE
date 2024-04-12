module hashFunctions_test

  use numPrecision
  use hashFunctions_func, only : FNV_1, knuthHash
  use pfUnit_mod

  implicit none

contains

  !!
  !! Test computation of FNV_1 hashes on short and long integer
  !!
!@Test
  subroutine testFNV1()
    integer(shortInt) :: hashVal1, hashResult1

    ! Test for short int
    call FNV_1('All hail gfortran!', hashVal1)
    hashResult1 = 64933584_shortInt !int(z'03deced0', shortInt)
#line 21 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertEqual(hashResult1, hashVal1, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 21) )
  if (anyExceptions()) return
#line 22 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"

    call FNV_1('MY HASH IS BROKEN', hashVal1)
    hashResult1 = -2038433774_shortInt !int(z'867ff812', shortInt)
#line 25 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertEqual(hashResult1, hashVal1, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 25) )
  if (anyExceptions()) return
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"

  end subroutine testFNV1

  !!
  !! Test knuth hash
  !!
!@Test
  subroutine testKnuthHash()
    integer(shortInt), dimension(100) :: keys
    integer(shortInt), dimension(100) :: hashes
    integer(shortInt)                 :: i, m

    ! Generate some random keys of shortInts
    keys(1) = 765875
    do i=2,100
      keys(i) = keys(i-1) * 22695477
    end do

    ! Case 1 range 2**m m = 4
    m = 4
    hashes = knuthHash(keys, m)
#line 47 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertLessThanOrEqual(0, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
#line 48 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertGreaterThan(2**m, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"

    ! Case 2 range m = 30
    m = 30
    hashes = knuthHash(keys, m)
#line 53 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertLessThanOrEqual(0, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
#line 54 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertGreaterThan(2**m, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"

    ! Case 3 range m = 37 (will be clipped to m = 31)
    m = 37
    hashes = knuthHash(keys, m)
#line 59 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertLessThanOrEqual(0, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 59) )
  if (anyExceptions()) return
#line 60 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
#line 60 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertGreaterThan(huge(m), hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"

    ! Case 4 -ve m (will be clipped to m=1
    m = -71
    hashes = knuthHash(keys, m)
#line 65 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertLessThanOrEqual(0, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
#line 66 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"
  call assertGreaterThan(2, hashes, &
 & location=SourceLocation( &
 & 'hashFunctions_test.f90', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/home/mhgk4/SCONE/SharedModules/Tests/hashFunctions_test.f90"

  end subroutine testKnuthHash


end module hashFunctions_test

module WraphashFunctions_test
   use pFUnit_mod
   use hashFunctions_test
   implicit none
   private

contains


end module WraphashFunctions_test

function hashFunctions_test_suite() result(suite)
   use pFUnit_mod
   use hashFunctions_test
   use WraphashFunctions_test
   type (TestSuite) :: suite

   suite = newTestSuite('hashFunctions_test_suite')

   call suite%addTest(newTestMethod('testFNV1', testFNV1))

   call suite%addTest(newTestMethod('testKnuthHash', testKnuthHash))


end function hashFunctions_test_suite

