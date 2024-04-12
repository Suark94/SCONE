module universe_test
  use numPrecision
  use charMap_class,  only : charMap
  use universe_inter, only : charToFill
  use pFUnit_mod
  implicit none


  !!
  !! Note that universe is abstract thus it cannot be tested by itself
  !!
  !! Tests for universe non-overridable procedures are in cellUniverse_test
  !!

contains

  !!
  !! Test charToFill
  !!
  !! Since universe is abstract it cannot be tested by itself (only via its
  !! subclasses)
  !!
!@Test
  subroutine test_charToFill()
    type(charMap)      :: mats
    character(nameLen) :: name
    character(100), parameter :: Here = 'parentScope'

    ! Load some material names and their (fake) matIdxs
    name = 'mat13'
    call mats % add(name, 13)

    name = 'city17'
    call mats % add(name, 17)

    name = 'mat47'
    call mats % add(name, 47)

    ! Test material conversion
    name = 'mat13'
#line 41 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"
  call assertEqual(13, charToFill(name, mats, Here), &
 & location=SourceLocation( &
 & 'universe_test.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"

    name = 'city17'
#line 44 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"
  call assertEqual(17, charToFill(name, mats, Here), &
 & location=SourceLocation( &
 & 'universe_test.f90', &
 & 44) )
  if (anyExceptions()) return
#line 45 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"

    ! Test universe ID conversion
    ! By convention returns -ve uniID
    ! NO SPACES IN THE CHAR!!!
    name = 'u<87>'
#line 50 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"
  call assertEqual(-87, charToFill(name, mats, Here), &
 & location=SourceLocation( &
 & 'universe_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"

    name = 'u<133>'
#line 53 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"
  call assertEqual(-133, charToFill(name, mats, Here), &
 & location=SourceLocation( &
 & 'universe_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"

    name = 'u<5>'
#line 56 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"
  call assertEqual(-5, charToFill(name, mats, Here), &
 & location=SourceLocation( &
 & 'universe_test.f90', &
 & 56) )
  if (anyExceptions()) return
#line 57 "/home/mhgk4/SCONE/Geometry/Universes/Tests/universe_test.f90"

  end subroutine test_charToFill



end module universe_test

module Wrapuniverse_test
   use pFUnit_mod
   use universe_test
   implicit none
   private

contains


end module Wrapuniverse_test

function universe_test_suite() result(suite)
   use pFUnit_mod
   use universe_test
   use Wrapuniverse_test
   type (TestSuite) :: suite

   suite = newTestSuite('universe_test_suite')

   call suite%addTest(newTestMethod('test_charToFill', test_charToFill))


end function universe_test_suite

