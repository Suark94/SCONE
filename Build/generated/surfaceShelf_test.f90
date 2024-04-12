module surfaceShelf_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use surface_inter,      only : surface
  use surfaceShelf_class, only : surfaceShelf
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SHELF_DEF = &
    " surf1 {id 18; type sphere; origin (1.0 1.0 1.0); radius 17.4;} &
    & surf2 {id 23; type xPlane; x0 3.1;} &
    & surf3 {id 1; type zPlane; z0  0.3;} "

  ! Variables
  type(surfaceShelf) :: shelf


contains

  !!
  !! Build the Shelf
  !!
!@Before
  subroutine setUp()
    type(dictionary) :: dict

    ! Try killing uninitialised
    call shelf % kill()

    call charToDict(dict, SHELF_DEF)
    call shelf % init(dict)

  end subroutine setUp

  !!
  !! Clean after SHelf
  !!
!@After
  subroutine cleanUp()

    call shelf % kill()

  end subroutine cleanUp

  !!
  !! Test
  !!
!@Test
  subroutine testShelf()
    class(surface), pointer :: ptr
    integer(shortInt)       :: idx

    ! 1st Surface
    idx = shelf % getIdx(18)
    ptr => shelf % getPtr(idx)
#line 60 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual(18, ptr % id(), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
#line 61 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual('sphere', ptr % myType(), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
#line 62 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual(18, shelf % getId(idx), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"

    ! 2nd Surface
    idx = shelf % getIdx(1)
    ptr => shelf % getPtr(idx)
#line 67 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual(1, ptr % id(), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
#line 68 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual('zPlane', ptr % myType(), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
#line 69 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual(1, shelf % getId(idx), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"

    ! Test size
#line 72 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"
  call assertEqual(3, shelf % getSize(), &
 & location=SourceLocation( &
 & 'surfaceShelf_test.f90', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/surfaceShelf_test.f90"

  end subroutine testShelf


end module surfaceShelf_test

module WrapsurfaceShelf_test
   use pFUnit_mod
   use surfaceShelf_test
   implicit none
   private

contains


end module WrapsurfaceShelf_test

function surfaceShelf_test_suite() result(suite)
   use pFUnit_mod
   use surfaceShelf_test
   use WrapsurfaceShelf_test
   type (TestSuite) :: suite

   suite = newTestSuite('surfaceShelf_test_suite')

   call suite%addTest(newTestMethod('testShelf', testShelf, setUp, cleanUp))


end function surfaceShelf_test_suite

