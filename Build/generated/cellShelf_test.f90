module cellShelf_test
  use numPrecision
  use universalVariables, only : OUTSIDE_MAT, VOID_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use cell_inter,         only : cell
  use pFUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SURFS_DEF = "&
  & surf1 {id 7; type sphere; origin (0.0 0.0 0.0); radius 2.0;} &
  & surf2 {id 8; type zPlane; z0 -1.3;} &
  & surf3 {id 9; type zPlane; z0 0.0;} &
  & surf4 {id 10; type zPlane; z0 1.0;}"

  character(*), parameter :: CELLS_DEF = &
  " cell1 {id 4; type simpleCell; surfaces (7); filltype outside;} &
  & cell2 {id 2; type simpleCell; surfaces (-7 10); filltype uni; universe 8;} &
  & cell3 {id 1; type simpleCell; surfaces (-7 -10 9); filltype mat; material void;} &
  & cell4 {id 5; type simpleCell; surfaces (-7 -9 8); filltype mat; material fuel;} &
  & cell5 {id 8; type simpleCell; surfaces (-7 -8); filltype mat; material pecorino;}"


  ! Variables
  type(charMap)      :: mats
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells

contains

  !!
  !! Build the cell
  !!
!@Before
  subroutine setUp()
    type(dictionary) :: dict
    character(nameLen) :: name

    ! Surface shelf
    call charToDict(dict, SURFS_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Material map
    name = 'void'
    call mats % add(name, VOID_MAT)
    name = 'fuel'
    call mats % add(name, 1)
    name = 'pecorino'
    call mats % add(name, 2)

    ! Build cells
    call charToDict(dict, CELLS_DEF)
    call cells % init(dict, surfs, mats)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
!@After
  subroutine cleanUp()

    call surfs % kill()
    call cells % kill()
    call mats % kill()

  end subroutine cleanUp

  !!
  !! Test shelf
  !!
!@Test
  subroutine test_get()
    integer(shortInt) :: idx
    class(cell), pointer :: ptr

    ! Cell ID 1
    idx = cells % getIdx(1)
    ptr => cells % getPtr(idx)
#line 86 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(1, ptr % id(), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
#line 87 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(1, cells % getID(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

    ! Cell ID 5
    idx = cells % getIdx(5)
    ptr => cells % getPtr(idx)
#line 92 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(5, ptr % id(), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
#line 93 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(5, cells % getID(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 93) )
  if (anyExceptions()) return
#line 94 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

    ! Cell ID 8
    idx = cells % getIdx(8)
    ptr => cells % getPtr(idx)
#line 98 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(8, ptr % id(), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
#line 99 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(8, cells % getID(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

    ! Test size
#line 102 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(5, cells % getSize(), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

  end subroutine test_get

  !!
  !! Test cell filling
  !!
!@Test
  subroutine test_cell_fill()
    integer(shortInt) :: idx

    ! Cell ID 4
    idx = cells % getIdx(4)
#line 115 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(OUTSIDE_MAT, cells % getFill(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

    ! Cell ID 2 -> universe
    idx = cells % getIdx(2)
#line 119 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(-8, cells % getFill(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

    ! Cell ID 1 -> void
    idx = cells % getIdx(1)
#line 123 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(VOID_MAT, cells % getFill(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

    ! Cell ID 8 -> Pecorino fill
    idx = cells % getIdx(8)
#line 127 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"
  call assertEqual(2, cells % getFill(idx), &
 & location=SourceLocation( &
 & 'cellShelf_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/Geometry/Cells/Tests/cellShelf_test.f90"

  end subroutine test_cell_fill

end module cellShelf_test

module WrapcellShelf_test
   use pFUnit_mod
   use cellShelf_test
   implicit none
   private

contains


end module WrapcellShelf_test

function cellShelf_test_suite() result(suite)
   use pFUnit_mod
   use cellShelf_test
   use WrapcellShelf_test
   type (TestSuite) :: suite

   suite = newTestSuite('cellShelf_test_suite')

   call suite%addTest(newTestMethod('test_get', test_get, setUp, cleanUp))

   call suite%addTest(newTestMethod('test_cell_fill', test_cell_fill, setUp, cleanUp))


end function cellShelf_test_suite

