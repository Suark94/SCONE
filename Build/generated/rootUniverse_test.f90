module rootUniverse_test

  use numPrecision
  use universalVariables, only : OUTSIDE_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use rootUniverse_class, only : rootUniverse
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = &
  "surf1 { id 4; type sphere; origin (0.0 5.0 0.0); radius 0.5;}&
  &surf2 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}"

  character(*), parameter :: UNI_DEF = &
  "id 1; type rootUniverse; border 1; fill u<17>;"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(rootUniverse) :: uni

contains

  !!
  !! Setup environment
  !!
!@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary) :: dict

    ! Build surfaces and MATS

    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)
    call dict % kill()

    ! Set index
    call uni % setIdx(8)

    ! Verify fill
#line 55 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual([-17, OUTSIDE_MAT], fill, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

  end subroutine setUp

  !!
  !! Clean environment
  !!
!@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality (of generic universe)
  !!
!@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! Get id
#line 80 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(1, uni % id(), &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

    ! Set ID
    call uni % setId(7)
#line 84 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(7, uni % id(), &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

    ! Test boundary surface
#line 87 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(surfs % getIdx(1), uni % border(), &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
!@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r ,dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Enter inside
    r = [1.0_defReal, -1.0_defReal, 1.0_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni % enter(new, r, dir)

#line 106 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(r,   new % r, TOL, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 107 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(dir, new % dir, TOL, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 108 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 109 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 110 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(1, new % localID, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

    ! Enter outside
    r = [2.0_defReal, -2.0_defReal, 1.0_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni % enter(new, r, dir)

#line 118 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(r,   new % r, TOL, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 119 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(dir, new % dir, TOL, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 120 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 121 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 122 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(2, new % localID, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
!@Test
  subroutine test_distance()
    real(defReal)     :: d, ref
    integer(shortInt) :: surfIdx
    type(coord)       :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Distance from inside -> only relevant
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 1

    call uni % distance(d, surfIdx, pos)

    ref = 1.0_defReal
#line 145 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(ref, d, ref * TOL, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
#line 146 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(surfs % getIdx(1), surfIdx, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
!@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! Cross into outside
    pos % r = [2.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 1

    idx = surfs % getIdx(1)
    call uni % cross(pos, idx)

#line 167 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual(2, pos % localID, &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 167) )
  if (anyExceptions()) return
#line 168 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
!@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Inside
    pos % r = [1.5_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 1

#line 184 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos), &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

    ! Outside
    pos % r = [2.5_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % localID = 2

#line 192 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos), &
 & location=SourceLocation( &
 & 'rootUniverse_test.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/Geometry/Universes/Tests/rootUniverse_test.f90"

  end subroutine test_cellOffset



end module rootUniverse_test

module WraprootUniverse_test
   use pFUnit_mod
   use rootUniverse_test
   implicit none
   private

contains


end module WraprootUniverse_test

function rootUniverse_test_suite() result(suite)
   use pFUnit_mod
   use rootUniverse_test
   use WraprootUniverse_test
   type (TestSuite) :: suite

   suite = newTestSuite('rootUniverse_test_suite')

   call suite%addTest(newTestMethod('test_misc', test_misc, setUp, clean))

   call suite%addTest(newTestMethod('test_enter', test_enter, setUp, clean))

   call suite%addTest(newTestMethod('test_distance', test_distance, setUp, clean))

   call suite%addTest(newTestMethod('test_cross', test_cross, setUp, clean))

   call suite%addTest(newTestMethod('test_cellOffset', test_cellOffset, setUp, clean))


end function rootUniverse_test_suite

