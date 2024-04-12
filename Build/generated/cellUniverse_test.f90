module cellUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables, only : UNDEF_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use cellUniverse_class, only : cellUniverse
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = &
  " surf1 { id 1; type sphere; origin (0.0 0.0 0.0); radius 2;}&
  & surf2 { id 2; type sphere; origin (4.0 0.0 0.0); radius 1;}"

  character(*), parameter :: CELL_DEF = &
  " cell1 {id 1; type simpleCell; surfaces (-1); filltype uni; universe 3;} &
  & cell2 {id 2; type simpleCell; surfaces (1 2); filltype uni; universe 4;}"


  !
  ! Note that rotation is such that following axis transformation applies:
  !   x -> z
  !   y -> -y
  !   z -> x
  !
  character(*), parameter :: UNI_DEF = &
  "id 1; type cellUniverse; origin (2.0 0.0 0.0); rotation (90.0 90.0 90.0); cells (1 2);"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(cellUniverse) :: uni

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

    call charToDict(dict, CELL_DEF)
    call cells % init(dict, surfs, mats)
    call dict % kill()

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)
    call dict % kill()

    ! Set index
    call uni % setIdx(8)

    ! Verify fill
#line 71 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual([-3, -4, UNDEF_MAT], fill, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

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
#line 96 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(1, uni % id(), &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

    ! Set ID
    call uni % setId(7)
#line 100 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(7, uni % id(), &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 100) )
  if (anyExceptions()) return
#line 101 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
!@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref, r, dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** Enter into local cell 1
    r = [0.0_defReal, 0.0_defReal, 3.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    u_ref = [ONE, ZERO, ZERO]
#line 122 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 123 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 124 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 124) )
  if (anyExceptions()) return
#line 125 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 125 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(1, new % localID, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 125) )
  if (anyExceptions()) return
#line 126 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 126 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(cells % getIdx(1), new % cellIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

    ! ** Enter into local cell 2
    r = [2.0_defReal, 0.0_defReal, 1.0_defReal]
    dir = [ZERO, ONE, ZERO]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = [-1.0_defReal, 0.0_defReal, 2.0_defReal]
    u_ref = [ZERO, -ONE, ZERO]
#line 137 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 138 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 139 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 140 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(2, new % localID, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 141 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(cells % getIdx(2), new % cellIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

    ! ** Enter into the UNDEFINED cell
    r = [0.0_defReal, 0.0_defReal, 6.5_defReal]
    dir = [ONE, ZERO, ZERO]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = [4.5_defReal, 0.0_defReal, 0.0_defReal]
    u_ref = [ZERO, ZERO, ONE]
#line 152 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 153 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 154 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 155 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(3, new % localID, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 156 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

    ! Verify rotation settings in coord
    ! * Do it only once
#line 160 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertTrue(new % isRotated, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 161 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual([ZERO, ZERO,  ONE], new % rotMat(1,:), TOL, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 162 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual([ZERO, -ONE, ZERO], new % rotMat(2,:), TOL, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 163 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual([ONE , ZERO, ZERO], new % rotMat(3,:), TOL, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"


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

    ! ** In local cell 1 distance to boundary
    pos % r = [-1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 8
    pos % cellIdx = cells % getIdx(1)
    pos % localId = 1

    call uni % distance(d, surfIdx, pos)

    ref = 3.0_defReal
#line 188 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 189 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(surfs % getIdx(1), surfIdx , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"


    ! ** In local cell 2 distance to surface 2
    pos % r = [7.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [-ONE, ZERO, ZERO]
    pos % cellIdx = cells % getIdx(2)
    pos % localId = 2

    call uni % distance(d, surfIdx, pos)

    ref = 2.0_defReal
#line 201 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 202 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(surfs % getIdx(2), surfIdx , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

    ! ** In local cell 2 distance to infinity
    ! surfIdx must be set to 0
    pos % dir = [ONE, ZERO, ZERO]
    call uni % distance(d, surfIdx, pos)

#line 209 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(INF, d, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 210 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(0, surfIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 210) )
  if (anyExceptions()) return
#line 211 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
!@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! Cross from cell 1 to cell 2
    pos % r   = [0.0_defReal, 2.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = cells % getIdx(1)
    pos % localId = 1

    idx = surfs % getIdx(1)
    call uni % cross(pos, idx)

#line 232 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(2, pos % localId, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 232) )
  if (anyExceptions()) return
#line 233 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
#line 233 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual(cells % getIdx(2), pos % cellIdx, &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 233) )
  if (anyExceptions()) return
#line 234 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
!@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Cell 1
    pos % r   = [0.0_defReal, 1.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = cells % getIdx(1)
    pos % localId = 1

#line 251 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 251) )
  if (anyExceptions()) return
#line 252 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

    ! Cell 2
    pos % r   = [-7.0_defReal, 2.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = cells % getIdx(2)
    pos % localId = 2

#line 260 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) , &
 & location=SourceLocation( &
 & 'cellUniverse_test.f90', &
 & 260) )
  if (anyExceptions()) return
#line 261 "/home/mhgk4/SCONE/Geometry/Universes/Tests/cellUniverse_test.f90"

  end subroutine test_cellOffset



end module cellUniverse_test

module WrapcellUniverse_test
   use pFUnit_mod
   use cellUniverse_test
   implicit none
   private

contains


end module WrapcellUniverse_test

function cellUniverse_test_suite() result(suite)
   use pFUnit_mod
   use cellUniverse_test
   use WrapcellUniverse_test
   type (TestSuite) :: suite

   suite = newTestSuite('cellUniverse_test_suite')

   call suite%addTest(newTestMethod('test_misc', test_misc, setUp, clean))

   call suite%addTest(newTestMethod('test_enter', test_enter, setUp, clean))

   call suite%addTest(newTestMethod('test_distance', test_distance, setUp, clean))

   call suite%addTest(newTestMethod('test_cross', test_cross, setUp, clean))

   call suite%addTest(newTestMethod('test_cellOffset', test_cellOffset, setUp, clean))


end function cellUniverse_test_suite

