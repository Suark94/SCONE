module pinUniverse_test

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use pinUniverse_class,  only : pinUniverse, MOVING_IN, MOVING_OUT
  use pfUnit_mod
  implicit none

  ! Parameters
  character(*), parameter :: UNI_DEF = &
  "id 7; type pinUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (2.5 1.5 0.0); fills (u<7> u<14> void);"

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(pinUniverse)  :: uni


contains

  !!
  !! Set-up test environment
  !!
!@Before
  subroutine setup()
    character(nameLen)                           :: name
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)                             :: dict

    ! Load void material
    name = 'void'
    call mats % add(name, 13)

    ! Build universe
    call charToDict(dict, UNI_DEF)
    call uni % init(fill, dict, cells, surfs, mats)

    ! Set index
    call uni % setIdx(3)

    ! Verify fill array
#line 50 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual([-14, -7, 13], fill, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"


  end subroutine setup

  !!
  !! Clean after test
  !!
!@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality
  !!
!@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! Get id
#line 76 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(7, uni % id(), &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! Set ID
    call uni % setId(7)
#line 80 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(7, uni % id(), &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

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
    r = [0.0_defReal, 1.0_defReal, 0.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 102 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 103 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 103) )
  if (anyExceptions()) return
#line 104 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 104 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 105 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(1, new % localID, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 106 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! ** Enter into local cell 2
    r = [2.3_defReal, 0.0_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 117 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 118 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 119 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 120 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(2, new % localID, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 121 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! ** Enter into local cell 3
    r = [2.6_defReal, 0.0_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 132 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 133 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 133) )
  if (anyExceptions()) return
#line 134 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 134 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 135 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(3, new % localID, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 135) )
  if (anyExceptions()) return
#line 136 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 136 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! VERIFY THAT ROTATION IS NOT SET (all angles were 0.0)
#line 139 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertFalse(new % isRotated, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

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
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni % distance(d, surfIdx, pos)

    ref = 0.5_defReal
#line 163 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 164 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(MOVING_OUT, surfIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! ** In outermost cell moving away
    pos % r = [2.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 3

    call uni % distance(d, surfIdx, pos)
#line 172 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(INF, d, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 172) )
  if (anyExceptions()) return
#line 173 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
    ! Surface momento is undefined -> No crossing

    ! In ordinary cell in-between
    pos % r = [0.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]
    pos % localId = 2

    call uni % distance(d, surfIdx, pos)
    ref = 0.1_defReal
#line 182 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 182) )
  if (anyExceptions()) return
#line 183 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 183 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(MOVING_IN, surfIdx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
!@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx
    real(defReal) :: eps

    ! Cross from cell 1 to cell 2
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = 0
    pos % localId = 1

    idx = MOVING_OUT
    call uni % cross(pos, idx)

#line 207 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(2, pos % localId, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 207) )
  if (anyExceptions()) return
#line 208 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! Cross form cell 2 to cell 1
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal+eps, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]

    idx = MOVING_IN
    call uni % cross(pos, idx)

#line 217 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(1, pos % localId, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

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
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 1

#line 235 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) , &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    ! Cell 3
    pos % r   = [-7.0_defReal, 2.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 3

#line 244 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni % cellOffset(pos) , &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 244) )
  if (anyExceptions()) return
#line 245 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

  end subroutine test_cellOffset

  !!
  !! Test surface transitions
  !!
  !! Check that there is no problem with distance calculations
  !! if particle is placed very close to an annulus surface (within SURF_TOL)
  !!
!@Test
  subroutine test_edgeCases()
    type(coord)       :: pos
    integer(shortInt) :: idx, localID, cellIdx
    real(defReal)     :: eps, d
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! At boundary between cell 1 and 2
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ONE, -0.00001_defReal, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx = 8
    pos % cellIdx = 0

    ! Should find particle in cell 1
    ! And return very small distance -> MOVING OUT
    call uni % findCell(localID, cellIDx, pos % r, pos % dir)
#line 272 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(1, localID, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 272) )
  if (anyExceptions()) return
#line 273 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

    pos % localID = 1
    call uni % distance(d, idx, pos)

#line 277 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(ZERO, d, 1.0E-3_defReal, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 277) )
  if (anyExceptions()) return
#line 278 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
#line 278 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"
  call assertEqual(MOVING_OUT, idx, &
 & location=SourceLocation( &
 & 'pinUniverse_test.f90', &
 & 278) )
  if (anyExceptions()) return
#line 279 "/home/mhgk4/SCONE/Geometry/Universes/Tests/pinUniverse_test.f90"

  end subroutine test_edgeCases


end module pinUniverse_test

module WrappinUniverse_test
   use pFUnit_mod
   use pinUniverse_test
   implicit none
   private

contains


end module WrappinUniverse_test

function pinUniverse_test_suite() result(suite)
   use pFUnit_mod
   use pinUniverse_test
   use WrappinUniverse_test
   type (TestSuite) :: suite

   suite = newTestSuite('pinUniverse_test_suite')

   call suite%addTest(newTestMethod('test_misc', test_misc, setup, clean))

   call suite%addTest(newTestMethod('test_enter', test_enter, setup, clean))

   call suite%addTest(newTestMethod('test_distance', test_distance, setup, clean))

   call suite%addTest(newTestMethod('test_cross', test_cross, setup, clean))

   call suite%addTest(newTestMethod('test_cellOffset', test_cellOffset, setup, clean))

   call suite%addTest(newTestMethod('test_edgeCases', test_edgeCases, setup, clean))


end function pinUniverse_test_suite

