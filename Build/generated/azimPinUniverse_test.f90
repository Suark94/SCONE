module azimPinUniverse_test

  use numPrecision
  use universalVariables,    only : INF, SURF_TOL
  use dictionary_class,      only : dictionary
  use dictParser_func,       only : charToDict
  use charMap_class,         only : charMap
  use coord_class,           only : coord
  use surfaceShelf_class,    only : surfaceShelf
  use cellShelf_class,       only : cellShelf
  use azimPinUniverse_class, only : azimPinUniverse, MOVING_IN, MOVING_OUT, MOVING_CLOCK, &
                                    MOVING_ANTI, MOVING_CLOCK_FORWARD, MOVING_CLOCK_BACK
  use pfUnit_mod
  implicit none

  ! Parameters
  character(*), parameter :: UNI_DEF1 = &
  "id 7; type azimPinUniverse; naz 4; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (2.5 1.5 0.0); fills (u<7> u<14> void);"
  character(*), parameter :: UNI_DEF2 = &
  "id 8; type azimPinUniverse; naz 8; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
  &radii (4.0 2.5 1.5 0.0); fills (u<20> u<7> u<14> void);"

  ! Variables
  type(surfaceShelf)    :: surfs
  type(cellShelf)       :: cells
  type(charMap)         :: mats
  type(azimPinUniverse) :: uni1, uni2


contains

  !!
  !! Set-up test enviroment
  !!
!@Before
  subroutine setup()
    character(nameLen)                           :: name
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)                             :: dict
    integer(shortInt), dimension(:), allocatable :: fillArray

    ! Load void material
    name = 'void'
    call mats % add(name, 13)

    ! Build universe
    call charToDict(dict, UNI_DEF1)
    call uni1 % init(fill, dict, cells, surfs, mats)
    
    ! Set index
    call uni1 % setIdx(3)
    
    ! Verify fill array
#line 55 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual([-14, -14, -14, -14, -7, -7, -7, -7, 13, 13, 13, 13], fill, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
    
    ! Build second universe
    fillArray = [-14, -14, -14, -14, -14, -14, -14, -14,-7, -7, -7, -7, -7, -7, -7, -7, &
            -20, -20, -20, -20, -20, -20, -20, -20, 13, 13, 13, 13, 13, 13, 13, 13]
    call charToDict(dict, UNI_DEF2)
    call uni2 % init(fill, dict, cells, surfs, mats)
    call uni2 % setIdx(26)
#line 63 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(fillArray, fill, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

  end subroutine setup

  !!
  !! Clean after test
  !!
!@After
  subroutine clean()

    call surfs % kill()
    call cells % kill()
    call mats % kill()
    call uni1 % kill()
    call uni2 % kill()

  end subroutine clean

  !!
  !! Test miscellaneous functionality
  !!
!@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! Get id
#line 89 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(7, uni1 % id(), &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! Set ID
    call uni1 % setId(7)
#line 93 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(7, uni1 % id(), &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 93) )
  if (anyExceptions()) return
#line 94 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

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
    r = [1.0_defReal, 0.0_defReal, 0.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 115 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 116 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 116) )
  if (anyExceptions()) return
#line 117 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 117 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 118 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(1, new % localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 119 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
    
    ! ** Enter into local cell 2
    r = [0.0_defReal, 1.0_defReal, 0.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 130 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 131 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 132 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 133 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(2, new % localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 133) )
  if (anyExceptions()) return
#line 134 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 134 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! ** Enter into local cell 8
    r = [0.0_defReal, -2.3_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 145 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 146 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 147 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 147) )
  if (anyExceptions()) return
#line 148 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 148 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(8, new % localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 149 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! ** Enter into local cell 11
    r = [-2.6_defReal, 0.0_defReal, -980.0_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    ! Verify location
    r_ref = r
    u_ref = dir
#line 160 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 161 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 162 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 163 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(11, new % localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 164 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! VERIFY THAT ROTATION IS NOT SET (all angles were 0.0)
#line 167 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertFalse(new % isRotated, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 167) )
  if (anyExceptions()) return
#line 168 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! In universe 2, enter local cell 21
    r = [-3.0_defReal, 0.1_defReal, 32.0_defReal]
    call uni2 % enter(new, r, dir)
    
    ! Verify location
    r_ref = r
    u_ref = dir
#line 176 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 176) )
  if (anyExceptions()) return
#line 177 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 177 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 178 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(26, new % uniIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 179 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(21, new % localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 180 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
!@Test
  subroutine test_distance()
    real(defReal)     :: d, ref
    integer(shortInt) :: surfIdx
    type(coord)       :: pos
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! ** In local cell 1 distance to radial boundary
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % distance(d, surfIdx, pos)

    ref = 0.5_defReal
#line 204 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 204) )
  if (anyExceptions()) return
#line 205 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 205 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_OUT, surfIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 205) )
  if (anyExceptions()) return
#line 206 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! ** In local cell 1 distance to anti-clockwise boundary
    pos % r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [-SQRT2_2, SQRT2_2, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2_2
#line 217 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 218 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_ANTI, surfIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! ** In local cell 3 distance to clockwise boundary
    pos % r = [-1.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [SQRT2_2, SQRT2_2, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 3

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2_2
#line 230 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 231 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_CLOCK, surfIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 231) )
  if (anyExceptions()) return
#line 232 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! ** In local cell 4 distance to clockwise boundary
    ! ** Moves back around!
    pos % r = [0.0_defReal, -1.0_defReal, 0.0_defReal]
    pos % dir = [SQRT2_2, SQRT2_2, ZERO]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 4

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2_2
#line 244 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 244) )
  if (anyExceptions()) return
#line 245 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 245 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_CLOCK_FORWARD, surfIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! ** In outermost cell moving away
    pos % r = [2.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localId = 9

    call uni1 % distance(d, surfIdx, pos)
#line 253 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(INF, d, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 253) )
  if (anyExceptions()) return
#line 254 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
    ! Surface momento is undefined -> No crossing

    ! In ordinary cell in-between
    pos % r = [0.0_defReal, 1.6_defReal, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]
    pos % localId = 5

    call uni1 % distance(d, surfIdx, pos)
    ref = 0.1_defReal
#line 263 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ref, d, ref * tol, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 263) )
  if (anyExceptions()) return
#line 264 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 264 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_IN, surfIdx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 264) )
  if (anyExceptions()) return
#line 265 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! Universe 2


  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
!@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx
    real(defReal) :: eps

    ! Cross from cell 2 to cell 6
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 8
    pos % cellIdx = 0
    pos % localId = 2

    idx = MOVING_OUT
    call uni1 % cross(pos, idx)

#line 291 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(6, pos % localId, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! Cross from cell 6 to cell 2
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal+eps, 0.0_defReal]
    pos % dir = [ZERO, -ONE, ZERO]

    idx = MOVING_IN
    call uni1 % cross(pos, idx)

#line 301 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(2, pos % localId, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 301) )
  if (anyExceptions()) return
#line 302 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! Cross from cell 1 to cell 4
    eps = HALF * SURF_TOL
    pos % r   = [SQRT2_2, -SQRT2_2, ZERO]
    pos % dir = [-SQRT2_2, -SQRT2_2, ZERO]
    pos % r   = pos % r - eps * pos % dir
    pos % localID = 1

    idx = MOVING_CLOCK_BACK
    call uni1 % cross(pos, idx)

#line 313 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(4, pos % localId, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 313) )
  if (anyExceptions()) return
#line 314 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! Universe 2
    ! 

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
!@Test
  subroutine test_cellOffset()
    type(coord)       :: pos

    ! Cell 2
    pos % r   = [0.0_defReal, 1.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 2

#line 334 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni1 % cellOffset(pos) , &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 334) )
  if (anyExceptions()) return
#line 335 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! Cell 11
    pos % r   = [-7.0_defReal, 0.0_defReal, 0.0_defReal]
    pos % dir = [ZERO, ONE, ZERO]
    pos % uniIdx = 3
    pos % cellIdx = 0
    pos % localId = 11

#line 343 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual([ZERO, ZERO, ZERO], uni1 % cellOffset(pos) , &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 343) )
  if (anyExceptions()) return
#line 344 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

  end subroutine test_cellOffset

  !!
  !! Test surface transitions
  !!
  !! Check that there is no problem with distance calculations
  !! if particle is placed very close to an annulus or plane surface 
  !! (within SURF_TOL)
  !!
!@Test
  subroutine test_edgeCases()
    type(coord)       :: pos
    integer(shortInt) :: idx, localID, cellIdx
    real(defReal)     :: eps, d
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! At boundary between cell 2 and 6
    eps = HALF * SURF_TOL
    pos % r   = [0.0_defReal, 1.5_defReal-eps, 0.0_defReal]
    pos % dir = [ONE, -0.00001_defReal, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx = 8
    pos % cellIdx = 0

    ! Should find particle in cell 2
    ! And return very small distance -> MOVING OUT
    call uni1 % findCell(localID, cellIdx, pos % r, pos % dir)
#line 372 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(2, localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 372) )
  if (anyExceptions()) return
#line 373 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    pos % localID = 2
    call uni1 % distance(d, idx, pos)

#line 377 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ZERO, d, 1.0E-3_defReal, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 377) )
  if (anyExceptions()) return
#line 378 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 378 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_OUT, idx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 378) )
  if (anyExceptions()) return
#line 379 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    ! At boundary between cell 4 and 1
    eps = 1.1*SURF_TOL
    pos % r   = [SQRT2_2, -SQRT2_2, ZERO]
    pos % dir = [-SQRT2_2, -SQRT2_2, ZERO]
    pos % r   = pos % r + eps * pos % dir
    pos % dir = -pos % dir
    
    ! Should find particle in cell 4
    ! And return very small distance -> MOVING_CLOCK_FORWARD
    call uni1 % findCell(localID, cellIDx, pos % r, pos % dir)
#line 390 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(4, localID, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 390) )
  if (anyExceptions()) return
#line 391 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"

    pos % localID = 4
    call uni1 % distance(d, idx, pos)

#line 395 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(ZERO, d, 1.0E-3_defReal, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 395) )
  if (anyExceptions()) return
#line 396 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
#line 396 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
  call assertEqual(MOVING_CLOCK_FORWARD, idx, &
 & location=SourceLocation( &
 & 'azimPinUniverse_test.f90', &
 & 396) )
  if (anyExceptions()) return
#line 397 "/home/mhgk4/SCONE/Geometry/Universes/Tests/azimPinUniverse_test.f90"
    
  end subroutine test_edgeCases


end module azimPinUniverse_test

module WrapazimPinUniverse_test
   use pFUnit_mod
   use azimPinUniverse_test
   implicit none
   private

contains


end module WrapazimPinUniverse_test

function azimPinUniverse_test_suite() result(suite)
   use pFUnit_mod
   use azimPinUniverse_test
   use WrapazimPinUniverse_test
   type (TestSuite) :: suite

   suite = newTestSuite('azimPinUniverse_test_suite')

   call suite%addTest(newTestMethod('test_misc', test_misc, setup, clean))

   call suite%addTest(newTestMethod('test_enter', test_enter, setup, clean))

   call suite%addTest(newTestMethod('test_distance', test_distance, setup, clean))

   call suite%addTest(newTestMethod('test_cross', test_cross, setup, clean))

   call suite%addTest(newTestMethod('test_cellOffset', test_cellOffset, setup, clean))

   call suite%addTest(newTestMethod('test_edgeCases', test_edgeCases, setup, clean))


end function azimPinUniverse_test_suite

