module latUniverse_test

  use numPrecision
  use genericProcedures
  use universalVariables, only : UNDEF_MAT
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use coord_class,        only : coord
  use surfaceShelf_class, only : surfaceShelf
  use cellShelf_class,    only : cellShelf
  use latUniverse_class,  only : latUniverse
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: UNI1_DEF = &
  "id 1; type latUniverse; origin (0.0 0.0 0.0); rotation (0.0 0.0 0.0); &
   pitch (1.0 2.0 3.0); shape (3 2 2); padMat void; &
  &map ( 3 4 5 &
  &      7 4 8 &
  &            &
  &      1 2 3 &
  &      4 5 6); "

  character(*), parameter :: UNI2_DEF = &
  "id 2; type latUniverse; pitch (1.0 2.0 0.0); shape (2 1 0); padMat u<1>; &
  &map (1 2); "

  ! Variables
  type(surfaceShelf) :: surfs
  type(cellShelf)    :: cells
  type(charMap)      :: mats
  type(latUniverse)  :: uni1
  type(latUniverse)  :: uni2

contains

  !!
  !! Setup environment
  !!
!@Before
  subroutine setUp()
    integer(shortInt), dimension(:), allocatable :: fill
    type(dictionary)   :: dict
    character(nameLen) :: name
    integer(shortInt), dimension(:), allocatable :: ref

    ! Add materials
    name = 'void'
    call mats % add(name, 3)

    ! Build universe 1
    call charToDict(dict, UNI1_DEF)
    call uni1 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni1 % setIdx(8)

    ! Verify fill vector
    ref = [-4, -5, -6, -1, -2, -3, -7, -4, -8, -3, -4, -5, 3]
#line 62 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, fill, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Build universe 2
    call charToDict(dict, UNI2_DEF)
    call uni2 % init(fill, dict, cells, surfs, mats)
    call dict % kill()
    call uni2 % setIdx(3)

    ref = [-1, -2, -1]
#line 71 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, fill, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

  end subroutine setUp

  !!
  !! Clean environment
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
  !! Test miscellaneous functionality (of generic universe)
  !!
!@Test
  subroutine test_misc()
    real(defReal), dimension(3,3) :: mat

    ! * Single universe is fine here
    ! Get id
#line 98 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(1, uni1 % id(), &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Set ID
    call uni1 % setId(7)
#line 102 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(7, uni1 % id(), &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

  end subroutine test_misc

  !!
  !! Test entering a universe
  !!
!@Test
  subroutine test_enter()
    type(coord) :: new
    real(defReal), dimension(3) :: r_ref, u_ref, r ,dir
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Enter inside -> Away from surface
    r = [1.0_defReal, 1.0_defReal, 0.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
#line 124 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 124) )
  if (anyExceptions()) return
#line 125 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 125 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 125) )
  if (anyExceptions()) return
#line 126 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 126 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 127 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(12, new % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 128 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Enter outside
    r = [1.6_defReal, 0.5_defReal, 0.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
#line 138 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 139 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 140 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 141 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(13, new % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 142 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Enter in a corner
    r = [-0.5_defReal, 0.0_defReal, 0.0_defReal ]
    dir = [-ONE, ONE, -ONE]
    dir = dir / norm2(dir)

    call uni1 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
#line 153 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 154 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 155 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(8, new % uniIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 156 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(4, new % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 157 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! ** 2D Universe
    ! Enter inside -> Away from surface
    r = [0.5_defReal, 0.5_defReal, 13.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni2 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
#line 168 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 168) )
  if (anyExceptions()) return
#line 169 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 169 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 169) )
  if (anyExceptions()) return
#line 170 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 170 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 171 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(2, new % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 172 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 172) )
  if (anyExceptions()) return
#line 173 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Enter outside
    r = [1.6_defReal, 0.5_defReal, 0.5_defReal ]
    dir = [ZERO, ZERO, ONE]

    call uni2 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
#line 182 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 182) )
  if (anyExceptions()) return
#line 183 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 183 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 184 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 185 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(3, new % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 186 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 186) )
  if (anyExceptions()) return
#line 187 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Enter on a face
    r = [0.0_defReal, 0.0_defReal, 0.0_defReal ]
    dir = [-ONE, ONE, -ONE]
    dir = dir / norm2(dir)

    call uni2 % enter(new, r, dir)

    r_ref = r
    u_ref = dir
#line 197 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(r_ref, new % r, TOL , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 197) )
  if (anyExceptions()) return
#line 198 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 198 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(u_ref, new % dir, TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 198) )
  if (anyExceptions()) return
#line 199 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 199 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(3, new % uniIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 199) )
  if (anyExceptions()) return
#line 200 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 200 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(1, new % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 200) )
  if (anyExceptions()) return
#line 201 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 201 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(0, new % cellIdx, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"


  end subroutine test_enter

  !!
  !! Test distance calculation
  !!
!@Test
  subroutine test_distance()
    real(defReal)     :: d, ref, eps
    integer(shortInt) :: surfIdx
    type(coord)       :: pos
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D universe
    ! Well inside a cell
    pos % r = [0.0_defReal, 0.1_defReal, 0.5_defReal]
    pos % dir = [ZERO, -ZERO, ONE]
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 11

    call uni1 % distance(d, surfIdx, pos)

    ref = 2.5_defReal
#line 227 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 227) )
  if (anyExceptions()) return
#line 228 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 228 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(-6, surfIdx , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 228) )
  if (anyExceptions()) return
#line 229 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! From outside -> miss
    pos % r = [-4.0_defReal, 0.1_defReal, 0.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 13

    call uni1 % distance(d, surfIdx, pos)

#line 240 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(INF, d, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 241 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(-7, surfIdx , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! After a surface undershoot
    eps = HALF * SURF_TOL
    pos % r = [-1.0_defReal, 0.0_defReal-eps, -0.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 4

    call uni1 % distance(d, surfIdx, pos)

    ref = SQRT2 * HALF
#line 255 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, d, ref * TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 255) )
  if (anyExceptions()) return
#line 256 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 256 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(-2, surfIdx , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 256) )
  if (anyExceptions()) return
#line 257 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! After overshoot via a corner
    pos % r = [-0.5_defReal+eps, 0.0_defReal+eps, -0.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 4

    call uni1 % distance(d, surfIdx, pos)
#line 267 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ZERO, d,  TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 267) )
  if (anyExceptions()) return
#line 268 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 268 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(-2, surfIdx , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 268) )
  if (anyExceptions()) return
#line 269 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    !** 2D universe
    ! Well inside a cell -> Vertical
    pos % r = [0.5_defReal, 0.6_defReal, 0.5_defReal]
    pos % dir = [ZERO, ZERO, ONE]
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 2

    call uni2 % distance(d, surfIdx, pos)

#line 280 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(INF, d, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 280) )
  if (anyExceptions()) return
#line 281 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Well inside a cell -> Shallow hit
    pos % r = [0.5_defReal, 0.6_defReal, 0.5_defReal]
    pos % dir = [ZERO, 0.01_defReal, ONE]
    pos % dir = pos % dir / norm2(pos % dir)

    call uni2 % distance(d, surfIdx, pos)

    ref = sqrt(40.0_defReal**2 + 0.4_defReal**2)
#line 290 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 290) )
  if (anyExceptions()) return
#line 291 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 291 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(-4, surfIdx , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! From outside -> Hit
    pos % r = [-1.5_defReal, 0.6_defReal, 0.5_defReal]
    pos % dir = [ONE, ZERO, ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 3

    call uni2 % distance(d, surfIdx, pos)

    ref = HALF * SQRT2
#line 304 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 304) )
  if (anyExceptions()) return
#line 305 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
#line 305 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(-7, surfIdx , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 305) )
  if (anyExceptions()) return
#line 306 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

  end subroutine test_distance

  !!
  !! Test cell-to cell crossing
  !!
!@Test
  subroutine test_cross()
    type(coord)       :: pos
    integer(shortInt) :: idx

    ! *** 3D Lattice
    ! Cross inside
    pos % r = [-1.0_defReal, 0.0_defReal, -0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 1

    call uni1 % cross(pos, -4)

#line 328 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(4, pos % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 328) )
  if (anyExceptions()) return
#line 329 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Cross from outside
    pos % r = [1.0_defReal, 2.0_defReal, -0.5_defReal]
    pos % dir = [ONE, -ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 13

    call uni1 % cross(pos, -7)

#line 338 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(6, pos % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 338) )
  if (anyExceptions()) return
#line 339 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Cross to outside
    pos % r   = [1.5_defReal, 1.0_defReal, -1.0_defReal]
    pos % dir = [ONE, ZERO, ZERO]
    pos % localID = 6

    call uni1 % cross(pos, -2)

#line 347 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(13, pos % localID , &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 347) )
  if (anyExceptions()) return
#line 348 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"


    ! *** 2D Lattice
    pos % r = [0.0_defReal, 0.0_defReal, 16.5_defReal]
    pos % dir = [ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 1

    call uni2 % cross(pos, -2)

#line 360 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(2, pos % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 360) )
  if (anyExceptions()) return
#line 361 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Cross from outside
    pos % r = [-1.0_defReal, -0.5_defReal, -78.5_defReal]
    pos % dir = [ONE, ONE, ZERO]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 3

    call uni2 % cross(pos, -7)

#line 370 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(1, pos % localID, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 370) )
  if (anyExceptions()) return
#line 371 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

  end subroutine test_cross

  !!
  !! Test cell offset
  !!
!@Test
  subroutine test_cellOffset()
    type(coord)                 :: pos
    real(defReal), dimension(3) :: ref
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! ** 3D lattice
    ! Inside
    pos % r = [0.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 8
    pos % cellIdx = 0
    pos % localId = 11

    ref = [0.0_defReal, 1.0_defReal, 1.5_defReal]
#line 393 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, uni1 % cellOffset(pos), TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 393) )
  if (anyExceptions()) return
#line 394 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Outside
    pos % r = [-7.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 13

    ref = ZERO
#line 402 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, uni1 % cellOffset(pos), TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 402) )
  if (anyExceptions()) return
#line 403 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! ** 2D Lattice
    pos % r = [0.5_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % uniIdx  = 3
    pos % cellIdx = 0
    pos % localId = 2

    ref = [0.5_defReal, 0.0_defReal, 0.0_defReal]
#line 413 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, uni2 % cellOffset(pos), TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 413) )
  if (anyExceptions()) return
#line 414 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

    ! Outside
    pos % r = [-7.0_defReal, 0.0_defReal, 0.5_defReal]
    pos % dir = [-ONE, ONE, -ONE]
    pos % dir = pos % dir / norm2(pos % dir)
    pos % localId = 3

    ref = ZERO
#line 422 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"
  call assertEqual(ref, uni2 % cellOffset(pos), TOL, &
 & location=SourceLocation( &
 & 'latUniverse_test.f90', &
 & 422) )
  if (anyExceptions()) return
#line 423 "/home/mhgk4/SCONE/Geometry/Universes/Tests/latUniverse_test.f90"

  end subroutine test_cellOffset


end module latUniverse_test

module WraplatUniverse_test
   use pFUnit_mod
   use latUniverse_test
   implicit none
   private

contains


end module WraplatUniverse_test

function latUniverse_test_suite() result(suite)
   use pFUnit_mod
   use latUniverse_test
   use WraplatUniverse_test
   type (TestSuite) :: suite

   suite = newTestSuite('latUniverse_test_suite')

   call suite%addTest(newTestMethod('test_misc', test_misc, setUp, clean))

   call suite%addTest(newTestMethod('test_enter', test_enter, setUp, clean))

   call suite%addTest(newTestMethod('test_distance', test_distance, setUp, clean))

   call suite%addTest(newTestMethod('test_cross', test_cross, setUp, clean))

   call suite%addTest(newTestMethod('test_cellOffset', test_cellOffset, setUp, clean))


end function latUniverse_test_suite

