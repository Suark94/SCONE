module box_test

  use numPrecision
  use universalVariables
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use box_class,          only : box
  use pfUnit_mod

  implicit none

  character(*), parameter :: BOX_DEF = "id 7; origin (1.0 2.0 1.0); halfwidth (1.0 2.0 3.0); "
  type(box)               :: surf

contains

  !!
  !! Build the surface
  !!
!@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, BOX_DEF)
    call surf % init(dict)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
!@After
  subroutine cleanUp()

    call surf % kill()

  end subroutine cleanUp

  !!
  !! Test Misc functionality
  !!
!@Test
  subroutine testMisc()
    real(defReal), dimension(6) :: aabb, ref
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! Test ID
#line 48 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(7, surf % id(), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Change ID
    call surf % setID(1)
#line 52 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(1, surf % id(), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 52) )
  if (anyExceptions()) return
#line 53 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Bounding box
    ref = [0.0_defReal, 0.0_defReal, -2.0_defReal, 2.0_defReal, 4.0_defReal, 4.0_defReal]
    aabb = surf % boundingBox()
#line 57 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Name
#line 60 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual('box', surf % myType(), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
!@Test
  subroutine testBC()
    real(defReal), dimension(3) :: r, u, r_ref, u_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Assign boundary conditions
    call surf % setBC([VACUUM_BC, REFLECTIVE_BC, PERIODIC_BC, &
                      PERIODIC_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! ** Explicit BC
    ! Vacuum surface
    r = [0.0_defReal, 1.0_defReal, -1.0_defReal]
    u = [ONE, ZERO, ZERO]
    r_ref = r
    u_ref = u
    call surf % explicitBC(r, u)
#line 83 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 84 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Reflective surface
    r = [1.0_defReal, 1.0_defReal, -2.0_defReal]
    u = [ZERO, ZERO, -ONE]
    r_ref = r
    u_ref = -u
    call surf % explicitBC(r, u)
#line 92 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 93 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 93) )
  if (anyExceptions()) return
#line 94 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Periodic surface
    r = [1.0_defReal, 0.0_defReal, -1.0_defReal]
    u = [ZERO, -ONE, ZERO]
    r_ref = [1.0_defReal, 4.0_defReal, -1.0_defReal]
    u_ref = u
    call surf % explicitBC(r, u)
#line 101 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 102 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Try the corner
    r = [2.0_defReal, 0.0_defReal, 4.0_defReal]
    u = [ONE, -ONE, ONE]
    u = u / norm2(u)

    r_ref = [2.0_defReal, 4.0_defReal, 4.0_defReal]
    u_ref = [-ONE, -ONE, -ONE]
    u_ref = u_ref / norm2(u_ref)

    call surf % explicitBC(r, u)
#line 114 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 115 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! ** Transform BC
    !
    r = [4.5_defReal, -10.0_defReal, 4.0_defReal]
    u = [ONE, ONE, ONE]
    u = u / norm2(u)

    r_ref = [-0.5_defReal, 2.0_defReal, 4.0_defReal]
    u_ref = [-ONE, ONE, -ONE]
    u_ref = u_ref / norm2(u_ref)

    call surf % transformBC(r, u)
#line 128 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 129 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 129) )
  if (anyExceptions()) return
#line 130 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Try the corner
    r = [2.0_defReal, 0.0_defReal, 4.0_defReal]
    u = [ONE, -ONE, ONE]
    u = u / norm2(u)

    r_ref = [2.0_defReal, 4.0_defReal, 4.0_defReal]
    u_ref = [-ONE, -ONE, -ONE]
    u_ref = u_ref / norm2(u_ref)

    call surf % transformBC(r, u)
#line 141 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(r_ref, r, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 142 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(u_ref, u, TOL, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

  end subroutine testBC

  !!
  !! Test halfspace membership
  !!
!@Test
  subroutine testHalfspace()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: eps

    ! ** Well inside the box
    r = [0.5_defReal, 1.0_defReal, 3.6_defReal]
    u = [ZERO, ZERO, ONE]
#line 157 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertFalse(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Difrent octant
    r = [1.5_defReal, 1.0_defReal, 0.5_defReal]
#line 161 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertFalse(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 161) )
  if (anyExceptions()) return
#line 162 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! ** Well outside the box
    ! Make sure point is at one of the planes outside
    r = [-0.5_defReal, 0.0_defReal, 3.6_defReal]
#line 166 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Diffrent octant
    r = [2.0_defReal, 5.0_defReal, 0.5_defReal]
#line 170 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! ** Proximity of a surface
    r = [1.5_defReal, 4.0_defReal, 0.5_defReal]
    u = [ZERO, ONE, ZERO]
#line 175 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Inside, within surface tolerance
    eps = -HALF * SURF_TOL
#line 179 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Inside, a bit outside surface tolerance
    eps = -1.00001_defReal * SURF_TOL
#line 183 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertFalse(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Parallel to the surface
    u2 = [ONE, ZERO, ZERO]
    eps = HALF * SURF_TOL

#line 189 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r + eps*u, u2), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 190 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertFalse(surf % halfspace(r - eps*u, u2), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 190) )
  if (anyExceptions()) return
#line 191 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
!@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! ** Outside the box
    r = [-2.0_defReal, 0.001_defReal, 1.0_defReal]

    ! Direct impact
    u = [ONE, ZERO, ZERO]
    ref = 2.0_defReal
#line 209 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Moving away
#line 212 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(INF, surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 212) )
  if (anyExceptions()) return
#line 213 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Miss
    u = [ONE, -0.002_defReal, ONE]
    u = u /norm2(u)
#line 217 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(INF, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Oblique Hit
    u = [ONE, ZERO, ONE]
    u = u /norm2(u)
    ref = 2.0_defReal * SQRT2
#line 223 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 223) )
  if (anyExceptions()) return
#line 224 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Corner skim [0.0, 0.0, 4.0]
    ! Use dirty values
    r = [-2.7_defReal, -0.3_defReal, 1.0_defReal/3.0_defReal]
    u = [2.7_defReal, 0.3_defReal, 3.0_defReal + 2.0_defReal/3.0_defReal]
    u = u /norm2(u)
#line 230 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(INF, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Parallel
    r = [-2.0_defReal, 0.000_defReal, 1.0_defReal]
    u = [ONE, ZERO, ZERO]
#line 235 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(INF, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! ** At the surface
    r = [ONE, TWO, -2.0_defReal]
    u = [ZERO, ZERO, ONE]
    ref = 6.0_defReal
#line 241 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 242 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(INF, surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! Outside within surface tolerance
    r = [ONE, TWO, -2.0_defReal - HALF * SURF_TOL]
    ref = 6.0_defReal + HALF * SURF_TOL
#line 247 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 247) )
  if (anyExceptions()) return
#line 248 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ! ** Inside
    r = [ONE, TWO, -1.0_defReal]
    u = [ZERO, ONE, ONE]
    u = u/norm2(u)

    ref = TWO * SQRT2
#line 255 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 255) )
  if (anyExceptions()) return
#line 256 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

    ref = SQRT2
#line 258 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertEqual(ref, surf % distance(r, -u), TOL*ref, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 258) )
  if (anyExceptions()) return
#line 259 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

  end subroutine testDistance

  !!
  !! Test Edge Cases
  !!
  !! Test unlikley cases to make sure that halfspace + distance
  !! procedures allow particle to escape
  !!
!@Test
  subroutine testEdgeCases()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: eps, d
    logical(defBool)            :: hs

    ! ** Corner
    ! * Particle is almost at the corner
    !   Either it is outside or can escape with a short movment in next step
    !
    ! Currently does not work for Y-position in plane (r(2) == ZERO)
    ! See `distance` documentation in box_class
    ! This comment applies to all cases in this test.
    !
    eps =  5.0_defReal * epsilon(eps)
    r = [2.0_defReal-eps, eps, 4.0_defReal-eps]
    u = [HALF, ZERO, -ONE]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
#line 289 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 289) )
  if (anyExceptions()) return
#line 290 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 290 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 290) )
  if (anyExceptions()) return
#line 291 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
    end if

    ! Try in other direction
    u = [-ONE, ZERO, HALF]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
#line 299 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 300 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 300) )
  if (anyExceptions()) return
#line 301 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
    end if

    ! Try asymertic corner
    r = [2.0_defReal-TWO*eps, eps, 4.0_defReal-eps]
    u = [HALF, ZERO, -ONE]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
#line 310 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 310) )
  if (anyExceptions()) return
#line 311 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 311 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 311) )
  if (anyExceptions()) return
#line 312 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
    end if

    ! Asymetric corner position
    ! Try other direction
    r = [2.0_defReal-eps, eps, 4.0_defReal-TWO*eps]
    u = [-ONE, ZERO, HALF]
    u = u/norm2(u)
    hs = surf % halfspace(r, u)
    if (.not.hs) then ! Perform small movment
      d = surf % distance(r, u)
#line 322 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue( abs(d) < 1.0E-6, &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 322) )
  if (anyExceptions()) return
#line 323 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
#line 323 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertTrue(surf % halfspace(r + d*u, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 323) )
  if (anyExceptions()) return
#line 324 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
    end if

  end subroutine testEdgeCases

  !!
  !! Test problematic cases
  !!
!@Test
  subroutine test_problems()
    type(box)            :: b_surf
    type(dictionary)     :: dict
    real(defReal), dimension(3) :: r, u

    call dict % init(5)
    call dict % store('type','box')
    call dict % store('id', 7)
    call dict % store('origin', [ZERO, ZERO, ZERO])
    call dict % store('halfwidth', [8.00_defReal, 1.26_defReal, 1.26_defReal])
    call b_surf % init(dict)

    r = [-7.59_defReal, 1.26_defReal, ZERO]
    u = [ZERO, -ONE, ZERO]
#line 346 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"
  call assertFalse(b_surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'box_test.f90', &
 & 346) )
  if (anyExceptions()) return
#line 347 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/box_test.f90"

  end subroutine test_problems

end module box_test

module Wrapbox_test
   use pFUnit_mod
   use box_test
   implicit none
   private

contains


end module Wrapbox_test

function box_test_suite() result(suite)
   use pFUnit_mod
   use box_test
   use Wrapbox_test
   type (TestSuite) :: suite

   suite = newTestSuite('box_test_suite')

   call suite%addTest(newTestMethod('testMisc', testMisc, setUp, cleanUp))

   call suite%addTest(newTestMethod('testBC', testBC, setUp, cleanUp))

   call suite%addTest(newTestMethod('testHalfspace', testHalfspace, setUp, cleanUp))

   call suite%addTest(newTestMethod('testDistance', testDistance, setUp, cleanUp))

   call suite%addTest(newTestMethod('testEdgeCases', testEdgeCases, setUp, cleanUp))

   call suite%addTest(newTestMethod('test_problems', test_problems, setUp, cleanUp))


end function box_test_suite

