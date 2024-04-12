module plane_test

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC, REFLECTIVE_BC, INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use plane_class,       only : plane
  use pfUnit_mod

  implicit none

  character(*), parameter :: PLANE_DEF = "id 7; coeffs (1.0 1.0 1.0 3.0);"
  type(plane)             :: surf

contains

  !!
  !! Build the surface
  !!
!@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, PLANE_DEF)
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
#line 48 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(7, surf % id(), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Change ID
    call surf % setID(1)
#line 52 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(1, surf % id(), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 52) )
  if (anyExceptions()) return
#line 53 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Bounding box
    ref = [-INF, -INF, -INF, INF, INF, INF]
    aabb = surf % boundingBox()
#line 57 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Name
#line 60 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual('plane', surf % myType(), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

  end subroutine testMisc

  !!
  !! Test boundary conditions
  !!
!@Test
  subroutine testBC()
    real(defReal), dimension(3) :: r, u, r_pre, u_pre

    ! Set Boundary conditions
    ! Should ignore extra entries
    call surf % setBC([VACUUM_BC, REFLECTIVE_BC, REFLECTIVE_BC])

    ! Apply BCs
    r = [1.0_defReal, 1.0_defReal, 1.0_defReal]
    u = [ZERO, ONE, ZERO]
    r_pre = r
    u_pre = u

    ! Explicit
    call surf % explicitBC(r, u)
#line 83 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
#line 84 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Transform
    call surf % transformBC(r, u)
#line 88 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
#line 89 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

  end subroutine testBC

  !!
  !! Test halfspace membership
  !!
!@Test
  subroutine testHalfspace()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: eps


    r = [ONE, ONE, ONE]
    u = [ONE, ONE, ONE]
    u = u /norm2(u)

    ! At the surface
#line 107 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertTrue(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! In within SURF_TOL
    eps = -0.5_defReal * SURF_TOL * sqrt(3.0_defReal)
#line 111 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertTrue(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! In outside SURF_TOL
    eps = -1.00001_defReal * SURF_TOL * sqrt(3.0_defReal)
#line 115 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertFalse(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Well Inside
    eps = -ONE
#line 119 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertFalse(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Well Outside
    eps = ONE
#line 123 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertTrue(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! Tangent particle should use position
    u2 = [-ONE, ZERO, ONE]
    u2 = u2/norm2(u2)
    eps = HALF * SURF_TOL

#line 130 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertTrue( surf % halfspace(r + eps*u, u2), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
#line 131 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertFalse( surf % halfspace(r - eps*u, u2), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
!@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u, u2
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    r = [-ONE, ZERO, ZERO]
    u = [ONE, ZERO, ZERO]

    ! ** Inside the surface
    ref = 4.0_defReal
#line 149 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! ** Exactly at the surface
    r(1) = 3.0_defReal
#line 153 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(INF, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! ** Inside within surface Tolerance
    r(1) = 3.0_defReal - SURF_TOL
#line 157 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(INF, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! ** Inside close to the surface tolerance
    r(1) = 3.0_defReal - sqrt(3.0_defReal) * 1.1_defReal * SURF_TOL
    ref = sqrt(3.0_defReal) * 1.1_defReal * SURF_TOL
    ! Use very liberal tolerance becouse the distance is veary small and sensitive to FP precision
#line 163 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(ref, surf % distance(r, u), 1.0E-2 * ref, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! ** Outside the surface
    r(1) = 4.0_defReal

    ! +ve Direction
#line 169 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(INF, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 169) )
  if (anyExceptions()) return
#line 170 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! -ve Direction
    ref = 1.0_defReal
#line 173 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(ref, surf % distance(r, -u), TOL * ref, &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 173) )
  if (anyExceptions()) return
#line 174 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

    ! ** Parallel to plane
    r = [-ONE, ZERO, ZERO]
    u2 = [-ONE, ZERO, ONE]
    u2 = u2/norm2(u2)
#line 179 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"
  call assertEqual(INF, surf % distance(r, u2), &
 & location=SourceLocation( &
 & 'plane_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/plane_test.f90"

  end subroutine testDistance

end module plane_test

module Wrapplane_test
   use pFUnit_mod
   use plane_test
   implicit none
   private

contains


end module Wrapplane_test

function plane_test_suite() result(suite)
   use pFUnit_mod
   use plane_test
   use Wrapplane_test
   type (TestSuite) :: suite

   suite = newTestSuite('plane_test_suite')

   call suite%addTest(newTestMethod('testMisc', testMisc, setUp, cleanUp))

   call suite%addTest(newTestMethod('testBC', testBC, setUp, cleanUp))

   call suite%addTest(newTestMethod('testHalfspace', testHalfspace, setUp, cleanUp))

   call suite%addTest(newTestMethod('testDistance', testDistance, setUp, cleanUp))


end function plane_test_suite

