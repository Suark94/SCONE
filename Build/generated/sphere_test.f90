module sphere_test

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC, REFLECTIVE_BC, INF
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use sphere_class,       only : sphere
  use pfUnit_mod

  implicit none

  character(*), parameter :: SPH_DEF = "id 7; origin (1.0 2.0 1.0); radius 2.0; "
  type(sphere)            :: surf

contains

  !!
  !! Build the surface
  !!
!@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, SPH_DEF)
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
#line 48 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(7, surf % id(), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Change ID
    call surf % setID(1)
#line 52 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(1, surf % id(), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 52) )
  if (anyExceptions()) return
#line 53 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Bounding box
    ref = [-1.0_defReal, 0.0_defReal, -1.0_defReal, 3.0_defReal, 4.0_defReal, 3.0_defReal]
    aabb = surf % boundingBox()
#line 57 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, aabb, TOL, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Name
#line 60 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual('sphere', surf % myType(), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

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
    r = [1.0_defReal, 0.0_defReal, 1.0_defReal]
    u = [ZERO, ONE, ZERO]
    r_pre = r
    u_pre = u

    ! Explicit
    call surf % explicitBC(r, u)
#line 83 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
#line 84 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Transform
    call surf % transformBC(r, u)
#line 88 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(r_pre, r, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
#line 89 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(u_pre, u, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

  end subroutine testBC

  !!
  !! Test halfspace membership
  !!
!@Test
  subroutine testHalfspace()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: eps

    ! Choose point At x-axis going at the surface
    r = [-1.0_defReal, 2.0_defReal, 1.0_defReal]
    u = [ONE, ZERO, ZERO]

    ! At the surface
#line 106 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertFalse(surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Out within SURF_TOL
    eps = -0.5_defReal * SURF_TOL
#line 110 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertFalse(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Out outside SURF_TOL
    eps = -1.00001_defReal * SURF_TOL
#line 114 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertTrue(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Well Outside
    eps = -ONE
#line 118 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertTrue(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Well within
    eps = ONE
#line 122 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertFalse(surf % halfspace(r + eps*u, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Tangent particle whould be outside
    u = [ZERO, ONE, ZERO]
#line 126 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertTrue( surf % halfspace(r, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

  end subroutine testHalfspace

  !!
  !! Test distance calculations
  !!
!@Test
  subroutine testDistance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref
    real(defReal), parameter :: TOL = 1.0E-7

    ! **Outside the sphere
    r = [-2.0_defReal, 2.0_defReal, 1.0_defReal]

    ! Perpendicular impact
    u = [ONE, ZERO, ZERO]
    ref = 1.0_defReal
#line 145 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Oblique impact
    ref = 1.5_defReal
    ! Calculate ange from Cosine theorem
    u(1) = (ref**2 + 3.0_defReal**2 - 2.0_defReal**2) / (TWO * ref * 3.0_defReal)
    u(2) = sqrt(ONE - u(1)**2)
    u = u / norm2(u)
#line 153 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, surf % distance(r, u), TOL * ref, &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! **Exactly at the surface
    r  = [-1.0_defReal, 2.0_defReal, 1.0_defReal]

    ! Tangent particle
    u = [ZERO, ONE, ZERO]
#line 160 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(INF, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! Particle going inside
    u = [ONE, ZERO, ZERO]
    ref = 4.0_defReal
#line 165 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! ** Outside in Surface tolerance
    r = r - [ONE, ZERO, ZERO] * HALF * SURF_TOL
    ref = 4.0_defReal + HALF * SURF_TOL
#line 170 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! **Inside the surface
    ! +ve direction
    r = [-0.5_defReal, 2.0_defReal, 1.0_defReal]
    u = [ONE, ZERO, ZERO]
    ref = 3.5_defReal
#line 177 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, surf % distance(r, u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

    ! -ve direction
    ref = 0.5_defReal
#line 181 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"
  call assertEqual(ref, surf % distance(r, -u), &
 & location=SourceLocation( &
 & 'sphere_test.f90', &
 & 181) )
  if (anyExceptions()) return
#line 182 "/home/mhgk4/SCONE/Geometry/Surfaces/Tests/sphere_test.f90"

  end subroutine testDistance

end module sphere_test

module Wrapsphere_test
   use pFUnit_mod
   use sphere_test
   implicit none
   private

contains


end module Wrapsphere_test

function sphere_test_suite() result(suite)
   use pFUnit_mod
   use sphere_test
   use Wrapsphere_test
   type (TestSuite) :: suite

   suite = newTestSuite('sphere_test_suite')

   call suite%addTest(newTestMethod('testMisc', testMisc, setUp, cleanUp))

   call suite%addTest(newTestMethod('testBC', testBC, setUp, cleanUp))

   call suite%addTest(newTestMethod('testHalfspace', testHalfspace, setUp, cleanUp))

   call suite%addTest(newTestMethod('testDistance', testDistance, setUp, cleanUp))


end function sphere_test_suite

