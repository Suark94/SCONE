module coord_test

  use numPrecision
  use coord_class, only : coord, coordList
  use pfUnit_mod

  implicit none

  ! Variables
  type(coordList) :: coords


contains

  !!
  !! Set up test envoroment
  !!
  !! Coords Placed at 3 levels
  !!   Level 1 -> global
  !!   Level 2 -> Rotated universe
  !!   Level 3 -> Translated universe
  !!
!@Before
  subroutine set_up()
    real(defReal), dimension(3,3) :: mat

    ! Set Nesting
    coords % nesting = 3
    coords % matIdx = 2
    coords % uniqueId = 7

    ! Set Level 1
    coords % lvl(1) % r   = [1.0_defReal, 0.0_defReal, -1.0_defReal]
    coords % lvl(1) % dir = [ZERO, ONE, ZERO]
    coords % lvl(1) % uniIdx    = 1
    coords % lvl(1) % uniRootID = 1
    coords % lvl(1) % localID   = 1
    coords % lvl(1) % cellIdx   = 1

    ! Set Level 2
    ! Rotation Y -> Z; Z -> -Y
    mat = ZERO
    mat(1,1) = ONE
    mat(3,2) = -ONE
    mat(2,3) = ONE
    coords % lvl(2) % r   = [1.0_defReal, 0.0_defReal, -1.0_defReal]
    coords % lvl(2) % dir = [ZERO, ZERO, -ONE]
    coords % lvl(2) % uniIdx    = 2
    coords % lvl(2) % uniRootID = 6
    coords % lvl(2) % localID   = 3
    coords % lvl(2) % isRotated = .true.
    coords % lvl(2) % rotMat    = mat
    coords % lvl(2) % cellIdx   = 3

    ! Set Level 3
    ! Translation to origin
    coords % lvl(3) % r = ZERO
    coords % lvl(3) % dir = [ZERO, ZERO, -ONE]
    coords % lvl(3) % uniIdx    = 4
    coords % lvl(3) % uniRootID = 12
    coords % lvl(3) % localID   = 2
    coords % lvl(3) % cellIdx   = 0

  end subroutine set_up

  !!
  !! Clean test enviroment
  !!
!@After
  subroutine clean_up()

    call coords % kill()

  end subroutine clean_up

  !!
  !! Test state changing procedures
  !!
!@Test
  subroutine test_changing_state()

    ! Test State
#line 83 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertTrue(coords % isPlaced(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 84 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % isAbove(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 85 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % isUninitialised(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

    ! Change to above
    call coords % takeAboveGeom()
#line 89 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % isPlaced(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 90 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertTrue(coords % isAbove(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 91 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % isUninitialised(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 91) )
  if (anyExceptions()) return
#line 92 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

    ! Chenge to uninitialised
    call coords % kill()
#line 95 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % isPlaced(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 96 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % isAbove(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 97 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertTrue(coords % isUninitialised(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 97) )
  if (anyExceptions()) return
#line 98 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

  end subroutine test_changing_state

  !!
  !! Test nesting level changes & cell inquiry
  !!
!@Test
  subroutine test_nesting_level()

    ! Move deeper
    call coords % addLevel()
#line 109 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(4, coords % nesting, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 110 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(0, coords % cell(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

    ! Move to higher level
    call coords % decreaseLevel(2)
#line 114 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(2, coords % nesting, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 115 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(3, coords % cell(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

  end subroutine test_nesting_level

  !!
  !! Test rotation
  !!
  !! Verifies only the deflection by mu !
  !!
!@Test
  subroutine test_rotation()
    real(defReal)               :: mu, phi
    real(defReal), dimension(3) :: u1, u2, u3
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    mu = 0.3_defReal
    phi = 2.1_defReal

    ! Save pre-rotation direction
    u1 = coords % lvl(1) % dir
    u2 = coords % lvl(2) % dir
    u3 = coords % lvl(3) % dir

    call coords % rotate(mu, phi)

    ! Verify deflection
#line 141 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(mu, dot_product(u1, coords % lvl(1) % dir), TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 142 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(mu, dot_product(u2, coords % lvl(2) % dir), TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 143 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(mu, dot_product(u3, coords % lvl(3) % dir), TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 143) )
  if (anyExceptions()) return
#line 144 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

  end subroutine test_rotation

  !!
  !! Test direction assigment
  !!
!@Test
  subroutine test_direction_assigment()
    real(defReal), dimension(3) :: u1, u2, u3
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Save pre-rotation direction
    u1 = coords % lvl(1) % dir
    u2 = coords % lvl(2) % dir
    u3 = coords % lvl(3) % dir

    ! Invert direction
    call coords % assignDirection(-coords % lvl(1) % dir)

    ! Verify
#line 164 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(-u1, coords % lvl(1) % dir, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 165 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(-u2, coords % lvl(2) % dir, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 166 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(-u3, coords % lvl(3) % dir, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

  end subroutine test_direction_assigment

  !!
  !! Test Movment
  !!
!@Test
  subroutine test_movement()
    real(defReal), dimension(3) :: u1, u2, u3, r1, r2, r3
    real(defReal)               :: d
    real(defReal), parameter    :: TOL = 1.0E-7_defReal

    ! Move local
    d = 0.3_defReal
    r1 = coords % lvl(1) % r
    r2 = coords % lvl(2) % r
    r3 = coords % lvl(3) % r
    u1 = coords % lvl(1) % dir
    u2 = coords % lvl(2) % dir
    u3 = coords % lvl(3) % dir

    call coords % moveLocal(d, 3)

    ! Verify
#line 191 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(r1 + d*u1, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 192 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(r2 + d*u2, coords % lvl(2) % r, TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 193 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(r3 + d*u3, coords % lvl(3) % r, TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 194 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertTrue(coords % isPlaced(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 194) )
  if (anyExceptions()) return
#line 195 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

    ! Move Global
    d = -13.0_defReal
    r1 = coords % lvl(1) % r
    u1 = coords % lvl(1) % dir

    call coords % moveGlobal(d)

    ! Verify
#line 204 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertEqual(r1 + d*u1, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 204) )
  if (anyExceptions()) return
#line 205 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
#line 205 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertTrue(coords % isAbove(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 205) )
  if (anyExceptions()) return
#line 206 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

  end subroutine test_movement

  !!
  !! Test coord validation
  !!
!@Test
  subroutine test_coord_valid()

#line 215 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertTrue(coords % lvl(1:coords % nesting) % isValid(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 215) )
  if (anyExceptions()) return
#line 216 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"

    call coords % kill()

#line 219 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"
  call assertFalse(coords % lvl % isValid(), &
 & location=SourceLocation( &
 & 'coord_test.f90', &
 & 219) )
  if (anyExceptions()) return
#line 220 "/home/mhgk4/SCONE/Geometry/Tests/coord_test.f90"


  end subroutine test_coord_valid

end module coord_test

module Wrapcoord_test
   use pFUnit_mod
   use coord_test
   implicit none
   private

contains


end module Wrapcoord_test

function coord_test_suite() result(suite)
   use pFUnit_mod
   use coord_test
   use Wrapcoord_test
   type (TestSuite) :: suite

   suite = newTestSuite('coord_test_suite')

   call suite%addTest(newTestMethod('test_changing_state', test_changing_state, set_up, clean_up))

   call suite%addTest(newTestMethod('test_nesting_level', test_nesting_level, set_up, clean_up))

   call suite%addTest(newTestMethod('test_rotation', test_rotation, set_up, clean_up))

   call suite%addTest(newTestMethod('test_direction_assigment', test_direction_assigment, set_up, clean_up))

   call suite%addTest(newTestMethod('test_movement', test_movement, set_up, clean_up))

   call suite%addTest(newTestMethod('test_coord_valid', test_coord_valid, set_up, clean_up))


end function coord_test_suite

