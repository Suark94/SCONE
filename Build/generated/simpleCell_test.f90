module simpleCell_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use surfaceShelf_class, only : surfaceShelf
  use simpleCell_class,   only : simpleCell
  use pFUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SURF_DEF = "&
  & surf1 { id 13; type sphere; origin (0.0 0.0 0.0); radius 2.0;} &
  & surf2 { id 4;  type xPlane; x0 0.0;} &
  & surf3 { id 99; type yPlane; y0 0.0;}"

  ! Note that fill is not really needed to build a cell. It is used by cellShelf only
  character(*), parameter :: CELL_DEF = "&
  & id 2; type simpleCell; surfaces (-13 4 99 ); filltype outside; "


  ! Variables
  type(surfaceShelf) :: surfs
  type(simpleCell)   :: cell

contains

  !!
  !! Build the cell
  !!
!@Before
  subroutine setUp()
    type(dictionary) :: dict

    call charToDict(dict, SURF_DEF)
    call surfs % init(dict)
    call dict % kill()
    call charToDict(dict, CELL_DEF)
    call cell % init(dict, surfs)

  end subroutine setUp

  !!
  !! Clean after tests
  !!
!@After
  subroutine cleanUp()

    call surfs % kill()
    call cell % kill()

  end subroutine cleanUp

  !!
  !! Test Miscellaneous functionality
  !!
!@Test
  subroutine test_misc()

    ! Test Id
#line 62 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(2, cell % id(), &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

    call cell % setID(7)
#line 65 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(7, cell % id(), &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

  end subroutine test_misc


  !!
  !! Test inside/outside determination
  !!
!@Test
  subroutine test_inside()
    real(defReal), dimension(3) :: r, u

    ! Few points inside
    r = [1.0_defReal, 0.3_defReal, 0.1_defReal]
    u = [ONE, ZERO, ZERO]
#line 80 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertTrue( cell % inside(r, u), &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

    r = [0.3_defReal, 1.3_defReal, 0.4_defReal]
#line 83 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertTrue( cell % inside(r, u), &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

    ! Few points outside
    r = [-0.1_defReal, 1.2_defReal, 0.1_defReal]
#line 87 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertFalse( cell % inside(r, u), &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

    r = [1.8_defReal, 0.8_defReal, 0.8_defReal]
#line 90 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertFalse( cell % inside(r, u), &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

  end subroutine test_inside

  !!
  !! Test distance calculations
  !!
!@Test
  subroutine test_distance()
    real(defReal), dimension(3) :: r, u
    real(defReal)               :: ref, d
    integer(shortInt)           :: idx, idx_ref
    real(defReal), parameter :: TOL = 1.0E-6

    ! Point inside
    ! X-Plane hit
    r = [0.3_defReal, 0.4_defReal, 0.0_defReal]
    u = [-ONE, ZERO, ZERO]
    ref = 0.3_defReal
    idx_ref = surfs % getIdx(4)
    call cell % distance(d, idx, r, u)
#line 111 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
#line 112 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(idx_ref, idx, &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

    ! Y-Plane hit
    u = [ZERO, -ONE, ZERO]
    ref = 0.4_defReal
    idx_ref = surfs % getIdx(99)
    call cell % distance(d, idx, r, u)
#line 119 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
#line 120 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(idx_ref, idx, &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

    ! Sphere hit
    u = [ONE, ZERO, ZERO]
    ref = 1.659591794_defReal
    idx_ref = surfs % getIdx(13)
    call cell % distance(d, idx, r, u)
#line 127 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(ref, d, TOL * ref, &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
#line 128 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"
  call assertEqual(idx_ref, idx, &
 & location=SourceLocation( &
 & 'simpleCell_test.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/Geometry/Cells/Tests/simpleCell_test.f90"

  end subroutine test_distance

end module simpleCell_test

module WrapsimpleCell_test
   use pFUnit_mod
   use simpleCell_test
   implicit none
   private

contains


end module WrapsimpleCell_test

function simpleCell_test_suite() result(suite)
   use pFUnit_mod
   use simpleCell_test
   use WrapsimpleCell_test
   type (TestSuite) :: suite

   suite = newTestSuite('simpleCell_test_suite')

   call suite%addTest(newTestMethod('test_misc', test_misc, setUp, cleanUp))

   call suite%addTest(newTestMethod('test_inside', test_inside, setUp, cleanUp))

   call suite%addTest(newTestMethod('test_distance', test_distance, setUp, cleanUp))


end function simpleCell_test_suite

