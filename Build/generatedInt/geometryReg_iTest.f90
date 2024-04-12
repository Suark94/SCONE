module geometryReg_iTest

  use numPrecision
  use genericProcedures, only : numToChar
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : charToDict
  use geometry_inter,    only : geometry
  use field_inter,       only : field
  use materialMenu_mod,  only : mm_init => init, mm_kill => kill
  use geometryReg_mod,   only : gr_geomIdx => geomIdx, gr_geomPtr => geomPtr,&
                                gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr, gr_kill => kill
  use geometryFactory_func, only : new_geometry
  use fieldFactory_func,    only : new_field
  use pfUnit_mod
  implicit none

  ! Material definition
  character(*), parameter :: MAT_DEF = &
  " fuel  {temp 8; composition {}} &
  & water {temp 2; composition{}}"

  ! Geometry definition
  character(*), parameter :: GEOM_DEF = &
  " type geometryStd; boundary (0); graph {type shrunk;} &
  & surfaces { surf1 {id 1; type sphere; origin (0.0 0.0 0.0); radius 2.0; }} &
  & cells {} &
  & universes { root {id 1; type rootUniverse; border 1; fill fuel;}}"

  ! Field definition
  character(*), parameter :: FIELD_DEF = &
  "type uniformScalarField; value -1.333;"




contains

  !!
  !! Prepare test enviroment
  !!
!@Before
  subroutine set_up()
    type(dictionary) :: dict

    ! Load material definitions
    call charToDict(dict, MAT_DEF)
    call mm_init(dict)

  end subroutine set_up

  !!
  !! Clean test enviroment
  !!
!@After
  subroutine clean_up()

    call mm_kill()

  end subroutine clean_up

  !!
  !! Test geometry functionality
  !!
!@Test
  subroutine test_geometry()
    type(dictionary)         :: dict
    character(nameLen)       :: name
    integer(shortInt)        :: i
    class(geometry), pointer :: ptr

    ! Get geometry definition
    call charToDict(dict, GEOM_DEF)

    ! Build 10 instances of geometry
    do i = 1, 10
      name = 'geom'//numToChar(i)
      call new_geometry(dict, name, silent=.true.)
    end do

    ! Get some indexes
    name = 'geom7'
#line 82 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertEqual(7, gr_geomIdx(name), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    name = 'geom2'
#line 85 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertEqual(2, gr_geomIdx(name), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    ! Get some pointers -> Check if associated
    ptr => gr_geomPtr(3)
#line 89 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertTrue(associated(ptr), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    ptr => gr_geomPtr(10)
#line 92 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertTrue(associated(ptr), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    ! Kill registry
    call gr_kill()

  end subroutine test_geometry

  !!
  !! Test field functionality
  !!
!@Test
  subroutine test_field()
    type(dictionary)      :: dict
    character(nameLen)    :: name
    integer(shortInt)     :: i
    class(field), pointer :: ptr

    ! Get field definition
    call charToDict(dict, FIELD_DEF)

    ! Build 10 instances of fields
    do i = 1, 10
      name = 'field'//numToChar(i)
      call new_field(dict, name)
    end do

    ! Get some indexes
    name = 'field7'
#line 120 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertEqual(7, gr_fieldIdx(name), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    name = 'field2'
#line 123 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertEqual(2, gr_fieldIdx(name), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    ! Get some pointers -> Check if associated
    ptr => gr_fieldPtr(3)
#line 127 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertTrue(associated(ptr), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    ptr => gr_fieldPtr(10)
#line 130 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"
  call assertTrue(associated(ptr), &
 & location=SourceLocation( &
 & 'geometryReg_iTest.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/Geometry/Tests/geometryReg_iTest.f90"

    ! Kill registry
    call gr_kill()

  end subroutine test_field


end module geometryReg_iTest

module WrapgeometryReg_iTest
   use pFUnit_mod
   use geometryReg_iTest
   implicit none
   private

contains


end module WrapgeometryReg_iTest

function geometryReg_iTest_suite() result(suite)
   use pFUnit_mod
   use geometryReg_iTest
   use WrapgeometryReg_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('geometryReg_iTest_suite')

   call suite%addTest(newTestMethod('test_geometry', test_geometry, set_up, clean_up))

   call suite%addTest(newTestMethod('test_field', test_field, set_up, clean_up))


end function geometryReg_iTest_suite

