module uniformScalarField_test

  use numPrecision
  use dictionary_class,         only : dictionary
  use particle_class,           only : particle
  use field_inter,              only : field
  use scalarField_inter,        only : scalarField, scalarField_CptrCast
  use uniformScalarField_class, only : uniformScalarField, uniformScalarField_TptrCast
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test Uniform Scalar Field
  !!
!@Test
  subroutine test_uniformScalarField()
    type(uniformScalarField), target  :: fieldT
    class(field), pointer             :: ref
    class(scalarField), pointer       :: ptr
    type(uniformScalarField), pointer :: ptr2
    type(dictionary)                  :: dict
    type(particle)                    :: p
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()

    ptr => scalarField_CptrCast(ref)
    ptr2 => uniformScalarField_TptrCast(ref)
#line 33 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
  call assertFalse(associated(ptr), &
 & location=SourceLocation( &
 & 'uniformScalarField_test.f90', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
#line 34 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
  call assertFalse(associated(ptr2), &
 & location=SourceLocation( &
 & 'uniformScalarField_test.f90', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"

    ! Test valid pointers
    ref => fieldT

    ptr => scalarField_CptrCast(ref)
    ptr2 => uniformScalarField_TptrCast(ref)

#line 42 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
  call assertTrue(associated(ptr, fieldT), &
 & location=SourceLocation( &
 & 'uniformScalarField_test.f90', &
 & 42) )
  if (anyExceptions()) return
#line 43 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
#line 43 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
  call assertTrue(associated(ptr2, fieldT), &
 & location=SourceLocation( &
 & 'uniformScalarField_test.f90', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"

    ! Initialise field
    call dict % init(2)
    call dict % store('type', 'uniformVectorField')
    call dict % store('value', 9.6_defReal)

    call fieldT % init(dict)

    ! Check value
#line 53 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"
  call assertEqual(9.6_defReal, fieldT % at(p), TOL, &
 & location=SourceLocation( &
 & 'uniformScalarField_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/Geometry/Fields/ScalarFields/Tests/uniformScalarField_test.f90"

    ! Kill
    call fieldT % kill()

  end subroutine test_uniformScalarField

end module uniformScalarField_test

module WrapuniformScalarField_test
   use pFUnit_mod
   use uniformScalarField_test
   implicit none
   private

contains


end module WrapuniformScalarField_test

function uniformScalarField_test_suite() result(suite)
   use pFUnit_mod
   use uniformScalarField_test
   use WrapuniformScalarField_test
   type (TestSuite) :: suite

   suite = newTestSuite('uniformScalarField_test_suite')

   call suite%addTest(newTestMethod('test_uniformScalarField', test_uniformScalarField))


end function uniformScalarField_test_suite

