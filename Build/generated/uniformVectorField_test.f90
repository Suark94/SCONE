module uniformVectorField_test

  use numPrecision
  use dictionary_class,         only : dictionary
  use particle_class,           only : particle
  use field_inter,              only : field
  use vectorField_inter,        only : vectorField, vectorField_CptrCast
  use uniformVectorField_class, only : uniformVectorField, uniformVectorField_TptrCast
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test Uniform Scalar Field
  !!
!@Test
  subroutine test_uniformVectorField()
    type(uniformVectorField), target  :: fieldT
    class(field), pointer             :: ref
    class(vectorField), pointer       :: ptr
    type(uniformVectorField), pointer :: ptr2
    type(dictionary)                  :: dict
    type(particle)                    :: p
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()

    ptr => vectorField_CptrCast(ref)
    ptr2 => uniformVectorField_TptrCast(ref)
#line 33 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
  call assertFalse(associated(ptr), &
 & location=SourceLocation( &
 & 'uniformVectorField_test.f90', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
#line 34 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
  call assertFalse(associated(ptr2), &
 & location=SourceLocation( &
 & 'uniformVectorField_test.f90', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"

    ! Test valid pointers
    ref => fieldT

    ptr => vectorField_CptrCast(ref)
    ptr2 => uniformVectorField_TptrCast(ref)

#line 42 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
  call assertTrue(associated(ptr, fieldT), &
 & location=SourceLocation( &
 & 'uniformVectorField_test.f90', &
 & 42) )
  if (anyExceptions()) return
#line 43 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
#line 43 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
  call assertTrue(associated(ptr2, fieldT), &
 & location=SourceLocation( &
 & 'uniformVectorField_test.f90', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"

    ! Initialise field
    call dict % init(2)
    call dict % store('type', 'uniformVectorField')
    call dict % store('value', [9.6_defReal, -8.0_defReal, 9.7_defReal])

    call fieldT % init(dict)

    ! Check value
#line 53 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"
  call assertEqual([9.6_defReal, -8.0_defReal, 9.7_defReal], fieldT % at(p), TOL, &
 & location=SourceLocation( &
 & 'uniformVectorField_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/Geometry/Fields/VectorFields/Tests/uniformVectorField_test.f90"

    ! Kill
    call fieldT % kill()

  end subroutine test_uniformVectorField

end module uniformVectorField_test

module WrapuniformVectorField_test
   use pFUnit_mod
   use uniformVectorField_test
   implicit none
   private

contains


end module WrapuniformVectorField_test

function uniformVectorField_test_suite() result(suite)
   use pFUnit_mod
   use uniformVectorField_test
   use WrapuniformVectorField_test
   type (TestSuite) :: suite

   suite = newTestSuite('uniformVectorField_test_suite')

   call suite%addTest(newTestMethod('test_uniformVectorField', test_uniformVectorField))


end function uniformVectorField_test_suite

