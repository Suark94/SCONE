module conversions_test
  use numPrecision
  use genericProcedures, only : charToInt
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test char to int conversion
  !!
!@Test
  subroutine testCharToInt()
    logical(defBool)  :: flag
    integer(shortInt) :: i
    ! Easy cases
#line 18 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertEqual(2, charToInt('   2 '), &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 18) )
  if (anyExceptions()) return
#line 19 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
#line 19 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertEqual(-1, charToInt('  -1'), &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 19) )
  if (anyExceptions()) return
#line 20 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"

    ! Edge cases
#line 22 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertEqual(7,charToInt('7                   A'), &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 22) )
  if (anyExceptions()) return
#line 23 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"


    ! Cases that should fail
    i = charToInt('7.0', error = flag)
#line 27 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertTrue(flag, &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 27) )
  if (anyExceptions()) return
#line 28 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"

    i = charToInt('7E-7', error = flag)
#line 30 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertTrue(flag, &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"

    i = charToInt('NaN', error = flag)
#line 33 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertTrue(flag, &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"

    i = charToInt('7 is Not A Number For Fortran Rex', error = flag)
#line 36 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertTrue(flag, &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"

    i = charToInt('253.', error = flag)
#line 39 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"
  call assertTrue(flag, &
 & location=SourceLocation( &
 & 'conversions_test.f90', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/home/mhgk4/SCONE/SharedModules/Tests/conversions_test.f90"

  end subroutine testCharToInt

    
end module conversions_test

module Wrapconversions_test
   use pFUnit_mod
   use conversions_test
   implicit none
   private

contains


end module Wrapconversions_test

function conversions_test_suite() result(suite)
   use pFUnit_mod
   use conversions_test
   use Wrapconversions_test
   type (TestSuite) :: suite

   suite = newTestSuite('conversions_test_suite')

   call suite%addTest(newTestMethod('testCharToInt', testCharToInt))


end function conversions_test_suite

