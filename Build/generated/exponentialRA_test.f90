module exponentialRA_test

  use numPrecision
  use exponentialRA_func, only : exponential
  use pfUnit_mod

  implicit none

contains

  !!
  !! Test exponential rational approximation on a few values
  !! ExponentialRA evaluated (1 - exp(-x))
  !! This is compared against the analytic equivalent
  !!
!@Test
  subroutine testExponentialRA()
    real(defFlt)                :: x, res, resRA
    real(defFlt), parameter     :: tol = 1E-5

    x = 0.5
    res   = ONE - exp(-x)
    resRA = exponential(x) 

#line 25 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"
  call assertEqual(res, resRA, tol, &
 & location=SourceLocation( &
 & 'exponentialRA_test.f90', &
 & 25) )
  if (anyExceptions()) return
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"

    x = 0.2_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

#line 31 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"
  call assertEqual(res, resRA, tol, &
 & location=SourceLocation( &
 & 'exponentialRA_test.f90', &
 & 31) )
  if (anyExceptions()) return
#line 32 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"

    x = 0.03_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

#line 37 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"
  call assertEqual(res, resRA, tol, &
 & location=SourceLocation( &
 & 'exponentialRA_test.f90', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"

    x = 3.0_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

#line 43 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"
  call assertEqual(res, resRA, tol, &
 & location=SourceLocation( &
 & 'exponentialRA_test.f90', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"


    x = 0.0001_defFlt
    res   = ONE - exp(-x)
    resRA = exponential(x) 

#line 50 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"
  call assertEqual(res, resRA, tol, &
 & location=SourceLocation( &
 & 'exponentialRA_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/SharedModules/Tests/exponentialRA_test.f90"
  end subroutine testExponentialRA


end module exponentialRA_test

module WrapexponentialRA_test
   use pFUnit_mod
   use exponentialRA_test
   implicit none
   private

contains


end module WrapexponentialRA_test

function exponentialRA_test_suite() result(suite)
   use pFUnit_mod
   use exponentialRA_test
   use WrapexponentialRA_test
   type (TestSuite) :: suite

   suite = newTestSuite('exponentialRA_test_suite')

   call suite%addTest(newTestMethod('testExponentialRA', testExponentialRA))


end function exponentialRA_test_suite

