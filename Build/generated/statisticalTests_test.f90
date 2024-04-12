module statisticalTests_test

  use numPrecision
  use statisticalTests_func, only : twoSampleKS
  use pfUnit_mod

  implicit none

contains

  !!
  !! Test Kolmogorov Smirnov Test on a simple case
  !! Taken from: http://influentialpoints.com/Training/kolmogorov-smirnov_test.htm
  !!
!@Test
  subroutine testTwoSamplesKS()
    real(defReal), dimension(11) :: S1 = [51.0, 71.0, 42.0, 37.0, 51.0, 78.0, 51.0, 49.0, 56.0, 47.0, 58.0]
    real(defReal), dimension(5)  :: S2 = [45.0, 87.0, 123.0, 120.0, 70.0]
    real(defReal), dimension(5)  :: W  = ONE
    real(defReal)                :: p, D

    p = twoSampleKS(S1, S2, W, D)

#line 24 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"
  call assertEqual(0.6182_defReal, D, 0.00005_defReal, &
 & location=SourceLocation( &
 & 'statisticalTests_test.f90', &
 & 24) )
  if (anyExceptions()) return
#line 25 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"
#line 25 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"
  call assertEqual(0.07226_defReal, p, 0.00005_defReal , &
 & location=SourceLocation( &
 & 'statisticalTests_test.f90', &
 & 25) )
  if (anyExceptions()) return
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"

  end subroutine testTwoSamplesKS

  !!
  !! Test Kolmogorov Smirnov Test on a simple case
  !! Modified for weights -> double weight of first element of S2
  !! Taken from: http://influentialpoints.com/Training/kolmogorov-smirnov_test.htm
  !!
!@Test
  subroutine testTwoSamplesKS_withWgt()
    real(defReal), dimension(11) :: S1 = [51.0, 71.0, 42.0, 37.0, 51.0, 78.0, 51.0, 49.0, 56.0, 47.0, 58.0]
    real(defReal), dimension(5)  :: S2 = [45.0, 87.0, 123.0, 120.0, 70.0]
    real(defReal), dimension(5)  :: W  = ONE
    real(defReal)                :: p, D

    W(1) = TWO
    p = twoSampleKS(S1, S2, W, D)

#line 44 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"
  call assertEqual(0.48485_defReal, D, 0.00005_defReal, &
 & location=SourceLocation( &
 & 'statisticalTests_test.f90', &
 & 44) )
  if (anyExceptions()) return
#line 45 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"
#line 45 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"
  call assertEqual(0.22280_defReal, p, 0.00005_defReal , &
 & location=SourceLocation( &
 & 'statisticalTests_test.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/SharedModules/Tests/statisticalTests_test.f90"

  end subroutine testTwoSamplesKS_withWgt


end module statisticalTests_test

module WrapstatisticalTests_test
   use pFUnit_mod
   use statisticalTests_test
   implicit none
   private

contains


end module WrapstatisticalTests_test

function statisticalTests_test_suite() result(suite)
   use pFUnit_mod
   use statisticalTests_test
   use WrapstatisticalTests_test
   type (TestSuite) :: suite

   suite = newTestSuite('statisticalTests_test_suite')

   call suite%addTest(newTestMethod('testTwoSamplesKS', testTwoSamplesKS))

   call suite%addTest(newTestMethod('testTwoSamplesKS_withWgt', testTwoSamplesKS_withWgt))


end function statisticalTests_test_suite

