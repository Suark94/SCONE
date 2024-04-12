module linearAlgebra_test
  use numPrecision
  use linearAlgebra_func, only : kill_linearAlgebra, eig, solve, solveAdjointProblem
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test eigenvalue calculation
  !!
!@Test
  subroutine testEigenvalue()
    real(defReal), dimension(3,3) :: mat33, v33
    real(defReal), dimension(4,4) :: mat44, v44
    real(defReal), dimension(3)   :: k3
    real(defReal), dimension(4)   :: k4
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Calculate eigenvectors of 3x3 Matrix
    mat33(1,:) = [ 2.0_defReal,  4.0_defReal, -2.0_defReal]
    mat33(2,:) = [ 2.0_defReal, -1.0_defReal, -2.0_defReal]
    mat33(3,:) = [-4.0_defReal, -2.0_defReal, -5.0_defReal]

    call eig(k3, v33, mat33)

    ! Check eigenvalues
#line 29 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual( 5.0_defReal, k3(1), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 29) )
  if (anyExceptions()) return
#line 30 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 30 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual(-3.0_defReal, k3(2), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 31 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual(-6.0_defReal, k3(3), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 31) )
  if (anyExceptions()) return
#line 32 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"

    ! Check eigenvectors
#line 34 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([-0.81649658_defReal,-0.40824829_defReal, 0.40824829_defReal], v33(:,1), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 35 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([-0.53452248_defReal, 0.80178372_defReal, 0.26726124_defReal], v33(:,2), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 35) )
  if (anyExceptions()) return
#line 36 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 36 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([ 0.05842062_defReal, 0.35052374_defReal, 0.93472998_defReal], v33(:,3), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"

    ! Calculate eigenvectors of 4x4 Matrix
    mat44(1,:) = [16.0_defReal,  2.0_defReal,  3.0_defReal, 13.0_defReal]
    mat44(2,:) = [ 5.0_defReal, 11.0_defReal, 10.0_defReal,  8.0_defReal]
    mat44(3,:) = [ 9.0_defReal,  7.0_defReal,  6.0_defReal, 12.0_defReal]
    mat44(4,:) = [ 4.0_defReal, 14.0_defReal, 15.0_defReal,  3.0_defReal]

    call eig(k4, v44, mat44)

    ! Check eigenvalues
#line 47 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual(34.51699899_defReal, k4(1), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 48 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual( 8.79398337_defReal, k4(2), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 49 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual(-7.43156094_defReal, k4(3), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 50 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual( 0.12057858_defReal, k4(4), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"

    ! Check eigenvectors
#line 53 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([-0.49616808_defReal, -0.49173210_defReal, -0.49480677_defReal, -0.51689675_defReal], v44(:,1), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 54 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([-0.83035254_defReal,  0.39013716_defReal,  0.01092828_defReal,  0.39772885_defReal], v44(:,2), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 55 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([ 0.38548337_defReal, -0.00693389_defReal,  0.45994233_defReal, -0.79987971_defReal], v44(:,3), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
#line 56 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([-0.25724392_defReal, -0.66999974_defReal,  0.64227835_defReal,  0.26908070_defReal], v44(:,4), TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 56) )
  if (anyExceptions()) return
#line 57 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"

    ! Kill workspace
    call kill_linearAlgebra()

  end subroutine testEigenvalue

  !!
  !! Test linear solve of Ax=b
  !!
!@Test
  subroutine testLinSolve()
    real(defReal),dimension(3,3) :: A
    real(defReal), dimension(3)  :: b
    real(defReal), dimension(3)  :: x
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set linear system
    A(1,:) = [8.0_defReal, 1.0_defReal, 6.0_defReal]
    A(2,:) = [3.0_defReal, 5.0_defReal, 7.0_defReal]
    A(3,:) = [4.0_defReal, 9.0_defReal, 2.0_defReal]

    b = [ 3.0_defReal, 2.0_defReal, 1.0_defReal]

    ! Solve equation
    call solve(A,x,b)

    ! Verify results
#line 84 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([0.21666667_defReal, -0.03333333_defReal, 0.21666667_defReal],x,TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"

  end subroutine testLinSolve

  !!
  !! Test solution of a generalised adjoint system
  !!
!@Test
  subroutine testAdjointSolve()
    real(defReal),dimension(2,2) :: A
    real(defReal), dimension(2)  :: s
    real(defReal), dimension(2)  :: x
    real(defReal), dimension(2)  :: f
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set linear system
    A(1,:) = [0.9_defReal,  0.1_defReal]
    A(2,:) = [0.1_defReal, 0.9_defReal]

    s = [ONE, -ONE] / sqrt(TWO)
    f = [ONE, ONE] / sqrt(TWO)

    call solveAdjointProblem(A, x, s, f)

    ! Verify solution
#line 109 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"
  call assertEqual([ONE/sqrt(TWO), -ONE/sqrt(TWO)], x, TOL, &
 & location=SourceLocation( &
 & 'linearAlgebra_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/LinearAlgebra/Tests/linearAlgebra_test.f90"

  end subroutine testAdjointSolve

end module linearAlgebra_test

module WraplinearAlgebra_test
   use pFUnit_mod
   use linearAlgebra_test
   implicit none
   private

contains


end module WraplinearAlgebra_test

function linearAlgebra_test_suite() result(suite)
   use pFUnit_mod
   use linearAlgebra_test
   use WraplinearAlgebra_test
   type (TestSuite) :: suite

   suite = newTestSuite('linearAlgebra_test_suite')

   call suite%addTest(newTestMethod('testEigenvalue', testEigenvalue))

   call suite%addTest(newTestMethod('testLinSolve', testLinSolve))

   call suite%addTest(newTestMethod('testAdjointSolve', testAdjointSolve))


end function linearAlgebra_test_suite

