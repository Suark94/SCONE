module fissionCE_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use fissionCE_class,              only : fissionCE, fissionCE_TptrCast
  use aceCard_class,                only : aceCard
  use pFUnit_mod
  implicit none

contains

  !!
  !! Integration test of elasticScattering reaction
  !! Tests:
  !!   -> building from ACE
  !!   -> pointer Casting
  !!   -> probability of Scattering
  !!
  !! Does NOT verify correctness of the sampling procedure
  !!
!@Test
  subroutine testFissionCE()
    type(fissionCE), target               :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(fissionCE),pointer               :: fissionPtr
    type(aceCard)                         :: ACE
    type(RNG)                             :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    fissionPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
#line 41 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertTrue(associated(unCorrPtr, reaction), &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"

    ! Elastic Scattering type cast
    fissionPtr => fissionCE_TptrCast(reaction)
#line 45 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertTrue(associated(fissionPtr, reaction), &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/92233JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_FISSION)

    ! Test trivial functionality
#line 54 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertFalse(reaction % inCMframe(), &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"

    ! Test neutron release
#line 57 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertEqual(2.65431_defReal, reaction % release(1.6_defReal), TOL, &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
#line 58 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertEqual(5.147534_defReal, reaction % release(17.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 58) )
  if (anyExceptions()) return
#line 59 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
#line 59 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertEqual(2.48098_defReal, reaction % releasePrompt(0.6E-6_defReal), TOL, &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 59) )
  if (anyExceptions()) return
#line 60 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
#line 60 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertEqual(6.73E-3_defReal, reaction % releaseDelayed(0.6E-6_defReal), TOL, &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
#line 61 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"
  call assertEqual(0.0041725_defReal, reaction % releaseDelayed(17.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'fissionCE_iTest.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/fissionCE_iTest.f90"

    ! Test probability density
    !@assertEqual(0.1618843E-01_defReal, reaction % probOf(0.7_defReal, 2.0_defReal, 0.1404_defReal, 2.0_defReal), TOL)
    !@assertEqual(0.1586432E-01, reaction % probOf(0.7_defReal, 2.0_defReal, 2.48077_defReal, 14.0_defReal), TOL)

    ! Clean
    call reaction % kill()

  end subroutine testFissionCE

end module fissionCE_iTest

module WrapfissionCE_iTest
   use pFUnit_mod
   use fissionCE_iTest
   implicit none
   private

contains


end module WrapfissionCE_iTest

function fissionCE_iTest_suite() result(suite)
   use pFUnit_mod
   use fissionCE_iTest
   use WrapfissionCE_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('fissionCE_iTest_suite')

   call suite%addTest(newTestMethod('testFissionCE', testFissionCE))


end function fissionCE_iTest_suite

