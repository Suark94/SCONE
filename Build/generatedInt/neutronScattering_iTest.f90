module neutronScattering_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use neutronScatter_class,         only : neutronScatter, neutronScatter_TptrCast
  use aceCard_class,                only : aceCard

  use pFUnit_mod
  implicit none


contains

  !!
  !! Integration test of elasticScattering reaction
  !! Tests:
  !!   -> building from ACE
  !!   -> pointer Casting
  !!
  !! There were some problems to establish same normalisation wrt JANIS
  !! Leave probability test for now!
  !!
  !! Does NOT verify correctness of the sampling procedure
  !!
!@Test
  subroutine testNeutronScatteringReaction()
    type(neutronScatter), target          :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(neutronScatter),pointer          :: scatterPtr
    type(aceCard)                         :: ACE
    type(RNG)                             :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    scatterPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
#line 45 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
  call assertTrue(associated(unCorrPtr, reaction), &
 & location=SourceLocation( &
 & 'neutronScattering_iTest.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"

    ! Elastic Scattering type cast
    scatterPtr => neutronScatter_TptrCast(reaction)
#line 49 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
  call assertTrue(associated(scatterPtr, reaction), &
 & location=SourceLocation( &
 & 'neutronScattering_iTest.f90', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/8016JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_2N)

    ! Verify simple functionality
#line 58 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
  call assertTrue(reaction % inCMFrame(), &
 & location=SourceLocation( &
 & 'neutronScattering_iTest.f90', &
 & 58) )
  if (anyExceptions()) return
#line 59 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
#line 59 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
  call assertEqual(TWO,  reaction % release(7.0_defReal), &
 & location=SourceLocation( &
 & 'neutronScattering_iTest.f90', &
 & 59) )
  if (anyExceptions()) return
#line 60 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
#line 60 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
  call assertEqual(TWO,  reaction % releasePrompt(1.0E-6_defReal), &
 & location=SourceLocation( &
 & 'neutronScattering_iTest.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
#line 61 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"
  call assertEqual(ZERO, reaction % releaseDelayed(1.5E-3_defReal), &
 & location=SourceLocation( &
 & 'neutronScattering_iTest.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/neutronScattering_iTest.f90"

    ! Clean
    call reaction % kill()

  end subroutine testNeutronScatteringReaction


end module neutronScattering_iTest

module WrapneutronScattering_iTest
   use pFUnit_mod
   use neutronScattering_iTest
   implicit none
   private

contains


end module WrapneutronScattering_iTest

function neutronScattering_iTest_suite() result(suite)
   use pFUnit_mod
   use neutronScattering_iTest
   use WrapneutronScattering_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('neutronScattering_iTest_suite')

   call suite%addTest(newTestMethod('testNeutronScatteringReaction', testNeutronScatteringReaction))


end function neutronScattering_iTest_suite

