module pureCapture_test

  use numPrecision
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use pureCapture_class,            only : pureCapture, pureCapture_TptrCast
  use dictDeck_class,               only : dictDeck

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
  subroutine testPureCaptureReaction()
    type(pureCapture), target             :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(pureCapture),pointer             :: pureCapPtr
    type(dictDeck)                        :: fakeDeck
    type(RNG)                             :: rand
    real(defReal)                         :: E_out, mu, phi
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    pureCapPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
#line 44 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertTrue(associated(unCorrPtr, reaction), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 44) )
  if (anyExceptions()) return
#line 45 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"

    ! Elastic Scattering type cast
    pureCapPtr => pureCapture_TptrCast(reaction)
#line 48 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertTrue(associated(pureCapPtr, reaction), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"

    ! Initialise -> not needed really
    call pureCapPtr % init(fakeDeck, 2)

    ! Test functionality
    ! Release & Delay rate sampling
#line 55 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, pureCapPtr % release(1.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 56 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, pureCapPtr % releasePrompt(2.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 56) )
  if (anyExceptions()) return
#line 57 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 57 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, pureCapPtr % releaseDelayed(2.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 58 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertFalse(pureCapPtr % inCMframe(), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 58) )
  if (anyExceptions()) return
#line 59 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"

    ! Verify probability distribution
#line 61 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ONE, pureCapPtr % probOf(ONE, ZERO, 1.0_defReal, 1.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 62 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ONE, pureCapPtr % probOf(ONE, TWO_PI, 1.0_defReal, 1.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 63 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, pureCapPtr % probOf(ONE, ZERO, 1.1_defReal, 1.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 64 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, pureCapPtr % probOf(ONE, 1.0_defReal, 1.0_defReal, 1.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 64) )
  if (anyExceptions()) return
#line 65 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 65 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, pureCapPtr % probOf(0.0_defReal, ZERO, 1.0_defReal, 1.0_defReal), &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"

    ! Sample outgoing
    call pureCapPtr % sampleOut(mu, phi, E_out, 8.0_defReal, rand)
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ONE, mu, TOL, &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(ZERO, phi, TOL, &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"
  call assertEqual(8.0_defReal, E_out, TOL, &
 & location=SourceLocation( &
 & 'pureCapture_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/pureCapture_test.f90"

    ! Clean
    call reaction % kill()

  end subroutine testPureCaptureReaction


end module pureCapture_test

module WrappureCapture_test
   use pFUnit_mod
   use pureCapture_test
   implicit none
   private

contains


end module WrappureCapture_test

function pureCapture_test_suite() result(suite)
   use pFUnit_mod
   use pureCapture_test
   use WrappureCapture_test
   type (TestSuite) :: suite

   suite = newTestSuite('pureCapture_test_suite')

   call suite%addTest(newTestMethod('testPureCaptureReaction', testPureCaptureReaction))


end function pureCapture_test_suite

