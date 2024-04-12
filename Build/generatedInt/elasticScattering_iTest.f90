module elasticScattering_iTest

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use reactionHandle_inter,         only : reactionHandle
  use uncorrelatedReactionCE_inter, only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use elasticNeutronScatter_class,  only : elasticNeutronScatter, elasticNeutronScatter_TptrCast
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
  subroutine testElasticNeutronScatteringReaction()
    type(elasticNeutronScatter), target   :: reaction
    class(reactionHandle),pointer         :: handlePtr
    class(uncorrelatedReactionCE),pointer :: unCorrPtr
    type(elasticNeutronScatter),pointer   :: elasticScatterPtr
    type(aceCard)                         :: ACE
    type(RNG)                             :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    unCorrPtr => null()
    elasticScatterPtr => null()

    ! Uncorrelated Reaction cast
    unCorrPtr => uncorrelatedReactionCE_CptrCast(reaction)
#line 41 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertTrue(associated(unCorrPtr, reaction), &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! Elastic Scattering type cast
    elasticScatterPtr => elasticNeutronScatter_TptrCast(reaction)
#line 45 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertTrue(associated(elasticScatterPtr, reaction), &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/8016JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_N_ELASTIC)

    ! Verify simple functionality
#line 54 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertTrue(reaction % inCMFrame(), &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 55 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ONE,  reaction % release(7.0_defReal), &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 56 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ONE,  reaction % releasePrompt(1.0E-6_defReal), &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 56) )
  if (anyExceptions()) return
#line 57 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 57 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % releaseDelayed(1.5E-3_defReal), &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! Varify probability distributions
    ! Isotropic range
#line 61 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ONE/TWO/TWO_PI, reaction % probOf(0.5_defReal, 2.0_defReal, 0.6E-6_defReal, 0.6E-6_defReal), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! High energy
    ! It is based on regression really but verified against JANIS polynomial representation
    ! with error of about 1.0E-3 absolute, which is accaptable
    ! It was immposible to access polynomial representation directly, and there might be additional
    ! Issues due to tabular representation of the data.
#line 68 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(0.4273399E-01_defReal, reaction % probOf(0.7_defReal, 2.0_defReal, 0.36_defReal, 0.36_defReal), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(0.1593347_defReal, reaction % probOf(0.75911_defReal, 2.0_defReal, 3.94_defReal, 3.94_defReal), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(0.7943654E-01_defReal, reaction % probOf(-0.3_defReal, 2.0_defReal, 8.9201_defReal, 8.9201_defReal), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! Test invalid angle ranges
#line 73 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % probOf(1.1_defReal, 2.0_defReal, ONE, ONE), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 74 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % probOf(0.7_defReal, -2.0_defReal, ONE, ONE), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 75 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % probOf(0.7_defReal, 2.0_defReal, -ONE, -ONE), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 75) )
  if (anyExceptions()) return
#line 76 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! Clean
    call reaction % kill()

  end subroutine testElasticNeutronScatteringReaction

  !!
  !! Test elasticScattering for a data with LOCB == 0
  !! Use Ta-126 from JEFF 3.1.1
  !!
  !! Does not test sampling
  !!
!@Test
  subroutine testElasticNeutronScattering_isotropic()
    type(elasticNeutronScatter)  :: reaction
    type(aceCard)                :: ACE
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/52126JEF311.ace', 1)

    ! Build reaction object
    call reaction % init(ACE, N_N_ELASTIC)

    ! Test probability distribution
#line 101 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ONE/TWO/TWO_PI, reaction % probOf(0.5_defReal, 2.0_defReal, 6.7_defReal, 6.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

    ! Test invalid angle aranges
#line 104 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % probOf(1.1_defReal, 2.0_defReal, ONE, ONE), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 105 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % probOf(0.7_defReal, -2.0_defReal, ONE, ONE), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
#line 106 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"
  call assertEqual(ZERO, reaction % probOf(0.7_defReal, 2.0_defReal, -ONE, -ONE), TOL, &
 & location=SourceLocation( &
 & 'elasticScattering_iTest.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/NuclearData/Reactions/Tests/elasticScattering_iTest.f90"

  end subroutine testElasticNeutronScattering_isotropic

end module elasticScattering_iTest

module WrapelasticScattering_iTest
   use pFUnit_mod
   use elasticScattering_iTest
   implicit none
   private

contains


end module WrapelasticScattering_iTest

function elasticScattering_iTest_suite() result(suite)
   use pFUnit_mod
   use elasticScattering_iTest
   use WrapelasticScattering_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('elasticScattering_iTest_suite')

   call suite%addTest(newTestMethod('testElasticNeutronScatteringReaction', testElasticNeutronScatteringReaction))

   call suite%addTest(newTestMethod('testElasticNeutronScattering_isotropic', testElasticNeutronScattering_isotropic))


end function elasticScattering_iTest_suite

