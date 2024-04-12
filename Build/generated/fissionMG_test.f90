module fissionMG_test

  use numPrecision
  use endfConstants
  use RNG_class,            only : RNG
  use dictionary_class,     only : dictionary
  use dictDeck_class,       only : dictDeck
  use reactionHandle_inter, only : reactionHandle
  use reactionMG_inter,     only : reactionMG, reactionMG_CptrCast
  use fissionMG_class,      only : fissionMG, fissionMG_TptrCast
  use pfUnit_mod

  implicit none

  !!
  !! Test data
  !!
  real(defReal),dimension(*),parameter :: nu = [2.3_defReal, 2.0_defReal, 1.3_defReal]
  real(defReal),dimension(*),parameter :: chi = [0.333333_defReal, 0.333333_defReal, 0.333334_defReal]



contains

  !!
  !! Tests pointer casting and basic (deterministic) functionality
  !! Does NOT test sampling procedure!
  !!
!@Test
  subroutine fissionMG_Build_And_Functionality()
    type(fissionMG), target       :: reaction
    class(reactionHandle),pointer :: handlePtr
    class(reactionMG),pointer     :: mgPtr
    type(fissionMG),pointer       :: fissPtr
    type(dictionary),target       :: dictT
    type(dictionary),pointer      :: dictPtr
    type(dictDeck)                :: data
    type(RNG)                     :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    mgPtr => null()
    fissPtr => null()

    ! Uncorrelated Reaction cast
    mgPtr => reactionMG_CptrCast(reaction)
#line 48 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertTrue(associated(mgPtr, reaction), &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"

    ! Elastic Scattering type cast
    fissPtr => fissionMG_TptrCast(reaction)
#line 52 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertTrue(associated(fissPtr, reaction), &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 52) )
  if (anyExceptions()) return
#line 53 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"

    ! Build dictionary for input
    call dictT % init(3)
    call dictT % store('numberOfGroups',3)
    call dictT % store('chi', chi)
    call dictT % store('nu',nu)

    ! Build data Deck and initialise
    data % dict => dictT
    call reaction % init(data, macroFission)

    ! Test Misc functionality
#line 65 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertEqual(ZERO, reaction % releaseDelayed(1), TOL, &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
#line 66 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertEqual(ZERO, reaction % sampleDelayRate(2, rand),TOL , &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"

    ! Test Release
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertEqual(ZERO, reaction % releasePrompt(-2), TOL, &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertEqual(ZERO, reaction % release(6), TOL, &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertEqual(2.0_defReal, reaction % release(2), TOL, &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"
  call assertEqual(1.3_defReal,reaction % releasePrompt(3), TOL, &
 & location=SourceLocation( &
 & 'fissionMG_test.f90', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/fissionMG_test.f90"

    ! Clean
    call dictT % kill()
    call reaction % kill()

  end subroutine fissionMG_Build_And_Functionality
    
end module fissionMG_test

module WrapfissionMG_test
   use pFUnit_mod
   use fissionMG_test
   implicit none
   private

contains


end module WrapfissionMG_test

function fissionMG_test_suite() result(suite)
   use pFUnit_mod
   use fissionMG_test
   use WrapfissionMG_test
   type (TestSuite) :: suite

   suite = newTestSuite('fissionMG_test_suite')

   call suite%addTest(newTestMethod('fissionMG_Build_And_Functionality', fissionMG_Build_And_Functionality))


end function fissionMG_test_suite

