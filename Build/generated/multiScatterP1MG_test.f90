module multiScatterP1MG_test

  use numPrecision
  use endfConstants
  use RNG_class,              only : RNG
  use dictionary_class,       only : dictionary
  use dictDeck_class,         only : dictDeck
  use reactionHandle_inter,   only : reactionHandle
  use reactionMG_inter,       only : reactionMG, reactionMG_CptrCast
  use multiScatterP1MG_class, only : multiScatterP1MG, multiScatterP1MG_TptrCast
  use pfUnit_mod

  implicit none

  !! Test parameters
  real(defReal),dimension(4),parameter :: P0   = [1.3_defReal, 0.7_defReal, 0.3_defReal, 4.0_defReal ]
  real(defReal),dimension(4),parameter :: P1   = [0.5_defReal, 0.1_defReal, ZERO, ZERO ]
  real(defReal),dimension(4),parameter :: prod = [1.1_defReal, 1.05_defReal, ONE, ONE]

contains

  !!
  !! Tests pointer casting and basic (deterministic) functionality
  !! Does NOT test sampling procedure!
  !!
  !! NOTE: Mostly copy of test for multiScatterMG. Significant functionality is
  !!   shared, but it is not necessarly so. Thus test is duplicated.
  !!
!@Test
  subroutine multiplicative_p1_scattering_test()
    type(multiScatterP1MG), target   :: reaction
    class(reactionHandle), pointer   :: handlePtr
    class(reactionMG), pointer       :: mgPtr
    type(multiScatterP1MG), pointer  :: multiTPtr
    type(dictionary),target          :: dictT
    type(dictDeck)                   :: data
    type(RNG)                        :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    mgPtr => null()
    multiTPtr => null()

    ! Uncorrelated Reaction class cast
    mgPtr => reactionMG_CptrCast(reaction)
#line 47 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertTrue(associated(mgPtr, reaction), &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Multiplicative Scattering type cast
    multiTPtr => multiScatterP1MG_TptrCast(reaction)
#line 51 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertTrue(associated(multiTPtr, reaction), &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 51) )
  if (anyExceptions()) return
#line 52 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Initialise multiplicative scattering
    call dictT % init(4)
    call dictT % store('numberOfGroups',2)
    call dictT % store('P0',P0)
    call dictT % store('P1',P1)
    call dictT % store('scatteringMultiplicity', prod)
    data % dict => dictT
    call reaction % init(data, anyScatter)

    ! Test misc functionality
#line 63 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % releaseDelayed(1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 64 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % sampleDelayRate(1, rand), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 64) )
  if (anyExceptions()) return
#line 65 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Test average release
#line 67 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % releasePrompt(-1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 68 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % release(-1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % releasePrompt(172), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % release(42), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ONE, reaction % release(2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(1.0825_defReal, reaction % releasePrompt(1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Test getting total scattering Xss
#line 75 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(2.0_defReal, reaction % scatterXS(1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 75) )
  if (anyExceptions()) return
#line 76 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 76 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(4.3_defReal, reaction % scatterXS(2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Test getting Group-To-Group production
#line 79 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(1.1_defReal, reaction % production(1,1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 79) )
  if (anyExceptions()) return
#line 80 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 80 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ONE, reaction % production(2,1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 81 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(1.05_defReal, reaction % production(1,2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Test that normalisation of P1 coefficients is OK
#line 84 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % P1(2,2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 85 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(ZERO, reaction % P1(1,2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 86 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(1.1538461538_defReal, reaction % P1(1,1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
#line 87 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"
  call assertEqual(0.4285714287_defReal, reaction % P1(2,1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterP1MG_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterP1MG_test.f90"

    ! Clean memory
    call dictT % kill()
    call reaction % kill()

  end subroutine multiplicative_p1_scattering_test

end module multiScatterP1MG_test

module WrapmultiScatterP1MG_test
   use pFUnit_mod
   use multiScatterP1MG_test
   implicit none
   private

contains


end module WrapmultiScatterP1MG_test

function multiScatterP1MG_test_suite() result(suite)
   use pFUnit_mod
   use multiScatterP1MG_test
   use WrapmultiScatterP1MG_test
   type (TestSuite) :: suite

   suite = newTestSuite('multiScatterP1MG_test_suite')

   call suite%addTest(newTestMethod('multiplicative_p1_scattering_test', multiplicative_p1_scattering_test))


end function multiScatterP1MG_test_suite

