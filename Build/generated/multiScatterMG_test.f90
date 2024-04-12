module multiScatterMG_test

  use numPrecision
  use endfConstants
  use RNG_class,            only : RNG
  use dictionary_class,     only : dictionary
  use dictDeck_class,       only : dictDeck
  use reactionHandle_inter, only : reactionHandle
  use reactionMG_inter,     only : reactionMG, reactionMG_CptrCast
  use multiScatterMG_class, only : multiScatterMG, multiScatterMG_TptrCast, multiScatterMG_CptrCast
  use pfUnit_mod

  implicit none

  !! Test parameters
  real(defReal),dimension(4),parameter :: P0   = [1.3_defReal, 0.7_defReal, 0.3_defReal, 4.0_defReal ]
  real(defReal),dimension(4),parameter :: prod = [1.1_defReal, 1.05_defReal, ONE, ONE]


contains

  !!
  !! Tests pointer casting and basic (deterministic) functionality
  !! Does NOT test sampling procedure!
  !!
!@Test
  subroutine multiplicative_scattering_test()
    type(multiScatterMG), target   :: reaction
    class(reactionHandle), pointer :: handlePtr
    class(reactionMG), pointer     :: mgPtr
    class(multiScatterMG), pointer :: multiCPtr
    type(multiScatterMG), pointer  :: multiTPtr
    type(dictionary),target       :: dictT
    type(dictDeck)                :: data
    type(RNG)                     :: rand
    real(defReal),parameter :: TOL = 1.0E-6_defReal

    ! Set pointers
    handlePtr => reaction
    mgPtr => null()
    multiCPtr => null()
    multiTPtr => null()

    ! Uncorrelated Reaction class cast
    mgPtr => reactionMG_CptrCast(reaction)
#line 46 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertTrue(associated(mgPtr, reaction), &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 46) )
  if (anyExceptions()) return
#line 47 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Multiplicative Scattering class cast
    multiCPtr => multiScatterMG_CptrCast(reaction)
#line 50 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertTrue(associated(multiCPtr, reaction), &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Multiplicative Scattering type cast
    multiTPtr => multiScatterMG_TptrCast(reaction)
#line 54 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertTrue(associated(multiTPtr, reaction), &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Initialise multiplicative scattering
    call dictT % init(3)
    call dictT % store('numberOfGroups',2)
    call dictT % store('P0',P0)
    call dictT % store('scatteringMultiplicity', prod)
    data % dict => dictT
    call reaction % init(data, anyScatter)

    ! Test misc functionality
#line 65 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ZERO, reaction % releaseDelayed(1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 66 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ZERO, reaction % sampleDelayRate(1, rand), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Test average release
#line 69 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ZERO, reaction % releasePrompt(-1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 70 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ZERO, reaction % release(-1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 70) )
  if (anyExceptions()) return
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 71 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ZERO, reaction % releasePrompt(172), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 72 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ZERO, reaction % release(42), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 73 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ONE, reaction % release(2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 74 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(1.0825_defReal, reaction % releasePrompt(1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Test getting total scattering Xss
#line 77 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(2.0_defReal, reaction % scatterXS(1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 78 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(4.3_defReal, reaction % scatterXS(2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Test getting Group-To-Group production
#line 81 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(1.1_defReal, reaction % production(1,1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 82 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(ONE, reaction % production(2,1), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
#line 83 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"
  call assertEqual(1.05_defReal, reaction % production(1,2), TOL, &
 & location=SourceLocation( &
 & 'multiScatterMG_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/NuclearData/Reactions/reactionMG/Tests/multiScatterMG_test.f90"

    ! Clean memory
    call dictT % kill()
    call reaction % kill()

  end subroutine multiplicative_scattering_test


end module multiScatterMG_test

module WrapmultiScatterMG_test
   use pFUnit_mod
   use multiScatterMG_test
   implicit none
   private

contains


end module WrapmultiScatterMG_test

function multiScatterMG_test_suite() result(suite)
   use pFUnit_mod
   use multiScatterMG_test
   use WrapmultiScatterMG_test
   type (TestSuite) :: suite

   suite = newTestSuite('multiScatterMG_test_suite')

   call suite%addTest(newTestMethod('multiplicative_scattering_test', multiplicative_scattering_test))


end function multiScatterMG_test_suite

