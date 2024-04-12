module urrProbabilityTables_iTest

  use numPrecision
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use particle_class,           only : particle
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_CptrCast
  use neutronXSPackages_class,  only : neutronMicroXSs
  use materialMenu_mod,         only : mm_init => init
  use ceNeutronCache_mod,       only : zaidCache, nuclideCache
  use pFUnit_mod

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = &
  & " uo2  { temp 1;           &
  &        composition {       &
  &        92235.03 1.0E-3;    &
  &        8016.03  2.0E-3;    &
  &        }                   &
  &      }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; ures 1 ;"

contains

  !!
  !! Test the use of probability tables
  !!
!@Test
  subroutine test_urrProbabilityTables()
    type(aceNeutronDatabase), target  :: data
    class(nuclearDatabase), pointer   :: ptr
    type(dictionary)                  :: matDict
    type(dictionary)                  :: dataDict
    class(aceNeutronNuclide), pointer :: U235, O16
    real(defReal), dimension(3)       :: val
    real(defReal), dimension(2)       :: eBounds
    class(ceNeutronNuclide), pointer  :: nuc
    type(particle)                    :: p
    type(neutronMicroXSs)             :: microXSs
    real(defReal), parameter          :: TOL = 1.0E-6

    ! Prepare dictionaries
    call charToDict(matDict, MAT_INPUT_STR)
    call charToDict(dataDict, ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)
    call data % activate([1])

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get nuclides
    U235  => aceNeutronNuclide_CptrCast( data % getNuclide(1))
    O16   => aceNeutronNuclide_CptrCast( data % getNuclide(2))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test probability tables

#line 73 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertTrue(U235 % hasProbTab, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 74 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertFalse(O16 % hasProbTab, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test energy bounds

    eBounds = U235 % probTab % getEbounds()

#line 81 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(2.25E-3_defReal, eBounds(1), TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 82 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(2.5E-2_defReal,  eBounds(2), TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"

#line 84 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(O16 % urrE(1), ZERO, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 85 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(O16 % urrE(2), ZERO, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test sampling from tables

    call U235 % probTab % sampleXSs(9.1E-3_defReal, 0.347_defReal, val)

#line 92 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(0.98499622_defReal, val(1), TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 93 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(0.83939802_defReal, val(2), TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 93) )
  if (anyExceptions()) return
#line 94 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 94 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(0.8515398_defReal,  val(3), TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 94) )
  if (anyExceptions()) return
#line 95 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting XSs

    ! U-235
    nuc  => ceNeutronNuclide_CptrCast( data % getNuclide(1))
    zaidCache(1) % E = 9.1E-3_defReal
    zaidCache(1) % xi = 0.347_defReal
    nuclideCache(1) % E_tot = ONE
    nuclideCache(1) % needsUrr = .true.

    call nuc % getMicroXSs(microXSs, 9.1E-3_defReal, p % pRNG)

#line 108 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(ONE, 15.317184903738868_defReal/ microXSs % total,            TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 109 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(ONE, 11.662135262310867_defReal/ microXSs % elasticScatter,   TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 110 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(ONE, 0.5743300000E-5_defReal   / microXSs % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 111 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(ONE, 0.999051523404001_defReal / microXSs % capture,          TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 112 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(ONE, 2.655992374724002_defReal / microXSs % fission,          TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
#line 113 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"
  call assertEqual(ONE, 6.462838469906821_defReal / microXSs % nuFission,        TOL, &
 & location=SourceLocation( &
 & 'urrProbabilityTables_iTest.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/urrProbabilityTables_iTest.f90"

  end subroutine test_urrProbabilityTables


end module urrProbabilityTables_iTest

module WrapurrProbabilityTables_iTest
   use pFUnit_mod
   use urrProbabilityTables_iTest
   implicit none
   private

contains


end module WrapurrProbabilityTables_iTest

function urrProbabilityTables_iTest_suite() result(suite)
   use pFUnit_mod
   use urrProbabilityTables_iTest
   use WrapurrProbabilityTables_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('urrProbabilityTables_iTest_suite')

   call suite%addTest(newTestMethod('test_urrProbabilityTables', test_urrProbabilityTables))


end function urrProbabilityTables_iTest_suite

