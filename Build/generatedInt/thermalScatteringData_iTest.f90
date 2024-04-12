module thermalScatteringData_iTest

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
  use ceNeutronCache_mod,       only : nuclideCache
  use pFUnit_mod

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR =        &
  & "water { temp 1;                               &
  &       moder {1001.03 h-h2o.49; }               &
  &       composition {                            &
  &       1001.03  2.0E-3;                         &
  &       8016.03  1.0E-3;                         &
  &                   }                            &
  &      }                                         &
  &  graphite { temp 1;                                &
  &          moder {6012.06  grph30.46;}               &
  &          composition {                             &
  &          6012.06 2.0E-3;                           &
  &                       }                            &
  &            }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; "

contains

  !!
  !! Test the use of thermal scattering libraries
  !!
!@Test
  subroutine test_thermalScatteringData()
    type(aceNeutronDatabase), target  :: data
    class(nuclearDatabase), pointer   :: ptr
    type(dictionary)                  :: matDict
    type(dictionary)                  :: dataDict
    class(aceNeutronNuclide), pointer :: H1, O16, C12
    real(defReal)                     :: val
    real(defReal), dimension(2)       :: eBounds
    class(ceNeutronNuclide), pointer  :: nuc
    type(particle)                    :: p
    type(neutronMicroXSs)             :: microXSs
    integer(shortInt)                 :: Nin
    logical(defBool)                  :: gotIt
    real(defReal), parameter          :: TOL = 1.0E-6

    ! Prepare dictionaries
    call charToDict(matDict, MAT_INPUT_STR)
    call charToDict(dataDict, ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)
    call data % activate(([1,2]))

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get nuclides
    O16 => aceNeutronNuclide_CptrCast( data % getNuclide(1))
    H1  => aceNeutronNuclide_CptrCast( data % getNuclide(2))
    C12 => aceNeutronNuclide_CptrCast( data % getNuclide(3))

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test scattering tables

#line 83 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertTrue(H1 % hasThData, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 84 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertFalse(O16 % hasThData, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

#line 86 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertFalse(H1 % thData % hasElastic, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

#line 88 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertTrue(C12 % hasThData, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 89 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertTrue(C12 % thData % hasElastic, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 90 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertTrue(C12 % thData % isCoherent, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test energy bounds
    eBounds = H1 % thData % getEbounds('inelastic')

#line 96 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(1.000E-11_defReal, eBounds(1), TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 97 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(1.000E-5_defReal,  eBounds(2), TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 97) )
  if (anyExceptions()) return
#line 98 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

#line 99 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(O16 % SabInel(1), ZERO, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 100 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(O16 % SabInel(2), ZERO, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 100) )
  if (anyExceptions()) return
#line 101 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

    eBounds = C12 % thData % getEbounds('elastic')

#line 104 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(1.000E-11_defReal, eBounds(1), TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 105 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(4.9000E-06,  eBounds(2), TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test sampling from tables

    val = H1 % thData % getInelXS(1.8E-6_defReal)
#line 111 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(21.018654322_defReal, val, TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

    val = H1 % thData % getElXS(1.8E-6_defReal)
#line 114 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ZERO, val, TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test Getting material XSs
    ! water

    call p % build([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO], 1.0E-6_defReal, ONE)

    ! Total XS of water
    p % E = 1.8E-6_defReal
#line 124 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ONE, data % getTotalMatXS(p , 1) / 0.0459700882_defReal , TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 124) )
  if (anyExceptions()) return
#line 125 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting XSs
    ! H-1
    nuc  => ceNeutronNuclide_CptrCast( data % getNuclide(2))
    nuclideCache(2) % E_tot = ONE
    nuclideCache(2) % needsSabInel = .true.

    call nuc % getMicroXSs(microXSs, 1.8E-6_defReal, p % pRNG)

#line 135 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ONE, 21.05810233858_defReal/ microXSs % total,          TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 135) )
  if (anyExceptions()) return
#line 136 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 136 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ONE, 21.01865432_defReal / microXSs % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 137 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ONE, 3.94480160E-002_defReal / microXSs % capture,      TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 138 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ZERO, microXSs % elasticScatter,   TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 139 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ZERO, microXSs % fission,          TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
#line 140 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"
  call assertEqual(ZERO, microXSs % nuFission,        TOL, &
 & location=SourceLocation( &
 & 'thermalScatteringData_iTest.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/thermalScatteringData_iTest.f90"

  end subroutine test_thermalScatteringData


end module thermalScatteringData_iTest

module WrapthermalScatteringData_iTest
   use pFUnit_mod
   use thermalScatteringData_iTest
   implicit none
   private

contains


end module WrapthermalScatteringData_iTest

function thermalScatteringData_iTest_suite() result(suite)
   use pFUnit_mod
   use thermalScatteringData_iTest
   use WrapthermalScatteringData_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('thermalScatteringData_iTest_suite')

   call suite%addTest(newTestMethod('test_thermalScatteringData', test_thermalScatteringData))


end function thermalScatteringData_iTest_suite

