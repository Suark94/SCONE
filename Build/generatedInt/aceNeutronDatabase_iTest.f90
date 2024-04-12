module aceNeutronDatabase_iTest

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,        only : linFind
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use charMap_class,            only : charMap
  use particle_class,           only : particle
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use materialHandle_inter,     only : materialHandle
  use nuclideHandle_inter,      only : nuclideHandle
  use reactionHandle_inter,     only : reactionHandle
  use ceNeutronMaterial_class,  only : ceNeutronMaterial, ceNeutronMaterial_TptrCast
  use ceNeutronNuclide_inter,   only : ceNeutronNuclide, ceNeutronNuclide_CptrCast
  use aceNeutronNuclide_class,  only : aceNeutronNuclide, aceNeutronNuclide_TptrCast
  use neutronXSPackages_class,  only : neutronMicroXSs, neutronMacroXSs
  use materialMenu_mod,         only : mm_init => init, mm_kill => kill
  use pFUnit_mod

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = &
  & "water { temp 273;           &
  &       composition {          &
  &       1001.03 5.028E-02;     &
  &       8016.03 2.505E-02;     &
  &                   }          &
  &        }                     &
  &  uo2  { temp 1;              &
  &        composition {         &
  &        92233.03 2.286E-02;   &
  &        8016.03  4.572E-02;   &
  &                    }         &
  &       }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = &
  & "aceLibrary ./IntegrationTestFiles/testLib; "

contains

  !!
  !! One big monster test to avoid expensive set up each test
  !!
!@Test
  subroutine test_aceNeutronDatabase()
    type(aceNeutronDatabase), target :: data
    class(nuclearDatabase), pointer  :: ptr
    type(dictionary)                 :: matDict
    type(dictionary)                 :: dataDict
    type(ceNeutronMaterial),pointer  :: mat
    class(ceNeutronNuclide), pointer :: nuc
    type(aceNeutronNuclide), pointer :: nuc2
    class(reactionHandle), pointer   :: reac
    type(charMap), pointer           :: matNames
    character(nameLen)               :: name
    type(particle)                   :: p
    type(neutronMicroXSs)            :: microXSs
    type(neutronMacroXSs)            :: macroXSs
    real(defReal)                    :: t1, t2
    integer(shortInt)                :: i, H1, O16, U233
    real(defReal), parameter         :: TOL = 1.0E-6

    ! Prepare dictionaries
    call charToDict(matDict, MAT_INPUT_STR)
    call charToDict(dataDict, ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)
    call data % activate([1,2])

    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !! Perform tests
    !!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    !!<><><><><><><><><><><><><><><><><><><><><><
    !! Test getting material
    mat => null()

    ! Get invalid materials
#line 89 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( ceNeutronMaterial_TptrCast( data % getMaterial(0))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 90 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( ceNeutronMaterial_TptrCast( data % getMaterial(-4))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 91 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( ceNeutronMaterial_TptrCast( data % getMaterial(3))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 91) )
  if (anyExceptions()) return
#line 92 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Get water
    mat => ceNeutronMaterial_TptrCast( data % getMaterial(1))
#line 95 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated(mat), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Make sure densities are present
#line 98 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(targetNotFound /= linFind(mat % dens, 5.028E-02_defReal, TOL), "H-1 dens is absent", &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 99 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(targetNotFound /= linFind(mat % dens, 2.505E-02_defReal, TOL), "O-16 dens is absent", &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

#line 101 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(1, mat % matIdx, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 102 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse( mat % isFissile(), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 103 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated(mat % data, data), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 103) )
  if (anyExceptions()) return
#line 104 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Get UO2
    mat => ceNeutronMaterial_TptrCast( data % getMaterial(2))
#line 107 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated(mat), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

#line 109 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(targetNotFound /= linFind(mat % dens, 2.286E-02_defReal, TOL), "U-233 dens is absent", &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 110 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(targetNotFound /= linFind(mat % dens, 4.572E-02_defReal, TOL), "O-16 dens is absent", &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

#line 112 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(2, mat % matIdx, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 113 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue( mat % isFissile(), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 114 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated(mat % data, data), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting nuclides
    nuc => null()

    ! Get invalid Nuclides
#line 121 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( ceNeutronNuclide_CptrCast( data % getNuclide(0))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 122 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( ceNeutronNuclide_CptrCast( data % getNuclide(4))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 123 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( ceNeutronNuclide_CptrCast( data % getNuclide(-3))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Get Valid Nuclides
#line 126 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated( ceNeutronNuclide_CptrCast( data % getNuclide(1))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 127 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated( ceNeutronNuclide_CptrCast( data % getNuclide(2))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 128 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated( ceNeutronNuclide_CptrCast( data % getNuclide(3))), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting reactions
    reac => null()


    ! Nuclides can have diffrent indexes if Hashes do not work correctly
    ! Need to explicitly find which index correcponds to which nuclide
    ! Usually will be the following
    ! Nuclides 1 -> O-16
    !          2 -> U-233
    !          3 -> H-1
    do i=1,3
      nuc2 => aceNeutronNuclide_TptrCast( data % getNuclide(i))
      select case(trim(adjustl(nuc2 % ZAID)))
        case('1001.03c')
          H1 = i
        case('8016.03c')
          O16 = i
        case('92233.03c')
          U233 = i
      end select
    end do

    ! Get Invalid Reaction
#line 154 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( data % getReaction(N_Nl(3), H1)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 155 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( data % getReaction(macroEscatter, O16)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 156 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( data % getReaction(N_total, U233)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 157 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( data % getReaction(N_N_elastic, 4)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 158 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( data % getReaction(N_N_elastic, -2)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 159 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertFalse(associated( data % getReaction(N_N_elastic, 0)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Get Valid Reaction
#line 162 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated( data % getReaction(N_fission, U233)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 163 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated( data % getReaction(N_N_elastic, O16)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 163) )
  if (anyExceptions()) return
#line 164 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 164 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue(associated( data % getReaction(N_NL(3), U233)), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 164) )
  if (anyExceptions()) return
#line 165 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Get material names dictionary
    matNames => data % matNamesMap()
    name = 'water'
#line 170 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue( 0 /= matNames % getOrDefault(name, 0), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    name = 'uo2'
#line 173 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertTrue( 0 /= matNames % getOrDefault(name, 0), &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 173) )
  if (anyExceptions()) return
#line 174 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"


    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test getting nuclide XSs
     !
    ! SET UP as a regression test! Take values with the grain of salt!
    !

    ! H-1
    nuc  => ceNeutronNuclide_CptrCast( data % getNuclide(H1))
#line 184 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 20.765855864000002_defReal/ nuc % getTotalXS(1.1E-6_defReal, p % pRNG), TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    call nuc % getMicroXSs(microXSs, 5.6E-3_defReal, p % pRNG)

    ! Absent XSs
#line 189 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, microXSs % fission, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 190 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, microXSs % nuFission, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 190) )
  if (anyExceptions()) return
#line 191 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 191 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, microXSs % inelasticScatter, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Present XSs
#line 194 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 19.731020820000000_defReal   / microXSs % total,          TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 194) )
  if (anyExceptions()) return
#line 195 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 195 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 19.730326000000000_defReal   / microXSs % elasticScatter, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 195) )
  if (anyExceptions()) return
#line 196 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 196 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 6.948036800000000e-04_defReal/ microXSs % capture,        TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 196) )
  if (anyExceptions()) return
#line 197 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><><><><><>
    ! Test Getting material XSs
    !

    call p % build([ZERO, ZERO, ZERO], [ONE, ZERO, ZERO], 1.0E-6_defReal, ONE)

    ! Total XS of water
    p % E = 1.1E-6_defReal
#line 206 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, data % getTotalMatXS(p , 1)/1.1406745607419302_defReal , TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 206) )
  if (anyExceptions()) return
#line 207 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    p % E = 19.9_defReal
#line 209 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, data % getTransMatXS(p , 1)/6.539039844E-02_defReal , TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"


    ! Total XS of UO2
    p % E = 1.1E-6_defReal
#line 214 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, data % getTotalMatXS(p , 2)/4.4149556129495560_defReal , TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 214) )
  if (anyExceptions()) return
#line 215 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    p % E = 19.9_defReal
#line 217 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, data % getTransMatXS(p , 2)/0.21869599644_defReal , TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Majorant
    p % E = 1.1E-6_defReal
#line 221 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, data % getMajorantXS(p) /4.4149556129495560_defReal , TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    p % E = 19.9_defReal
#line 224 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, data % getMajorantXS(p)/0.21869599644_defReal , TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 224) )
  if (anyExceptions()) return
#line 225 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    !<><><><><><><><><><><><><><><><><><><><>
    ! Test getting Macroscopic XSs
    !
    ! Water
    mat => ceNeutronMaterial_TptrCast( data % getMaterial(1))
    call mat % getMacroXSs(macroXss, 3.6E-1_defReal, p % pRNG)

    ! Absent XSs
#line 234 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, macroXSs % fission, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 234) )
  if (anyExceptions()) return
#line 235 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 235 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, macroXSs % nuFission, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

#line 237 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 0.466713100775700_defReal/ macroXSs     % total, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 238 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 0.466710902790000_defReal/ macroXSs     % elasticScatter, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 238) )
  if (anyExceptions()) return
#line 239 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 239 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, macroXSs % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 239) )
  if (anyExceptions()) return
#line 240 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 240 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, 2.198066842597500e-06_defReal/ macroXSs % capture, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Water with some inelastic collisions
    call mat % getMacroXSs(macroXss, 6.525_defReal, p % pRNG)

#line 245 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(ONE, macroXSs % inelasticScatter/1.903667536E-04_defReal, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"


    !<><><><><><><><><><><><><><><><><><><><>
    ! Test getting energy bounds
    !
    call data % energyBounds(t1,t2)
#line 252 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(1.0E-11_defReal, t1, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 252) )
  if (anyExceptions()) return
#line 253 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
#line 253 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"
  call assertEqual(20.0_defReal,    t2, TOL, &
 & location=SourceLocation( &
 & 'aceNeutronDatabase_iTest.f90', &
 & 253) )
  if (anyExceptions()) return
#line 254 "/home/mhgk4/SCONE/NuclearData/ceNeutronData/aceDatabase/Tests/aceNeutronDatabase_iTest.f90"

    ! Clean everything
    call data % kill()
    call mm_kill()

  end subroutine test_aceNeutronDatabase

end module aceNeutronDatabase_iTest

module WrapaceNeutronDatabase_iTest
   use pFUnit_mod
   use aceNeutronDatabase_iTest
   implicit none
   private

contains


end module WrapaceNeutronDatabase_iTest

function aceNeutronDatabase_iTest_suite() result(suite)
   use pFUnit_mod
   use aceNeutronDatabase_iTest
   use WrapaceNeutronDatabase_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('aceNeutronDatabase_iTest_suite')

   call suite%addTest(newTestMethod('test_aceNeutronDatabase', test_aceNeutronDatabase))


end function aceNeutronDatabase_iTest_suite

