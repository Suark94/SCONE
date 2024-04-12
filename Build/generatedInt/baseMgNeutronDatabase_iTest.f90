module baseMgNeutronDatabase_iTest

  use numPrecision
  use endfConstants
  use pFUnit_mod
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use particle_class,     only : particle

  ! Nuclear Data Objects & Interfaces
  use baseMgNeutronDatabase_class, only : baseMgNeutronDatabase, baseMgNeutronDatabase_CptrCast, &
                                          baseMgNeutronDatabase_TptrCast
  use baseMgNeutronMaterial_class, only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast, &
                                          baseMgNeutronMaterial_TptrCast
  use fissionMG_class,             only : fissionMG, fissionMG_TptrCast
  use multiScatterMG_class,        only : multiScatterMG, multiScatterMG_CptrCast, &
                                          multiScatterMG_TptrCast
  use multiScatterP1MG_class,      only : multiScatterP1MG, multiScatterP1MG_TptrCast
  use materialMenu_mod,            only : mm_init => init, mm_kill => kill
  use nuclearDatabase_inter,       only : nuclearDatabase
  use materialHandle_inter,        only : materialHandle
  use nuclideHandle_inter,         only : nuclideHandle
  use neutronXsPackages_class,     only : neutronMacroXSs
  use reactionHandle_inter,        only : reactionHandle



  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = "   &
  mat1 { temp 273;                              &
         composition {                          &
         1001.03 5.028E-02;                     &
         8016.03 2.505E-02;                     &
         }                                      &
         xsFile ./IntegrationTestFiles/mgMat1;  &
       }                                        &
  mat2  { temp 1;                               &
          composition {                         &
          92233.03 2.286E-02;                   &
          8016.03  4.572E-02;                   &
          }                                     &
          xsFile ./IntegrationTestFiles/mgMat2; &
        }"


contains

  !!
  !! Monster test to build and verify data in baseMgNeutronDatabase with P0 scattering
  !!
!@Test
  subroutine testBaseMgNeutronDatabaseWithP0()
    type(baseMgNeutronDatabase), target  :: database
    class(nuclearDatabase), pointer      :: data_ptr
    type(dictionary)                     :: databaseDef
    type(dictionary)                     :: matMenuDict
    type(particle)                       :: p
    type(neutronMacroXSs)                :: xss
    type(baseMgNeutronMaterial),pointer  :: mat
    class(baseMgNeutronMaterial),pointer :: matClass
    class(reactionHandle), pointer       :: reac
    real(defReal),parameter :: TOL = 1.0E-6_defReal


    data_ptr => database

    ! Load materialMenu
    call charToDict(matMenuDict, MAT_INPUT_STR)
    call mm_init(matMenuDict )

    ! Build database
    call databaseDef % init(1)
    call databaseDef % store('PN','P0')
    call database % init(databaseDef, data_ptr, silent = .true.)
    call database % activate([1])

    ! Varify number of groups
#line 80 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(4, database % nGroups(), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting Transport XS
    p % G = 1
#line 84 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(2.1_defReal, database % getTransMatXS(p, 1), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting Total XS
    p % G = 1
#line 88 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(3.1_defReal, database % getTotalMatXS(p, 2), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    p % G = 3
#line 91 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(6.0_defReal, database % getTotalMatXS(p, 1), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 91) )
  if (anyExceptions()) return
#line 92 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting Majorant
    p % G = 1
#line 95 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(2.1_defReal, database % getMajorantXS(p), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"


    ! Get a material and verify macroXSS
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(2))
#line 100 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(mat), "Type Ptr Cast has failed", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 100) )
  if (anyExceptions()) return
#line 101 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
    call mat % getMacroXSs(xss, 1, p % pRNG)

    ! Check that is fissile
#line 104 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(mat % isFissile(), "Is not fissile but should", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 104) )
  if (anyExceptions()) return
#line 105 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

#line 106 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(3.1_defReal, xss % total, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 107 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, xss % elasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 108 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(1.1_defReal, xss % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 109 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(1.0_defReal, xss % capture, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 110 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(1.0_defReal, xss % fission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 111 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(2.3_defReal, xss % nuFission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    matClass => baseMgNeutronMaterial_CptrCast(database % getMaterial(1))
#line 114 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(matClass), "Type Ptr Cast has failed", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
    call matClass % getMacroXSs(xss, 4, p % pRNG)

#line 117 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(matClass % isFissile(), "Is fissile but should not", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

#line 119 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(7.1_defReal, xss % total, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 120 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, xss % elasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 121 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(3.1_defReal, xss % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 122 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(4.0_defReal, xss % capture, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 123 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(0.0_defReal, xss % fission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 124 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(0.0_defReal, xss % nuFission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 124) )
  if (anyExceptions()) return
#line 125 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get some invalid Materials
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(0))
#line 128 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(mat), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(-2))
#line 131 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(mat), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(3))
#line 134 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(mat), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get Fission reaction and verify type
    reac => fissionMG_TptrCast(database % getReaction(macroFission, 1))
#line 138 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), "Pointer for the mission reaction is not null", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => fissionMG_TptrCast(database % getReaction(macroFission, 2))
#line 141 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(reac), "Pointer fission reaction is wrong type or null", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get Scattering reaction and verify type
    reac => multiScatterMG_TptrCast(database % getReaction(macroIEScatter, 1))
#line 145 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(reac), "Wrong type of scattering reaction", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get some invalid reactions
    reac => database % getReaction(anyScatter, 0)
#line 149 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => database % getReaction(anyScatter, -1)
#line 152 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => database % getReaction(anyScatter, 3)
#line 155 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => database % getReaction(anyCapture, 1)
#line 158 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! **** Note that anyFission is not present !
    reac => database % getReaction(anyFission, 2)
#line 162 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 162) )
  if (anyExceptions()) return
#line 163 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting nuclide
#line 165 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(database % getNuclide(1)), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Clean up
    call database % kill()
    call mm_kill()
    call matMenuDict % kill()
    call databaseDef % kill()

  end subroutine testBaseMgNeutronDatabaseWithP0

  !!
  !! Monster test to build and verify data in baseMgNeutronDatabase with P1 scattering
  !! *Copy and pasted from the above with only the type of scattering changed
  !!
!@Test
  subroutine testBaseMgNeutronDatabaseWithP1()
    type(baseMgNeutronDatabase), target  :: database
    class(nuclearDatabase), pointer      :: data_ptr
    type(dictionary)                     :: databaseDef
    type(dictionary)                     :: matMenuDict
    type(particle)                       :: p
    type(neutronMacroXSs)                :: xss
    type(baseMgNeutronMaterial),pointer  :: mat
    class(baseMgNeutronMaterial),pointer :: matClass
    class(reactionHandle), pointer       :: reac
    real(defReal),parameter :: TOL = 1.0E-6_defReal


    data_ptr => database

    ! Load materialMenu
    call charToDict(matMenuDict, MAT_INPUT_STR)
    call mm_init(matMenuDict )

    ! Build database
    call databaseDef % init(1)
    call databaseDef % store('PN','P1')
    call database % init(databaseDef, data_ptr, silent = .true.)
    call database % activate([1])

    ! Varify number of groups
#line 206 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(4, database % nGroups(), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 206) )
  if (anyExceptions()) return
#line 207 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting Transport XS
    p % G = 1
#line 210 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(2.1_defReal, database % getTransMatXS(p, 1), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 210) )
  if (anyExceptions()) return
#line 211 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting Total XS
    p % G = 1
#line 214 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(3.1_defReal, database % getTotalMatXS(p, 2), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 214) )
  if (anyExceptions()) return
#line 215 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    p % G = 3
#line 217 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(6.0_defReal, database % getTotalMatXS(p, 1), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 217) )
  if (anyExceptions()) return
#line 218 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting Majorant
    p % G = 1
#line 221 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(2.1_defReal, database % getMajorantXS(p), TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"


    ! Get a material and verify macroXSS
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(2))
#line 226 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(mat), "Type Ptr Cast has failed", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
    call mat % getMacroXSs(xss, 1, p % pRNG)

    ! Check that is fissile
#line 230 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(mat % isFissile(), "Is not fissile but should", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

#line 232 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(3.1_defReal, xss % total, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 232) )
  if (anyExceptions()) return
#line 233 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 233 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, xss % elasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 233) )
  if (anyExceptions()) return
#line 234 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 234 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(1.1_defReal, xss % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 234) )
  if (anyExceptions()) return
#line 235 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 235 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(1.0_defReal, xss % capture, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 236 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(1.0_defReal, xss % fission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 236) )
  if (anyExceptions()) return
#line 237 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 237 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(2.3_defReal, xss % nuFission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    matClass => baseMgNeutronMaterial_CptrCast(database % getMaterial(1))
#line 240 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(matClass), "Type Ptr Cast has failed", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
    call matClass % getMacroXSs(xss, 4, p % pRNG)

#line 243 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(matClass % isFissile(), "Is fissile but should not", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 243) )
  if (anyExceptions()) return
#line 244 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

#line 245 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(7.1_defReal, xss % total, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 246 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(ZERO, xss % elasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 246) )
  if (anyExceptions()) return
#line 247 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 247 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(3.1_defReal, xss % inelasticScatter, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 247) )
  if (anyExceptions()) return
#line 248 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 248 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(4.0_defReal, xss % capture, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 248) )
  if (anyExceptions()) return
#line 249 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 249 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(0.0_defReal, xss % fission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 249) )
  if (anyExceptions()) return
#line 250 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
#line 250 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertEqual(0.0_defReal, xss % nuFission, TOL, &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 250) )
  if (anyExceptions()) return
#line 251 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get some invalid Materials
    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(0))
#line 254 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(mat), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 254) )
  if (anyExceptions()) return
#line 255 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(-2))
#line 257 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(mat), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 257) )
  if (anyExceptions()) return
#line 258 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    mat => baseMgNeutronMaterial_TptrCast(database % getMaterial(3))
#line 260 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(mat), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 260) )
  if (anyExceptions()) return
#line 261 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get Fission reaction and verify type
    reac => fissionMG_TptrCast(database % getReaction(macroFission, 1))
#line 264 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), "Pointer for the mission reaction is not null", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 264) )
  if (anyExceptions()) return
#line 265 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => fissionMG_TptrCast(database % getReaction(macroFission, 2))
#line 267 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(reac), "Pointer fission reaction is wrong type or null", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 267) )
  if (anyExceptions()) return
#line 268 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get Scattering reaction and verify type
    reac => multiScatterP1MG_TptrCast(database % getReaction(macroIEScatter, 1))
#line 271 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertTrue(associated(reac), "Wrong type of scattering reaction", &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 271) )
  if (anyExceptions()) return
#line 272 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Get some invalid reactions
    reac => database % getReaction(anyScatter, 0)
#line 275 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 275) )
  if (anyExceptions()) return
#line 276 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => database % getReaction(anyScatter, -1)
#line 278 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 278) )
  if (anyExceptions()) return
#line 279 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => database % getReaction(anyScatter, 3)
#line 281 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 281) )
  if (anyExceptions()) return
#line 282 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    reac => database % getReaction(anyCapture, 1)
#line 284 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 284) )
  if (anyExceptions()) return
#line 285 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! **** Note that anyFission is not present !
    reac => database % getReaction(anyFission, 2)
#line 288 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(reac), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 288) )
  if (anyExceptions()) return
#line 289 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Test getting nuclide
#line 291 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"
  call assertFalse(associated(database % getNuclide(1)), &
 & location=SourceLocation( &
 & 'baseMgNeutronDatabase_iTest.f90', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/home/mhgk4/SCONE/NuclearData/mgNeutronData/baseMgNeutron/Tests/baseMgNeutronDatabase_iTest.f90"

    ! Clean up
    call database % kill()
    call mm_kill()
    call matMenuDict % kill()
    call databaseDef % kill()

  end subroutine testBaseMgNeutronDatabaseWithP1



end module baseMgNeutronDatabase_iTest

module WrapbaseMgNeutronDatabase_iTest
   use pFUnit_mod
   use baseMgNeutronDatabase_iTest
   implicit none
   private

contains


end module WrapbaseMgNeutronDatabase_iTest

function baseMgNeutronDatabase_iTest_suite() result(suite)
   use pFUnit_mod
   use baseMgNeutronDatabase_iTest
   use WrapbaseMgNeutronDatabase_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('baseMgNeutronDatabase_iTest_suite')

   call suite%addTest(newTestMethod('testBaseMgNeutronDatabaseWithP0', testBaseMgNeutronDatabaseWithP0))

   call suite%addTest(newTestMethod('testBaseMgNeutronDatabaseWithP1', testBaseMgNeutronDatabaseWithP1))


end function baseMgNeutronDatabase_iTest_suite

