module materialMenu_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict

  use materialMenu_mod,   only : init_menu => init, kill_menu => kill, nameMap, materialDefs, &
                                              display, nMat, getMatPtr, materialItem
  use pFUnit_mod

  implicit none

  character(*),parameter :: INPUT_STR = "     &
  mat1 { temp 273;                            &
         composition {                        &
         1001.03 12;                          &
         8016.07 0.00654;                     &
         }                                    &
         xsPath ./A_PATH;                     &
       }                                      &
  mat2 { temp 1;                              &
          composition {                       &
          100253.00 7.0E+4;                   &
          }                                   &
        }                                     "

contains

  !!
  !! Test the initialisation of matarialMenu
  !!
!@Test
  subroutine testMaterialMenu()
    type(dictionary)           :: matDict
    type(dictionary)           :: emptyDict
    type(materialItem),pointer :: matPtr
    integer(shortInt)          :: i1, i2, i
    character(nameLen)         :: name
    real(defReal), parameter   :: TOL = 1.0E-6_defReal

    ! Build from empty and see if anything crashes
    call emptyDict % init(1)
    call init_menu(emptyDict)
    call kill_menu()

    ! Build some real definitions
    call charToDict(matDict, INPUT_STR)
    call init_menu(matDict)

    ! Check that materials are present in the nameMap
    name = 'mat1'
    i1 = nameMap % getOrDefault(name,0)
#line 53 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertTrue(i1 == 1 .or. i1 == 2, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    name = 'mat2'
    i2 = nameMap % getOrDefault(name,0)
#line 57 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertTrue(i2 == 1 .or. i2 == 2, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    ! Verify definitions

    ! mat 1
#line 62 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(273.0_defReal, materialDefs(i1) % T, TOL, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 63 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(i1, materialDefs(i1) % matIdx, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 63) )
  if (anyExceptions()) return
#line 64 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    call materialDefs(i1) % extraInfo % get(name,'xsPath')
#line 66 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual('./A_PATH', trim(name), &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    ! Check individual compositions
    do i=1,2
      if(materialDefs(i1) % nuclides(i) % Z == 1) then
#line 71 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(1, materialDefs(i1) % nuclides(i) % A, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 72 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(3, materialDefs(i1) % nuclides(i) % T, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 72) )
  if (anyExceptions()) return
#line 73 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 73 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(12.0_defReal, materialDefs(i1) % dens(i), TOL*12.0_defReal, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 73) )
  if (anyExceptions()) return
#line 74 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
      else if(materialDefs(i1) % nuclides(i) % Z == 8) then
#line 75 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(16, materialDefs(i1) % nuclides(i) % A, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 75) )
  if (anyExceptions()) return
#line 76 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 76 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(7, materialDefs(i1) % nuclides(i) % T, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 77 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(0.00654_defReal, materialDefs(i1) % dens(i), TOL*0.00654_defReal, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
      else
#line 79 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertTrue(.false., 'Error when reading Mat 1 compositions', &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 79) )
  if (anyExceptions()) return
#line 80 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
      end if
    end do

    ! mat 2
#line 84 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(1.0_defReal, materialDefs(i2) % T, TOL, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 85 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(i2, materialDefs(i2) % matIdx, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 86 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(7.0E+4_defReal, materialDefs(i2) % dens(1), TOL*7.0E+4_defReal, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 87 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(100, materialDefs(i2) % nuclides(1) % Z, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 88 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(253, materialDefs(i2) % nuclides(1) % A, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
#line 89 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(0,   materialDefs(i2) % nuclides(1) % T, &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    ! Get number of materials
#line 92 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual(2, nMat(), &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 92) )
  if (anyExceptions()) return
#line 93 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    ! Get pointer to material1
    matPtr => getMatPtr(1)
#line 96 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual('mat1', trim(matPtr % name), &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    ! Test writing nuclide info back to definition string
#line 99 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual('1001.03', matPtr % nuclides(1) % toChar(), &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    matPtr => getMatPtr(2)
#line 102 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"
  call assertEqual('100253.00', matPtr % nuclides(1) % toChar(), &
 & location=SourceLocation( &
 & 'materialMenu_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/NuclearData/Tests/materialMenu_test.f90"

    ! Clean
    call kill_menu()

  end subroutine testMaterialMenu


end module materialMenu_test

module WrapmaterialMenu_test
   use pFUnit_mod
   use materialMenu_test
   implicit none
   private

contains


end module WrapmaterialMenu_test

function materialMenu_test_suite() result(suite)
   use pFUnit_mod
   use materialMenu_test
   use WrapmaterialMenu_test
   type (TestSuite) :: suite

   suite = newTestSuite('materialMenu_test_suite')

   call suite%addTest(newTestMethod('testMaterialMenu', testMaterialMenu))


end function materialMenu_test_suite

