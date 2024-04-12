module aceLibraryRead_iTest

  use numPrecision
  use aceLibrary_mod, only : load, new_neutronACE
  use aceCard_class,  only : aceCard
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test reading ACE Library
  !!   Make sure that a libery file is read without any problems and that all ACE cards
  !!   can be read as well
  !!
!@Test
  subroutine testReadingACELibrary()
    type(aceCard)      :: ACE
    character(nameLen) :: ZAID
    real(defReal), parameter :: TOL = 1.0E-6


    ! Load library
    call load('./IntegrationTestFiles/testLib')

    ! Load U-233
    ZAID = '92233.03'
    call new_neutronACE(ACE,ZAID)

    ! Verify
#line 32 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual('92233.03c', ACE % ZAID, &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 32) )
  if (anyExceptions()) return
#line 33 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
#line 33 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual(231.038000_defReal, ACE % AW, TOL * 231.038000_defReal, &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
#line 34 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual(25, ACE % numMT(), &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 34) )
  if (anyExceptions()) return
#line 35 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"

    ! Load H-1 at 300 K
    ZAID = '1001.03'
    call new_neutronACE(ACE,ZAID)
#line 39 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual('1001.03c', ACE % ZAID, &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
#line 40 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual(0.999170_defReal, ACE % AW, TOL, &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
#line 41 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual(2.5852E-08 , ACE % TZ ,TOL * 2.5852E-08, &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
#line 42 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"
  call assertEqual(1, ACE % numMT(), &
 & location=SourceLocation( &
 & 'aceLibraryRead_iTest.f90', &
 & 42) )
  if (anyExceptions()) return
#line 43 "/home/mhgk4/SCONE/NuclearData/Tests/aceLibraryRead_iTest.f90"


  end subroutine testReadingACELibrary


end module aceLibraryRead_iTest

module WrapaceLibraryRead_iTest
   use pFUnit_mod
   use aceLibraryRead_iTest
   implicit none
   private

contains


end module WrapaceLibraryRead_iTest

function aceLibraryRead_iTest_suite() result(suite)
   use pFUnit_mod
   use aceLibraryRead_iTest
   use WrapaceLibraryRead_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('aceLibraryRead_iTest_suite')

   call suite%addTest(newTestMethod('testReadingACELibrary', testReadingACELibrary))


end function aceLibraryRead_iTest_suite

