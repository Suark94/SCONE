module aceCard_iTest

  use numPrecision
  use aceCard_class, only : aceCard
  use pFUnit_mod

  implicit none


contains


  !!
  !! Test Reloading aceCard
  !!
  !! Test Correcness of fission-related data after reloading
  !!
!@Test
  subroutine testAceReloading_fissionData()
    type(aceCard) :: ACE
    character(*),parameter :: path1 = './IntegrationTestFiles/91231JEF311.ace'
    character(*),parameter :: path2 = './IntegrationTestFiles/91232JEF311.ace'
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! Load 1st Card with Delayed Fission Neutrons
    call ACE % readFromFile(path1, 1)

    ! Test HEADER contents
    ! Mass
#line 30 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual(229.05_defReal,     ACE % AW, TOL* ACE % AW , &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
    ! Temperarture
#line 32 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual(2.5852E-08_defReal, ACE % TZ, TOL * ACE % TZ, &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 32) )
  if (anyExceptions()) return
#line 33 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 33 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual('91231.03c', ACE % ZAID, &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 33) )
  if (anyExceptions()) return
#line 34 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"

    ! Test Fission Data
#line 36 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual(8, ACE % precursorGroups(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 37 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % isFissile(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 38 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % hasNuPrompt(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 38) )
  if (anyExceptions()) return
#line 39 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 39 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % hasNuDelayed(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 39) )
  if (anyExceptions()) return
#line 40 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 40 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % hasNuTotal(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"

    ! Load 2nd Card without Delayed Fission Neutrons
    call ACE % readFromFile(path2, 1)

    ! Test HEADER Contents
    ! Mass
#line 47 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual(230.045_defReal,    ACE % AW, TOL* ACE % AW , &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
    ! Temperature
#line 49 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual(2.5852E-08_defReal, ACE % TZ, TOL * ACE % TZ, &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 49) )
  if (anyExceptions()) return
#line 50 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 50 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual('91232.03c', ACE % ZAID, &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"

    ! Test Fission Data
#line 53 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertEqual(0, ACE % precursorGroups(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 54 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % isFissile(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 54) )
  if (anyExceptions()) return
#line 55 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 55 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % hasNuPrompt(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 55) )
  if (anyExceptions()) return
#line 56 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 56 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertFalse(ACE % hasNuDelayed(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 56) )
  if (anyExceptions()) return
#line 57 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
#line 57 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"
  call assertTrue(ACE % hasNuTotal(), &
 & location=SourceLocation( &
 & 'aceCard_iTest.f90', &
 & 57) )
  if (anyExceptions()) return
#line 58 "/home/mhgk4/SCONE/NuclearData/DataDecks/Tests/aceCard_iTest.f90"

  end subroutine testAceReloading_fissionData



end module aceCard_iTest

module WrapaceCard_iTest
   use pFUnit_mod
   use aceCard_iTest
   implicit none
   private

contains


end module WrapaceCard_iTest

function aceCard_iTest_suite() result(suite)
   use pFUnit_mod
   use aceCard_iTest
   use WrapaceCard_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('aceCard_iTest_suite')

   call suite%addTest(newTestMethod('testAceReloading_fissionData', testAceReloading_fissionData))


end function aceCard_iTest_suite

