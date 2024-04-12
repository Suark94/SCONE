module dictParser_iTest
  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : fileToDict
  use pFUnit_mod
  implicit none

contains

  !!
  !! Test Reading a Dictionary from a File
  !!
!@Test
  subroutine testFromFile()
    type(dictionary) :: dict
    integer(shortInt)  :: tempInt
    real(defReal)      :: tempReal
    character(nameLen) :: tempChar
    class(dictionary), pointer :: dictPtr
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    real(defReal), dimension(:), allocatable      :: tempRealArray
    character(nameLen), dimension(:), allocatable :: tempCharArray

    call fileToDict(dict,'./IntegrationTestFiles/testDictionary')

    ! Verify integer values
    call dict % get(tempInt, 'myInt')
    call dict % get(tempIntArray, 'intArray')

#line 30 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
  call assertEqual(7, tempInt, &
 & location=SourceLocation( &
 & 'dictParser_iTest.f90', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
#line 31 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
  call assertEqual([1, 2, 4, 5], tempIntArray, &
 & location=SourceLocation( &
 & 'dictParser_iTest.f90', &
 & 31) )
  if (anyExceptions()) return
#line 32 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"

    ! Verify real values
    call dict % get(tempReal, 'myReal')
    call dict % get(tempRealArray, 'realArray')

#line 37 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
  call assertEqual(1.3_defReal, tempReal, &
 & location=SourceLocation( &
 & 'dictParser_iTest.f90', &
 & 37) )
  if (anyExceptions()) return
#line 38 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
#line 38 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
  call assertEqual([1.0_defReal, 2.2_defReal, 3.5_defReal], tempRealArray, &
 & location=SourceLocation( &
 & 'dictParser_iTest.f90', &
 & 38) )
  if (anyExceptions()) return
#line 39 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"

    ! Verify nested dictionary
    dictPtr => dict % getDictPtr('subDict')
    call dictPtr % get(tempInt, 'myInt')
    call dictPtr % get(tempReal, 'myReal')

#line 45 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
  call assertEqual(3, tempInt, &
 & location=SourceLocation( &
 & 'dictParser_iTest.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
#line 46 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"
  call assertEqual(3.2_defReal, tempReal, &
 & location=SourceLocation( &
 & 'dictParser_iTest.f90', &
 & 46) )
  if (anyExceptions()) return
#line 47 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_iTest.f90"

  end subroutine testFromFile


end module dictParser_iTest

module WrapdictParser_iTest
   use pFUnit_mod
   use dictParser_iTest
   implicit none
   private

contains


end module WrapdictParser_iTest

function dictParser_iTest_suite() result(suite)
   use pFUnit_mod
   use dictParser_iTest
   use WrapdictParser_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('dictParser_iTest_suite')

   call suite%addTest(newTestMethod('testFromFile', testFromFile))


end function dictParser_iTest_suite

