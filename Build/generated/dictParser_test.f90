module dictParser_test
  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use pFUnit_mod
  implicit none

contains

  !!
  !! Sets up test_dictionary object we can use in a number of tests
  !!
!@Test
  subroutine testFromChar()
    type(dictionary)   :: dict
    integer(shortInt)  :: tempInt
    real(defReal)      :: tempReal
    character(nameLen) :: tempChar
    class(dictionary), pointer :: dictPtr
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    real(defReal), dimension(:), allocatable      :: tempRealArray
    character(nameLen), dimension(:), allocatable :: tempCharArray
    character(*),parameter :: tape = " myInt 7;                                 &
                                       myChar my;                               &
                                       myReal 1.3;                              &
                                       weirdFloat 1E-11;                        &
                                       intArray (1 2 4 5);                      &
                                       realArray (1.1 2.2 3.4 1E-11);           &
                                       charArray (One element );                &
                                       subDict { myInt 3; myReal 3.2; }"


    ! Create dictionary
    call charToDict(dict, tape)

    ! Verify integer values
    call dict % get(tempInt, 'myInt')
    call dict % get(tempIntArray, 'intArray')

#line 40 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual(7, tempInt, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 40) )
  if (anyExceptions()) return
#line 41 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
#line 41 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual([1, 2, 4, 5], tempIntArray, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 41) )
  if (anyExceptions()) return
#line 42 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"

    ! Verify real values
    call dict % get(tempReal, 'myReal')
    call dict % get(tempRealArray, 'realArray')

#line 47 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual(1.3_defReal, tempReal, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
#line 48 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual([1.1_defReal, 2.2_defReal, 3.4_defReal, 1.0E-11_defReal], tempRealArray, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"

    ! Verify a problematic value
    ! It may not be parsed correclty with wrong fromat settings
    call dict % get(tempReal, 'weirdFloat')
#line 53 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual(1.0E-11_defReal, tempReal, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 53) )
  if (anyExceptions()) return
#line 54 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"

    ! Verify nested dictionary
    dictPtr => dict % getDictPtr('subDict')
    call dictPtr % get(tempInt, 'myInt')
    call dictPtr % get(tempReal, 'myReal')

#line 60 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual(3, tempInt, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 60) )
  if (anyExceptions()) return
#line 61 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
#line 61 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"
  call assertEqual(3.2_defReal, tempReal, &
 & location=SourceLocation( &
 & 'dictParser_test.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/DataStructures/Tests/dictParser_test.f90"

  end subroutine testFromChar


end module dictParser_test

module WrapdictParser_test
   use pFUnit_mod
   use dictParser_test
   implicit none
   private

contains


end module WrapdictParser_test

function dictParser_test_suite() result(suite)
   use pFUnit_mod
   use dictParser_test
   use WrapdictParser_test
   type (TestSuite) :: suite

   suite = newTestSuite('dictParser_test_suite')

   call suite%addTest(newTestMethod('testFromChar', testFromChar))


end function dictParser_test_suite

