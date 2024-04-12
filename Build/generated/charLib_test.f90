module charLib_test

  use numPrecision
  use charLib_func, only : splitChar
  use pFUnit_mod

  implicit none


contains

  !!
  !! Test char splitting
  !!
!@Test
  subroutine testSplitChar()
    character(:), allocatable :: line
    integer(shortInt), dimension(:,:),allocatable :: SE

    ! CASE 1: Standard Use Case
    line = " This is a first attempt at splitting"

    SE = splitChar(line, ' ')

#line 25 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([2, 7, 10, 12, 18, 26, 29], SE(1,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 25) )
  if (anyExceptions()) return
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
#line 26 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([5, 8, 10, 16, 24, 27, 37], SE(2,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 26) )
  if (anyExceptions()) return
#line 27 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"

    deallocate(line)
    deallocate(SE)

    ! CASE 2: Long delimiter sequences
    line ="...........12345......123..."
    SE = splitChar(line,'.')

#line 35 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([12, 23], SE(1,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 35) )
  if (anyExceptions()) return
#line 36 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
#line 36 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([16, 25], SE(2,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 36) )
  if (anyExceptions()) return
#line 37 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
    deallocate(line)
    deallocate(SE)

    ! CASE 3: Empty String
    line = ""
    SE = splitChar(line,'.')
#line 43 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([0], SE(1,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 43) )
  if (anyExceptions()) return
#line 44 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
#line 44 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([0], SE(2,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 44) )
  if (anyExceptions()) return
#line 45 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
    deallocate(line)
    deallocate(SE)

    ! CASE 4: Delimiter only string
    line ="#################################"
    SE = splitChar(line,'#')
#line 51 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([0], SE(1,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 51) )
  if (anyExceptions()) return
#line 52 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
#line 52 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
  call assertEqual([0], SE(2,:), &
 & location=SourceLocation( &
 & 'charLib_test.f90', &
 & 52) )
  if (anyExceptions()) return
#line 53 "/home/mhgk4/SCONE/SharedModules/Tests/charLib_test.f90"
    deallocate(line)
    deallocate(SE)


  end subroutine testSplitChar


    
end module charLib_test

module WrapcharLib_test
   use pFUnit_mod
   use charLib_test
   implicit none
   private

contains


end module WrapcharLib_test

function charLib_test_suite() result(suite)
   use pFUnit_mod
   use charLib_test
   use WrapcharLib_test
   type (TestSuite) :: suite

   suite = newTestSuite('charLib_test_suite')

   call suite%addTest(newTestMethod('testSplitChar', testSplitChar))


end function charLib_test_suite

