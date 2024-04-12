module timer_test
  use numPrecision
  use timer_mod, only : registerTimer, timerStart, timerStop, timerReset, timerTime
  use pfUnit_mod

  implicit none


contains

  !!
  !! Test logic for dynamic space for timers
  !! Make sure it does not segment
  !!
!@Test
  subroutine testRegisterTimer()
    integer(shortInt) :: i, j

    do i = 1,100
      j = registerTimer('myName')
    end do

#line 23 "/home/mhgk4/SCONE/SharedModules/Tests/timer_test.f90"
  call assertEqual(100, j, &
 & location=SourceLocation( &
 & 'timer_test.f90', &
 & 23) )
  if (anyExceptions()) return
#line 24 "/home/mhgk4/SCONE/SharedModules/Tests/timer_test.f90"

  end subroutine testRegisterTimer

end module timer_test

module Wraptimer_test
   use pFUnit_mod
   use timer_test
   implicit none
   private

contains


end module Wraptimer_test

function timer_test_suite() result(suite)
   use pFUnit_mod
   use timer_test
   use Wraptimer_test
   type (TestSuite) :: suite

   suite = newTestSuite('timer_test_suite')

   call suite%addTest(newTestMethod('testRegisterTimer', testRegisterTimer))


end function timer_test_suite

