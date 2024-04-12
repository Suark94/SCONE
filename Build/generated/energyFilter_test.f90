module energyFilter_test

  use numPrecision
  use energyFilter_class, only : energyFilter
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_energyFilter
    private
    type(energyFilter) :: filter
    type(energyFilter) :: mgFilter
    type(energyFilter) :: cemgFilter
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyFilter

  !! Parameters
  real(defReal), parameter    :: E_MIN = 1.27E-6_defReal
  real(defReal), parameter    :: E_MAX = 1.34_defReal
  real(defReal), parameter    :: E_Delta = 0.999_defReal ! Rel. diff. Must be < 1 for correct tests !
  integer(shortInt),parameter :: G_TOP = 3
  integer(shortInt),parameter :: G_LOW = 7

contains

  !!
  !! Sets up test_energyFilter object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_energyFilter), intent(inout) :: this
    type(dictionary)                        :: tempDict

    ! Build CE filter
    call tempDict % init(4)
    call tempDict % store('Emin', E_MIN)
    call tempDict % store('Emax', E_MAX)
    call this % filter % init(tempDict)

    ! Build MG-CE filter
    call tempDict % store('Gtop', G_TOP)
    call tempDict % store('Glow', G_LOW)
    call this % cemgFilter % init(tempDict)

    ! Build MG filter
    call tempDict % kill()
    call tempDict % init(2)
    call tempDict % store('Gtop', G_TOP)
    call tempDict % store('Glow', G_LOW)
    call this % mgFilter % init(tempDict)

  end subroutine setUp

  !!
  !! Kills test_energyFilter object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_energyFilter), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the CE filter
  !!
!@Test
  subroutine testCEFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    real(defReal)                           :: testE
    logical(defBool)                        :: filterRes

    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    state % E = testE
    filterRes = this % filter % isPass(state)
#line 84 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Emin -> TRUE
    testE = E_MIN
    state % E = testE
    filterRes = this % filter % isPass(state)
#line 90 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    state % E = testE
    filterRes = this % filter % isPass(state)
#line 96 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! E max -> TRUE
    testE = E_MAX
    state % E = testE
    filterRes = this % filter % isPass(state)
#line 102 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    state % E = testE
    filterRes = this % filter % isPass(state)
#line 108 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 108) )
  if (anyExceptions()) return
#line 109 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

  end subroutine testCEFilter

  !!
  !! Test correct behaviour of the CEMG filter
  !!
!@Test
  subroutine testMGFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    integer(shortInt)                       :: testG
    logical(defBool)                        :: filterRes

    state % isMG = .true.
    ! Below G_TOP -> FALE
    testG = G_TOP - 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 127 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Exactly G_TOP -> TRUE
    testG = G_TOP
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 133 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 133) )
  if (anyExceptions()) return
#line 134 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Above G_TOP -> TRUE
    testG = G_TOP + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 139 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Above G_LOW -> FALSE
    testG = G_LOW + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 145 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Exactly G_LOW -> TRUE
    testG = G_LOW
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 151 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTRUE(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

  end subroutine testMGFilter

  !!
  !! Test correct behaviour of the CEMG filter
  !!
!@Test
  subroutine testCEMGFilter(this)
    class(test_energyFilter), intent(inout) :: this
    type(particleState)                     :: state
    real(defReal)                           :: testE
    integer(shortInt)                       :: testG
    logical(defBool)                        :: filterRes

    state % isMG = .false.
    ! Below Emin -> FALSE
    testE = E_MIN * E_Delta
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
#line 171 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Emin -> TRUE
    testE = E_MIN
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
#line 177 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Below E max -> TRUE
    testE = E_MAX * E_Delta
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
#line 183 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! E max -> TRUE
    testE = E_MAX
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
#line 189 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Above Emax -> FALSE
    testE = E_MAX / E_Delta
    state % E = testE
    filterRes = this % cemgFilter % isPass(state)
#line 195 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 195) )
  if (anyExceptions()) return
#line 196 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    state % isMG = .true.
    ! Below G_TOP -> FALE
    testG = G_TOP - 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 202 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Exactly G_TOP -> TRUE
    testG = G_TOP
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 208 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 208) )
  if (anyExceptions()) return
#line 209 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Above G_TOP -> TRUE
    testG = G_TOP + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 214 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTrue(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 214) )
  if (anyExceptions()) return
#line 215 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Above G_LOW -> FALSE
    testG = G_LOW + 1
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 220 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertFalse(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 220) )
  if (anyExceptions()) return
#line 221 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

    ! Exactly G_LOW -> TRUE
    testG = G_LOW
    state % G = testG
    filterRes = this % cemgFilter % isPass(state)
#line 226 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"
  call assertTRUE(filterRes, &
 & location=SourceLocation( &
 & 'energyFilter_test.f90', &
 & 226) )
  if (anyExceptions()) return
#line 227 "/home/mhgk4/SCONE/Tallies/TallyFilters/Tests/energyFilter_test.f90"

  end subroutine testCEMGFilter



end module energyFilter_test

module WrapenergyFilter_test
   use pFUnit_mod
   use energyFilter_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_energyFilter) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use energyFilter_test
        class (test_energyFilter), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod) result(aTest)
#ifdef INTEL_13
      use pfunit_mod, only: testCase
#endif
      type (WrapUserTestCase) :: aTest
#ifdef INTEL_13
      target :: aTest
      class (WrapUserTestCase), pointer :: p
#endif
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
   end function makeCustomTest

end module WrapenergyFilter_test

function energyFilter_test_suite() result(suite)
   use pFUnit_mod
   use energyFilter_test
   use WrapenergyFilter_test
   type (TestSuite) :: suite

   suite = newTestSuite('energyFilter_test_suite')

   call suite%addTest(makeCustomTest('testCEFilter', testCEFilter))

   call suite%addTest(makeCustomTest('testMGFilter', testMGFilter))

   call suite%addTest(makeCustomTest('testCEMGFilter', testCEMGFilter))


end function energyFilter_test_suite

