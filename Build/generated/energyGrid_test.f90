module energyGrid_test
  use numPrecision
  use universalVariables
  use pfUnit_mod
  use energyGrid_class, only : energyGrid

  implicit none

!@TestCase
  type, extends(TestCase) :: test_energyGrid
    type(energyGrid) :: linGrid
    type(energyGrid) :: logGrid
    type(energyGrid) :: unstructGrid
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_energyGrid

  !! Patameters
  real(defReal), parameter :: FP_TOL = 50*epsilon(1.0_defReal)

contains

  !!
  !! Initialise test case for structured linear energyGrid
  !!
  subroutine setUp(this)
    class(test_energyGrid), intent(inout) :: this
    character(nameLen)                    :: type

    ! Simple linear energy grid
    type = 'lin'
    call this % linGrid % init(ZERO, 20.0_defReal, 20, type)

    ! Simple logarithmic (equilethargic) grid
    type = 'log'
    call this % logGrid % init(1.0E-9_defReal, 1.0_defReal, 9, type)

    ! Simple unstructured grid CASMO-8
    call this % unstructGrid % init([10.00000_defReal, &
                                     0.821000_defReal, &
                                     0.005530_defReal, &
                                     4.00e-06_defReal, &
                                     6.25e-07_defReal, &
                                     2.80e-07_defReal, &
                                     1.40e-07_defReal, &
                                     5.80e-08_defReal, &
                                     1.00e-11_defReal])

  end subroutine setUp

  !!
  !! Deconstruct test case for structured linear energyGrid
  !!
  subroutine tearDown(this)
    class(test_energyGrid), intent(inout) :: this

    call this % linGrid % kill()
    call this % logGrid % kill()
    call this % unstructGrid % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Proper Tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test searching linear grid
  !!
!@Test
  subroutine testSearchLinEnergyGrid(this)
    class(test_energyGrid), intent(inout) :: this

    ! Valid searches
#line 76 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1, this % linGrid % search(19.5_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 77 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(20, this % linGrid % search(0.5_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 78 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(11, this % linGrid % search(9.9_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

    ! Invalid searches
#line 81 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % linGrid % search(20.5_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 82 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % linGrid % search(-0.1_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 83 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % linGrid % search(-2.1_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testSearchLinEnergyGrid

  !!
  !! Test searching logarithmic grid
  !!
!@Test
  subroutine testSearchLogEnergyGrid(this)
    class(test_energyGrid), intent(inout) :: this

    ! Valid searches
#line 95 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1, this % logGrid % search(0.7_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 96 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(6, this % logGrid % search(1.1E-6_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 96) )
  if (anyExceptions()) return
#line 97 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 97 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(9, this % logGrid % search(2.0E-9_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 97) )
  if (anyExceptions()) return
#line 98 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

    ! Invalid searches
#line 100 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % logGrid % search(1.1_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 100) )
  if (anyExceptions()) return
#line 101 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 101 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % logGrid % search(-0.1_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 101) )
  if (anyExceptions()) return
#line 102 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 102 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % logGrid % search(1.0E-11_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 102) )
  if (anyExceptions()) return
#line 103 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testSearchLogEnergyGrid

  !!
  !! Test searching unstructured grid
  !!
!@Test
  subroutine testSearchUnstructEnergyGrid(this)
    class(test_energyGrid), intent(inout) :: this

    ! Valid searches
#line 114 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1, this % unstructGrid % search(1.7_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 115 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(3, this % unstructGrid % search(5.1E-5_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 116 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(8, this % unstructGrid % search(2.0E-9_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 116) )
  if (anyExceptions()) return
#line 117 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

    ! Invalid searches
#line 119 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % unstructGrid % search(11.0_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 120 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % unstructGrid % search(-0.1_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 121 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(valueOutsideArray, this % unstructGrid % search(0.9E-11_defReal), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testSearchUnstructEnergyGrid

  !!
  !! Test getting bin values from linear grid
  !!
!@Test
  subroutine testBinLin(this)
    class(test_energyGrid), intent(inout) :: this
    real(defReal),parameter :: TOL =1.0E-9

    ! Valid Bins
#line 134 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(20.0_defReal, this % linGrid % bin(1), TOL * 20.0_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 135 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(0.0_defReal, this % linGrid % bin(21), TOL, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 135) )
  if (anyExceptions()) return
#line 136 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 136 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(17.0_defReal, this % linGrid % bin(4), TOL * 17.0_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 136) )
  if (anyExceptions()) return
#line 137 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

    ! Invalid Bins
#line 139 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(-huge(ONE), this % linGrid % bin(22), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 140 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(-huge(ONE), this % linGrid % bin(0), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testBinLin

  !!
  !! Test getting bin values from logarithmic grid
  !!
!@Test
  subroutine testBinLog(this)
    class(test_energyGrid), intent(inout) :: this
    real(defReal),parameter :: TOL =1.0E-9

    ! Valid bins
#line 153 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1.0_defReal, this % logGrid % bin(1), TOL * 1.0_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 154 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1.0E-8_defReal, this % logGrid % bin(9), TOL * 1.0E-8_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 155 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1.0E-9_defReal, this % logGrid % bin(10), TOL * 1.0E-9_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 156 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1.0E-3_defReal, this % logGrid % bin(4), TOL * 1.0E-3_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

    ! Invalid Bins
#line 159 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(-huge(ONE), this % logGrid % bin(11), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 160 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(-huge(ONE), this % logGrid % bin(0), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testBinLog

  !!
  !! Test getting bin values from logarithmic grid
  !!
!@Test
  subroutine testBinUnstruct(this)
    class(test_energyGrid), intent(inout) :: this
    real(defReal),parameter :: TOL =1.0E-9

    ! Valid bins
#line 173 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(10.00000_defReal, this % unstructGrid % bin(1), TOL * 10.00000_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 173) )
  if (anyExceptions()) return
#line 174 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 174 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(5.80e-08_defReal, this % unstructGrid % bin(8), TOL * 5.80e-08_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 175 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(1.00e-11_defReal, this % unstructGrid % bin(9), TOL * 1.00e-11_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 176 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(4.00e-06_defReal, this % unstructGrid % bin(4), TOL * 4.00e-06_defReal, &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 176) )
  if (anyExceptions()) return
#line 177 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

    ! Invalid Bins
#line 179 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(-huge(ONE), this % unstructGrid % bin(10), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 179) )
  if (anyExceptions()) return
#line 180 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 180 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(-huge(ONE), this % unstructGrid % bin(0), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testBinUnstruct

  !!
  !! Test getting size (number of bins/groups)
  !!
!@Test
  subroutine testGetSize(this)
    class(test_energyGrid), intent(inout) :: this

#line 191 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(20, this % linGrid % getSize(), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 192 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(9, this % logGrid % getSize(), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
#line 193 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"
  call assertEqual(8, this % unstructGrid % getSize(), &
 & location=SourceLocation( &
 & 'energyGrid_test.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/SharedModules/Tests/energyGrid_test.f90"

  end subroutine testGetSize

end module energyGrid_test

module WrapenergyGrid_test
   use pFUnit_mod
   use energyGrid_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_energyGrid) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use energyGrid_test
        class (test_energyGrid), intent(inout) :: this
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

end module WrapenergyGrid_test

function energyGrid_test_suite() result(suite)
   use pFUnit_mod
   use energyGrid_test
   use WrapenergyGrid_test
   type (TestSuite) :: suite

   suite = newTestSuite('energyGrid_test_suite')

   call suite%addTest(makeCustomTest('testSearchLinEnergyGrid', testSearchLinEnergyGrid))

   call suite%addTest(makeCustomTest('testSearchLogEnergyGrid', testSearchLogEnergyGrid))

   call suite%addTest(makeCustomTest('testSearchUnstructEnergyGrid', testSearchUnstructEnergyGrid))

   call suite%addTest(makeCustomTest('testBinLin', testBinLin))

   call suite%addTest(makeCustomTest('testBinLog', testBinLog))

   call suite%addTest(makeCustomTest('testBinUnstruct', testBinUnstruct))

   call suite%addTest(makeCustomTest('testGetSize', testGetSize))


end function energyGrid_test_suite

