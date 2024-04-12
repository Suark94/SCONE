module grid_test
  use numPrecision
  use universalVariables
  use pfUnit_mod
  use grid_class, only : grid

  implicit none

!@TestCase
  type, extends(TestCase) :: test_grid
    type(grid) :: linGrid
    type(grid) :: logGrid
    type(grid) :: unstructGrid
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_grid

  !! Patameters
  real(defReal), parameter :: FP_TOL = 50*epsilon(1.0_defReal)

contains

  !!
  !! Initialise test case for structured linear grid
  !!
  subroutine setUp(this)
    class(test_grid), intent(inout) :: this
    character(nameLen)              :: type

    ! Typical PWR Assembly Grid - Linear Grid
    type = 'lin'
    call this % linGrid % init(-10.71_defReal, 10.71_defReal, 17, type)

    ! Equilethargic 70 Group Spectrum - Log Grid
    type = 'log'
    call this % logGrid % init(1.0E-11_defReal, 20.0_defReal, 70, type)

    ! CASMO 8 Group Structure - Unstructured grid
    call this % unstructGrid % init([1.00e-11_defReal, &
                                     5.80e-08_defReal, &
                                     1.40e-07_defReal, &
                                     2.80e-07_defReal, &
                                     6.25e-07_defReal, &
                                     4.00e-06_defReal, &
                                     0.005530_defReal, &
                                     0.821000_defReal, &
                                     10.00000_defReal])

  end subroutine setUp

  !!
  !! Deconstruct test case for structured linear grid
  !!
  subroutine tearDown(this)
    class(test_grid), intent(inout) :: this

    call this % linGrid % kill()
    call this % logGrid % kill()
    call this % unstructGrid % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on linear structured grid
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test valid quary about the bin boundary
!!
!@test
  subroutine testValidBinQuery_lin(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(5)               :: bin
    integer(shortInt),dimension(5),parameter :: idx = [1,4,7,14,18]

    bin = this % linGrid % bin(idx)
#line 78 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(-10.71_defReal, bin(1), abs(FP_TOL*bin(1)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 79 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(-6.93_defReal,  bin(2), abs(FP_TOL*bin(2)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 79) )
  if (anyExceptions()) return
#line 80 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 80 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(-3.15_defReal,  bin(3), abs(FP_TOL*bin(3)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 81 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(5.67_defReal,   bin(4), abs(FP_TOL*bin(4)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 82 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(10.71_defReal,  bin(5), abs(FP_TOL*bin(5)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testValidBinQuery_lin

!!
!! Test invalid quary about the bin boundary
!!
!@test
  subroutine testInvalidBinQuery_lin(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(3)               :: bin
    integer(shortInt),dimension(3),parameter :: idx = [-1,0, 20]
    real(defReal),parameter                  :: INVALID_BIN = -huge(bin)

    bin = this % linGrid % bin(idx)
#line 97 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(1), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 97) )
  if (anyExceptions()) return
#line 98 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 98 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(2), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 99 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(3), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testInvalidBinQuery_lin

!!
!! Test search
!!
!@test
  subroutine testSearch_lin(this)
    class(test_grid), intent(inout) :: this
    integer(shortInt),dimension(5)  :: idx
    real(defReal),dimension(5)      :: keys

    keys = [-10.71_defReal, 3.13245_defReal, -8.96_defReal, -20.0_defReal, 10.72_defReal]
    idx = this % linGrid % search(keys)

#line 115 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(1,  idx(1), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 116 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(11, idx(2), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 116) )
  if (anyExceptions()) return
#line 117 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 117 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(2,  idx(3), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

    ! Invalid Searches
#line 120 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(valueOutsideArray, idx(4), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 120) )
  if (anyExceptions()) return
#line 121 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 121 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(valueOutsideArray, idx(5), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testSearch_lin

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on logarithmic structured grid
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test valid quary about the bin boundary
!!
!@test
  subroutine testValidBinQuery_log(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(4)               :: bin
    integer(shortInt),dimension(4),parameter :: idx = [1,35,70,71]

    bin = this % logGrid % bin(idx)
#line 139 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(1.0E-11_defReal,               bin(1), abs(FP_TOL*bin(1)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 140 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(9.435957972757912E-6_defReal,  bin(2), abs(FP_TOL*bin(2)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 141 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(13.344459739056777_defReal,    bin(3), abs(FP_TOL*bin(3)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 142 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(20.0_defReal,                  bin(4), abs(FP_TOL*bin(4)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"


  end subroutine testValidBinQuery_log

!!
!! Test invalid quary about the bin boundary
!!
!@test
  subroutine testInvalidBinQuery_log(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(3)               :: bin
    integer(shortInt),dimension(3),parameter :: idx = [-1,0, 72]
    real(defReal),parameter                  :: INVALID_BIN = -huge(bin)

    bin = this % logGrid % bin(idx)
#line 158 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(1), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 159 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(2), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 160 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(3), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testInvalidBinQuery_log

!!
!! Test search
!!
!@test
  subroutine testSearch_log(this)
    class(test_grid), intent(inout) :: this
    integer(shortInt),dimension(5)  :: idx
    real(defReal),dimension(5)      :: keys

    keys = [1.0E-11_defReal, 6.7E-4_defReal, 1.0_defReal, 9.0E-12_defReal, 21.0_defReal]
    idx = this % logGrid % search(keys)

#line 176 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(1,  idx(1), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 176) )
  if (anyExceptions()) return
#line 177 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 177 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(45, idx(2), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 178 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(63, idx(3), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

    ! Invalid Searches
#line 181 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(valueOutsideArray, idx(4), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 181) )
  if (anyExceptions()) return
#line 182 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 182 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(valueOutsideArray, idx(5), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 182) )
  if (anyExceptions()) return
#line 183 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testSearch_log

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on unstructured grid
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test valid quary about the bin boundary
!!
!@test
  subroutine testValidBinQuery_unstruct(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(4)               :: bin
    integer(shortInt),dimension(4),parameter :: idx = [1,3,7,9]

    bin = this % unstructGrid % bin(idx)
#line 200 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(1.0E-11_defReal,  bin(1), abs(FP_TOL*bin(1)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 200) )
  if (anyExceptions()) return
#line 201 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 201 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(1.4E-7_defReal,   bin(2), abs(FP_TOL*bin(2)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 202 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(0.00553_defReal,  bin(3), abs(FP_TOL*bin(3)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 202) )
  if (anyExceptions()) return
#line 203 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 203 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(10.0_defReal,     bin(4), abs(FP_TOL*bin(4)) , &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 203) )
  if (anyExceptions()) return
#line 204 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testValidBinQuery_unstruct

!!
!! Test invalid quary about the bin boundary
!!
!@test
  subroutine testInvalidBinQuery_unstruct(this)
    class(test_grid), intent(inout)          :: this
    real(defReal),dimension(3)               :: bin
    integer(shortInt),dimension(3),parameter :: idx = [-1,0, 72]
    real(defReal),parameter                  :: INVALID_BIN = -huge(bin)

    bin = this % unstructGrid % bin(idx)
#line 218 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(1), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 219 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(2), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 219) )
  if (anyExceptions()) return
#line 220 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 220 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(INVALID_BIN, bin(3), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 220) )
  if (anyExceptions()) return
#line 221 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testInvalidBinQuery_unstruct

!!
!! Test search
!!
!@test
  subroutine testSearch_unstruct(this)
    class(test_grid), intent(inout) :: this
    integer(shortInt),dimension(5)  :: idx
    real(defReal),dimension(5)      :: keys

    keys = [1.0E-11_defReal, 6.7E-4_defReal, 1.0_defReal, 9.0E-12_defReal, 21.0_defReal]
    idx = this % unstructGrid % search(keys)

#line 236 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(1,  idx(1), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 236) )
  if (anyExceptions()) return
#line 237 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 237 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(6,  idx(2), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 237) )
  if (anyExceptions()) return
#line 238 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 238 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(8,  idx(3), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 238) )
  if (anyExceptions()) return
#line 239 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

    ! Invalid Searches
#line 241 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(valueOutsideArray, idx(4), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 242 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(valueOutsideArray, idx(5), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 242) )
  if (anyExceptions()) return
#line 243 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testSearch_unstruct

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!Tests on all grids
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test getting number of bins (size)
  !!
!@Test
  subroutine testGetSize(this)
    class(test_grid), intent(inout) :: this

#line 257 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(17, this % linGrid % getSize(), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 257) )
  if (anyExceptions()) return
#line 258 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 258 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(70, this % logGrid % getSize(), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 258) )
  if (anyExceptions()) return
#line 259 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
#line 259 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"
  call assertEqual(8, this % unstructGrid % getSize(), &
 & location=SourceLocation( &
 & 'grid_test.f90', &
 & 259) )
  if (anyExceptions()) return
#line 260 "/home/mhgk4/SCONE/SharedModules/Tests/grid_test.f90"

  end subroutine testGetSize
end module grid_test

module Wrapgrid_test
   use pFUnit_mod
   use grid_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_grid) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use grid_test
        class (test_grid), intent(inout) :: this
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

end module Wrapgrid_test

function grid_test_suite() result(suite)
   use pFUnit_mod
   use grid_test
   use Wrapgrid_test
   type (TestSuite) :: suite

   suite = newTestSuite('grid_test_suite')

   call suite%addTest(makeCustomTest('testValidBinQuery_lin', testValidBinQuery_lin))

   call suite%addTest(makeCustomTest('testInvalidBinQuery_lin', testInvalidBinQuery_lin))

   call suite%addTest(makeCustomTest('testSearch_lin', testSearch_lin))

   call suite%addTest(makeCustomTest('testValidBinQuery_log', testValidBinQuery_log))

   call suite%addTest(makeCustomTest('testInvalidBinQuery_log', testInvalidBinQuery_log))

   call suite%addTest(makeCustomTest('testSearch_log', testSearch_log))

   call suite%addTest(makeCustomTest('testValidBinQuery_unstruct', testValidBinQuery_unstruct))

   call suite%addTest(makeCustomTest('testInvalidBinQuery_unstruct', testInvalidBinQuery_unstruct))

   call suite%addTest(makeCustomTest('testSearch_unstruct', testSearch_unstruct))

   call suite%addTest(makeCustomTest('testGetSize', testGetSize))


end function grid_test_suite

