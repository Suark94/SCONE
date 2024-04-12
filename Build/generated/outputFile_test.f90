module outputFile_test

  use numPrecision
  use outputFile_class, only : outputFile
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_outputFile
    private
    type(outputFile) :: outFile
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_outputFile


contains

  !!
  !! Sets up test_outputFile object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_outputFile), intent(inout) :: this

    call this % outFile % init('dummyPrinter', fatalErrors = .false.)

  end subroutine setUp

  !!
  !! Kills test_outputFile object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_outputFile), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that block errors are caught
  !!
!@Test
  subroutine testBlockLogic(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen)                    :: myBlock, myBlock2, nextBlock, myArray

    myBlock = 'myBlock'
    myBlock2 = 'myBlock2'
    nextBlock = 'nextBlock'
    myArray = 'myArray'

    ! Test for too early block closure
    call this % outFile % endBlock()
#line 58 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 58) )
  if (anyExceptions()) return
#line 59 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Test for to many block closures
    call this % outFile % startBlock(myBlock)
    call this % outFile % startBlock(myBlock2)
    call this % outFile % endBlock()
    call this % outFile % endBlock()
    call this % outFile % endBlock()
#line 67 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Try to close block from an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2,3])
    call this % outFile % endBlock()
#line 74 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Try to start block from an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2,3])
    call this % outFile % startBlock(nextBlock)
#line 81 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

  end subroutine testBlockLogic

  !!
  !! Test that Array Errors are caught
  !!
!@Test
  subroutine testArrayLogic(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen), dimension(2) :: charArray
    character(nameLen)               :: name
    character(nameLen)               :: myBlock, myArray, myArray2
    integer(shortInt), dimension(:), allocatable :: temp_int

    myBlock = 'myBlock'
    myArray = 'myArray'
    myArray2 = 'myArray2'

    ! Start array in array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2,3])
    call this % outFile % startArray(myArray2,[2])
#line 105 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Change from result to value
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    call this % outFile % addValue(0.0_defReal)
#line 113 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Change type value in value array Real -> Int
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % addValue(0)
#line 121 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Change type value in value array Real -> Char
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addValue(0.0_defReal)
    name = 'jgvjj'
    call this % outFile % addValue(name)
#line 130 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Change type value in value array Int -> Real
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[3])
    call this % outFile % addValue(0)
    call this % outFile % addValue(0.0_defReal)
#line 138 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Array overflow defReal
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % addValue(0.0_defReal)
#line 147 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 147) )
  if (anyExceptions()) return
#line 148 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Array overflow shortInt
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0)
    call this % outFile % addValue(0)
    call this % outFile % addValue(0)
#line 156 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Array overflow longInt
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0_longInt)
    call this % outFile % addValue(0_longInt)
    call this % outFile % addValue(0_longInt)
#line 165 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Array overflow char
    name ='lkm'
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(name)
    call this % outFile % addValue(name)
    call this % outFile % addValue(name)
#line 175 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Array overflow Result
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
    call this % outFile % addResult(0.0_defReal, 1.0_defReal)
#line 184 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Close array with too little entries provided
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    call this % outFile % endArray()
#line 192 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Printing shortInt value inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(0, name)
#line 201 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Printing longInt value inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(0_longInt, name)
#line 210 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 210) )
  if (anyExceptions()) return
#line 211 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Printing defReal value inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(0.0_defReal, name)
#line 219 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 219) )
  if (anyExceptions()) return
#line 220 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Printing character value inside an array
    charArray(1) = 'sth'
    charArray(2) ='sth else'
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printValue(charArray(1), name)
#line 230 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Printing result inside an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % startArray(myArray,[2])
    call this % outFile % addValue(0.0_defReal)
    name ='myVal'
    call this % outFile % printResult(0.0_defReal,1.0_defReal, name)
#line 239 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 239) )
  if (anyExceptions()) return
#line 240 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Add defReal value without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(0.0_defReal)
#line 245 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 245) )
  if (anyExceptions()) return
#line 246 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Add shortInt value without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(0)
#line 251 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 251) )
  if (anyExceptions()) return
#line 252 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Add longInt value without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(0_longInt)
#line 257 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 257) )
  if (anyExceptions()) return
#line 258 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Add char value without starting an array
    name ='char value'
    call this % outFile % startBlock(myBlock)
    call this % outFile % addValue(name)
#line 264 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 264) )
  if (anyExceptions()) return
#line 265 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Add result without starting an array
    call this % outFile % startBlock(myBlock)
    call this % outFile % addResult(0.0_defReal, 0.0_defReal)
#line 270 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 270) )
  if (anyExceptions()) return
#line 271 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Create a degenerate arrays
    allocate(temp_int(0))
    call this % outFile % startArray(name, temp_int)
    call this % outFile % endArray()
#line 277 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 277) )
  if (anyExceptions()) return
#line 278 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    call this % outFile % startArray(name, [2, 0, 7])
    call this % outFile % endArray()
#line 282 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 282) )
  if (anyExceptions()) return
#line 283 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

  end subroutine testArrayLogic

  !!
  !! Test Repeated Names
  !!
!@Test
  subroutine testRepeatedNames(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen)                    :: name

    name = 'myKey'

    ! Print ordinary values
    call this % outFile % printResult(1.0_defReal, 0.5_defReal, name)
    call this % outFile % printResult(1.0_defReal, 1.5_defReal, name)
#line 300 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 300) )
  if (anyExceptions()) return
#line 301 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    call this % outFile % printValue(1.0_defReal, name)
    call this % outFile % printValue(2.0_defReal, name)
#line 305 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 305) )
  if (anyExceptions()) return
#line 306 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    call this % outFile % printValue(1, name)
    call this % outFile % printValue(2, name)
#line 310 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 310) )
  if (anyExceptions()) return
#line 311 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    call this % outFile % printValue(1_longInt, name)
    call this % outFile % printValue(2_longInt, name)
#line 315 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 315) )
  if (anyExceptions()) return
#line 316 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    call this % outFile % printValue("char", name)
    call this % outFile % printValue("charizard", name)
#line 320 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 320) )
  if (anyExceptions()) return
#line 321 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    ! Test with blocks & arrays
    call this % outFile % startBlock(name)
    call this % outFile % endBlock()
    call this % outFile % startBlock(name)
    call this % outFile % endBlock()
#line 328 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 328) )
  if (anyExceptions()) return
#line 329 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

    call this % outFile % startBlock(name)
    call this % outFile % endBlock()
    call this % outFile % startArray(name, [1])
    call this % outFile % addValue(1)
    call this % outFile % endArray()
#line 336 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertFalse( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 336) )
  if (anyExceptions()) return
#line 337 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

  end subroutine testRepeatedNames

  !!
  !! Test deep nested blocks
  !!
!@Test
  subroutine testNestedBlocks(this)
    class(test_outputFile), intent(inout) :: this
    character(nameLen)                    :: name
    integer(shortInt)                     :: i

    name = 'myKey'

    do i = 1, 30
      call this % outFile % startBlock(name)
    end do
    do i = 1, 30
      call this % outFile % endBlock()
    end do

#line 359 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
  call assertTrue( this % outFile % isValid(), &
 & location=SourceLocation( &
 & 'outputFile_test.f90', &
 & 359) )
  if (anyExceptions()) return
#line 360 "/home/mhgk4/SCONE/UserInterface/fileOutput/Tests/outputFile_test.f90"
    call this % outFile % reset()

  end subroutine testNestedBlocks


end module outputFile_test

module WrapoutputFile_test
   use pFUnit_mod
   use outputFile_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_outputFile) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use outputFile_test
        class (test_outputFile), intent(inout) :: this
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

end module WrapoutputFile_test

function outputFile_test_suite() result(suite)
   use pFUnit_mod
   use outputFile_test
   use WrapoutputFile_test
   type (TestSuite) :: suite

   suite = newTestSuite('outputFile_test_suite')

   call suite%addTest(makeCustomTest('testBlockLogic', testBlockLogic))

   call suite%addTest(makeCustomTest('testArrayLogic', testArrayLogic))

   call suite%addTest(makeCustomTest('testRepeatedNames', testRepeatedNames))

   call suite%addTest(makeCustomTest('testNestedBlocks', testNestedBlocks))


end function outputFile_test_suite

