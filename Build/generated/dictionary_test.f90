module dictionary_test
  use numPrecision
  use dictionary_class, only : dictionary
  use pFUnit_mod
  implicit none

!@TestCase
  type, extends(TestCase) :: test_dictionary
    type(dictionary) :: dict
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_dictionary

  !! Parameters
  real(defReal),parameter                  :: realVal     = 3.3_defReal
  integer(shortInt),parameter              :: boolVal     = 1 
  integer(shortInt), parameter             :: intVal      = 1_shortInt
  character(nameLen),parameter             :: charNameLen = 'GoFortran_DownWithCpp'
  character(pathLen), parameter            :: charPathLen ='/home/KyloRen/VaderFanFic'
  real(defReal), dimension(2), parameter   :: realArray = [-1.0E-17_defReal, 14.7_defReal]
  integer(shortInt),dimension(3),parameter :: intArray =[-6475_shortInt, 13_shortInt, 14_shortInt]
  character(nameLen), dimension(1), parameter :: charNameLenArray = ['TK-421']
  character(pathLen), dimension(2), parameter :: charPathLenArray = ['C:\User\Tarkin\DeathStarPlans              ', &
                                                                     '/home/Dodonna/Articles/whyRebelsUseUNIX.odt']


contains

  !!
  !! Sets up test_dictionary object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_dictionary), intent(inout) :: this
    type(dictionary) :: tempDict

    call this % dict % init(1)
    call this % dict % store('myReal', realVal)
    call this % dict % store('myInt', intVal )
    call this % dict % store('myCharNameLen', charNameLen)
    call this % dict % store('myCharPathLen', charPathLen)
    call this % dict % store('realArray', realArray)
    call this % dict % store('intArray', intArray)
    call this % dict % store('charNameLenArray', charNameLenArray)
    call this % dict % store('charPathLenArray', charPathLenArray)
    call this % dict % store('myBool',boolVal)

    tempDict = this % dict
    call this % dict % store('nestedDict', tempDict)

  end subroutine setUp

  !!
  !! SKills test_dictionary object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_dictionary), intent(inout) :: this

    call this % dict % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!  TESTS PROPER BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test extracting Real value
!!
!@test
  subroutine testGettingReal(this)
    class(test_dictionary), intent(inout)    :: this
    real(defReal)                            :: tempReal

    call this % dict % get(tempReal,'myReal')
#line 76 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(realVal, tempReal, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 76) )
  if (anyExceptions()) return
#line 77 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal,'myReal',7.0_defReal)
#line 79 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(realVal, tempReal, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 79) )
  if (anyExceptions()) return
#line 80 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal,'invalid',7.0_defReal)
#line 82 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(7.0_defReal, tempReal, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 82) )
  if (anyExceptions()) return
#line 83 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Int as real
    call this % dict % get(tempReal,'myInt')
#line 86 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(real(intVal,defReal), tempReal, 'Retrival of int into real', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 86) )
  if (anyExceptions()) return
#line 87 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal,'myInt',7.0_defReal)
#line 89 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(real(intVal,defReal), tempReal, 'Get or Default Retrival of Int for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingReal

!!
!! Test extracting Real Array
!!
!@test
  subroutine testGettingRealArray(this)
    class(test_dictionary), intent(inout)    :: this
    real(defReal),dimension(:),allocatable   :: tempReal
    real(defReal),dimension(:),pointer       :: tempReal_ptr => null()

    call this % dict % get(tempReal,'realArray')
#line 103 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(realArray, tempReal, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 103) )
  if (anyExceptions()) return
#line 104 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal,'realArray', [7.0_defReal])
#line 106 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(realArray, tempReal, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal,'myRealNon',[7.0_defReal])
#line 109 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual([7.0_defReal], tempReal, 'Get or Default Retrival Failed for Absent Keyword[ptr]', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! With pointer attribute
    call this % dict % get(tempReal_ptr,'realArray')
#line 113 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(realArray, tempReal_ptr, 'Ordinary Retrival Failed[ptr]', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal_ptr,'realArray', [7.0_defReal])
#line 116 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(realArray, tempReal_ptr, 'Get or Default Retrival Failed for Present Keyword[ptr]', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 116) )
  if (anyExceptions()) return
#line 117 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal_ptr,'myRealNon',[7.0_defReal])
#line 119 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual([7.0_defReal], tempReal_ptr, 'Get or Default Retrival Failed for Absent Keyword[ptr]', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Retrival of int array into real array
    call this % dict % get(tempReal,'intArray')
#line 123 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(real(intArray,defReal), tempReal, 'Ordinary Retrival of int Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal,'intArray', [7.0_defReal])
#line 126 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(real(intArray,defReal), tempReal, 'Get or Default int Retrival for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 126) )
  if (anyExceptions()) return
#line 127 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % get(tempReal_ptr,'intArray')
#line 129 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(real(intArray,defReal), tempReal_ptr, 'Ordinary Retrival of int Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 129) )
  if (anyExceptions()) return
#line 130 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(tempReal_ptr,'intArray', [7.0_defReal])
#line 132 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(real(intArray,defReal), tempReal_ptr, 'Get or Default int Retrival for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingRealArray

!!
!! Test extracting Integer Value
!!
!@test
  subroutine testGettingInt(this)
    class(test_dictionary), intent(inout)    :: this
    integer(shortInt)                        :: temp

    call this % dict % get(temp,'myInt')
#line 145 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(intVal, temp, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'myInt',7_shortInt)
#line 148 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(intVal, temp, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', 7_shortInt)
#line 151 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(7_shortInt, temp, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingInt


!!
!! Test extracting logical Value
!!
!@test
  subroutine testGettingBool(this)
    class(test_dictionary), intent(inout)    :: this
    logical(defBool)                         :: temp

    call this % dict % get(temp,'myBool')
#line 165 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(temp, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'myBool',.false.)
#line 168 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue( temp, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 168) )
  if (anyExceptions()) return
#line 169 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', .false.)
#line 171 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertFalse(temp, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingBool

!!
!! Test extracting Integer Array
!!
!@test
  subroutine testGettingIntArray(this)
    class(test_dictionary), intent(inout)      :: this
    integer(shortInt),dimension(:),allocatable :: temp
    integer(shortInt),dimension(:),pointer     :: temp_ptr => null()

    call this % dict % get(temp,'intArray')
#line 185 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(intArray, temp, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'intArray',[7_shortInt])
#line 188 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(intArray, temp, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 188) )
  if (anyExceptions()) return
#line 189 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', [7_shortInt])
#line 191 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual([7_shortInt], temp, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! With pointer attribute
    call this % dict % get(temp_ptr,'intArray')
#line 195 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(intArray, temp_ptr, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 195) )
  if (anyExceptions()) return
#line 196 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp_ptr,'intArray',[7_shortInt])
#line 198 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(intArray, temp_ptr, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 198) )
  if (anyExceptions()) return
#line 199 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp_ptr,'invalid', [7_shortInt])
#line 201 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual([7_shortInt], temp_ptr, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 201) )
  if (anyExceptions()) return
#line 202 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingIntArray

!!
!! Test extracting nameLen long character
!!
!@test
  subroutine testGettingNameLenChar(this)
    class(test_dictionary), intent(inout) :: this
    character(nameLen)                    :: temp
    character(nameLen),parameter :: default = 'Mes Que Nada'

    call this % dict % get(temp,'myCharNameLen')
#line 215 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(charNameLen, temp, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 215) )
  if (anyExceptions()) return
#line 216 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'myCharNameLen',default)
#line 218 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(charNameLen, temp, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 218) )
  if (anyExceptions()) return
#line 219 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', default)
#line 221 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(default , temp, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingNameLenChar

!!
!! Test extracting pathLen long character
!!
!@test
  subroutine testGettingPathLenChar(this)
    class(test_dictionary), intent(inout) :: this
    character(pathLen)                    :: temp
    character(pathLen),parameter  :: default = 'Mes Que Nada'

    call this % dict % get(temp,'myCharPathLen')
#line 235 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(charPathLen, temp, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 235) )
  if (anyExceptions()) return
#line 236 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'myCharPathLen', default)
#line 238 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(charPathLen, temp, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 238) )
  if (anyExceptions()) return
#line 239 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', default)
#line 241 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(default , temp, 'Get or Default Retrival Failed for Absent Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 241) )
  if (anyExceptions()) return
#line 242 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingPathLenChar

!!
!! Test extracting nameLen long array of Chars
!!
!@test
  subroutine testGettingNameLenCharArray(this)
    class(test_dictionary), intent(inout)        :: this
    character(nameLen),dimension(:),allocatable  :: temp
    character(nameLen),dimension(:),pointer      :: temp_ptr => null()
    character(nameLen),dimension(1),parameter    :: default = ['Brasil, meu Brasil Brasileiro']
    logical(defBool)                             :: isSame

    call this % dict % get(temp,'charNameLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charNameLenArray == temp)
#line 260 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 260) )
  if (anyExceptions()) return
#line 261 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'charNameLenArray', default)
    isSame = all(charNameLenArray == temp)
#line 264 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 264) )
  if (anyExceptions()) return
#line 265 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', default)
    isSame = all(default == temp)
#line 268 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 268) )
  if (anyExceptions()) return
#line 269 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    !* With  pointer attribute
    call this % dict % get(temp_ptr,'charNameLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charNameLenArray == temp_ptr)
#line 275 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 275) )
  if (anyExceptions()) return
#line 276 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp_ptr,'charNameLenArray', default)
    isSame = all(charNameLenArray == temp_ptr)
#line 279 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 279) )
  if (anyExceptions()) return
#line 280 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp_ptr,'invalid', default)
    isSame = all(default == temp_ptr)
#line 283 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 283) )
  if (anyExceptions()) return
#line 284 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingNameLenCharArray

!!
!! Test extracting pathLen long array of Chars
!!
!@test
  subroutine testGettingPathLenCharArray(this)
    class(test_dictionary), intent(inout)        :: this
    character(pathLen),dimension(:),allocatable  :: temp
    character(pathLen),dimension(:),pointer      :: temp_ptr => null()
    character(pathLen),dimension(1),parameter    :: default = ['Brasil, meu Brasil Brasileiro']
    logical(defBool)                             :: isSame

    call this % dict % get(temp,'charPathLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charPathLenArray == temp)
#line 302 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 302) )
  if (anyExceptions()) return
#line 303 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'charPathLenArray', default)
    isSame = all(charPathLenArray == temp)
#line 306 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 306) )
  if (anyExceptions()) return
#line 307 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp,'invalid', default)
    isSame = all(default == temp)
#line 310 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 310) )
  if (anyExceptions()) return
#line 311 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    !* With pointer attribute
    call this % dict % get(temp_ptr,'charPathLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charPathLenArray == temp_ptr)
#line 317 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Ordinary Retrival Failed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 317) )
  if (anyExceptions()) return
#line 318 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp_ptr,'charPathLenArray', default)
    isSame = all(charPathLenArray == temp_ptr)
#line 321 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 321) )
  if (anyExceptions()) return
#line 322 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    call this % dict % getOrDefault(temp_ptr,'invalid', default)
    isSame = all(default == temp_ptr)
#line 325 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 325) )
  if (anyExceptions()) return
#line 326 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingPathLenCharArray

!@test
  subroutine testGettingNestedDictionary(this)
    class(test_dictionary), intent(inout) :: this
    type(dictionary)                      :: temp
    real(defReal)                         :: tempReal
    integer(shortInt)                     :: tempInt
    character(nameLen)                    :: tempCharNameLen
    character(pathLen)                    :: tempCharPathLen
    real(defReal),dimension(:),allocatable        :: tempRealArray
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    character(nameLen), dimension(:), allocatable :: tempCharArrayNameLen
    character(pathLen), dimension(:), allocatable :: tempCharArrayPathLen
    logical(defBool)                      :: isSame

    call this % dict % get(temp,'nestedDict')

    ! Get all contents of the dictionary
    call temp % get(tempReal, 'myReal')
    call temp % get(tempint, 'myInt')
    call temp % get(tempCharNameLen, 'myCharNameLen')
    call temp % get(tempCharPathLen, 'myCharPathLen')
    call temp % get(tempRealArray, 'realArray')
    call temp % get(tempIntArray, 'intArray')
    call temp % get(tempCharArrayNameLen, 'charNameLenArray')
    call temp % get(tempCharArrayPathLen, 'charPathLenArray')

    ! Verify that content was not deformed
    isSame = tempReal == realVal
    isSame = isSame .and. tempInt == intVal
    isSame = isSame .and. tempCharNameLen == charNameLen
    isSame = isSame .and. tempCharPathLen == charPathLen
    isSame = isSame .and. all(tempIntArray == intArray)
    isSame = isSame .and. all(tempRealArray == realArray)
    isSame = isSame .and. all(tempCharArrayNameLen == charNameLenArray)
    isSame = isSame .and. all(tempCharArrayPathLen == charPathLenArray)

#line 365 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Contents of nested dictionary were changed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 365) )
  if (anyExceptions()) return
#line 366 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGettingNestedDictionary

  !!
  !!
  !!
!@Test
  subroutine testPointerPassing(this)
    class(test_dictionary), intent(inout)        :: this
    type(dictionary)                      :: temp
    real(defReal)                         :: tempReal
    integer(shortInt)                     :: tempInt
    character(nameLen)                    :: tempCharNameLen
    character(pathLen)                    :: tempCharPathLen
    real(defReal),dimension(:),allocatable        :: tempRealArray
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    character(nameLen), dimension(:), allocatable :: tempCharArrayNameLen
    character(pathLen), dimension(:), allocatable :: tempCharArrayPathLen
    logical(defBool)                      :: isSame

    ! Retrieve pointer in a subroutine
    ! It is to test whether pointer to a nested dictionary
    ! going out of scope upon subroutine termination
    ! causes deallocation of memory.
    !
    call getAndFinalPointer(this % dict)

    call this % dict % get(temp,'nestedDict')

    ! Get all contents of the dictionary
    call temp % get(tempReal, 'myReal')
    call temp % get(tempint, 'myInt')
    call temp % get(tempCharNameLen, 'myCharNameLen')
    call temp % get(tempCharPathLen, 'myCharPathLen')
    call temp % get(tempRealArray, 'realArray')
    call temp % get(tempIntArray, 'intArray')
    call temp % get(tempCharArrayNameLen, 'charNameLenArray')
    call temp % get(tempCharArrayPathLen, 'charPathLenArray')

    ! Verify that content was not deformed
    isSame = tempReal == realVal
    isSame = isSame .and. tempInt == intVal
    isSame = isSame .and. tempCharNameLen == charNameLen
    isSame = isSame .and. tempCharPathLen == charPathLen
    isSame = isSame .and. all(tempIntArray == intArray)
    isSame = isSame .and. all(tempRealArray == realArray)
    isSame = isSame .and. all(tempCharArrayNameLen == charNameLenArray)
    isSame = isSame .and. all(tempCharArrayPathLen == charPathLenArray)

#line 415 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isSame, 'Contents of nested dictionary were changed', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 415) )
  if (anyExceptions()) return
#line 416 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testPointerPassing

  !!
  !!
  !!
  subroutine getAndFinalPointer(dict)
    class(dictionary), intent(in) :: dict
    class(dictionary),pointer     :: ptr
    logical(defBool)              :: myBool

    ptr => dict % getDictPtr('nestedDict')
    myBool = ptr % isPresent('bla')

  end subroutine getAndFinalPointer

!!
!! Test keys retrival
!!
!@test
  subroutine testKeys(this)
    class(test_dictionary), intent(inout)        :: this
    character(nameLen),dimension(:),allocatable  :: tempAll
    character(nameLen),dimension(:),allocatable  :: tempReal
    character(nameLen),dimension(:),allocatable  :: tempInt
    character(nameLen),dimension(:),allocatable  :: tempChar
    character(nameLen),dimension(:),allocatable  :: tempRealArray
    character(nameLen),dimension(:),allocatable  :: tempIntArray
    character(nameLen),dimension(:),allocatable  :: tempCharArray
    character(nameLen),dimension(:),allocatable  :: tempDict
    character(nameLen)                           :: keyword
    logical(defBool) :: isValid

    ! Obtain all keys arrays
    call this % dict % keys(tempAll)
    call this % dict % keys(tempReal,'real')
    call this % dict % keys(tempInt,'int')
    call this % dict % keys(tempChar,'char')
    call this % dict % keys(tempDict,'dict')

    call this % dict % keys(tempRealArray,'realArray')
    call this % dict % keys(tempIntArray,'intArray')
    call this % dict % keys(tempCharArray,'charArray')

    ! Verify keys for real
    keyword = 'myReal'
    isValid = contains(tempAll, keyword) .and. contains(tempReal, keyword)
#line 463 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid, 'Keywords failed for real', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 463) )
  if (anyExceptions()) return
#line 464 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Verify keys for int
    keyword = 'myInt'
    isValid = contains(tempAll, keyword) .and. contains(tempInt, keyword)
#line 468 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid,' Keywords failed for int', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 468) )
  if (anyExceptions()) return
#line 469 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Verify keys for char
    keyword = 'myCharNameLen'
    isValid = contains(tempAll, keyword) .and. contains(tempChar, keyword)
    keyword = 'myCharPathLen'
    isValid = isValid .and. contains(tempAll, keyword) .and. contains(tempChar, keyword)
#line 475 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid,' Keywords failed for char', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 475) )
  if (anyExceptions()) return
#line 476 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Verify keys for dict
    keyword = 'nestedDict'
    isValid = contains(tempAll, keyword) .and. contains(tempDict, keyword)
#line 480 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid,' Keywords failed for nested dictionary', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 480) )
  if (anyExceptions()) return
#line 481 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Verify keys for realArray
    keyword = 'realArray'
    isValid = contains(tempAll, keyword) .and. contains(tempRealArray, keyword)
#line 485 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid, 'Keywords failed for real array ', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 485) )
  if (anyExceptions()) return
#line 486 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Verify keys for int
    keyword = 'intArray'
    isValid = contains(tempAll, keyword) .and. contains(tempIntArray, keyword)
#line 490 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid,' Keywords failed for int array', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 490) )
  if (anyExceptions()) return
#line 491 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Verify keys for char
    keyword = 'charNameLenArray'
    isValid = contains(tempAll, keyword) .and. contains(tempCharArray, keyword)
    keyword = 'charPathLenArray'
    isValid = isValid .and. contains(tempAll, keyword) .and. contains(tempCharArray, keyword)
#line 497 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isValid,' Keywords failed for char array', &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 497) )
  if (anyExceptions()) return
#line 498 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"


  contains
    function contains(array,keyword) result(doesIt)
      character(nameLen), dimension(:) :: array
      character(*)                     :: keyword
      logical(defBool)                 :: doesIt

      doesIt = count(array == keyword) == 1

    end function contains
  end subroutine testKeys

  !!
  !! Test isPresent function of a dictionary
  !!
!@test
  subroutine testIsPresent(this)
    class(test_dictionary), intent(inout) :: this
    logical(defBool)                      :: isPresent

    isPresent = this % dict % isPresent('nestedDict')
#line 520 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertTrue(isPresent, &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 520) )
  if (anyExceptions()) return
#line 521 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    isPresent = this % dict % isPresent('invalid')
#line 523 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertFalse(isPresent, &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 523) )
  if (anyExceptions()) return
#line 524 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testIsPresent

  !!
  !! Test getSize function of a dictionary
  !!
!@test
  subroutine testGetSize(this)
    class(test_dictionary), intent(inout) :: this

    ! Get size of scalar
#line 535 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(1, this % dict % getSize('myInt'), &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 535) )
  if (anyExceptions()) return
#line 536 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Get size of int Array
#line 538 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(3, this % dict % getSize('intArray'), &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 538) )
  if (anyExceptions()) return
#line 539 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Get size of realArray
#line 541 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(2, this % dict % getSize('realArray'), &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 541) )
  if (anyExceptions()) return
#line 542 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Get size of word Array
#line 544 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(1, this % dict % getSize('charNameLenArray'), &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 544) )
  if (anyExceptions()) return
#line 545 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

    ! Get length of the dictionary (number of entries)
#line 547 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"
  call assertEqual(10, this % dict % length(), &
 & location=SourceLocation( &
 & 'dictionary_test.f90', &
 & 547) )
  if (anyExceptions()) return
#line 548 "/home/mhgk4/SCONE/DataStructures/Tests/dictionary_test.f90"

  end subroutine testGetSize


end module dictionary_test

module Wrapdictionary_test
   use pFUnit_mod
   use dictionary_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_dictionary) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use dictionary_test
        class (test_dictionary), intent(inout) :: this
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

end module Wrapdictionary_test

function dictionary_test_suite() result(suite)
   use pFUnit_mod
   use dictionary_test
   use Wrapdictionary_test
   type (TestSuite) :: suite

   suite = newTestSuite('dictionary_test_suite')

   call suite%addTest(makeCustomTest('testGettingReal', testGettingReal))

   call suite%addTest(makeCustomTest('testGettingRealArray', testGettingRealArray))

   call suite%addTest(makeCustomTest('testGettingInt', testGettingInt))

   call suite%addTest(makeCustomTest('testGettingBool', testGettingBool))

   call suite%addTest(makeCustomTest('testGettingIntArray', testGettingIntArray))

   call suite%addTest(makeCustomTest('testGettingNameLenChar', testGettingNameLenChar))

   call suite%addTest(makeCustomTest('testGettingPathLenChar', testGettingPathLenChar))

   call suite%addTest(makeCustomTest('testGettingNameLenCharArray', testGettingNameLenCharArray))

   call suite%addTest(makeCustomTest('testGettingPathLenCharArray', testGettingPathLenCharArray))

   call suite%addTest(makeCustomTest('testGettingNestedDictionary', testGettingNestedDictionary))

   call suite%addTest(makeCustomTest('testPointerPassing', testPointerPassing))

   call suite%addTest(makeCustomTest('testKeys', testKeys))

   call suite%addTest(makeCustomTest('testIsPresent', testIsPresent))

   call suite%addTest(makeCustomTest('testGetSize', testGetSize))


end function dictionary_test_suite

