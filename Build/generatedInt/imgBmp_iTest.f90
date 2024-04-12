module imgBmp_iTest

  use numPrecision
  use imgBmp_func, only : intToByte, imgBmp
  use pFUnit_mod

  implicit none


contains

  !!
  !! Test little endiness byte printng
  !!
!@Test
  subroutine test_byte_print()
    character(1) :: ref1
    character(2) :: ref2
    character(3) :: ref3
    character(4) :: ref4

    ! 1 Byte value
    ! 134
    ref1 = char(int(z'86'))

    ! 2 Byte value -> Construct byte by byte Necessary to preserve endianess
    ! 259
    ref2(1:1) = char(int(z'03'))
    ref2(2:2) = char(int(z'01'))

    ! 3 Byte value -
    ! 776655
    ref3(1:1) = char(int(z'CF'))
    ref3(2:2) = char(int(z'D9'))
    ref3(3:3) = char(int(z'0B'))

    ! 4 Byte Value
    ! 133316666
    ref4(1:1) = char(int(z'3A'))
    ref4(2:2) = char(int(z'40'))
    ref4(3:3) = char(int(z'F2'))
    ref4(4:4) = char(int(z'07'))

    ! Verify
#line 45 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
  call assertEqual(ref1, intToByte(134, 1), &
 & location=SourceLocation( &
 & 'imgBmp_iTest.f90', &
 & 45) )
  if (anyExceptions()) return
#line 46 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
#line 46 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
  call assertEqual(ref2, intToByte(259, 2), &
 & location=SourceLocation( &
 & 'imgBmp_iTest.f90', &
 & 46) )
  if (anyExceptions()) return
#line 47 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
#line 47 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
  call assertEqual(ref3, intToByte(776655, 3), &
 & location=SourceLocation( &
 & 'imgBmp_iTest.f90', &
 & 47) )
  if (anyExceptions()) return
#line 48 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
#line 48 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
  call assertEqual(ref4, intToByte(133316666, 4), &
 & location=SourceLocation( &
 & 'imgBmp_iTest.f90', &
 & 48) )
  if (anyExceptions()) return
#line 49 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"

  end subroutine test_byte_print

  !!
  !! Test creation of a BMP image by comparing the output to the reference BMP file
  !!
!@Test
  subroutine test_image()
    integer(shortInt), dimension(5,5) :: img
    integer(shortInt)                 :: black, green, red, yellow, white, blue
    integer(shortInt)                 :: file, i
    character(1)                      :: ref
    character(:), allocatable         :: image
    character(*), parameter           :: path = './IntegrationTestFiles/sample.bmp'

    ! Write colors
    black  = int(z'000000', shortInt)
    white  = int(z'FFFFFF', shortInt)
    yellow = int(z'FFFF00', shortInt)
    green  = int(z'00FF00', shortInt)
    red    = int(z'0000FF', shortInt)
    blue   = int(z'FF0000', shortInt)

    img = black
    img(1,1) = white
    img(5,1) = yellow
    img(3,3) = green
    img(1,5) = red
    img(5,5) = blue

    ! Create Image as a character string
    image = imgBmp(img)

    ! Compare byte by byte with the reference
    open(newunit=file, file=path, access='stream', status='old', action='read')

    do i = 1, len(image)
      read(file) ref
#line 87 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
  call assertEqual(ref, image(i:i), &
 & location=SourceLocation( &
 & 'imgBmp_iTest.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/UserInterface/Graphics/Tests/imgBmp_iTest.f90"
    end do

    close(file)

  end subroutine test_image

end module imgBmp_iTest

module WrapimgBmp_iTest
   use pFUnit_mod
   use imgBmp_iTest
   implicit none
   private

contains


end module WrapimgBmp_iTest

function imgBmp_iTest_suite() result(suite)
   use pFUnit_mod
   use imgBmp_iTest
   use WrapimgBmp_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('imgBmp_iTest_suite')

   call suite%addTest(newTestMethod('test_byte_print', test_byte_print))

   call suite%addTest(newTestMethod('test_image', test_image))


end function imgBmp_iTest_suite

