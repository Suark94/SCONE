module tabularPdf_test

  use numPrecision
  use endfConstants
  use tabularPdf_class, only : tabularPdf
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_tabularPdf
    private
    type(tabularPdf) :: pdf_lin
    type(tabularPdf) :: pdf_lin_CDF
    type(tabularPdf) :: pdf_hist
    type(tabularPdf) :: pdf_hist_CDF
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_tabularPdf


contains

  !!
  !! Sets up test_tabularPdf object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_tabularPdf), intent(inout) :: this
    integer(shortInt),parameter :: R = defReal
    real(defReal),dimension(4), parameter :: Grid = [ 1.0_R, 2.0_R, 3.0_R, 4.0_R]
    real(defReal),dimension(4), parameter :: PDF_L  = [ 0.0_R, 0.5_R, 0.5_R, 0.0_R]
    real(defReal),dimension(4), parameter :: CDF_L  = [ 0.0_R, 0.25_R, 0.75_R, 1.0_R]
    real(defReal),dimension(4), parameter :: PDF_H  = [ 0.3_R, 0.5_R, 0.2_R, 0.2_R]
    real(defReal),dimension(4), parameter :: CDF_H  = [ 0.0_R, 0.3_R, 0.8_R, 1.0_R]

    ! Initialise Linear interpolation table with and without provided CDF
    ! Symmetric trapezoidal PDF between 1-4
    call this % pdf_lin % init(Grid, PDF_L, tabPdfLinLin)
    call this % pdf_lin_CDF % init(Grid, PDF_L, CDF_L, tabPdfLinLin)

    ! Initialise histogram interpolation table with and without provided CDF
    ! Assymetric histogram with 3 bins (0.3, 0.5, 0.2) between 1-4
    call this % pdf_hist % init(Grid, PDF_H, tabPdfHistogram)
    call this % pdf_hist_CDF % init(Grid, PDF_H, CDF_H, tabPdfHistogram)

  end subroutine setUp

  !!
  !! Kills test_tabularPdf object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_tabularPdf), intent(inout) :: this

    call this % pdf_lin % kill()
    call this % pdf_lin_CDF % kill()
    call this % pdf_hist % kill()
    call this % pdf_hist_CDF % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !!
  !!
!@Test
  subroutine testReadingTable(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9

    ! Test linear pdf
#line 74 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.5_defReal, this % pdf_lin % probabilityOf(2.1_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 75 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.5_defReal, this % pdf_lin_CDF % probabilityOf(2.1_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 75) )
  if (anyExceptions()) return
#line 76 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 77 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.35_defReal, this % pdf_lin % probabilityOf(1.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 78 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.35_defReal, this % pdf_lin_CDF % probabilityOf(1.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 80 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.5_defReal, this % pdf_lin % probabilityOf(3.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 81 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.5_defReal, this % pdf_lin_CDF % probabilityOf(3.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 81) )
  if (anyExceptions()) return
#line 82 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

    ! Test Histogram PDF
#line 84 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.5_defReal, this % pdf_hist % probabilityOf(2.1_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 85 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.5_defReal, this % pdf_hist_CDF % probabilityOf(2.1_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 87 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.3_defReal, this % pdf_hist % probabilityOf(1.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 87) )
  if (anyExceptions()) return
#line 88 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 88 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.3_defReal, this % pdf_hist_CDF % probabilityOf(1.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 88) )
  if (anyExceptions()) return
#line 89 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 90 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.2_defReal, this % pdf_hist % probabilityOf(3.01_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 91 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(0.2_defReal, this % pdf_hist_CDF % probabilityOf(3.01_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 91) )
  if (anyExceptions()) return
#line 92 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

  end subroutine testReadingTable

  !!
  !! Test getting bounds
  !!
!@Test
  subroutine testGettingBounds(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9
    real(defReal)                         :: bottom, top

    call this % pdf_lin % bounds(bottom, top)
#line 105 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, bottom, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 106 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, top, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

    call this % pdf_lin_CDF % bounds(bottom, top)
#line 109 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, bottom, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 109) )
  if (anyExceptions()) return
#line 110 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 110 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, top, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

    call this % pdf_hist % bounds(bottom, top)
#line 113 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, bottom, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 114 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, top, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

    call this % pdf_hist_CDF % bounds(bottom, top)
#line 117 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, bottom, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 118 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, top, TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

  end subroutine testGettingBounds

  !!
  !! Test inversion of CDF
  !!
!@Test
  subroutine testSample(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9

    ! Linear PDFs
#line 131 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.5_defReal, this % pdf_lin % sample(0.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 131) )
  if (anyExceptions()) return
#line 132 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 132 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.5_defReal, this % pdf_lin_CDF % sample(0.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 132) )
  if (anyExceptions()) return
#line 133 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 134 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.9_defReal, this % pdf_lin % sample(0.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 135 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.9_defReal, this % pdf_lin_CDF % sample(0.7_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 135) )
  if (anyExceptions()) return
#line 136 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 137 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3.367544468_defReal, this % pdf_lin % sample(0.9_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 138 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3.367544468_defReal, this % pdf_lin_CDF % sample(0.9_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 140 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.0_defReal, this % pdf_lin % sample(0.25_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 141 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.0_defReal, this % pdf_lin_CDF % sample(0.25_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 143 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, this % pdf_lin % sample(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 143) )
  if (anyExceptions()) return
#line 144 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 144 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, this % pdf_lin_CDF % sample(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 144) )
  if (anyExceptions()) return
#line 145 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 146 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, this % pdf_lin % sample(1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 147 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, this % pdf_lin_CDF % sample(1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 147) )
  if (anyExceptions()) return
#line 148 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

    ! Histograms PDFs
#line 150 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.4_defReal, this % pdf_hist % sample(0.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 150) )
  if (anyExceptions()) return
#line 151 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 151 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.4_defReal, this % pdf_hist_CDF % sample(0.5_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 153 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3.0_defReal, this % pdf_hist % sample(0.8_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 154 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3.0_defReal, this % pdf_hist_CDF % sample(0.8_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 156 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, this % pdf_hist % sample(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 157 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(1.0_defReal, this % pdf_hist_CDF % sample(0.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 159 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, this % pdf_hist % sample(1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 160 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, this % pdf_hist_CDF % sample(1.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

  end subroutine testSample

  !!
  !! Test inversion of CDF with bin
  !!
!@Test
  subroutine testSampleWithBin(this)
    class(test_tabularPdf), intent(inout) :: this
    real(defReal),parameter               :: TOL = 1.0E-9
    integer(shortInt)                     :: bin

    ! Linear PDFs
#line 174 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.5_defReal, this % pdf_lin % sample(0.5_defReal, bin), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 175 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2, bin, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 177 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3.367544468_defReal, this % pdf_lin_CDF % sample(0.9_defReal, bin), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 177) )
  if (anyExceptions()) return
#line 178 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 178 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3, bin, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 178) )
  if (anyExceptions()) return
#line 179 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

    ! Histograms PDFs
#line 181 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2.4_defReal, this % pdf_hist % sample(0.5_defReal, bin), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 181) )
  if (anyExceptions()) return
#line 182 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 182 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(2, bin, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 182) )
  if (anyExceptions()) return
#line 183 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

#line 184 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(4.0_defReal, this % pdf_hist_CDF % sample(1.0_defReal, bin), TOL, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
#line 185 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"
  call assertEqual(3, bin, &
 & location=SourceLocation( &
 & 'tabularPdf_test.f90', &
 & 185) )
  if (anyExceptions()) return
#line 186 "/home/mhgk4/SCONE/NuclearData/NuclearDataStructures/Tests/tabularPdf_test.f90"

  end subroutine testSampleWithBin


end module tabularPdf_test

module WraptabularPdf_test
   use pFUnit_mod
   use tabularPdf_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_tabularPdf) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use tabularPdf_test
        class (test_tabularPdf), intent(inout) :: this
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

end module WraptabularPdf_test

function tabularPdf_test_suite() result(suite)
   use pFUnit_mod
   use tabularPdf_test
   use WraptabularPdf_test
   type (TestSuite) :: suite

   suite = newTestSuite('tabularPdf_test_suite')

   call suite%addTest(makeCustomTest('testReadingTable', testReadingTable))

   call suite%addTest(makeCustomTest('testGettingBounds', testGettingBounds))

   call suite%addTest(makeCustomTest('testSample', testSample))

   call suite%addTest(makeCustomTest('testSampleWithBin', testSampleWithBin))


end function tabularPdf_test_suite

