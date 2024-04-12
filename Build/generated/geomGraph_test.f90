module geomGraph_test

  use numPrecision
  use universalVariables, only : OUTSIDE_MAT
  use dictionary_class,   only : dictionary
  use intMap_class,       only : intMap
  use uniFills_class,     only : uniFills
  use geomGraph_class,    only : geomGraph
  use pFUnit_mod

  implicit none


  ! Variables
  type(uniFills) :: geom

contains

  !!
  !! Set up fills
  !!
  !! Filling structure (->) contained materials (|) nested universes.
  !! See definitions for order in terms of localIDs
  !!
  !!  3 -> 0
  !!  | 7
  !!  | | 1001 -> 1 4
  !!  | | 1001 -> 1 4
  !!  | | 1003 -> 1 2
  !!
!@Before
  subroutine set_up()
    integer(shortInt), dimension(:), allocatable :: fill
    type(intMap)                                 :: map

    call geom % init(7)

    ! Root uni id 3
    fill = [-7, OUTSIDE_MAT]
    call geom % addUniverse(1, 3, fill)
    call map % add(3, 1)

    ! Uni id 7
    fill = [-1001, -1001, -1003]
    call geom % addUniverse(2, 2, fill)
    call map % add(7, 2)

    ! Uni id 1001
    fill = [1, 4]
    call geom % addUniverse(3, 1001, fill)
    call map % add(1001, 3)

    ! Uni id 1002
    fill = [1, 3]
    call geom % addUniverse(4, 1002, fill)
    call map % add(1002, 4)

    ! Uni id 1003
    fill = [1, 2]
    call geom % addUniverse(5, 1003, fill)
    call map % add(1003, 5)

    ! Uni id 200
    fill = [-1001, -200]
    call geom % addUniverse(6, 200, fill)
    call map % add(200, 6)

    ! Uni id 201
    fill = [-1001, -200]
    call geom % addUniverse(7, 201, fill)
    call map % add(201, 7)

    ! Finish build
    call geom % finishBuild(map)
    call geom % setRoot(1)
  end subroutine set_up

  !!
  !! Clean up test enviroment
  !!
!@After
  subroutine clean_up()

    call geom % kill()

  end subroutine clean_up

  !!
  !! Test geometry graphs build in 'shrunk' version
  !!
!@Test
  subroutine test_shrunk()
    type(geomGraph)   :: graph
    integer(shortInt) :: i, idx, id
    type(dictionary)  :: dict

    ! Create input dictionary
    call dict % init(1)
    call dict % store('type', 'shrunk')

    call graph % kill() ! Should work fine for uninitialised as well
    call graph % init(geom, dict)

    ! Verify location array
    ! Test will fail if traverse order is change (even if structure is correct)!
#line 106 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual([-2, 0, -3, -3, -5, 1, 4, 1, 2], graph % array % idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 107 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual([ 3, 0,  6,  6,  8, 1, 2, 3, 4], graph % array % id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Verify number of unique cells
#line 110 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(4, graph % uniqueCells, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 110) )
  if (anyExceptions()) return
#line 111 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Verify used materials
#line 113 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual([1, 2, 4], graph % usedMats, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 113) )
  if (anyExceptions()) return
#line 114 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Test gtting content
    ! Universe
    call graph % getFill(idx, id, 1, 1)
#line 118 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(-2, idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 119 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(3, id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 119) )
  if (anyExceptions()) return
#line 120 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    call graph % getFill(idx, id, 3, 3)
#line 122 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(-5, idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 122) )
  if (anyExceptions()) return
#line 123 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 123 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(8, id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 123) )
  if (anyExceptions()) return
#line 124 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Material
    call graph % getFill(idx, id, 8, 2)
#line 127 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(2, idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 127) )
  if (anyExceptions()) return
#line 128 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 128 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(4, id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 128) )
  if (anyExceptions()) return
#line 129 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Clean up
    call graph % kill()

  end subroutine test_shrunk

  !!
  !! Test extended setting with instances copying
  !!
!@Test
  subroutine test_extended()
    type(geomGraph)   :: graph
    integer(shortInt) :: i, idx, id
    type(dictionary)  :: dict

    ! Create input dictionary
    call dict % init(1)
    call dict % store('type', 'extended')

    call graph % kill() ! Should work fine for uninitialised as well
    call graph % init(geom, dict)

    ! Verify location array
    ! Test will fail if traverse order is change (even if structure is correct)!
#line 153 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual([-2, 0, -3, -3, -5, 1, 4, 1, 4, 1, 2], graph % array % idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 154 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual([ 3, 0,  6,  8, 10, 1, 2, 3, 4, 5, 6], graph % array % id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 154) )
  if (anyExceptions()) return
#line 155 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Verify number of unique cells
#line 157 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(6, graph % uniqueCells, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Verify used materials
#line 160 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual([1, 2, 4], graph % usedMats, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 160) )
  if (anyExceptions()) return
#line 161 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Test gtting content
    ! Universe
    call graph % getFill(idx, id, 1, 1)
#line 165 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(-2, idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 165) )
  if (anyExceptions()) return
#line 166 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 166 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(3, id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 166) )
  if (anyExceptions()) return
#line 167 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    call graph % getFill(idx, id, 3, 3)
#line 169 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(-5, idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 169) )
  if (anyExceptions()) return
#line 170 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 170 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(10, id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Material
    call graph % getFill(idx, id, 8, 2)
#line 174 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(4, idx, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
#line 175 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"
  call assertEqual(4, id, &
 & location=SourceLocation( &
 & 'geomGraph_test.f90', &
 & 175) )
  if (anyExceptions()) return
#line 176 "/home/mhgk4/SCONE/Geometry/Tests/geomGraph_test.f90"

    ! Clean up
    call graph % kill()

  end subroutine test_extended


end module geomGraph_test

module WrapgeomGraph_test
   use pFUnit_mod
   use geomGraph_test
   implicit none
   private

contains


end module WrapgeomGraph_test

function geomGraph_test_suite() result(suite)
   use pFUnit_mod
   use geomGraph_test
   use WrapgeomGraph_test
   type (TestSuite) :: suite

   suite = newTestSuite('geomGraph_test_suite')

   call suite%addTest(newTestMethod('test_shrunk', test_shrunk, set_up, clean_up))

   call suite%addTest(newTestMethod('test_extended', test_extended, set_up, clean_up))


end function geomGraph_test_suite

