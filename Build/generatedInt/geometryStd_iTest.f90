module geometryStd_iTest

  use numPrecision
  use universalVariables
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : charToDict
  use charMap_class,     only : charMap
  use dictParser_func,   only : fileToDict
  use coord_class,       only : coordList
  use geometry_inter,    only : geometry
  use geometryStd_class, only : geometryStd
  use visualiser_class,  only : visualiser
  use pFUnit_mod

  implicit none


contains

  !!
  !! Geometry integration test -> Simple 2x2 lattice
  !!
!@Test
  subroutine test_lattice_geom()
    class(geometryStd), target, allocatable   :: geom
    character(*), parameter     :: path = './IntegrationTestFiles/Geometry/test_lat'
    type(charMap)               :: mats
    integer(shortInt)           :: i, idx, matIdx, uniqueID, event
    type(dictionary)            :: dict
    real(defReal), dimension(3) :: r, u, r_ref, u_ref
    type(dictionary),pointer    :: tempDict
    character(nameLen)          :: name
    type(coordList)             :: coords
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt), dimension(10,10)           :: img
    real(defReal), dimension(6)                   :: aabb
    real(defReal)                                 :: maxDist
    class(geometry), pointer                      :: geomP
    type(visualiser)                              :: viz
    type(dictionary)                              :: vizDict
    real(defReal), parameter :: TOL = 1.0E-7_defReal
    
    ! Load dictionary
    call fileToDict(dict, path)

    ! Load materials
    tempDict => dict % getDictPtr('nuclearData')
    tempDict => tempDict % getDictPtr('materials')
    call tempDict % keys(keys, 'dict')
    do i = 1, size(keys)
      call mats % add(keys(i), i)
    end do
    
    ! Build geometry
    allocate(geom)
    call geom % init(dict, mats, silent=.true.)
    
    ! Get material at few locations
    name = 'water'
    idx = mats % get(name)
    call geom % whatIsAt(matIdx, uniqueID, [ZERO, ZERO, ZERO])
#line 62 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 62) )
  if (anyExceptions()) return
#line 63 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    name = 'mox43'
    idx = mats % get(name)
    r = [0.63_defReal, -0.09_defReal, 0.0_defReal]
    u = [ZERO, -ONE, ZERO]
    call geom % whatIsAt(matIdx, uniqueID, r, u)
#line 69 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 69) )
  if (anyExceptions()) return
#line 70 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Place coordinates
    r = [0.1_defReal, 0.1_defReal, 0.0_defReal]
    u = [ZERO, ZERO, ONE]
    call coords % init(r, u)
    call geom % placeCoord(coords)
    
    ! Verify positions
#line 78 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 79 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r, coords % lvl(2) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 79) )
  if (anyExceptions()) return
#line 80 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 80 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r - [0.63_defReal, 0.63_defReal, 0.0_defReal], coords % lvl(3) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Verify directions
#line 83 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 84 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u, coords % lvl(2) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 84) )
  if (anyExceptions()) return
#line 85 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 85 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u, coords % lvl(3) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 85) )
  if (anyExceptions()) return
#line 86 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Construct visualiser and verify slice plotting
    geomP => geom
    call charToDict(vizDict, ' ')
    call viz % init(geomP, vizDict)
    
    ! Slice plot -> Material
    call viz % slicePlot(img, [ZERO, ZERO, ZERO], 'z', 'material')
    
    ! Verify some pixels
    name = 'water'
    idx = mats % get(name)
#line 98 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, img(1, 1), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 98) )
  if (anyExceptions()) return
#line 99 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 99 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, img(2, 6), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 99) )
  if (anyExceptions()) return
#line 100 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    name = 'mox43'
    idx = mats % get(name)
#line 103 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, img(3, 7), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 103) )
  if (anyExceptions()) return
#line 104 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    name = 'uox'
    idx = mats % get(name)
#line 107 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, img(3, 3), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 107) )
  if (anyExceptions()) return
#line 108 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Slice plot -> UniqueID
    r = [-0.63_defReal, -0.63_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'z', 'uniqueID', [1.26_defReal, 1.26_defReal])

    ! Verify some pixels
    ! Note that this test depends on universe layout order in geomGraph
    ! If it changes this test fill fail
#line 116 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(2, img(5,5), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 116) )
  if (anyExceptions()) return
#line 117 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 117 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(3, img(1,1), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Verify bounds
    aabb = [-1.26_defReal, -1.26_defReal, 0.0_defReal, 1.26_defReal, 1.26_defReal, 0.0_defReal]
#line 121 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(aabb, geom % bounds(), TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    !*** Test teleport movement
    r = [ZERO, ZERO, ZERO]
    u = [-ONE, -TWO, ZERO]
    u = u/norm2(u)
    call coords % init(r, u)

    call geom % teleport(coords, 3.0_defReal)

    r_ref = [-1.1783592_defReal, -0.1632816_defReal, ZERO]
    u_ref = [ONE, -TWO, ZERO]
    u_ref = u_ref / norm2(u_ref)
    name = 'water'
    idx = mats % get(name)

#line 137 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r_ref, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 138 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u_ref, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 139 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, coords % matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 139) )
  if (anyExceptions()) return
#line 140 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
    
    !*** Test global movement
    r = [ZERO, ZERO, ZERO]
    u = [ZERO, -ONE, ZERO]
    call coords % init(r, u)

    ! Collosion movement
    maxDist = 1.0_defReal
    call geom % moveGlobal(coords, maxDist, event)

    r_ref = [ZERO, -1.0_defReal, ZERO]
    u_ref = u
    name = 'water'
    idx = mats % get(name)

#line 155 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r_ref, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 155) )
  if (anyExceptions()) return
#line 156 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 156 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u_ref, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 157 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(COLL_EV, event, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 158 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, coords % matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 158) )
  if (anyExceptions()) return
#line 159 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 159 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(1.0_defReal, maxDist, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 159) )
  if (anyExceptions()) return
#line 160 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Boundary Hit
    maxDist = 1.0_defReal
    call geom % moveGlobal(coords, maxDist, event)

    r_ref = [ZERO, 1.26_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

#line 170 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r_ref, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 171 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u_ref, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 171) )
  if (anyExceptions()) return
#line 172 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 172 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(BOUNDARY_EV, event, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 172) )
  if (anyExceptions()) return
#line 173 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 173 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, coords % matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 173) )
  if (anyExceptions()) return
#line 174 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 174 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(0.26_defReal, maxDist, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 174) )
  if (anyExceptions()) return
#line 175 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    !*** Normal Movment (easy case)
    r = [-0.63_defReal, -0.63_defReal, 0.0_defReal]
    u = [ZERO, -ONE, ZERO]
    call coords % init(r, u)
    call geom % placeCoord(coords)

    ! Local cell crossing
    maxDist = 1.0_defReal
    call geom % move(coords, maxDist, event)

    r_ref = [-0.63_defReal, -1.13_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

#line 191 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r_ref, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 191) )
  if (anyExceptions()) return
#line 192 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 192 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u_ref, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 192) )
  if (anyExceptions()) return
#line 193 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 193 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(CROSS_EV, event, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 194 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, coords % matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 194) )
  if (anyExceptions()) return
#line 195 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 195 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(0.5_defReal, maxDist, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 195) )
  if (anyExceptions()) return
#line 196 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Boundary Hit
    maxDist = 1.0_defReal
    call geom % move(coords, maxDist, event)

    r_ref = [-0.63_defReal, 1.26_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

#line 206 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r_ref, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 206) )
  if (anyExceptions()) return
#line 207 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 207 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u_ref, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 207) )
  if (anyExceptions()) return
#line 208 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 208 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(BOUNDARY_EV, event, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 208) )
  if (anyExceptions()) return
#line 209 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 209 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, coords % matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 210 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(0.13_defReal, maxDist, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 210) )
  if (anyExceptions()) return
#line 211 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Collision
    maxDist = 0.08_defReal
    call geom % move(coords, maxDist, event)

    r_ref = [-0.63_defReal, 1.18_defReal, ZERO]
    u_ref = u_ref
    name = 'water'
    idx = mats % get(name)

#line 221 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(r_ref, coords % lvl(1) % r, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 222 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(u_ref, coords % lvl(1) % dir, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 223 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(COLL_EV, event, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 223) )
  if (anyExceptions()) return
#line 224 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 224 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idx, coords % matIdx, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 224) )
  if (anyExceptions()) return
#line 225 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 225 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(0.08_defReal, maxDist, TOL, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 225) )
  if (anyExceptions()) return
#line 226 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Kill geometry
    call geom % kill()

  end subroutine test_lattice_geom

  !!
  !! Test geometry with tilted cylinder
  !!
!@Test
  subroutine test_tilted_cylinder()
    class(geometryStd), target, allocatable   :: geom
    character(*), parameter     :: path = './IntegrationTestFiles/Geometry/test_cyl'
    type(charMap)               :: mats
    integer(shortInt)           :: idxW, idxF, i
    type(dictionary)            :: dict
    type(dictionary),pointer    :: tempDict
    character(nameLen)          :: name
    character(nameLen), dimension(:), allocatable :: keys
    integer(shortInt), dimension(20,20)    :: img
    integer(shortInt), dimension(20,20,20) :: img3
    class(geometry), pointer    :: geomP
    type(visualiser)            :: viz
    type(dictionary)            :: vizDict
    real(defReal), dimension(3) :: r

    ! Load dictionary
    call fileToDict(dict, path)

    ! Load materials
    tempDict => dict % getDictPtr('nuclearData')
    tempDict => tempDict % getDictPtr('materials')
    call tempDict % keys(keys, 'dict')
    do i = 1, size(keys)
      call mats % add(keys(i), i)
    end do

    ! Build geometry
    allocate(geom)
    call geom % init(dict, mats, silent=.true.)

    ! Get fuel and water index
    name = 'water'
    idxW = mats % get(name)

    name = 'mox43'
    idxF = mats % get(name)
    
    ! Construct visualiser and verify slice plotting
    geomP => geom
    call charToDict(vizDict, ' ')
    call viz % init(geomP, vizDict)

    !*** Test slice normal to x & y
    ! X-axis at 1.0
    r = [1.0_defReal, 0.0_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'x', 'material')

    ! Test some pixels
#line 285 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxW, img(8, 11), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 285) )
  if (anyExceptions()) return
#line 286 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 286 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxW, img(17, 3), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 286) )
  if (anyExceptions()) return
#line 287 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 287 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxF, img(10, 10), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 287) )
  if (anyExceptions()) return
#line 288 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 288 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxF, img(18, 1), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 288) )
  if (anyExceptions()) return
#line 289 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
    
    ! Y-axis at 3.0
    r = [0.0_defReal, 3.0_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'y', 'material')

#line 294 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxW, img(15, 1), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 294) )
  if (anyExceptions()) return
#line 295 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 295 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxW, img(13, 4), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 295) )
  if (anyExceptions()) return
#line 296 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 296 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxF, img(13, 3), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 296) )
  if (anyExceptions()) return
#line 297 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
#line 297 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxF, img(14, 2), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 297) )
  if (anyExceptions()) return
#line 298 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    !*** Test voxel plot
    ! Full plot
    ! Value of r is irrelevant
    call viz % voxelPlot(img3, r, 'material')

    ! Checksome against 2D plot
    r = [0.0_defReal, 2.75_defReal, 0.0_defReal]
    call viz % slicePlot(img, r, 'y', 'material')

#line 308 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(img, img3(:,16,:), &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 308) )
  if (anyExceptions()) return
#line 309 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

    ! Small box all inside fuel
    r = [ 1.0_defReal, 0.0_defReal, 0.0_defReal]
    call viz % voxelPlot(img3, r, 'material', [0.5_defReal, 0.5_defReal, 0.5_defReal])

#line 314 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"
  call assertEqual(idxF, img3, &
 & location=SourceLocation( &
 & 'geometryStd_iTest.f90', &
 & 314) )
  if (anyExceptions()) return
#line 315 "/home/mhgk4/SCONE/Geometry/Tests/geometryStd_iTest.f90"

  end subroutine test_tilted_cylinder


end module geometryStd_iTest

module WrapgeometryStd_iTest
   use pFUnit_mod
   use geometryStd_iTest
   implicit none
   private

contains


end module WrapgeometryStd_iTest

function geometryStd_iTest_suite() result(suite)
   use pFUnit_mod
   use geometryStd_iTest
   use WrapgeometryStd_iTest
   type (TestSuite) :: suite

   suite = newTestSuite('geometryStd_iTest_suite')

   call suite%addTest(newTestMethod('test_lattice_geom', test_lattice_geom))

   call suite%addTest(newTestMethod('test_tilted_cylinder', test_tilted_cylinder))


end function geometryStd_iTest_suite

