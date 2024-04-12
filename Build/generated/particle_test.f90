module particle_test
  use numPrecision
  use particle_class, only : particle, particleState, P_NEUTRON, P_PHOTON, verifyType
  use pFUnit_mod

  implicit none

!@TestCase
  type, extends(TestCase) :: test_particle
    type(particle) :: p_MG
    type(particle) :: p_CE
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_particle

  ! Test paramethers
  real(defReal),dimension(3),parameter :: r0 = [-1.3_defReal, 97.7_defReal, -9.0_defReal]
  real(defReal),dimension(3),parameter :: u0 = [1.0_defReal, 0.0_defReal, 0.0_defReal]
  real(defReal),parameter              :: t0 = 13.0_defReal
  real(defReal),parameter              :: w0 = 0.7_defReal
  real(defReal),parameter              :: E0 = 0.00034_defReal
  integer(shortInt),parameter          :: G0 = 3

  real(defReal),dimension(3),parameter :: lev1_offset = [0.3_defReal, -13.5_defReal, 1.0_defReal]
  real(defReal),dimension(3),parameter :: lev2_offset = [-3.4_defReal, 3.0_defReal, -8.0_defReal]

  integer(shortInt),parameter          :: lev1_uni     = 5
  integer(shortInt),parameter          :: lev1_uniRoot = 67
  integer(shortInt),parameter          :: lev2_uni     = 3
  integer(shortInt),parameter          :: lev2_uniRoot = 32

contains

  !!
  !! Set-up particle for testing
  !!
  subroutine setUp(this)
    class(test_particle), intent(inout) :: this

    ! Build MG particle
    call this % p_MG % build(r0, u0, G0, w0, t0)
    this % p_MG % type = P_NEUTRON

    ! Build CE particle
    call this % p_CE % build(r0, u0, E0, w0)
    this % p_CE % type = P_PHOTON

    ! Add 2 levels to the particles
    call this % p_MG % coords % addLevel()
    this % p_MG % coords % lvl(2) % r = r0 - lev1_offset
    this % p_MG % coords % lvl(2) % dir = u0

    call this % p_CE % coords % addLevel()
    this % p_CE % coords % lvl(2) % r = r0 - lev1_offset
    this % p_CE % coords % lvl(2) % dir = u0

    call this % p_MG % coords % addLevel()
    this % p_MG % coords % lvl(3) % r = r0 - lev1_offset - lev2_offset
    this % p_MG % coords % lvl(3) % dir = u0

    call this % p_CE % coords % addLevel()
    this % p_CE % coords % lvl(3) % r = r0 - lev1_offset - lev2_offset
    this % p_CE % coords % lvl(3) % dir = u0

    ! Set MatIdx
    this % p_MG % coords % matIdx = 7
    this % p_CE % coords % matIdx = 7

    ! Set uniqueID
    this % p_MG % coords % uniqueID = 34
    this % p_CE % coords % uniqueID = 34

    ! NOTE: THIS MAY BREAK AT SOME POINT
    ! HAND SET SOME COORD PARAMETERS
    ! BREAKS ENCAPSULATION (It is broken anyway)

    ! Level 1
    this % p_MG % coords % lvl(1) % uniIdx    = 1
    this % p_MG % coords % lvl(1) % uniRootID = 0
    this % p_MG % coords % lvl(1) % localID   = 1
    this % p_MG % coords % lvl(1) % cellIdx   = 3

    this % p_CE % coords % lvl(1) % uniIdx    = 1
    this % p_CE % coords % lvl(1) % uniRootID = 0
    this % p_CE % coords % lvl(1) % localID   = 1
    this % p_CE % coords % lvl(1) % cellIdx   = 3

    ! Level 2
    this % p_MG % coords % lvl(2) % uniIdx    = lev1_uni
    this % p_MG % coords % lvl(2) % uniRootID = lev1_uniRoot
    this % p_MG % coords % lvl(2) % localID   = 4
    this % p_MG % coords % lvl(2) % cellIdx   = 2

    this % p_CE % coords % lvl(2) % uniIdx    = lev1_uni
    this % p_CE % coords % lvl(2) % uniRootID = lev1_uniRoot
    this % p_CE % coords % lvl(2) % localID   = 4
    this % p_CE % coords % lvl(2) % cellIdx   = 2

    ! Level 3
    this % p_MG % coords % lvl(3) % uniIdx    = lev2_uni
    this % p_MG % coords % lvl(3) % uniRootID = lev2_uniRoot
    this % p_MG % coords % lvl(3) % localID   = 2
    this % p_MG % coords % lvl(3) % cellIdx   = 8

    this % p_CE % coords % lvl(3) % uniIdx    = lev2_uni
    this % p_CE % coords % lvl(3) % uniRootID = lev2_uniRoot
    this % p_CE % coords % lvl(3) % localID   = 2
    this % p_CE % coords % lvl(3) % cellIdx   = 8

  end subroutine setUp

  !!
  !! Kill particle after test
  !!
  subroutine tearDown(this)
    class(test_particle), intent(inout) :: this
    ! DO NOTHING
  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test initialisation correctness
  !!
!@test
  subroutine correctInitialisation(this)
    class(test_particle), intent(inout) :: this

    ! Test energies
#line 133 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(E0, this % p_CE % E, 'CE energy value', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 133) )
  if (anyExceptions()) return
#line 134 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 134 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(G0, this % p_MG % G, 'MG group number', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test time
#line 137 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(ZERO, this % p_CE % time ,'Default time initialisation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 138 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(t0, this % p_MG % time ,'Time initialisation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 138) )
  if (anyExceptions()) return
#line 139 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test weight
#line 141 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, this % p_CE % w, 'CE weight initialisation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 142 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, this % p_CE % w0, 'CE initial weigth initialisation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 142) )
  if (anyExceptions()) return
#line 143 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

#line 144 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, this % p_MG % w, 'MG weight initialisation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 144) )
  if (anyExceptions()) return
#line 145 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 145 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, this % p_MG % w0, 'MG initial weigth initialisation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test isMG flag
#line 148 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertFalse( this % p_CE % isMG, 'isMG flag for CE particle', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 148) )
  if (anyExceptions()) return
#line 149 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 149 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue( this % p_MG % isMG, 'isMG flag for MG particle', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 149) )
  if (anyExceptions()) return
#line 150 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test isDead flag
#line 152 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertFalse( this % p_CE % isDead, 'isDead flag CE particle', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 152) )
  if (anyExceptions()) return
#line 153 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 153 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertFalse( this % p_MG % isDead, 'isDead flag MG particle', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 153) )
  if (anyExceptions()) return
#line 154 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test timeMax default value
#line 156 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(ZERO, this % p_CE % timeMax, 'timeMax initialises to 0 by default', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 156) )
  if (anyExceptions()) return
#line 157 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 157 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(ZERO, this % p_MG % timeMax, 'timeMax initialises to 0 by default', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 157) )
  if (anyExceptions()) return
#line 158 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine correctInitialisation

  !!
  !! Test access to positions
  !!
!@test
  subroutine testPositionAccess(this)
    class(test_particle), intent(inout) :: this
    real(defReal),dimension(3)          :: r_global
    real(defReal),dimension(3)          :: r_lev1
    real(defReal),dimension(3)          :: r_lev2
    real(defReal),dimension(3)          :: r_lev3

    ! Get position at all levels
    r_global = this % p_CE % rGlobal()

    r_lev1 = this % p_CE % rLocal(1)
    r_lev2 = this % p_CE % rLocal(2)
    r_lev3 = this % p_CE % rLocal()

    ! Verify correctness
#line 180 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0,r_global,'Global Position Wrong', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 181 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(r_global == r_lev1, 'Level 1 r is not equal global r', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 181) )
  if (anyExceptions()) return
#line 182 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 182 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0-lev1_offset, r_lev2, 'Level 2 position is wrong', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 182) )
  if (anyExceptions()) return
#line 183 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 183 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0-lev1_offset-lev2_offset, r_lev3, 'Level 3 position is wrong', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 183) )
  if (anyExceptions()) return
#line 184 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testPositionAccess

  !!
  !! Test access to directions
  !!
!@test
  subroutine testDirectionAccess(this)
    class(test_particle), intent(inout) :: this
    real(defReal),dimension(3)          :: u_global
    real(defReal),dimension(3)          :: u_lev1
    real(defReal),dimension(3)          :: u_lev2
    real(defReal),dimension(3)          :: u_lev3

    ! Get direction at all levels
    u_global = this % p_CE % dirGlobal()

    u_lev1 = this % p_CE % dirLocal(1)
    u_lev2 = this % p_CE % dirLocal(2)
    u_lev3 = this % p_CE % dirLocal()

    ! Verify correctness
#line 206 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(u0,u_global,'Global Direction Wrong', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 206) )
  if (anyExceptions()) return
#line 207 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 207 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(u_global == u_lev1, 'Level 1 dir is not equal global dir', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 207) )
  if (anyExceptions()) return
#line 208 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 208 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(u0, u_lev2, 'Level 2 Direction is wrong', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 208) )
  if (anyExceptions()) return
#line 209 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 209 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(u0, u_lev3, 'Level 3 Direction is wrong', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 209) )
  if (anyExceptions()) return
#line 210 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testDirectionAccess

  !!
  !! Tests access procedures to extra parameters. TODO: add uniqueID
  !!
!@test
  subroutine testMiscAccess(this)
    class(test_particle), intent(inout) :: this
    integer(shortInt)                   :: matIdx, cellIdx, uniIdx

    ! Verify nesting level
#line 222 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(3, this % p_CE % nesting(), 'Nesting Level', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 222) )
  if (anyExceptions()) return
#line 223 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Level 3
    matIdx  = this % p_CE % matIdx()
    cellIdx = this % p_CE % getCellIdx()
    uniIdx  = this % p_CE % getUniIdx()

    ! Verify correctness
#line 230 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(7, matIdx, 'Material Index', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 230) )
  if (anyExceptions()) return
#line 231 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 231 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(8, cellIdx, 'Cell Index. Deepest level.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 231) )
  if (anyExceptions()) return
#line 232 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 232 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(3, uniIdx, 'Universe Index. Deepest level.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 232) )
  if (anyExceptions()) return
#line 233 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Level 2
    cellIdx = this % p_CE % getCellIdx(2)
    uniIdx  = this % p_CE % getUniIdx(2)

    ! Verify correctness
#line 239 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(2, cellIdx, 'Cell Index. Level 2.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 239) )
  if (anyExceptions()) return
#line 240 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 240 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(5, uniIdx, 'Universe Index. Level 2.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 240) )
  if (anyExceptions()) return
#line 241 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Level 1
    cellIdx = this % p_CE % getCellIdx(1)
    uniIdx  = this % p_CE % getUniIdx(1)

    ! Verify correctness
#line 247 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(3, cellIdx, 'Cell Index. Level 1.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 247) )
  if (anyExceptions()) return
#line 248 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 248 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(1, uniIdx, 'Universe Index. Level 1.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 248) )
  if (anyExceptions()) return
#line 249 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testMiscAccess

  !!
  !! Test material setting outside geometry
  !!
!@test
  subroutine testSetMatIdx(this)
    class(test_particle), intent(inout) :: this

    call this % p_CE % setMatIdx(3)

#line 261 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(3, this % p_CE % matIdx(), &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 261) )
  if (anyExceptions()) return
#line 262 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testSetMatIdx


  !!
  !! Test movement procedures
  !!
!@test
  subroutine testMovementProcedures(this)
    class(test_particle), intent(inout) :: this
    real(defReal)                       :: dist = 1.0_defReal
    real(defReal),dimension(3)          :: r0_lvl1 ,r0_lvl2, r0_lvl3
    real(defReal),dimension(3)          :: r_lvl1, r_lvl2, r_lvl3

    ! Move local on lowest level
    call this % p_CE % moveLocal(dist,3)

    r_lvl1 = this % p_CE % rLocal(1)
    r_lvl2 = this % p_CE % rLocal(2)
    r_lvl3 = this % p_CE % rLocal()

    ! Calculate expected position
    r0_lvl1 = r0 + u0 * dist
    r0_lvl2 = r0 - lev1_offset + u0 * dist
    r0_lvl3 = r0 - lev1_offset - lev2_offset + u0 * dist

    ! Verify position
#line 289 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0_lvl1, r_lvl1, 'Level 1 position. Local Movement.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 289) )
  if (anyExceptions()) return
#line 290 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 290 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0_lvl2, r_lvl2, 'Level 2 position. Local Movement.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 290) )
  if (anyExceptions()) return
#line 291 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 291 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0_lvl3, r_lvl3, 'Level 3 position. Local Movement.', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 291) )
  if (anyExceptions()) return
#line 292 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 292 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(this % p_CE % coords % isPlaced(), 'Particle is placed in the geometry', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 292) )
  if (anyExceptions()) return
#line 293 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Move on global level
    call this % p_CE % moveGlobal(dist)
    r_lvl1  = this % p_CE % rGlobal()
    r0_lvl1 = r0_lvl1 + u0 * dist

#line 299 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0_lvl1, r_lvl1, 'Global position after global movement', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 299) )
  if (anyExceptions()) return
#line 300 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 300 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(this % p_CE % coords % isAbove(), 'Particle is above geometry', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 300) )
  if (anyExceptions()) return
#line 301 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Teleport on global level
    r0_lvl1 = [-1.0_defReal, 6.0_defReal, -14.6868_defReal]
    call this % p_CE % teleport(r0_lvl1)
    r_lvl1 = this % p_CE % rGlobal()

#line 307 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0_lvl1, r_lvl1, 'Global position after global teleport', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 307) )
  if (anyExceptions()) return
#line 308 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 308 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(this % p_CE % coords % isAbove(), 'Particle is above geometry', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 308) )
  if (anyExceptions()) return
#line 309 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  end subroutine testMovementProcedures

  !!
  !!
  !!
!@test
  subroutine testRotationProcedures(this)
    class(test_particle), intent(inout) :: this
    real(defReal), dimension(3)         :: dir
    real(defReal)                       :: tol

    dir = [4.0_defReal, -5.0_defReal, 1.0_defReal]
    dir = dir / norm2(dir)

    ! Test point procedure
    call this % p_CE % point(dir)

#line 326 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirGlobal(), 'Global direction after pointing', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 326) )
  if (anyExceptions()) return
#line 327 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 327 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirLocal(1),' Local direction after pointing. Level 1', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 327) )
  if (anyExceptions()) return
#line 328 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 328 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirLocal(2),' Local direction after pointing. Level 2', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 328) )
  if (anyExceptions()) return
#line 329 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 329 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirLocal(),' Local direction after pointing. Level 3', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 329) )
  if (anyExceptions()) return
#line 330 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 330 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(this % p_CE % coords % isPlaced(), 'Particle is still placed in geometry', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 330) )
  if (anyExceptions()) return
#line 331 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test rotation
    call this % p_CE % rotate(0.3_defReal, 1.3_defReal)

    ! Rotation was performed with an independent, verified MATLAB implementation
    dir = [0.927517049363521_defReal, 0.312003100719400_defReal,   -0.205830484334726_defReal]
    tol = 50.0 * epsilon(dir)
#line 338 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirGlobal(), tol, 'Global direction after rotation', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 338) )
  if (anyExceptions()) return
#line 339 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 339 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirLocal(1), tol,'Local direction after rotation. Level 1', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 339) )
  if (anyExceptions()) return
#line 340 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 340 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirLocal(2), tol,'Local direction after rotation. Level 2', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 340) )
  if (anyExceptions()) return
#line 341 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 341 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(dir, this % p_CE % dirLocal(), tol,'Local direction after rotation. Level 3', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 341) )
  if (anyExceptions()) return
#line 342 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 342 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(this % p_CE % coords % isPlaced(),'Particle is still placed in geometry', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 342) )
  if (anyExceptions()) return
#line 343 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Test taking above geometry
    call this % p_CE % takeAboveGeom()
#line 346 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(this % p_CE % coords % isAbove(), 'Particle is above geometry', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 346) )
  if (anyExceptions()) return
#line 347 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testRotationProcedures

  !!
  !! Test that assignments between particle and particleState are correct
  !!
!@test
  subroutine testParticleStateAssignments(this)
    class(test_particle), intent(inout) :: this
    type(particleState)                 :: mg_coord
    type(particleState)                 :: ce_coord
    type(particle)                      :: mg_p
    type(particle)                      :: ce_p

    ! Copy particle to phaseCoords
    mg_coord = this % p_MG
    ce_coord = this % p_CE

    ! Copy phaseCoords to particle
    mg_p = mg_coord
    ce_p = ce_coord

    ! Verify correctness (Use only MG execept where it differs from CE: energy and isMG flag)

    ! Compare global position and direction
#line 372 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0, mg_coord % r, 'Global position. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 372) )
  if (anyExceptions()) return
#line 373 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 373 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(r0, mg_p % rGlobal(), 'Global position. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 373) )
  if (anyExceptions()) return
#line 374 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

#line 375 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(u0, mg_coord % dir, 'Global direction, PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 375) )
  if (anyExceptions()) return
#line 376 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 376 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(u0, mg_p % dirGLobal(), 'Global direction. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 376) )
  if (anyExceptions()) return
#line 377 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Compare weight
#line 379 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, mg_coord % wgt, 'Weight. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 379) )
  if (anyExceptions()) return
#line 380 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 380 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, mg_p % w, ' Weight. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 380) )
  if (anyExceptions()) return
#line 381 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 381 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(w0, mg_p % w0, 'Starting weight. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 381) )
  if (anyExceptions()) return
#line 382 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Compare time
#line 384 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(t0, mg_coord % time, 'Time. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 384) )
  if (anyExceptions()) return
#line 385 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 385 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(t0, mg_p % time, 'Time. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 385) )
  if (anyExceptions()) return
#line 386 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Compare energy
#line 388 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(G0, mg_coord % G, 'Energy group. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 388) )
  if (anyExceptions()) return
#line 389 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 389 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(E0, ce_coord % E, 'Energy. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 389) )
  if (anyExceptions()) return
#line 390 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 390 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(mg_coord % isMG, 'isMG flag for MG. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 390) )
  if (anyExceptions()) return
#line 391 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 391 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertFalse(ce_coord % isMG, 'isMG flag for CE. PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 391) )
  if (anyExceptions()) return
#line 392 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

#line 393 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(G0, mg_p % G, 'Energy group. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 393) )
  if (anyExceptions()) return
#line 394 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 394 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(E0, ce_p % E, 'Energy. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 394) )
  if (anyExceptions()) return
#line 395 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 395 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(mg_p % isMG, 'isMG flag for MG. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 395) )
  if (anyExceptions()) return
#line 396 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 396 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertFalse(ce_p % isMG, 'isMG flag for CE. Particle from PhaseCoord', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 396) )
  if (anyExceptions()) return
#line 397 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Compare Type
#line 399 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(P_PHOTON, ce_p % type, 'Type of particle has changed', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 399) )
  if (anyExceptions()) return
#line 400 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 400 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(P_NEUTRON, mg_p % type, 'Type of particle has changed', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 400) )
  if (anyExceptions()) return
#line 401 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Verify material, cell IDXs and unique ID
#line 403 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(7, mg_coord % matIdx, 'Material index', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 403) )
  if (anyExceptions()) return
#line 404 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 404 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(8, mg_coord% cellIdx, 'Cell index', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 404) )
  if (anyExceptions()) return
#line 405 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 405 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertEqual(34, mg_coord % uniqueID, 'Unique ID', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 405) )
  if (anyExceptions()) return
#line 406 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  end subroutine testParticleStateAssignments

!  !!
!  !! Test that assignment between particle and particleState is correct
!  !!   DOES NOT RETEST ASSIGNMENT OF phaseCoords parts
!  !!
!@test
!  subroutine testParticleStateAssignment(this)
!    class(test_particle), intent(inout) :: this
!    type(particleState)                 :: state
!
!    state = this % p_MG
!
!    ! Verify
!
!
!  end subroutine testParticleStateAssignment

  !!
  !! Test state saving procedures
  !!
!@test
  subroutine testStateSaving(this)
    class(test_particle), intent(inout) :: this
    type(particleState)                 :: stateRef
    logical(defBool)                    :: isCorrect

    ! **** Do CE particle first*****************************
    ! Save reference
    stateRef = this % p_CE

    ! Check prehistory saving
    call this % p_CE % savePreHistory()
    isCorrect = this % p_CE % preHistory == stateRef
#line 440 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'preHistory check CE', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 440) )
  if (anyExceptions()) return
#line 441 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Check pretransition saving
    call this % p_CE % savePreTransition()
    isCorrect = this % p_CE % preTransition == stateRef
#line 445 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'preTransition check CE', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 445) )
  if (anyExceptions()) return
#line 446 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Check prePath saving
    call this % p_CE % savePrePath()
    isCorrect = this % p_CE % prePath == stateRef
#line 450 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'prePath check CE', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 450) )
  if (anyExceptions()) return
#line 451 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Check precollision saving
    call this % p_CE % savePreCollision()
    isCorrect = this % p_CE % preCollision == stateRef
#line 455 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'preCollision check CE', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 455) )
  if (anyExceptions()) return
#line 456 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! **** Do MG particle **********************************
    ! Save reference
    stateRef = this % p_MG

    ! Check prehistory saving
    call this % p_MG % savePreHistory()
    isCorrect = this % p_MG % preHistory == stateRef
#line 464 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'preHistory check MG', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 464) )
  if (anyExceptions()) return
#line 465 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Check pretransition saving
    call this % p_MG % savePreTransition()
    isCorrect = this % p_MG % preTransition == stateRef
#line 469 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'preTransition check MG', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 469) )
  if (anyExceptions()) return
#line 470 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Check prePath saving
    call this % p_MG % savePrePath()
    isCorrect = this % p_MG % prePath == stateRef
#line 474 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'prePath check MG', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 474) )
  if (anyExceptions()) return
#line 475 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

    ! Check precollision saving
    call this % p_MG % savePreCollision()
    isCorrect = this % p_MG % preCollision == stateRef
#line 479 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue(isCorrect,'preCollision check MG', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 479) )
  if (anyExceptions()) return
#line 480 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testStateSaving

  !!
  !! Test type verification
  !!
!@Test
  subroutine testParticleTypeVerification(this)
    class(test_particle), intent(inout) :: this

#line 490 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue( verifyType(P_NEUTRON), 'Particle Neutron', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 490) )
  if (anyExceptions()) return
#line 491 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 491 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertTrue( verifyType(P_PHOTON), 'Particle Photon', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 491) )
  if (anyExceptions()) return
#line 492 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
#line 492 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"
  call assertFalse( verifyType(-876864), 'Invalid particle type parameter', &
 & location=SourceLocation( &
 & 'particle_test.f90', &
 & 492) )
  if (anyExceptions()) return
#line 493 "/home/mhgk4/SCONE/ParticleObjects/Tests/particle_test.f90"

  end subroutine testParticleTypeVerification
end module particle_test

module Wrapparticle_test
   use pFUnit_mod
   use particle_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_particle) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use particle_test
        class (test_particle), intent(inout) :: this
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

end module Wrapparticle_test

function particle_test_suite() result(suite)
   use pFUnit_mod
   use particle_test
   use Wrapparticle_test
   type (TestSuite) :: suite

   suite = newTestSuite('particle_test_suite')

   call suite%addTest(makeCustomTest('correctInitialisation', correctInitialisation))

   call suite%addTest(makeCustomTest('testPositionAccess', testPositionAccess))

   call suite%addTest(makeCustomTest('testDirectionAccess', testDirectionAccess))

   call suite%addTest(makeCustomTest('testMiscAccess', testMiscAccess))

   call suite%addTest(makeCustomTest('testSetMatIdx', testSetMatIdx))

   call suite%addTest(makeCustomTest('testMovementProcedures', testMovementProcedures))

   call suite%addTest(makeCustomTest('testRotationProcedures', testRotationProcedures))

   call suite%addTest(makeCustomTest('testParticleStateAssignments', testParticleStateAssignments))

   call suite%addTest(makeCustomTest('testStateSaving', testStateSaving))

   call suite%addTest(makeCustomTest('testParticleTypeVerification', testParticleTypeVerification))


end function particle_test_suite

