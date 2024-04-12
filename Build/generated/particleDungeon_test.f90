module particleDungeon_test
  use numPrecision
  use RNG_class,             only : RNG
  use particle_class,        only : particle, particleState
  use particleDungeon_class, only : particleDungeon
  use pFUnit_mod

  implicit none


contains

  !!
  !! Test stack like access. Test is a dummy use case
  !!
!@Test
  subroutine testStackInterface()
    type(particle)        :: p
    type(particleState)   :: phase
    type(particleDungeon) :: dungeon
    integer(shortInt)     :: i

    ! Is empty result for uninitialised dungeon
#line 24 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertTrue(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 24) )
  if (anyExceptions()) return
#line 25 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Initialise
    call dungeon % init(10)

    ! Is empty initialised. No particles stored
#line 30 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertTrue(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 30) )
  if (anyExceptions()) return
#line 31 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 31 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(0, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 31) )
  if (anyExceptions()) return
#line 32 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Store some particles with energy
    p % isMG = .false.
    p % E  = 8.6_defReal

    do i = 1,5
      call dungeon % detain(p)
    end do

    ! Store some phase coordinates with energy
    phase % isMG = .false.
    phase % E = 3.4_defReal

    do i = 1,4
      call dungeon % detain(phase)
    end do

    ! Verify size
#line 50 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertFalse(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 50) )
  if (anyExceptions()) return
#line 51 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 51 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(9, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 51) )
  if (anyExceptions()) return
#line 52 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Remove particles
    do i = 1,4
      call dungeon % release(p)
#line 56 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(3.4_defReal, p % E, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 56) )
  if (anyExceptions()) return
#line 57 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
    end do

     do i = 1,5
      call dungeon % release(p)
#line 61 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(8.6_defReal, p % E, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 61) )
  if (anyExceptions()) return
#line 62 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
    end do

    ! Verify size
#line 65 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertTrue(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 65) )
  if (anyExceptions()) return
#line 66 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 66 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(0, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 66) )
  if (anyExceptions()) return
#line 67 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Clean
    call dungeon % kill()
  end subroutine testStackInterface

  !!
  !! Test array like interface
  !!
!@Test
  subroutine testArrayInterface
    type(particle)        :: p
    type(particleState)   :: phase
    type(particleDungeon) :: dungeon
    integer(shortInt)     :: i

    ! Is empty result for uninitialised dungeon
#line 83 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertTrue(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Initialise to fixed size
    call dungeon % setSize(2)

    ! Is filled with random particles
#line 89 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertFalse(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 89) )
  if (anyExceptions()) return
#line 90 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 90 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(2, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 90) )
  if (anyExceptions()) return
#line 91 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Extend size
    call dungeon % setSize(9)
#line 94 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertFalse(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 94) )
  if (anyExceptions()) return
#line 95 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 95 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(9, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 95) )
  if (anyExceptions()) return
#line 96 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Store some particles with energy
    p % isMG = .false.
    p % E  = 8.6_defReal

    do i = 1,5
      call dungeon % replace(p, i)
    end do

    ! Store some phase coordinates with energy
    phase % isMG = .false.
    phase % E = 3.4_defReal

    do i = 1,4
      call dungeon % replace(phase, 5 + i)
    end do

    ! Verify size
#line 114 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertFalse(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 114) )
  if (anyExceptions()) return
#line 115 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 115 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(9, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 115) )
  if (anyExceptions()) return
#line 116 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Raplace by particle and phaseCoords
    p % E = 13.0_defReal
    phase % E = 17.0_defReal
    call dungeon % replace(p, 9)
    call dungeon % replace(phase, 1)

    ! Verify particles by copies
    call dungeon % copy(p, 9)
#line 125 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(13.0_defReal, p % E, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 125) )
  if (anyExceptions()) return
#line 126 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    do i=8,6,-1
      call dungeon % copy(p, i)
#line 129 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(3.4_defReal, p % E, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 129) )
  if (anyExceptions()) return
#line 130 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
    end do

    do i=5,2,-1
      call dungeon % copy(p, i)
#line 134 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(8.6_defReal, p % E, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 134) )
  if (anyExceptions()) return
#line 135 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
    end do
    call dungeon % copy(p, 1)
#line 137 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(17.0_defReal, p % E, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 137) )
  if (anyExceptions()) return
#line 138 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Verify that population has not changed
#line 140 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertFalse(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 140) )
  if (anyExceptions()) return
#line 141 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 141 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(9, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 141) )
  if (anyExceptions()) return
#line 142 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Shrink size
    call dungeon % setSize(2)
#line 145 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertFalse(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 145) )
  if (anyExceptions()) return
#line 146 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 146 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(2, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 146) )
  if (anyExceptions()) return
#line 147 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Clean population
    call dungeon % cleanPop()
#line 150 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertTrue(dungeon % isEmpty(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 150) )
  if (anyExceptions()) return
#line 151 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
#line 151 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(0, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 151) )
  if (anyExceptions()) return
#line 152 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Clean memory
    call dungeon % kill()

  end subroutine testArrayInterface

  !!
  !! Test weight normalisation and inquiry
  !!
!@Test
  subroutine testWeightNorm()
    type(particle)        :: p
    type(particleDungeon) :: dungeon
    integer(shortInt)     :: i
    real(defReal), parameter :: TOL = 1.0E-9

    ! Initialise
    call dungeon % init(10)
#line 170 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(ZERO, dungeon % popWeight(), TOL, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 170) )
  if (anyExceptions()) return
#line 171 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Store some particles with non-uniform weight
    do i = 1,5
      p % w = 0.5_defReal + i * 0.1_defReal
      call dungeon % detain(p)

    end do

    ! Verify total weight
#line 180 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(4.0_defReal, dungeon % popWeight(), TOL, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 180) )
  if (anyExceptions()) return
#line 181 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Normalise weight
    call dungeon % normWeight(2.0_defReal)
#line 184 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(2.0_defReal, dungeon % popWeight(), TOL, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 184) )
  if (anyExceptions()) return
#line 185 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Get particles and compare weight
    do i = 5,1,-1
      call dungeon % release(p)
#line 189 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual( 0.25_defReal + i * 0.05_defReal, p % w, TOL, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 189) )
  if (anyExceptions()) return
#line 190 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
    end do

    ! Verify weight
#line 193 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(0.0_defReal, dungeon % popWeight(), TOL, &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 193) )
  if (anyExceptions()) return
#line 194 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Clean
    call dungeon % kill()
  end subroutine testWeightNorm

  !!
  !! Test normalisation of population to smaller number
  !! Particles with non-uniform weight
  !!  NOTE: Weight preservation is disabled for now
  !!
!@Test
  subroutine testNormPopDown()
    type(particleDungeon)    :: dungeon
    type(particle)           :: p
    type(RNG)                :: pRNG
    integer(shortInt)        :: i
    real(defReal), parameter :: TOL = 1.0E-9


    ! Initialise
    call dungeon % init(10)
    call pRNG % init(7865856_longInt)

    ! Store some particles with non-uniform weight
    do i = 1,10
      p % w = 0.5_defReal + i * 0.1_defReal
      call dungeon % detain(p)
    end do

    ! Normalise population
    call dungeon % normSize(5, pRNG)

    ! Verify size
#line 227 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(5, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 227) )
  if (anyExceptions()) return
#line 228 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Verify weight *** DISABLED
    !@assertEqual(6.05_defReal, dungeon % popWeight(), TOL)

    ! Clean memory
    call dungeon % kill()

  end subroutine testNormPopDown

  !!
  !! Test normalisation of population to smaller number
  !! Particles with non-uniform weight
  !!  NOTE: Weight preservation is disabled for now
  !!
!@Test
  subroutine testNormPopUp()
    type(particleDungeon)    :: dungeon
    type(particle)           :: p
    type(RNG)                :: pRNG
    integer(shortInt)        :: i
    real(defReal), parameter :: TOL = 1.0E-9

    ! Initialise
    call dungeon % init(20)
    call pRNG % init(435468_longInt)

    ! Store some particles with non-uniform weight
    do i = 1,10
      p % w = 0.5_defReal + i * 0.1_defReal
      call dungeon % detain(p)
    end do

    ! Normalise population
    call dungeon % normSize(15, pRNG)

    ! Verify size
#line 264 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"
  call assertEqual(15, dungeon % popSize(), &
 & location=SourceLocation( &
 & 'particleDungeon_test.f90', &
 & 264) )
  if (anyExceptions()) return
#line 265 "/home/mhgk4/SCONE/ParticleObjects/Tests/particleDungeon_test.f90"

    ! Verify weight *** DISABLED
    !@assertEqual(18.15_defReal, dungeon % popWeight(), TOL)

    ! Clean memory
    call dungeon % kill()

  end subroutine testNormPopUp


end module particleDungeon_test

module WrapparticleDungeon_test
   use pFUnit_mod
   use particleDungeon_test
   implicit none
   private

contains


end module WrapparticleDungeon_test

function particleDungeon_test_suite() result(suite)
   use pFUnit_mod
   use particleDungeon_test
   use WrapparticleDungeon_test
   type (TestSuite) :: suite

   suite = newTestSuite('particleDungeon_test_suite')

   call suite%addTest(newTestMethod('testStackInterface', testStackInterface))

   call suite%addTest(newTestMethod('testArrayInterface', testArrayInterface))

   call suite%addTest(newTestMethod('testWeightNorm', testWeightNorm))

   call suite%addTest(newTestMethod('testNormPopDown', testNormPopDown))

   call suite%addTest(newTestMethod('testNormPopUp', testNormPopUp))


end function particleDungeon_test_suite

