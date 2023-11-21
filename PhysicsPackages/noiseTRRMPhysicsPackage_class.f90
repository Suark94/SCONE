module noiseTRRMPhysicsPackage_class

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential, exponentialComplex
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile
  use rng_class,                      only : RNG
  use physicsPackage_inter,           only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryStd_class,              only : geometryStd
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx, gr_kill    => kill

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat            => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init         => init, &
                                             ndReg_getMatNames  => getMatNames, &
                                             ndReg_activate     => activate, &
                                             ndReg_kill         => kill, &
                                             ndReg_getNeutronMG => getNeutronMG
  use materialHandle_inter,           only : materialHandle
  use mgNeutronDatabase_inter,        only : mgNeutronDatabase
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  use fissionMG_class,                only : fissionMG, fissionMG_TptrCast

  ! Visualisation
  use visualiser_class,               only : visualiser

  ! Tally map for fission rate
  use tallyMap_inter,                 only : tallyMap
  use tallyMapFactory_func,           only : new_tallyMap

  ! Random ray - or a standard particle
  ! Also particleState for easier output
  use particle_class,                 only : ray => particle, particleState

  ! For locks
  use omp_lib

  implicit none
  private

  ! Parameter for when to skip a tiny volume
  real(defReal), parameter :: volume_tolerance = 1.0E-12

  complex(defFlt), parameter :: imag = cmplx(0.0_defFlt, 1.0_defFlt, defFlt)

  ! Parameter for entropy vector length: must be divisible by 2
  integer(shortInt), parameter :: entLength = 400

  !!
  !! Physics package to perform The Random Ray Method (TRRM) neutron noise calculations
  !!
  !! Solves both the eigenvalue and noise equations simultaneously such that one need 
  !! not store the angular flux for use in the noise solve. This will result in a slower
  !! solution (as opposed to running a noise solver by itself) but should give a large
  !! reduction in memory.
  !!
  !! Necessitates the use of complex variables for the scalar flux noise vector.
  !!
  !! The solve is performed for a fixed noise frequency
  !!
  !! Tracks rays across the geometry, attenuating their flux. After some dead length,
  !! rays begin scoring to estimates of the fluxes and volume. Each ray has a
  !! uniform termination length, after which it is stopped and the next ray is tracked.
  !! Once all rays have been tracked, a cycle concludes and fluxes, sources, and keff
  !! are updated.
  !!
  !! Both inactive and active cycles occur, as in Monte Carlo. These can be terminated
  !! after a specified number of iterations or on reaching some chosen convergence
  !! criterion (Shannon entropy-based for regular flux, phase-based for noise).
  !!
  !! Calculates relative volume of different materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised such that the total domain
  !! volume is 1.0.
  !!
  !! IMPORTANT N.B.: Geometry type must be extended! Won't run if shrunk.
  !! This is because spatial discretisation is determined by the number of unique cells in the
  !! geometry.
  !! Also, this is obviously for multi-group calculations only.
  !!
  !! Must provide a neutron noise source.
  !!
  !! Sample Input Dictionary:
  !!   PP {
  !!     type noiseTRRMPhysicsPackage;
  !!     dead 10;              // Dead length where rays do not score to scalar fluxes
  !!     termination 100;      // Length a ray travels before it is terminated
  !!     rays 1000;            // Number of rays to sample per iteration
  !!     inactive 100;         // Number of convergence cycles (would use accum and std otherwise)
  !!     active 200;           // Number of scoring cycles (would use eps otherwise)
  !!     #seed 86868;#         // Optional RNG seed
  !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
  !!     #fissionMap {<map>}#  // Optionally output fission rates according to a given map
  !!     #fluxMap {<map>}#     // Optionally output one-group fluxes according to a given map
  !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
  !!
  !!     noise {
  !!       omega 0.01;           // Noise frequency in Hz
  !!       mat   <name>;         // Name of noisy material
  !!       .....NEED TO MAKE LOTS OF DECISIONS HERE
  !!     }
  !!
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition>}
  !!   }
  !!
  !! Private Members
  !!   geom        -> Pointer to the geometry.
  !!   geomIdx     -> Index of the geometry in geometry Registry.
  !!   top         -> Top co-ordinates of the geometry bounding box.
  !!   bottom      -> Bottom co-ordinates of the geometry bounding box.
  !!   rand        -> Random number generator.
  !!   timerMain   -> Index of the timer defined to measure calculation time.
  !!   mgData      -> MG database. Calculation obviously cannot be run in CE.
  !!   nG          -> Number of energy groups, kept for convenience.
  !!   nCells      -> Number of unique cells in the geometry, kept for convenience.
  !!   nMat        -> Number of unique materials in the geometry - for convenience.
  !!   lengthPerIt -> Distance all rays travel in a single iteration - for convenience.
  !!
  !!   omega       -> Frequency of the noise source
  !!   nP          -> Number of precursors
  !!
  !!   termination -> Distance a ray can travel before it is terminated
  !!   dead        -> Distance a ray must travel before it becomes active
  !!   pop         -> Number of rays to track per cycle
  !!   inactive    -> Number of inactive cycles to perform
  !!   active      -> Number of active cycles to perform
  !!   cache       -> Logical check whether to use distance caching
  !!   outputFile  -> Output file name
  !!   outputFormat-> Output file format
  !!   plotResults -> Plot results?
  !!   printFluxes -> Print fluxes?
  !!   printVolume -> Print volumes?
  !!   printCells  -> Print cell positions?
  !!   viz         -> Output visualiser
  !!   mapFission  -> Output fission rates across a given map?
  !!   resultsMap  -> The map across which to output fission rate results
  !!   mapFission  -> Output 1G flux across a given map?
  !!   resultsMap  -> The map across which to output 1G flux results
  !!
  !!   sigmaT      -> Local total cross section vector
  !!   nuSigmaF    -> Local nuSigmaF vector
  !!   sigmaS      -> Local flattened scattering matrix
  !!   chi         -> Local chi effective vector
  !!   chiP        -> Local prompt chi vector
  !!   chiD        -> Local delayed chi vector
  !!   lambda      -> Local precursor decay constant vector
  !!   speed       -> Local neutron speed vector
  !!   dSigmaT     -> Local perturbed SigmaT vector (in freq domain)
  !!   dNuSigmaF   -> Local perturbed NuSigmaF vector (in freq domain)
  !!   dSigmaS     -> Local perturbed SigmaS vector (in freq domain)
  !!
  !!   keff        -> Estimated value of keff
  !!   keffScore   -> Vector holding cumulative keff score and keff^2 score
  !!   scalarFlux  -> Array of scalar flux values of length = nG * nCells
  !!   prevFlux    -> Array of previous scalar flux values of length = nG * nCells
  !!   fluxScores  -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!   source      -> Array of neutron source values of length = nG * nCells
  !!
  !!   scalarNoise -> Array of scalar neutron noise values of length = nG * nCells
  !!   prevNoise   -> Array of previous scalar neutron noise values of length = nG * nCells
  !!   noiseScores -> Array of scalar neutron noise values and squared values to be reported
  !!                  in results, dimension = [nG * nCells, 2]
  !!   sourceNoise -> Array of neutron noise source values of length = nG * nCells
  !!   fixedSourceNoise -> Array of fixed neutron noise source values of length = nG * nCells
  !!
  !!   volume      -> Array of stochastically estimated cell volumes of length = nCells
  !!   cellHit     -> Array tracking whether given cells have been hit during tracking
  !!   cellFound   -> Array tracking whether a cell was ever found
  !!   cellPos     -> Array of cell positions, populated once they are found
  !!
  !!   entropyFluxVec  -> vector containing entropy values over consecutive iterations
  !!   entropyNoiseVec -> vector containing 'entropy' noise values over consecutive iterations
  !!
  !!   entropy         -> entropy values in individual cells corresponding to entropyMap
  !!   entropyNoise    -> 'entropy' noise values in individual cells corresponding to entropyMap
  !!   entropyMap      -> map over which to calculate entropies
  !!   entMapSize      -> size of the entropy map  
  !!
  !!   locks       -> OpenMP locks used when writing to scalar flux during transport
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: noiseTRRMPhysicsPackage
    private
    ! Components
    class(geometryStd), pointer           :: geom
    integer(shortInt)                     :: geomIdx     = 0
    real(defReal), dimension(3)           :: top         = ZERO
    real(defReal), dimension(3)           :: bottom      = ZERO
    type(RNG)                             :: rand
    class(baseMgNeutronDatabase), pointer :: mgData      => null()
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    integer(shortInt)                     :: nMat        = 0
    real(defReal)                         :: lengthPerIt = ZERO
    
    integer(shortInt)                     :: nP = 0
    real(defFlt)                          :: omega = 0.0_defFlt

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .false.
    logical(defBool)   :: doEntropy   = .false.
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    logical(defBool)   :: plotResults = .false.
    logical(defBool)   :: printFlux   = .false.
    logical(defBool)   :: printVolume = .false.
    logical(defBool)   :: printCells  = .false.
    type(visualiser)   :: viz
    logical(defBool)   :: mapFission  = .false.
    class(tallyMap), allocatable :: resultsMap
    logical(defBool)   :: mapFlux     = .false.
    class(tallyMap), allocatable :: fluxMap

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:), allocatable     :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi
    real(defFlt), dimension(:), allocatable     :: chiP
    real(defFlt), dimension(:), allocatable     :: chiD
    real(defFlt), dimension(:), allocatable     :: lambda
    real(defFlt), dimension(:), allocatable     :: beta
    real(defFlt), dimension(:), allocatable     :: speed
    
    complex(defFlt), dimension(:), allocatable  :: dSigmaT
    complex(defFlt), dimension(:), allocatable  :: dNuSigmaF
    complex(defFlt), dimension(:), allocatable  :: dSigmaS

    ! Results space
    real(defReal)                                :: keff
    real(defReal), dimension(2)                  :: keffScore
    real(defFlt), dimension(:), allocatable      :: scalarFlux
    real(defFlt), dimension(:), allocatable      :: prevFlux
    real(defReal), dimension(:,:), allocatable   :: fluxScores
    real(defFlt), dimension(:), allocatable      :: source
    real(defReal), dimension(:), allocatable     :: volume
    real(defReal), dimension(:), allocatable     :: volumeTracks
   
    ! Noise results
    complex(defFlt), dimension(:), allocatable    :: scalarNoise
    complex(defFlt), dimension(:), allocatable    :: prevNoise
    complex(defReal), dimension(:,:), allocatable :: noiseScores
    complex(defFlt), dimension(:), allocatable    :: sourceNoise
    complex(defFlt), dimension(:), allocatable    :: fixedSourceNoise

    ! Tracking cell properites
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos
    integer(longInt)                             :: intersectionsTotal = 0

    ! Entropies
    real(defReal), dimension(entLength) :: entropyFluxVec  = ZERO
    real(defReal), dimension(entLength) :: entropyNoiseVec = ZERO

    real(defReal), dimension(:,:), allocatable :: entropy
    real(defReal), dimension(:,:), allocatable :: entropyNoise
    class(tallyMap), allocatable :: entropyMap
    integer(shortInt) :: entMapSize = -1


    ! OMP locks
    integer(kind=omp_lock_kind), dimension(:), allocatable :: locks

    ! Timer bins
    integer(shortInt) :: timerMain
    integer(shortInt) :: timerTransport
    real (defReal)    :: time_transport = ZERO
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: initialiseNoise
    procedure, private :: cycles
    procedure, private :: initialiseRay
    procedure, private :: transportSweep
    procedure, private :: sourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: calculateEntropy
    procedure, private :: calculateEntropyMap
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings

  end type noiseTRRMPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
    integer(shortInt)                             :: seed_temp, i, g, g1, m, nP, nPNew, p, idx0
    integer(longInt)                              :: seed
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    class(dictionary),pointer                     :: tempDict, graphDict, noiseDict
    class(mgNeutronDatabase),pointer              :: db
    character(nameLen)                            :: geomName, graphType, nucData
    class(geometry), pointer                      :: geom
    type(outputFile)                              :: test_out
    class(baseMgNeutronMaterial), pointer         :: mat
    class(materialHandle), pointer                :: matPtr
    class(fissionMG), pointer                     :: fiss
    character(100), parameter :: Here = 'init (noiseTRRMPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)
    
    ! Load settings
    call dict % get( nucData, 'XSdata')
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    call dict % get(self % active, 'active')
    
    ! Use entropy to measure convergence?
    call dict % getOrDefault(self % doEntropy, 'entropy', .false.)
    if (.not. self % doEntropy) call dict % get(self % inactive, 'inactive')

    ! Perform distance caching?
    call dict % getOrDefault(self % cache, 'cache', .false.)

    ! Print fluxes?
    call dict % getOrDefault(self % printFlux, 'printFlux', .false.)

    ! Print volumes?
    call dict % getOrDefault(self % printVolume, 'printVolume', .false.)

    ! Print cell positions?
    call dict % getOrDefault(self % printCells, 'printCells', .false.)

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Check settings
    if (self % termination <= ZERO) call fatalError(Here, 'Ray termination distance (termination) is less than or equal to zero.')
    if (self % pop < 1) call fatalError(Here, 'Must have 1 or more rays (pop).')

    ! Dead length can be less than zero but will be reset to zero if so
    if (self % dead < ZERO) then
      self % dead = ZERO
      print *,'Warning: Dead length of rays (dead) was negative. This has been set to 0 instead.'
    end if

    ! Ensure termination length is longer than dead length
    if (self % termination <= self % dead) call fatalError(Here,&
            'Ray termination length must be greater than ray dead length')

    ! Check whether there is a map for outputting fission rates
    ! If so, read and initialise the map to be used
    if (dict % isPresent('fissionMap')) then
      self % mapFission = .true.
      tempDict => dict % getDictPtr('fissionMap')
      call new_tallyMap(self % resultsMap, tempDict)
    else
      self % mapFission = .false.
    end if
    
    ! Check whether there is a map for outputting one-group fluxes
    ! If so, read and initialise the map to be used
    if (dict % isPresent('fluxMap')) then
      self % mapFlux = .true.
      tempDict => dict % getDictPtr('fluxMap')
      call new_tallyMap(self % fluxMap, tempDict)
    else
      self % mapFlux = .false.
    end if
    
    ! Check whether there is a map for outputting one-group fluxes
    ! If so, read and initialise the map to be used
    if (dict % isPresent('entropyMap')) then
      self % doEntropy = .true.
      tempDict => dict % getDictPtr('entropyMap')
      call new_tallyMap(self % entropyMap, tempDict)
      self % entMapSize = product(self % entropyMap % binArrayShape())
      allocate(self % entropy(entLength, self % entMapSize))
      allocate(self % entropyNoise(entLength, self % entMapSize))
      self % entropy = ZERO
      self % entropyNoise = ZERO
    end if


    ! Register timer
    self % timerMain = registerTimer('simulationTime')
    self % timerTransport = registerTimer('transportTime')

    ! Initialise RNG
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')
    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % rand % init(seed)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'randomRayGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    geom    => gr_geomPtr(self % geomIdx)

    ! Ensure geometry is geometryStd
    select type(geom)
      type is (geometryStd)
        self % geom => geom
      class default
        call fatalError(Here,'Unrecognised geometry type')
    end select

    ! Ensure that geometry graph is extended
    graphDict => tempDict % getDictPtr('graph')
    call graphDict % get(graphType,'type')
    if (graphType /= 'extended') call fatalError(Here,&
            'Geometry graph type must be "extended" for random ray calculations.')

    ! Activatee nuclear data
    call ndReg_activate(P_NEUTRON_MG, nucData, self % geom % activeMats())

    ! Ensure that nuclear data is multi-group
    db => ndReg_getNeutronMG()
    if (.not. associated(db)) call fatalError(Here,&
            'No MG nuclear database was constructed')

    ! Ensure nuclear data is baseMgNeutronDatabase
    select type(db)
      type is (baseMgNeutronDatabase)
        self % mgData => db
      class default
        call fatalError(Here,'Unrecognised MG database type')
    end select

    ! Store number of energy groups for convenience
    self % nG = self % mgData % nGroups()

    ! Get lower and upper corner of bounding box
    associate (aabb => self % geom % bounds())
      self % bottom = aabb(1:3)
      self % top    = aabb(4:6)
    end associate
    
    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      print *, "Constructing visualisation"
      call self % viz % makeViz()
      call self % viz % kill()
    endif
    
    ! Check for results plotting and initialise VTK
    call dict % getOrDefault(self % plotResults,'plot',.false.)
    if (self % plotResults) then
      ! Initialise a visualiser to be used when results are available
      print *, "Initialising results visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      print *, "Constructing geometry visualisation"
      call self % viz % initVTK()
    end if

    ! Store number of cells in geometry for convenience
    self % nCells = self % geom % numberOfCells()

    ! Allocate results space
    allocate(self % scalarFlux(self % nCells * self % nG))
    allocate(self % prevFlux(self % nCells * self % nG))
    allocate(self % fluxScores(self % nCells * self % nG, 2))
    allocate(self % source(self % nCells * self % nG))
    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))
    allocate(self % scalarNoise(self % nCells * self % nG))
    allocate(self % prevNoise(self % nCells * self % nG))
    allocate(self % noiseScores(self % nCells * self % nG, 2))
    allocate(self % sourceNoise(self % nCells * self % nG))
    allocate(self % fixedSourceNoise(self % nCells * self % nG))
    
    ! Set active length traveled per iteration
    self % lengthPerIt = (self % termination - self % dead) * self % pop
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
      call omp_init_lock(self % locks(i))
    end do

    ! Initialise local nuclear data
    ! TODO: clean nuclear database afterwards! It is no longer used
    !       and takes up memory.
    print *,'Setting up XS arrays'
    self % nMat = mm_nMat()
    allocate(self % sigmaT(self % nMat * self % nG))
    allocate(self % nuSigmaF(self % nMat * self % nG))
    allocate(self % chi(self % nMat * self % nG))
    allocate(self % speed(self % nMat * self % nG))
    allocate(self % sigmaS(self % nMat * self % nG * self % nG))

    ! Keep track of number of precursor groups
    ! Error if not uniform across all materials
    nP = -1
    fiss => null()
    do m = 1, self % nMat
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      do g = 1, self % nG
        self % sigmaT(self % nG * (m - 1) + g) = real(mat % getTotalXS(g, self % rand),defFlt)
        self % nuSigmaF(self % nG * (m - 1) + g) = real(mat % getNuFissionXS(g, self % rand),defFlt)
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, self % rand),defFlt)
        self % speed(self % nG * (m - 1) + g) = real(mat % getSpeed(g, self % rand),defFlt)
        do g1 = 1, self % nG
          self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)  = &
                  real(mat % getScatterXS(g1, g, self % rand), defFlt)
        end do
      end do
      ! This is going to be very fudgy...
      ! Determine number of precusors
      fiss => fissionMG_TptrCast(self % mgData % getReaction(macroFission,m))
      if (associated(fiss)) then
        if (allocated(fiss % delayData)) then
          nPNew = size(fiss % delayData,1)
          if (nP < 0) then
            nP = nPNew
          elseif (nP /= nPNew) then
            call fatalError(Here,'Different numbers of precusor groups for different materials: '//&
                    numToChar(nP)//' and '//numToChar(nPNew))
          end if
        end if
      end if
      fiss => null()
    end do

    ! Also fudgy
    ! Prepare delayed neutron data
    ! Assumes chi in the input is the cumulative chi - back calculate chi prompt
    if (nP /= -1) then
      self % nP = nP
      allocate(self % chiP(self % nMat * self % nG))
      allocate(self % chiD(self % nMat * self % nG * nP))
      allocate(self % beta(self % nMat * nP))
      allocate(self % lambda(self % nMat * nP))
      self % chiD = 0.0_defFlt
      self % chiP = 0.0_defFlt
      !self % chiP = self % chi
      do m = 1, self % nMat
        fiss => fissionMG_TptrCast(self % mgData % getReaction(macroFission,m))
        if (associated(fiss)) then
          if (allocated(fiss % delayData)) then
            do p = 1, self % nP
              self % lambda(self % nP * (m - 1) + p) = real(fiss % delayData(p,1),defFlt)
              self % beta(self % nP * (m - 1) + p) = real(fiss % delayData(p,2),defFlt)
              do g = 1, self % nG
                self % chiD(self % nP * self % nG * (m - 1) + self % nP * (g - 1) + p) = &
                        real(fiss % chiDelayed(p,g),defFlt)
              end do
            end do
          else
            idx0 = self % nP * self % nG * (m - 1)
            do i = 1, self % nP * self % nG
              self % chiD(idx0 + i) = 0.0_defFlt
            end do
            idx0 = self % nP * (m - 1)
            do i = 1, self % nP
              self % beta(idx0 + i) = 0.0_defFlt
              self % lambda(idx0 + i) = 0.0_defFlt
            end do
          end if
        else
          idx0 = self % nP * self % nG * (m - 1)
          do i = 1, self % nP * self % nG
            self % chiD(idx0 + i) = 0.0_defFlt
          end do
          idx0 = self % nP * (m - 1)
          do i = 1, self % nP
            self % beta(idx0 + i) = 0.0_defFlt
            self % lambda(idx0 + i) = 0.0_defFlt
          end do
        end if
        fiss => null()

        ! Construct chi_prompt
        idx0 = self % nG * (m - 1)
        do g = 1, self % nG
          self % chiP(idx0 + g) = self % chi(idx0 + g) 
        end do
        do g = 1, self % nG
          do p = 1, self % nP
            self % chiP(self % nG * (m - 1) + g) = self % chiP(self % nG * (m - 1) + g) &
                    - self % beta(self % nP * (m - 1) + p) * &
                    self % chiD(self % nP * self % nG * (m - 1) + self % nP * (g - 1) + p)
          end do
          self % chiP(self % nG * (m - 1) + g) = self % chiP(self % nG * (m - 1) + g) / &
                  (1.0_defFlt - sum(self % beta((self % nP * (m - 1) + 1):(self % nP * (m - 1) + self % nP))))
        end do
      end do
    else
      self % nP = 1
      allocate(self % chiP(self % nMat * self % nG))
      allocate(self % chiD(self % nMat * self % nG * self % nP))
      allocate(self % beta(self % nMat * self % nP))
      allocate(self % lambda(self % nMat * self % nP))
      self % chiP   = self % chi
      self % chiD   = 0.0_defFlt
      self % beta   = 0.0_defFlt
      self % lambda = 0.0_defFlt

    end if

    ! Initialise noise source
    print *,'Initialising noise source'
    noiseDict => dict % getDictPtr('noise')
    call self % initialiseNoise(noiseDict)

  end subroutine init

  !!
  !! Initialises the noise source to be used in the simulation
  !! Takes a dictionary containing noise frequency, names of materials in the geometry,
  !! cross sections which vary, and their relative amplitudes in each group
  !! and places these in the appropriate elements of the dXS vectors
  !! Assumes everything is an absorber of variable strength
  !!
  !! Example dictionary:
  !!
  !! noise {
  !!   omega 0.1; // omega = 2*pi*freq
  !!
  !!   materials {
  !!
  !!    #matname1 {#
  !!
  !!       fixedReal (1.0 23.7 0.0 ...);
  !!       fixedImag (123123 0.7 3.2 ...);
  !!
  !!      reactions{
  !!         capture (0.05 0.01 0.02 ... ); // vector of nG elements describing group-wise capture variation
  !!         scatter (0.01 0.01 0.01 ... ); // vector of nG elements describing group-wise scatter variation
  !!         fission (0.1  0.005 0.2 ... ); // vector of nG elements describing group-wise fission variation
  !!         total   (0.1  0.005 0.2 ... ); // vector of nG elements describing group-wise total variation
  !!       }
  !!    # } #
  !!
  !!   #matname2 {#
  !!   # ... #
  !!   # } #
  !!
  !!   }
  !!
  !! }
  !!
  subroutine initialiseNoise(self, dict)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
    class(dictionary),pointer                     :: matDict, myMatDict, reacDict
    character(nameLen),dimension(:), allocatable  :: names, reactions
    real(defReal), dimension(:), allocatable      :: temp, realPart, imagPart
    real(defReal)                                 :: omega
    integer(shortInt)                             :: g, g1, m, i, r, cIdx
    integer(shortInt), save                       :: idx, gL, matIdx
    logical(defBool)                              :: found, foundIdx, totalPresent, isSine
    character(nameLen)                            :: localName
    character(nameLen), save                      :: local
    class(baseMgNeutronMaterial), pointer         :: mat
    class(materialHandle), pointer                :: matPtr
    character(100), parameter :: Here = 'initialiseNoise (noiseTRRMPhysicsPackage_class.f90)'
    !$omp threadprivate(matIdx, local, idx, gL)

    ! Allocate dXSs
    allocate(self % dSigmaT(self % nMat * self % nG))
    allocate(self % dNuSigmaF(self % nMat * self % nG))
    allocate(self % dSigmaS(self % nMat * self % nG * self % nG))
    self % dSigmaT   = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    self % dNuSigmaF = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    self % dSigmaS   = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    self % fixedSourceNoise = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    
    ! Read frequency for solve
    call dict % get(omega, 'omega')
    self % omega = real(omega, defFlt)

    matDict => dict % getDictPtr('materials')
    call matDict % keys(names)

    ! Cycle through materials in the dictionary
    do i = 1, size(names)

      ! Check this material exists in the problem and find its index, m
      found = .false.
      do m = 1, self % nMat
        localName = mm_matName(m)
        if (localName == names(i)) found = .true.
        if (found) exit
      end do
      
      if (.not. found) call fatalError(Here,'The material '//trim(names(i))//' does not correspond to '//&
              'any material found in the geometry.')
      
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)

      myMatDict => matDict % getDictPtr(names(i))
      
      ! Read fixed noise sources
      if (myMatDict % isPresent('fixedReal') .OR. myMatDict % isPresent('fixedImag')) then
        
        if (myMatDict % isPresent('fixedReal')) then
          call myMatDict % get(realPart,'fixedReal')
        else
          allocate(realPart(self % nG))
          realPart = 0.0_defFlt
        end if
        if (myMatDict % isPresent('fixedImag')) then
          call myMatDict % get(imagPart,'fixedImag')
        else
          allocate(imagPart(self % nG))
          imagPart = 0.0_defFlt
        end if
        
        ! Ensure correct number of energy groups
        if (size(realPart) /= self % nG) call fatalError(Here,'Real source '//trim(names(i))//&
                ' has '//numToChar(size(realPart))//' groups rather than '//numToChar(self % nG))
        if (size(imagPart) /= self % nG) call fatalError(Here,'Imaginary source '//trim(names(i))//&
                ' has '//numToChar(size(imagPart))//' groups rather than '//numToChar(self % nG))
      
        ! Make sure that the source corresponds to a material present in the geometry
        foundIdx = .false.
        !$omp parallel do schedule(static) 
        do cIdx = 1, self % nCells

          matIdx    = self % geom % geom % graph % getMatFromUID(cIdx)
          local = mm_matName(matIdx)

          if (local == names(i)) then

            foundIdx = .true.
            do gL = 1, self % nG
            
              idx = (cIdx - 1) * self % nG + gL
              self % fixedSourceNoise(idx)%re = real(realPart(gL),defFlt)
              self % fixedSourceNoise(idx)%im = real(imagPart(gL),defFlt)

            end do

          end if

        end do
        !$omp end parallel do

        if (.not. found) call fatalError(Here,'The source '//trim(names(i))//' does not correspond to '//&
                'any material found in the geometry.')

      end if

      ! Read XS noise sources
      if (myMatDict % isPresent('reactions')) then
      
        reacDict => myMatDict % getDictPtr('reactions') 
        call reacDict % keys(reactions)

        if (size(reactions) == 0) call fatalError(Here,'No reactions specified for the oscillating material')
      
        ! Go through all noisy reactions, setting up the dXS terms
        totalPresent = .false.
        do r = 1, size(reactions)
        
          call reacDict % get(temp, reactions(r))

          ! Ensure correct number of energy groups
          if (size(temp) /= self % nG) call fatalError(Here,'Noise in reaction '//reactions(r)//&
                  ' has '//numToChar(size(temp))//' groups rather than '//numToChar(self % nG))
        
          ! Perform unique operations for each permitted reaction type
          select case(reactions(r))
            case('capture')
              do g = 1, self % nG
                if(.not. totalPresent) self % dSigmaT(self % nG * (m - 1) + g)%re = self % dSigmaT(self % nG * (m - 1) + g)%re + &
                      real(temp(g) * mat % getCaptureXS(g, self % rand),defFlt)
              end do

            case('fission')
              do g = 1, self % nG
                self % dNuSigmaF(self % nG * (m - 1) + g)%re = real(temp(g) * mat % getNuFissionXS(g, self % rand),defFlt)
                if (.not. totalPresent) self % dSigmaT(self % nG * (m - 1) + g)%re = self % dSigmaT(self % nG * (m - 1) + g)%re + &
                      real(temp(g) * mat % getFissionXS(g, self % rand),defFlt)
              end do

            case('scatter')
              do g = 1, self % nG
                do g1 = 1, self % nG
                  if(.not. totalPresent) self % dSigmaT(self % nG * (m - 1) + g)%re = self % dSigmaT(self % nG * (m - 1) + g)%re + &
                        real(temp(g) * mat % getScatterXS(g, g1, self % rand),defFlt)
                  self % dSigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)%re = &
                        real(temp(g) * mat % getScatterXS(g, g1, self % rand),defFlt)
                end do
              end do

            ! Total overwrites any changes in dSigmaT due to other perturbations
            ! Essentially allows easier configuration of dSigmaC perturbations (which might otherwise be the balancing perturbation)
            case('total')
              totalPresent = .true.
              do g = 1, self % nG
                self % dSigmaT(self % nG * (m - 1) + g)%re = real(temp(g) * self % sigmaT(self % nG * (m - 1) + g),defFlt)
              end do

            case default
              call fatalError(Here,"Unrecognised reaction type. Must be capture, fission, scatter, or total")
          end select

        end do

      end if

    end do
    
    ! Apply scaling to make all noise perturbations due to 'variable strength absorbers'
    ! -i corresponds to sine, 1 corresponds to cosine
    ! For now, uniform across all materials
    call dict % get(isSine,'sine')
    if (isSine) then
      self % dSigmaT%im   = -sign(1.0_defFlt,self % omega) * real(PI * self % dSigmaT%re,defFlt)
      self % dNuSigmaF%im = -sign(1.0_defFlt,self % omega) * real(PI * self % dNuSigmaF%re, defFlt)
      self % dSigmaS%im   = -sign(1.0_defFlt,self % omega) * real(PI * self % dSigmaS%re, defFlt)
      self % dSigmaT%re   = 0.0_defFlt
      self % dNuSigmaF%re = 0.0_defFlt
      self % dSigmaS%re   = 0.0_defFlt
    else
      self % dSigmaT%re   = real(PI * self % dSigmaT%re, defFlt)
      self % dNuSigmaF%re = real(PI * self % dNuSigmaF%re, defFlt)
      self % dSigmaS%re   = real(PI * self % dSigmaS%re, defFlt)
      self % dSigmaT%im   = 0.0_defFlt
      self % dNuSigmaF%im = 0.0_defFlt
      self % dSigmaS%im   = 0.0_defFlt
    end if

  end subroutine initialiseNoise

          
  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self

    call self % printSettings()
    call self % cycles()
    call self % printResults()

  end subroutine run

  !!
  !! Perform cycles of The Random Ray Method for a combined noise/eigenvalue solve.
  !! This subroutine looks just like the eigenvalue only version.
  !!
  !! Randomly places the ray starting point and direction uniformly.
  !! Rays are tracked until they reach some specified termination length.
  !! During tracking, fluxes are attenuated (and adjusted according to BCs),
  !! scoring to fluxes and volume estimates when the ray has surpassed its 
  !! specified dead length.
  !!
  !! Inactive and active iterations occur, terminating subject either to 
  !! given criteria or when a fixed number of iterations has passed.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  subroutine cycles(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defFlt)                                  :: hitRate, ONE_KEFF
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: stoppingCriterion, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    !$omp threadprivate(pRNG, r, ints)

    print *,'Beginning calculation'

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise fluxes 
    self % keff       = ONE
    self % scalarFlux = 0.0_defFlt
    self % prevFlux   = 1.0_defFlt
    self % fluxScores = ZERO
    self % keffScore  = ZERO
    self % source     = 0.0_defFlt
    
    self % scalarNoise = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    self % prevNoise   = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    self % noiseScores = cmplx(ZERO,ZERO,defReal)
    self % sourceNoise = cmplx(0.0_defFlt,0.0_defFlt,defFlt)

    ! Initialise other results
    self % cellHit      = 0
    self % volume       = ZERO
    self % volumeTracks = ZERO
    self % intersectionsTotal  = 0

    ! Initialise cell information
    self % cellFound = .false.
    self % cellPos = -INFINITY

    ! Stopping criterion is initially on flux convergence or number of convergence iterations.
    ! Will be replaced by RMS error in flux or number of scoring iterations afterwards.
    itInac = 0
    itAct  = 0
    isActive = .false.
    stoppingCriterion = .true.
    
    ! Power iteration
    do while( stoppingCriterion )
      
      if (isActive) then
        itAct = itAct + 1
      else
        itInac = itInac + 1
      end if
      it = itInac + itAct
      
      ONE_KEFF = real(ONE / self % keff,defFlt)
      !$omp parallel do schedule(static)
      do i = 1, self % nCells
        call self % sourceUpdateKernel(i, ONE_KEFF)
      end do
      !$omp end parallel do
      
      ! Reset and start transport timer
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
      intersections = 0
      
      !$omp parallel do schedule(dynamic) reduction(+: intersections)
      do i = 1, self % pop

        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 

        ! Set ray attributes
        call self % initialiseRay(r)

        ! Transport ray until termination criterion met
        call self % transportSweep(r,ints)
        intersections = intersections + ints

      end do
      !$omp end parallel do

      self % intersectionsTotal = self % intersectionsTotal + intersections
      
      call timerStop(self % timerTransport)

      ! Update RNG on master thread
      call self % rand % stride(self % pop + 1)

      ! Normalise flux estimate and combines with source
      call self % normaliseFluxAndVolume(it)

      ! Calculate new k
      call self % calculateKeff()

      ! Accumulate flux scores
      if (isActive) call self % accumulateFluxAndKeffScores()

      ! Calculate proportion of cells that were hit
      hitRate = real(sum(self % cellHit),defFlt) / self % nCells
      self % cellHit = 0
      
      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        stoppingCriterion = (itAct < self % active)
      elseif (self % doEntropy) then
        if (allocated(self % entropyMap)) then
          call self % calculateEntropyMap(isActive, it)
        else
          call self % calculateEntropy(isActive, it)
        end if
      else
        isActive = (itInac >= self % inactive)
      end if

      ! Set previous iteration fluxes to current fluxes
      ! and zero current fluxes
      call self % resetFluxes()

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)
      transport_T = timerTime(self % timerTransport)
      self % time_transport = self % time_transport + transport_T

      ! Predict time to end
      if (self % doEntropy) then
        end_T = real(self % active + itInac, defReal) * elapsed_T / it
      else
        end_T = real(self % active + self % inactive, defReal) * elapsed_T / it
      end if
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(it)
      print *
      if (self % doEntropy) then
        if (isActive) then
          print *, 'Iteration: ', numToChar(it),' of ', numToChar(self % active + itInac)
        else
          print *, 'Iteration: ', numToChar(it)
        end if
      else
        print *, 'Iteration: ', numToChar(it), ' of ', numToChar(self % active + self % inactive)
      end if
      if(isActive) then
        print *,'Active iterations'
      else
        print *,'Inactive iterations'
        if (self % doEntropy) print *,'Noise phase sum: ',trim(numToChar(self % entropyNoiseVec(1)))
      end if
      print *, 'Cell hit rate: ', trim(numToChar(real(hitRate,defReal)))
      print *, 'keff: ', trim(numToChar(real(self % keff,defReal)))
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      print *, 'Time per integration (ns): ', &
              trim(numToChar(transport_T*10**9/(self % nG * intersections)))

    end do

    ! Finalise flux scores
    call self % finaliseFluxAndKeffScores(itAct)

  end subroutine cycles

  !!
  !! Initialises rays: samples initial position and direction,
  !! and performs the build operation
  !!
  subroutine initialiseRay(self, r)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, x
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (randomRayPhysicsPackage_class.f90)'

    i = 0
    mu = TWO * r % pRNG % get() - ONE
    phi = TWO_PI * r % pRNG % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    rejection : do
      rand3(1) = r % pRNG % get()
      rand3(2) = r % pRNG % get()
      rand3(3) = r % pRNG % get()
      x = self % bottom + (self % top - self % bottom) * rand3

      ! Exit if point is inside the geometry
      call self % geom % whatIsAt(matIdx, cIdx, x, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u, 1, ONE)
    call self % geom % placeCoord(r % coords)

    if (.not. self % cellFound(cIdx)) then
      !$omp critical 
      self % cellFound(cIdx) = .true.
      self % cellPos(cIdx,:) = x
      !$omp end critical
    end if

  end subroutine initialiseRay

  !!
  !! Moves ray through geometry, updating angular flux (static and noise)
  !! and scoring scalar flux (static and noise) and volume.
  !! Records the number of intersections/ray movements.
  !!
  subroutine transportSweep(self, r, ints)
    class(noiseTRRMPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, idx, event, matIdx0, baseIdx, baseIdxNG, mIdx
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, fluxVec, fluxMinSour, omegaV
    complex(defFlt), dimension(self % nG)                 :: attCmplx, deltaN, fluxNVec, iVOmega
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec, speedVec
    complex(defFlt), pointer, dimension(:)                :: scalarNVec, sourceNVec, dTotVec
    real(defReal), dimension(3)                           :: r0, mu0
    real(defFlt)                                          :: lenFlt

    
    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % source(idx)
      fluxNVec(g) = self % sourceNoise(idx)
    end do

    ints = 0
    matIdx0 = 0
    totalLength = ZERO
    activeRay = .false.
    do while (totalLength < self % termination)

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        
        ! Cache total cross sections, speed, and frequency quantities
        mIdx = (matIdx - 1) * self % nG
        totVec => self % sigmaT((mIdx + 1):(mIdx + self % nG))
        dTotVec => self % dSigmaT((mIdx + 1):(mIdx + self % nG))
        speedVec => self % speed((mIdx + 1):(mIdx + self % nG))
    
        !$omp simd
        do g = 1, self % nG
          omegaV(g) = self % omega/speedVec(g)
          iVOmega(g) = imag * speedVec(g) / self % omega
        end do
      end if

      ! Remember co-ordinates to set new cell's position
      if (.not. self % cellFound(cIdx)) then
        r0 = r % rGlobal()
        mu0 = r % dirGlobal()
      end if
          
      ! Set maximum flight distance and ensure ray is active
      if (totalLength >= self % dead) then
        length = self % termination - totalLength 
        activeRay = .true.
      else
        length = self % dead - totalLength
      end if

      ! Move ray
      ! Use distance caching or standard ray tracing
      ! Distance caching seems a little bit more unstable
      ! due to FP error accumulation, but is faster.
      ! This can be fixed by resetting the cache after X number
      ! of distance calculations.
      if (self % cache) then
        if (mod(ints,100_longInt) == 0)  cache % lvl = 0
        call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
      else
        call self % geom % moveRay_noCache(r % coords, length, event, hitVacuum)
      end if
      totalLength = totalLength + length
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = r0 + HALF * length * mu0
        !$omp end critical
      end if

      ints = ints + 1

      ! Get mesh cell quantities for calculation
      baseIdx = (cIdx - 1) * self % nG + 1
      baseIdxNG = baseIdx + self % nG
      sourceVec => self % source(baseIdx:baseIdxNG)
      sourceNVec => self % sourceNoise(baseIdx:baseIdxNG)
      
      lenFlt = real(length,defFlt)

      ! Calculate real and complex attentuations
      !$omp simd
      do g = 1, self % nG
        attenuate(g) = exponential(totVec(g) * lenFlt)
        attCmplx(g) = 1.0_defFlt - (1.0_defFlt - attenuate(g)) * cmplx(cos(omegaV(g)*lenFlt),-sin(omegaV(g)*lenFlt),defFlt)
      end do
      
      ! Calculate angular flux update for static equation
      !$omp simd
      do g = 1, self % nG
        fluxMinSour(g) = fluxVec(g) - sourceVec(g) 
        delta(g) = fluxMinSour(g) * attenuate(g)
        fluxVec(g) = fluxVec(g) - delta(g)
      end do
      
      ! Calculate angular flux update for noise equation
      !$omp simd
      do g = 1, self % nG
        deltaN(g) = (fluxNVec(g) - sourceNVec(g) - iVOmega(g) * dTotVec(g) * fluxMinSour(g)) * attCmplx(g)
        fluxNVec(g) = fluxNVec(g) - deltaN(g) - iVOmega(g) * dTotVec(g) * delta(g)
      end do
      
      ! Accumulate to scalar flux
      if (activeRay) then
      
        call OMP_set_lock(self % locks(cIdx))
      
        scalarVec => self % scalarFlux(baseIdx:baseIdxNG)
        scalarNVec => self % scalarNoise(baseIdx:baseIdxNG)
        
        !$omp simd
        do g = 1, self % nG
          scalarVec(g) = scalarVec(g) + delta(g) 
          scalarNVec(g) = scalarNVec(g) + deltaN(g)
        end do
        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
        call OMP_unset_lock(self % locks(cIdx))

        if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1
      
      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, self % nG
          fluxVec(g) = 0.0_defFlt
          fluxNVec(g) = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
        end do
      end if

    end do

  end subroutine transportSweep

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    real(defReal), save                           :: vol
    complex(defFlt), save                         :: iOmegaV, dTotal, iVOmega
    real(defFlt), save                            :: total
    integer(shortInt), save                       :: g, matIdx, idx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(total, dTotal, iOmegaV, iVOmega, idx, g, matIdx, vol)

    norm = real(ONE / self % lengthPerIt, defFlt)
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = self % volume(cIdx)
      
      do g = 1, self % nG

        total = self % sigmaT((matIdx - 1) * self % nG + g)
        dTotal = self % dSigmaT((matIdx - 1) * self % nG + g)
        iOmegaV = imag * self % omega / self % speed((matIdx - 1) * self % nG + g)
        iVOmega = imag * self % speed((matIdx - 1) * self % nG + g) / self % omega

        idx   = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
          ! Normalise scalar flux by number of tracks and total XS
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm / ( total * real(vol,defFlt))
          
          ! Normalise scalar noise by number of tracks and total complex XS
          self % scalarNoise(idx) = self % scalarNoise(idx) * norm / ( (total + iOmegaV) * real(vol,defFlt))
          
          ! Increment noise by scaled scalar flux
          self % scalarNoise(idx) = self % scalarNoise(idx) + iVOmega * self % scalarFlux(idx) * dTotal
        end if

        ! Increment both by appropriate source term
        self % scalarFlux(idx) = self % scalarFlux(idx) + self % source(idx) 
        self % scalarNoise(idx) = self % scalarNoise(idx) + self % sourceNoise(idx)

      end do
      
    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF)
    class(noiseTRRMPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    real(defFlt)                                          :: scatter, fission, ONE_MIN_B
    complex(defFlt)                                       :: scatterNoise, fissionNoise, spectrum, iOmega, lambdaBetaChi
    real(defFlt), dimension(:), pointer                   :: nuFission, total, chi, scatterXS 
    real(defFlt), dimension(:), pointer                   :: fluxVec, scatterVec
    real(defFlt), dimension(:), pointer                   :: beta, lambda, chiP, chiD, chiDVec, speed
    complex(defFlt), dimension(:), pointer                :: dNuFission, dTotal, dScatterXS, noiseVec, dScatterVec
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, delayIdx, scatIdx, p, &
                                                             matIdx0, matIdxNG, sMatIdx0, sMatIdxNG, pIdx0, pIdxNG

    iOmega = imag * self % omega
    
    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 

    ! Obtain XSs
    delayIdx = (matIdx - 1) * self % nP
    beta => self % beta((delayIdx + 1):(delayIdx + self % nP))
    ONE_MIN_B = 1.0_defFlt - sum(beta)
    lambda => self % lambda((delayIdx + 1):(delayIdx + self % nP))
    
    matIdx0 = (matIdx - 1) * self % nG + 1
    matIdxNG = matIdx * self % nG
    
    speed => self % speed(matIdx0:matIdxNG)

    total => self % sigmaT(matIdx0:matIdxNG)
    dTotal => self % dSigmaT(matIdx0:matIdxNG)
    
    nuFission => self % nuSigmaF(matIdx0:matIdxNG)
    dNuFission => self % dNuSigmaF(matIdx0:matIdxNG)
    
    sMatIdx0 = (matIdx - 1) * self % nG * self % nG + 1
    sMatIdxNG = matIdx * self % nG * self % nG 
    
    scatterXS => self % sigmaS(sMatIdx0:sMatIdxNG)
    dScatterXS => self % dSigmaS(sMatIdx0:sMatIdxNG)
    
    chi => self % chi(matIdx0:matIdxNG)
    chiP => self % chiP(matIdx0:matIdxNG)

    pIdx0 = delayIdx * self % nG + 1
    pIdxNG = delayIdx * self % nG + self % nG * self % nP
    
    chiD => self % chiD(pIdx0:pIdxNG)

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx + 1):(baseIdx + self % nG))
    noiseVec => self % prevNoise((baseIdx + 1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    fissionNoise = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    !$omp simd reduction(+:fission, fissionNoise)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
      fissionNoise = fissionNoise + noiseVec(gIn) * nuFission(gIn) + fluxVec(gIn) * dNuFission(gIn)
    end do
    fission = fission * ONE_KEFF
    fissionNoise = fissionNoise * ONE_KEFF

    do g = 1, self % nG

      scatIdx = self % nG * (g - 1)
      scatterVec => scatterXS((scatIdx + 1):(scatIdx + self % nG))
      dScatterVec => dScatterXS((scatIdx + 1):(scatIdx + self % nG))

      ! Calculate scattering source
      scatter = 0.0_defFlt
      scatterNoise = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
      !$omp simd reduction(+:scatter, scatterNoise)
      do gIn = 1, self % nG
        scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
        scatterNoise = scatterNoise + fluxVec(gIn) * dScatterVec(gIn) + noiseVec(gIn) * scatterVec(gIn)
      end do

      ! Output index
      idx = baseIdx + g

      ! Calculate standard source
      self % source(idx) = chi(g) * fission + scatter
      self % source(idx) = self % source(idx) / total(g) 

      ! Calculate complex fission spectrum for noise
      chiDVec => chiD((self % nP * (g - 1) + 1):(self % nP * (g - 1) + self % nP))
      lambdaBetaChi = 0.0_defFlt
      !$omp simd reduction(+:lambdaBetaChi)
      do p = 1, self % nP
        lambdaBetaChi = lambdaBetaChi + lambda(p) * beta(p) * chiDVec(p) / (iOmega + lambda(p))
      end do
      spectrum = ONE_MIN_B * chiP(g) + lambdaBetaChi
      
      ! Calculate noise source
      self % sourceNoise(idx) = spectrum * fissionNoise + scatterNoise &
              - dTotal(g) * self % source(idx) + self % fixedSourceNoise(idx)
      self % sourceNoise(idx) = self % sourceNoise(idx) / (total(g) + iOmega / speed(g))

    end do

  end subroutine sourceUpdateKernel

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    real(defReal)                                 :: fissionRate, prevFissionRate
    real(defReal), save                           :: fissLocal, prevFissLocal, vol
    integer(shortInt), save                       :: matIdx, g, idx, mIdx
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, fissLocal, prevFissLocal, matIdx, g, idx, mIdx, vol)

    fissionRate     = ZERO
    prevFissionRate = ZERO
    !$omp parallel do schedule(static) reduction(+: fissionRate, prevFissionRate)
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)
      if (.not. mat % isFissile()) cycle

      vol = self % volume(cIdx)

      if (vol <= volume_tolerance) cycle

      fissLocal = ZERO
      prevFissLocal = ZERO
      mIdx = (matIdx - 1) * self % nG
      do g = 1, self % nG
        
        ! Source index
        idx = self % nG * (cIdx - 1) + g
        fissLocal     = fissLocal     + self % scalarFlux(idx) * self % nuSigmaF(mIdx + g)
        prevFissLocal = prevFissLocal + self % prevFlux(idx) * self % nuSigmaF(mIdx + g)

      end do

      fissionRate     = fissionRate     + fissLocal * vol
      prevFissionRate = prevFissionRate + prevFissLocal * vol

    end do
    !$omp end parallel do

    ! Update k
    self % keff = self % keff * fissionRate / prevFissionRate

  end subroutine calculateKeff
  
  !!
  !! Calculate the entropy
  !!
  subroutine calculateEntropy(self, isActive, it)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    logical(defBool), intent(inout)               :: isActive
    integer(shortInt), intent(in)                 :: it
    real(defReal)                                 :: sumFlux, entropyFlux, entropyNoise, &
                                                     meanF1, meanF2, meanN1, meanN2, sdF1, sdN1
    real(defReal), save                           :: fluxLocal, noiseAmpLocal, vol, pF, pN
    integer(shortInt), save                       :: g, idx
    integer(shortInt)                             :: cIdx, i
    !$omp threadprivate(fluxLocal, noiseAmpLocal, g, idx, vol)

    ! Calculate full flux and noise amplitude
    sumFlux = ZERO
    !$omp parallel do schedule(static) reduction(+:sumFlux)
    do cIdx = 1, self % nCells
      
      vol = self % volume(cIdx)
      if (vol <= volume_tolerance) cycle

      fluxLocal = ZERO
      noiseAmpLocal = ZERO
      do g = 1, self % nG
        
        ! Flux index
        idx = self % nG * (cIdx - 1) + g
        fluxLocal = fluxLocal + self % scalarFlux(idx)

      end do

      sumFlux = sumFlux + fluxLocal * vol

    end do
    !$omp end parallel do
    
    ! Calculate the entropy
    entropyFlux  = ZERO
    entropyNoise = ZERO
    !$omp parallel do schedule(static) reduction(+:entropyFlux, entropyNoise)
    do cIdx = 1, self % nCells

      vol = self % volume(cIdx)
      if (vol <= volume_tolerance) cycle

      fluxLocal = ZERO
      noiseAmpLocal = ZERO
      do g = 1, self % nG
        
        ! Flux index
        idx = self % nG * (cIdx - 1) + g
        fluxLocal     = fluxLocal     + self % scalarFlux(idx)
        noiseAmpLocal = noiseAmpLocal + atan2(self % scalarNoise(idx) % im,self % scalarNoise(idx)%re)

      end do
      
      pF = fluxLocal * vol / sumFlux
      if (pF > 0) entropyFlux  = entropyFlux  - pF * log(pF)
      pN = noiseAmpLocal * vol
      entropyNoise = entropyNoise + pN 

    end do
    !$omp end parallel do

    ! Juggle entropy values
    self % entropyFluxVec(2:entLength) = self % entropyFluxVec(1:(entLength - 1))
    self % entropyFluxVec(1) = entropyFlux
    self % entropyNoiseVec(2:entLength) = self % entropyNoiseVec(1:(entLength - 1))
    self % entropyNoiseVec(1) = entropyNoise

    ! If 200 iterations have occurred, check for convergence 
    ! Compare mean + 1 SD of the first 100 values with the last 100 values
    if (it >= entLength) then

      meanF1 = ZERO
      meanF2 = ZERO
      meanN1 = ZERO
      meanN2 = ZERO
      sdF1 = ZERO
      sdN1 = ZERO
      !$omp parallel do schedule(static) reduction(+:meanF1,sdF1,meanF2,meanN1,sdN1,meanN2) 
      do i = 1, entLength/2
        meanF1 = meanF1 + self % entropyFluxVec(i)
        sdF1   = sdF1 + self % entropyFluxVec(i) * self % entropyFluxVec(i)
        meanF2 = meanF2 + self % entropyFluxVec(i + entLength/2)
        meanN1 = meanN1 + self % entropyNoiseVec(i)
        sdN1   = sdN1 + self % entropyNoiseVec(i) * self % entropyNoiseVec(i)
        meanN2 = meanN2 + self % entropyNoiseVec(i + entLength/2)
      end do
      !$omp end parallel do

      meanF1 = meanF1 / (entLength/2)
      meanF2 = meanF2 / (entLength/2)
      sdF1   = sdF1 / (entLength/2)
      sdF1   = sqrt((sdF1 - meanF1 * meanF1) / (entLength/2 -1))
      meanN1 = meanN1 / (entLength/2)
      meanN2 = meanN2 / (entLength/2)
      sdN1   = sdN1 / (entLength/2)
      sdN1   = sqrt((sdN1 - meanN1 * meanN1) / (entLength/2 -1))

      if ( abs(meanF1 - meanF2) < sdF1 .and. abs(meanN1 - meanN2) < sdN1) isActive = .true.

    end if

  end subroutine calculateEntropy
  
  !!
  !! Calculate the entropy based no a specified entropy map
  !! Hopefully a bit harder to satisfy...
  !!
  subroutine calculateEntropyMap(self, isActive, it)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    logical(defBool), intent(inout)               :: isActive
    integer(shortInt), intent(in)                 :: it
    real(defReal), save                           :: fluxLocal, noiseAmpLocal, vol, pF, pN
    real(defReal), dimension(self % entMapSize)   :: sumFlux, entropyFlux, entropyNoise, &
                                                     meanF1, meanF2, meanN1, meanN2, sdF1, sdN1
    integer(shortInt), save                       :: g, idx, j
    integer(shortInt)                             :: cIdx, i
    type(particleState), save                     :: s
    logical(defBool)                              :: passFlux, passNoise
    !$omp threadprivate(fluxLocal, noiseAmpLocal, g, idx, vol, s, j, pF, pN)

    ! Calculate full flux and noise amplitude in every entropy mesh element
    sumFlux = ZERO
    !$omp parallel do schedule(static) reduction(+:sumFlux)
    do cIdx = 1, self % nCells
      
      vol = self % volume(cIdx)
      if (vol <= volume_tolerance) cycle

      fluxLocal = ZERO
      noiseAmpLocal = ZERO
      do g = 1, self % nG
        
        ! Flux index
        idx = self % nG * (cIdx - 1) + g
        fluxLocal = fluxLocal + self % scalarFlux(idx)

      end do

      ! Map to the local entropy bin
      ! Fudge a particle state to search tally map
      s % r = self % cellPos(cIdx,:)
      idx = self % entropyMap % map(s)
      if (idx > 0) then
        sumFlux(idx) = sumFlux(idx) + fluxLocal * vol
      end if

    end do
    !$omp end parallel do

    ! Calculate the entropy
    entropyFlux  = ZERO
    entropyNoise = ZERO
    !$omp parallel do schedule(static) reduction(+:entropyFlux, entropyNoise)
    do cIdx = 1, self % nCells

      vol = self % volume(cIdx)
      if (vol <= volume_tolerance) cycle

      fluxLocal = ZERO
      noiseAmpLocal = ZERO
      do g = 1, self % nG
        
        ! Flux index
        idx = self % nG * (cIdx - 1) + g
        fluxLocal     = fluxLocal     + self % scalarFlux(idx)
        noiseAmpLocal = noiseAmpLocal + atan2(self % scalarNoise(idx) % im,self % scalarNoise(idx)%re)

      end do

      ! Map to the local entropy bin
      ! Fudge a particle state to search tally map
      s % r = self % cellPos(cIdx,:)
      idx = self % entropyMap % map(s)
      
      if (idx > 0) then
        pF = fluxLocal * vol / sumFlux(idx)
        if (pF > 0) entropyFlux(idx)  = entropyFlux(idx)  - pF * log(pF)
        pN = noiseAmpLocal * vol 
        entropyNoise(idx) = entropyNoise(idx) + pN 
      end if 

    end do
    !$omp end parallel do

    ! Juggle entropy values
    self % entropy(2:entLength,:) = self % entropy(1:(entLength - 1),:)
    self % entropy(1,:) = entropyFlux
    self % entropyNoise(2:entLength,:) = self % entropyNoise(1:(entLength - 1),:)
    self % entropyNoise(1,:) = entropyNoise

    ! If 200 iterations have occurred, check for convergence 
    ! Compare mean + 1 SD of the first 100 values with the last 100 values
    if (it >= entLength) then

      meanF1 = ZERO
      meanF2 = ZERO
      meanN1 = ZERO
      meanN2 = ZERO
      sdF1 = ZERO
      sdN1 = ZERO
      !$omp parallel do schedule(static) reduction(+:meanF1,sdF1,meanF2,meanN1,sdN1,meanN2) 
      do i = 1, self % entMapSize
        do j = 1, entLength/2
          meanF1(i) = meanF1(i) + self % entropy(j,i)
          sdF1(i)   = sdF1(i) + self % entropy(j,i) * self % entropy(j,i)
          meanF2(i) = meanF2(i) + self % entropy(j + entLength/2,i)
          meanN1(i) = meanN1(i) + self % entropyNoise(j,i)
          sdN1(i)   = sdN1(i) + self % entropyNoise(j,i) * self % entropyNoise(j,i)
          meanN2(i) = meanN2(i) + self % entropyNoise(j + entLength/2,i)
        end do
      end do
      !$omp end parallel do

      meanF1 = meanF1 / (entLength/2)
      meanF2 = meanF2 / (entLength/2)
      sdF1   = sdF1 / (entLength/2)
      sdF1   = sqrt((sdF1 - meanF1 * meanF1) / (entLength/2 - 1))
      meanN1 = meanN1 / (entLength/2)
      meanN2 = meanN2 / (entLength/2)
      sdN1   = sdN1 / (entLength/2)
      sdN1   = sqrt((sdN1 - meanN1 * meanN1) / (entLength/2 -1))

      ! Do 68% of the values fall within 1SD?
      passFlux = (real(count(abs(meanF1 - meanF2) < sdF1),defReal)/self % entMapSize > 0.68)
      passNoise = (real(count(abs(meanN1 - meanN2) < sdN1),defReal)/self % entMapSize > 0.68)
      if (passFlux .and. passNoise) isActive = .true.

    end if

  end subroutine calculateEntropyMap
  
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxes(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt)                             :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % prevNoise(idx) = self % scalarNoise(idx)
      self % scalarFlux(idx) = 0.0_defFlt
      self % scalarNoise(idx) = cmplx(0.0_defFlt,0.0_defFlt,defFlt)
    end do
    !$omp end parallel do

  end subroutine resetFluxes

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    real(defReal), save                           :: flux
    complex(defReal), save                        :: noise
    integer(shortInt)                             :: idx
    !$omp threadprivate(flux, noise)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = real(self % scalarFlux(idx),defReal)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux*flux
      noise = cmplx(self % scalarNoise(idx)%re,self % scalarNoise(idx)%im,defReal)
      self % noiseScores(idx,1) = self % noiseScores(idx, 1) + noise
      self % noiseScores(idx,2) = self % noiseScores(idx, 2) + noise%re*noise%re + noise%im*noise%im
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) + self % keff
    self % keffScore(2) = self % keffScore(2) + self % keff * self % keff

  end subroutine accumulateFluxAndKeffScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxAndKeffScores(self,it)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt)                             :: idx
    real(defReal)                                 :: N1, Nm1

    if (it /= 1) then
      Nm1 = ONE/(it - 1)
    else
      Nm1 = ONE
    end if
    N1 = ONE/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) * N1
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) * N1
      self % fluxScores(idx,2) = Nm1 *(self % fluxScores(idx,2) - &
            self % fluxScores(idx,1) * self % fluxScores(idx,1)) 
      if (self % fluxScores(idx,2) <= ZERO) then
        self % fluxScores(idx,2) = ZERO
      else
        self % fluxScores(idx,2) = sqrt(self % fluxScores(idx,2))
      end if

      self % noiseScores(idx,1) = self % noiseScores(idx, 1) * N1
      self % noiseScores(idx,2) = self % noiseScores(idx, 2) * N1
      self % noiseScores(idx,2) = Nm1 *(self % noiseScores(idx,2) - &
            abs(self % noiseScores(idx,1)) * abs(self % noiseScores(idx,1))) 
      if (abs(self % noiseScores(idx,2)) <= ZERO) then
        self % noiseScores(idx,2) = ZERO
      else
        self % noiseScores(idx,2) = sqrt(self % noiseScores(idx,2))
      end if
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) * N1
    self % keffScore(2) = self % keffScore(2) * N1
    self % keffScore(2) = sqrt(Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    integer(shortInt)                             :: cIdx, g1
    integer(shortInt), save                       :: idx, matIdx, i, g
    real(defFlt), save                            :: SigmaF
    real(defReal), save                           :: vol
    type(particleState), save                     :: s
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal), dimension(:), allocatable      :: groupFlux, fiss, fissSTD, totVol
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(idx, matIdx, i, mat, matPtr, vol, s, SigmaF, g)

    call out % init(self % outputFormat)
    
    name = 'seed'
    call out % printValue(self % rand % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % active,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Total_Transport_Time'
    call out % printValue(self % time_transport,name)
    
    name = 'Time_Per_Integration'
    call out % printValue(self % time_transport/(self % intersectionsTotal * self % nG),name)
    
    name = 'Clock_Time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print keff
    name = 'keff'
    call out % startBlock(name)
    call out % printResult(self % keffScore(1), self % keffScore(2), name)
    call out % endBlock()

    ! Print cell volumes
    if (self % printVolume) then
      name = 'volume'
      call out % startBlock(name)
      resArrayShape = [size(self % volume)]
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        call out % addResult(self % volume(cIdx), ZERO)
      end do
      call out % endArray()
      call out % endBlock()
    end if

    ! Print cell positions
    if (self % printCells) then
      name = 'position'
      call out % startBlock(name)
      resArrayShape = [size(self % cellPos)]
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        call out % addResult(self % cellPos(cIdx,1), ZERO)
        call out % addResult(self % cellPos(cIdx,2), ZERO)
        call out % addResult(self % cellPos(cIdx,3), ZERO)
      end do
      call out % endArray()
      call out % endBlock()
    end if

    ! Print fluxes
    if (self % printFlux) then
      resArrayShape = [size(self % volume)]
      do g = 1, self % nG
        name = 'flux_g'//numToChar(g)
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g
          call out % addResult(self % fluxScores(idx,1), self % fluxScores(idx,2))
        end do
        call out % endArray()
        call out % endBlock()
      end do
      resArrayShape = [size(self % volume)]
      do g = 1, self % nG
        name = 'realNoise_g'//numToChar(g)
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g
          call out % addResult(self % noiseScores(idx,1)%re, abs(self % noiseScores(idx,2)))
        end do
        call out % endArray()
        call out % endBlock()
      end do
      resArrayShape = [size(self % volume)]
      do g = 1, self % nG
        name = 'imagNoise_g'//numToChar(g)
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g
          call out % addResult(self % noiseScores(idx,1)%im, abs(self % noiseScores(idx,2)))
        end do
        call out % endArray()
        call out % endBlock()
      end do
    end if

    ! Send fission rates to map output
    if (self % mapFission) then
      resArrayShape = self % resultsMap % binArrayShape()
      allocate(fiss(self % resultsMap % bins(0)))
      allocate(fissSTD(self % resultsMap % bins(0)))
      fiss    = ZERO
      fissSTD = ZERO

      ! Find whether cells are in map and sum their contributions
      !$omp parallel do reduction(+: fiss, fissSTD)
      do cIdx = 1, self % nCells
        
        ! Identify material
        matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
        matPtr => self % mgData % getMaterial(matIdx)
        mat    => baseMgNeutronMaterial_CptrCast(matPtr)
        vol    =  self % volume(cIdx)

        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % resultsMap % map(s)

        if (i > 0) then
          do g = 1, self % nG
            SigmaF = real(mat % getFissionXS(g, self % rand),defFlt)
            idx = (cIdx - 1)* self % nG + g
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1) * SigmaF
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    vol * vol * self % fluxScores(idx,2)*self % fluxScores(idx,2) * SigmaF * SigmaF
          end do
        end if

      end do
      !$omp end parallel do

      do i = 1,size(fissSTD)
        fissSTD(i) = sqrt(fissSTD(i))
        if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
      end do

      name = 'fissionRate'
      call out % startBlock(name)
      call out % startArray(name, resArrayShape)
      ! Add all map elements to results
      do idx = 1, self % resultsMap % bins(0)
        call out % addResult(fiss(idx), fissSTD(idx))
      end do
      call out % endArray()
      ! Output tally map
      call self % resultsMap % print(out)
      call out % endBlock()
      
      deallocate(fiss)
      deallocate(fissSTD)
    end if

    ! Send fluxes to map output
    if (self % mapFlux) then

      !Forward fluxes
      resArrayShape = self % fluxMap % binArrayShape()
      allocate(fiss(self % fluxMap % bins(0)))
      allocate(fissSTD(self % fluxMap % bins(0)))
      allocate(totVol(self % fluxMap % bins(0)))

      name = 'flux'
      call out % startBlock(name)
      
      ! Find whether cells are in map and sum their contributions
      do g1 = 1, self % nG
        
        fiss    = ZERO
        fissSTD = ZERO

        !$omp parallel do reduction(+: fiss, fissSTD)
        do cIdx = 1, self % nCells
        
          vol    =  self % volume(cIdx)
          if (vol < volume_tolerance) cycle

          ! Fudge a particle state to search tally map
          s % r = self % cellPos(cIdx,:)
          i = self % fluxMap % map(s)

          if (i > 0) then
            idx = (cIdx - 1)* self % nG + g1
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1)
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    self % fluxScores(idx,2)*self % fluxScores(idx,2) * vol * vol
          end if

        end do
        !$omp end parallel do

        do i = 1,size(fissSTD)
          fissSTD(i) = sqrt(fissSTD(i))
          if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
        end do

        name = 'g'//numToChar(g1)
        call out % startArray(name, resArrayShape)
        ! Add all map elements to results
        do idx = 1, self % fluxMap % bins(0)
          call out % addResult(fiss(idx), fissSTD(idx))
        end do
        call out % endArray()

      end do
      
      ! Output tally map
      call self % fluxMap % print(out)
      call out % endBlock()
        
      name = 'noiseReal'
      call out % startBlock(name)

      ! Real part of noise
      do g1 = 1, self % nG
        
        fiss    = ZERO
        fissSTD = ZERO

        !$omp parallel do reduction(+: fiss, fissSTD)
        do cIdx = 1, self % nCells
        
          vol    =  self % volume(cIdx)
          if (vol < volume_tolerance) cycle

          ! Fudge a particle state to search tally map
          s % r = self % cellPos(cIdx,:)
          i = self % fluxMap % map(s)

          if (i > 0) then
            idx = (cIdx - 1)* self % nG + g1
            fiss(i) = fiss(i) + vol * self % noiseScores(idx,1)%re
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    abs(self % noiseScores(idx,2))*abs(self % noiseScores(idx,2)) * vol * vol
          end if

        end do
        !$omp end parallel do

        do i = 1,size(fissSTD)
          fissSTD(i) = sqrt(fissSTD(i))
          if (abs(fiss(i)) > 0) fissSTD(i) = fissSTD(i) / abs(fiss(i))
        end do

        name = 'g'//numToChar(g1)
        call out % startArray(name, resArrayShape)
        ! Add all map elements to results
        do idx = 1, self % fluxMap % bins(0)
          call out % addResult(fiss(idx), fissSTD(idx))
        end do
        call out % endArray()

      end do
        
      ! Output tally map
      call self % fluxMap % print(out)
      call out % endBlock()

      ! Imaginary part of noise
      name = 'noiseImag'
      call out % startBlock(name)
      
      do g1 = 1, self % nG
        
        fiss    = ZERO
        fissSTD = ZERO
        totVol  = ZERO

        !$omp parallel do reduction(+: fiss, fissSTD, totVol)
        do cIdx = 1, self % nCells
        
          vol    =  self % volume(cIdx)
          if (vol < volume_tolerance) cycle

          ! Fudge a particle state to search tally map
          s % r = self % cellPos(cIdx,:)
          i = self % fluxMap % map(s)

          if (i > 0) then
            idx = (cIdx - 1)* self % nG + g1
            totVol(i) = totVol(i) + vol
            fiss(i) = fiss(i) + vol * self % noiseScores(idx,1)%im
            ! Is this correct? DEFINITELY FIX THIS
            !fissSTD(i) = ZERO
            fissSTD(i) = fissSTD(i) + &
                    abs(self % noiseScores(idx,2))*abs(self % noiseScores(idx,2)) * vol * vol
          end if

        end do
        !$omp end parallel do

        !do i = 1,size(fiss)
        !  fiss(i) = fiss(i) / totVol(i)
        !end do

        do i = 1,size(fissSTD)
          fissSTD(i) = sqrt(fissSTD(i))
          if (abs(fiss(i)) > 0) fissSTD(i) = fissSTD(i) / abs(fiss(i))
        end do

        name = 'g'//numToChar(g1)
        call out % startArray(name, resArrayShape)
        ! Add all map elements to results
        do idx = 1, self % fluxMap % bins(0)
          call out % addResult(fiss(idx), fissSTD(idx))
        end do
        call out % endArray()

      end do
        
      ! Output tally map
      call self % fluxMap % print(out)
      call out % endBlock()
      
      deallocate(fiss)
      deallocate(fissSTD)
    end if
    
    call out % writeToFile(self % outputFile)

    ! Send all fluxes and stds to VTK
    if (self % plotResults) then
      allocate(groupFlux(self % nCells))
      do g1 = 1, self % nG
        name = 'flux_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'std_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,2) /self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'noiseAmp_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = abs(self % noiseScores(idx,1))
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'stdAmp_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          if (abs(self % noiseScores(idx,1)) > ZERO) then
            groupFlux(cIdx) = abs(self % noiseScores(idx,2) /self % noiseScores(idx,1))
          else
            groupFlux(cIdx) = ZERO
          end if
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'noisePhase_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = atan(self % noiseScores(idx,1)%im,self % noiseScores(idx,1)%re)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'realNoise_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % noiseScores(idx,1)%re
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'imagNoise_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % noiseScores(idx,1)%im
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'imagSource_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = real(self % sourceNoise(idx)%im,defReal)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'realSource_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = real(self % sourceNoise(idx)%re,defReal)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      call self % viz % finaliseVTK
    end if


  end subroutine printResults

  !!
  !! Print settings of the random ray calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(noiseTRRMPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RANDOM RAY NOISE CALCULATION /\/\"
    print *, "Using "//numToChar(self % inactive)// " iterations for "&
              //"the inactive cycles"
    print *, "Using "//numToChar(self % active)// " iterations for "&
              //"the active cycles"
    print * 
    print *, "Rays per cycle: "// numToChar(self % pop)
    print *, "Ray dead length: "//numToChar(self % dead)
    print *, "Ray termination length: "//numToChar(self % termination)
    print *, "Initial RNG Seed:   "// numToChar(self % rand % getSeed())
    print *
    print *, "Number of cells in the geometry: "// numToChar(self % nCells)
    print *, "Number of energy groups: "// numToChar(self % nG)
    if (self % cache) print *, "Accelerated with distance caching"
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(noiseTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean Nuclear Data, Geometry and visualisation
    call gr_kill()
    call ndreg_kill()
    call self % viz % kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    self % timerMain = 0
    self % timerTransport = 0

    self % top       = ZERO
    self % bottom    = ZERO
    self % mgData    => null()
    self % nG        = 0
    self % nP        = 0
    self % omega     = -1.0_defFlt

    if(allocated(self % locks)) then
      do i = 1, self % nCells
        call OMP_destroy_lock(self % locks(i))
      end do
      deallocate(self % locks)
    end if
    self % nCells    = 0
    self % nMat      = 0
    if(allocated(self % sigmaT)) deallocate(self % sigmaT)
    if(allocated(self % sigmaS)) deallocate(self % sigmaS)
    if(allocated(self % nusigmaF)) deallocate(self % nuSigmaF)
    if(allocated(self % dSigmaT)) deallocate(self % dSigmaT)
    if(allocated(self % dSigmaS)) deallocate(self % dSigmaS)
    if(allocated(self % dNusigmaF)) deallocate(self % dNuSigmaF)
    if(allocated(self % chi)) deallocate(self % chi)
    if(allocated(self % chiP)) deallocate(self % chiP)
    if(allocated(self % chiD)) deallocate(self % chiD)
    if(allocated(self % lambda)) deallocate(self % lambda)
    if(allocated(self % beta)) deallocate(self % beta)
    if(allocated(self % speed)) deallocate(self % speed)

    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % inactive    = 0
    self % active      = 0
    self % cache       = .false.
    self % mapFission  = .false.
    self % mapFlux     = .false.
    self % plotResults = .false.
    self % printFlux   = .false.
    self % printVolume = .false.
    self % printCells  = .false.

    self % keff        = ZERO
    self % keffScore   = ZERO
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % scalarNoise)) deallocate(self % scalarNoise)
    if(allocated(self % prevNoise)) deallocate(self % prevNoise)
    if(allocated(self % noiseScores)) deallocate(self % noiseScores)
    if(allocated(self % sourceNoise)) deallocate(self % sourceNoise)
    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % cellHit)) deallocate(self % cellHit)

    if(allocated(self % resultsMap)) then
      call self % resultsMap % kill()
      deallocate(self % resultsMap)
    end if

    if(allocated(self % fluxMap)) then
      call self % resultsMap % kill()
      deallocate(self % fluxMap)
    end if

  end subroutine kill

end module noiseTRRMPhysicsPackage_class
