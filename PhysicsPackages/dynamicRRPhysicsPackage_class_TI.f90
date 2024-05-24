module dynamicRRPhysicsPackage_class_TI

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential
  use exponentialRA_TD_func,          only : exponential_TD
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile
  use rng_class,                      only : RNG
  use physicsPackage_inter,           only : physicsPackage
  use transient_class,                only : transientData

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryStd_class,              only : geometryStd
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use geometryFactory_func,           only : new_geometry

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
!--> DNP 230918
  use fissionMG_class,                only : fissionMG, fissionMG_TptrCast
!<-- DNP 230918
  
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

  !!
  !! Physics package to perform The Random Ray Method (TRRM) time-dependent calculations 
  !! using the time-implicit approach
  !!
  !! Tracks rays across the geometry, attenuating their flux. After some dead length,
  !! rays begin scoring to estimates of the scalar flux and volume. Each ray has a
  !! uniform termination length, after which it is stopped and the next ray is tracked.
  !! Once all rays have been tracked, a cycle concludes and fluxes, sources, and keff
  !! are updated.
!--> MK 240205
  !!
  !! This procedure is repeated for each time step, thus converging the spatial flux shape for 
  !! each time step individually. The segment-avergae flux for the subseqeuent time step is set 
  !! equal to an isotopric approximation using the scalar flux of the respective cell.
  !! 
  !! The first time step corresponds to the steady-state (enforced by diving nuSigmaf by keff).
  !! The initial flux is scaled such that the sum of the steady-state equals a specified input value.
  !!
  !! Both inactive and active cycles occur, as in Monte Carlo. These are terminated
  !! after a specified number of iterations and on reaching the Shannon entropy convergence
  !! criterion.
!<-- MK 240205
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
  !! Sample Input Dictionary:
  !!   PP {
  !!     type dynamicRRPhysicsPackage_TI;
  !!     dead 10;              // Dead length where rays do not score to scalar fluxes
  !!     termination 100;      // Length a ray travels before it is terminated
  !!     rays 1000;            // Number of rays to sample per iteration
  !!     inactive 100;         // Number of convergence cycles (would use accum and std otherwise)
  !!     active 200;           // Number of scoring cycles (would use eps otherwise)
!--> MK 240205
  !!     tstp 0.01;            // Time step size
  !!     nT 200;               // Number of time steps
  !!     tOut 100;             // Number of printed time steps for output files
  !!     #initVal 1;#          // Optional normalisation value for total flux in steady state
  !!     #window 100;#         // Optional width of moving window over which flux values are averaged for calculating RMS error
  !!     #reduceLines 1;#      // Optionally reduce printed output lines during run
  !!     #isotropicDeriv 1;#   // Optionally use isotropic approximation for entire flux time derivative
!<-- MK 240205
  !!     #seed 86868;#         // Optional RNG seed
  !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
  !!     #fissionMap {<map>}#  // Optionally output fission rates according to a given map
  !!     #fluxMap {<map>}#     // Optionally output one-group fluxes according to a given map
  !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
  !!
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition>}
  !!     transient {<Transient data definition>}
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
  !!   termination -> Distance a ray can travel before it is terminated
  !!   dead        -> Distance a ray must travel before it becomes active
  !!   pop         -> Number of rays to track per cycle
  !!   inactive    -> Number of inactive cycles to perform
  !!   active      -> Number of active cycles to perform
  !!   cache       -> Logical check whether to use distance caching
  !!   outputFile  -> Output file name
  !!   outputFormat-> Output file format
  !!   plotResults -> Plot results?
  !!   printFlux   -> Print fluxes?
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
  !!   chi         -> Local chi vector
  !!
  !!   keff        -> Estimated value of keff
  !!   keffScore   -> Vector holding cumulative keff score and keff^2 score
  !!   scalarFlux  -> Array of scalar flux values of length = nG * nCells
  !!   prevFlux    -> Array of previous scalar flux values of length = nG * nCells 
  !!   fluxScore   -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!   source      -> Array of neutron source values of length = nG * nCells
  !!   volume      -> Array of stochastically estimated cell volumes of length = nCells
  !!   cellHit     -> Array tracking whether given cells have been hit during tracking
  !!   cellFound   -> Array tracking whether a cell was ever found
  !!   cellPos     -> Array of cell positions, populated once they are found
  !!
  !!   locks       -> OpenMP locks used when writing to scalar flux during transport
!--> MK 230601
  !!   vel         -> Local neutron velocity vector
  !!   nu          -> Local average number of neutrons released per fission event
  !!   chiP        -> Local prompt chi vector 
  !!   chiD        -> Local delayed chi vector
  !!   lambda      -> Local DNP decay constants
  !!   beta        -> Local DNP group fractions
  !!   omega 0, omegaN, omega_1, omega_2 --> Auxiliary values stored for DNP derivative solution
  !!
  !!   nT          -> Number of time steps
  !!   nP          -> Number of DNP groups
  !!   tstp        -> Time step size
  !!   initVal     -> Initialisation value for scalar flux
  !!   window      -> Width of moving window for Shannon entropy
  !!   reduceLines -> Reduce output lines that are printed while running?
  !!   isotropic   -> Use isotropic approximation for entire flux time derivative?
  !!   totalIt     -> Total number of iterations in transient calculation to be reported in results
  !!   nTOut       -> Number of time steps for which results are printed
  !!   outInterval -> Interval of time steps at the end of which result is printed (Print every 'outInterval' time steps)
  !!
  !!   fluxHist    -> Array of accumulated and averaged scalar flux values of timestep n-1,
  !!                  dimension = [nG * nCells]
  !!   prevFission_1 -> Array of accumulated and averaged fission source values * keff of timestep n-1 to be used 
  !!                    for DNP calculation, dimension = [nCells]
  !!   prevFission_2 -> Array of accumulated and averaged fission source values * keff of timestep n-2 to be used 
  !!                    for DNP calculation, dimension = [nCells]
  !!   fissionScore  -> Array of accumulated and (eventually) averaged fission source values * keff of current timestep n to be used
  !!                    for DNP calculation, dimension = [nCells]
  !!   fissionSource -> Array of fission source values * keff , dimension = [nCells]
  !!   dnpScores   -> Array of accumulated and (eventually) averaged DNP values of current timestep n,
  !!                  dimension = [nP * nCells]
  !!   dnp         -> Array of DNP values of current timestep n, dimension = [nP * nCells]
  !!   prevDnp     -> Array of accumulated and averaged DNP values of timestep n-1,
  !!                  dimension = [nP * nCells]
  !!
  !!   scalarFluxTD -> same quantities as above but for time-dependent calculation
  !!   prevFluxTD  
  !!   fluxScoreTD
  !!   sourceTD 
  !!   fluxHistTD 
  !!   prevFission_TD1 
  !!   prevFission_TD2 
  !!   fissionScoreTD
  !!   fissionSourceTD
  !!   dnpScoresTD 
  !!   dnpTD
  !!   prevDnpTD
  !!
!<-- MK 230601
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: dynamicRRPhysicsPackage_TI
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
!---> DNP 230918
    integer(shortInt)                     :: nP = 0
!<-- DNP 230918

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .false.
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
!--> MK 230718
    integer(shortInt)  :: nT          = 0
    integer(shortInt)  :: outInterval = 0
    integer(shortInt)  :: window      = 0
    real(defFlt)       :: tstp        = ZERO
    real(defFlt)       :: initVal     = ZERO
    integer(shortInt)  :: totalIt     = 0        
    integer(shortInt)  :: totalAct    = 0         
    integer(shortInt)  :: totalInact  = 0         
    integer(shortInt)  :: nTOut       = 0
    logical(defBool)   :: reduceLines = .false.
    logical(defBool)   :: isotropic = .false.
!<-- MK 230718   
!--> MK 240523
    real(defReal)      :: AccVol = ZERO
    real(defFlt)       :: sarePrec = ZERO
    real(defFlt)       :: rmsPrec = ZERO
!<-- MK 240523

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:), allocatable     :: sigmaS
!--> MK 240517
    real(defFlt), dimension(:), allocatable     :: sigmaF
!<-- MK 240517
    real(defFlt), dimension(:), allocatable     :: chi
    real(defFlt), dimension(:), allocatable     :: nu

    real(defFlt), dimension(:), allocatable     :: chiP
    real(defFlt), dimension(:), allocatable     :: chiD
    real(defFlt), dimension(:), allocatable     :: lambda
    real(defFlt), dimension(:), allocatable     :: beta

    real(defFlt), dimension(:), allocatable     :: omega0
    real(defFlt), dimension(:), allocatable     :: omegaN
    real(defFlt), dimension(:), allocatable     :: omega_1
    real(defFlt), dimension(:), allocatable     :: omega_2

    real(defFlt), dimension(:), allocatable     :: vel

    ! Results space
    real(defFlt)                                :: keff
    real(defReal), dimension(2)                 :: keffScore  
    
    real(defFlt), dimension(:), allocatable     :: scalarFluxTD
    real(defFlt), dimension(:), allocatable     :: prevFluxTD
    real(defReal), dimension(:,:), allocatable  :: fluxScoresTD
    real(defFlt), dimension(:), allocatable     :: fluxHistTD
    real(defFlt), dimension(:), allocatable     :: sourceTD
    real(defFlt), dimension(:), allocatable     :: prevFissionTD_1
    real(defFlt), dimension(:), allocatable     :: prevFissionTD_2
!--> MK 240517
    real(defReal), dimension(:,:), allocatable  :: fissionRateScore
    real(defFlt), dimension(:), allocatable     :: fissionRate
!<-- MK 240517
    real(defReal), dimension(:), allocatable    :: fissionScoreTD
    real(defFlt), dimension(:), allocatable     :: fissionSourceTD
    real(defReal), dimension(:), allocatable    :: dnpScoresTD
    real(defFlt), dimension(:), allocatable     :: dnpTD
    real(defFlt), dimension(:), allocatable     :: prevDnpTD

    real(defFlt), dimension(:), allocatable     :: scalarFlux
    real(defFlt), dimension(:), allocatable     :: prevFlux
    real(defReal), dimension(:), allocatable    :: fluxScores
    real(defFlt), dimension(:), allocatable     :: fluxHist
    real(defFlt), dimension(:), allocatable     :: source
    real(defFlt), dimension(:), allocatable     :: prevFission_1
    real(defFlt), dimension(:), allocatable     :: prevFission_2
    real(defReal), dimension(:), allocatable    :: fissionScore
    real(defFlt), dimension(:), allocatable     :: fissionSource
    real(defReal), dimension(:), allocatable    :: dnpScores
    real(defFlt), dimension(:), allocatable     :: dnp
    real(defFlt), dimension(:), allocatable     :: prevDnp

    real(defReal), dimension(:), allocatable    :: volume
    real(defReal), dimension(:), allocatable    :: volumeTracks 

    ! Tracking cell properites
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos

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
    procedure, private :: cycles
    procedure, private :: initialiseRay
    procedure, private :: transportSweep_isotropic
    procedure, private :: transportSweep
    procedure, private :: calculateEntropy
    procedure, private :: sourceUpdateKernel_TD
    procedure, private :: sourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume_TD
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings
    procedure, private :: printTransientFlux
    procedure, private :: printFissionRate
    procedure, private :: updatePreviousQuantities
!--> MK 240517
    procedure, private :: calculateFissionRate
!<-- MK 240517

  end type dynamicRRPhysicsPackage_TI

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
!--> MK 230718
    integer(shortInt)                             :: seed_temp, i, g, g1, m, t, nP, nPNew, p, idx0, m2
    real(defReal)                                 :: tstpReal, initReal, outTemp, precReal
    type(transientData)                           :: variations      
    integer(shortInt)                             :: tStart, tEnd, baseIdx0, baseIdx1, baseIdx, tCnt
    real(defFlt)                                  :: magnitude
    character(nameLen)                            :: variationType, xsName
    real(defFlt), dimension(:), allocatable       :: xsType
    real(defReal)                                 :: kappa0, kappa1, kappa2
    real(defFlt), dimension(:), allocatable       :: origXS, startXS          
!<-- MK 230718
    integer(longInt)                              :: seed
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    class(dictionary),pointer                     :: tempDict, graphDict
    class(mgNeutronDatabase),pointer              :: db
    character(nameLen)                            :: geomName, graphType, nucData
    class(geometry), pointer                      :: geom
    type(outputFile)                              :: test_out
    class(baseMgNeutronMaterial), pointer         :: mat, mat2
    class(materialHandle), pointer                :: matPtr, matPtr2
!--> DNP 230918
    class(fissionMG), pointer                     :: fiss
    real(defFlt)                                  :: omega, timeXS, velInv
!<-- DNP 230918
!--> MK 240412 3D transient
    real(defFlt)                                  :: rodPartial, rodCellsInserted
    integer(shortInt)                             :: rodFull
!<-- MK 240412 3D transient
    character(100), parameter :: Here = 'init (dynamicRRPhysicsPackage_class_TI.f90)'

    call cpu_time(self % CPU_time_start)
    
    ! Load settings
    call dict % get( nucData, 'XSdata')
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    call dict % get(self % active, 'active')
!--> MK 230718
    call dict % get(self % nT, 'nT')
    call dict % getOrDefault(self % window, 'window', 200)
    call dict % getOrDefault(self % inactive, 'inactive', self % window)
    
    if (self % window > self % inactive) call fatalError(Here,&
            'Moving window must be smaller then number of inactive cycles!')
 
    ! Load number of timesteps printed for output and determine output interval
    call dict % get(self % nTOut, 'tOut')
    outTemp = real(self % nT / self % nTOut, defFlt)
    self % outInterval = NINT(outTemp)

    call dict % get(tstpReal, 'tstp')
    self % tstp = real(tstpReal, defFlt)
    call dict % getOrDefault(initReal, 'initVal', 1.0_defReal)
    self % initVal = real(initReal, defFlt)
    
    if (self % nT > 1) then
        self % nT = self % nT + 1 
    end if
    
!--> MK 240523
    call dict % get(precReal, 'sarePrec')
    self % sarePrec = real(precReal, defFlt)
    
    call dict % get(precReal, 'rmsPrec')
    self % rmsPrec = real(precReal, defFlt)
!<-- MK 240523
    
    ! Reduce number of ouput lines?
    call dict % getOrDefault(self % reduceLines, 'reduceLines', .false.)
    
    ! Use isotropic approximation for entire time derivative?
    call dict % getOrDefault(self % isotropic, 'isotropicDeriv', .false.)
!<-- MK 230718   
    
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
    call new_geometry(tempDict, geomName)
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
    allocate(self % scalarFluxTD(self % nCells * self % nG))
    allocate(self % prevFluxTD(self % nCells * self % nG))
    allocate(self % fluxScoresTD(self % nCells * self % nG, 2))
    allocate(self % sourceTD(self % nCells * self % nG))
    allocate(self % fluxHistTD(self % nCells * self % nG))
    allocate(self % prevFissionTD_1(self % nCells))
    allocate(self % prevFissionTD_2(self % nCells))
    allocate(self % fissionScoreTD(self % nCells))
    allocate(self % fissionSourceTD(self % nCells))
!--> MK 240517
    allocate(self % fissionRateScore(self % nCells, 2))
    allocate(self % fissionRate(self % nCells))
!<-- MK 240517

    allocate(self % scalarFlux(self % nCells * self % nG))
    allocate(self % prevFlux(self % nCells * self % nG))
    allocate(self % fluxScores(self % nCells * self % nG))
    allocate(self % source(self % nCells * self % nG))
    allocate(self % fluxHist(self % nCells * self % nG))
    allocate(self % prevFission_1(self % nCells))
    allocate(self % prevFission_2(self % nCells)) 
    allocate(self % fissionScore(self % nCells))
    allocate(self % fissionSource(self % nCells))

    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))
    
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
    self % nMat = mm_nMat()
!--> MK 230307
    allocate(self % sigmaT(self % nMat * self % nG * self % nT))
    allocate(self % nuSigmaF(self % nMat * self % nG * self % nT))
!--> MK 240517
    allocate(self % SigmaF(self % nMat * self % nG * self % nT))
!<-- MK 240517
    allocate(self % chi(self % nMat * self % nG))
    allocate(self % nu(self % nMat * self % nG))
    allocate(self % sigmaS(self % nMat * self % nG * self % nG * self % nT))
    allocate(self % vel(self % nMat * self % nG * self % nT))
!<-- MK 230307

!--> MK 230307
    ! Keep track of number of precursor groups
    ! Error if not uniform across all materials
    nP = -1
    fiss => null()
    do m = 1, self % nMat
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      do g = 1, self % nG
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, self % rand),defFlt)
!--> MK 240411
        self % nu(self % nG * (m - 1) + g) = real(mat % getNu(g, self % rand),defFlt)
!<-- MK 240411
        do t = 1, self % nT
            self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = real(mat % getVel(g, self % rand),defFlt)
            self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = real(mat % getTotalXS(g, self % rand),defFlt) 
!             self % nuSigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                 real(mat % getNuFissionXS(g, self % rand),defFlt) 
!--> MK 240517
            self % sigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
                real(mat % getFissionXS(g, self % rand),defFlt) 
!<-- MK 240517
            do g1 = 1, self % nG
                self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
                self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) 
            end do
        end do
      end do
!<-- MK 230307
      
!--> DNP 230918
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
!<-- DNP 230918

!--> MK 240228
    ! Read transient data           !TODO This is probably very ugly
    call variations % init(dict % getDictPtr('transient'))
    
    do i = 1, size(variations % xsArr)
        m = variations % matArr(i)
        matPtr  => self % mgData % getMaterial(m)
        mat     => baseMgNeutronMaterial_CptrCast(matPtr)
        
        tCnt = 0
        tStart = variations % timeStartArr(i) + 1 !To account for steady-state at beginning
        tEnd = variations % timeEndArr(i) + 1
        if (tEnd .le. tStart) tEnd = self % nT  !or flag error
        
        magnitude = variations % magnitudeArr(i)        !magnitude in terms of original XS
        variationType = variations % typeArr(i)
        xsName = variations % xsArr(i)
        
        select case(xsName)
            case('fission')
!--> MK 240517
                allocate (xsType(size(self % sigmaF)))        !Do they need to be deallocated
                xsType = self % sigmaF
!<-- MK 240517
            case('total')
                allocate (xsType(size(self % sigmaT)))
                xsType = self % sigmaT
            case('scatter')
                allocate (xsType(size(self % sigmaS)))
                xsType = self % sigmaS
        end select
        
        baseIdx0 = self % nG * self % nT * (m - 1)
        baseIdx1 = self % nG * self % nT * (m - 1) + (tStart - 2) * self % nG
        
        if (xsName == 'scatter') then
            baseIdx0 = self % nG * baseIdx0
            baseIdx1 = self % nG * baseIdx1
            
            allocate (origXS(self % nG * self % nG))
            allocate (startXS(self % nG * self % nG))
            origXS = xsType(baseIdx0 + 1 : baseIdx0 + self % nG * self % nG)
            startXS = xsType(baseIdx1 + 1 : baseIdx1 + self % nG * self % nG)
        else
            allocate (origXS(self % nG))
            allocate (startXS(self % nG))        
            origXS = xsType(baseIdx0 + 1 : baseIdx0 + self % nG) 
            startXS = xsType(baseIdx1 + 1 : baseIdx1 + self % nG)
        end if

        print *
        write(*, '(A,A,A,A,A,E16.9)') 'XS (', trim(xsName), ' ', trim(mm_matName(m)) , ') in Group 1 before variation:' &
            , xsType(baseIdx1 + 1)
    
        do t = tStart, tEnd
        
            baseIdx = self % nG * self % nT * (m - 1) + (t - 1) * self % nG
            
            if (xsName == 'scatter') baseIdx = self % nG * baseIdx

            select case(variationType)
            
                case('step')
                    do g = 1, self % nG
                        if (xsName == 'scatter') then
                            do g1 = 1, self % nG
                                xsType(baseIdx + self % nG * (g - 1) + g1) = startXS(self % nG * (g - 1) + g1) + &
                                    magnitude * origXS(self % nG * (g - 1) + g1)
                            end do                
                        else
                            xsType(baseIdx + g) = startXS(g) + magnitude * origXS(g)                       
                        end if
                    end do
                    
                case('ramp')
                    tCnt = tCnt + 1
                    timeXS = tCnt * self % tstp
                    do g = 1, self % nG
                        if (xsName == 'scatter') then
                            do g1 = 1, self % nG
                                xsType(baseIdx + self % nG * (g - 1) + g1) = startXS(self % nG * (g - 1) + g1) + & 
                                    magnitude * timeXS * origXS(self % nG * (g - 1) + g1)
                            end do                
                        else
                            xsType(baseIdx + g) = startXS(g) + magnitude * timeXS * origXS(g)
                        end if
                    end do
            
            end select
            
             if (t == tStart) write(*, '(A,A,A,A,A,F9.6,A,E16.9)') 'XS (', trim(xsName), ' ', trim(mm_matName(m)),&
                ') in Group 1 at beginning of variation (t = ', (tStart - 1)*self % tstp, ' s): ', xsType(baseIdx + 1)
             if (t == tEnd) write(*, '(A,A,A,A,A,F9.6,A,E16.9)') 'XS (', trim(xsName), ' ', trim(mm_matName(m)),&
                ') in Group 1 at end of variation (t = ', (tEnd - 1)*self % tstp, ' s): ', xsType(baseIdx + 1)        
        end do
        
        select case(xsName)
            case('fission')
!--> MK 240517
                self % sigmaF = xsType
!<-- MK 240517
            case('total')
                self % sigmaT = xsType
            case('scatter')
                self % sigmaS = xsType
        end select
                
        if(allocated(origXS)) deallocate(origXS)
        if(allocated(startXS)) deallocate(startXS)
        if(allocated(xsType)) deallocate(xsType)

    end do
    
    call variations % kill()
!<-- MK 240228  
    
! !--> MK 230718
!     do m = 1, self % nMat       !TODO This is only provisional
!       matPtr  => self % mgData % getMaterial(m)
!       mat     => baseMgNeutronMaterial_CptrCast(matPtr)
!       do g = 1, self % nG
!         do t = 2, self % nT
! !             self % nuSigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
! !             real(mat % getNuFissionXS(g, self % rand)*1.001,defFlt)
! !--> MK 240517
!             self % sigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!             real(mat % getFissionXS(g, self % rand)*1.001,defFlt)
! !<-- MK 240517
!         end do
!       end do
!     end do
! !<-- MK 230718

! !--> MK 230718 Moderator density change
!     omega = 0.8
!     m = 1       !TODO This is only provisional
!       matPtr  => self % mgData % getMaterial(m)
!       mat     => baseMgNeutronMaterial_CptrCast(matPtr)
!       do g = 1, self % nG
!        
!         do t = 2, self % nT
!             timeXS = (t-1) * self % tstp
!             if (timeXS <= 1) then
!                 self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat % getTotalXS(g, self % rand),defFlt) * (1 - (1 - omega) * timeXS)
!                 do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) &
!                         * (1 - (1 - omega) * timeXS)
!                 end do
!             end if
!             
!             if (timeXS > 1 .AND. timeXS <= 2) then
!                 self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat % getTotalXS(g, self % rand),defFlt) * (omega + (1 - omega) * (timeXS - 1))
!                 do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) &
!                         * (omega + (1 - omega) * (timeXS - 1))
!                 end do
!             end if
! 
!         end do
!       end do
! !<-- MK 230718

! !--> MK 240404 Control rod insertion !TODO This is only provisional
!     omega = 0.01
!     m = 2       ! rod inside
!     matPtr  => self % mgData % getMaterial(m)
!     mat     => baseMgNeutronMaterial_CptrCast(matPtr)
!     
!     m2 = 1      ! CR
!     matPtr2  => self % mgData % getMaterial(m2)
!     mat2     => baseMgNeutronMaterial_CptrCast(matPtr2)
!     
!     do g = 1, self % nG
!     
!         do t = 2, self % nT
!             timeXS = (t-1) * self % tstp
!             if (timeXS <= 1) then
!                 self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat % getTotalXS(g, self % rand),defFlt) + omega * &
!                     (real(mat2 % getTotalXS(g, self % rand) - mat % getTotalXS(g, self % rand),defFlt)) * timeXS
!                     
!                 velInv = real(ONE / mat % getVel(g, self % rand),defFlt) + omega * &
!                     (real((ONE / mat2 % getVel(g, self % rand)) - (ONE / mat % getVel(g, self % rand)),defFlt)) * timeXS    !TODO
!                     
!                 do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) &
!                         + omega * (real(mat2 % getScatterXS(g1, g, self % rand) - mat % getScatterXS(g1, g, self % rand), defFlt)) &
!                         * timeXS
!                 end do
!             end if
!             
!             if (timeXS > 1 .AND. timeXS <= 2) then
!                 self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat % getTotalXS(g, self % rand),defFlt) + omega * &
!                     (real(mat2 % getTotalXS(g, self % rand) - mat % getTotalXS(g, self % rand),defFlt)) * (2 - timeXS)
!                     
!                 velInv = real(ONE / mat % getVel(g, self % rand),defFlt) + omega * &        !TODO
!                     (real((ONE / mat2 % getVel(g, self % rand)) - (ONE / mat % getVel(g, self % rand)),defFlt)) * (2 - timeXS)
!                     
!                 do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) &
!                         + omega * (real(mat2 % getScatterXS(g1, g, self % rand) - mat % getScatterXS(g1, g, self % rand), defFlt)) &
!                         * (2 - timeXS)
!                 end do
!             end if
!             
!             self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = 1.0_defFlt / velInv !TODO
! 
!         end do
!     end do
! !<-- MK 240404

! !--> MK 240412 3D Control rod insertion !TODO This is only provisional
! !     omega = 0.01
! !     m = 2       ! rod inside
! !     matPtr  => self % mgData % getMaterial(m)
! !     mat     => baseMgNeutronMaterial_CptrCast(matPtr)
!     
!     m2 = 8      ! CR
!     matPtr2  => self % mgData % getMaterial(m2)
!     mat2     => baseMgNeutronMaterial_CptrCast(matPtr2)
!     
!     do g = 1, self % nG
!     
!         do t = 2, self % nT
!             timeXS = (t-1) * self % tstp
!             
!             ! Insertion
!             if (timeXS <= 2) then
!                 rodCellsInserted = 17 * timeXS
!                 rodFull = FLOOR(rodCellsInserted)
!                 rodPartial = rodCellsInserted - rodFull
!                 
! !                 if (g == 1) then
! !                     print *, timeXS
! !                     print *, rodCellsInserted
! !                     print *, rodFull
! !                     print *, rodPartial
! !                     print *
! !                 end if
!                 
!                 if (rodCellsInserted > 34) then
!                     print *, "Number of inserted rod cells: ", rodCellsInserted, " at t = ", timeXS
!                     rodPartial = 0.0_defFlt
!                 end if
!                 
!                 do m = 8 + 1, 8 + rodFull !matIdx for bank 1 starts at 8
!                   self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat2 % getTotalXS(g, self % rand),defFlt)
!                   self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = real(mat2 % getVel(g, self % rand),defFlt)
!                     
!                   do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat2 % getScatterXS(g1, g, self % rand),defFlt)
!                   end do
!                 end do
!                 
!                 m = 8 + rodFull + 1
!                 self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat % getTotalXS(g, self % rand),defFlt) + rodPartial * &
!                     (real(mat2 % getTotalXS(g, self % rand) - mat % getTotalXS(g, self % rand),defFlt)) * timeXS
!                     
!                 velInv = real(ONE / mat % getVel(g, self % rand),defFlt) + rodPartial * &
!                     (real((ONE / mat2 % getVel(g, self % rand)) - (ONE / mat % getVel(g, self % rand)),defFlt)) * timeXS
!                 self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = 1.0_defFlt / velInv
!                     
!                 do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) &
!                         + rodPartial * &
!                         (real(mat2 % getScatterXS(g1, g, self % rand) - mat % getScatterXS(g1, g, self % rand), defFlt)) * timeXS
!                 end do
! !                 
! !                 if (g == 1) then
! !                     print *, m
! !                     print *, self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g)
! !                     print *, rodPartial
! !                     print *
! !                 end if
!                 
!             end if
!             
!             ! Extraction
!             if (timeXS > 2 .AND. timeXS <= 4) then
!                 rodCellsInserted = -17 * timeXS + 34
!                 rodFull = FLOOR(rodCellsInserted)
!                 rodPartial = rodCellsInserted - rodFull
!                 
!                 if (rodCellsInserted < 0) then
!                     print *, "Number of inserted rod cells: ", rodCellsInserted, " at t = ", timeXS
!                     rodFull = 0
!                     rodPartial = 0.0_defFlt
!                 end if
!                 
!                 do m = 8 + 1, 8 + rodFull !matIdx for bank 1 starts at 8
!                   self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat2 % getTotalXS(g, self % rand),defFlt)
!                   self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = real(mat2 % getVel(g, self % rand),defFlt)
!                     
!                   do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat2 % getScatterXS(g1, g, self % rand),defFlt) 
!                   end do
!                 end do
!                 
!                 m = 8 + rodFull + 1
!                 self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
!                     real(mat % getTotalXS(g, self % rand),defFlt) + rodPartial * &
!                     (real(mat2 % getTotalXS(g, self % rand) - mat % getTotalXS(g, self % rand),defFlt)) * timeXS
!                     
!                 velInv = real(ONE / mat % getVel(g, self % rand),defFlt) + rodPartial * &
!                     (real((ONE / mat2 % getVel(g, self % rand)) - (ONE / mat % getVel(g, self % rand)),defFlt)) * timeXS
!                 self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = 1.0_defFlt / velInv
!                     
!                 do g1 = 1, self % nG
!                     self % sigmaS(self % nG * self % nG * self % nT * (m - 1) + (t - 1) * self % nG * self % nG + &
!                         self % nG * (g - 1) + g1) = real(mat % getScatterXS(g1, g, self % rand), defFlt) &
!                         + rodPartial * &
!                         (real(mat2 % getScatterXS(g1, g, self % rand) - mat % getScatterXS(g1, g, self % rand), defFlt)) * timeXS
!                 end do
!             end if
!             
! !             self % vel(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = 1.0_defFlt / velInv !TODO
! 
!         end do
!     end do
! !<-- MK 240412

!--> MK 240517
    do m = 1, self % nMat
      do g = 1, self % nG
        do t = 1, self % nT
            self % nuSigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
                self % nu(self % nG * (m - 1) + g) * self % sigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g)
        end do
      end do
    end do
!<-- MK 240517

!--> DNP 230918
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
            self % chiD(idx0 + 1 : idx0 + self % nP * self % nG) = 0.0_defFlt
            
            idx0 = self % nP * (m - 1)
            self % beta(idx0 + 1 : idx0 + self % nP) = 0.0_defFlt
            self % lambda(idx0 + 1 : idx0 + self % nP) = 1.0_defFlt
          end if
        else
          idx0 = self % nP * self % nG * (m - 1)
          self % chiD(idx0 + 1: idx0 + self % nP * self % nG) = 0.0_defFlt
          
          idx0 = self % nP * (m - 1)
          self % beta(idx0 + 1 : idx0 + self % nP) = 0.0_defFlt
          self % lambda(idx0 + 1 : idx0 + self % nP) = 1.0_defFlt
        end if
        fiss => null()

        ! Construct chi_prompt
        idx0 = self % nG * (m - 1)
        self % chiP(idx0 + 1: idx0 + self % nG) = &
                self % chi(idx0 + 1 : idx0 + self % nG) 
        do g = 1, self % nG
          do p = 1, self % nP
            self % chiP(self % nG * (m - 1) + g) = self % chiP(self % nG * (m - 1) + g) &
                    - self % beta(self % nP * (m - 1) + p) * &
                    self % chiD(self % nP * self % nG * (m - 1) + self % nP * (g - 1) + p)
          end do
          self % chiP(self % nG * (m - 1) + g) = self % chiP(self % nG * (m - 1) + g) / &
                  (1.0_defFlt - sum(self % beta(self % nP * (m - 1) + 1: self % nP * (m - 1) + self % nP)))
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
      self % lambda = 1.0_defFlt
    end if
    
    allocate(self % dnpScoresTD(self % nCells * self % nP))
    allocate(self % prevDnpTD(self % nCells * self % nP))
    allocate(self % dnpTD(self % nCells * self % nP))
    allocate(self % dnpScores(self % nCells * self % nP))
    allocate(self % prevDnp(self % nCells * self % nP))
    allocate(self % dnp(self % nCells * self % nP))
!<-- DNP 230918
!--> DNP 2nd 240311
    allocate(self % omega0(self % nMat * self % nP))
    allocate(self % omegaN(self % nMat * self % nP))
    allocate(self % omega_1(self % nMat * self % nP))
    allocate(self % omega_2(self % nMat * self % nP))

    if (nP /= -1) then
        do m = 1, self % nMat
            idx0 = self % nP * (m - 1)
            !$omp simd
            do p = 1, self % nP 
                kappa0 = exponential(self % lambda(idx0 + p) * self % tstp)     ! exponential(tau) = 1 - exp(-tau)
                kappa1 = 1 - (kappa0/(self % lambda(idx0 + p) * self % tstp))
                kappa2 = 1 - (2 * kappa1/(self % lambda(idx0 + p) * self % tstp))
                
                self % omega0(idx0 + p) = - (exponential(self % lambda(idx0 + p) * self % tstp) - 1)
                self % omegaN(idx0 + p) = real((kappa1 + kappa2) / 2, defFlt)
                self % omega_1(idx0 + p) = real(kappa0 - kappa2, defFlt)
                self % omega_2(idx0 + p) = real((kappa2 - kappa1) / 2, defFlt)
            end do 
        end do
    else
        self % omega0 = 0.0_defFlt
        self % omegaN = 0.0_defFlt
        self % omega_1 = 0.0_defFlt
        self % omega_2 = 0.0_defFlt
    end if
!<-- DNP 2nd 240311

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self

    call self % printSettings()
    call self % cycles() 
    call self % printResults()

  end subroutine run

!--> MK 230823       
  !!
  !! Perform cycles of The Random Ray Method.
  !!
  !! Randomly places the ray starting point and direction uniformly.
  !! Rays are tracked until they reach some specified termination length.
  !! During tracking, fluxes are attenuated (and adjusted according to BCs),
  !! scoring to fluxes and volume estimates when the ray has surpassed its 
  !! specified dead length.
  !!
  !! Inactive and active iterations occur, terminating subject either to 
  !! given criteria or when a fixed number of iterations has been passed.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  subroutine cycles(self)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defFlt)                                  :: hitRate, ONE_KEFF
    real(defReal), dimension(self % window)        :: entropyVec
    real(defReal)                                 :: elapsed_T, transport_T, scaleFac
    logical(defBool)                              :: stoppingCriterion, isActive
    integer(shortInt)                             :: i, itInac, itAct, it, t, nFiss
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    logical(defBool)                              :: convBool, entropyBool
    integer(shortInt)                             :: tCnt, outCnt
!--> MK 240517
    real(defFlt)                                  :: sare
    integer(shortInt)                             :: idx
    real(defReal)                                 :: N1, Nm1, stdSum
    real(defReal), save                           :: fluxMean, fluxStd
!<-- MK 240517
!--> MK 240520
    real(defFlt)                                  :: RMS, RMS0
    real(defFlt), save                            :: sumValue
    integer(shortInt)                             :: cIdx, cnt, popTemp
    real(defFlt), save                            :: sumMean1, sumMean2, err 
    integer(shortInt), save                       :: g, baseIdx
    real(defReal), dimension(self % nCells)       :: fluxNew, fluxOld
    !$omp threadprivate(sumMean1, sumMean2, err, pRNG, r, ints, fluxMean, fluxStd, sumValue, g, baseIdx)
!<-- MK 240520

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)
    
    ! Initialise scores and fluxes
    scaleFac          = ZERO   
    self % keff       = 1.0_defFlt     
    self % keffScore  = ZERO        
   
    self % fluxScoresTD = ZERO
    self % fissionRateScore = ZERO
    self % prevFluxTD = 1.0_defFlt  
    self % fluxHistTD = ZERO   
    self % fissionScoreTD = ZERO
    self % dnpScoresTD = ZERO  
    self % prevDnpTD = 0.0_defFlt
    self % prevFissionTD_1 = 0.0_defFlt
    self % prevFissionTD_2 = 0.0_defFlt
    
    self % fluxScores = ZERO
    self % prevFlux = 1.0_defFlt  
    self % fluxHist = ZERO    
    self % fissionScore = ZERO
    self % dnpScores = ZERO
    self % prevDnp = 0.0_defFlt
    self % prevFission_1 = 0.0_defFlt
    self % prevFission_2 = 0.0_defFlt
    
    ! Initialise counters and volumes
    tCnt = 0
    outCnt = 0
    self % totalIt = 0
    self % totalAct = 0
    self % totalInact = 0
    self % volumeTracks = ZERO  ! To score volumes over time steps
    self % volume       = ZERO  !
!--> MK 240523   
    self % AccVol = ZERO
    RMS = 0.0_defFlt
!<-- MK 240523   
    
    ! Time-stepping
    do t = 1, self % nT
    
!--> MK 240520
    cnt = 0
    fluxNew = 0.0_defFlt
    fluxOld = 0.0_defFlt
!<-- MK 240520
        
        tCnt = tCnt + 1     

        ! Initialise fluxes and sources
        self % scalarFluxTD = 0.0_defFlt
        self % scalarFlux   = 0.0_defFlt
        self % sourceTD     = 0.0_defFlt
        self % source       = 0.0_defFlt

        ! Initialise other results
        self % cellHit      = 0
        
        ! Initialise cell information
        self % cellFound = .false.
        self % cellPos = -INFINITY

        ! Stopping criterion is initially on flux convergence or number of convergence iterations.
        ! Will be replaced by RMS error in flux or number of scoring iterations afterwards.
        itInac = 0
        itAct  = 0
        isActive = .false.
        stoppingCriterion = .true.        
        convBool = .true.
        entropyVec = 0.0_defFlt
        sare = 1.0_defFlt
        
        ! Power iteration
        do while( stoppingCriterion )
            self % dnpTD = 0.0_defFlt
            self % dnp = 0.0_defFlt
            
            if (isActive) then
                itAct = itAct + 1
                self % totalAct = self % totalAct + 1
                self % totalInact = self % totalInact + 1
            else
                itInac = itInac + 1
            end if
            it = itInac + itAct
                    
            self % totalIt = self % totalIt + 1      ! For metrics, also volume accumlation MK 231011
                    
            ONE_KEFF = 1.0_defFlt / self % keff

            !$omp parallel do schedule(static)
            do i = 1, self % nCells
                call self % sourceUpdateKernel_TD(i, ONE_KEFF, t)
                call self % sourceUpdateKernel(i, ONE_KEFF, t)
            end do
            !$omp end parallel do
            
            ! Reset and start transport timer
            call timerReset(self % timerTransport)
            call timerStart(self % timerTransport)
            intersections = 0
            
!--> MK 240523
            ! Adjustment of ray population size over convergence 
            if (.not. isActive) then
!             if ((t > 1) .and. (.not. isActive)) then
!               self % pop  = INT(650* 0.1 * (1 + it/1000))
!               self % pop  = INT(650* ((0.9 * RMS)/(self % convPrec - 1) + (0.1 - 0.9/(self % convPrec - 1))))
!               popTemp  = INT(650* ((0.9/log(2E-4)) * log(RMS) + 0.1))
                if (it > 2+2) then
                    popTemp = INT(650* ((0.9/log(2E-4/RMS0)) * log(RMS/RMS0) + 0.1))
                    if (popTemp > self % pop) self % pop = popTemp
                else
                    self % pop = INT(650* 0.1)
                    RMS0 = RMS
                end if 
            end if
!<-- MK 240523
            
            !$omp parallel do schedule(runtime) reduction(+: intersections) 
            do i = 1, self % pop
            
                ! Set seed
                pRNG = self % rand
                call pRNG % stride(i)
                r % pRNG => pRNG 

                ! Set ray attributes
                call self % initialiseRay(r)
                
                ! Transport ray until termination criterion met
                if (self % isotropic) then
                    call self % transportSweep_isotropic(r,ints,t)
                else
                    call self % transportSweep(r,ints,t)
                end if
                intersections = intersections + ints

            end do
            !$omp end parallel do      
            
            call timerStop(self % timerTransport)

            ! Update RNG on master thread
            call self % rand % stride(self % pop + 1)
            
            ! Normalise flux estimate and combines with source
            call self % normaliseFluxAndVolume_TD(t)
            call self % normaliseFluxAndVolume()
            
            ! Calculate proportion of cells that were hit
            hitRate = real(sum(self % cellHit),defFlt) / self % nCells
            self % cellHit = 0
            
            ! Calculate keff
            call self % calculateKeff(t)

            ! Scale initial flux to specified value
            if (t == 1) then
                scaleFac = (self % initVal / sum(self % scalarFluxTD))
                self % scalarFluxTD = self % scalarFluxTD * real(scaleFac, defFlt)
            end if 

            ! Scale unperturbed flux for keff recalcualtion. This is always done, also in first time step.
            scaleFac = (self % initVal / sum(self % scalarFlux))
            self % scalarFlux = self % scalarFlux * real(scaleFac, defFlt)         

!--> MK 240517
            ! Source Average Relative Error convergence critertion for active cycles
            if (isActive) then
!--> MK 240517
                !$omp parallel do schedule(static)
                do i = 1, self % nCells
                    call self % calculateFissionRate(i, t)
                end do
                !$omp end parallel do
!<-- MK 240517

                call self % accumulateFluxAndKeffScores(t)
            
                if (itAct /= 1) then
                    Nm1 = 1.0_defReal/(itAct - 1)
                else
                    Nm1 = 1.0_defReal
                end if
                N1 = 1.0_defReal/itAct
                
                stdSum = 0.0_defFlt
                nFiss = 0
                
                !$omp parallel do schedule(static) reduction(+: stdSum, nFiss) 
                do idx = 1, self % nCells
                    fluxMean = self % fissionRateScore(idx, 1) * N1
                    fluxStd = self % fissionRateScore(idx, 2) * N1
                    fluxStd = Nm1 *(fluxStd - fluxMean * fluxMean) 
                    
                    if (fluxStd <= ZERO) then
                        fluxStd = ZERO
                    else
                        nFiss = nFiss + 1
                        fluxStd = sqrt(fluxStd)
                        fluxStd = fluxStd / fluxMean
                    end if      
                    
                    stdSum = stdSum + fluxStd
                    
                end do
                !$omp end parallel do
                
                sare = real(stdSum/nFiss,defFlt)
            end if
!<-- MK 240517

!--> MK 240520
            ! RMS convergence criterion for inactive cycles
!             if ((t > 1) .and. (.not. isActive)) then
            if (.not. isActive) then
            
                fluxOld  = fluxNew
                
                !$omp parallel do schedule(static)
                do cIdx = 1, self % nCells
                    
                    sumValue = 0.0_defFlt
                    
                    baseIdx = self % nG * (cIdx - 1)
                    do g = 1, self % nG
                        sumValue = sumValue + self % scalarFlux(baseIdx + g)
                    end do
                    
                    fluxNew(cIdx) = fluxNew(cIdx) + sumValue
                    
                end do 
                !$omp end parallel do

                if (it > 2) then
                    RMS = 0.0_defFlt
                    
                    !$omp parallel do schedule(static) reduction(+: RMS) 
                    do cIdx = 1, self % nCells
                        
                        sumMean1 = real(fluxOld(cIdx) / (it - 1) ,defFlt)
                        sumMean2 = real(fluxNew(cIdx) / it ,defFlt)
                        
                        err = (sumMean2 - sumMean1) / sumMean2 

                        RMS = RMS + err*err
                    end do  
                    !$omp end parallel do
                            
                    !Stop execution if RMS error is Nan
                    if (RMS /= RMS) then
                            print *, RMS, err, sumMean1, sumMean2, sumValue
                            call fatalError('Convergence','RMS is NaN')
                    end if
                    
                    RMS = sqrt(RMS/self % nCells)
                end if
            end if 
            
            stoppingCriterion = (sare > self % sarePrec .or. itAct < 10)
!<-- MK 240520

!--> MK 240520
            ! Check convergence of active cycles
!              if (.not. isActive .and. t == 1) call self % calculateEntropy(entropyVec, entropyBool)
! 
!             ! Check convergence of inactive cycles
!             if ( t == 1 .AND. convBool .AND. (it >= self % inactive)) then   
!                 isActive = .true.
!                 convBool = .false.
!             end if
            
!             if (t > 1 .AND. convBool .AND. RMS < 2E-4) then
            if (convBool .AND. RMS < self % rmsPrec) then
                cnt = cnt + 1
                
                if(cnt == 10) then      ! Convergence criterion should be met for 10 consecutive cycles
                    isActive = .true.
                    convBool = .false.
                end if
            else
                cnt = 0                 !TODO get rid of cnt?
            end if
!<-- MK 240530

            ! Set previous iteration flux to scalar flux and zero scalar flux
            call self % resetFluxes()

            ! Calculate times
            call timerStop(self % timerMain)
            elapsed_T = timerTime(self % timerMain)
            transport_T = timerTime(self % timerTransport)
            self % time_transport = self % time_transport + transport_T

            ! Display progress
            call printFishLineR(it)
            print *
            print *, 'Time step: ', numToChar(t), ' of ', numToChar(self % nT)
            if(isActive) then
                print *, 'Iteration: ', numToChar(itAct)
                if (.not. self % reduceLines) print *,'Active iterations'
                if (.not. self % reduceLines) print *, 'SARE: ', trim(numToChar(real(sare,defReal)))
            else
                print *, 'Iteration: ', numToChar(it)
                if (.not. self % reduceLines) print *,'Inactive iterations'
                if (.not. self % reduceLines) print *, 'RMS: ', trim(numToChar(real(RMS,defReal)))
            end if
            if (.not. self % reduceLines) print *, 'Ray population: ', self % pop
            if (.not. self % reduceLines) print *, 'Cell hit rate: ', trim(numToChar(real(hitRate,defReal)))
            if (.not. self % reduceLines) print *, 'keff: ', trim(numToChar(real(self % keff,defReal)))
            print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
            print *, 'Time per integration (ns): ', &
                    trim(numToChar(transport_T*10**9/(self % nG * intersections)))
        end do
        
        ! Finalise flux scores
        call self % finaliseFluxAndKeffScores(itAct, t) 
        
        ! Print transient output
        if (tCnt == self % outInterval .OR. t == 1 .OR. t == self % nT) then
            outCnt = outCnt + 1
            if (self % printFlux) call self % printTransientFlux(outCnt,t) 
            if (self % mapFission) call self % printFissionRate(outCnt,t) 
            tCnt = 0
        end if
        
        ! Update quantities for previous time steps needed for time derivatives
        call self % updatePreviousQuantities(t)
        
    end do
    
    ! Save number of time steps that are actually printed
    self % nTOut = outCnt

  end subroutine cycles
!<-- MK 230823 
  
!--> MK 240129
  !!
  !! Calculates the Shannon entropy over all cells and averages them over the halves of a travelling window
  !! Checks whether the mean values of either window half lie within one standard deviation of each other 
  !!
  subroutine calculateEntropy(self, entropyVec, entropyBool)
      class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
      real(defReal), dimension(self % window), intent(inout)   :: entropyVec
      logical(defBool), intent(out)                 :: entropyBool
      real(defReal)                                  :: entropyMean1, entropyMean2, entropyStd1, entropyStd2
      real(defFlt)                                  :: sumFlux
      real(defReal)                                  :: newEntropy
      integer(shortInt)                             :: cIdx, halfwindow
      real(defFlt), save                            :: P
      !$omp threadprivate(P)

!     integer(shortInt), intent(in)                         :: t
!       real(defFlt), dimension(self % nCells)   :: fissTemp
!       integer(shortInt), save                            :: matIdx, g, idx
!       real(defFlt), save                                  :: vol
!       
!       !$omp threadprivate(P, matIdx, g, idx, vol)
      
      halfwindow = INT(self % window * 0.5)

! !--> MK 240521
! if (t>1) then
!     fissTemp = 0.0_defFlt
!     
!     !$omp parallel do reduction(+: fissTemp)
!     do cIdx = 1, self % nCells
!         
!         ! Identify material
!         matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
!        vol = real(self % volume(cIdx),defFlt)
! 
!         if (vol < volume_tolerance) cycle
! 
!         if (matIdx == 1) then    
!           do g = 1, self % nG
!             idx = self % nG * (cIdx - 1) + g 
!             fissTemp(cIdx) = fissTemp(cIdx) + self % scalarFlux(idx)
!           end do
!         end if
!     end do
!     !$omp end parallel do   
! else
!     fissTemp = self % fissionRate
! end if
! !<-- MK 240521
      
!--> MK 240517      
      ! Calculate Shannon entropy
!       sumFlux = sum(fissTemp)

      sumFlux = sum(self % fissionRate)
      
      entropyVec(1:(self % window-1)) = entropyVec(2:self % window)
      newEntropy = ZERO
      
      !$omp parallel do schedule(static) reduction(+: newEntropy) 
      do cIdx = 1, self % nCells
        P = self % fissionRate(cIdx) / sumFlux
!         P = fissTemp(cIdx) / sumFlux

        if (P <= ZERO) then
            newEntropy = newEntropy
        else
            newEntropy = newEntropy - real(P * (log(P) / log(2.0)), defReal) 
        end if
      end do
      !$omp end parallel do
!<-- MK 240517
      
      ! Calculate average and standard deviation over either half of travelling window
      entropyVec(self % window) = newEntropy
      entropyMean1 = sum(entropyVec(1 : halfwindow)) / halfwindow
      entropyMean2 = sum(entropyVec((halfwindow+1) : self % window)) / halfwindow
      
      entropyStd1 = (1.0_defFlt / (halfwindow-1)) * sum((entropyVec(1 : halfwindow) - entropyMean1) &
            * (entropyVec(1 : halfwindow) - entropyMean1))
      entropyStd2 = (1.0_defFlt / (halfwindow-1)) * sum((entropyVec((halfwindow+1) : self % window) - entropyMean2) &
            * (entropyVec((halfwindow+1) : self % window) - entropyMean2))
            
      if (entropyStd1 <= ZERO) then
        entropyStd1 = ZERO
      else
        entropyStd1 = sqrt(entropyStd1)
      end if      
      
      if (entropyStd2 <= ZERO) then
        entropyStd2 = ZERO
      else
        entropyStd2 = sqrt(entropyStd2)
      end if
      
!       print *
!       print *, ABS(entropyMean1 - entropyMean2), entropyStd2, self % window
! !       print *, entropyMean1, entropyMean2
!       print *, newEntropy
! !       print *, ABS(entropyMean1 - entropyMean2)/entropyMean2
!       print *
      
      ! Check whether Shannon entropy lies within standard deviation
      entropyBool = (ABS(entropyMean1 - entropyMean2) < entropyStd2 .AND. ABS(entropyMean1 - entropyMean2) < entropyStd1)
      
  end subroutine calculateEntropy
!<-- MK 240129

  !!
  !! Initialises rays: samples initial position and direction,
  !! and performs the build operation
  !!
  subroutine initialiseRay(self, r)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, x
    integer(shortInt)                             :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (dynamicRRPhysicsPackage_class_TI.f90)'

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
  
!--> MK 240129
  !!
  !! Transient transport sweep
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep_isotropic(self, r, ints, t)
    class(dynamicRRPhysicsPackage_TI), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(shortInt), intent(in)                         :: t
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, event, matIdx0, baseIdx, baseIdx0
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt), dimension(self % nG)                    :: attenuate, deltaTD, fluxVec, fluxVecSteady, deltaSteady
    real(defFlt), pointer, dimension(:)                   :: totVec, scalarVec, sourceVec, scalarVecSteady, &
        sourceVecSteady, totVecSteady
    real(defFlt)                                          :: lenFlt
    real(defReal), dimension(3)                           :: r0, mu0

    ! Set initial angular flux to angle average of cell source (space)
    cIdx = r % coords % uniqueID

    baseIdx = (cIdx - 1) * self % nG
    sourceVec => self % sourceTD(baseIdx + 1 : baseIdx + self % nG)

    !$omp simd
    do g = 1, self % nG
        fluxVec(g) = sourceVec(g)
    end do

! --> MK240319 Recalculate
    sourceVecSteady => self % source(baseIdx + 1 : baseIdx + self % nG)
    !$omp simd
    do g = 1, self % nG
        fluxVecSteady(g) = sourceVecSteady(g)
    end do
! <-- MK240319 Recalculate
    
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
        
        ! Cache total cross section
        baseIdx0 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG
        totVec => self % sigmaT((baseIdx0 + 1) : (baseIdx0 + self % nG))
        
        baseIdx0 = (matIdx - 1) * self % nG * self % nT
        totVecSteady => self % sigmaT((baseIdx0 + 1) : (baseIdx0 + self % nG))
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
      lenFlt = real(length, defFlt)
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = r0 + HALF * length * mu0
        !$omp end critical
      end if
      
      ints = ints + 1
      
      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % sourceTD(baseIdx + 1 : baseIdx + self % nG)

      ! Calculate time-dependent flux
      !$omp simd
      do g = 1, self % nG
        attenuate(g) = exponential(totVec(g) * lenFlt)
        deltaTD(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
        fluxVec(g) = fluxVec(g) - deltaTD(g)
      end do
      
! --> MK240319 Recalculate    
      sourceVecSteady => self % source(baseIdx + 1 : baseIdx + self % nG) 

      !$omp simd            !I guess this could also be made to be combined with loop above
      do g = 1, self % nG
        attenuate(g) = exponential(totVecSteady(g) * lenFlt)
        deltaSteady(g) = (fluxVecSteady(g) - sourceVecSteady(g)) * attenuate(g)
        fluxVecSteady(g) = fluxVecSteady(g) - deltaSteady(g)
      end do
! <-- MK240319 Recalculate  

      ! Accumulate to scalar flux
      if (activeRay) then
        scalarVec => self % scalarFluxTD(baseIdx + 1 : baseIdx + self % nG)
        
        call OMP_set_lock(self % locks(cIdx))
        !$omp simd
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + deltaTD(g)
        enddo
        
! --> MK240319 Recalculate 
        scalarVecSteady => self % scalarFlux(baseIdx + 1 : baseIdx + self % nG)
        
        !call OMP_set_lock(self % locks(cIdx))
        !$omp simd
        do g = 1, self % nG
            scalarVecSteady(g) = scalarVecSteady(g) + deltaSteady(g)
        enddo
! <-- MK240319 Recalculate 

        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + lenFlt
        call OMP_unset_lock(self % locks(cIdx))

        if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1
      end if
    
      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, size(fluxVec)
            fluxVec(g) = 0.0_defFlt
        end do
! --> MK240319 Recalculate 
        !$omp simd
        do g = 1, size(fluxVecSteady)
            fluxVecSteady(g) = 0.0_defFlt
        end do 
! <-- MK240319 Recalculate 
      end if
    end do

  end subroutine transportSweep_isotropic
!<-- MK 240129

!--> MK 240129
  !!
  !! Transient transport sweep
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep(self, r, ints, t)
    class(dynamicRRPhysicsPackage_TI), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(shortInt), intent(in)                         :: t
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, event, matIdx0, baseIdx, baseIdx0
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, flux_avg, flux_avgprv, denom, deltaTD, fluxVec, &
        fluxVecSteady, deltaSteady
    real(defFlt), pointer, dimension(:)                   :: totVec, velVec, scalarVec, sourceVec, scalarVecSteady, &
        sourceVecSteady, totVecSteady, velVecSteady
    real(defFlt)                                          :: lenFlt
    real(defReal), dimension(3)                           :: r0, mu0

    ! Set initial angular flux to angle average of cell source (space)
    cIdx = r % coords % uniqueID

    baseIdx = (cIdx - 1) * self % nG
    sourceVec => self % sourceTD(baseIdx + 1 : baseIdx + self % nG)

    !$omp simd
    do g = 1, self % nG
        fluxVec(g) = sourceVec(g)
    end do

! --> MK240319 Recalculate
    sourceVecSteady => self % source(baseIdx + 1 : baseIdx + self % nG)
    !$omp simd
    do g = 1, self % nG
        fluxVecSteady(g) = sourceVecSteady(g)
    end do
! <-- MK240319 Recalculate
    
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
        
        ! Cache total cross section and neutron velocity
!         baseIdx0 = (matIdx - 1) * self % nG       !TODO
!         velVec => self % vel((baseIdx0 + 1) : (baseIdx0 + self % nG))
        
        baseIdx0 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG
        totVec => self % sigmaT((baseIdx0 + 1) : (baseIdx0 + self % nG))
        velVec => self % vel((baseIdx0 + 1) : (baseIdx0 + self % nG))
        
        baseIdx0 = (matIdx - 1) * self % nG * self % nT
        totVecSteady => self % sigmaT((baseIdx0 + 1) : (baseIdx0 + self % nG))
        velVecSteady => self % vel((baseIdx0 + 1) : (baseIdx0 + self % nG))
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
      lenFlt = real(length, defFlt)
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .true.
        self % cellPos(cIdx,:) = r0 + HALF * length * mu0
        !$omp end critical
      end if
      
      ints = ints + 1
      
      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % sourceTD(baseIdx + 1 : baseIdx + self % nG)
      scalarVec => self % fluxHistTD(baseIdx + 1 : baseIdx + self % nG)

      ! Isotropic approximation
      !$omp simd
      do g = 1, self % nG
        flux_avgprv(g) = real(scalarVec(g) * ONE_FOUR_PI, defFlt)
      end do
      
      if (t==1) then
        ! Calculate steady-state flux
        !$omp simd
        do g = 1, self % nG
            attenuate(g) = exponential(totVec(g) * lenFlt)
            deltaTD(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
            fluxVec(g) = fluxVec(g) - deltaTD(g)
        end do
      else
        ! Calculate time-dependent flux
        !$omp simd                                                    
        do g = 1, self % nG
            attenuate(g) = exponential_TD(totVec(g) * lenFlt)
            denom(g) = 1 / (totVec(g) * velVec(g) * self % tstp)
            flux_avg(g) = (attenuate(g) * (fluxVec(g) - sourceVec(g) - flux_avgprv(g) * denom(g)) &
                + sourceVec(g) + flux_avgprv(g) * denom(g)) / (1 + denom(g) * (1 - attenuate(g)))                         
            delta(g) = (fluxVec(g) - sourceVec(g) + (flux_avg(g) - flux_avgprv(g)) * denom(g)) * attenuate(g) * totVec(g) * lenFlt     ! * totVec(g) * lenFlt to account for changed exponetial function
            fluxVec(g) = fluxVec(g) - delta(g)
            deltaTD(g) = delta(g) - lenFlt * (flux_avg(g) - flux_avgprv(g))/(velVec(g) * self % tstp)
        end do
      end if
      
! --> MK240319 Recalculate    
      sourceVecSteady => self % source(baseIdx + 1 : baseIdx + self % nG) 
      scalarVecSteady => self % fluxHist(baseIdx + 1 : baseIdx + self % nG)
     
      ! Isotropic approximation
      !$omp simd
      do g = 1, self % nG
        flux_avgprv(g) = real(scalarVecSteady(g) * ONE_FOUR_PI, defFlt)
      end do
      
      if (t==1) then
        ! Calculate steady-state flux
        !$omp simd
        do g = 1, self % nG
            attenuate(g) = exponential(totVecSteady(g) * lenFlt)
            deltaSteady(g) = (fluxVecSteady(g) - sourceVecSteady(g)) * attenuate(g)
            fluxVecSteady(g) = fluxVecSteady(g) - deltaSteady(g)
        end do
      else
        !Calculate time-dependent flux
        !$omp simd                                                    
        do g = 1, self % nG
            attenuate(g) = exponential_TD(totVecSteady(g) * lenFlt)
            denom(g) = 1 / (totVecSteady(g) * velVecSteady(g) * self % tstp)
            flux_avg(g) = (attenuate(g) * (fluxVecSteady(g) - sourceVecSteady(g) - flux_avgprv(g) * denom(g)) &
                + sourceVecSteady(g) + flux_avgprv(g) * denom(g)) / (1 + denom(g) * (1 - attenuate(g)))                         
            delta(g) = (fluxVecSteady(g) - sourceVecSteady(g) + (flux_avg(g) - flux_avgprv(g)) * denom(g)) * attenuate(g) &
                * totVecSteady(g) * lenFlt
            fluxVecSteady(g) = fluxVecSteady(g) - delta(g)
            deltaSteady(g) = delta(g) - lenFlt * (flux_avg(g) - flux_avgprv(g))/(velVecSteady(g) * self % tstp)
        end do
      end if
! <-- MK240319 Recalculate  

      ! Accumulate to scalar flux
      if (activeRay) then
        scalarVec => self % scalarFluxTD(baseIdx + 1 : baseIdx + self % nG)
        
        call OMP_set_lock(self % locks(cIdx))
        !$omp simd
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + deltaTD(g)
        enddo
        
! --> MK240319 Recalculate 
        scalarVecSteady => self % scalarFlux(baseIdx + 1 : baseIdx + self % nG)
        
        !$omp simd
        do g = 1, self % nG
            scalarVecSteady(g) = scalarVecSteady(g) + deltaSteady(g)
        enddo
! <-- MK240319 Recalculate 

        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + lenFlt
        call OMP_unset_lock(self % locks(cIdx))

        if (self % cellHit(cIdx) == 0) self % cellHit(cIdx) = 1
      end if
    
      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, size(fluxVec)
            fluxVec(g) = 0.0_defFlt
        end do
! --> MK240319 Recalculate 
        !$omp simd
        do g = 1, size(fluxVecSteady)
            fluxVecSteady(g) = 0.0_defFlt
        end do 
! <-- MK240319 Recalculate 
      end if
    end do

  end subroutine transportSweep
!<-- MK 240129

!--> MK 230718
  !!
  !! Normalise transient flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume_TD(self, t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: t
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    integer(shortInt), save                       :: g, idx, matIdx
    integer(shortInt)                             :: cIdx
    real(defFlt), save                            :: total, vol
    !$omp threadprivate(total, vol, idx, g, matIdx)

!--> MK 240523
!     norm = real(ONE / self % lengthPerIt, defFlt)   
!     normVol = ONE / (self % lengthPerIt * self % totalIt)        ! To score volumes over time steps
!     
    norm = real(ONE / ((self % termination - self % dead) * self % pop), defFlt)    !TODO I guess lengthPerIt could suimply be updated somewhere as well?
    self % AccVol = self % AccVol + ((self % termination - self % dead) * self % pop)
    normVol = ONE / self % AccVol
!<-- MK 240523

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = real(self % volume(cIdx),defFlt)

      do g = 1, self % nG
        total = self % sigmaT((matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG + g)
        idx = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
          self % scalarFluxTD(idx) = real(self % scalarFluxTD(idx) * FOUR_PI, defFlt) * norm  / (total * vol)
        end if
        
        self % scalarFluxTD(idx) = self % scalarFluxTD(idx) + real(self % sourceTD(idx) * FOUR_PI, defFlt)
        
        if (self % scalarFluxTD(idx) < 0) self % scalarFluxTD(idx) = 0.0_defFlt
      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume_TD
!<-- MK 230718

  !!
  !! Normalise unperturbed flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    real(defFlt)                                  :: norm
    real(defFlt), save                            :: total, vol
    integer(shortInt), save                       :: g, matIdx, idx
    integer(shortInt)                             :: cIdx
    !$omp threadprivate(total, vol, idx, g, matIdx)

!--> MK 240523
!     norm = real(ONE / self % lengthPerIt, defFlt)
    norm = real(ONE / ((self % termination - self % dead) * self % pop), defFlt)
!<-- MK 240523

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Read volume
      vol = real(self % volume(cIdx),defFlt)
      
      do g = 1, self % nG

        total = self % sigmaT((matIdx - 1) * self % nG * self % nT + g)
        idx   = self % nG * (cIdx - 1) + g

        if (vol > volume_tolerance) then
           self % scalarFlux(idx) = real(self % scalarFlux(idx) * FOUR_PI, defFlt) * norm / ( total * vol)
        end if

        self % scalarFlux(idx) =  self % scalarflux(idx) + real(self % source(idx) * FOUR_PI, defFlt)
        if (self % scalarFlux(idx) < 0) self % scalarFlux(idx) = 0.0_defFlt     ! Fudged

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume

!--> MK 230822   
  !!
  !! Kernel to update transient sources given a cell index
  !!
  subroutine sourceUpdateKernel_TD(self, cIdx, ONE_KEFF, t)
    class(dynamicRRPhysicsPackage_TI), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    integer(shortInt), intent(in)                         :: t
    real(defFlt)                                          :: scatter, fission, ONE_MIN_B
    real(defFlt), dimension(:), pointer                   :: nuFission, total, scatterXS, velocity
    real(defFlt), dimension(:), pointer                   :: beta, lambda, chiP, chiD, chiDVec, dnpPrevVec, dnpVec
!--> DNP 2nd 240311
    real(defFlt), dimension(:), pointer                   :: omega0, omegaN, omega_1, omega_2
    real(defFlt)                                          :: omegaAcc, dnpAcc, dnpAcc_1, dnpAcc_2, sourceDnp, sourceDnpPrev, &
        fissionVec_1, fissionVec_2
!<-- DNP 2nd 240311
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, matIdx1, matIdx2, dnpIdx, p
    integer(shortInt)                                     :: baseIdx0, baseIdx1
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec
    real(defFlt), dimension(self % nG)                    :: fluxTimeDeriv

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= VOID_MAT) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % sourceTD(idx) = 0.0_defFlt
      end do
      return
    end if
    
    ! Define indeces
    matIdx1 = (matIdx - 1) * self % nG 
    baseIdx0 = matIdx1 * self % nP
    baseIdx1 = (matIdx - 1) * self % nP
    matIdx2 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG
    
    ! Obtain data for DNPs
    beta => self % beta((baseIdx1 + 1) : (baseIdx1 + self % nP))
    ONE_MIN_B = 1.0_defFlt - sum(beta)
    lambda => self % lambda((baseIdx1 + 1) : (baseIdx1 + self % nP))

    omega0 => self % omega0((baseIdx1 + 1) : (baseIdx1 + self % nP))
    omega_1 => self % omega_1((baseIdx1 + 1) : (baseIdx1 + self % nP))
    omega_2 => self % omega_2((baseIdx1 + 1) : (baseIdx1 + self % nP))
    omegaN => self % omegaN((baseIdx1 + 1) : (baseIdx1 + self % nP))

    dnpIdx = self % nP * (cIdx - 1)
    dnpPrevVec => self % prevDnpTD((dnpIdx+1) : (dnpIdx+self % nP))
    dnpVec => self % dnpTD((dnpIdx+1) : (dnpIdx+self % nP))
    
    ! Obtain XSs and other data
    chiP => self % chiP((matIdx1 + 1) : (matIdx1 + self % nG))
    chiD => self % chiD((baseIdx0 + 1) : (baseIdx0 + self % nG * self % nP))
!     velocity => self % vel((matIdx1 + 1) : (matIdx1 + self % nG))     !TODO
        
    velocity => self % vel((matIdx2 + 1) : (matIdx2 + self % nG))
    total => self % sigmaT((matIdx2 + 1) : (matIdx2 + self % nG))                        
    scatterXS => self % sigmaS((matIdx2 * self % nG + 1) : (matIdx2 * self % nG + self % nG*self % nG))
    nuFission => self % nuSigmaF((matIdx2 + 1) : (matIdx2 + self % nG))

    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % prevFluxTD((baseIdx+1) : (baseIdx+self % nG))
    
    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)                                 
    do gIn = 1, self % nG                                         
        fission = fission + fluxVec(gIn) * nuFission(gIn) 
    end do
    
    ! Calculate DNPs
    if (t==1) then
        !$omp simd
        do p = 1, self % nP 
            dnpVec(p) = (beta(p) / lambda(p)) * ONE_KEFF * fission
        end do
        
        fissionVec_1 = ONE_KEFF * fission
        fissionVec_2 = ONE_KEFF * fission
        dnpPrevVec => self % dnpTD((dnpIdx+1) : (dnpIdx+self % nP))
    else
        fissionVec_1 = self % prevFissionTD_1(cIdx)
        fissionVec_2 = self % prevFissionTD_2(cIdx)
        
        !$omp simd
        do p = 1, self % nP 
            dnpVec(p) = omega0(p) * dnpPrevVec(p) + (beta(p) / lambda(p)) * (ONE_KEFF * fission * omegaN(p) + &
                fissionVec_1 * omega_1(p) + fissionVec_2 * omega_2(p))
        end do
    end if       

    ! Calculate time derivative source term
    if (self % isotropic) then
        if (t==1) then       
            !$omp simd
            do gIn = 1, self % nG                                         
                fluxTimeDeriv(gIn) = 0.0_defFlt 
            end do
        else
            !$omp simd
            do gIn = 1, self % nG                                         
                fluxTimeDeriv(gIn) = (fluxVec(gIn) - real(self % fluxHistTD(baseIdx + gIn), defFlt)) / (velocity(gIn) * self % tstp) 
            end do
        end if
    else
        !$omp simd
        do gIn = 1, self % nG                                         
            fluxTimeDeriv(gIn) = 0.0_defFlt 
        end do
    end if

    do g = 1, self % nG

        ! Calculate DNP source
        omegaAcc = 0.0_defFlt
        dnpAcc = 0.0_defFlt
        dnpAcc_1 = 0.0_defFlt
        dnpAcc_2 = 0.0_defFlt
        baseIdx0 = self % nP * (g - 1)
        chiDVec => chiD((baseIdx0 + 1):(baseIdx0 + self % nP))
        
        !$omp simd reduction(+:omegaAcc,dnpAcc,dnpAcc_1,dnpAcc_2)  
        do p = 1, self % nP 
            omegaAcc = omegaAcc + chiDVec(p) * beta(p) * real(omegaN(p), defFlt)
            dnpAcc = dnpAcc + chiDVec(p) * lambda(p) * real(omega0(p), defFlt) * dnpPrevVec(p)
            dnpAcc_1 = dnpAcc_1 + chiDVec(p) * beta(p) * real(omega_1(p), defFlt)
            dnpAcc_2 = dnpAcc_2 + chiDVec(p) * beta(p) * real(omega_2(p), defFlt)
        end do
        
        sourceDnpPrev = dnpAcc + fissionVec_1 * dnpAcc_1 + fissionVec_2 * dnpAcc_2  
        sourceDnp = omegaAcc * fission * ONE_KEFF + sourceDnpPrev

        ! Calculate scattering source
        baseIdx0 = self % nG * (g - 1)
        scatterVec => scatterXS((baseIdx0 + 1) : (baseIdx0 + self % nG))    

        scatter = 0.0_defFlt
        
        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter)
        do gIn = 1, self % nG
            scatter = scatter + fluxVec(gIn) * scatterVec(gIn)  
        end do

        ! Output index
        idx = baseIdx + g

        ! Add all contributions
        self % sourceTD(idx) = ONE_MIN_B * chiP(g) * fission * ONE_KEFF + scatter + sourceDnp - fluxTimeDeriv(g)
        self % sourceTD(idx) = self % sourceTD(idx) / (real(total(g) * FOUR_PI, defFlt))
    end do
    
    ! Save fission source for following time steps
    self % fissionSourceTD(cIdx) = ONE_KEFF * fission

  end subroutine sourceUpdateKernel_TD
!<-- MK 230822   

!--> MK 230822   
  !!
  !! Kernel to update unperturbed sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF, t)
    class(dynamicRRPhysicsPackage_TI), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    integer(shortInt), intent(in)                         :: t
    real(defFlt)                                          :: scatter, fission, ONE_MIN_B
    real(defFlt), dimension(:), pointer                   :: nuFission, total, scatterXS, velocity
    real(defFlt), dimension(:), pointer                   :: beta, lambda, chiP, chiD, chiDVec, dnpPrevVec, dnpVec
    real(defFlt), dimension(:), pointer                   :: omega0, omegaN, omega_1, omega_2
    real(defFlt)                                          :: omegaAcc, dnpAcc, dnpAcc_1, dnpAcc_2, sourceDnp, sourceDnpPrev, &
        fissionVec_1, fissionVec_2
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, matIdx1, matIdx2, dnpIdx, p
    integer(shortInt)                                     :: baseIdx0, baseIdx1
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec
    real(defFlt), dimension(self % nG)                    :: fluxTimeDeriv

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx >= VOID_MAT) then
      baseIdx = self % ng * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
      end do
      return
    end if
    
    ! Define indeces
    matIdx1 = (matIdx - 1) * self % nG 
    baseIdx0 = matIdx1 * self % nP
    baseIdx1 = (matIdx - 1) * self % nP
    matIdx2 = (matIdx - 1) * self % nG * self % nT
    
    ! Obtain data for DNPs
    beta => self % beta((baseIdx1 + 1) : (baseIdx1 + self % nP))
    ONE_MIN_B = 1.0_defFlt - sum(beta)
    lambda => self % lambda((baseIdx1 + 1) : (baseIdx1 + self % nP))

    omega0 => self % omega0((baseIdx1 + 1) : (baseIdx1 + self % nP))
    omega_1 => self % omega_1((baseIdx1 + 1) : (baseIdx1 + self % nP))
    omega_2 => self % omega_2((baseIdx1 + 1) : (baseIdx1 + self % nP))
    omegaN => self % omegaN((baseIdx1 + 1) : (baseIdx1 + self % nP))
    
    dnpIdx = self % nP * (cIdx - 1)
    dnpPrevVec => self % prevDnp((dnpIdx+1) : (dnpIdx+self % nP))
    dnpVec => self % dnp((dnpIdx+1) : (dnpIdx+self % nP))

    ! Obtain XSs and other data
    chiP => self % chiP((matIdx1 + 1) : (matIdx1 + self % nG))
    chiD => self % chiD((baseIdx0 + 1) : (baseIdx0 + self % nG * self % nP))
!     velocity => self % vel((matIdx1 + 1) : (matIdx1 + self % nG)) !TODO
 
    velocity => self % vel((matIdx2 + 1) : (matIdx2 + self % nG))
    total => self % sigmaT((matIdx2 + 1) : (matIdx2 + self % nG))                        
    scatterXS => self % sigmaS((matIdx2 * self % nG + 1) : (matIdx2 * self % nG + self % nG*self % nG))
    nuFission => self % nuSigmaF((matIdx2 + 1) : (matIdx2 + self % nG)) 
    
    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx+1) : (baseIdx+self % nG))
    
    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)                                 
    do gIn = 1, self % nG                                         
        fission = fission + fluxVec(gIn) * nuFission(gIn) 
    end do
    
    ! Calculate DNPs
    if (t==1) then      
        !$omp simd
        do p = 1, self % nP 
            dnpVec(p) = (beta(p) / lambda(p)) * ONE_KEFF * fission
        end do
        
        fissionVec_1 = ONE_KEFF * fission
        fissionVec_2 = ONE_KEFF * fission
        dnpPrevVec => self % dnp((dnpIdx+1) : (dnpIdx+self % nP))
    else
        fissionVec_1 = self % prevFission_1(cIdx)
        fissionVec_2 = self % prevFission_2(cIdx)
        
        !$omp simd
        do p = 1, self % nP 
            dnpVec(p) = omega0(p) * dnpPrevVec(p) + (beta(p) / lambda(p)) * (ONE_KEFF * fission * omegaN(p) + &
                fissionVec_1 * omega_1(p) + fissionVec_2 * omega_2(p))
        end do
    end if   

    ! Calculate time derivative source term
    if (self % isotropic) then
        if (t==1) then       
            !$omp simd
            do gIn = 1, self % nG                                         
                fluxTimeDeriv(gIn) = 0.0_defFlt 
            end do
        else
            !$omp simd
            do gIn = 1, self % nG                                         
                fluxTimeDeriv(gIn) = (fluxVec(gIn) - real(self % fluxHist(baseIdx + gIn), defFlt)) / (velocity(gIn) * self % tstp) 
            end do
        end if
    else
        !$omp simd
        do gIn = 1, self % nG                                         
            fluxTimeDeriv(gIn) = 0.0_defFlt 
        end do
    end if

    do g = 1, self % nG
        ! Calculate DNP source
        omegaAcc = 0.0_defFlt
        dnpAcc = 0.0_defFlt
        dnpAcc_1 = 0.0_defFlt
        dnpAcc_2 = 0.0_defFlt
        baseIdx0 = self % nP * (g - 1)
        chiDVec => chiD((baseIdx0 + 1):(baseIdx0 + self % nP))
        
        !$omp simd reduction(+:omegaAcc,dnpAcc,dnpAcc_1,dnpAcc_2)  
        do p = 1, self % nP 
            omegaAcc = omegaAcc + chiDVec(p) * beta(p) * real(omegaN(p), defFlt)
            dnpAcc = dnpAcc + chiDVec(p) * lambda(p) * real(omega0(p), defFlt) * dnpPrevVec(p)
            dnpAcc_1 = dnpAcc_1 + chiDVec(p) * beta(p) * real(omega_1(p), defFlt)
            dnpAcc_2 = dnpAcc_2 + chiDVec(p) * beta(p) * real(omega_2(p), defFlt)
        end do
        
        sourceDnpPrev = dnpAcc + fissionVec_1 * dnpAcc_1 + fissionVec_2 * dnpAcc_2  
        sourceDnp = omegaAcc * fission * ONE_KEFF + sourceDnpPrev

        baseIdx0 = self % nG * (g - 1)
        scatterVec => scatterXS((baseIdx0 + 1) : (baseIdx0 + self % nG))    

        ! Calculate scattering source
        scatter = 0.0_defFlt
        
        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter)
        do gIn = 1, self % nG
            scatter = scatter + fluxVec(gIn) * scatterVec(gIn)  
        end do

        ! Output index
        idx = baseIdx + g

        ! Add all contributions
        self % source(idx) = ONE_MIN_B * chiP(g) * fission * ONE_KEFF + scatter + sourceDnp - fluxTimeDeriv(g)
        self % source(idx) = self % source(idx) / (real(total(g) * FOUR_PI, defFlt))
    end do
    
    ! Save fission source for following time steps
    self % fissionSource(cIdx) = ONE_KEFF * fission

  end subroutine sourceUpdateKernel
!<-- MK 230822   

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self, t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: t
    real(defFlt)                                  :: fissionRate, prevFissionRate
    real(defFlt), save                            :: fissLocal, prevFissLocal, vol
    integer(shortInt), save                       :: matIdx, g, idx, mIdx
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, fissLocal, prevFissLocal, matIdx, g, idx, mIdx, vol)

    fissionRate     = 0.0_defFlt
    prevFissionRate = 0.0_defFlt
    !$omp parallel do schedule(static) reduction(+: fissionRate, prevFissionRate)
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      if (matIdx >= VOID_MAT) cycle
      
      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)
      if (.not. mat % isFissile()) cycle

      vol = real(self % volume(cIdx), defFlt)

      if (vol <= volume_tolerance) cycle

      fissLocal = 0.0_defFlt
      prevFissLocal = 0.0_defFlt

      mIdx = (matIdx - 1) * self % nG * self % nT
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
  !! Sets prevFlux to scalarFluxTD and zeros scalarFluxTD
  !!
  subroutine resetFluxes(self)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt)                                :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFluxTD)
        ! Reset transient fluxes
        self % prevFluxTD(idx) = self % scalarFluxTD(idx)
        self % scalarFluxTD(idx) = 0.0_defFlt
        
        ! Reset steady-state fluxes
        self % prevFlux(idx) = self % scalarFlux(idx)
        self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do    

  end subroutine resetFluxes
  
!--> MK 240517   
  !!
  !! Calculate fission rate
  !!
  subroutine calculateFissionRate(self, cIdx, t)
    class(dynamicRRPhysicsPackage_TI), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    integer(shortInt), intent(in)                         :: t
    real(defFlt)                                          :: fRateTemp
    real(defFlt), dimension(:), pointer                   :: fissionXS
    integer(shortInt)                                     :: matIdx, gIn, baseIdx, matIdx2
    real(defFlt), pointer, dimension(:)                   :: fluxVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Define indeces
    matIdx2 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG
    fissionXS => self % sigmaF((matIdx2 + 1) : (matIdx2 + self % nG)) 

    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % scalarFluxTD((baseIdx+1) : (baseIdx+self % nG))

    fRateTemp = 0.0_defFlt
    !$omp simd reduction(+:fRateTemp)                                 
    do gIn = 1, self % nG                                         
        fRateTemp = fRateTemp + fluxVec(gIn) * fissionXS(gIn)
    end do
 
    self % fissionRate(cIdx) = fRateTemp * real(self % volume(cIdx),defFlt)

  end subroutine calculateFissionRate
!<-- MK 240517   

  
!--> MK 240402  
  !!
  !! Updates quantities of previous timesteps that are needed for time derivatives
  !!
  subroutine updatePreviousQuantities(self, t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: t
    integer(shortInt)                             :: cIdx, pIdx, g
    integer(shortInt), save                       :: idx0!, g
    !$omp threadprivate(idx0)
!     !$omp threadprivate(idx0, g)

    ! Set values for quantities of previous time steps when entering transient calculation
    if (t == 1) then
        ! Fission rate in previous time step
        self % prevFissionTD_1 = real(self % fissionScoreTD,defFlt) !--> DNP 2nd 240311
        self % prevFission_1(:) = real(self % fissionScore,defFlt)   
    end if
           
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
        idx0 = (cIdx - 1) * self % nG
        
        ! Update previous scalar flux, also to initialise source for next time step
        !$omp simd
        do g = 1, self % nG
        
            self % fluxHistTD(idx0 + g) = real(self % fluxScoresTD(idx0 + g, 1),defFlt)
            self % prevFluxTD(idx0 + g) = real(self % fluxScoresTD(idx0 + g, 1), defFlt)
            self % fluxScoresTD(idx0 + g, :) = ZERO
            
            self % fluxHist(idx0 + g) = real(self % fluxScores(idx0 + g),defFlt)
            self % prevFlux(idx0 + g) = real(self % fluxScores(idx0 + g), defFlt) 
            self % fluxScores(idx0 + g) = 0! ZERO
        end do

        ! Update fission sources of previous time steps
        self % prevFissionTD_2(cIdx) = self % prevFissionTD_1(cIdx)
!--> MK 240517
         self % prevFissionTD_1(cIdx) = real(self % fissionScoreTD(cIdx),defFlt)
         self % fissionScoreTD(cIdx) = ZERO
         self % fissionRateScore(cIdx, 1) = ZERO       !TODO necessary?
         self % fissionRateScore(cIdx, 2) = ZERO       !TODO necessary?
!<-- MK 240517

        self % prevFission_2(cIdx) = self % prevFission_1(cIdx)
        self % prevFission_1(cIdx) = real(self % fissionScore(cIdx),defFlt)
        self % fissionScore(cIdx) = ZERO
        
    end do
    !$omp end parallel do
    
    ! Update previous precursors
    !$omp parallel do schedule(static)
    do pIdx = 1, self % nP * self % nCells
        self % prevDnpTD(pIdx) = real(self % dnpScoresTD(pIdx),defFlt)
        self % dnpScoresTD(pIdx) = ZERO

        self % prevDnp(pIdx) = real(self % dnpScores(pIdx),defFlt)
        self % dnpScores(pIdx) = ZERO
    end do
    !$omp end parallel do

  end subroutine updatePreviousQuantities
!<-- MK 240402   

!--> MK 230822   
  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self, t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: t
!--> MK 240517
!     real(defReal), save                           :: flux
    real(defReal), save                           :: flux, source
!<-- MK 240517
    integer(shortInt)                             :: idx
    !$omp threadprivate(flux, source)

    !$omp parallel do schedule(static)
    do idx = 1, self % nCells * self % nG
        ! Transient scalar flux scores 
        flux = real(self % scalarFluxTD(idx), defReal)        
        self % fluxScoresTD(idx,1) = self % fluxScoresTD(idx, 1) + flux
        self % fluxScoresTD(idx,2) = self % fluxScoresTD(idx, 2) + flux*flux
        
        ! Steady-state scalar flux scores
        self % fluxScores(idx) = self % fluxScores(idx) + self % scalarFlux(idx)
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, self % nP * self % nCells
        ! Transient DNP scores 
        self % dnpScoresTD(idx) = self % dnpScoresTD(idx) + self % dnpTD(idx)
        
        ! Steady-state DNP scores 
        self % dnpScores(idx) = self % dnpScores(idx) + self % dnp(idx)       
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, self % nCells
        ! Transient fission source scores
!--> MK 240517
        source = real(self % fissionRate(idx), defReal) 
        self % fissionRateScore(idx,1) = self % fissionRateScore(idx,1) + source
        self % fissionRateScore(idx,2) = self % fissionRateScore(idx, 2) + source*source

        self % fissionScoreTD(idx) = self % fissionScoreTD(idx) + real(self % fissionSourceTD(idx), defReal) 
!<-- MK 240517
        
        ! Steady-state fission source scores
        self % fissionScore(idx) = self % fissionScore(idx) + self % fissionSource(idx)
    end do
    !$omp end parallel do
    
    ! keff scores
    if (t == 1) then
        self % keffScore(1) = self % keffScore(1) + self % keff
        self % keffScore(2) = self % keffScore(2) + self % keff * self % keff
    end if
        
  end subroutine accumulateFluxAndKeffScores
!<-- MK 230822  

  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxAndKeffScores(self,it, t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt), intent(in)                 :: t
    integer(shortInt)                             :: idx
    real(defReal)                                 :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it
    
    !$omp parallel do schedule(static)
    do idx = 1, self % nCells * self % nG
        ! Transient scalar flux scores 
        self % fluxScoresTD(idx,1) = self % fluxScoresTD(idx, 1) * N1
        self % fluxScoresTD(idx,2) = self % fluxScoresTD(idx, 2) * N1
        self % fluxScoresTD(idx,2) = Nm1 *(self % fluxScoresTD(idx,2) - &
                self % fluxScoresTD(idx,1) * self % fluxScoresTD(idx,1)) 
        if (self % fluxScoresTD(idx,2) <= ZERO) then
            self % fluxScoresTD(idx,2) = ZERO
        else
            self % fluxScoresTD(idx,2) = sqrt(self % fluxScoresTD(idx,2))
        end if
        
        ! Steady-state scalar flux score
        self % fluxScores(idx) = self % fluxScores(idx) * N1           
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, self % nP * self % nCells
        ! Transient DNP scores 
        self % dnpScoresTD(idx) = self % dnpScoresTD(idx) * N1
        
        ! Steady-state DNP scores 
        self % dnpScores(idx) = self % dnpScores(idx) * N1
    end do
    !$omp end parallel do

    !$omp parallel do schedule(static)
    do idx = 1, self % nCells
!--> MK 240517
         ! Transient fission source scores
        self % fissionScoreTD(idx) = self % fissionScoreTD(idx) * N1
        
        self % fissionRateScore(idx,1) = self % fissionRateScore(idx,1) * N1
        self % fissionRateScore(idx,2) = self % fissionRateScore(idx,2) * N1
        self % fissionRateScore(idx,2) = Nm1 *(self % fissionRateScore(idx,2) - &  
                self % fissionRateScore(idx,1) * self % fissionRateScore(idx,1)) 
        if (self % fissionRateScore(idx,2) <= ZERO) then
            self % fissionRateScore(idx,2) = ZERO
        else
            self % fissionRateScore(idx,2) = sqrt(self % fissionRateScore(idx,2))
        end if
!<-- MK 240517

        ! Steady-state fission source scores
        self % fissionScore(idx) = self % fissionScore(idx) * N1
    end do
    !$omp end parallel do

    ! keff scores
    if (t == 1) then
        self % keffScore(1) = self % keffScore(1) * N1
        self % keffScore(2) = self % keffScore(2) * N1
        self % keffScore(2) = sqrt(Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 
    end if

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    integer(shortInt)                             :: cIdx, g1
!--> MK 230307
    integer(shortInt), save                       :: idx, matIdx, i, g
!<-- MK 230307
    real(defReal), save                           :: vol, SigmaF
    type(particleState), save                     :: s
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal), dimension(:), allocatable      :: groupFlux, fiss, fissSTD
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(idx, matIdx, i, mat, matPtr, vol, s, SigmaF, g)

    call out % init(self % outputFormat)
    
    name = 'seed'
    call out % printValue(self % rand % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

!     name = 'Inactive_Cycles'
!     call out % printValue(self % inactive,name)
! 
!     name = 'Active_Cycles'
!     call out % printValue(self % active,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Total_Transport_Time'
    call out % printValue(self % time_transport,name)
    
    name = 'Clock_Time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print keff
    name = 'keff'
    call out % startBlock(name)
    call out % printResult(real(self % keffScore(1),defReal), real(self % keffScore(2),defReal), name)
    call out % endBlock()
    
!--> MK 230718   
    ! Print number of time steps
    name = 'Time_Steps'
    call out % printValue(self % nT,name)
    
    ! Print time step width
    name = 'Time_Step_Width'
    call out % printValue(real(self % tstp,defReal),name)
    
    ! Print time step width
    name = 'Cell_Number'
    call out % printValue(self % nCells,name)
    
    ! Print time step width
    name = 'Time_Steps_Printed'
    call out % printValue(self % nTOut,name)
    
    ! Print cycles when converged
    name = 'Total_Cycles'
    call out % printValue(self % totalIt,name)
    
    name = 'Total_Active_Cycles'
    call out % printValue(self % totalAct,name)
    
    name = 'Total_Inactive_Cycles'
    call out % printValue(self % totalInact,name)
!<-- MK 230718 

    ! Print cell volumes
    if (self % printVolume) then
      name = 'volume'
      call out % startBlock(name)
      resArrayShape = [size(self % volume)]
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        call out % addResult(real(self % volume(cIdx),defReal), ZERO)
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
!--> MK 240205
            idx = (cIdx - 1)* self % nG + g
!<-- MK 240205 
            fiss(i) = fiss(i) + vol * self % fluxScoresTD(idx,1) * SigmaF
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    vol * vol * self % fluxScoresTD(idx,2)*self % fluxScoresTD(idx,2) * SigmaF * SigmaF
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
        call out % addResult(real(fiss(idx),defReal), real(fissSTD(idx),defReal))
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
      resArrayShape = self % fluxMap % binArrayShape()
      allocate(fiss(self % fluxMap % bins(0)))
      allocate(fissSTD(self % fluxMap % bins(0)))
      fiss    = ZERO
      fissSTD = ZERO

      ! Find whether cells are in map and sum their contributions
      !$omp parallel do reduction(+: fiss, fissSTD)
      do cIdx = 1, self % nCells
        
        vol    =  self % volume(cIdx)
        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % fluxMap % map(s)

        if (i > 0) then
          do g = 1, self % nG
!--> MK 240205
            idx = (cIdx - 1)* self % nG + g
!<-- MK 240205 
            fiss(i) = fiss(i) + vol * self % fluxScoresTD(idx,1)
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    self % fluxScoresTD(idx,2)*self % fluxScoresTD(idx,2) * vol * vol
          end do
        end if

      end do
      !$omp end parallel do

      do i = 1,size(fissSTD)
        fissSTD(i) = sqrt(fissSTD(i))
        if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
      end do

      name = 'flux1G'
      call out % startBlock(name)
      call out % startArray(name, resArrayShape)
      ! Add all map elements to results
      do idx = 1, self % fluxMap % bins(0)
        call out % addResult(real(fiss(idx),defReal), real(fissSTD(idx),defReal))
      end do
      call out % endArray()
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
!--> MK 240205
            idx = (cIdx - 1)* self % nG + g
!<-- MK 240205 
            groupFlux(cIdx) = self % fluxScoresTD(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(real(groupFlux,defReal),name)
      end do
      do g1 = 1, self % nG
        name = 'std_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScoresTD(idx,2) /self % fluxScoresTD(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(real(groupFlux,defReal),name)
      end do
      call self % viz % finaliseVTK
    end if

  end subroutine printResults
  
!--> MK 231003
  !!
  !! Print transient flux at specified ouput time in output file (output = array)
  !!
  subroutine printTransientFlux(self,nOut,t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: nOut
    integer(shortInt), intent(in)                 :: t
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    character(nameLen)                            :: outputName
    integer(shortInt)                             :: idx, g, cIdx
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal)                                 :: tOut

    outputName = trim(self % outputFile)//'_t'//numToChar(nOut)
        
    call out % init(self % outputFormat)
    
    name = 't'
    
    tOut = (t - 1) * self % tstp
    call out % printValue(tOut,name)
    
    ! Print transient fluxes
    resArrayShape = [size(self % volume)]
    do g = 1, self % nG
        name = 'flux_g'//numToChar(g)//'t'
        call out % startBlock(name)
        call out % startArray(name, resArrayShape)
        do cIdx = 1, self % nCells 
            idx = self % nG * (cIdx - 1) + g 
            call out % addResult(real(self % fluxScoresTD(idx,1),defReal), real(self % fluxScoresTD(idx,2),defReal))           
        end do
        call out % endArray()
        call out % endBlock()
    end do
    
    call out % writeToFile(outputName)

  end subroutine printTransientFlux
!<-- MK 231003 

!--> MK 231005
  !!
  !! Write file for fission rate disctribution at specified ouput time (output = map)
  !!
  subroutine printFissionRate(self,nOut,t)
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt), intent(in)                 :: nOut
    integer(shortInt), intent(in)                 :: t
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    character(nameLen)                            :: outputName
    integer(shortInt), save                       :: idx, matIdx, i, g, matIdx0
    integer(shortInt)                             :: cIdx
    real(defReal)                                 :: tOut
    real(defReal), save                           :: vol, SigmaF
    type(particleState), save                     :: s
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal), dimension(:), allocatable      :: fiss, fissSTD
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(idx, matIdx, i, mat, matPtr, vol, s, SigmaF, g, matIdx0)

    outputName = trim(self % outputFile)//'_fission_t'//numToChar(nOut)
    
    call out % init(self % outputFormat)
    
    name = 't'
    
    tOut = (t - 1) * self % tstp
    call out % printValue(tOut,name)
    
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
        mat    => baseMgNeutronMaterial_CptrCast(matPtr)    !TODO is mat necessary?
        vol    =  self % volume(cIdx)

        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % resultsMap % map(s)

        if ((i > 0) .AND. (matIdx /= 8)) then                               !TODO matIdx /= 7 for C5G7-TD3, matIdx /= 8 for C5G7-TD1
!--> MK 240517
!           do g = 1, self % nG
!             matIdx0 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG + g    !TODO This works but should be checked again in the future
!             if (self % nuSigmaF(matIdx0) /= 0.0_defFlt) then
!                 SigmaF = self % nuSigmaF(matIdx0) / self % nu(self % nG * (matIdx - 1) + g)
!             else
!                 SigmaF = ZERO
!             end if
!    
!             idx = self % nG * (cIdx - 1) + g 
!             fiss(i) = fiss(i) + vol * self % fluxScoresTD(idx,1) * SigmaF
!             ! Is this correct? Also neglects uncertainty in volume - assumed small. Covariance neglected
!             fissSTD(i) = fissSTD(i) + &
!                    vol * vol * self % fluxScoresTD(idx,2)*self % fluxScoresTD(idx,2) * SigmaF * SigmaF
              fiss(i) = fiss(i) + self % fissionRateScore(cIdx,1)
              fissSTD(i) = fissSTD(i) + self % fissionRateScore(cIdx,2) * self % fissionRateScore(cIdx,2)
!               fiss(i) = vol * self % fissionRateScore(idx,1)
!               fissSTD(i) = vol * vol * self % fissionRateScore(idx,2)
!           end do
!<-- MK 240517
        end if

    end do
    !$omp end parallel do

!--> MK 240517
    do i = 1,size(fissSTD)
        fissSTD(i) = sqrt(fissSTD(i))
        if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
    end do
!<-- MK 240517

    name = 'fissionRate'
    call out % startBlock(name)
    call out % startArray(name, resArrayShape)
    ! Add all map elements to results
    do idx = 1, self % resultsMap % bins(0)
        call out % addResult(real(fiss(idx),defReal), real(fissSTD(idx),defReal))
    end do
    call out % endArray()
    ! Output tally map
    call self % resultsMap % print(out)
    call out % endBlock()
      
    deallocate(fiss)
    deallocate(fissSTD)    
      
    call out % writeToFile(outputName)

  end subroutine printFissionRate
!<-- MK 231005   
    
  !!
  !! Print settings of the random ray calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(dynamicRRPhysicsPackage_TI), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RANDOM RAY TIME-DEPENDENT CALCULATION /\/\"
    print *, "Using "//numToChar(self % inactive)// " iterations for "&
              //"the inactive cycles"
    print *, "Using "//numToChar(self % active)// " iterations for "&
              //"the active cycles"
    print * 
    print *, "Rays per cycle: "// numToChar(self % pop)
    print *, "Ray dead length: "//numToChar(self % dead)
    print *, "Ray termination length: "//numToChar(self % termination)
!--> MK 230302   
    print *, "Number of timesteps: "//numToChar(self % nT)
    print *, "Total simulated time: "//numToChar(self % nT * real(self % tstp, defReal))
!<-- MK 230302  
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
    class(dynamicRRPhysicsPackage_TI), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean Nuclear Data, Geometry and visualisation
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
    
    if(allocated(self % locks)) then
      do i = 1, self % nCells
        call OMP_destroy_lock(self % locks(i))
      end do
      deallocate(self % locks)
    end if
    self % nCells    = 0
    self % nMat      = 0
    
    ! Clean nuclear data
    if(allocated(self % sigmaT)) deallocate(self % sigmaT)
    if(allocated(self % sigmaS)) deallocate(self % sigmaS)
    if(allocated(self % nusigmaF)) deallocate(self % nuSigmaF)
!--> MK 240517
    if(allocated(self % sigmaF)) deallocate(self % sigmaF)
!<-- MK 240517
    if(allocated(self % chi)) deallocate(self % chi)
    if(allocated(self % nu)) deallocate(self % nu)
!--> MK 240403
    if(allocated(self % vel)) deallocate(self % vel)
    if(allocated(self % chiP)) deallocate(self % chiP)
    if(allocated(self % chiD)) deallocate(self % chiD)
    if(allocated(self % lambda)) deallocate(self % lambda)
    if(allocated(self % beta)) deallocate(self % beta)
    
    if(allocated(self % omega0)) deallocate(self % omega0)
    if(allocated(self % omegaN)) deallocate(self % omegaN)
    if(allocated(self % omega_1)) deallocate(self % omega_1)
    if(allocated(self % omega_2)) deallocate(self % omega_2)
!<-- MK 240403

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
!--> MK 230718   
    self % nT          = 0
    self % outInterval = 0
    self % window      = 0
    self % tstp        = ZERO
    self % initVal     = ZERO
    self % totalIt     = 0
    self % totalAct    = 0
    self % totalInact  = 0
    self % nTOut       = 0
    self % reduceLines = .false.
    self % isotropic = .false.
!<-- MK 230718
!--> MK 240523   
    self % AccVol = ZERO
    self % sarePrec = ZERO
    self % rmsPrec = ZERO
!<-- MK 240523  

    ! Clean scores
    self % keff        = ZERO
    self % keffScore   = ZERO
    if(allocated(self % scalarFluxTD)) deallocate(self % scalarFluxTD)
    if(allocated(self % prevFluxTD)) deallocate(self % prevFluxTD)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScoresTD)) deallocate(self % fluxScoresTD)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % sourceTD)) deallocate(self % sourceTD)
    
!--> MK 240403
    if(allocated(self % dnpScoresTD)) deallocate(self % dnpScoresTD)
    if(allocated(self % dnpTD)) deallocate(self % dnpTD)
    if(allocated(self % prevDnpTD)) deallocate(self % prevDnpTD)
    if(allocated(self % dnpScores)) deallocate(self % dnpScores)
    if(allocated(self % dnp)) deallocate(self % dnp)
    if(allocated(self % prevDnp)) deallocate(self % prevDnp)
    if(allocated(self % fluxHistTD)) deallocate(self % fluxHistTD)
    if(allocated(self % fluxHist)) deallocate(self % fluxHist)

    if(allocated(self % prevFissionTD_1)) deallocate(self % prevFissionTD_1)
    if(allocated(self % prevFissionTD_2)) deallocate(self % prevFissionTD_2)
    if(allocated(self % fissionScoreTD)) deallocate(self % fissionScoreTD)
    if(allocated(self % fissionSourceTD)) deallocate(self % fissionSourceTD)
    if(allocated(self % prevFission_1)) deallocate(self % prevFission_1)
    if(allocated(self % prevFission_2)) deallocate(self % prevFission_2)
    if(allocated(self % fissionScore)) deallocate(self % fissionScore)
    if(allocated(self % fissionSource)) deallocate(self % fissionSource)
    if(allocated(self % fissionRateScore)) deallocate(self % fissionRateScore)
    if(allocated(self % fissionRate)) deallocate(self % fissionRate)
!<-- MK 240403
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

end module dynamicRRPhysicsPackage_class_TI
