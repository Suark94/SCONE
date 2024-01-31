module randomRayPhysicsPackage_class

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_TD_func,          only : exponential_TD
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
  real(defFlt), parameter :: volume_tolerance = 1.0E-12

  !!
  !! Physics package to perform The Random Ray Method (TRRM) eigenvalue calculations
  !!
  !! Tracks rays across the geometry, attenuating their flux. After some dead length,
  !! rays begin scoring to estimates of the scalar flux and volume. Each ray has a
  !! uniform termination length, after which it is stopped and the next ray is tracked.
  !! Once all rays have been tracked, a cycle concludes and fluxes, sources, and keff
  !! are updated.
  !!
  !! Both inactive and active cycles occur, as in Monte Carlo. These can be terminated
  !! after a specified number of iterations or on reaching some chosen convergence
  !! criterion (though the latter hasn't been implemented yet).
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
  !!     type randomRayPhysicsPackage;
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
  !!   nT          -> Number of time steps
  !!   tstp        -> Time step size
  !!   initVal     -> Initialisation value for scalar flux
  !!   chiP
  !!   chiD
  !!   lambda
  !!   beta
  !!   window      -> Width of moving window for Shannon entropy
  !!   reduceLines -> reduce output lines that are printed while running
  !!   itConv      -> total number of iterations in transient calculation
!<-- MK 230601
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: randomRayPhysicsPackage
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
    integer(shortInt)  :: nT         = 0
    integer(shortInt)  :: nTOut      = 0
    integer(shortInt)  :: window     = 0
    real(defFlt)       :: tstp       = ZERO
    real(defFlt)       :: initVal    = ZERO
    integer(shortInt)  :: itConv     = 0         
    integer(shortInt)  :: nPrintOut  = 0
    logical(defBool)   :: reduceLines = .false.
!<-- MK 230718   

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable     :: sigmaT
    real(defFlt), dimension(:), allocatable     :: nuSigmaF
    real(defFlt), dimension(:), allocatable     :: sigmaS
    real(defFlt), dimension(:), allocatable     :: chi
!--> DNP 230918
    real(defFlt), dimension(:), allocatable     :: chiP
    real(defFlt), dimension(:), allocatable     :: chiD
    real(defFlt), dimension(:), allocatable     :: lambda
    real(defFlt), dimension(:), allocatable     :: beta
    real(defFlt), dimension(:), allocatable     :: DNPScores
    real(defFlt), dimension(:), allocatable     :: DNP
    real(defFlt), dimension(:), allocatable     :: prevDNP
!<-- DNP 230918
!--> MK 230823
    real(defFlt), dimension(:), allocatable     :: vel
    real(defFlt), dimension(:), allocatable     :: fluxHist
!<-- MK 230823

    ! Results space
    real(defFlt)                                :: keff
    real(defFlt), dimension(2)                  :: keffScore
    real(defFlt), dimension(:), allocatable     :: scalarFlux
    real(defFlt), dimension(:), allocatable     :: prevFlux
    real(defFlt), dimension(:,:), allocatable   :: fluxScores
    real(defFlt), dimension(:), allocatable     :: source
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
    procedure, private :: transportSweep_steady
    procedure, private :: transportSweep_transient
    procedure, private :: calculateEntropy
    procedure, private :: sourceUpdateKernel
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings
    procedure, private :: printTransientFlux
    procedure, private :: printFissionRate

  end type randomRayPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(randomRayPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
!--> MK 230718
    integer(shortInt)                             :: seed_temp, i, g, g1, m, t, nP, nPNew, p, idx0
    real(defReal)                                 :: tstpReal, initReal, outTemp
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
    class(baseMgNeutronMaterial), pointer         :: mat
    class(materialHandle), pointer                :: matPtr
!--> DNP 230918
    class(fissionMG), pointer                     :: fiss
    real(defFlt)                                  :: omega, timeXS
!<-- DNP 230918
    character(100), parameter :: Here = 'init (randomRayPhysicsPackage_class.f90)'

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
 
    call dict % get(self % nPrintOut, 'tOut')
    outTemp = real(self % nT / self % nPrintOut, defFlt)
    self % nTOut = NINT(outTemp)

    call dict % get(tstpReal, 'tstp')
    self % tstp = real(tstpReal, defFlt)
    call dict % get(initReal, 'initVal')
    self % initVal = real(initReal, defFlt)
    
    if (self % nT > 1) then
        self % nT = self % nT + 1 
    end if
    
    call dict % getOrDefault(self % reduceLines, 'reduceLines', .false.)
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
!--> MK 230823
    allocate(self % scalarFlux(self % nCells * self % nG))
    allocate(self % prevFlux(self % nCells * self % nG))
    allocate(self % fluxScores(self % nCells * self % nG, 2))
    allocate(self % fluxHist(self % nCells * self % nG))
    allocate(self % source(self % nCells * self % nG))
!<-- MK 230823
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
    allocate(self % chi(self % nMat * self % nG))
    allocate(self % sigmaS(self % nMat * self % nG * self % nG * self % nT))
    allocate(self % vel(self % nMat * self % nG))
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
        if (self % nT > 1) then 
            self % vel(self % nG * (m - 1) + g) = real(mat % getVel(g, self % rand),defFlt)
        end if
        do t = 1, self % nT
            self % sigmaT(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = real(mat % getTotalXS(g, self % rand),defFlt) 
            self % nuSigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
            real(mat % getNuFissionXS(g, self % rand),defFlt) 
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
    
!--> MK 230718
    do m = 1, self % nMat       !TODO This is only provisional
      matPtr  => self % mgData % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      do g = 1, self % nG
        do t = 2, self % nT
            self % nuSigmaF(self % nG * self % nT * (m - 1) + (t - 1) * self % nG + g) = &
            real(mat % getNuFissionXS(g, self % rand)*1.001,defFlt) 
        end do
      end do
    end do
!<-- MK 230718

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
    
    allocate(self % DNPScores(self % nCells * self % nP))
    allocate(self % prevDNP(self % nCells * self % nP))
    allocate(self % DNP(self % nCells * self % nP))
!<-- DNP 230918

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(randomRayPhysicsPackage), intent(inout) :: self

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
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defFlt)                                  :: hitRate, ONE_KEFF, scaleFac
    real(defFlt), dimension(self % window)        :: conv
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: stoppingCriterion, isActive
    integer(shortInt)                             :: i, itInac, itAct, it, t
    integer(longInt), save                        :: ints
    integer(longInt)                              :: intersections
    logical(defBool)                              :: convBool, entropyBool
    integer(shortInt)                             :: cIdx, pIdx, tCnt, outCnt, actCycles
    integer(shortInt), save                       :: idx0, g
    !$omp threadprivate(pRNG, r, ints, idx0,g)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)
    
    self % keff       = 1.0_defFlt
    scaleFac          = 0.0_defFlt
    self % keffScore  = 0.0_defFlt
    self % DNPScores = 0.0_defFlt
    self % prevDNP = 0.0_defFlt
    self % fluxHist = 0.0_defFlt
    tCnt = 0
    outCnt = 0
    self % itConv = 0
    self % volumeTracks = 0.0_defReal       ! To score volumes over time steps
    self % volume       = 0.0_defReal       !
    
    ! Time-stepping
    do t = 1, self % nT
    
        tCnt = tCnt + 1

        ! Initialise fluxes 
        self % scalarFlux = 0.0_defFlt
        
        if (t==1) then
            self % prevFlux = 1.0_defFlt  
        else
            self % prevFlux = self % fluxScores(:,1)
        end if
        
        self % fluxScores = 0.0_defFlt
        self % source     = 0.0_defFlt

        ! Initialise other results
        self % cellHit      = 0
!         self % volume       = 0.0_defReal   ! To score volumes over cycles only
!         self % volumeTracks = 0.0_defReal   !   
        
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
        conv = 0.0_defFlt
        
        ! Power iteration
        do while( stoppingCriterion )
            self % DNP = 0.0_defFlt
            
            if (isActive) then
                itAct = itAct + 1
            else
                itInac = itInac + 1
            end if
            it = itInac + itAct
                    
            self % itConv = self % itConv + 1      ! For metrics, also volume accumlation MK 231011
                    
            if (t == 1) then
                ONE_KEFF = 1.0_defFlt / self % keff
            else
                ONE_KEFF = 1.0_defFlt / self % keffScore(1)
            end if

            !$omp parallel do schedule(static)
            do i = 1, self % nCells
                call self % sourceUpdateKernel(i, ONE_KEFF, t)
            end do
            !$omp end parallel do
            
            ! Reset and start transport timer
            call timerReset(self % timerTransport)
            call timerStart(self % timerTransport)
            intersections = 0
            
            if (t == 1) then
                !$omp parallel do schedule(runtime) reduction(+: intersections) 
                do i = 1, self % pop
                
                    ! Set seed
                    pRNG = self % rand
                    call pRNG % stride(i)
                    r % pRNG => pRNG 

                    ! Set ray attributes
                    call self % initialiseRay(r)

                    ! Transport ray until termination criterion met
                    call self % transportSweep_steady(r,ints)
                    intersections = intersections + ints

                end do
                !$omp end parallel do
            else
                !$omp parallel do schedule(runtime) reduction(+: intersections) 
                do i = 1, self % pop
                
                    ! Set seed
                    pRNG = self % rand
                    call pRNG % stride(i)
                    r % pRNG => pRNG 

                    ! Set ray attributes
                    call self % initialiseRay(r)

                    ! Transport ray until termination criterion met
                    call self % transportSweep_transient(r,ints,t)
                    intersections = intersections + ints

                end do
                !$omp end parallel do      
            end if
            
            call timerStop(self % timerTransport)

            ! Update RNG on master thread
            call self % rand % stride(self % pop + 1)
            
            ! Normalise flux estimate and combines with source
            call self % normaliseFluxAndVolume(it, t)
            
            ! Calculate proportion of cells that were hit
            hitRate = real(sum(self % cellHit),defFlt) / self % nCells
            self % cellHit = 0

            if (t == 1) then
                call self % calculateKeff()  
                
                scaleFac = (self % initVal / sum(self % scalarFlux))
                self % scalarFlux = self % scalarFlux * scaleFac
            end if

            ! Accumulate flux scores
            if (isActive) call self % accumulateFluxAndKeffScores(t)
                
            call self % calculateEntropy(conv, entropyBool)

            if (t == 1) then
                stoppingCriterion = (itAct < self % active * 4)
                actCycles = self % active * 4
            else
                stoppingCriterion = (itAct < self % active)
                actCycles = self % active
            end if

            if (convBool .AND. (it >= self % inactive) .AND. entropyBool) then   
                    isActive = .true.
                    convBool = .false.
            end if

            ! Set previous iteration flux to scalar flux
            ! and zero scalar flux
            call self % resetFluxes()

            ! Calculate times
            call timerStop(self % timerMain)
            elapsed_T = timerTime(self % timerMain)
            transport_T = timerTime(self % timerTransport)
            self % time_transport = self % time_transport + transport_T

            ! Predict time to end
            end_T = real(self % active + self % inactive, defReal) * elapsed_T / it
            T_toEnd = max(ZERO, end_T - elapsed_T)

            ! Display progress
            call printFishLineR(it)
            print *
            if(isActive) then
                print *, 'Iteration: ', numToChar(itAct), ' of ', numToChar(actCycles)
                if (.not. self % reduceLines) print *,'Active iterations'
            else
                print *, 'Iteration: ', numToChar(it), ' of ', numToChar(self % inactive)
                if (.not. self % reduceLines) print *,'Inactive iterations'
            end if
            print *, 'Time step: ', numToChar(t), ' of ', numToChar(self % nT)
            if (.not. self % reduceLines) print *, 'Cell hit rate: ', trim(numToChar(real(hitRate,defReal)))
            if (t ==1) then
            if (.not. self % reduceLines) print *, 'keff: ', trim(numToChar(real(self % keff,defReal)))
            else
            if (.not. self % reduceLines) print *, 'keff: ', trim(numToChar(real(self % keffScore(1),defReal)))
            end if
            print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
            print *, 'Time per integration (ns): ', &
                    trim(numToChar(transport_T*10**9/(self % nG * intersections)))
            !print *, 'Volume: ', self % volume(10)
        end do
        
        ! Finalise flux scores
        call self % finaliseFluxAndKeffScores(itAct, t)  
        
        ! Update previous scalar flux
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
            idx0 = (cIdx - 1) * self % nG
            do g = 1, self % nG
                self % fluxHist(idx0 + g) = self % fluxScores(idx0 + g, 1)
            end do
        end do
        !$omp end parallel do
        
        ! Update previous precursors
        !$omp parallel do schedule(static)
        do pidx = 1, self % nP * self % nCells
        self % prevDNP(pidx) = self % DNPScores(pidx)
        self % DNPScores(pidx) = ZERO
        end do
        !$omp end parallel do

        ! Print transient output
        if (tCnt == self % nTOut .OR. t == 1 .OR. t == self % nT) then
            outCnt = outCnt + 1
            if (self % printFlux) call self % printTransientFlux(outCnt,t) 
            if (self % mapFission) call self % printFissionRate(outCnt,t) 
            tCnt = 0
        end if
    
    end do
    
    self % nPrintOut = outCnt

  end subroutine cycles
!<-- MK 230823 
  
!--> MK 240129
  subroutine calculateEntropy(self, conv, entropyBool)
      class(randomRayPhysicsPackage), intent(inout) :: self
      real(defFlt), dimension(self % window), intent(inout)   :: conv
      logical(defBool), intent(out)                 :: entropyBool
      real(defFlt)                                  :: conv_mean1, conv_mean2, conv_std2
      real(defFlt)                                  :: sumFlux
      real(defFlt)                                  :: convTemp
      integer(shortInt)                             :: cIdx, halfwindow
      integer(shortInt), save                       :: idx
      real(defFlt), save                            :: P
      !$omp threadprivate(P, idx)
      
      halfwindow = INT(self % window * 0.5)
      
    ! Convergence - Shannon entropy
      sumFlux = sum(self % scalarFlux)
      conv(1:(self % window-1)) = conv(2:self % window)
      convTemp = 0.0_defFlt
      
      !$omp parallel do schedule(static) reduction(+: convTemp) 
      do cIdx = 1, self % nCells
        idx = self % nG * (cIdx - 1)
        P = sum(self % scalarFlux(idx + 1 : idx + self % nG)) / sumFlux
        if (P <= ZERO) then
            convTemp = convTemp
        else
            convTemp = convTemp - P * log(P)
        end if
      end do
      !$omp end parallel do
      
       conv(self % window) = convTemp
       conv_mean1 = sum(conv(1 : halfwindow)) / halfwindow
       conv_mean2 = sum(conv((halfwindow+1) : self % window)) / halfwindow
       conv_std2 = (1.0_defFlt / (halfwindow-1)) * sum((conv((halfwindow+1) : self % window) - conv_mean2) &
            * (conv((halfwindow+1) : self % window) - conv_mean2))
      if (conv_std2 <= ZERO) then
        conv_std2 = ZERO
      else
        conv_std2 = sqrt(conv_std2)
      end if
      
      entropyBool = (ABS(conv_mean1 - conv_mean2) < conv_std2)
      
  end subroutine calculateEntropy
!<-- MK 240129

  !!
  !! Initialises rays: samples initial position and direction,
  !! and performs the build operation
  !!
  subroutine initialiseRay(self, r)
    class(randomRayPhysicsPackage), intent(inout) :: self
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

!--> MK 240129
  !!
  !! Steady-state transport sweep
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep_steady(self, r, ints)
    class(randomRayPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, event, matIdx0, baseIdx, baseIdx0
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, fluxVec
    real(defFlt), pointer, dimension(:)                   :: totVec, velVec, scalarVec, sourceVec
    real(defFlt)                                          :: lenFlt
    real(defReal), dimension(3)                           :: r0, mu0

    ! Set initial angular flux to angle average of cell source (space)
    cIdx = r % coords % uniqueID

    baseIdx = (cIdx - 1) * self % nG
    sourceVec => self % source(baseIdx + 1 : baseIdx + self % nG)

    !$omp simd
    do g = 1, self % nG
        fluxVec(g) = sourceVec(g)
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
        
        ! Cache total cross section
        baseIdx0 = (matIdx - 1) * self % nG 
        velVec => self % vel((baseIdx0 + 1) : (baseIdx0 + self % nG))
        
        baseIdx0 = (matIdx - 1) * self % nG * self % nT
        totVec => self % sigmaT((baseIdx0 + 1) : (baseIdx0 + self % nG))
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
      
      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % source(baseIdx + 1 : baseIdx + self % nG)
      
      ints = ints + 1

        ! Calculate steady-state flux
        !$omp simd
        do g = 1, self % nG
            attenuate(g) = exponential_TD(totVec(g) * lenFlt)               ! exponential_TD(tau) = (1 - exp(-tau)/tau)
            delta(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g) * totVec(g) * lenFlt
            fluxVec(g) = fluxVec(g) - delta(g)
        end do

      ! Accumulate to scalar flux
      if (activeRay) then
        scalarVec => self % scalarFlux(baseIdx + 1 : baseIdx + self % nG)
        
        call OMP_set_lock(self % locks(cIdx))
        !$omp simd
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + delta(g)
        enddo

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
      end if
    end do

  end subroutine transportSweep_steady
!<-- MK 240129

!--> MK 240129
  !!
  !! Transient transport sweep
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep_transient(self, r, ints, t)
    class(randomRayPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                              :: r
    integer(shortInt), intent(in)                         :: t
    integer(longInt), intent(out)                         :: ints
    integer(shortInt)                                     :: matIdx, g, cIdx, event, matIdx0, baseIdx, baseIdx0
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt), dimension(self % nG)                    :: attenuate, delta, flux_avg, flux_avgprv, denom, deltaAcc, fluxVec
    real(defFlt), pointer, dimension(:)                   :: totVec, velVec, scalarVec, sourceVec
    real(defFlt)                                          :: lenFlt
    real(defReal), dimension(3)                           :: r0, mu0

    ! Set initial angular flux to angle average of cell source (space)
    cIdx = r % coords % uniqueID

    baseIdx = (cIdx - 1) * self % nG
    sourceVec => self % source(baseIdx + 1 : baseIdx + self % nG)

    !$omp simd
    do g = 1, self % nG
        fluxVec(g) = sourceVec(g)
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
        
        ! Cache total cross section
        baseIdx0 = (matIdx - 1) * self % nG 
        velVec => self % vel((baseIdx0 + 1) : (baseIdx0 + self % nG))
        
        baseIdx0 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG
        totVec => self % sigmaT((baseIdx0 + 1) : (baseIdx0 + self % nG))
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
      
      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % source(baseIdx + 1 : baseIdx + self % nG)
      
      ints = ints + 1

      ! Isotropic approximation
      scalarVec => self % fluxHist(baseIdx + 1 : baseIdx + self % nG)
        
      !$omp simd
      do g = 1, self % nG
        flux_avgprv(g) = real(scalarVec(g) * ONE_FOUR_PI, defFlt)
      end do
        
      ! Transient calculations
      !$omp simd                                                    
      do g = 1, self % nG
        attenuate(g) = exponential_TD(totVec(g) * lenFlt)
        denom(g) = 1 / (totVec(g) * velVec(g) * self % tstp)
        flux_avg(g) = (attenuate(g) * (fluxVec(g) - sourceVec(g) - flux_avgprv(g) * denom(g)) &
            + sourceVec(g) + flux_avgprv(g) * denom(g)) / (1 + denom(g) * (1 - attenuate(g)))                         
        delta(g) = (fluxVec(g) - sourceVec(g) + (flux_avg(g) - flux_avgprv(g)) * denom(g)) * attenuate(g) * totVec(g) * lenFlt     ! * totVec(g) * lenFlt to account for changed exponetial function
        fluxVec(g) = fluxVec(g) - delta(g)
        deltaAcc(g) = delta(g) - lenFlt * (flux_avg(g) - flux_avgprv(g))/(velVec(g) * self % tstp)
        flux_avgprv(g) = flux_avg(g) 
      end do

      ! Accumulate to scalar flux
      if (activeRay) then
        scalarVec => self % scalarFlux(baseIdx + 1 : baseIdx + self % nG)
        
        call OMP_set_lock(self % locks(cIdx))
        !$omp simd
        do g = 1, self % nG
            scalarVec(g) = scalarVec(g) + deltaAcc(g)
        enddo

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
      end if
    end do

  end subroutine transportSweep_transient
!<-- MK 240129

!--> MK 230718
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it, t)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt), intent(in)                 :: t
    real(defFlt)                                  :: norm
    real(defReal)                                 :: normVol
    integer(shortInt), save                       :: g, idx, matIdx
    integer(shortInt)                             :: cIdx
    real(defFlt), save                            :: total, vol
    !$omp threadprivate(total, vol, idx, g, matIdx)

    norm = real(ONE / self % lengthPerIt, defFlt)   
    normVol = ONE / (self % lengthPerIt * self % itConv)        ! To score volumes over time steps
!    normVol = ONE / (self % lengthPerIt * it)                  ! To score volumes over cycles only

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
          self % scalarFlux(idx) = real(self % scalarFlux(idx) * FOUR_PI, defFlt) * norm  / (total * vol)
        end if
        
        self % scalarFlux(idx) = self % scalarFlux(idx) + real(self % source(idx) * FOUR_PI, defFlt)
        
        if (self % scalarFlux(idx) < 0) self % scalarFlux(idx) = 0.0_defFlt
      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
!<-- MK 230718

!--> MK 230822   
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx, ONE_KEFF, t)
    class(randomRayPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                         :: cIdx
    real(defFlt), intent(in)                              :: ONE_KEFF
    integer(shortInt), intent(in)                         :: t
    real(defFlt)                                          :: scatter, fission, DNP_sum, ONE_MIN_B
    real(defFlt), dimension(:), pointer                   :: nuFission, total, scatterXS
    real(defFlt), dimension(:), pointer                   :: beta, lambda, chiP, chiD, chiDVec, DNP_prev, DNPVec
    integer(shortInt)                                     :: matIdx, g, gIn, baseIdx, idx, matIdx1, matIdx2, DNPIdx, p
    integer(shortInt)                                     :: baseIdx0, baseIdx1
    real(defFlt), pointer, dimension(:)                   :: fluxVec, scatterVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    matIdx1 = (matIdx - 1) * self % nG 
    baseIdx0 = matIdx1 * self % nP
    baseIdx1 = (matIdx - 1) * self % nP
    
    !chi => self % chi(matIdx1 + (1):(self % nG))                 
    chiP => self % chiP((matIdx1 + 1) : (matIdx1 + self % nG))
    chiD => self % chiD((baseIdx0 + 1) : (baseIdx0 + self % nG * self % nP))
    
    beta => self % beta((baseIdx1 + 1) : (baseIdx1 + self % nP))
    ONE_MIN_B = 1.0_defFlt - sum(beta)
    lambda => self % lambda((baseIdx1 + 1) : (baseIdx1 + self % nP))

    ! Obtain XSs
    matIdx2 = (matIdx - 1) * self % nG * self % nT + (t - 1) * self % nG
        
    total => self % sigmaT((matIdx2 + 1) : (matIdx2 + self % nG))                        
    scatterXS => self % sigmaS((matIdx2 * self % nG + 1) : (matIdx2 * self % nG + self % nG*self % nG))
    nuFission => self % nuSigmaF((matIdx2 + 1) : (matIdx2 + self % nG)) 
    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx+1) : (baseIdx+self % nG))
    
    DNPIdx = self % nP * (cIdx - 1)
    DNP_prev => self % prevDNP((DNPIdx+1) : (DNPIdx+self % nP))
    DNPVec => self % DNP((DNPIdx+1) : (DNPIdx+self % nP))
    
    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)                                 
    do gIn = 1, self % nG                                         
        fission = fission + fluxVec(gIn) * nuFission(gIn) 
    end do
    
!--> DNP 230918
    ! Calculate DNPs
    if (t==1) then
        !$omp simd
        do p = 1, self % nP 
            DNPVec(p) = (beta(p) / lambda(p)) * ONE_KEFF * fission
        end do
    else
        !$omp simd
        do p = 1, self % nP 
            DNPVec(p) = (beta(p) * ONE_KEFF * fission * self % tstp + DNP_prev(p))/(1 + lambda(p) * self % tstp)
        end do
    end if
!<-- DNP 230918

    do g = 1, self % nG
!--> DNP 230918
        ! Calculate DNP source
        DNP_sum = 0.0_defFlt
        baseIdx0 = self % nP * (g - 1)
        chiDVec => chiD((baseIdx0 + 1) : (baseIdx0 + self % nP))
            
        !$omp simd reduction(+:DNP_sum) 
        do p = 1, self % nP 
            DNP_sum = DNP_sum + chiDVec(p) * lambda(p) * DNPVec(p)
        end do
!<-- DNP 230918
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

        self % source(idx) = ONE_MIN_B * chiP(g) * fission * ONE_KEFF + scatter + DNP_sum
        self % source(idx) = self % source(idx) / (real(total(g) * FOUR_PI, defFlt))
    end do

  end subroutine sourceUpdateKernel
!<-- MK 230822   

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
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
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxes(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt)                             :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = ZERO
    end do
    !$omp end parallel do
    
  end subroutine resetFluxes

!--> MK 230822   
  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self, t)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: t
    real(defFlt), save                            :: flux, sourceflux, DNP_val
    integer(shortInt)                             :: idx
    !$omp threadprivate(flux, sourceflux, DNP_val)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarFlux(idx)
      sourceflux = self % source(idx)
      
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux*flux
    end do
    !$omp end parallel do
    
    !$omp parallel do schedule(static)
    do idx = 1, self % nP * self % nCells
      DNP_val = self % DNP(idx)
      self % DNPScores(idx) = self % DNPScores(idx) + DNP_val
    end do
    !$omp end parallel do
    
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
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt), intent(in)                 :: t
    integer(shortInt)                             :: idx
    real(defFlt)                                  :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defFLt/(it - 1)
    else
      Nm1 = 1.0_defFlt
    end if
    N1 = 1.0_defFlt/it
    
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
    end do
    !$omp end parallel do
    
!--> DNP 230919
    !$omp parallel do schedule(static)
    do idx = 1, self % nP * self % nCells
      self % DNPScores(idx) = self % DNPScores(idx) * N1
    end do
    !$omp end parallel do
!<-- DNP 230919

!--> MK 230822 
    if (t == 1) then
        self % keffScore(1) = self % keffScore(1) * N1
        self % keffScore(2) = self % keffScore(2) * N1
        self % keffScore(2) = sqrt(Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 
    end if
!<-- MK 230822 

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    integer(shortInt)                             :: cIdx, g1
!--> MK 230307
    integer(shortInt), save                       :: idx, matIdx, i, g! , t
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

    name = 'Inactive_Cycles'
    call out % printValue(self % inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % active,name)

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
    call out % printValue(self % nPrintOut,name)
    
    ! Print cycles when converged
    name = 'Total_Cycles'
    call out % printValue(self % itConv,name)
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
      
!     ! Print steady-state fluxes
!     if (self % printFlux) then
!       resArrayShape = [size(self % volume)]
!       do g = 1, self % nG
!         name = 'flux_g'//numToChar(g)
!         call out % startBlock(name)
!         call out % startArray(name, resArrayShape)
!         do cIdx = 1, self % nCells
! !--> MK 230302   
!           idx = (cIdx - 1)* self % nG + g
! !<-- MK 230302  
!           call out % addResult(real(self % fluxScores(idx,1),defReal), real(self % fluxScores(idx,2),defReal))
!         end do
!         call out % endArray()
!         call out % endBlock()
!       end do
!     end if   

! !--> MK 230307
!     ! Print transient fluxes
!     if (self % printFlux) then
!       resArrayShape = [size(self % volume) * self % nT]
!       do g = 1, self % nG
!         name = 'flux_g'//numToChar(g)//'t'
!         call out % startBlock(name)
!         call out % startArray(name, resArrayShape)
!         do t = 1, self % nT
!             do cIdx = 1, self % nCells 
!                 idx = (t - 1) * self % nG * self % nCells + self % nG * (cIdx - 1) + g 
!                 call out % addResult(real(self % fluxScores(idx,1),defReal), real(self % fluxScores(idx,2),defReal))           
!             end do
!         end do
!         call out % endArray()
!         call out % endBlock()
!       end do
!     end if 
! !<-- MK 230307 

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
            idx = (cIdx - 1)* self % nG + g
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1)
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    self % fluxScores(idx,2)*self % fluxScores(idx,2) * vol * vol
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
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(real(groupFlux,defReal),name)
      end do
      do g1 = 1, self % nG
        name = 'std_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,2) /self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(real(groupFlux,defReal),name)
      end do
      call self % viz % finaliseVTK
    end if

  end subroutine printResults
  
!--> MK 231003
  !!
  !! Print transient flux
  !!
  subroutine printTransientFlux(self,nOut,t)
    class(randomRayPhysicsPackage), intent(inout) :: self
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
            call out % addResult(real(self % fluxScores(idx,1),defReal), real(self % fluxScores(idx,2),defReal))           
        end do
        call out % endArray()
        call out % endBlock()
    end do
    
    call out % writeToFile(outputName)

  end subroutine printTransientFlux
!<-- MK 231003 

!--> MK 231005
  !!
  !! Print transient flux
  !!
  subroutine printFissionRate(self,nOut,t)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: nOut
    integer(shortInt), intent(in)                 :: t
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    character(nameLen)                            :: outputName
    integer(shortInt), save                       :: idx, matIdx, i, g !, t
    integer(shortInt)                             :: cIdx
    real(defReal)                                 :: tOut
    real(defReal), save                           :: vol, SigmaF
    type(particleState), save                     :: s
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal), dimension(:), allocatable      :: fiss, fissSTD
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(idx, matIdx, i, mat, matPtr, vol, s, SigmaF, g)

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
        mat    => baseMgNeutronMaterial_CptrCast(matPtr)
        vol    =  self % volume(cIdx)

        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % resultsMap % map(s)

        if ((i > 0) .AND. (matIdx /= 7)) then                               !TODO matIdx /= 7 only for C5G7 - TD3
          do g = 1, self % nG
            SigmaF = real(mat % getFissionXS(g, self % rand),defFlt)        !TODO if Sigma_f changes, it also has to change here! Maybe read Sigma_f with other Sigmas?
            idx = self % nG * (cIdx - 1) + g 
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1) * SigmaF
            ! Is this correct? Also neglects uncertainty in volume - assumed small. Covariance neglected
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
    class(randomRayPhysicsPackage), intent(in) :: self

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
    class(randomRayPhysicsPackage), intent(inout) :: self
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
    if(allocated(self % chi)) deallocate(self % chi)
!--> MK 230302
    if(allocated(self % vel)) deallocate(self % vel)
    if(allocated(self % chiP)) deallocate(self % chiP)
    if(allocated(self % chiD)) deallocate(self % chiD)
    if(allocated(self % lambda)) deallocate(self % lambda)
    if(allocated(self % beta)) deallocate(self % beta)
!<-- MK 230302

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
    self % nTOut       = 0
    self % window      = 0
    self % tstp        = ZERO
    self % initVal     = ZERO
    self % itConv      = 0
    self % nPrintOut   = 0
    self % reduceLines = .false.
!<-- MK 230718

    self % keff        = ZERO
    self % keffScore   = ZERO
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
!--> DNP 230919
    if(allocated(self % DNPScores)) deallocate(self % DNPScores)
    if(allocated(self % DNP)) deallocate(self % DNP)
    if(allocated(self % prevDNP)) deallocate(self % prevDNP)
    if(allocated(self % fluxHist)) deallocate(self % fluxHist)
!<-- DNP 230919
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

end module randomRayPhysicsPackage_class
