module transient_class

  use numPrecision
  use universalVariables,  only : NOT_FOUND
  use genericProcedures,   only : fatalError, numToChar
  use dictionary_class,    only : dictionary
!  use intMap_class,        only : intMap
!  use transient_inter,       only : variation
!  use surfaceFactory_func, only : new_surface_ptr
  use materialMenu_mod,     only : mm_matIdx => matIdx

  implicit none
  private
  
  !!
  !! Storage space for surfaces defined in the geometry
  !!
  !! Sample dictionary input:
  !! transient {
  !!     variation1 {
  !!     type xsVariation;
  !!     crosssection absorption;
  !!     material reflector;
  !!     transientType{type ramp; magnitude 0.1; timeStart 1, timeEnd 10;}
  !!     }
  !! }
  !!
  !! Private Members:
  !!   surfaces -> Array to store pointers to polymorphic surfaces
  !!   idMap -> Map between surface ID and corresponding index
  !!
  !! Interface:
  !!   init    -> Initialise from a dictionary
  !!   getPtr  -> Return pointer to a surface given its index
  !!   getIdx  -> Return index of a surface given its id
  !!   getID   -> Return id of a surface given its idx
  !!   getSize    -> Return the number of surfaces (max surfIdx)
  !!   kill    -> Return to uninitialised state
  !!
  !! NOTE: Becouse surfaces are stored as pointers, calling `kill` is crutial to prevent
  !!   memory leaks. TODO: Add `final` procedure here ?
  !!
  type, public :: transientData
    character(nameLen), dimension(:), allocatable   :: xsArr
    integer(shortInt), dimension(:), allocatable    :: matArr
    character(nameLen), dimension(:), allocatable   :: typeArr
    real(defFlt), dimension(:), allocatable         :: magnitudeArr    
    integer(shortInt), dimension(:), allocatable    :: timeStartArr
    integer(shortInt), dimension(:), allocatable    :: timeEndArr    

  contains
    procedure :: init
    procedure :: kill
  end type transientData

contains

  !!
  !! Load surfaces into shelf
  !!
  !! Args:
  !!   dict [in] -> Dictionary with subdictionaries that contain surface definitions
  !!
  !! Errors:
  !!   fatalError if there are clashes in surface ID
  !!
  subroutine init(self, dict)
    class(transientData), intent(inout)           :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: names
    class(dictionary),pointer                     :: tempDict
    character(nameLen)                            :: matTemp
    real(defReal)                                 :: magnitudeReal
    integer(shortInt)                             :: i
    character(nameLen), dimension(*), parameter :: AVAILABLE_XS = ['fission',&
                                                                   'total  ',&
                                                                   'scatter']
    character(nameLen), dimension(*), parameter :: AVAILABLE_TYPE = ['step',&
                                                                     'ramp']
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter :: Here = 'init (transient_class.f90)'

    ! Get all keys for subdictionaries
    call dict % keys(names, 'dict')

    ! Allocate space
    allocate (self % xsArr(size(names)))
    allocate (self % matArr(size(names)))
    allocate (self % typeArr(size(names)))
    allocate (self % magnitudeArr(size(names)))
    allocate (self % timeStartArr(size(names)))
    allocate (self % timeEndArr(size(names)))
    
    ! Build surfaces
    do i = 1, size(names)
        tempDict => dict % getDictPtr(names(i))
        
        call tempDict % get(self % xsArr(i), 'crosssection')
        IF (.not. (ANY((AVAILABLE_XS) == self % xsArr(i)))) THEN
            call fatalError(Here,'Unrecognised type of a surface: '//trim(self % xsArr(i)))
            print '(A)' , ' AVAILABLE XS: '
            print '(A)' , AVAILABLE_XS
        END IF
        
        call tempDict % get(matTemp, 'material')
        self % matArr(i) = mm_matIdx(matTemp)
        if (self % matArr(i) == NOT_FOUND) call fatalError(Here, 'Unknown material: '//trim(matTemp))
        
        call tempDict % get(self % typeArr(i), 'transientType') !TODO implement error messages if incomplete / wrong data
        IF (.not. (ANY((AVAILABLE_TYPE) == self % typeArr(i)))) THEN
            call fatalError(Here,'Unrecognised variation type: '//trim(self % typeArr(i)))
            print '(A)' , ' AVAILABLE_TYPE: '
            print '(A)' , AVAILABLE_TYPE
        END IF
        
        call tempDict % get(magnitudeReal, 'magnitude')
        self % magnitudeArr(i) = real(magnitudeReal, defFlt)
            
        call tempDict % get(self % timeStartArr(i), 'timeStart') 
        call tempDict % getOrDefault(self % timeEndArr(i), 'timeEnd', NOT_PRESENT)
        IF (self % timeStartArr(i) < 1) THEN
            call fatalError(Here,'Start time is: '//numToChar(self % timeStartArr(i))//&
                ' but must be larger than first timestep = 1')
        END IF
        IF ((self % timeEndArr(i) < self % timeStartArr(i)) .and. .not. (self % timeEndArr(i) == NOT_PRESENT))  THEN
            call fatalError(Here,'End time is: '//numToChar(self % timeEndArr(i))//' but must be larger than start time: '&
                //numToChar(self % timeStartArr(i)))
        END IF

        print *, 'Transient in ', trim(matTemp), ' of type ', trim(self % typeArr(i)), ' with magnitude ', numToChar(magnitudeReal) 

    end do

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(transientData), intent(inout) :: self

    if(allocated(self % xsArr)) deallocate(self % xsArr)
    if(allocated(self % matArr)) deallocate(self % matArr)
    if(allocated(self % typeArr)) deallocate(self % typeArr)
    if(allocated(self % magnitudeArr)) deallocate(self % magnitudeArr)
    if(allocated(self % timeStartArr)) deallocate(self % timeStartArr)
    if(allocated(self % timeEndArr)) deallocate(self % timeEndArr)

  end subroutine kill

end module transient_class
