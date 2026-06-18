#include "w3macros.h"
!/ ------------------------------------------------------------------- /
MODULE W3NMLPRNQMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |           R. Gorman               |
  !/                  |                                   |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         11-Jun-2026 |
  !/                  +-----------------------------------+
  !/
  !/    For updates see subroutines.
  !/
  !  1. Purpose :
  !
  !     Manages namelists from configuration file ww3_prnq.nml for ww3_prnq program
  !
  !/ ------------------------------------------------------------------- /

  ! module defaults
  IMPLICIT NONE

  PUBLIC

  ! grid structure
  TYPE NML_GRID_T
    CHARACTER(256)              :: COORD
    CHARACTER(256)              :: CLOS
    INTEGER                     :: LVLREFT
    INTEGER                     :: LVLMAX
  END TYPE NML_GRID_T

  ! rect structure
  TYPE NML_RECT_T
    INTEGER                     :: NX
    INTEGER                     :: NY
    REAL                        :: SX
    REAL                        :: SY
    REAL                        :: SF
    REAL                        :: X0
    REAL                        :: Y0
    REAL                        :: SF0
  END TYPE NML_RECT_T

  ! file structure
  TYPE NML_FILE_T
    CHARACTER(256)              :: FILENAME
    INTEGER                     :: OUTTYPE
    CHARACTER(256)              :: OUTFILEC
    CHARACTER(256)              :: OUTFILEQ
  END TYPE NML_FILE_T

  ! outvar structure
  TYPE NML_OUTVAR_COUNT_T
    INTEGER                     :: NVAR
  END TYPE NML_OUTVAR_COUNT_T

  TYPE NML_OUTVAR_T
    CHARACTER(256)              :: OUTVARNAME
  END TYPE NML_OUTVAR_T
  
  ! subgrid structure
  TYPE NML_SUBGRID_COUNT_T
    INTEGER                     :: NSUBGRID
    INTEGER                     :: RANK_XYT
  END TYPE NML_SUBGRID_COUNT_T

  TYPE NML_SUBGRID_T
    CHARACTER(256)              :: FILENAME
    INTEGER                     :: IPSG
    CHARACTER(256)              :: XYTDIMNAME(3)
    CHARACTER(256)              :: XYTVARNAME(3)
    INTEGER                     :: ITFIRST
    INTEGER                     :: ITLAST
    INTEGER                     :: ITSTEP
  END TYPE NML_SUBGRID_T

  TYPE NML_SUBGRID_VAR_T
    CHARACTER(256)              :: VARNAME
  END TYPE NML_SUBGRID_VAR_T

  ! miscellaneous
  CHARACTER(256)                :: MSG
  INTEGER                       :: NDSN



CONTAINS
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3NMLPRNQ (NDSI, INFILE, NML_GRID, NML_RECT, &
       NML_FILE, NML_OUTVAR_COUNT, NML_OUTVAR,            &
       NML_SUBGRID_COUNT, NML_SUBGRID, NML_SUBGRID_VAR, IERR)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !
    !  1. Purpose :
    !
    !     Reads all the namelists to define inputs for ww3_prnq
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NDSI              Int.
    !      INFILE            Char.
    !      NML_GRID          Type
    !      NML_RECT          Type
    !      NML_FILE          Type
    !      NML_OUTVAR_COUNT  Type
    !      NML_OUTVAR        Type
    !      NML_SUBGRID_COUNT Type
    !      NML_SUBGRID       Type
    !      NML_SUBGRID_VAR   Type
    !      IERR              Int.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !      READ_GRID_NML
    !      REPORT_GRID_NML
    !      READ_RECT_NML
    !      REPORT_RECT_NML
    !      READ_FILE_NML
    !      REPORT_FILE_NML
    !      READ_OUTVAR_NML
    !      REPORT_OUTVAR_NML
    !      READ_SUBGRID_NML
    !      REPORT_SUBGRID_NML
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      WW3_PRNQ  Prog   N/A     Preprocessor for quadtree data.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE W3ODATMD, ONLY: NDSE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                         :: NDSI
    CHARACTER*(*), INTENT(IN)                   :: INFILE
    TYPE(NML_GRID_T), INTENT(INOUT)             :: NML_GRID
    TYPE(NML_RECT_T), INTENT(INOUT)             :: NML_RECT
    TYPE(NML_FILE_T), INTENT(INOUT)             :: NML_FILE
    TYPE(NML_OUTVAR_COUNT_T), INTENT(INOUT)     :: NML_OUTVAR_COUNT
    TYPE(NML_OUTVAR_T), ALLOCATABLE, INTENT(INOUT) :: NML_OUTVAR(:)
    TYPE(NML_SUBGRID_COUNT_T), INTENT(INOUT)    :: NML_SUBGRID_COUNT
    TYPE(NML_SUBGRID_T), ALLOCATABLE, INTENT(INOUT) :: NML_SUBGRID(:)
    TYPE(NML_SUBGRID_VAR_T), ALLOCATABLE, INTENT(INOUT) :: NML_SUBGRID_VAR(:,:)
    INTEGER, INTENT(OUT)                        :: IERR
#ifdef W3_S
    INTEGER, SAVE                             :: IENT = 0
#endif

    IERR = 0
#ifdef W3_S
    CALL STRACE (IENT, 'W3NMLPRNQ')
#endif

    ! open namelist log file
    NDSN = 3
    OPEN (NDSN, file=TRIM(INFILE)//'.log', form='formatted', iostat=IERR)
    IF (IERR.NE.0) THEN
      WRITE (NDSE,'(A)') 'ERROR: open full nml file '//TRIM(INFILE)//'.log failed'
      RETURN
    END IF

    ! open input file
    OPEN (NDSI, file=TRIM(INFILE), form='formatted', status='old', iostat=IERR)
    IF (IERR.NE.0) THEN
      WRITE (NDSE,'(A)') 'ERROR: open input file '//TRIM(INFILE)//' failed'
      RETURN
    END IF

    ! read grid namelist
    CALL READ_GRID_NML (NDSI, NML_GRID)
    CALL REPORT_GRID_NML (NML_GRID)

    ! read rect namelist
    CALL READ_RECT_NML (NDSI, NML_RECT)
    CALL REPORT_RECT_NML (NML_RECT)

    ! read file namelist
    CALL READ_FILE_NML (NDSI, NML_FILE)
    CALL REPORT_FILE_NML (NML_FILE)

    ! read outvar namelist
    CALL READ_OUTVAR_NML (NDSI, NML_OUTVAR_COUNT, NML_OUTVAR)
    CALL REPORT_OUTVAR_NML (NML_OUTVAR_COUNT, NML_OUTVAR)

    ! read subgrid namelist
    CALL READ_SUBGRID_NML (NDSI, NML_OUTVAR_COUNT%NVAR,        &
                           NML_SUBGRID_COUNT, NML_SUBGRID,     &
                           NML_SUBGRID_VAR)
    CALL REPORT_SUBGRID_NML (NML_OUTVAR_COUNT%NVAR, NML_SUBGRID_COUNT, &
                             NML_SUBGRID, NML_SUBGRID_VAR)

    ! close namelist files
    CLOSE (NDSI)
    CLOSE (NDSN)

  END SUBROUTINE W3NMLPRNQ


  !/ ------------------------------------------------------------------- /


  SUBROUTINE READ_GRID_NML (NDSI, NML_GRID)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NDSI             Int.
    !      NML_GRID         Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLPRNQ Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  :: NDSI
    TYPE(NML_GRID_T), INTENT(INOUT)      :: NML_GRID

    ! locals
    INTEGER                   :: IERR
    TYPE(NML_GRID_T) :: GRID
    NAMELIST /GRID_NML/ GRID
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

    IERR = 0
#ifdef W3_S
    CALL STRACE (IENT, 'READ_GRID_NML')
#endif

    ! set default values for grid structure
    GRID%COORD      = 'unset'
    GRID%CLOS       = 'unset'
    GRID%LVLREFT    = 0
    GRID%LVLMAX     = 0

    ! read grid namelist
    REWIND (NDSI)
    READ (NDSI, nml=GRID_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.NE.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_GRID_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (4)
    END IF

    ! save namelist
    NML_GRID = GRID

  END SUBROUTINE READ_GRID_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE READ_RECT_NML (NDSI, NML_RECT)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NDSI             Int.
    !      NML_RECT         Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLPRNQ Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  :: NDSI
    TYPE(NML_RECT_T), INTENT(INOUT)      :: NML_RECT

    ! locals
    INTEGER                   :: IERR
    TYPE(NML_RECT_T) :: RECT
    NAMELIST /RECT_NML/ RECT
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

    IERR = 0
#ifdef W3_S
    CALL STRACE (IENT, 'READ_RECT_NML')
#endif

    ! set default values for rect structure
    RECT%NX         = 0
    RECT%NY         = 0
    RECT%SX         = 0.
    RECT%SY         = 0.
    RECT%SF         = 1.
    RECT%X0         = 0.
    RECT%Y0         = 0.
    RECT%SF0        = 1.

    ! read rect namelist
    REWIND (NDSI)
    READ (NDSI, nml=RECT_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_RECT_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (5)
    END IF

    ! save namelist
    NML_RECT = RECT

  END SUBROUTINE READ_RECT_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE READ_FILE_NML (NDSI, NML_FILE)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         15-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NDSI             Int.
    !      NML_FILE         Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLPRNQ Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                  :: NDSI
    TYPE(NML_FILE_T), INTENT(INOUT)      :: NML_FILE

    ! locals
    INTEGER                   :: IERR
    TYPE(NML_FILE_T) :: FILES
    NAMELIST /FILE_NML/ FILES
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

    IERR = 0
#ifdef W3_S
    CALL STRACE (IENT, 'READ_FILE_NML')
#endif

    ! set default values for file structure
    FILES%FILENAME  = 'unset'
    FILES%OUTTYPE  =  1
    FILES%OUTFILEC  = 'unset'
    FILES%OUTFILEQ  = 'unset'

    ! read file namelist
    REWIND (NDSI)
    READ (NDSI, nml=FILE_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_FILE_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (6)
    END IF

    ! save namelist
    NML_FILE = FILES

  END SUBROUTINE READ_FILE_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE READ_OUTVAR_NML (NDSI, NML_OUTVAR_COUNT, NML_OUTVAR)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R, Gorman               |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NDSI             Int.
    !      NML_OUTVAR_COUNT Type.
    !      NML_OUTVAR       Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLPRNQ Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: NDSI
    TYPE(NML_OUTVAR_COUNT_T), INTENT(INOUT)   :: NML_OUTVAR_COUNT
    TYPE(NML_OUTVAR_T), ALLOCATABLE, INTENT(INOUT)   :: NML_OUTVAR(:)

    ! locals
    INTEGER                   :: IERR, I
    TYPE(NML_OUTVAR_COUNT_T) :: OUTVAR_COUNT
    NAMELIST /OUTVAR_COUNT_NML/ OUTVAR_COUNT
    TYPE(NML_OUTVAR_T), ALLOCATABLE :: OUTVAR(:)
    NAMELIST /OUTVAR_NML/ OUTVAR
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

    IERR = 0
#ifdef W3_S
    CALL STRACE (IENT, 'READ_OUTVAR_NML')
#endif

    ! set default values for outvar count structure
    OUTVAR_COUNT%NVAR    = 3

    ! read outvar count namelist
    REWIND (NDSI)
    READ (NDSI, nml=OUTVAR_COUNT_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_OUTVAR_COUNT_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (14)
    END IF

    ! allocate the total count of output variables
    ALLOCATE(OUTVAR(OUTVAR_COUNT%NVAR))
    ALLOCATE(NML_OUTVAR(OUTVAR_COUNT%NVAR))

    ! set default values for outvar structure
    IF (OUTVAR_COUNT%NVAR .NE. 0 ) THEN
      DO I=1,OUTVAR_COUNT%NVAR
        IF (I.EQ.1) THEN 
          OUTVAR(I)%OUTVARNAME  = "bed_elevation"
        ELSEIF (I.EQ.2) THEN 
          OUTVAR(I)%OUTVARNAME  = "obs"
        ELSEIF (I.EQ.3) THEN 
          OUTVAR(2)%OUTVARNAME  = "obs_x"
          OUTVAR(3)%OUTVARNAME  = "obs_y"
        ELSE 
          OUTVAR(I)%OUTVARNAME  = "unset"
        END IF
      END DO
    END IF

    ! read outvar namelist
    REWIND (NDSI)
    READ (NDSI, nml=OUTVAR_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_OUTVAR_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (15)
    END IF

    ! save namelist
    NML_OUTVAR_COUNT = OUTVAR_COUNT
    NML_OUTVAR = OUTVAR

  END SUBROUTINE READ_OUTVAR_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE READ_SUBGRID_NML (NDSI, NVAR, NML_SUBGRID_COUNT, &
                               NML_SUBGRID, NML_SUBGRID_VAR)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R, Gorman               |
    !/                  |                                   |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NDSI              Int.
    !      NVAR              Int.
    !      NML_SUBGRID_COUNT Type.
    !      NML_SUBGRID       Type.
    !      NML_SUBGRID_VAR   Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLPRNQ Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                      :: NDSI
    INTEGER, INTENT(IN)                      :: NVAR
    TYPE(NML_SUBGRID_COUNT_T), INTENT(INOUT)   :: NML_SUBGRID_COUNT
    TYPE(NML_SUBGRID_T), ALLOCATABLE, INTENT(INOUT)   :: NML_SUBGRID(:)
    TYPE(NML_SUBGRID_VAR_T), ALLOCATABLE, INTENT(INOUT)   :: NML_SUBGRID_VAR(:,:)

    ! locals
    INTEGER                   :: IERR, I, IVAR
    TYPE(NML_SUBGRID_COUNT_T) :: SUBGRID_COUNT
    NAMELIST /SUBGRID_COUNT_NML/ SUBGRID_COUNT
    TYPE(NML_SUBGRID_T), ALLOCATABLE :: SUBGRID(:)
    NAMELIST /SUBGRID_NML/ SUBGRID
    TYPE(NML_SUBGRID_VAR_T), ALLOCATABLE :: SUBGRID_VAR(:,:)
    NAMELIST /SUBGRID_VAR_NML/ SUBGRID_VAR
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

    IERR = 0
#ifdef W3_S
    CALL STRACE (IENT, 'READ_SUBGRID_NML')
#endif

    ! set default values for subgrid count structure
    SUBGRID_COUNT%NSUBGRID    = 0

    ! read outvar count namelist
    REWIND (NDSI)
    READ (NDSI, nml=SUBGRID_COUNT_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_SUBGRID_COUNT_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (14)
    END IF

    ! allocate the total count of subgrids
    ALLOCATE(SUBGRID(SUBGRID_COUNT%NSUBGRID))
    ALLOCATE(NML_SUBGRID(SUBGRID_COUNT%NSUBGRID))
    ALLOCATE(SUBGRID_VAR(SUBGRID_COUNT%NSUBGRID,NVAR))
    ALLOCATE(NML_SUBGRID_VAR(SUBGRID_COUNT%NSUBGRID,NVAR))

    ! set default values for outvar structure
    IF (SUBGRID_COUNT%NSUBGRID .NE. 0 ) THEN
      DO I=1,SUBGRID_COUNT%NSUBGRID
        SUBGRID(I)%FILENAME  = 'unset'
        SUBGRID(I)%XYTDIMNAME(1) = 'longitude'
        SUBGRID(I)%XYTDIMNAME(2) = 'latitude'
        SUBGRID(I)%XYTDIMNAME(3) = 'time'
        SUBGRID(I)%XYTVARNAME(1) = 'longitude'
        SUBGRID(I)%XYTVARNAME(2) = 'latitude'
        SUBGRID(I)%XYTVARNAME(3) = 'time'
        SUBGRID(I)%IPSG  = 1
        SUBGRID(I)%ITFIRST  = 1
        SUBGRID(I)%ITLAST = 1
        SUBGRID(I)%ITSTEP = 1
        DO IVAR=1,NVAR
          SUBGRID_VAR(I,IVAR)%VARNAME = 'unset'
        END DO
      END DO
    END IF
      

    ! read subgrid namelist
    REWIND (NDSI)
    READ (NDSI, nml=SUBGRID_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_SUBGRID_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (15)
    END IF

    ! read subgrid_var namelist
    REWIND (NDSI)
    READ (NDSI, nml=SUBGRID_VAR_NML, iostat=IERR, iomsg=MSG)
    IF (IERR.GT.0) THEN
      WRITE (NDSE,'(A,/A)') &
           'ERROR: READ_SUBGRID_VAR_NML: namelist read error', &
           'ERROR: '//TRIM(MSG)
      CALL EXTCDE (15)
    END IF

    ! save namelist
    NML_SUBGRID_COUNT = SUBGRID_COUNT
    NML_SUBGRID = SUBGRID
    NML_SUBGRID_VAR = SUBGRID_VAR

  END SUBROUTINE READ_SUBGRID_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE REPORT_GRID_NML (NML_GRID)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NML_GRID  Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLGRID Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    TYPE(NML_GRID_T), INTENT(IN) :: NML_GRID
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

#ifdef W3_S
    CALL STRACE (IENT, 'REPORT_GRID_NML')
#endif

    WRITE (MSG,'(A)') 'GRID % '
    WRITE (NDSN,'(A)')
    WRITE (NDSN,10) TRIM(MSG),'COORD      = ', TRIM(NML_GRID%COORD)
    WRITE (NDSN,10) TRIM(MSG),'CLOS       = ', TRIM(NML_GRID%CLOS)
    WRITE (NDSN,11) TRIM(MSG),'LVLREFT    = ', NML_GRID%LVLREFT
    WRITE (NDSN,11) TRIM(MSG),'LVLMAX     = ', NML_GRID%LVLMAX

10  FORMAT (A,2X,A,A)
11  FORMAT (A,2X,A,I8)
13  FORMAT (A,2X,A,L1)
14  FORMAT (A,2X,A,F8.2)

  END SUBROUTINE REPORT_GRID_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE REPORT_RECT_NML (NML_RECT)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NML_RECT  Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLGRID Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    TYPE(NML_RECT_T), INTENT(IN) :: NML_RECT
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

#ifdef W3_S
    CALL STRACE (IENT, 'REPORT_RECT_NML')
#endif

    WRITE (MSG,'(A)') 'RECT % '
    WRITE (NDSN,'(A)')
    WRITE (NDSN,11) TRIM(MSG),'NX         = ', NML_RECT%NX
    WRITE (NDSN,11) TRIM(MSG),'NY         = ', NML_RECT%NY
    WRITE (NDSN,14) TRIM(MSG),'SX         = ', NML_RECT%SX
    WRITE (NDSN,14) TRIM(MSG),'SY         = ', NML_RECT%SY
    WRITE (NDSN,14) TRIM(MSG),'SF         = ', NML_RECT%SF
    WRITE (NDSN,14) TRIM(MSG),'X0         = ', NML_RECT%X0
    WRITE (NDSN,14) TRIM(MSG),'Y0         = ', NML_RECT%Y0
    WRITE (NDSN,14) TRIM(MSG),'SF0        = ', NML_RECT%SF0

10  FORMAT (A,2X,A,A)
11  FORMAT (A,2X,A,I8)
13  FORMAT (A,2X,A,L1)
14  FORMAT (A,2X,A,F12.2)

  END SUBROUTINE REPORT_RECT_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE REPORT_FILE_NML (NML_FILE)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           R. Gorman               |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         12-Jun-2026 |
    !/                  +-----------------------------------+
    !/
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NML_FILE  Type.
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLGRID Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    TYPE(NML_FILE_T), INTENT(IN) :: NML_FILE
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

#ifdef W3_S
    CALL STRACE (IENT, 'REPORT_FILE_NML')
#endif

    WRITE (MSG,'(A)') 'FILES % '
    WRITE (NDSN,'(A)')
    WRITE (NDSN,10) TRIM(MSG),'FILENAME   = ', NML_FILE%FILENAME
    WRITE (NDSN,11) TRIM(MSG),'OUTTYPE    = ', NML_FILE%OUTTYPE
    WRITE (NDSN,10) TRIM(MSG),'OUTFILEC   = ', NML_FILE%OUTFILEC
    WRITE (NDSN,10) TRIM(MSG),'OUTFILEQ   = ', NML_FILE%OUTFILEQ

10  FORMAT (A,2X,A,A)
11  FORMAT (A,2X,A,I8)
13  FORMAT (A,2X,A,L1)
14  FORMAT (A,2X,A,F12.2)

  END SUBROUTINE REPORT_FILE_NML

  !/ ------------------------------------------------------------------- /

  !/ ------------------------------------------------------------------- /

  SUBROUTINE REPORT_OUTVAR_NML (NML_OUTVAR_COUNT, NML_OUTVAR)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           M. Accensi              |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         15-May-2018 |
    !/                  +-----------------------------------+
    !/
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NML_OUTVAR_COUNT  Type
    !      NML_OUTVAR        Type
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLGRID Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    TYPE(NML_OUTVAR_COUNT_T), INTENT(IN) :: NML_OUTVAR_COUNT
    TYPE(NML_OUTVAR_T), INTENT(IN) :: NML_OUTVAR(NML_OUTVAR_COUNT%NVAR)

    ! locals
    INTEGER              :: I
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

#ifdef W3_S
    CALL STRACE (IENT, 'REPORT_OUTVAR_NML')
#endif

    WRITE (MSG,'(A)') 'OUTVAR_COUNT % '
    WRITE (NDSN,'(A)')
    WRITE (NDSN,11) TRIM(MSG),'NVAR       = ', NML_OUTVAR_COUNT%NVAR

    IF (NML_OUTVAR_COUNT%NVAR .NE. 0) THEN
      DO I=1,NML_OUTVAR_COUNT%NVAR
        WRITE (MSG,'(A,I8,A)') 'OUTVAR(',I,') % '
        WRITE (NDSN,'(A)')
        WRITE (NDSN,10) TRIM(MSG),'OUTVARNAME   = ', NML_OUTVAR(I)%OUTVARNAME
        WRITE (NDSN,'(A)')
      END DO
    END IF

10  FORMAT (A,2X,A,A)
11  FORMAT (A,2X,A,I8)
13  FORMAT (A,2X,A,L1)
14  FORMAT (A,2X,A,F8.2)

  END SUBROUTINE REPORT_OUTVAR_NML

  !/ ------------------------------------------------------------------- /



  !/ ------------------------------------------------------------------- /

  SUBROUTINE REPORT_SUBGRID_NML (NVAR, NML_SUBGRID_COUNT, NML_SUBGRID, &
                                 NML_SUBGRID_VAR)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |           M. Accensi              |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         15-May-2018 |
    !/                  +-----------------------------------+
    !/
    !/
    !  1. Purpose :
    !
    !
    !  2. Method :
    !
    !     See source term routines.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !      NML_SUBGRID        Type
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD SUBROUTINE tracing.
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      TYPE  Module   Description
    !     ----------------------------------------------------------------
    !      W3NMLGRID Subr.   N/A    Namelist configuration routine.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /

#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif

    IMPLICIT NONE

    INTEGER, INTENT(IN)                   :: NVAR
    TYPE(NML_SUBGRID_COUNT_T), INTENT(IN) :: NML_SUBGRID_COUNT
    TYPE(NML_SUBGRID_T), INTENT(IN)       :: NML_SUBGRID(NML_SUBGRID_COUNT%NSUBGRID)
    TYPE(NML_SUBGRID_VAR_T), INTENT(IN)   ::                  &
                         NML_SUBGRID_VAR(NML_SUBGRID_COUNT%NSUBGRID, NVAR)

    ! locals
    INTEGER              :: I, K, NSUBGRID
#ifdef W3_S
    INTEGER, SAVE                           :: IENT = 0
#endif

#ifdef W3_S
    CALL STRACE (IENT, 'REPORT_SUBGRID_NML')
#endif

    NSUBGRID = NML_SUBGRID_COUNT%NSUBGRID
    WRITE (MSG,'(A)') 'SUBGRID_COUNT % '
    WRITE (NDSN,'(A)')
    WRITE (NDSN,11) TRIM(MSG),'NSUBGRID       = ', NSUBGRID

    IF (NSUBGRID .NE. 0) THEN
      DO I=1,NSUBGRID
        WRITE (MSG,'(A,I8,A)') 'SUBGRID(',I,') % '
        WRITE (NDSN,'(A)')
        WRITE (NDSN,10) TRIM(MSG),'FILENAME      = ', NML_SUBGRID(I)%FILENAME
        WRITE (NDSN,11) TRIM(MSG),'IPSG          = ', NML_SUBGRID(I)%IPSG
        WRITE (NDSN,10) TRIM(MSG),'XYTDIMNAME(1) = ', NML_SUBGRID(I)%XYTDIMNAME(1)
        WRITE (NDSN,10) TRIM(MSG),'XYTDIMNAME(2) = ', NML_SUBGRID(I)%XYTDIMNAME(2)
        WRITE (NDSN,10) TRIM(MSG),'XYTDIMNAME(3) = ', NML_SUBGRID(I)%XYTDIMNAME(3)
        WRITE (NDSN,10) TRIM(MSG),'XYTVARNAME(1) = ', NML_SUBGRID(I)%XYTVARNAME(1)
        WRITE (NDSN,10) TRIM(MSG),'XYTVARNAME(2) = ', NML_SUBGRID(I)%XYTVARNAME(2)
        WRITE (NDSN,10) TRIM(MSG),'XYTVARNAME(3) = ', NML_SUBGRID(I)%XYTVARNAME(3)
        WRITE (NDSN,11) TRIM(MSG),'ITFIRST       = ', NML_SUBGRID(I)%ITFIRST
        WRITE (NDSN,11) TRIM(MSG),'ITLAST        = ', NML_SUBGRID(I)%ITLAST
        WRITE (NDSN,11) TRIM(MSG),'ITSTEP        = ', NML_SUBGRID(I)%ITSTEP
        DO K=1,NVAR
          WRITE (NDSN,'(A,A,I4,A)') TRIM(MSG),'VARNAME(',K,') = ',   &
                        NML_SUBGRID_VAR(I,K)%VARNAME
        END DO
        WRITE (NDSN,'(A)')
      END DO
    END IF

10  FORMAT (A,2X,A,A)
11  FORMAT (A,2X,A,I8)
13  FORMAT (A,2X,A,L1)
14  FORMAT (A,2X,A,F8.2)

  END SUBROUTINE REPORT_SUBGRID_NML

  !/ ------------------------------------------------------------------- /



END MODULE W3NMLPRNQMD

!/ ------------------------------------------------------------------- /
