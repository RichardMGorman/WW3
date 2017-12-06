
!/ ------------------------------------------------------------------- /
      PROGRAM W3PRNQ
!/
!         Richard Gorman, NIWA
!         November, 2017:   Origination
!
!/
!  1. Purpose :
!
!     Prepares a file containing data distributed on a quadtree structure, 
!     suitable for adaptive model applications, e.g. Wavewatch, from input 
!     data (in separate netcdf files) on multiple regular grids.
!
!  2. Method :
!
!     Input is taken from a set of netcdf files, each containing data on 
!     a (different) regular grid, i.e. with position (Cartesian X and Y,
!     or latitude and longitude) on rank 2 [NX x NY] or rank 1 (N[X] 
!     and [NY] arrays. Data are interpolated to a quadtree structure,
!     relative to a user-specified reference grid, that provides the
!     maximum refinement level available within each region covered by 
!     an input grid.
!
!     Two output file types are supported:
!         1. Bathymetry and (optionally) subgrid obstruction data
!            output in an ASCII format, along with two separate ASCII
!            files (cell data and quad data) containing full specifications 
!            for the quadtree structure. In this case data are specified
!            at quadtree cell centres at the highest refinement level 
!            locally available from the input grids, but also 
!            "pre-averaged" to all lower refinement levels.
!            These files may be used as input for ww3_grid to set up
!            quadtree-adaptive Wavewatch simulations.
!         2. Multiple time-dependent fields in netcdf format. In this case,
!            the output fields will be at quadtree cell centres at the 
!            highest refinement level locally available from the inputs,
!            but not "pre-averaged" to lower refinement levels, i.e.
!            similar to output from a quadtree-adaptive model. However
!            the single output file will only use one invariant quadtree 
!            structure. Such files can be used as input to ww3_prnc
!            to prepare inputs (e.g. wind, current, water level, ice, ...)
!            for quadtree-adaptive Wavewatch simulations.
!            Note that no time interpolation or checking is done: timing
!            information is taken from the first input file, and the user
!            should specify time indices (first, last, stride) to extract
!            from each input file that ensures consistent timing.
!            
!
!  3. Parameters :
!     N/A
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     GRIDINT    Subr. LOCAL    Compute interpolation indices & weights
!     SGPROPS    Subr. LOCAL    Compute some regular grid properties
!     OPNNCFILR  Subr. LOCAL    Open a regular-grid netcdf file for input
!     OPNNCFILW  Subr. LOCAL    Open a quadtree netcdf file for output
!     CPATTS     Subr. LOCAL    Copy attributes between netcdf files
!     NEXTLN     Subr. W3SERVMD Read past comments up to the next input line
!     EXTCDE     Subr. W3SERVMD Handle fatal errors
!     QA_ALLOC   Subr. qa_utils Allocate a quadtree structure
!     QA_RECT    Subr. qa_utils Set up a quadtree representation (level 0)
!     QA_IOQT    Subr. qa_utils Unformatted IO of a quadtree structure
!     QA_MLG2QT  Subr. qa_utils Derive a quadtree structure and compute 
!                               values for all cells in it for quantities
!                               (e.g. depth, transmission arrays, winds,
!                               currents, ice, ...) provided on a set of
!                               regular input grids 
!     QA_REMUNDEF Subr.qa_utils Remove undefined cells
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
! 10. Source code :
!
!!
       USE NETCDF
       USE QA_UTILS
       USE W3SERVMD, ONLY: NEXTLN, EXTCDE
       IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
       TYPE GRIDSET
          REAL, ALLOCATABLE :: ARRIN(:,:,:)
       END TYPE
       TYPE(GRIDSET), ALLOCATABLE :: SG(:)
       TYPE(QA_TREE)        :: QTREE
       INTEGER              :: UNDEF_TYPE, VALID_TYPE
       INTEGER              :: NX, NY
       INTEGER              :: NX0, NY0
       INTEGER              :: NXREF, NYREF
       INTEGER              :: NXP, NYP
       INTEGER              :: LVLREF, LVLREFT, NSUBGRID
       INTEGER              :: OUTTYPE
       INTEGER              :: NVAR, IVAR, ISG, ISG2, NVMAP
       INTEGER              :: LVLMXB, NINDBG
       INTEGER              :: LVLT, LVLMAX, LV0, LVLDIF
       INTEGER              :: NCUT
       INTEGER              :: I, J, K, II, JJ, ICELL, IST
       INTEGER              :: NSEAMXB, NQMXB
       INTEGER              :: NCIDO
       INTEGER              :: NDSI, NDSO, NDSE, NDSA
       INTEGER              :: ICLOSE
       INTEGER              :: IDUM, IERR, IPAR, IRET
       INTEGER              :: ITMAX, ITOUT, NINDML
       INTEGER              :: RANK_XYT
       INTEGER              :: DIMLEN_OUT(3)
       INTEGER              :: MAP_TYPES(3)
       INTEGER              :: NSVIDO(5)
       INTEGER              :: ISTART2(2), ICOUNT2(2)
       INTEGER              :: ISTART3(3), ICOUNT3(3)
       INTEGER              :: ISCALE
       INTEGER              :: NREDO, IQ, ICELL2
       INTEGER, ALLOCATABLE :: NSUM(:)
       INTEGER, ALLOCATABLE :: NXYTF(:)
       INTEGER, ALLOCATABLE :: IPSG(:)
       INTEGER, ALLOCATABLE :: ITSG(:)
       INTEGER, ALLOCATABLE :: ITFIRST(:), ITLAST(:), ITSTEP(:) 
       INTEGER, ALLOCATABLE :: ATYPE(:) 
       INTEGER, ALLOCATABLE :: NCID(:)
       INTEGER, ALLOCATABLE :: NXF(:), NYF(:), NTF(:)
       INTEGER, ALLOCATABLE :: NXYTDIMID(:,:)
       INTEGER, ALLOCATABLE :: NXYTVARID(:,:)
       INTEGER, ALLOCATABLE :: NVID(:,:)
       INTEGER, ALLOCATABLE :: NXSG(:), NYSG(:)
       INTEGER, ALLOCATABLE :: IQint(:,:)
       INTEGER, ALLOCATABLE :: LVLREL(:)
       INTEGER, ALLOCATABLE :: NEWCELL(:)
       INTEGER, ALLOCATABLE :: MAPML(:)
       INTEGER, ALLOCATABLE :: VNAME(:)
       INTEGER, ALLOCATABLE :: Iint(:,:,:), Jint(:,:,:)
       INTEGER, ALLOCATABLE :: NVIDO(:)
       INTEGER, ALLOCATABLE :: IDDXYT(:)
       INTEGER, ALLOCATABLE :: IDVXYT(:)
       INTEGER, ALLOCATABLE :: IDVAR(:)
       INTEGER, ALLOCATABLE :: RANKF(:)
       INTEGER, ALLOCATABLE :: NXYTVRANK(:,:)
       LOGICAL              :: FLAGLL
       LOGICAL              :: STOPNOW
       LOGICAL, ALLOCATABLE :: HAS_VALID_RANGE(:,:)
       LOGICAL, ALLOCATABLE :: HAS_FILLVAL(:,:)
       REAL                 :: FACTOR
       REAL                 :: SX, SY, X0, Y0, VSC
       REAL                 :: SXP, SYP
       REAL                 :: X0P, Y0P
       REAL                 :: WT
       REAL                 :: SCALE
       REAL                 :: TIME(1)
       REAL                 :: XYPARS(6)
       REAL                 :: X0VAL, Y0VAL
       REAL                 :: ARRVAL
       REAL                 :: VRANGE(2)
       REAL                 :: UNDEF_VAL
       REAL, ALLOCATABLE    :: X0SG(:), Y0SG(:)
       REAL, ALLOCATABLE    :: SXSG(:), SYSG(:)
       REAL, ALLOCATABLE    :: XIN(:,:), YIN(:,:)
       REAL, ALLOCATABLE    :: XSG(:,:), YSG(:,:)
       REAL, ALLOCATABLE    :: ARRQ(:,:)
       REAL, ALLOCATABLE    :: ARANGE(:,:) 
       REAL, ALLOCATABLE    :: Wint(:,:,:)
       REAL, ALLOCATABLE    :: ARMAP(:,:,:)
       REAL, ALLOCATABLE    :: ARRIN(:,:)
       REAL, ALLOCATABLE    :: UNDEF(:)
       REAL, ALLOCATABLE    :: WQint(:,:)
       REAL, ALLOCATABLE    :: X(:), Y(:)
       REAL, ALLOCATABLE    :: XBOUND(:,:), YBOUND(:,:)
       REAL, ALLOCATABLE    :: WTSUM(:)
       REAL, ALLOCATABLE    :: ASUM(:)
       REAL, ALLOCATABLE    :: VALID_RANGE(:,:,:)
       REAL, ALLOCATABLE    :: FILLVAL(:,:)
       CHARACTER            :: COMSTR*1
       CHARACTER            :: CSTRG*4
       CHARACTER*16         :: FMT
       CHARACTER*80         :: DIMNAME_OUT(3)
       CHARACTER*80         :: SVARNAME_OUT(5)
       CHARACTER*80         :: OUTFILE, OUTFILEC, OUTFILEQ
       CHARACTER*80, ALLOCATABLE :: OUTVARNAME(:)
       CHARACTER*80, ALLOCATABLE :: FNAMES(:)
       CHARACTER*80, ALLOCATABLE :: XYTDIMNAME(:,:), XYTVARNAME(:,:)
       CHARACTER*80, ALLOCATABLE :: VARNAME(:,:)
       CHARACTER*80, ALLOCATABLE :: VN(:)
       CHARACTER*80, ALLOCATABLE :: XYTDN(:)
       CHARACTER*80, ALLOCATABLE :: XYTVN(:)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 0. Constants
!
       NDSA = 11
       NDSI = 10
       NDSO = 6
       NDSE = 6
       UNDEF_TYPE = 0
       VALID_TYPE = 1
!
! 1. User input
!
      OPEN (NDSI,FILE='ww3_prnq.inp',STATUS='OLD',                     &
            ERR=2000,IOSTAT=IERR)
      WRITE (NDSO,900)
      WRITE (NDSO,3000)
      READ (NDSI,'(A)',END=2001,ERR=2002) COMSTR
      IF (COMSTR.EQ.' ') COMSTR = '$'
      WRITE (NDSO,901) COMSTR
!
! 1.a Type of grid (lat/lon?) & closure type
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) FLAGLL, CSTRG
!
       IF ( FLAGLL ) THEN
          FACTOR = 1.
          WRITE (NDSO,3001) 'spherical'
       ELSE
          FACTOR = 1.E-3
          WRITE (NDSO,3001) 'Cartesian'
       END IF
       ICLOSE = 0
       SELECT CASE (TRIM(CSTRG))
          CASE ('NONE')
             ICLOSE = 0
             WRITE (NDSO,3002) 'none'
          CASE ('SMPL', 'TRUE', 't', 'T','.tru','.TRU')
             ICLOSE = 1
             WRITE (NDSO,3002) 'simple'
       END SELECT
!
! 1.b  reference grid dimensions
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) NX, NY
       NX     = MAX ( 3 , NX )
       NY     = MAX ( 3 , NY )
       WRITE (NDSO,3003) NX, NY
!
! 1.c  reference grid spacings
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) SX, SY, VSC
       VSC    = MAX ( 1.E-7 , VSC )
       SX     = SX / VSC
       SY     = SY / VSC
       SX     = MAX ( 1.E-7 , SX )
       SY     = MAX ( 1.E-7 , SY )
       IF ( ICLOSE.EQ.1 ) SX = 360. / REAL(NX)
!
! 1.d  reference grid origin
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) X0, Y0, VSC
       VSC    = MAX ( 1.E-7 , VSC )
       X0     = X0 / VSC
       Y0     = Y0 / VSC
       IF ( FLAGLL ) THEN
              WRITE (NDSO,3004) FACTOR*SX, FACTOR*SY,         &
                     FACTOR*X0, FACTOR*(X0+REAL(NX-1)*SX),    &
                     FACTOR*Y0, FACTOR*(Y0+REAL(NY-1)*SY)
       ELSE
              WRITE (NDSO,3005) FACTOR*SX, FACTOR*SY,         &
                     FACTOR*X0, FACTOR*(X0+REAL(NX-1)*SX),    &
                     FACTOR*Y0, FACTOR*(Y0+REAL(NY-1)*SY)
       END IF
!
! 1.e  Max. number of levels to coarsen and refine the reference grid,
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) LVLREFT, LVLMAX
       WRITE (NDSO,3006)  LVLREFT, LVLMAX
!
! 1.f  Data output file
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) OUTFILE
       WRITE (NDSO,3007) TRIM(OUTFILE)
!
! 1.g  Output type: 1. ASCII depth [ & obstruction] data,
!                      no time dependence
!                   2. other data (netcdf format)
!                      may be time dependent
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) OUTTYPE
!
! 1.h  Number of variables to process
!         
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) NVAR
       IF ( OUTTYPE.EQ.1 ) NVAR = MIN(NVAR,3)
!
! 1.i  Depending on OUTTYPE, either ...
!         
       ALLOCATE ( OUTVARNAME(NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'OUTVARNAME alloc. error', IERR
       IF ( OUTTYPE.EQ.1 ) THEN
!         ... additional output files:
!         1.i.1 Quadtree cell file
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) OUTFILEC
!         1.i.2 Quadtree quad file
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) OUTFILEQ
          !
          WRITE (NDSO,3008) TRIM(OUTFILEC), TRIM(OUTFILEQ)
          !
          ! Names for output variables are fixed:
          IF (NVAR.GE.1) OUTVARNAME(1) = "bed_elevation"
          IF (NVAR.EQ.2) OUTVARNAME(2) = "obs"
          IF (NVAR.GE.3) THEN
             OUTVARNAME(2) = "obs_x"
             OUTVARNAME(3) = "obs_y"
          END IF
          !
          ! Don't look for a time dimension or variable, just X&Y:
          RANK_XYT = 2
       ELSE
          WRITE (NDSO,3009) 
!         ... or output variables
          DO IVAR=1,NVAR
             CALL NEXTLN ( COMSTR , NDSI , NDSE )
             READ (NDSI,*,END=2001,ERR=2002) OUTVARNAME(IVAR)
          END DO
          !
          ! Look for X, Y and time dimensions and variables:
          RANK_XYT = 3
       END IF
!
       DO IVAR=1,NVAR
          WRITE (NDSO,3010) IVAR, TRIM(OUTVARNAME(IVAR))
       END DO
!
! 1.j Read subgrid information         
!      number of input grids/files
!
       CALL NEXTLN ( COMSTR , NDSI , NDSE )
       READ (NDSI,*,END=2001,ERR=2002) NSUBGRID
       WRITE (NDSO,3011)  NSUBGRID
!
       ALLOCATE ( FNAMES(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'FNAMES alloc. error', IERR
       ALLOCATE ( IPSG(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'IPSG alloc. error', IERR
       ALLOCATE ( XYTDIMNAME(NSUBGRID,RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XYTDIMNAME alloc. error', IERR
       ALLOCATE ( XYTVARNAME(NSUBGRID,RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XYTVARNAME alloc. error', IERR
       ALLOCATE ( VARNAME(NSUBGRID,NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'VARNAME alloc. error', IERR
       ALLOCATE ( ITFIRST(NSUBGRID), ITLAST(NSUBGRID),                 &
                  ITSTEP(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ITFIRST, etc. alloc. error', IERR
!
       WRITE (NDSO,3011) 
!
       DO ISG=1,NSUBGRID
          !
          ! 1.j.1  Subgrid index (ignored), 
          !        Index of the subgrid which is the parent of this one
          !          (i.e. a more extensive, lower resolution subgrid
          !           containing it)
          !        Name of the netcdf file
          CALL NEXTLN ( COMSTR , NDSI , NDSE )
          READ (NDSI,*,END=2001,ERR=2002) IDUM, IPSG(ISG), FNAMES(ISG)
          WRITE (NDSO,3012)  ISG, IPSG(ISG), TRIM(FNAMES(ISG))
          !
          ! 1.j.2  Name of the X/lon, Y/lat [, time] dimensions and
          !        variables in this file 
          DO K=1,RANK_XYT
             CALL NEXTLN ( COMSTR , NDSI , NDSE )
             READ (NDSI,*,END=2001,ERR=2002) XYTDIMNAME(ISG,K)
             CALL NEXTLN ( COMSTR , NDSI , NDSE )
             READ (NDSI,*,END=2001,ERR=2002) XYTVARNAME(ISG,K)
             WRITE (NDSO,3013)  ISG, K, TRIM(XYTDIMNAME(ISG,K)),       &
                                TRIM(XYTVARNAME(ISG,K))
          END DO
          !
          IF ( OUTTYPE.EQ.2 ) THEN
          ! 1.j.3  Time indices to use: first, last, stride
          !        (ignored for OUTTYPE = 1)
             CALL NEXTLN ( COMSTR , NDSI , NDSE )
             READ (NDSI,*,END=2001,ERR=2002) ITFIRST(ISG),             &
                                             ITLAST(ISG), ITSTEP(ISG)
             WRITE (NDSO,3014)  ISG, ITFIRST(ISG), ITLAST(ISG),        &
                                     ITSTEP(ISG)
          END IF
          !
          ! 1.j.4  Name, in this file, corresponding to each 
          !        requested output variable
          DO IVAR=1,NVAR
             CALL NEXTLN ( COMSTR , NDSI , NDSE )
             READ (NDSI,*,END=2001,ERR=2002) VARNAME(ISG,IVAR)
             WRITE (NDSO,3015)  ISG, IVAR, TRIM(VARNAME(ISG,IVAR))
          END DO
       END DO  ! loop over subgrids 
!
! 2. Some preliminaries:
!
! 2.a determine the number of cells [NX0,NY0] in the level-0 grid, and the 
!     refinement level LVLREF of the reference grid, which may be less than 
!     requested, so that NXREF = NX0*2**LVLREF, NYREF = NY0*2**LVLREF, :
!
       NXREF = NX
       NYREF = NY
       IF ( LVLREFT.LT.0 ) LVLREFT = 1000 
       NX0 = NX
       NY0 = NY
       LVLREF = 0
       DO WHILE( MOD(NX0,2).EQ.0 .AND. MOD(NY0,2).EQ.0 .AND.       &
                 LVLREF.LT.LVLREFT )
          NX0 = NX0/2
          NY0 = NY0/2
          LVLREF = LVLREF + 1
       END DO
       WRITE (NDSO,3006)  LVLREF, LVLMAX
!
! 2.b Settings for the mapping from regular subgrids to the quadtree 
!  
!    2.b.1 Number of variables to be mapped from regular
!          grids to the quadtree structure:
!
       IF ( OUTTYPE.EQ.1 ) THEN
          ! mapping is applied to the output variables directly
          NVMAP = NVAR
       ELSE
          ! mapping is applied to the interpolation indices and weights
          NVMAP = 13
       END IF
!
!    2.b.2 Variable type for mapping/averaging to lower refinement level,
!          and valid range:
!
       ALLOCATE ( ATYPE(NVMAP), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ATYPE alloc. error', IERR
       ALLOCATE ( ARANGE(NVMAP,2), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARANGE alloc. error', IERR
!
       IF ( OUTTYPE.EQ.1 ) THEN
          ! simple averaging is applied to bed elevation 
          IF (NVAR.GE.1) ATYPE(1) = 1
          ! special averaging is applied to X and Y obstruction
          IF (NVAR.GE.2) ATYPE(2) = 3
          IF (NVAR.GE.3) ATYPE(3) = 4    
       ELSE
          ! rather than averaging, values from a single (wet) cell
          ! are copied to the lower refinement level. This is just to 
          ! flag to lower levels whether there is at least one higher
          ! -level wet cell that can be averaged, though the resulting 
          ! indices and weights won't be the best way to do that
          ATYPE = 5
       END IF
!
! 3. Second loop over subgrids to determine array parameters
!

       ALLOCATE ( NCID(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NSUBGRID alloc. error', IERR
       ALLOCATE ( NXF(NSUBGRID), NYF(NSUBGRID), NTF(NSUBGRID),         &
                                                    stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NXF, NYF, NTF alloc. error', IERR
       ALLOCATE ( NXYTDIMID(NSUBGRID,RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NXYTDIMID alloc. error', IERR
       ALLOCATE ( NXYTVARID(NSUBGRID,RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NXYTVARID alloc. error', IERR
       ALLOCATE ( NXYTVRANK(NSUBGRID,RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NXYTVRANK alloc. error', IERR
       ALLOCATE ( NXYTF(RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NXYTF alloc. error', IERR
       ALLOCATE ( RANKF(RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'RANKF alloc. error', IERR
       ALLOCATE ( NVID(NSUBGRID,NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NVID alloc. error', IERR
       ALLOCATE ( NXSG(NSUBGRID), NYSG(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NXSG, NYSG alloc. error', IERR
       ALLOCATE ( X0SG(NSUBGRID), Y0SG(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'X0SG, Y0SG alloc. error', IERR
       ALLOCATE ( SXSG(NSUBGRID), SYSG(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'SXSG, SYSG alloc. error', IERR
       ALLOCATE ( LVLREL(NSUBGRID), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'LVLREL alloc. error', IERR
!
       !write(*,*) 'Allocating with NVAR, RANK_XYT = ', NVAR, RANK_XYT
       ALLOCATE ( VN(NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'VN alloc. error', IERR
       ALLOCATE ( XYTDN(RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XYTDN alloc. error', IERR
       ALLOCATE ( XYTVN(RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XYTVN alloc. error', IERR
       ALLOCATE ( IDDXYT(RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'IDDXYT alloc. error', IERR
       ALLOCATE ( IDVXYT(RANK_XYT), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'IDVXYT alloc. error', IERR
       ALLOCATE ( IDVAR(NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'IDVAR alloc. error', IERR
!
       NXF = 0
       NYF = 0
       NTF = 0
       LVLMXB = 0
       NINDBG = 0
       WRITE (NDSO,3011)  NSUBGRID
       DO ISG=1,NSUBGRID
          !
          ! Open the netcdf file, and read the dimensions, and
          ! ids for time, x, y and the input variables
          !
          XYTDN = XYTDIMNAME(ISG,:)
          XYTVN = XYTVARNAME(ISG,:)
          VN = VARNAME(ISG,:)
          CALL OPNNCFILR( FNAMES(ISG), RANK_XYT, RANK_XYT, NVAR,       &
                          XYTDN, XYTVN, VN, NCID(ISG), NXYTF,          &
                          IDDXYT, IDVXYT, RANKF, IDVAR, NDSE )
          NXYTDIMID(ISG,:) = IDDXYT
          NXYTVARID(ISG,:) = IDVXYT
          NVID(ISG,:) = IDVAR
          IF (RANK_XYT.GE.1) THEN
             NXF(ISG) = NXYTF(1)
             NXYTVRANK(ISG,1) = RANKF(1)
          END IF
          IF (RANK_XYT.GE.2) THEN
             NYF(ISG) = NXYTF(2)
             NXYTVRANK(ISG,2) = RANKF(2)
          END IF
          IF (RANK_XYT.GE.3) THEN
             NTF(ISG) = NXYTF(3)
             NXYTVRANK(ISG,3) = RANKF(3)
          END IF
          !
          WRITE (NDSO,4001)  ISG, NXF(ISG), NYF(ISG), NTF(ISG)
          !
          ! Allocate spatial arrays for the input subgrid
          !
          ALLOCATE ( XIN(NXF(ISG),NYF(ISG)), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XIN alloc. error', IERR
          ALLOCATE ( YIN(NXF(ISG),NYF(ISG)), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'YIN alloc. error', IERR
          !
          ! Read in spatial arrays from the netcdf file
          !
          !write(*,*) 'Reading X from grid, NCID = ', NCID(ISG),        &
          !               NXYTVARID(ISG,1)
          IF ( NXYTVRANK(ISG,1).EQ.1 ) THEN
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,1),         &
                          XIN(:,1), start=(/1/), count=(/NXF(ISG)/) )
             CALL CHECK_ERR(IRET,NDSE)
             DO J=2,NYF(ISG)
                XIN(:,J) = XIN(:,1)
             END DO
          ELSEIF ( NXYTVRANK(ISG,1).EQ.2 ) THEN
             ISTART2(1) = 1
             ISTART2(2) = 1
             ICOUNT2(1) = NXF(ISG)
             ICOUNT2(2) = NYF(ISG)
             !write(*,*) 'ISTART2 = ', ISTART2
             !write(*,*) 'ICOUNT2 = ', ICOUNT2
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,1), XIN,    &
                               start=ISTART2, count=ICOUNT2 )
             CALL CHECK_ERR(IRET,NDSE)
          ELSE
             WRITE(NDSE,*) 'unable to read X of rank ',NXYTVRANK(ISG,1)
             CALL EXTCDE(1)
          END IF
          !write(*,*) 'Reading Y from grid, NCID = ', NCID(ISG),        &
          !            NXYTVARID(ISG,2)
          IF ( NXYTVRANK(ISG,2).EQ.1 ) THEN
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,2),         &
                          YIN(1,:), start=(/1/), count=(/NYF(ISG)/) )
             CALL CHECK_ERR(IRET,NDSE)
             DO I=2,NXF(ISG)
                YIN(I,:) = YIN(1,:)
             END DO
          ELSEIF ( NXYTVRANK(ISG,2).EQ.2 ) THEN
             ISTART2(1) = 1
             ISTART2(2) = 1
             ICOUNT2(1) = NXF(ISG)
             ICOUNT2(2) = NYF(ISG)
             !write(*,*) 'ISTART2 = ', ISTART2
             !write(*,*) 'ICOUNT2 = ', ICOUNT2
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,2), YIN,    &
                               start=ISTART2, count=ICOUNT2 )
             CALL CHECK_ERR(IRET,NDSE)
          ELSE
             WRITE(NDSE,*) 'unable to read Y of rank ',NXYTVRANK(ISG,2)
             CALL EXTCDE(1)
          END IF
          !
          ! Evaluate some properties of the input grid:
          !     XYPARS(1:2) = minimum cell X,Y-extent
          !     XYPARS(3:4) = min. & max cell X value
          !     XYPARS(5:6) = min. & max cell Y value
          !
          CALL SGPROPS( NXF(ISG), NYF(ISG), XIN, YIN, XYPARS )
          !
          WRITE (NDSO,4002)  ISG, (XYPARS(K), K=1,6)
          !
          ! Deallocate the input arrays, ready to reuse for the next grid
          !
          IF ( ALLOCATED(XIN) ) THEN
             DEALLOCATE ( XIN, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XIN dealloc. error', IERR
          END IF
          IF ( ALLOCATED(YIN) ) THEN
             DEALLOCATE ( YIN, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'YIN dealloc. error', IERR
          END IF
          !
          ! Select the matching quadtree subgrid, that covers the input grid
          ! at similar (or slightly higher) refinement level
          !
          DO LVLT=0,LVLMAX
             SXP = SX*2.**(LVLREF-LVLT)
             SYP = SY*2.**(LVLREF-LVLT)
             LVLREL(ISG) = LVLT - LVLREF
             IF ( SXP.LT.1.01*XYPARS(1) .AND. SYP.LT.1.01*XYPARS(2) ) EXIT
          END DO
          !
          WRITE (NDSO,4003)  ISG, LVLREL(ISG)
          !
          ! Update the maximum refinement level
          !
          LVLMXB = MAX(LVLMXB, LVLREF + LVLREL(ISG))
          !
          ! Size of the parent subgrid
          !
          IF ( ISG.EQ.1 ) THEN
             !
             ! level 0 relative to the reference level
             !
             LV0 = -LVLREF
             ! 
             ! The first subgrid's parent is the full level-0 grid:
             ! Number of cells in the full parent subgrid
             ! 
             NXP = NX0
             NYP = NY0
             ! 
             ! step size on the parent subgrid relative to the reference grid:
             ! 
             SCALE = 2.**LVLREF
             ! 
             ! Origin of the parent subgrid:
             ! 
             X0P = X0 - 0.5*SX*(1.-SCALE)
             Y0P = Y0 - 0.5*SY*(1.-SCALE)
             ! 
             ! level of this subgrid relative to its parent:
             ! 
             LVLDIF = LVLREL(ISG) + LVLREF
             ! 
          ELSE
             ! 
             ! Index of this subgrid's parent  
             ! 
             IPAR = IPSG(ISG)
             ! 
             ! Number of cells in the full parent subgrid
             ! 
             NXP = NXSG(IPAR)
             NYP = NYSG(IPAR)
             ! 
             ! step size on the parent subgrid relative to the reference grid:
             ! 
             SCALE = 2.**(-LVLREL(IPAR))
             ! 
             ! Origin of the parent subgrid:
             ! 
             X0P = X0SG(IPAR)
             Y0P = Y0SG(IPAR)
             ! 
             ! level of this subgrid relative to its parent:
             ! 
             LVLDIF = LVLREL(ISG)-LVLREL(IPAR)
             ! 
          END IF
          !
          ! Step size of the parent subgrid:
          !
          SXP = SX*SCALE
          SYP = SY*SCALE
          !
          ! Origin of the full parent subgrid:
          ! ... restricted to lie within the input grid 
          ! trim cells off the west:
          !
          NCUT = CEILING( ( XYPARS(3) - X0P )/SXP )
          IF ( NCUT.GT.0 ) THEN
             X0P =  X0P + NCUT*SXP
             NXP = NXP - NCUT
          END IF
          !
          ! trim cells off the east:
          !
          NCUT = CEILING( ( X0P + (NXP-1)*SXP - XYPARS(4) )/SXP ) 
          IF ( NCUT.GT.0 ) NXP = NXP - NCUT
          !
          ! trim cells off the south:
          !
          NCUT = CEILING( ( XYPARS(5) - Y0P )/SYP )
          IF ( NCUT.GT.0 ) THEN
             Y0P =  Y0P + NCUT*SYP
             NYP = NYP - NCUT
          END IF
          !
          ! trim cells off the north:
          !
          NCUT = CEILING( ( Y0P + (NYP-1)*SYP - XYPARS(6) )/SYP )
          IF ( NCUT.GT.0 ) NYP = NYP - NCUT
          !
          ! step sizes, origin, number of cells in the the reduced subgrid:
          !
          ISCALE = 2**LVLDIF
          SCALE = 1./FLOAT(ISCALE)
          SXSG(ISG) = SCALE*SXP
          SYSG(ISG) = SCALE*SYP
          X0SG(ISG) = X0P - 0.5*SXP*(1.-SCALE)
          Y0SG(ISG) = Y0P - 0.5*SYP*(1.-SCALE)
          NXSG(ISG) = NXP*ISCALE
          NYSG(ISG) = NYP*ISCALE
          WRITE (NDSO,4004) ISG, 'nearest regular subgrid',            &
                            NXSG(ISG), NYSG(ISG),                      &
                            SXSG(ISG), SYSG(ISG), X0SG(ISG), Y0SG(ISG)
          WRITE (NDSO,4004) ISG, 'parent', NXP, NYP, SXP, SYP, X0P, Y0P
          !
          ! Cumulative number of cells in the (local) multilevel grid
          !
          IF ( ISG.EQ.1 ) THEN
             ! include cells from LVLREL(IPAR) to LVLREL(ISG)
             NINDBG = NXP*NYP*NINT( ( 4.**(LVLDIF+1) - 1.)/3. ) 
          ELSE
             ! include cells from LVLREL(IPAR)+1 to LVLREL(ISG)
             NINDBG = NINDBG + CEILING( NXP*NYP*(4.**LVLDIF - 1.)*4./3. )
          END IF
          !
       END DO  ! loop over subgrids 
!
! 4. Preparations for the quadtree structure
!     4.a  Allocate arrays
!
       IF ( OUTTYPE.EQ.1 ) THEN
          ALLOCATE ( ARRQ(NINDBG,NVAR), stat=IERR )
       ELSE
          ALLOCATE ( ARRQ(NINDBG,12), stat=IERR )
       END IF
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARRQ alloc. error', IERR
       ALLOCATE ( NEWCELL(NINDBG), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NEWCELL alloc. error', IERR
       NINDML = NINT( NX0*NY0*(4.**(LVLMXB+1) - 1)/3. )
       ALLOCATE ( MAPML(NINDML), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'MAPML alloc. error', IERR
       ALLOCATE ( NSUM(NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'NSUM alloc. error', IERR
       ALLOCATE ( ASUM(NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ASUM alloc. error', IERR
       ALLOCATE ( WTSUM(NVAR), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'WTSUM alloc. error', IERR
!
       WRITE (NDSO,4005) NINDBG, NINDML
!
! Initialise MAPML
!
       MAP_TYPES(1) =  -999   ! flag for cells with type yet to be determined
       MAP_TYPES(2) =  VALID_TYPE     ! flag for default wet cells 
       MAP_TYPES(3) =  2      ! flag for cells not to be refined
       MAPML = MAP_TYPES(1)
!
! Set up a quadtree structure for the level zero grid
!
       CALL QA_ALLOC ( QTREE, NINDBG, NINDBG, UNDEF_TYPE )
       CALL QA_RECT ( NX0, NY0, LVLREF, ICLOSE.EQ.1, VALID_TYPE,    &
                         UNDEF_TYPE, QTREE, ierr=IERR, ndse=NDSE )
       IF ( IERR .NE. 0 ) CALL EXTCDE ( IERR )
       !write(*,*) 'after QA_RECT:'
       !write(*,*) 'NQUAD:',QTREE%NQUAD
       !write(*,*) 'NCELL:',QTREE%NCELL
       !write(*,*) 'LVLREF:',QTREE%LVLREF
       !write(*,*) 'LVLMAX:',QTREE%LVLMAX
       !write(*,*) 'LVLHI:',QTREE%LVLHI
!
!  Allow cell indices to be retained when they are refined?
!
       IF ( OUTTYPE.EQ.1 ) THEN
          ! The bathymetry quadtree has pre-averaged data at multiple 
          ! levels
          QTREE%KEEP_REF = .TRUE.
       ELSE
          ! Other quadtrees only have data at the highest local refinemnt
          ! level
          QTREE%KEEP_REF = .FALSE.
       END IF
!  Allow for the quadtree to adapt?
       QTREE%DYNAMIC = .FALSE.
       QTREE%LVLMAX = LVLMXB
!
! 5. Third loop over subgrids to process the spatial data
!        
       ALLOCATE ( VNAME(NVAR), stat=IERR ) 
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'VNAME alloc. error', IERR
       ALLOCATE ( UNDEF(NVAR), stat=IERR ) 
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'UNDEF alloc. error', IERR
       ALLOCATE ( FILLVAL(NSUBGRID,NVAR), stat=IERR ) 
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'FILLVAL alloc. error', IERR
       ALLOCATE ( VALID_RANGE(NSUBGRID,NVAR,2), stat=IERR ) 
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'VALID_RANGE alloc. error', IERR
       ALLOCATE ( HAS_FILLVAL(NSUBGRID,NVAR), stat=IERR ) 
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'HAS_FILLVAL alloc. error', IERR
       ALLOCATE ( HAS_VALID_RANGE(NSUBGRID,NVAR), stat=IERR ) 
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'HAS_VALID_RANGE alloc. error', IERR
       HAS_FILLVAL = .FALSE.
       HAS_VALID_RANGE = .FALSE.
       DO ISG=1,NSUBGRID
          !
          ! Allocate spatial arrays for the input subgrid
          !
          ALLOCATE ( XIN(NXF(ISG),NYF(ISG)), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XIN alloc. error', IERR
          ALLOCATE ( YIN(NXF(ISG),NYF(ISG)), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'YIN alloc. error', IERR
          !
          ! Read in spatial arrays from the netcdf file
          !
          !write(*,*) 'Reading X from grid, NCID = ', NCID(ISG),        &
          !               NXYTVARID(ISG,1)
          IF ( NXYTVRANK(ISG,1).EQ.1 ) THEN
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,1),         &
                          XIN(:,1), start=(/1/), count=(/NXF(ISG)/) )
             CALL CHECK_ERR(IRET,NDSE)
             DO J=2,NYF(ISG)
                XIN(:,J) = XIN(:,1)
             END DO
          ELSEIF ( NXYTVRANK(ISG,1).EQ.2 ) THEN
             ISTART2(1) = 1
             ISTART2(2) = 1
             ICOUNT2(1) = NXF(ISG)
             ICOUNT2(2) = NYF(ISG)
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,1), XIN,    &
                               start=ISTART2, count=ICOUNT2 )
             CALL CHECK_ERR(IRET,NDSE)
          ELSE
             WRITE(NDSE,*) 'unable to read X of rank ',NXYTVRANK(ISG,1)
             CALL EXTCDE(1)
          END IF
          !write(*,*) 'Reading Y from grid, NCID = ', NCID(ISG),        &
          !            NXYTVARID(ISG,2)
          IF ( NXYTVRANK(ISG,2).EQ.1 ) THEN
             IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,2),         &
                          YIN(1,:), start=(/1/), count=(/NYF(ISG)/) )
             CALL CHECK_ERR(IRET,NDSE)
             DO I=2,NXF(ISG)
                YIN(I,:) = YIN(1,:)
             END DO
          ELSEIF ( NXYTVRANK(ISG,2).EQ.2 ) THEN
             ISTART2(1) = 1
             ISTART2(2) = 1
             ICOUNT2(1) = NXF(ISG)
             ICOUNT2(2) = NYF(ISG)
            IRET = nf90_get_var( NCID(ISG), NXYTVARID(ISG,2), YIN,    &
                               start=ISTART2, count=ICOUNT2 )
             CALL CHECK_ERR(IRET,NDSE)
          ELSE
             WRITE(NDSE,*) 'unable to read Y of rank ',NXYTVRANK(ISG,2)
             CALL EXTCDE(1)
          END IF
          !
          ! Assign spatial arrays for the single-level subgrid
          !
          ALLOCATE ( XSG(NXSG(ISG),NYSG(ISG)), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XSG alloc. error', IERR
          ALLOCATE ( YSG(NXSG(ISG),NYSG(ISG)), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'YSG alloc. error', IERR
          DO J=1,NYSG(ISG)
             DO I=1,NXSG(ISG)
                XSG(I,J) = X0SG(ISG) + SXSG(ISG)*(I-1)
                YSG(I,J) = Y0SG(ISG) + SYSG(ISG)*(J-1)
             END DO          
          END DO          
          !
          ! Compute indices and weights for interpolation from the input grid
          ! to the selected single-level regular subgrid
          ! Indices Iint, Jint and weights Wint are computed such that, 
          ! if an array A1 is defined on the input grid, it can be 
          ! interpolated to an array A2 defined on the regular subgrid
          !       with
          ! A2(i,j) = sum(k=1:4) { Wint(i,j,k)*A1(Iint(i,j,k),Jint(i,j,k)) }
          !
          ALLOCATE ( Iint(NXSG(ISG),NYSG(ISG),4), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'Iint alloc. error', IERR
          ALLOCATE ( Jint(NXSG(ISG),NYSG(ISG),4), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'Jint alloc. error', IERR
          ALLOCATE ( Wint(NXSG(ISG),NYSG(ISG),4), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'Wint alloc. error', IERR
          CALL gridint( NXF(ISG), NYF(ISG), XIN, YIN,                  &
                        NXSG(ISG), NYSG(ISG), XSG, YSG,                &
                        Iint, Jint, Wint )
          WRITE (NDSO,4006) ISG, NXF(ISG), NYF(ISG),                   &
                            NXSG(ISG), NYSG(ISG)
          !
          ! Deallocate the input arrays
          !
          IF ( ALLOCATED(XIN) ) THEN
             DEALLOCATE ( XIN, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XIN dealloc. error', IERR
          END IF
          IF ( ALLOCATED(YIN) ) THEN
             DEALLOCATE ( YIN, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'YIN dealloc. error', IERR
          END IF
          !
          ALLOCATE ( ARMAP(NXSG(ISG),NYSG(ISG),NVMAP), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARMAP alloc. error', IERR
          !
          ! Determine valid range and fill value for each variable
          DO IVAR=1,NVAR
             !  Read the fill value and valid range, 
             IRET = nf90_get_att( NCID(ISG), NVID(ISG,IVAR),           &
                                     '_FillValue', UNDEF_VAL )
             IF ( IRET /= nf90_noerr ) THEN
                ! No fill value in this input file
                IF ( ISG.GT.1 ) THEN
                   ! First option: use value from the first input file
                   UNDEF_VAL = FILLVAL(1,IVAR)
                ELSE IF ( OUTTYPE.EQ.1 ) THEN
                   ! Second option: use defaults for bathymetry, obstruction
                   IF ( IVAR.EQ.1 ) THEN
                      UNDEF_VAL = -999999.
                   ELSE
                      UNDEF_VAL = -9.
                   END IF
                ELSE 
                   ! Third option: arbitrary default
                   UNDEF_VAL = -999999.
                END IF
             ELSE
                HAS_FILLVAL(ISG,IVAR) = .TRUE.
             END IF 
             FILLVAL(ISG,IVAR) = UNDEF_VAL
             IRET = nf90_get_att( NCID(ISG), NVID(ISG,IVAR),     &
                                     'Valid_Range', VRANGE )
             IF ( IRET /= nf90_noerr ) THEN
                ! No valid range in this input file
                IF ( ISG.GT.1 ) THEN
                   ! First option: use value from the first input file
                   VRANGE = VALID_RANGE(1,IVAR,:)
                ELSE IF ( OUTTYPE.EQ.1 ) THEN
                   IF ( IVAR.EQ.1 ) THEN
                      VRANGE(1) = -99999.
                      VRANGE(2) = 99999.
                       !  ARANGE(IVAR,1) = -99999.
                       !  ARANGE(IVAR,2) = 99999.
                   ELSE
                      VRANGE(1) = -0.001
                      VRANGE(2) = 1.001
                      !   ARANGE(IVAR,1) = -0.001
                      !   ARANGE(IVAR,2) = 1.001
                   END IF
                ELSE
                   VRANGE(1) = -HUGE(1.)
                   VRANGE(2) = HUGE(1.)
                END IF 
             ELSE
                HAS_VALID_RANGE(ISG,IVAR) = .TRUE.
             END IF
             VALID_RANGE(ISG,IVAR,:) = VRANGE
          END DO
          !
          IF ( OUTTYPE.EQ. 1 ) THEN
             !
             ! Read in static arrays on the input grid
             !
             ALLOCATE ( ARRIN(NXF(ISG),NYF(ISG)), stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARRIN alloc. error', IERR
             ARMAP = 0.
             !
             ISTART2(1) = 1
             ISTART2(2) = 1
             ICOUNT2(1) = NXF(ISG)
             ICOUNT2(2) = NYF(ISG)
             DO IVAR=1,NVAR
                ! Use fill value and valid range, preferably specific to 
                ! this subgrid
                ! otherwise from the first file containing values
                UNDEF(IVAR) = FILLVAL(ISG,IVAR)
                IF ( .NOT.HAS_FILLVAL(ISG,IVAR) ) THEN
                   DO ISG2 = 1,NSUBGRID
                      IF ( HAS_FILLVAL(ISG2,IVAR) ) THEN
                        UNDEF(IVAR) = FILLVAL(ISG2,IVAR)
                        EXIT
                      END IF
                   END DO
                END IF
                ARANGE(IVAR,:) = VALID_RANGE(ISG,IVAR,:)
                IF ( .NOT.HAS_VALID_RANGE(ISG,IVAR) ) THEN
                   DO ISG2 = 1,NSUBGRID
                      IF ( HAS_VALID_RANGE(ISG2,IVAR) ) THEN
                         ARANGE(IVAR,:) = VALID_RANGE(ISG2,IVAR,:)
                         EXIT
                      END IF
                   END DO
                END IF
                !
                ! Read in the variable array for this subgrid
                !
                IRET = nf90_get_var( NCID(ISG), NVID(ISG,IVAR), ARRIN, &
                                     start=ISTART2, count=ICOUNT2 )
                CALL CHECK_ERR(IRET,NDSE)
                WRITE (NDSO,4007) ISG, IVAR, TRIM(VARNAME(ISG,IVAR)),  &
                                  'read'
                !
                ! Interpolate static arrays to the single level subgrid
                !
                DO J=1,NYSG(ISG)
                   DO I=1,NXSG(ISG)
                      NSUM(IVAR) = 0
                      ASUM(IVAR) = 0.
                      WTSUM(IVAR) = 0.
                      DO K=1,4
                         II = Iint(I,J,K)
                         JJ = Jint(I,J,K)
                         WT = Wint(I,J,K)
                         IF ( II.LE.0 .OR. JJ.LE.0 ) CYCLE
                         ARRVAL = ARRIN(II,JJ)
                         IF ( ARRVAL.LT.ARANGE(IVAR,1) .OR.            &
                              ARRVAL.GT.ARANGE(IVAR,2) ) CYCLE
                         IF ( ABS(ARRVAL-FILLVAL(ISG,IVAR)).LT.TINY(1.) &
                             .AND. HAS_FILLVAL(ISG,IVAR) ) CYCLE
                         WTSUM(IVAR) = WTSUM(IVAR) + WT
                         NSUM(IVAR) = NSUM(IVAR) + 1
                         ASUM(IVAR) = ASUM(IVAR) +  WT*ARRVAL
                      END DO ! K loop
                      IF ( NSUM(IVAR).EQ.0 .OR.                        &
                           WTSUM(IVAR).LT.1.E-9 ) THEN
                         ARMAP(I,J,IVAR) = UNDEF(IVAR)
                      ELSE
                         ARMAP(I,J,IVAR) = ASUM(IVAR)/WTSUM(IVAR)
                      END IF
                   END DO  ! I loop 
                END DO  ! J loop
                WRITE (NDSO,4007) ISG, IVAR, TRIM(VARNAME(ISG,IVAR)),  &
                                  'interpolated'
             END DO  ! end of loop over variables
          ELSE ! test on OUTTYPE
             ! Transfer data to input array for quadtree mapping
             !  1 subgrid number, 2:5 interpolation X indices,
             !  6:9 interpolation X indices, 10:13 weights 
             ARMAP(:,:,1) = REAL(ISG)
             ARMAP(:,:,2:5) = REAL(Iint)
             ARMAP(:,:,6:9) = REAL(Jint)
             ARMAP(:,:,10:13) = Wint
             ! Valid ranges for each data type
             ARANGE(1:9,1) = 0.5
             ARANGE(1,2) = NSUBGRID + 0.5
             ARANGE(2:5,2) = NXSG(ISG) + 0.5
             ARANGE(6:9,2) = NYSG(ISG) + 0.5
             ARANGE(10:13,1) = 0.
             ARANGE(10:13,2) = 1.
!
             WRITE (NDSO,4008) ISG
          END IF  ! test on OUTTYPE
!
! 5.h Add cells for this subgrid to the quadtree grid structure
!
          IF ( ISG.EQ.1 ) THEN
             LVLDIF = LVLREL(ISG)-LV0
          ELSE
             LVLDIF = LVLREL(ISG)-LVLREL(IPSG(ISG))
          END IF
          ! SW corner of sub-grid origin (1,1) in reference grid 
          ! coordinates:
          X0VAL = (X0SG(ISG) - 0.5*SXSG(ISG) - X0)/SX + 1.
          Y0VAL = (Y0SG(ISG) - 0.5*SYSG(ISG) - Y0)/SY + 1.
          !write(*,*) 'calling QA_MLGTQ:, ISG= ', ISG
          !write(*,*) 'LVLDIF:',LVLDIF
          !write(*,*) 'LVLREL:',LVLREL(ISG)
          !write(*,*) 'LVLREF:',LVLREF
          !write(*,*) 'NX0,NY0:',NX0,NY0
          !write(*,*) 'X0VAL,Y0VAL:',X0VAL,Y0VAL
          CALL QA_MLG2QT( LVLDIF, LVLREL(ISG), LVLREF,                 &
                    NX0, NY0, X0VAL, Y0VAL, ARMAP, QTREE,              &
                    MAPML, ARRQ, atype=ATYPE, arange=ARANGE,           &
                    map_types=MAP_TYPES, ierr=IERR, ndse=NDSE )
          IF ( IERR .NE. 0 ) CALL EXTCDE ( IERR )
          QTREE%LVLHI = MAX ( QTREE%LVLHI, LVLREF + LVLREL(ISG) )
          !write(*,*) 'after QA_MLGTQ:, ISG= ', ISG
          !write(*,*) 'NQUAD:',QTREE%NQUAD
          !write(*,*) 'NCELL:',QTREE%NCELL
          !write(*,*) 'LVLREF:',QTREE%LVLREF
          !write(*,*) 'LVLMAX:',QTREE%LVLMAX
          !write(*,*) 'LVLHI:',QTREE%LVLHI
          NSEAMXB = QTREE%NCELL
          NQMXB = QTREE%NQUAD
          WRITE (NDSO,4009) ISG, LVLREL(ISG), LVLDIF
          !call qa_ioqt(NDSO,QTREE,3,IERR, NDSE)
          !call qa_ioqt(NDSO,QTREE,4,IERR, NDSE)
          !
          IF ( OUTTYPE.EQ. 1 ) THEN
             IF ( ALLOCATED(ARRIN) ) THEN
                DEALLOCATE ( ARRIN, stat=IERR )
                IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARRIN dealloc. error', &
                                                IERR
             END IF
          END IF
          IF ( ALLOCATED(ARMAP) ) THEN
             DEALLOCATE ( ARMAP, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARMAP dealloc. error',    &
                                             IERR
          END IF
          IF ( ALLOCATED(XSG) ) THEN
             DEALLOCATE ( XSG, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XSG dealloc. error', IERR
          END IF
          IF ( ALLOCATED(YSG) ) THEN
             DEALLOCATE ( YSG, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'YSG dealloc. error', IERR
          END IF
          IF ( ALLOCATED(Iint) ) THEN
             DEALLOCATE ( Iint, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'Iint dealloc. error', IERR
          END IF
          IF ( ALLOCATED(Jint) ) THEN
             DEALLOCATE ( Jint, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'Jint dealloc. error', IERR
          END IF
          IF ( ALLOCATED(Wint) ) THEN
             DEALLOCATE ( Wint, stat=IERR )
             IF ( IERR.NE.0 ) WRITE(NDSE,*) 'Wint dealloc. error', IERR
          END IF
!
! 5.i Remove any undefined cells from the quadtree grid structure, and
!     make the correponding changes to the data array:
!
          CALL QA_REMUNDEF( QTREE, newcell=NEWCELL, ierr=IERR,         &
                            ndse=NDSE )
          IF ( IERR .NE. 0 ) CALL EXTCDE ( IERR )
          IF ( NSEAMXB.GT.QTREE%NCELL ) THEN
             DO ICELL=1,NSEAMXB
                IST = NEWCELL(ICELL)
                IF ( IST.GT.0 ) ARRQ(IST,:) = ARRQ(ICELL,:)
             END DO
             NSEAMXB = QTREE%NCELL
             NQMXB = QTREE%NQUAD
             !call qa_ioqt(NDSO,QTREE,3,IERR, NDSE)
             !call qa_ioqt(NDSO,QTREE,4,IERR, NDSE)
          END IF
       END DO  ! loop over subgrids        INTEGER     ::  IERR

!
! 5.j  Compute the cell centre coordinates, and bounds
!
       ALLOCATE( X(NSEAMXB), Y(NSEAMXB), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'X, Y alloc. error', IERR
       ALLOCATE( XBOUND(NSEAMXB,2), YBOUND(NSEAMXB,2), stat=IERR )
       IF ( IERR.NE.0 ) WRITE(NDSE,*) 'XBOUND alloc. error', IERR
       DO ICELL=1,NSEAMXB
          X(ICELL) = X0 + SX*( QTREE%XYVAL(ICELL,1) - 1. )
          Y(ICELL) = Y0 + SY*( QTREE%XYVAL(ICELL,2) - 1. )
          SCALE = 2.**(LVLREF - QTREE%INDLVL(ICELL) -1)
          XBOUND(ICELL,1) = X(ICELL) - SX*SCALE
          XBOUND(ICELL,2) = X(ICELL) + SX*SCALE
          YBOUND(ICELL,1) = Y(ICELL) - SY*SCALE
          YBOUND(ICELL,2) = Y(ICELL) + SY*SCALE
       END DO
!
!  6.
!
       IF ( OUTTYPE.EQ.1 ) THEN
!
!  6.a For static data, write the outputs
!
          WRITE (NDSO,4010) 'static ASCII'
          !
          ! Quadtree cell file
          !
          OPEN( UNIT=NDSA, FILE=OUTFILEC, FORM='FORMATTED' )
          CALL QA_IOQT(NDSA,QTREE,4,IERR, NDSE)
          IF ( IERR .NE. 0 ) CALL EXTCDE ( IERR )
          CLOSE(NDSA)
          WRITE (NDSO,4011) 'quadtree cell structure',                 &
                            TRIM(OUTFILEC)
          !
          ! Quadtree quad file
          !
          OPEN( UNIT=NDSA, FILE=OUTFILEQ, FORM='FORMATTED' )
          CALL QA_IOQT(NDSA,QTREE,3,IERR, NDSE)
          IF ( IERR .NE. 0 ) CALL EXTCDE ( IERR )
          CLOSE(NDSA)
          WRITE (NDSO,4011) 'quadtree quad structure',                 &
                            TRIM(OUTFILEQ)
          !
          ! data file
          !
          OPEN( UNIT=NDSA, FILE=OUTFILE, FORM='FORMATTED' )
          FMT = '(I10,1X,XXE13.5)'
          WRITE(FMT(9:10),'(I2)') NVAR+2
          WRITE(NDSA,5100) NVAR+7, 'No. of lines in the header'
          WRITE(NDSA,5100) NVAR+3, 'No. of columns in the file'
          WRITE(NDSA,5200) 1, 'Cell number'
          IF ( FLAGLL ) THEN
             WRITE(NDSA,5200) 2, 'longitude (deg)'
             WRITE(NDSA,5200) 3, 'latitude (deg)'
          ELSE
             WRITE(NDSA,5200) 2, 'X (m)'
             WRITE(NDSA,5200) 3, 'Y (m)'
          END IF
          DO IVAR=1,NVAR
             WRITE(NDSA,5200) IVAR+3, TRIM(OUTVARNAME(IVAR))
          END DO
          WRITE(NDSA,*) '--------------------------------------'
          DO ICELL=1,NSEAMXB
             WRITE(NDSA,FMT) ICELL, X(ICELL), Y(ICELL),                &
                             (ARRQ(ICELL,IVAR),IVAR=1,NVAR)
          END DO
          CLOSE(NDSA)
          !
          WRITE (NDSO,4011) 'quadtree output fields',                  &
                            TRIM(OUTFILE)
       ELSE ! test on OUTTYPE
!
!  6.b Non-stationary data:
!
!  Re-map the interpolation indices and weights
          ALLOCATE ( IQint(NSEAMXB,9), WQint(NSEAMXB,4), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'IQint alloc. error', IERR
          IQint = NINT(ARRQ(1:NSEAMXB,1:9))
          WQint = ARRQ(1:NSEAMXB,10:13)
!
!  Open the netcdf output file, and define dimensions, variables,
!  IDs
!
          DIMNAME_OUT(1) = 'node'
          DIMNAME_OUT(2) = 'nbound'
          DIMNAME_OUT(3) = 'time'
          DIMLEN_OUT(1) = NSEAMXB
          DIMLEN_OUT(2) = 2
          DIMLEN_OUT(3) = 0
          IF ( FLAGLL ) THEN
            SVARNAME_OUT(1) = 'longitude'
            SVARNAME_OUT(2) = 'latitude'
            SVARNAME_OUT(3) = 'time'
            SVARNAME_OUT(4) = 'longitude_bounds'
            SVARNAME_OUT(5) = 'latitude_bounds'
          ELSE
            SVARNAME_OUT(1) = 'x'
            SVARNAME_OUT(2) = 'y'
            SVARNAME_OUT(3) = 'time'
            SVARNAME_OUT(4) = 'x_bounds'
            SVARNAME_OUT(5) = 'y_bounds'
          END IF
          !
          CALL OPNNCFILW( OUTFILE, DIMNAME_OUT, DIMLEN_OUT,             &
                          SVARNAME_OUT, OUTVARNAME, NCIDO, NSVIDO,      &
                          NVIDO, NDSE )
          !
          WRITE (NDSO,4010) 'netcdf'
          WRITE (NDSO,4011) 'dimensions and header', TRIM(OUTFILE)
          !
          !  Copy attributes from the first input file
          !
          !  global
          CALL CPATTS( NCID(1), NF90_GLOBAL, NCIDO, NF90_GLOBAL, NDSE )
          !  X/lon
          CALL CPATTS( NCID(1), NXYTVARID(1,1), NCIDO, NSVIDO(1), NDSE )
          !  Y/lat
          CALL CPATTS( NCID(1), NXYTVARID(1,2), NCIDO, NSVIDO(2), NDSE )
          !  time
          CALL CPATTS( NCID(1), NXYTVARID(1,3), NCIDO, NSVIDO(3), NDSE )
          !  output variables
          !
          DO IVAR=1,NVAR
             ! Use the fill value from the first input file
             ISG2 = 1
             IF ( .NOT.HAS_FILLVAL(1,IVAR) ) THEN
                ! ... or any other file if it doesn't have one
                DO ISG = 1,NSUBGRID
                   IF ( HAS_FILLVAL(ISG,IVAR) ) THEN
                      ISG2 = ISG
                      EXIT
                   END IF
                END DO
             END IF
             UNDEF(IVAR) = FILLVAL(ISG2,IVAR)
             CALL CPATTS( NCID(ISG2), NVID(ISG2,IVAR), NCIDO,          &
                          NVIDO(IVAR), NDSE )
          END DO
          !
          !  exit define mode
          IRET = nf90_enddef(NCIDO)
          CALL CHECK_ERR(IRET,NDSE)
          !
          ! Write the spatial arrays
          IRET = nf90_put_var( NCIDO, NSVIDO(1), X, start=(/1/),       &
                                                    count=(/NSEAMXB/) )
          CALL CHECK_ERR(IRET,NDSE)
          IRET = nf90_put_var( NCIDO, NSVIDO(2), Y, start=(/1/),       &
                                                    count=(/NSEAMXB/) )
          CALL CHECK_ERR(IRET,NDSE)
          ISTART2(1) = 1
          ISTART2(2) = 1
          ICOUNT2(1) = NSEAMXB
          ICOUNT2(2) = 2
          IRET = nf90_put_var( NCIDO, NSVIDO(4), XBOUND,               &
                               start=ISTART2, count=ICOUNT2 )
          CALL CHECK_ERR(IRET,NDSE)
          IRET = nf90_put_var( NCIDO, NSVIDO(5), YBOUND,               &
                               start=ISTART2, count=ICOUNT2 )
          CALL CHECK_ERR(IRET,NDSE)
          !
          WRITE (NDSO,4011) 'spatial data', TRIM(OUTFILE)
          !
          ALLOCATE ( ARRQ(NSEAMXB,NVAR), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'ARRQ alloc. error', IERR
          ALLOCATE ( SG(NSUBGRID), stat=IERR )
          IF ( IERR.NE.0 ) WRITE(NDSE,*) 'SG alloc. error', IERR
          DO ISG=1,NSUBGRID
            ALLOCATE ( SG(ISG)%ARRIN(NXF(ISG),NYF(ISG),NVAR), stat=IERR )
            IF ( IERR.NE.0 ) WRITE(NDSE,*) 'SG%ARRIN alloc. error', IERR
          END DO
          !
          ! time loop
          !
          ITMAX = 10000
          STOPNOW = .FALSE.
          DO ITOUT=1,ITMAX
            ! time index on each subgrid
            DO ISG=1,NSUBGRID
              ITSG(ISG) = ITFIRST(ISG) + (ITOUT-1)*ITSTEP(ISG)
              IF ( ITSG(ISG).GT.ITLAST(ISG) .OR.                       &
                   ITSG(ISG).GT.NTF(ISG) ) STOPNOW = .TRUE.
            END DO
            IF ( STOPNOW ) EXIT
            !
            ! Read the time from the first subgrid
            IRET = nf90_get_var( NCID(1), NXYTVARID(1,3), TIME,        &
                                 start=(/ITSG(1)/), count=(/1/) )
            CALL CHECK_ERR(IRET,NDSE)
            WRITE (NDSO,4012) 'time', 'read', ITSG(1), 1
            WRITE (NDSO,4013) TIME
            !
            ! Write the time
            IRET = nf90_put_var( NCIDO, NSVIDO(3), TIME,               &
                                 start=(/ITOUT/), count=(/1/) )
            CALL CHECK_ERR(IRET,NDSE)
            WRITE (NDSO,4011) 'time', TRIM(OUTFILE)
            !
            ! loop through subgrids
            !
            ARRQ = 0.
            !
            DO ISG=1,NSUBGRID
              ! Read in arrays on the input grid
              !
              ISTART3(1) = 1
              ISTART3(2) = 1
              ISTART3(3) = ITSG(ISG)
              ICOUNT3(1) = NXF(ISG)
              ICOUNT3(2) = NYF(ISG)
              ICOUNT3(3) = 1
              DO IVAR=1,NVAR
                IRET = nf90_get_var( NCID(ISG), NVID(ISG,IVAR),        &
                                     SG(ISG)%ARRIN(:,:,IVAR),          &
                                     start=ISTART3, count=ICOUNT3 )
                CALL CHECK_ERR(IRET,NDSE)
                WRITE (NDSO,4012) TRIM(VARNAME(ISG,IVAR)), 'read',     &
                                  ITSG(ISG), ISG
              END DO
            END DO ! end of loop over subgrids
            !
            ! Interpolate arrays to the quadtree
            !
            NREDO = 0
            DO ICELL=1,NSEAMXB
               ! Is this a retained cell with children?
               IF ( QTREE%INDSUB(ICELL).EQ.0 .AND. LVLMXB.GT.0 ) THEN
                  ! Yes: leave it to be averaged later from child values 
                  NREDO = NREDO + 1
                  CYCLE
               END IF
               ! Otherwise: interpolate from an input grid
               ! Which one do we use?
               ISG2 = IQint(ICELL,1)
               NSUM = 0
               ASUM = 0.
               WTSUM = 0.
               DO K=1,4
                  ! Indices for that input grid:
                  II = IQint(ICELL,K+1)
                  JJ = IQint(ICELL,K+5)
                  WT = WQint(ICELL,K)
                  DO IVAR=1,NVAR
                     IF ( II.LE.0 .OR. JJ.LE.0 ) CYCLE
                     ARRVAL = SG(ISG2)%ARRIN(II,JJ,IVAR)
                     IF ( HAS_VALID_RANGE(ISG2,IVAR) .AND.             &
                          (ARRVAL.LT.VALID_RANGE(ISG2,IVAR,1) .OR.     &
                           ARRVAL.GT.VALID_RANGE(ISG2,IVAR,2)) ) CYCLE
                     IF ( HAS_FILLVAL(ISG2,IVAR) .AND.                 &
                          ABS(ARRVAL-FILLVAL(ISG2,IVAR)).LT.TINY(1.) ) &
                                                               CYCLE
                     ASUM(IVAR) = ASUM(IVAR) + WT*ARRVAL
                     NSUM(IVAR) = NSUM(IVAR) + 1
                     WTSUM(IVAR) = WTSUM(IVAR) + WT
                  END DO  ! end of loop over variables
               END DO ! K loop
               DO IVAR=1,NVAR
                  IF ( NSUM(IVAR).EQ.0 .OR.                            &
                       WTSUM(IVAR).LT.TINY(1.) ) THEN
                     ARRQ(ICELL,IVAR) = UNDEF(IVAR)
                  ELSE
                     ARRQ(ICELL,IVAR) = ASUM(IVAR)/WTSUM(IVAR)
                  END IF
               END DO 
            END DO  ! end loop over cells
            !
            ! Further passes to average from finer to coarser levels
            IF ( NREDO.GT.0 ) THEN
               ! Loop from second finest to coarsest levels
               DO LVLT=LVLMXB-1,0,-1
                  DO ICELL=1, NSEAMXB
                     ! Is this a retained level LVLT cell with children?
                     IF ( QTREE%INDLVL(ICELL).EQ.LVLT .AND.            &
                          QTREE%INDSUB(ICELL).EQ.0 ) THEN
                        NSUM = 0
                        ASUM = 0.
                        ! loop over the child cells
                        IQ = QTREE%INDQUAD(ICELL)
                        DO K = 1,4
                           ICELL2 = QTREE%QICELL(IQ,K)
                           IF ( ICELL2.GT.0 ) THEN
                              IF ( QTREE%CELL_TYPE(ICELL2).EQ.         &
                                   QTREE%UNDEF_TYPE ) CYCLE
                              DO IVAR=1,NVAR
                                 ARRVAL = ARRQ(ICELL2,IVAR)
                                 ASUM(IVAR) = ASUM(IVAR) + ARRVAL
                                 NSUM(IVAR) = NSUM(IVAR) + 1
                              END DO  ! end of loop over variables
                           END IF
                        END DO ! loop over subquad index
                        DO IVAR=1,NVAR
                           IF ( NSUM(IVAR).GT.0 )                      &
                              ARRQ(ICELL,IVAR) = ASUM(IVAR)/NSUM(IVAR)
                        END DO  ! end of loop over variables
                        NREDO = NREDO - 1
                     END IF ! Test on level and having children to average
               
                  END DO ! loop over all cells
                  !  Quit if we've (re-)computed all cells with children
                  IF ( NREDO.LE.0 ) EXIT
               END DO  ! loop over refinement levels
            END IF ! Test on NREDO 
            !
            ! Write the interpolated variables to file
            !
            ISTART2(1) = 1
            ISTART2(2) = ITOUT
            ICOUNT2(1) = NSEAMXB
            ICOUNT2(2) = 1
            DO IVAR=1,NVAR
               IRET = nf90_put_var( NCIDO, NVIDO(IVAR), ARRQ(:,IVAR),  &
                                    start=ISTART2, count=ICOUNT2 )
               CALL CHECK_ERR(IRET,NDSE)
               WRITE (NDSO,4014) TRIM(OUTVARNAME(IVAR)),               &
                                 TRIM(OUTFILE),  ITOUT
            END DO
          END DO ! end of time loop
          !
          ! Close the output file
          IRET = nf90_close(NCIDO)
          CALL CHECK_ERR(IRET,NDSE)
       END IF ! test on OUTTYPE
!
!  7. Wrap up: close the input files
!
       !DO ISG=1,NSUBGRID
       !   IRET = nf90_close(NCID(ISG))
       !   CALL CHECK_ERR(IRET,NDSE)
       !END DO
       WRITE (NDSO,4015)
!
      GOTO 2222
!
! Escape locations read errors :
!
 2000 CONTINUE
      WRITE (NDSE,1000) IERR
      CALL EXTCDE ( 60 )
!
 2001 CONTINUE
      WRITE (NDSE,1001)
      CALL EXTCDE ( 61 )
!
 2002 CONTINUE
      WRITE (NDSE,1002) IERR
      CALL EXTCDE ( 62 )
!
 2222 CONTINUE
  900 FORMAT (/15X,'    *** WAVEWATCH III Quadtree preprocessor *** '/ &
               15X,'==============================================='/)
  901 FORMAT ( '  Comment character is ''',A,''''/)
 3000 FORMAT (/'  The reference grid: '/                              &
               ' --------------------------------------------------')
 3001 FORMAT ( '       Coordinate system           : ',A)
 3002 FORMAT ( '       Index closure type          : ',A)
 3003 FORMAT ( '       Dimensions                  : ',I6,I8)
 3004 FORMAT (/'       Increments           (deg.) :',2F10.4/         &
               '       Longitude range      (deg.) :',2F10.4/         &
               '       Latitude range       (deg.) :',2F10.4)
 3005 FORMAT ( '       Increments             (km) :',2F8.2/          &
               '       X range                (km) :',2F8.2/          &
               '       Y range                (km) :',2F8.2)
 3006 FORMAT ( '       Max. levels coarsened       : ',I2/            &
               '       Max. levels refined         : ',I2)
 3007 FORMAT ( ' --------------------------------------------------'/ &
               '          Output file: ',A)
 3008 FORMAT ( '  ASCII format, with additional files:'/              &
               '   Quadtree cell file: ',A/                           &
               '   Quadtree quad file: ',A/                           &
               '  Output variables:'/                                 &
               ' --------------------------------------------------' )
 3009 FORMAT ( '  Netcdf format'/                                     &
               '  Output variables:'/                                 &
               ' --------------------------------------------------' )
 3010 FORMAT ( '  var. ',I2,': ',A)
 3011 FORMAT ( ' --------------------------------------------------'/ &
               '  Number of subgrids/input files: ',I2/               &
               ' --------------------------------------------------' )
 3012 FORMAT ( '  grid: ',I2,' parent grid: ',I2,' file: ',A)
 3013 FORMAT ( '  grid: ',I2,' dim/var: ',I2,/                        &
               '    dimension name: ',A/                              &
               '     variable name: ',A)
 3014 FORMAT ( '  grid: ',I2,' time indices used: ',/                 &
               '    First: ',I3,' Last: ',I3,' Stride: ',I3)
 3015 FORMAT ( '  grid: ',I2,' field variable: ',I2,/                 &
               '     name in input file: ',A/                         &
               ' --------------------------------------------------' )
 4001 FORMAT ( '  grid: ',I2,' NXF: ',I6,' NYF: ',I6,' NTF: ',I6)
 4002 FORMAT ( '  grid: ',I2/                                         &
               '   Min. X/Lon cell size       : ',F10.6/              &
               '   Min. Y/Lat cell size       : ',F10.6/              &
               '   X/Longitude range          : ',2F10.4/             &
               '   Y/Latitude range           : ',2F10.4)
 4003 FORMAT ( '  grid: ',I2,' relative refinement level: ',I2)
 4004 FORMAT ( '  grid: ',I2,' ',A,':'/                               &
               '    grid size                 : ',2I6/                &
               '    X/Lon cell size           : ',F10.6/              &
               '    Y/Lat cell size           : ',F10.6/              &
               '    X/Lon origin              : ',F10.4/              &
               '    Y/Lat origin              : ',F10.4)
 4005 FORMAT ( ' --------------------------------------------------'/ &
               '  Size of multigrid: [local, global]: ',2I10/         &
               ' --------------------------------------------------' )
 4006 FORMAT ( '  Indices and weights computed for interpolation'/    &
               '  from input grid: ',I2,' of size ',I6,' x ',I6/      &
               '  to the nearest regular subgrid  ',I6,' x ',I6)
 4007 FORMAT ( '  grid: ',I2,' variable ',I2,' : ',A,1X,A)
 4008 FORMAT ( '  grid: ',I2,' interpolation indices and weights'/    &
               '          appended to multilevel array')
 4009 FORMAT ( '  grid: ',I2,' interpolation indices and weights'/    &
               '          transferred to quadtree structure'/         &
               '          at level ',I2,' and ',I2,' lower levels')
 4010 FORMAT ( '-------------------------------------------------'/   &
               '  Outputs in ',A,' format: ')
 4011 FORMAT ( '    ',A,' written to ',A)
 4012 FORMAT ( '    ',A,' ',A,' for time index ',I4/                  &
               '          from input file ',I2)
 4013 FORMAT ( '  time =  ',F10.4)
 4014 FORMAT ( '    ',A,' written to ',A,' for time index ',I4)
 4015 FORMAT ( ' -------------------------------------------------'/  &
               '                  All done '/                         &
               ' -------------------------------------------------' )
!
 5100 FORMAT (I10,1X,A)
 5200 FORMAT ('Column ',I2,1X,A)
!
 1000 FORMAT (/' *** ERROR IN W3PRNQ : '/                             &
               '     ERROR IN OPENING INPUT FILE'/                    &
               '     IOSTAT =',I5/)
!
 1001 FORMAT (/' *** ERROR IN W3PRNQ : '/                             &
               '     PREMATURE END OF INPUT FILE'/)
!
 1002 FORMAT (/' *** ERROR IN W3PRNQ : '/                             &
               '     ERROR IN READING FROM INPUT FILE'/               &
               '     IOSTAT =',I5/)
!
!/
!/ ------------------------------------------------------------------- /
!/
       CONTAINS
!/
!/ ------------------------------------------------------------------- /
!/
      subroutine gridint(M1,N1,X1,Y1,M2,N2,X2,Y2,Iint,Jint,Wint)
!
!     Richard Gorman, NIWA
!
!       Compute indices and weights for bilinear interpolation from one 
!       curvilinear  grid to another.
!       Assuming that:
!       Input grid has cells at coordinates 
!           X1(i,j), Y1(i,j), i=1,...,M1, j=1,...N1
!       Output grid has cells at coordinates 
!           X2(i,j), Y2(i,j), i=1,...,M2, j=1,...N2
!       Then indices Iint, Jint and weights Wint are computed such that,
!       if an array A1 is defined on grid 1, it can be interpolated to 
!       an array A2 defined on grid 2, with
!         A2(i,j) = sum(k=1:4) { Wint(i,j,k)*A1(Iint(i,j,k),Jint(i,j,k)) }
!
      implicit none
      INTEGER, INTENT(IN)    ::   M1, N1, M2, N2
      REAL, INTENT(IN)       ::   X1(M1,N1), Y1(M1,N1)
      REAL, INTENT(IN)       ::   X2(M2,N2), Y2(M2,N2)
      INTEGER, INTENT(OUT)   ::   Iint(M2,N2,4), Jint(M2,N2,4)
      REAL, INTENT(OUT)      ::   Wint(M2,N2,4)
!
!      Local variables
      INTEGER :: i1, i2, j1, j2, k, i1p, j1p
      REAL    :: XA, YA, X00, X0p, Xp0, Xpp, Y00, Y0p, Yp0, Ypp
      REAL    :: c000p, c0p00, cppp0, cp0pp
      REAL    :: cp000, c00p0, c0ppp, cpp0p
      REAL    :: cx000p, cx0p00, cxppp0, cxp0pp
      REAL    :: cxp000, cx00p0, cx0ppp, cxpp0p
      REAL    :: TOL, Wtsum
      INTEGER :: Itried(M1,N1),Jtried(M1,N1)
      INTEGER :: inew, jnew, i2step, i2start, i2end
      LOGICAL :: outside, nomove
!
      TOL = 1.e-20
!
      Itried = 0
      Jtried = 0
      i1 = 1
      j1 = 1
      DO j2=1,N2
         if (mod(j2,2).eq.1) then
            i2start = 1
            i2end = M2
            i2step = 1
         else
            i2start = M2
            i2end = 1
            i2step = -1
         end if
         !WRITE(*,*) 'j2, i2start, i2end, i2step =',j2, i2start, i2end, i2step
         DO i2=i2start,i2end,i2step
!             Coordinates of the Interpolation point A:
            XA = X2(i2,j2)
            YA = Y2(i2,j2)
            !WRITE(*,*) 'i,j, X, Y =',i2,j2, XA, YA
            DO k=1,4
               Iint(i2,j2,k) = 0
               Jint(i2,j2,k) = 0
               Wint(i2,j2,k) = 0.0
            END DO
            outside = .TRUE.
            DO WHILE (outside)
               !WRITE(*,'(2I4,3X,2I4)') i2,j2, i1,j1
               Itried(i1,j1) = i2
               Jtried(i1,j1) = j2
               i1p = i1 + 1
               j1p = j1 + 1
!
!             Coordinates of the corners of a quadrilateral in the input grid:
               X00 = X1(i1,j1)
               Xp0 = X1(i1p,j1)
               X0p = X1(i1,j1p)
               Xpp = X1(i1p,j1p)
               Y00 = Y1(i1,j1)
               Yp0 = Y1(i1p,j1)
               Y0p = Y1(i1,j1p)
               Ypp = Y1(i1p,j1p)
!
!             Compute cross products:
!                cSSTT  = cross product of vector from corner SS to interp. point A
!                               with vector from corner SS to corner TT
!                         SS, TT = 00, 0p, p0, or pp
!             Check if signs of these show that A is on the wrong (outer) side of
!             the sides of the quadrilateral
!             left & right walls: step if same side of both
               c000p = (XA-X00)*(Y0p-Y00) - (YA-Y00)*(X0p-X00)
               cp0pp = (XA-Xp0)*(Ypp-Yp0) - (YA-Yp0)*(Xpp-Xp0)
               c0p00 = (XA-X0p)*(Y00-Y0p) - (YA-Y0p)*(X00-X0p)
               cppp0 = (XA-Xpp)*(Yp0-Ypp) - (YA-Ypp)*(Xp0-Xpp)
 !
               nomove = .FALSE.
               if ( c000p*cp0pp .GT. 0. ) then
                  if ( ABS(c000p) .lt. ABS(cp0pp) ) then
                     inew = i1-1
                  ELSE
                     inew = i1+1
                  END IF
                  IF( inew.GE.1 .and. inew.LT.M1 ) then
                     IF(Itried(inew,j1).eq.i2 .and. Jtried(inew,j1).eq.j2) THEN
                        nomove = .TRUE.
                     ELSE
                        i1 = inew
                        !WRITE(*,*) 'wall test: change i'
                        CYCLE
                     END IF
                  ELSE
                     nomove = .TRUE.
                  END IF
               END IF
!             top & bottom walls: step if same side of both
               cp000 = (XA-Xp0)*(Y00-Yp0) - (YA-Yp0)*(X00-Xp0)
               c00p0 = (XA-X00)*(Yp0-Y00) - (YA-Y00)*(Xp0-X00)
               c0ppp = (XA-X0p)*(Ypp-Y0p) - (YA-Y0p)*(Xpp-X0p)
               cpp0p = (XA-Xpp)*(Y0p-Ypp) - (YA-Ypp)*(X0p-Xpp)
               if ( cp000*cpp0p .GT. 0. ) then
                  if ( ABS(cp000) .lt. ABS(cpp0p) ) then
                     jnew = j1-1
                  ELSE
                     jnew = j1+1
                  END IF
                  IF( jnew.GE.1 .and. jnew.LT.N1 ) THEN
                     IF(Itried(i1,jnew).eq.i2 .and. Jtried(i1,jnew).eq.j2) THEN
                        nomove = .TRUE.
                     ELSE
                        j1 = jnew
                        !WRITE(*,*) 'wall test: change j'
                        CYCLE
                     END IF
                  ELSE
                     nomove = .TRUE.
                  END IF
               END IF
!
!             If a move is needed but unavailable, we may be stuck in a loop.
!             Try starting from other points to find one not yet tested
               IF (nomove) THEN
                  ploop: do inew = 1,M1-1,M1-2
                     do jnew = 1,N1-1,N1-2
                        IF(Itried(inew,jnew).ne.i2 .and. Jtried(inew,jnew).ne.j2) THEN
                           i1 = inew
                           j1 = jnew
                           !WRITE(*,*) 'spiral: i2, j2, i1, j1 = ', i2,j2,i1,j1
                           nomove = .FALSE.
                           EXIT ploop
                        END IF
                     end do
                  end do ploop
                  IF (nomove) THEN
!             If a move is needed but unavailable, point A is outside grid 1: give up
                     EXIT
                  ELSE
                     CYCLE
                  END IF
               END IF
!
!             A is inside (or on a wall of) the quadrilateral.
               outside = .FALSE.
               !WRITE(*,*) i2,j2,i1,j1
!             Compute more cross products:
!                cxSSTT  = cross product of vector from corner SS to its diagonally
!			        opposite corner
!                               with vector from corner SS to corner TT
!                         SS, TT = 00, 0p, p0, or pp
               cx000p = (Xpp-X00)*(Y0p-Y00) - (Ypp-Y00)*(X0p-X00)
               cx0p00 = (Xp0-X0p)*(Y00-Y0p) - (Yp0-Y0p)*(X00-X0p)
               cxppp0 = (X00-Xpp)*(Yp0-Ypp) - (Y00-Ypp)*(Xp0-Xpp)
               cxp0pp = (X0p-Xp0)*(Ypp-Yp0) - (Y0p-Yp0)*(Xpp-Xp0)
               cxp000 = (X0p-Xp0)*(Y00-Yp0) - (Y0p-Yp0)*(X00-Xp0)
               cx00p0 = (Xpp-X00)*(Yp0-Y00) - (Ypp-Y00)*(Xp0-X00)
               cx0ppp = (Xp0-X0p)*(Ypp-Y0p) - (Yp0-Y0p)*(Xpp-X0p)
               cxpp0p = (X00-Xpp)*(Y0p-Ypp) - (Y00-Ypp)*(X0p-Xpp)
!             Save the indices:
               Iint(i2,j2,1) = i1
               Jint(i2,j2,1) = j1
               Iint(i2,j2,2) = i1p
               Jint(i2,j2,2) = j1
               Iint(i2,j2,3) = i1
               Jint(i2,j2,3) = j1p
               Iint(i2,j2,4) = i1p
               Jint(i2,j2,4) = j1p
!             Compute the weights:
               IF (ABS(cxppp0).GT.tol .AND. ABS(cxpp0p).GT.tol) THEN
                  Wint(i2,j2,1) = cppp0*cpp0p/(cxppp0*cxpp0p)
               END IF
               IF (ABS(cx0ppp).GT.tol .AND. ABS(cx0p00).GT.tol) THEN
                  Wint(i2,j2,2) = c0ppp*c0p00/(cx0ppp*cx0p00)
               END IF
               IF (ABS(cxp0pp).GT.tol .AND. ABS(cxp000).GT.tol) THEN
                  Wint(i2,j2,3) = cp0pp*cp000/(cxp0pp*cxp000)
               END IF
               IF (ABS(cx000p).GT.tol .AND. ABS(cx00p0).GT.tol) THEN
                  Wint(i2,j2,4) = c000p*c00p0/(cx000p*cx00p0)
               END IF
!             Ensure all weights in the range 0 to 1:
               Wint(i2,j2,1) = MIN(MAX(Wint(i2,j2,1),0.0),1.0)
               Wint(i2,j2,2) = MIN(MAX(Wint(i2,j2,2),0.0),1.0)
               Wint(i2,j2,3) = MIN(MAX(Wint(i2,j2,3),0.0),1.0)
               Wint(i2,j2,4) = MIN(MAX(Wint(i2,j2,4),0.0),1.0)
               Wtsum = Wint(i2,j2,1) + Wint(i2,j2,2)                  &
                     + Wint(i2,j2,3) + Wint(i2,j2,4)
               IF(Wtsum.GT.TOL) then
                   Wint(i2,j2,1) = Wint(i2,j2,1)/Wtsum
                   Wint(i2,j2,2) = Wint(i2,j2,2)/Wtsum
                   Wint(i2,j2,3) = Wint(i2,j2,3)/Wtsum
                   Wint(i2,j2,4) = Wint(i2,j2,4)/Wtsum
               END if

            END DO
         END DO
      END DO

      end subroutine
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE SGPROPS( NX, NY, XIN, YIN, XYPARS )
!/
!         Richard Gorman, NIWA
!           November, 2017:   Origination
!
!/
!  1. Purpose :
!
!     Evaluate some properties of the input grid
!
!  3. Parameters :
!       Name       Type I/O Description
!       ------------------------------------------------------
!       NX, NY     Int.  I  Size of arrays
!       XIN        R.A.  I  Rank 2 array of x coordinates
!       YIN        R.A.  I  Rank 2 array of y coordinates
!       XYPARS     R.A.  I  Array of grid properties:
!                             XYPARS(1:2) = minimum cell X,Y-extent
!                             XYPARS(3:4) = min. & max cell X value
!                             XYPARS(5:6) = min. & max cell Y value
!       ------------------------------------------------------
!
!--------------------------------------------------------------------
!
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!
!     Subroutine Parameters
!
      INTEGER, INTENT(IN)  :: NX, NY
      REAL, INTENT(IN)  :: XIN(NX,NY)
      REAL, INTENT(IN)  :: YIN(NX,NY)
      REAL, INTENT(OUT) :: XYPARS(6)
!
!--------------------------------------------------------------------
!
!     Local Parameters
!
      INTEGER   :: IX, IY, IXP, IYP
      REAL      :: XSW, XNW, XSE, XNE
      REAL      :: YSW, YNW, YSE, YNE
      REAL      :: XMIN, XMAX, YMIN, YMAX
      REAL      :: XLEN, YLEN
!
!--------------------------------------------------------------------
!
      DO IY=1,NY-1
        IYP = IY+1
        DO IX=1,NX-1
          IXP = IX+1
          XSW = XIN(IX,IY)
          YSW = YIN(IX,IY)
          XSE = XIN(IXP,IY)
          YSE = YIN(IXP,IY)
          XNW = XIN(IX,IYP)
          YNW = YIN(IX,IYP)
          XNE = XIN(IXP,IYP)
          YNE = YIN(IXP,IYP)
          !
          ! min & max X & Y in this cell
          XMIN = MIN(XSW, XNW)
          XMAX = MAX(XSE, XNE)
          YMIN = MIN(YSW, YSE)
          YMAX = MAX(YNW, YNE)
          !
          ! cell extent in X and Y directions, using midpoints
          XLEN = 0.5*(XSE + XNE - XSW - XNW) 
          YLEN = 0.5*(YNE + YNW - YSE - YSW)
          !
          IF ( IX.EQ.1 .AND. IY.EQ.1 ) THEN
             XYPARS(1) = XLEN
             XYPARS(2) = YLEN
             XYPARS(3) = XMIN
             XYPARS(4) = XMAX
             XYPARS(5) = YMIN
             XYPARS(6) = YMAX
          ELSE
          !
             XYPARS(1) = MIN(XYPARS(1),XLEN)
             XYPARS(2) = MIN(XYPARS(2),YLEN)
             XYPARS(3) = MIN(XYPARS(3),XMIN)
             XYPARS(4) = MAX(XYPARS(4),XMAX)
             XYPARS(5) = MIN(XYPARS(5),YMIN)
             XYPARS(6) = MAX(XYPARS(6),YMAX)
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE SGPROPS
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE OPNNCFILR( FNAME, NXYTDIM, NXYTVAR, NVARIO,           &
                            XYTDIMNAME, XYTVARNAME, VARNAME,           &
                            NCID, XYTDIMLEN, XYTDIMID, XYTVARID,       &
                            XYTVRANK, VARID, NDSE )
!/
!         Richard Gorman, NIWA
!           November, 2017:   Origination
!
!/
!  1. Purpose :
!
!     Opens a netcdf file for reading, and identifies lengths and IDs
!     of specified dimensions, and IDs of specified variables
!
!
!  3. Parameters :
!       Name       Type I/O Description
!       ------------------------------------------------------
!       FNAME      Char  I  Name of the netcdf file
!       NXYTDIM    Int.  I  Number of target dimensions
!       NXYTVAR    Int.  I  Number of target spatio=temporal variables
!       NVARIO     Int.  I  Number of field variables
!       XYTDIMNAME C.A.  I  Array of target dimension names
!       XYTVARNAME C.A.  I  Array of target spatio-temporal variable names
!       VARNAME    C.A.  I  Array of target field variable names
!       NCID       Int   O  Netcdf file ID
!       XYTDIMLEN  I.A.  O  Array of dimension lengths
!       XYTDIMID   I.A.  O  Array of dimension IDs
!       XYTVARID   I.A.  O  Array of spatio-temporal variable IDs
!       XYTVRANK   I.A.  O  Array of spatio-temporal variable ranks
!       VARID      I.A.  O  Array of field variable IDs
!       NDSE       Int.  I  Unit number for error messages
!       ------------------------------------------------------
!
! ------------------------------------------------------------------- 
      USE NETCDF
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!
!     Subroutine Parameters
!
      CHARACTER*(*), INTENT(IN) :: FNAME
      INTEGER, INTENT(IN) :: NXYTDIM
      INTEGER, INTENT(IN) :: NXYTVAR
      INTEGER, INTENT(IN) :: NVARIO
      CHARACTER*(*), INTENT(IN) :: XYTDIMNAME(NXYTDIM)
      CHARACTER*(*), INTENT(IN) :: XYTVARNAME(NXYTVAR)
      CHARACTER*(*), INTENT(IN) :: VARNAME(NVARIO)
      INTEGER, INTENT(OUT) :: NCID
      INTEGER, INTENT(OUT) :: XYTDIMLEN(NXYTDIM)
      INTEGER, INTENT(OUT) :: XYTDIMID(NXYTDIM)
      INTEGER, INTENT(OUT) :: XYTVARID(NXYTVAR)
      INTEGER, INTENT(OUT) :: XYTVRANK(NXYTVAR)
      INTEGER, INTENT(OUT) :: VARID(NVARIO)
      INTEGER, INTENT(IN)  :: NDSE
!
!--------------------------------------------------------------------
!
!     Local Parameters
!
      INTEGER  ::  IRET
      INTEGER  ::  NDIMS, NVARS, NGATTS, UNLIMDIMID
      INTEGER  ::  ID, IO, DLEN
      CHARACTER*100 :: NAME
!
!--------------------------------------------------------------------
!
! Default outputs
!
      XYTDIMLEN = -1
      XYTDIMID = -1
      XYTVARID = -1
      VARID = -1
!
      !write(*,*) 'Opening file: ', trim(fname)
      IRET = nf90_open(FNAME, nf90_nowrite, NCID)
      CALL CHECK_ERR(IRET,NDSE)
      IRET = nf90_inquire(NCID, NDIMS, NVARS, NGATTS, UNLIMDIMID )
      CALL CHECK_ERR(IRET,NDSE)
      !write(*,*) 'NCID: ', NCID
      !write(*,*) 'NDIMS: ', NDIMS
      !write(*,*) 'NVARS: ', NVARS
      !write(*,*) 'NGATTS: ', NGATTS
!
      DO ID=1,NDIMS
         IRET = nf90_inquire_dimension(ncid, ID, NAME, DLEN )
         CALL CHECK_ERR(IRET,NDSE)
         !write(*,*) 'dimension  ', ID, trim(name), DLEN
         DO IO=1,NXYTDIM
            IF ( TRIM(NAME) .EQ. TRIM(XYTDIMNAME(IO)) ) THEN
               XYTDIMID(IO) = ID
               XYTDIMLEN(IO) = DLEN
               !write(*,*) 'match with target dim. ', IO
               EXIT
            END IF
         END DO
      END DO
!
      DO ID=1,NVARS
         IRET = nf90_inquire_variable(ncid, ID, NAME, ndims=NDIMS )
         CALL CHECK_ERR(IRET,NDSE)
         !write(*,*) 'variable  ', ID, trim(name)
         !write(*,*) 'NDIMS = ', NDIMS
         DO IO=1,NXYTVAR
            IF ( TRIM(NAME) .EQ. TRIM(XYTVARNAME(IO)) ) THEN
               XYTVARID(IO) = ID
               XYTVRANK(IO) = NDIMS
               !write(*,*) 'match with target XYT var. ', IO
               EXIT
            END IF
         END DO
         DO IO=1,NVARIO
            IF ( TRIM(NAME) .EQ. TRIM(VARNAME(IO)) ) THEN
               VARID(IO) = ID
               !write(*,*) 'match with target field var. ', IO
               EXIT
            END IF
         END DO
      END DO
      RETURN
      END SUBROUTINE OPNNCFILR
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE OPNNCFILW( FNAME, DIMNAME, DIMLEN, SVARNAME,       &
                            VARNAME, NCID, SVARID, VARID, NDSE )
!/
!         Richard Gorman, NIWA
!           November, 2017:   Origination
!
!/
!  1. Purpose :
!
!     Opens a netcdf file for writing, creates dimensions of specified 
!     names and lengths, and variables of specified names
!
!
!  3. Parameters :
!       Name       Type I/O Description
!       ------------------------------------------------------
!       FNAME      Char  I  Name of the netcdf file
!       DIMNAME    C.A.  I  Array of dimension names
!       DIMLEN     I.A.  I  Array of dimension lengths
!       SVARNAME   C.A.  I  Array of spatio-temporal variable names
!       VARNAME    C.A.  I  Array of field variable names
!       NCID       Int   O  Netcdf file ID
!       SVARID     I.A.  O  Array of spatio-temporal variable IDs
!       VARID      I.A.  O  Array of field variable IDs
!       NDSE       Int.  I  Unit number for error messages
!       ------------------------------------------------------
!
! ------------------------------------------------------------------- 
      USE NETCDF
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!
!     Subroutine Parameters
!
      CHARACTER*(*), INTENT(IN) :: FNAME
      CHARACTER*(*), INTENT(IN) :: DIMNAME(:)
      INTEGER, INTENT(IN)       :: DIMLEN(:)
      CHARACTER*(*), INTENT(IN) :: SVARNAME(:)
      CHARACTER*(*), INTENT(IN) :: VARNAME(:)
      INTEGER, INTENT(OUT)      :: NCID
      INTEGER, INTENT(OUT)      :: SVARID(:)
      INTEGER, INTENT(OUT)      :: VARID(:)
      INTEGER, INTENT(IN)       :: NDSE
!
!--------------------------------------------------------------------
!
!     Local Parameters
!
      INTEGER  ::  IRET, NVAR, NDIM, NSVAR
      INTEGER  ::  SDIMID, BDIMID, TDIMID
      INTEGER  ::  DLEN
      INTEGER  ::  ID
      CHARACTER*100 :: NAME
!
!--------------------------------------------------------------------
!
      NVAR = SIZE(VARNAME,1)
      NVAR = MIN(NVAR,SIZE(VARID,1))
      NDIM = SIZE(DIMNAME,1)
      NDIM = MIN(NDIM,SIZE(DIMLEN,1))
      NSVAR = SIZE(SVARNAME,1)
      NSVAR = MIN(NSVAR,SIZE(SVARID,1))
!
! Default outputs
!
      SVARID = -1
      VARID = -1
!
! Default DIM IDS
!
      SDIMID = -1
      BDIMID = -1
      TDIMID = -1
!
      IRET = nf90_open(FNAME, nf90_write, NCID)
      CALL CHECK_ERR(IRET,NDSE)
!
      DO ID=1,NDIM
         NAME = TRIM(DIMNAME(ID))
         DLEN = DIMLEN(ID)
         IF ( DLEN.LE.0 ) THEN 
            DLEN = NF90_UNLIMITED
            IRET = nf90_def_dim(ncid, NAME, DLEN, TDIMID )
         ELSEIF ( DLEN.LE.2 ) THEN 
            IRET = nf90_def_dim(ncid, NAME, DLEN, BDIMID )
         ELSE
            IRET = nf90_def_dim(ncid, NAME, DLEN, SDIMID )
         END IF
         CALL CHECK_ERR(IRET,NDSE)
      END DO
!
!  Define spatio-temporal variables
!
      DO ID=1,NSVAR
         IF ( ID.EQ.1 .OR. ID.EQ.2 ) THEN
           ! X/lon and Y/lat have the node dimension
           IRET = nf90_def_var( ncid, NAME, NF90_FLOAT, SDIMID,     &
                                SVARID(ID) )
         ELSE IF ( ID.EQ.3 ) THEN
           ! time has the time dimension
           IRET = nf90_def_var( ncid, NAME, NF90_FLOAT, TDIMID,     &
                                SVARID(ID) )
         ELSE IF ( ID.EQ.4 .OR. ID.EQ.5 ) THEN
           ! X/lon and Y/lat bounds have the node and nbounds dimensions
           IRET = nf90_def_var( ncid, NAME, NF90_FLOAT,             &
                                (/SDIMID,BDIMID/), SVARID(ID) )
         END IF
         CALL CHECK_ERR(IRET,NDSE)
      END DO
!
!  Define field variables: all have node and time dimensions
!
      DO ID=1,NVAR
         IRET = nf90_def_var( ncid, NAME, NF90_FLOAT,               &
                              (/SDIMID,TDIMID/), VARID(ID) )
         CALL CHECK_ERR(IRET,NDSE)
      END DO
      RETURN
      END SUBROUTINE OPNNCFILW
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE CPATTS( NCID1, VARID1, NCID2, VARID2, NDSE )
!/
!         Richard Gorman, NIWA
!           November, 2017:   Origination
!
!/
!  1. Purpose :
!
!     Copy attributes of variable(s) in netcdf file to corrsponding
!     variables in another file.
!
!  3. Parameters :
!       Name       Type I/O Description
!       ------------------------------------------------------
!       NCID1      Int   I  Source netcdf file ID
!       VARID1     I.A.  I  Array of source variable IDs
!       NCID2      Int   I  Destination netcdf file ID
!       VARID2     I.A.  I  Array of destination variable IDs
!       NDSE       Int.  I  Unit number for error messages
!       ------------------------------------------------------
!
! ------------------------------------------------------------------- 
      USE NETCDF
      IMPLICIT NONE
!
!--------------------------------------------------------------------
!
!     Subroutine Parameters
!
      INTEGER, INTENT(IN)      :: NCID1
      INTEGER, INTENT(IN)      :: VARID1(:)
      INTEGER, INTENT(IN)      :: NCID2
      INTEGER, INTENT(IN)      :: VARID2(:)
      INTEGER, INTENT(IN)      :: NDSE
!
!--------------------------------------------------------------------
!
!     Local Parameters
!
      INTEGER  ::  IRET, NVAR, NATTS
      INTEGER  ::  VID1, VID2
      INTEGER  ::  ID, IATT
      CHARACTER*100 :: ATTNAME
!
!--------------------------------------------------------------------
!
      NVAR = SIZE(VARID1,1)
      NVAR = MIN(NVAR,SIZE(VARID2,1))
!
!
      DO ID=1,NVAR
         VID1 = VARID1(ID)
         VID2 = VARID2(ID)
         IRET = nf90_inquire_variable(NCID1, VID1, nAtts=NATTS )
         CALL CHECK_ERR(IRET,NDSE)
         DO IATT=1,NATTS
            IRET = nf90_inq_ATTNAME(NCID1, VID1, IATT, ATTNAME )
            CALL CHECK_ERR(IRET,NDSE)
            IRET = nf90_copy_att( NCID1, VID1, ATTNAME, NCID2, VID2 )
            CALL CHECK_ERR(IRET,NDSE)
         END DO
      END DO
      RETURN
      END SUBROUTINE CPATTS
!/
!/ ------------------------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE CHECK_ERR(IRET, NDSE)
      USE NETCDF
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IRET
      INTEGER, INTENT(IN) :: NDSE
      IF ( IRET /= nf90_noerr ) THEN
         WRITE(NDSE,*) nf90_strerror(IRET)
         STOP
      END IF
      END SUBROUTINE CHECK_ERR
!/
!/ ------------------------------------------------------------------- /
!/

       END PROGRAM W3PRNQ
!/
!/ ------------------------------------------------------------------- /
!/
