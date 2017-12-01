!/
!/ ------------------------------------------------------------------- /
!/
      MODULE QA_UTILS
!/
!/
!  1. Purpose :
!
!     Routines for handling quadtree grids.
!
!  2. Method :
!
!
!      Name         Type  Scope    Description
!     ---------------------------------------------------------------
!     QA_ADGR       Subr. Public   Adapt the quadtree grid, based on a 
!                                  diagnostic variable
!     QA_ADVAR      Subr. Public   Recompute an array after adapting a 
!                                  quadtree grid.
!     QA_ALLOC      Subr. Public   Allocate a threaded quadtree structure
!     QA_BMLSTRUC   Subr. Public   Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC    Subr. Public   Create a threaded quadtree structure
!                                  from a rectangular grid
!     QA_COARSEN    Subr. Public   Coarsen selected quads of a
!                                  quadtree grid
!     QA_CPQT       Subr. Public   Copy contents of one quadtree structure 
!                                  to another
!     QA_D2WTS      Subr. Public   Compute weights for higher derivatives
!                                  using second-order neighbours
!     QA_DERIVC     Subr. Public   Evaluate 1st to 3rd order
!                                  derivatives
!     QA_DERIVWTS   Subr. Public   Compute weights for 1st to 3rd
!                                  order derivatives
!     QA_DERWTS     Subr. Public   Creates arrays of precomputed 
!                                  interpolation weights for all
!                                  possible arrangements of quadtree 
!                                  neighbours
!     QA_DISPREL    Subr. Public   Compute CG using the
!                                  dispersion relation
!     QA_FILLSPARE  Subr. Public   Renumber cells to eliminate spare 
!                                  cells 
!     QA_FINDNBR    Subr. Public   Compute the seapoint neighbour
!                                  indices
!     QA_INCR       Subr. Public   Increment a multiple-digit counter 
!                                  in arbitrary base
!     QA_INTERP     Subr. Public   Interpolate to a target location 
!                                  using weighted values from 
!                                  neighbouring cells.
!     QA_IOQT       Subr. Public   Unformatted IO of a quadtree structure
!     QA_LINGAUSS   Subr. Public   Solve linear equations by
!                                  Gaussian elimination
!     QA_MLG2QT     Subr. Public   Derive a quadtree structure and compute 
!                                  values for all cells in it for quantities
!                                  (e.g. depth, transmission arrays, winds,
!                                  currents, ice, ...) provided on a set of
!                                  regular input grids 
!     QA_MLGRID     Subr. Public   Precompute depth and transmission arrays, 
!                                  for all possible cells of a multilevel grid 
!                                  from inputs given on a set of regular input
!                                  grids
!     QA_NBDIST     Subr. Public   Identifies the relative distribution
!                                  of neighbours around each cell
!     QA_PRCELL     Subr. Public   Print out the quadtree cell
!                                  structure
!     QA_Q2NOKEEP   Subr. Public   Convert a quadtree from a form that retains
!                                  refined cells to a form that doesn't  
!     QA_Q2QMAP     Subr. Public   Mapping between cells in two quadtree 
!                                  structures
!     QA_QMEAN      Subr. Public   Compute the mean of a variable
!                                  over the cells of each quad
!     QA_QT2SMC     Subr. Public   Derive SMC cell and face arrays from a
!                                  quadtree strucrure
!     QA_QTCHECK    Subr. Public   Perform consistency checks on a
!                                  quadtree strucrure
!     QA_READWT     Subr. Public   Read in precomputed quadtree
!                                  weights
!     QA_RECT       Subr. Public   Set up a quadtree 
!                                  representation (level 0)
!     QA_REFINE     Subr. Public   Refine selected cells of a 
!                                  quadtree grid
!     QA_REFVAL     Subr. Public   Assign values of a field in refined 
!                                  cells using estimated gradients 
!     QA_REMUNDEF   Subr. Public   Remove undefined cells
!     QA_TSORDER    Subr. Public   Compute the order to carry out
!                                  time iterations for multiple
!                                  levels
!     QA_XY2CELL    Subr  Public   Given the (reference grid) coordinates
!                                  of a point, locate the cell (if any) in a 
!                                  threaded quadtree structure containing it.
!
!     ML_CHILD		Func. Public   Return the multilevel index
!                                  for a cell's first (SW) child
!     ML_INDEX      Func. Public   Return an index into multilevel arrays for 
!                                  a cell in a quadtree structure
!     ML_PARENT		Func. Public   Return the multilevel index
!                                  for a cell's parent
!     MULTILEVEL_INDEX  Func. Public   Return an index into
!                                      multilevel arrays for a
!                                      QA_TREE cell
!     Qsort         Subr. Public   Quicksort
!
!     ---------------------------------------------------------------
!
      PUBLIC
!
      TYPE QA_TREE
         INTEGER :: NQUAD
         INTEGER :: NCELL
         INTEGER :: NCELL_DEF
         INTEGER :: LVLREF
         INTEGER :: LVLMAX
         INTEGER :: LVLHI
         INTEGER :: NX0
         INTEGER :: NY0
         INTEGER :: UNDEF_TYPE
         LOGICAL :: KEEP_REF
         LOGICAL :: DYNAMIC
         INTEGER :: IWTORDER
         INTEGER, ALLOCATABLE :: INDQUAD(:)
         INTEGER, ALLOCATABLE :: INDSUB(:)
         INTEGER, ALLOCATABLE :: INDLVL(:)
         INTEGER, ALLOCATABLE :: NGBR(:,:)
         INTEGER, ALLOCATABLE :: CELL_TYPE(:)
         INTEGER, ALLOCATABLE :: INDML(:)
         INTEGER, ALLOCATABLE :: NCASE(:)
         REAL, ALLOCATABLE    :: XYVAL(:,:)
         INTEGER, ALLOCATABLE :: QICELL(:,:)
         INTEGER, ALLOCATABLE :: QLEVEL(:)
         INTEGER, ALLOCATABLE :: QPARENT(:)
         INTEGER, ALLOCATABLE :: QNBR(:,:)
         INTEGER, ALLOCATABLE :: QCHILD(:,:)
         INTEGER, ALLOCATABLE :: INDWT(:)
         INTEGER, ALLOCATABLE :: GNBR(:,:)
         INTEGER, ALLOCATABLE :: NGNBR(:)
      END TYPE QA_TREE
!
!      QA_TREE structure components:
!          NCELL   Int.  Number of cells (including undefined cells)
!          NCELL_DEF Int.  Number of active cells
!          NQUAD   Int.  Number of quads
!          LVLREF  Int.  Refinement level of reference grid
!          LVLMAX  Int.  Maximum refinement level allowed
!          LVLHI   Int.  Maximum refinement level present
!          NX0     Int.  X dimension of level-zero base grid
!          NY0     Int.  Y dimension of level-zero base grid
!          UNDEF_TYPE Int. value of flag for unused cell indices
!          KEEP_REF Log. retain cell indices for refined cells
!          DYNAMIC Log.  Dynamic, i.e. nonstationary, grid 
!          IWTORDER Int. Order of interpolation weights included
!          INDQUAD I.A.  Index of the quad containing each input cell
!          INDSUB  I.A.  Index of the sub-quad containing each cell:
!                               0 = none (home cell is at level zero)
!                               1,2,3,4 = SW, SE, NW, NE corner of quad
!          INDLVL  I.A.  Level of each cell
!          NGBR    I.A.  Sea-point indices for the neighbours of input cells
!                           NGBR(I,1:4) = index of primary (i.e. of equal or 
!                           lower level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!                           NGBR(I,5:8) = index of secondary (i.e. of higher 
!                           level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!          CELL_TYPE I.A. Flag for wet/dry/boundary/etc. cell
!          INDML   I.A.  Multilevel index of cells
!          NCASE   I.A.  Index identifying the relative configuration of 
!                        neighbouring cells (out of 2401 possibilities)
!                        (used if IWTORDER>=1)
!          XYVAL   R.A.  Coordinates of cells (w.r.t. 
!                          reference rectangular grid)
!          QICELL  I.A.  Sea-point indices for the cells of the quad
!          QLEVEL  I.A.  Level of the cells within each quad
!          QPARENT I.A.  Quad index for the parent of each quad
!          QNBR    I.A.  Quad indices for the neighbours of each quad
!          QCHILD  I.A.  Quad indices for the children of each quad
!        if IWTORDER>=2:
!          INDWT   I.A.  Index, for each cell, into the lookup table of
!                        (second order) weights
!          GNBR    I.A.  Grand neighbours of each cell
!          NGNBR   I.A.  Number of grand neighbours of each cell
!
      TYPE QA_SPARE
         INTEGER :: NPART
         INTEGER, ALLOCATABLE :: PARTNO(:)
         INTEGER, ALLOCATABLE :: JSVAL(:)
         INTEGER, ALLOCATABLE :: NCLOC(:)
         INTEGER, ALLOCATABLE :: ISVAL(:,:)
         INTEGER, ALLOCATABLE :: NSPARE(:)
         INTEGER, ALLOCATABLE :: JSSPARE(:,:)
      END TYPE QA_SPARE
!
!      QA_SPARE structure components:
!          NPART   Int.  Number of partitions
!          PARTNO  I.A.  Partition number for each cell
!          JSVAL   I.A.  Local cell index for each global index
!          NCLOC   I.A.  Actual number of cells in each partition
!          ISVAL   I.A.  Global cell index for each local index & partition
!          NSPARE  I.A.  Number of spare local cell indices
!          JSSPARE I.A.  List of spare local cell indices
!
      TYPE QA_WEIGHT
         INTEGER :: MCASE=2401
         INTEGER :: NNGBR=8
         INTEGER :: NCASEQUAD
         LOGICAL :: WTSREAD(4)
         INTEGER :: IQUADDIR(2401,4)
         REAL    :: XGRADWTC(2401,0:8)
         REAL    :: YGRADWTC(2401,0:8) 
         REAL    :: XGRADWTW(2401,0:8)
         REAL    :: XGRADWTE(2401,0:8)
         REAL    :: YGRADWTS(2401,0:8)
         REAL    :: YGRADWTN(2401,0:8) 
         REAL    :: LAMBDA0(2401,0:8)
         REAL    :: VALWTW(2401,0:8)
         REAL    :: VALWTE(2401,0:8)
         REAL    :: VALWTS(2401,0:8)
         REAL    :: VALWTN(2401,0:8) 
         REAL    :: XCOORDQUAD(2401,0:8)
         REAL    :: YCOORDQUAD(2401,0:8) 
         REAL    :: XXDER2WTC(2401,0:8) 
         REAL    :: XYDER2WTC(2401,0:8) 
         REAL    :: YYDER2WTC(2401,0:8) 
      END TYPE QA_WEIGHT
!
!      QA_WEIGHT structure components:
!          MCASE      Int. Allocated number of possible cases (2401)
!          NNGBR      Int. Maximum number of neighbouring cells (8)
!          WTSREAD    L.A. Flags for whether arrays have been read fully:
!                           WTSREAD(1) = .TRUE. if centre gradient weights read
!                           WTSREAD(2) = .TRUE. if wall value weights read
!                           WTSREAD(3) = .TRUE. if wall gradient weights read
!                           WTSREAD(4) = .TRUE. if centre 2nd der. weights read
!          NCASEQUAD  Int. Number of precomputed cell configurations
!          IQUADDIR   I.A. Code for type of neighbour(s) in
!                              [W,E,S,N] direction 
!                           = 1 for a neighbour of equal level
!                           = 2 for two neighbours of higher level
!                           = 3 for a neighbour of lower level,
!                               displaced in -ve direction
!                               (S for E,W nbrs, W for S,N nbrs)
!                           = 4 for a neighbour of lower level,
!                               displaced in +ve direction
!                               (N for E,W nbrs, E for S,N nbrs)
!                           = 5 for one neighbour of higher level,
!                               displaced in -ve direction
!                               (S for E,W nbrs, W for S,N nbrs)
!                           = 6 for one neighbour of higher level,
!                               displaced in +ve direction
!                               (N for E,W nbrs, E for S,N nbrs)
!                           = 7 for no neighbour
!  Weights at each cell of the precomputed configuration:
!          XGRADWTC   R.A. Weights for X-gradient at cell centre
!          YGRADWTC   R.A. Weights for Y-gradient at cell centre
!          XGRADWTW   R.A. Weights for X-gradient at west wall
!          XGRADWTE   R.A. Weights for X-gradient at east wall
!          YGRADWTS   R.A. Weights for Y-gradient at south wall
!          YGRADWTN   R.A. Weights for Y-gradient at north wall
!          LAMBDA0    R.A. Natural neighbour coordinates relative 
!                            to cell centre
!          VALWTW     R.A. Weights for value at west wall
!          VALWTE     R.A. Weights for value at east wall
!          VALWTS     R.A. Weights for value at south wall
!          VALWTN     R.A. Weights for value at north wall
!          XCOORDQUAD R.A. X-coordinates of target and neighbour 
!                           cell centres
!          YCOORDQUAD R.A. Y-coordinates of target and neighbour
!                           cell centres
!          XXDER2WTC  R.A   Weights for d2/dXdX at cell centre
!          XYDER2WTC  R.A   Weights for d2/dXdY at cell centre
!          YYDER2WTC  R.A   Weights for d2/dYdY at cell centre
!
      TYPE QA_WEIGHTS1
         INTEGER :: MCASE=2401
         INTEGER :: NNGBR=8
         INTEGER :: IQUADDIR(2401,4), MAXORDER(2401)
         REAL    :: DERWTC(2401,0:8,5)
         REAL    :: XYREL(2401,0:8,2)
         REAL    :: DETERM(2401)
         REAL    :: LAMBDA0(2401,0:8)
      END TYPE QA_WEIGHTS1
!
!      QA_WEIGHTS1 structure components:
!          MCASE    Int.    Maximum number of configurations of neighbours
!                           around each cell in a quadtree structure
!          NNGBR    Int.    Maximum number of neighbouring cells 
!                           (8 for quadtree)
!          IQUADDIR I.A     Code for type of neighbour(s) in
!                           each direction, for each possible case IC
!                           IQUADDIR(IC,P)
!                            = 1 for a neighbour of equal level
!                   		 = 2 for two neighbours of higher level
!                   		 = 3 for a neighbour of lower level,
!                                displaced in -ve direction
!                   		     (S for E,W nbrs, W for S,N nbrs)
!                   		 = 4 for a neighbour of lower level,
!                                displaced in +ve direction
!                   		     (N for E,W nbrs, E for S,N nbrs)
!                   		 = 5 for one neighbour of higher level,
!                                displaced in -ve direction
!                   		     (S for E,W nbrs, W for S,N nbrs)
!                   		 = 6 for one neighbour of higher level,
!                                displaced in +ve direction
!                   		     (N for E,W nbrs, E for S,N nbrs)
!                   		 = 7 for no neighbour
!                          in the [W,E,S,N] direction for P = [1,2,3,4]
! 	   MAXORDER I.A    Maximum order of derivative computed for each case.
!          DERWTC   R.A    Weights for derivatives at cell centre:
!                          DERWTC(IC,K,1) = wt for ngbr K in d/dX for case IC,
!                          DERWTC(IC,K,2) = wt for ngbr K in d/dY for case IC,
!                          DERWTC(IC,K,3) = wt for ngbr K in d2/dX2 for case IC,
!                          DERWTC(IC,K,4) = wt for ngbr K in d2/dY2 for case IC,
!                          DERWTC(IC,K,5) = wt for ngbr K in d2/dXdY for case IC
!          XYREL    R.A    X- and Y- coordinates of target and neighbour 
!                              cell centres (relative)
!          DETERM   R.A    Determinant of matrix used to derive weights,
!                          for each case.
!          LAMBDA0  R.A.   Natural neighbour coordinates relative 
!                            to cell centre
!
      TYPE QA_WEIGHTS2
         INTEGER :: NCSAV
         INTEGER, ALLOCATABLE :: IC2SAV(:,:)
         REAL, ALLOCATABLE    :: DERWTSAV(:,:,:)
      END TYPE QA_WEIGHTS2
!
!      QA_WEIGHTS2 structure components:
!         NCSAV    Int.    Actual number of configurations of 1st and 2nd order 
!                          neighbours that have been saved in the look up table 
!         IC2SAV   I.A.    Array of configurations of 1st and 2nd order ngbrs
!                          for which weights have been saved
!         DERWTSAV R.A.    Look up table of weights for derivatives at 
!                          cell centres:
!                           DERWTSAV(ISAV,NB,IDER) = weight for "type IDER"
!                           derivative at the cell centre, from nearby cell NB,
!                           i.e. for a variable A, it's various derivatives are
!                             D(ISEA) =   sum     DERWTSAV(ISAV,NB,IDER)*A(JSEA(NB))
!                                     {NB=0:MGNBR} 
!                                        where ISAV = INDWT(ISEA), 
!                                        JSEA(NB) = {ISEA          for NB=0
!					            {GNBR(ISEA,NB) for NB>0
!                           and "D" is d/dX for IDER=1,
!                                      d/dY          2,
!                                     d2/dx2         3,
!                                     d2/dY2         4,
!                                     d2/dXdY        5,
!                                     d3/dX3         6,
!                                     d3/dY3         7,
!                                     d3/dX2dY       8,
!                                     d3/dXdY2       9

      CONTAINS
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_ADGR( QTREE, QWTS, DIAGVAR, DVMAX, DVTOLFAC,      &
                          NCTARGET, NMOD, NFROM, IFROM, WFROM, ITO,   &
                          QSPARE, QTREE_AUX, IND_AUX, EXACT_AUX,      & 
                          GET_TYPE, NO_ADAPT_TYPE, NDSE, IERR, MAPML )
!/
!         Richard Gorman, NIWA
!           February, 2014:   Origination
!
!  1. Purpose :
!
!    Adapt the quadtree grid, based on a diagnostic variable DIAGVAR.
!    Cells with DIAGVAR greater than DVMAX are refined.
!    Quads with a value of DIAGVAR (averaged over the cells within)
!    less than DVMAX/DVTOLFAC are coarsened.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ---------------------------------------------------------------
!     QTREE   QA_TREE I/O  Quadtree structure, of which the following 
!                          components are used:
!          NCELL   Int.  Number of cells
!          NCELL_DEF Int.  Number of active cells
!          NQUAD   Int.  Number of quads
!          LVLREF  Int.  Refinement level of reference grid
!          LVLMAX  Int.  Maximum refinement level
!          UNDEF_TYPE Int. value of flag for unused cell indices
!          INDLVL  I.A.  Level of each cell
!          NGBR    I.A.  Cell indices for the neighbours of input cells:
!                           NGBR(I,1:4) = index of primary (i.e. of equal or 
!                           lower level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!                           NGBR(I,5:8) = index of secondary (i.e. of higher 
!                           level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!          CELL_TYPE I.A. Flag for wet/dry/boundary/etc. cell
!          NCASE   I.A.  Index identifying the relative configuration of 
!                        neighbouring cells (out of 2401 possibilities)
!          QICELL  I.A.  Sea-point indices for the cells of the quad
!          QLEVEL  I.A.  Level of the cells within each quad
!     QWTS  QA_WEIGHT1 I   first-order weights structure
!     DIAGVAR   R.A.   I   Diagnostic cell variable
!     DVMAX     Real  I/O  Refining threshold for the diagnostic variable.
!                            Input value not used if NCTARGET>0
!     DVTOLFAC  Real   I   Ratio of refining to coarsening threshold
!     NCTARGET  Int.   I   Target number of cells. Not used if <=0
!     NMOD      Int.   O   Number of cells for which new data needs to be 
!                          assigned
!     NFROM     I.A.   O   Number of cells from which old values are taken for 
!                          each new cell
!     IFROM     I.A.   O   Indices of cells from which old values are taken
!     WFROM     R.A.   O   Weights for data at cells for from which old values
!                          are taken
!     ITO       I.A.   O   Indices of cells for which new values are computed
!     QSPARE QA_SPARE I/O* Spare cell structure
!     QTREE_AUX QA_TREE.A. I* Auxiliary quadtree structure(s)
!     IND_AUX   I.A.  I/O* Index maps from QTREE to QTREE_AUX
!     EXACT_AUX L.A.  I/O* Flag whether IND_AUX represents an exact mapping
!                          (otherwise a mapping to a parent cell)
!     GET_TYPE  L.A.   I*  Flags to derive CELL_TYPE in QTREE from CELL_TYPE 
!                          in QTREE_AUX
!     NO_ADAPT_TYPE I.A. I* Values of CELL_TYPE for cells which shouldn't 
!                           be adapted
!     NDSE      Int.   I*  Unit number for error messages
!     IERR      Int.   O*  Return error code
!     ---------------------------------------------------------------
!                       * = optional
!
!     After this subroutine is called, variables defined on quadtree cells
!     need to be mapped as follows:
!         A(ITO(k))  = sum(i=1:NFROM(k)) WFROM(k,i)*A(IFROM(k,i))
!     e.g. by calling QA_ADVAR
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     QA_REFINE  Subr. qa_utils Refine selected cells of a quadtree grid
!     QA_COARSEN Subr. qa_utils Coarsen selected quads of a quadtree grid
!     QA_QMEAN   Subr. qa_utils Compute the mean of a variable over the
!                               cells of each quad
!     QA_FINDNBR Subr. qa_utils Compute the seapoint neighbour indices
!     QA_Q2QMAP  Subr. qa_utils Compute the mapping from one quadtree 
!                               structure to another
!     QA_NBDIST  Subr. qa_utils Identify neighbour distributions for
!                               mapping to precomputed first-order
!                               gradient weights
!     QSORT      Subr. qa_utils Quicksort
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
!/ ------------------------------------------------------------------- /
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(INOUT)            :: QTREE
      TYPE(QA_WEIGHTS1), INTENT(INOUT)        :: QWTS
      REAL, INTENT(IN)                        :: DIAGVAR(:)
      REAL, INTENT(INOUT)                     :: DVMAX
      REAL, INTENT(IN)                        :: DVTOLFAC
      INTEGER, INTENT(IN)                     :: NCTARGET
      INTEGER, INTENT(OUT)                    :: NMOD
      INTEGER, INTENT(OUT)                    :: NFROM(:)
      INTEGER, INTENT(OUT)                    :: IFROM(:,:)
      REAL, INTENT(OUT)                       :: WFROM(:,:)
      INTEGER, INTENT(OUT)                    :: ITO(:)
      TYPE(QA_SPARE), OPTIONAL, INTENT(INOUT) :: QSPARE
      TYPE(QA_TREE), OPTIONAL, INTENT(IN)     :: QTREE_AUX(:)
      INTEGER, OPTIONAL, INTENT(INOUT)        :: IND_AUX(:,:)
      LOGICAL, OPTIONAL, INTENT(INOUT)        :: EXACT_AUX(:,:)
      LOGICAL, OPTIONAL, INTENT(IN)           :: GET_TYPE(:)
      INTEGER, OPTIONAL, INTENT(IN)           :: NO_ADAPT_TYPE(:)
      INTEGER, OPTIONAL, INTENT(IN)           :: NDSE
      INTEGER, OPTIONAL, INTENT(OUT)          :: IERR
      INTEGER, OPTIONAL, INTENT(IN)           :: MAPML(:)
!/
!/ ------------------------------------------------------------------- /
!/ Local variables
!/
      INTEGER, PARAMETER   :: NNGBR=8
      REAL                 :: UNDEF
      REAL                 :: DVTOLLO, DVTOLHI
      INTEGER, ALLOCATABLE :: IDVSORT(:), IDVMSORT(:)
      INTEGER, ALLOCATABLE :: QCRS(:)
      REAL, ALLOCATABLE    :: DVSORT(:), DVMSORT(:)
      REAL, ALLOCATABLE    :: DVMEAN(:)
      INTEGER              :: ISEA, NSEA, NQUAD, MSEA, NSEA_DEF
      INTEGER              :: NCUNDEF, NQUNDEF, NQOK, NCRSMAX,        &
                              NREFMAX, NREF, NSHIFT
      INTEGER              :: NCNEW, NCRS, NC, IC, ICRS, NSUM, ISEA1, &
                              ISEA2, ISEA3, IERS, IUN, IQUAD, IREF
      INTEGER              :: IS2, ISHI, ISTEP, IQ, IQS, ISORT, II,   &
                              LV0, LVN, IK, INBR, IMODE, ISHIFT
      INTEGER              :: IAUX, NAUX_GRIDS
      INTEGER              :: NT_NO_ADAPT, NC_NO_ADAPT
      INTEGER              :: NINDML
      LOGICAL              :: GTYPE
      INTEGER, ALLOCATABLE :: ISEAOLDR(:), ISEANEWR(:,:)
      INTEGER, ALLOCATABLE :: ISEAOLDC(:,:), ISEANEWC(:)
      INTEGER, ALLOCATABLE :: IQREF(:), ISREF(:)
      INTEGER, ALLOCATABLE :: OLDSEA(:), NEWSEA(:)
      INTEGER, ALLOCATABLE :: MAPMLT(:)
      REAL                 :: WXNBR(4), WYNBR(4)
!
!     X,Y Coordinates of the four child cells in a quad:
!
      WXNBR = (/-0.25,+0.25,-0.25,+0.25/)
      WYNBR = (/-0.25,-0.25,+0.25,+0.25/)
!
!     Initialise outputs:
!
      NFROM = 0
      IFROM = 0
      WFROM = 0.
      ITO = 0
!
!     Defaults for optional inputs:
!
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(QTREE_AUX) .AND. PRESENT(IND_AUX) .AND.            &
           PRESENT(EXACT_AUX) ) THEN
         NAUX_GRIDS = MIN( SIZE(QTREE_AUX,1), SIZE(IND_AUX,2) )
         NAUX_GRIDS = MIN( NAUX_GRIDS, SIZE(EXACT_AUX,2) )
      ELSE
         NAUX_GRIDS = 0
      END IF
      NT_NO_ADAPT =0
      IF ( PRESENT(NO_ADAPT_TYPE) ) NT_NO_ADAPT = SIZE(NO_ADAPT_TYPE,1) 
      IF ( PRESENT(IERR) ) IERR = 0
      NINDML = 1
      IF ( PRESENT(MAPML) ) NINDML = SIZE(MAPML,1)
      ALLOCATE ( MAPMLT(NINDML) )
      IF ( PRESENT(MAPML) ) THEN
         MAPMLT = MAPML
      ELSE
         MAPMLT = 1
      END IF
!
!     Allocate arrays:
!
      MSEA = SIZE(DIAGVAR,1)
      NSEA = QTREE%NCELL
      NSEA_DEF = QTREE%NCELL_DEF
      NQUAD = QTREE%NQUAD
      !WRITE(*,*) 'QA_ADGR STARTING, NSEA, NSEA_DEF, NQUAD = ',        &
      !           NSEA, NSEA_DEF, NQUAD
      ALLOCATE ( ISEAOLDR(MSEA), ISEANEWR(MSEA,4) )
      ALLOCATE ( ISEAOLDC(MSEA,4), ISEANEWC(MSEA) )
      ALLOCATE ( IQREF(MSEA), ISREF(MSEA) )
      ALLOCATE ( OLDSEA(MSEA), NEWSEA(MSEA) )
!
! Find the number of valid cells, and sort the diagnostic variable array, 
! assuming that invalid values will be negative
!
      ALLOCATE (DVSORT(NSEA))
      ALLOCATE (IDVSORT(NSEA))
      DVSORT = DIAGVAR(1:NSEA)
      UNDEF = MIN(0.,MINVAL(DVSORT)) - 1.
      DO ISEA=1,NSEA
         IDVSORT(ISEA) = ISEA
         IF ( QTREE%CELL_TYPE(ISEA) .EQ. QTREE%UNDEF_TYPE )           &
            DVSORT(ISEA) = UNDEF
      END DO
      NCUNDEF = NSEA - NSEA_DEF
      CALL QSORT(DVSORT,IDVSORT,NSEA)
!
! Get the mean value of the diagnostic variable over each quad:
!
      ALLOCATE (DVMEAN(NQUAD))
      ALLOCATE (DVMSORT(NQUAD))
      ALLOCATE (IDVMSORT(NQUAD))
      ALLOCATE (QCRS(NQUAD))
      CALL QA_QMEAN ( QTREE, array=DIAGVAR, amean=DVMEAN,             &
                      a_undef=UNDEF )
      !WRITE(*,*) 'QUADMEAN done'
!
! Find the number of valid quads (NQOK)
!
      DVMSORT = DVMEAN(1:NQUAD)
      NQUNDEF = 0
      DO IQUAD=1,NQUAD
         IDVMSORT(IQUAD) = IQUAD
         IF ( ABS(DVMEAN(IQUAD)-UNDEF).LT.TINY(UNDEF) )               &
            NQUNDEF = NQUNDEF+1
      END DO
      NQOK = NQUAD - NQUNDEF
!
! Sort the quad means, assuming that invalid
! values will be negative
!
      CALL QSORT(DVMSORT,IDVMSORT,NQUAD)
!
! Find tolerances to adapt the grid by one of two methods.
!
      IF ( NCTARGET.GT.0 ) THEN
!
! Method 1: use a target number of cells, and adjust DVTOLLO and DVTOLHI
!           (in a fixed ratio) to achieve that
!
         ! Loop in ascending order of DVMSORT:
         ISHI = NSEA
         IS2 = 1
         ISTEP = -1
         DO IQS=NQUNDEF+1,NQUAD
            IQ = IDVMSORT(IQS)
            DVTOLLO = DVMSORT(IQS)
            ! If we coarsened all quads with DVMEAN less than this,
            ! the number of coarsened quads and the upper (refining)
            ! tolerance would be:
            NCRSMAX = IQS - NQUNDEF - 1
            DVTOLHI = DVTOLLO*DVTOLFAC
            ! Find the cell index (ISHI) corresponding to that refining 
            ! threshold by searching upward (except the first time)
            DO ISORT=ISHI,IS2,ISTEP
               IF ( DVSORT(ISORT).LE.DVTOLHI ) THEN
                  ISHI = ISORT
                  EXIT
               END IF
            END DO
            IS2=NSEA
            ISTEP=1
            ! If we refined all cells with DIAGVAR greater than this,
            ! the number of refined cells would be:
            NREFMAX = NSEA-ISHI
            !
            ! Total number of cells expected after adapting:
            NCNEW = NSEA_DEF + 3*(NREFMAX-NCRSMAX)
            !
            ! Stop when we've reached the target:
            IF ( NCNEW.LE.NCTARGET ) EXIT
         END DO
         DVMAX = DVTOLHI
      ELSE
!
! Method 2: use fixed tolerances DVTOLLO and DVTOLHI = DVTOLFAC*DVTOLLO
!           DVTOLFAC is fixed
         DVTOLHI = DVMAX
         DVTOLLO = DVTOLHI/DVTOLFAC
         NCRSMAX = NQOK
         NREFMAX = NSEA_DEF
      END IF
!
! Loop through quads in increasing order of quad-mean diagnostic variable
! Identify quads with estimated mean < tolerance
!
      NCRS = 0
      QCRS = 0
      DO IQS=NQUNDEF+1,NQUAD
         IQ = IDVMSORT(IQS)        
         IF ( DVMEAN(IQ).GT.DVTOLLO ) EXIT
         !
         ! Skip quads with children, and level-0 leaf cells
         IF (QTREE%QLEVEL(IQ).EQ.1 .AND. QTREE%QICELL(IQ,0).GT.0)     &
           CYCLE
         NC = 0
         DO II=1,4
            IF ( QTREE%QCHILD(IQ,II).GT.0 ) NC = NC + 1
         END DO
         IF ( NC.GT.0 ) CYCLE
         !
         ! Find the highest level of neighbouring cells
         ! and the number of "non-adaptable" cells in the quad
         LV0 = QTREE%QLEVEL(IQ)
         LVN = LV0
         NC_NO_ADAPT = 0
         DO II=0,4
            ISEA = QTREE%QICELL(IQ,II)
            IF ( ISEA.GT.0 ) THEN
               IF ( NT_NO_ADAPT.GT.0 ) THEN
                  IF ( ANY(NO_ADAPT_TYPE.EQ.QTREE%CELL_TYPE(ISEA)) )  &
                     NC_NO_ADAPT = NC_NO_ADAPT + 1
               END IF
               DO IK=1,NNGBR
                  INBR = QTREE%NGBR(ISEA,IK)
                  IF ( INBR.GT.0 ) LVN = MAX( QTREE%INDLVL(INBR),LVN )
               END DO 
            END IF
         END DO
         !
         ! Only coarsen if there are no finer neighbours, and no 
         ! "non-adaptable" cell types
         IF ( LVN.GT.LV0 .OR. NC_NO_ADAPT.GT.0 ) CYCLE
         !
         NCRS = NCRS + 1
         QCRS(NCRS) = IQ
         !WRITE(*,*) 'COARSEN QUAD ', IQ
         !
         IF ( NCTARGET.GT.0 ) THEN
            IF ( NCRS.GE.NCRSMAX ) EXIT
         ELSE
            IF ( DVMEAN(IQ).GT.DVTOLLO ) EXIT
         END IF
      END DO
      !WRITE(*,*) 'QUADS TO COARSEN = ', NCRS
!
!  Coarsen flagged quads
!     TO CHECK/FIX: use of QSPARE in QA_COARSEN
!
      IMODE = 1
      IF ( PRESENT(QSPARE) ) THEN
         CALL QA_COARSEN ( IMODE, QTREE, NCRS, QCRS, ISEAOLDC,        &
                           ISEANEWC, NSHIFT, qspare=QSPARE,           &
                           mapml=MAPMLT )
      ELSE
         CALL QA_COARSEN ( IMODE, QTREE, NCRS, QCRS, ISEAOLDC,        &
                           ISEANEWC, NSHIFT, mapml=MAPMLT )
      END IF
      NSEA = QTREE%NCELL
      NSEA_DEF = QTREE%NCELL_DEF
      NQUAD = QTREE%NQUAD
      !WRITE(*,*) 'QA_COARSEN DONE, NCRS, NSEA, NSEA_DEF, NQUAD = ',   &
      !           NCRS, NSEA, NSEA_DEF, NQUAD
      !QTREE%CELL_TYPE = 1
!
! Indices and weights for computing mean quantities from the coarsened quad 
! to apply to the single cells that replace them
!
      DO ICRS=1,NCRS
         ISEA2 = ISEANEWC(ICRS)
         ITO(ICRS) = ISEA2
         NSUM = 0
         DO II=1,4
            ISEA1 = ISEAOLDC(ICRS,II)
            IF ( ISEA1.GT.0 ) THEN
               NSUM = NSUM + 1
               IFROM(ICRS,NSUM) = ISEA1
               WFROM(ICRS,NSUM) = 1.
            END IF
         END DO
         NFROM(ICRS) = NSUM
         IF ( NSUM.GT.1 ) WFROM(ICRS,:) = WFROM(ICRS,:)/NSUM
      END DO
      !
      ! Cell renumbering during coarsening:
      DO ISEA = 1,NSEA
         NEWSEA(ISEA) = ISEA
         OLDSEA(ISEA) = ISEA
      END DO
      IF ( NSHIFT.GT.0 ) THEN
         DO ISHIFT=1,NSHIFT
            DO II=1,4
               ISEA1 = ISEAOLDC(NCRS+ISHIFT,II)
               IF ( ISEA1.GT.0 ) EXIT
            END DO
            ISEA2 = ISEANEWC(NCRS+ISHIFT)
            NEWSEA(ISEA1) = ISEA2
            OLDSEA(ISEA2) = ISEA1
            ITO(NCRS+ISHIFT) = ISEA2
            IFROM(NCRS+ISHIFT,1) = ISEA1
            WFROM(NCRS+ISHIFT,1) = 1.
            NFROM(NCRS+ISHIFT) = 1
         END DO
      END IF
!
! Identify cells with diagnostic variable > tolerance
!
      NREF = 0
      IQREF = 0
      ISREF = 0
      ISEAOLDR = 0
      NREFMAX = MIN(NREFMAX,(MSEA-NSEA)/3)
      DO ISORT=NSEA,1,-1
         ISEA1 = IDVSORT(ISORT)  ! original cell index
         ISEA = NEWSEA(ISEA1)    ! remapped cell index after coarsening
         IF ( ISEA.LE.0 .OR. ISEA.GT.NSEA ) CYCLE
         !  Don't refine if the cell is of a type that shouldn't be 
         !  adapted:
         IF ( QTREE%CELL_TYPE(ISEA).EQ.QTREE%UNDEF_TYPE ) CYCLE
         IF ( NT_NO_ADAPT.GT.0 ) THEN
            IF ( ANY(NO_ADAPT_TYPE.EQ.QTREE%CELL_TYPE(ISEA)) ) CYCLE
         END IF
         IF ( DIAGVAR(ISEA1).LT.DVTOLHI ) EXIT
         ! Don't refine cells that are already at the maximum level
         LV0 = QTREE%INDLVL(ISEA)
         IF ( LV0.EQ.QTREE%LVLMAX ) CYCLE
         !QA_ADGR
         ! Find the lowest level of neighbouring cells
         LVN = LV0
         DO IK=1,NNGBR
            INBR = QTREE%NGBR(ISEA,IK)
            IF ( INBR.GT.0 ) LVN = MIN( QTREE%INDLVL(INBR),LVN )
         END DO 
         !
         ! Only refine if there are no coarser neighbours
         IF ( LVN.LT.LV0 ) CYCLE
         NREF = NREF + 1
         ISEAOLDR(NREF) = ISEA
         !WRITE(*,*) 'REFINE CELL ', ISEA
         !QA_ADGR
         IF ( NCTARGET.GT.0 ) THEN
            IF ( NREF.GE.NREFMAX ) EXIT
         ELSE
            IF ( DIAGVAR(ISEA1).LT.DVTOLLO ) EXIT
         END IF
      END DO
      !WRITE(*,*) 'CELLS TO REFINE = ', NREF
!
! Refine flagged cells
!     TO CHECK/FIX: use of QSPARE in QA_REFINE
!
      IMODE = 2
      IF ( PRESENT(QSPARE) ) THEN
         CALL QA_REFINE ( IMODE, QTREE, NREF, IQREF, ISREF,           &
                          ISEAOLDR, ISEANEWR, ierr=IERS, ndse=IUN,    &
                          qspare=QSPARE, mapml=MAPMLT )
      ELSE
         CALL QA_REFINE ( IMODE, QTREE, NREF, IQREF, ISREF,           &
                          ISEAOLDR, ISEANEWR, mapml=MAPMLT,           &
                          ierr=IERS, ndse=IUN )
      END IF
      NSEA = QTREE%NCELL
      NSEA_DEF = QTREE%NCELL_DEF
      NQUAD = QTREE%NQUAD
      IF ( IERS.GT.0 ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_ADGR CALLING QA_REFINE'
         END IF
         IF ( PRESENT(IERR) ) IERR = IERS
         RETURN
      END IF
      !
      !WRITE(*,*) 'QA_REFINE DONE, NREF, NSEA, NSEA_DEF, NQUAD = ',    &
      !           NREF, NSEA, NSEA_DEF, NQUAD
      NMOD = NCRS + NSHIFT
      DO IREF=1,NREF
         ISEA = ISEAOLDR(IREF)
         ISEA1 = OLDSEA(ISEA)
! Weights for computing gradients
         IC = QTREE%NCASE(ISEA)
! Apply input & precomputed data to the new cells
         
         DO II=1,4
            ISEA2 = ISEANEWR(IREF,II)
            IF (ISEA2.GT.0) THEN
               NMOD = NMOD + 1
               ITO(NMOD) = ISEA2
               NSUM = 0
               DO INBR=0,NNGBR
                  IF ( INBR.EQ.0 ) THEN
                     ISEA3 = ISEA
                  ELSE
                     ISEA3 = QTREE%NGBR(ISEA,INBR)
                  END IF
                  IF ( ISEA3.GT.0 ) THEN
                     ISEA3 = OLDSEA(ISEA3)
                     IF ( ISEA3.GT.0 ) THEN
                        NSUM = NSUM + 1
                        IFROM(NMOD,NSUM) = ISEA3
                        WFROM(NMOD,NSUM) =                         &
                                WXNBR(II)*QWTS%DERWTC(IC,INBR,1)   &
                              + WYNBR(II)*QWTS%DERWTC(IC,INBR,2)
                     END IF
                  END IF
               END DO
               WFROM(NMOD,1) = WFROM(NMOD,1) + 1.
               NFROM(NMOD) = NSUM
            END IF
         END DO
      END DO
!
! Recompute mappings to auxiliary grids, if any:
!
      IF ( NAUX_GRIDS.GT.0 ) THEN
         DO IAUX = 1,NAUX_GRIDS
            IF ( PRESENT(GET_TYPE) .AND.                              &
                 IAUX.GE.SIZE(GET_TYPE,1) ) THEN
                GTYPE = GET_TYPE(IAUX)
            ELSE
                GTYPE = .FALSE.
            END IF
            CALL QA_Q2QMAP( QTREE, QTREE_AUX(IAUX), GTYPE,            &
                            IND_AUX(:,IAUX),                          &
                            exact12=EXACT_AUX(:,IAUX), ncalc=NMOD,    &
                            icalc=ITO(1:NMOD), ierr=IERS, ndse=IUN )
            IF ( IERS.GT.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE(IUN,*) 'ERROR IN QA_ADGR CALLING QA_Q2QMAP'
                  WRITE(IUN,*) 'FOR AUXILIARY GRID ', IAUX
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
         END DO
      END IF
!
! Recompute cell neighbours
!    TO BE CHECKED/FIXED: only recalculate neighbours for cells that 
!    have been modified. Will need to check "reverse neighbour" adjustment
!    within QA_FINDNBR. Meanwhile, recalculate neighbours for ALL cells:
!
      !CALL QA_FINDNBR ( QTREE, ierr=IERS, ndse=IUN, ncalc=NMOD,       &
      !                  icalc=ITO(1:NMOD) )
      CALL QA_FINDNBR ( QTREE, ierr=IERS, ndse=IUN )
      IF ( IERS.GT.0 ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_ADGR CALLING QA_FINDNBR'
         END IF
         IF ( PRESENT(IERR) ) IERR = IERS
         RETURN
      END IF
      !WRITE(*,*) 'FINDNBR done'
      !WRITE(*,*) 'NSEA, NQUAD = ', NSEA, NQUAD
!
! Recompute weights
      CALL QA_NBDIST( QTREE, ncalc=NMOD, icalc=ITO(1:NMOD) )
      NSEA = QTREE%NCELL
      NSEA_DEF = QTREE%NCELL_DEF
      NQUAD = QTREE%NQUAD
      !WRITE(*,*) 'QA_ADGR ENDING, NSEA, NSEA_DEF, NQUAD = ',          &
      !           NSEA, NSEA_DEF, NQUAD
!
      IF ( PRESENT(IERR) ) IERR = 0
      RETURN
      END SUBROUTINE QA_ADGR
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_ADVAR( NMOD, NFROM, IFROM, WFROM, ITO, A )
!/
!         Richard Gorman, NIWA
!           May, 2014:   Origination
!
!  1. Purpose :
!
!    Recompute an array after adapting a quadtree grid.
!
!  2. Method :
!         A(ITO(k))  = sum(i=1:NFROM(k)) WFROM(k,i)*A(IFROM(k,i))
!
!  3. Parameters :
!
!     Parameter list
!     ---------------------------------------------------------------
!     NMOD      Int.   I   Number of cells for which new data needs to be 
!                          assigned
!     NFROM     I.A.   I   Number of cells from which old values are taken for 
!                          each new cell
!     IFROM     I.A.   I   Indices of cells from which old values are taken
!     WFROM     R.A.   I   Weights for data at cells for from which old values
!                          are taken
!     ITO       I.A.   I   Indices of cells for which new values are computed
!     A         R.A.  I/O  Array to be recomputed
!     ---------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)                    :: NMOD
      INTEGER, INTENT(IN)                    :: NFROM(:)
      INTEGER, INTENT(IN)                    :: IFROM(:,:)
      REAL, INTENT(IN)                       :: WFROM(:,:)
      INTEGER, INTENT(IN)                    :: ITO(:)
      REAL, INTENT(INOUT)                    :: A(:)
!/
!/ ------------------------------------------------------------------- /
!/ Local variables
!/
      REAL    :: ANEW(NMOD)
      INTEGER :: IMOD, ISEA, II, ISEA2
!
!  First loop to calculate new values by weighted sum of old values
!
      ANEW = 0.
      DO IMOD=1,NMOD
         !ISEA = ITO(IMOD)
         DO II=1,NFROM(IMOD)
            ISEA2 = IFROM(IMOD,II)
            ANEW(IMOD) = ANEW(IMOD) + WFROM(IMOD,II)*A(ISEA2)
         END DO
      END DO
!
!  Second loop to replace old with new values
!      
      DO IMOD=1,NMOD
         ISEA = ITO(IMOD)
         A(ISEA) = ANEW(IMOD)
      END DO
      RETURN
      END SUBROUTINE QA_ADVAR
!/
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_ALLOC ( QTREE, MCELL, MQUAD, UNDEF_TYPE,          &
                            IWTORDER, QSPARE, MPART, MCLOC )
!/
!/       Richard Gorman, NIWA
!/         Jan-April, 2014:      Origination.
!
!  1. Purpose :
!
!      Allocate arrays for a threaded quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE   QA_TREE  I/O Quadtree structure, with the following components
!                           affected:
!          NCELL   Int.  Number of cells
!          NCELL_DEF Int.  Number of active cells
!          NQUAD   Int.  Number of quads
!          LVLREF  Int.  Refinement level of reference grid
!          LVLMAX  Int.  Maximum refinement level allowed
!          LVLHI   Int.  Maximum refinement level actually present
!          NX0     Int.  X dimension of level-zero base grid
!          NY0     Int.  Y dimension of level-zero base grid
!          UNDEF_TYPE Int. value of flag for unused cell indices
!          KEEP_REF Log. retain cell indices for refined cells
!          INDQUAD I.A.  Index of the quad containing each input cell
!          INDSUB  I.A.  Index of the sub-quad containing each cell:
!                               0 = none (home cell is at level zero)
!                               1,2,3,4 = SW, SE, NW, NE corner of quad
!          INDLVL  I.A.  Level of each cell
!          NGBR    I.A.  Sea-point indices for the neighbours of input cells
!                           NGBR(I,1:4) = index of primary (i.e. of equal or 
!                           lower level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!                           NGBR(I,5:8) = index of secondary (i.e. of higher 
!                           level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!          CELL_TYPE I.A. Flag for wet/dry/boundary/etc. cell
!          INDML   I.A.  Multilevel index of cells
!          NCASE   I.A.  Index identifying the relative configuration of 
!                        neighbouring cells (out of 2401 possibilities)
!          XYVAL   R.A.  Coordinates of cells (w.r.t. 
!                          reference rectangular grid)
!          QICELL  I.A.  Sea-point indices for the cells of the quad
!          QLEVEL  I.A.  Level of the cells within each quad
!          QPARENT I.A.  Quad index for the parent of each quad
!          QNBR    I.A.  Quad indices for the neighbours of each quad
!          QCHILD  I.A.  Quad indices for the children of each quad
!          INDWT   I.A.  Index, for each cell, into the lookup table of
!                           (2nd order) weights
!          GNBR    I.A.  Grand neighbours of each cell
!          NGNBR   I.A.  Number of grand neighbours of each cell
!       MCELL      Int.   I   Allocated max. number of cells
!       MQUAD      Int.   I   Allocated max. number of quads
!       UNDEF_TYPE Int.   I* Value of flag for unused cell indices
!       IWTORDER   Int.   I* Type of interpolation weights used: 
!                             0=none,1=1st order,2=2nd order 
!       QSPARE  QA_SPARE  O* Arrays for treating spare cells in quadtree 
!                            structure, with the following components used:
!          NPART   Int.   I   Number of partitions
!          PARTNO  I.A.  I/O  Partition number for each cell
!          JSVAL   I.A.  I/O  Local cell index for each global index
!          NCLOC   I.A.  I/O  Actual number of cells in each partition
!          ISVAL   I.A.  I/O  Global cell index for each local index & partition
!          NSPARE  I.A.  I/O  Number of spare local cell indices
!          JSSPARE I.A.  I/O  List of spare local cell indices
!       MPART   Int.     I*  Allocated max. number of partitions
!       MCLOC   Int.     I*  Allocated max. number of cells per partition
!     ----------------------------------------------------------------
!                         * = optional
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_BMLSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC   Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(INOUT)          :: QTREE
      INTEGER, INTENT(IN)                   :: MCELL
      INTEGER, INTENT(IN)                   :: MQUAD
      INTEGER, INTENT(IN), OPTIONAL         :: UNDEF_TYPE
      INTEGER, INTENT(IN), OPTIONAL         :: IWTORDER
      TYPE(QA_SPARE), OPTIONAL, INTENT(OUT) :: QSPARE
      INTEGER, OPTIONAL, INTENT(IN)         :: MPART
      INTEGER, OPTIONAL, INTENT(IN)         :: MCLOC
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
      INTEGER   :: NPART, UNDEF, IORD, MCL, IPART, ICELL, JCELL
      INTEGER   :: MC1, MC2, MN2
!/
      IF ( PRESENT(UNDEF_TYPE) ) THEN
         UNDEF = UNDEF_TYPE
      ELSE
         UNDEF = 0
      END IF
      IF ( PRESENT(IWTORDER) ) THEN
         IORD = IWTORDER
      ELSE
         IORD = 0
      END IF
!
!  QTREE scalars:
!    Initialise, as empty
!
      QTREE%NCELL = 0
      QTREE%NCELL_DEF = 0
      QTREE%NQUAD = 0
      QTREE%NX0 = 0
      QTREE%NY0 = 0
      QTREE%LVLREF = 0
      QTREE%LVLMAX = 0
      QTREE%LVLHI = 0
      QTREE%KEEP_REF = .FALSE.
      QTREE%DYNAMIC = .TRUE.
      QTREE%UNDEF_TYPE = UNDEF
      QTREE%IWTORDER = IORD
!
!  QTREE QUAD allocatable arrays:
!
      IF ( ALLOCATED (QTREE%QLEVEL) )  DEALLOCATE ( QTREE%QLEVEL )
      IF ( ALLOCATED (QTREE%QNBR) )    DEALLOCATE ( QTREE%QNBR )
      IF ( ALLOCATED (QTREE%QPARENT) ) DEALLOCATE ( QTREE%QPARENT )
      IF ( ALLOCATED (QTREE%QCHILD) )  DEALLOCATE ( QTREE%QCHILD )
      IF ( ALLOCATED (QTREE%QICELL) )  DEALLOCATE ( QTREE%QICELL )
!
      IF ( MQUAD.GT.0 ) THEN
!
!    Calling with MQUAD<=0 will just deallocate.
!    Otherwise, allocate with the specified size ...
!
         ALLOCATE ( QTREE%QLEVEL(MQUAD) )
         ALLOCATE ( QTREE%QNBR(MQUAD,4) )
         ALLOCATE ( QTREE%QPARENT(MQUAD) )
         ALLOCATE ( QTREE%QCHILD(MQUAD,4) )
         ALLOCATE ( QTREE%QICELL(MQUAD,0:4) )
!
!    ... and initialise as null
!
         QTREE%QLEVEL = 0
         QTREE%QNBR = 0
         QTREE%QPARENT = 0
         QTREE%QCHILD = 0
         QTREE%QICELL = 0
      END IF
!
!  QTREE CELL pointer arrays:
!
      IF ( ALLOCATED (QTREE%INDQUAD) )   DEALLOCATE (QTREE%INDQUAD)
      IF ( ALLOCATED (QTREE%INDSUB) )    DEALLOCATE (QTREE%INDSUB)
      IF ( ALLOCATED (QTREE%INDLVL) )    DEALLOCATE (QTREE%INDLVL)
      IF ( ALLOCATED (QTREE%NGBR) )      DEALLOCATE (QTREE%NGBR)
      IF ( ALLOCATED (QTREE%CELL_TYPE) ) DEALLOCATE (QTREE%CELL_TYPE)
      IF ( ALLOCATED (QTREE%INDML) )     DEALLOCATE (QTREE%INDML)
      IF ( ALLOCATED (QTREE%XYVAL) )     DEALLOCATE (QTREE%XYVAL)
!
      IF ( ALLOCATED (QTREE%NCASE) )     DEALLOCATE (QTREE%NCASE)
!
      IF ( ALLOCATED( QTREE%NGNBR) )     DEALLOCATE (QTREE%NGNBR)
      IF ( ALLOCATED( QTREE%GNBR) )      DEALLOCATE (QTREE%GNBR)
      IF ( ALLOCATED( QTREE%INDWT) )     DEALLOCATE (QTREE%INDWT)
!
      IF ( MCELL.GT.0 ) THEN
!
!    Calling with MCELL<=0 will just deallocate.
!    Otherwise, allocate with the specified size ...
!
         ALLOCATE ( QTREE%INDQUAD(MCELL) )
         ALLOCATE ( QTREE%INDSUB(MCELL) )
         ALLOCATE ( QTREE%INDLVL(MCELL) )
         ALLOCATE ( QTREE%NGBR(MCELL,8) )
         ALLOCATE ( QTREE%CELL_TYPE(MCELL) )
         ALLOCATE ( QTREE%INDML(MCELL) )
         ALLOCATE ( QTREE%XYVAL(MCELL,2) )
!
!    If weight arrays are needed, allocate them with the specified size
!    otherwise with a dummy size 1
!
         IF ( IORD.GE.1 ) THEN
            MC1 = MCELL
         ELSE
            MC1 = 1
         END IF
         ALLOCATE ( QTREE%NCASE(MC1) )
!
         IF ( IORD.GE.2 ) THEN
            MC2 = MCELL
            MN2 = MCELL
         ELSE
            MC2 = 1
            MN2 = 1
         END IF
         ALLOCATE ( QTREE%NGNBR(MC2) )
         ALLOCATE ( QTREE%GNBR(MC2,MN2) )
         ALLOCATE ( QTREE%INDWT(MC2) )
!
!    ... and initialise all quadtree cell arrays as null
!
         QTREE%INDQUAD = 0
         QTREE%INDSUB = 0
         QTREE%INDLVL = 0
         QTREE%NGBR = 0
         QTREE%CELL_TYPE = UNDEF
         QTREE%INDML = 0
         QTREE%XYVAL = 0.
!
         QTREE%NCASE = 0
!
         QTREE%NGNBR = 0
         QTREE%GNBR = 0
         QTREE%INDWT = 0
      END IF
!      
      IF ( .NOT.PRESENT(QSPARE) ) RETURN
!      
      IF ( PRESENT(MPART) ) THEN
         NPART = MPART
      ELSE
         NPART = 1
      END IF
      IF ( PRESENT(MCLOC) ) THEN
         MCL = MCLOC
      ELSE
         MCL = CEILING(FLOAT(MCELL)/FLOAT(NPART))
      END IF
!
      QSPARE%NPART = NPART
 
      IF ( ALLOCATED( QSPARE%PARTNO ) ) DEALLOCATE( QSPARE%PARTNO )
      IF ( ALLOCATED( QSPARE%JSVAL ) ) DEALLOCATE( QSPARE%JSVAL )
      IF ( ALLOCATED( QSPARE%NSPARE ) ) DEALLOCATE( QSPARE%NSPARE )
      IF ( ALLOCATED( QSPARE%NCLOC ) ) DEALLOCATE( QSPARE%NCLOC )
      IF ( ALLOCATED( QSPARE%ISVAL ) ) DEALLOCATE( QSPARE%ISVAL )
      IF ( ALLOCATED( QSPARE%JSSPARE ) ) DEALLOCATE( QSPARE%JSSPARE )
!
      IF ( MCELL.GT.0 ) THEN
         ALLOCATE( QSPARE%PARTNO(MCELL) )
         ALLOCATE( QSPARE%JSVAL(MCELL) )
      END IF
      IF ( NPART.GT.0 ) THEN
         ALLOCATE( QSPARE%NSPARE(NPART) )
         QSPARE%NSPARE = 0
         ALLOCATE( QSPARE%NCLOC(NPART) )
         IF ( MCL.GT.0 ) THEN
            ALLOCATE( QSPARE%ISVAL(MCL,NPART) )
            ALLOCATE( QSPARE%JSSPARE(MCL,NPART) )
            QSPARE%JSSPARE = 0
         END IF
      END IF
!
!  Initialise      
!
      IF ( MCELL.GT.0 .AND. NPART.EQ.1 ) THEN
         !
         ! Trivial partitioning scheme for 1 partition 
         QSPARE%PARTNO = 1
         QSPARE%NCLOC = MCELL
         !
         IPART = 1
         DO ICELL=1,MCELL
            QSPARE%JSVAL(ICELL) = ICELL
            QSPARE%ISVAL(ICELL,IPART) = ICELL
         END DO
      ELSEIF ( MCELL.GT.0 .AND. NPART.GT.1 ) THEN
         !
         ! Wavewatch partitioning scheme: may override later
         IPART = 0
         JCELL = 1
         QSPARE%NCLOC = 0
         DO ICELL=1,MCELL
            IPART = IPART + 1
            IF ( IPART.GT.NPART ) THEN
               IPART = 1
               JCELL = JCELL + 1
            END IF
            QSPARE%JSVAL(ICELL) = JCELL
            QSPARE%ISVAL(JCELL,IPART) = ICELL
            QSPARE%PARTNO(ICELL) = IPART
            QSPARE%NCLOC(IPART) = QSPARE%NCLOC(IPART) + 1
         END DO
      END IF
!
      RETURN
      END SUBROUTINE QA_ALLOC
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_BMLSTRUC( NX, NY, NMAXB, LVLREF, LVLMAX, MAPSTA,  &
                              SIG1, GLOBAL, INDBML, ZBML, ISWET,      &
                              NSEA, NQUAD, IERR, NDSE, QTREE, IBML,   &
                              ZLEV )
!/
!       Richard Gorman, NIWA
!         September, 2012:  Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!      Create a threaded quadtree structure based on bathymetry
!      defined on a multilevel grid structure, with a limited number of
!      cells produced by coarsening where a "CFL factor" = 
!      CG/cell width is greatest
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ---------------------------------------------------------------
!       NX, NY  Int.   I   Size of reference grid
!       NMAXB   Int.   I   Size of multilevel bathymetry grid
!       LVLREF  Int.   I   Refinement level of reference grid
!       LVLMAX  Int.   I   Maximum refinement level allowed
!       MAPSTA  I.A.   I   Flag for sea/boundary/excluded points 
!                          (on rectangular grid)
!       SIG1    Real   I   Radian frequency of lowest bin
!       GLOBAL  Log.   I   Flag for global grid (i.e. cyclic in X)
!       INDBML  I.A.   I   Index from bathymetry multilevel grid to 
!                          full multilevel grid
!       ZBML    R.A.   I   Elevation of the seabed on multilevel 
!                          bathymetry grid (m above datum)
!       NSEA    Int.   O   Actual number of cells
!       NQUAD   Int.   O   Actual number of quads
!       IERR    Int.   O*  Return flag = 1 for error, else 0
!       NDSE    Int.   I*  Unit number for error output (if >0)
!       QTREE   QA_TREE O* Quadtree structure, of which the following 
!                          components are affected:
!         NCELL   Int.       Number of cells (including undefined cells)
!         NCELL_DEF Int.     Number of active cells
!         NQUAD   Int.       Number of quads
!         LVLMAX  Int.       Maximum refinement level allowed
!         LVLHI   Int.       Maximum refinement level actually present
!         INDLVL  I.A.       Level of each cell
!         INDML   I.A.     
!         NGBR    I.A.       Neighbours (up to 8) of each cell
!         XYVAL   R.A.       Coordinates of cells (w.r.t.  
!                            reference grid)
!         CELL_TYPE  I.A.    Flag for sea/boundary/excluded points
!         INDQUAD I.A.       Index of the quad containing each input 
!                            cell
!         INDSUB  I.A.       Index of the sub-quad containing each sea
!                            point:
!                              0 = none (home cell is at level 0) 
!                              1,2,3,4 = SW, SE, NW, NE corner of quad
!         QICELL  I.A.       Sea-point indices for the cells of the quad
!         QLEVEL  I.A.       Level of the cells within each quad
!         QPARENT I.A.       Quad index for the parent of each quad
!         QNBR    I.A.       Quad indices for the neighbours of each quad
!         QCHILD  I.A.       Quad indices for the children of each quad
!       IBML    I.A.   O*  Index to multilevel bathy grid at each cell
!       ZLEV    R.A.   O*  Elevation of the seabed at each sea
!                            point (m above datum)
!     ---------------------------------------------------------------
!                       * = OPTIONAL
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     QA_ALLOC   Subr. qa_utils Allocate a quadtree structure
!     QA_RECT    Subr. qa_utils Set up a quadtree representation (level 0)
!     QA_REFINE  Subr. qa_utils Refine selected cells of a quadtree grid
!     QA_COARSEN Subr. qa_utils Coarsen selected quads of a quadtree grid
!     QA_QMEAN   Subr. qa_utils Compute the mean of a variable over 
!                               the cells of each quad
!     QA_FINDNBR Subr. qa_utils Compute the seapoint neighbour indices
!     ML_CHILD   Func. qa_utils Multilevel index of the SW child of a given cell
!     QA_DISPREL Subr. qa_utils Compute CG using the dispersion relation
!     ----------------------------------------------------------------
!
!  5. Called by :
!
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
      INTEGER, INTENT(IN)   :: NX, NY, NMAXB, LVLREF, LVLMAX
      REAL, INTENT(IN)      :: SIG1
      INTEGER, INTENT(IN)   :: MAPSTA(NY,NX)
      LOGICAL, INTENT(IN)   :: GLOBAL
      INTEGER, INTENT(IN)   :: INDBML(NMAXB)
      REAL, INTENT(IN)      :: ZBML(NMAXB)
      LOGICAL, INTENT(IN)   :: ISWET(NMAXB)
      INTEGER, INTENT(OUT)  :: NSEA, NQUAD
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
      TYPE(QA_TREE), INTENT(INOUT), OPTIONAL  :: QTREE
      INTEGER, INTENT(OUT), OPTIONAL  :: IBML(:)
      REAL, INTENT(OUT), OPTIONAL  :: ZLEV(:)
!
! Local variables.
!
!  Temporary version of the output QTREE, declared for the
!  full number of seapoints on the finest grid
      TYPE(QA_TREE)  :: QTREE_TEMP
!  Temporary versions of the optional output quad variables
      INTEGER              :: NSEA2
!
      INTEGER              :: MSEA, MQUAD, NSEA_DEF
      REAL                 :: CFL(NMAXB), ZLEVT(NMAXB)
      INTEGER              :: IBMLT(NMAXB)
      INTEGER              :: ISREF(NMAXB), ISEAOLDR(NMAXB),          &
                              ISEANEWR(NMAXB,4)
      INTEGER              :: ISEAOLDC(NMAXB,4), ISEANEWC(NMAXB)
      INTEGER              :: IQREF(NMAXB)
      INTEGER              :: NX0, NY0, ILEV
      INTEGER              :: IQ, ISEA, IMODE, NREF, ISUB,            &
                              NCRS, IX, IY, NSEATEM, MAPTMP
      REAL                 :: DELX0, DELY0, ZT
      INTEGER              :: IML, ISG, IMLSW, NINDML, NSHIFT
      INTEGER              :: IQUADMAX(NMAXB), IQUADMIN(NMAXB)
      REAL                 :: CFLMEAN(NMAXB)
      REAL                 :: RTG, SIGFAC, RTDEP, KH, SIGND, CGND, CG
      INTEGER, ALLOCATABLE :: ISSETG(:), MAPML(:)
      INTEGER              :: VALID_TYPE, UNDEF_TYPE
      INTEGER              :: N_AVE, IQ_AVE(1), EX_TYPE(1)
      INTEGER              :: IUN, IERS
!
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE      
      IF ( PRESENT(IERR) ) IERR = 0
!
!  Constants:
      RTG = 3.1321        ! sqrt(g)
      SIGFAC = SIG1/RTG
!
!  Initialise scalars:
      IMODE = 0
      NREF = 0
      CGND = 0.
!
      VALID_TYPE = 1
      UNDEF_TYPE = 0
!
      N_AVE = 0
      IQ_AVE = 0
      EX_TYPE = UNDEF_TYPE
!
!  Initialise arrays:
      ZLEVT = 0.
      CFL = 0.
      IQUADMAX = 0
      IQUADMIN = 0
      IQREF = 0
      ISREF = 0
      ISEAOLDC = 0
      ISEANEWC = 0
      ISEAOLDR = 0
      ISEANEWR = 0
      CFLMEAN = 0.
!
      NX0 = NINT(NX*(2.**(-LVLREF)))
      NY0 = NINT(NY*(2.**(-LVLREF)))
      DELX0 = 2.**LVLREF
      DELY0 = DELX0
!
!  Allocate temporary quadtree arrays
!
      CALL QA_ALLOC ( QTREE_TEMP, NMAXB, NMAXB, undef_type=UNDEF_TYPE )
!
!  Allow cell indices to be retained when they are refined
!
      QTREE_TEMP%KEEP_REF = .TRUE.
      QTREE_TEMP%LVLMAX = LVLMAX
!
! Create a mapping ISSETG from full multilevel grid indices to bathymetry 
! multigrid indices, from INDBML which goes the other way
!
      !write(*,*) 'QA_BMLSTRUC: NX, NY, NMAXB, LVLREF: ',              &
      !                         NX, NY, NMAXB, LVLREF
      NINDML = NINT( NX0*NY0*(4.**(LVLMAX+1) - 1)/3. )
      !write(*,*) 'QA_BMLSTRUC: NX0,NY0,LVLMAX,NINDML: ',              &
      !                        NX0,NY0,LVLMAX,NINDML
      ALLOCATE ( ISSETG(NINDML), MAPML(NINDML) )
      ISSETG = 0
      MAPML = UNDEF_TYPE
      DO ISG=1,NMAXB
         IML = INDBML(ISG)
         ISSETG(IML) = ISG
         IF ( ISWET(ISG) ) MAPML(IML) = VALID_TYPE
         !write(*,*) ISG, IML, ZBML(ISG), MAPML(IML)
      END DO
!
! Loop over increasingly fine levels
      DO ILEV=0,LVLREF
         IF (ILEV.EQ.0) THEN
!
! Set up a quadtree structure for the level zero grid
            CALL QA_RECT ( NX0, NY0, LVLREF, GLOBAL, VALID_TYPE,      &
                        UNDEF_TYPE, QTREE_TEMP, ierr=IERS, ndse=IUN )
            IF ( IERS.GT.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BMLSTRUC'
                  WRITE (IUN,*) ' WHILE CALLING QA_RECT'
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
            QTREE_TEMP%LVLMAX = LVLMAX
         ELSE
!
! Refine all quads by one level
! Note: this call to QA_REFINE won't use SPARE approach, 
            CALL QA_REFINE ( IMODE, QTREE_TEMP, NREF, IQREF, ISREF,   &
                             ISEAOLDR, ISEANEWR, ierr=IERS, ndse=IUN, &
                             mapml=MAPML )
            IF ( IERS.GT.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BMLSTRUC'
                  WRITE (IUN,*) ' WHILE CALLING QA_REFINE (1ST PASS)'
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
!               
         END IF
         NSEA = QTREE_TEMP%NCELL
         NSEA_DEF = QTREE_TEMP%NCELL_DEF
         NQUAD = QTREE_TEMP%NQUAD
         !write(*,*) 'QA_BMLSTRUC (1) ILEV,NSEA,NQUAD: ',ILEV,NSEA,NQUAD
!
      END DO
!
! Loop over increasingly fine levels
      DO ILEV=LVLREF+1,LVLMAX+1
!
! Loop through cells, only keeping those with positive depth
         NSEA2 = 0
         NREF = 0
         DO ISEA=1,NSEA
            IF ( QTREE_TEMP%CELL_TYPE(ISEA).EQ.UNDEF_TYPE ) CYCLE
            IQ = QTREE_TEMP%INDQUAD(ISEA)
            ISUB = QTREE_TEMP%INDSUB(ISEA)
            IML = QTREE_TEMP%INDML(ISEA)
            IF ( IML.LT.1 .OR. IML.GT.NINDML ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BMLSTRUC'
                  WRITE (IUN,*) ' INDML OUT OF RANGE'
                  WRITE (IUN,*) ' NX0:', NX0
                  WRITE (IUN,*) ' NY0:', NY0
                  WRITE (IUN,*) ' LVLREF:', LVLREF
                  WRITE (IUN,*) ' LVLMAX:', LVLMAX
                  WRITE (IUN,*) ' ILEV:', ILEV
                  WRITE (IUN,*) ' ISEA:', ISEA
                  WRITE (IUN,*) ' IML:', IML
               END IF
               IF ( PRESENT(IERR) ) IERR = 1
               RETURN
            END IF
            ISG = ISSETG(IML)
            IF (ISG.GT.0) THEN
               ZT = ZBML(ISG)
            ELSE
               ZT = 1.
            END IF
            !MAPTMP = MAPML(IML)
            IF ( ILEV.EQ.LVLREF+1 ) THEN
               IX = NINT( QTREE_TEMP%XYVAL(ISEA,1) )
               IY = NINT( QTREE_TEMP%XYVAL(ISEA,2) )
               MAPTMP = MAPSTA(IY,IX)
            ELSE
               MAPTMP = VALID_TYPE
            END IF
            IF ( ZT.GE.0. .OR. MAPTMP.EQ.UNDEF_TYPE ) THEN
               NSEATEM = 0
            ELSE
               NSEA2 = NSEA2 + 1
               NSEATEM = NSEA2
               QTREE_TEMP%XYVAL(NSEA2,1) = QTREE_TEMP%XYVAL(ISEA,1)
               QTREE_TEMP%XYVAL(NSEA2,2) = QTREE_TEMP%XYVAL(ISEA,2)
               ZLEVT(NSEA2) = ZT
               IBMLT(NSEA2) = ISG
               QTREE_TEMP%INDQUAD(NSEA2) = IQ
               QTREE_TEMP%INDSUB(NSEA2) = ISUB
               QTREE_TEMP%INDML(NSEA2) = IML
               QTREE_TEMP%CELL_TYPE(NSEA2) = MAPTMP
               MAPML(IML) = MAPTMP
!
!              Look for refinement if we aren't at LVLMAX
               IF ( ILEV.LE.LVLMAX ) THEN
!                 Compute (nondimensional) wavenumber from the nondimensional
!                 frequency omega*sqrt(ZLEV/g): with a first estimate,
                  RTDEP = SQRT(ABS(ZT))
                  SIGND = SIGFAC*RTDEP
                  IF (SIGND .LT.1) THEN
                     KH = SIGND
                  ELSE
                     KH = SIGND*SIGND
                  END IF
!
!                 Then with iteration
                  CALL QA_DISPREL(SIGND, KH, CGND, ierr=IERS, ndse=IUN)
                  IF ( IERS.GT.0 ) THEN
                     IF ( IUN.GT.0 ) THEN
                        WRITE (IUN,*) ' ERROR IN QA_BMLSTRUC'
                        WRITE (IUN,*) ' WHILE CALLING QA_DISPREL'
                     END IF
                     IF ( PRESENT(IERR) ) IERR = IERS
                     RETURN
                  END IF
                  CG = CGND*RTDEP*RTG
                  CFL(NSEA2) = CG*2.**QTREE_TEMP%QLEVEL(IQ)
                  !
                  !  See if this cell is "refinable"
                  !IF ( QTREE_TEMP%CELL_TYPE(NSEA2) .NE. UNDEF_TYPE ) THEN
                  IF ( ISUB.GT.0 .OR.                                  &
                       ( ALL(QTREE_TEMP%QICELL(IQ,1:4).EQ.0) .AND.     &
                         ALL(QTREE_TEMP%QCHILD(IQ,1:4).EQ.0) ) ) THEN
                     IMLSW = ML_CHILD ( IML, NX0, NY0 )
                     !write(*,*) 'ISEA, ISUB, IML IMLSW: ',            &
                     !              ISEA, ISUB, IML, IMLSW
                     IF ( IMLSW.GT.0 .AND. IMLSW.LE.NINDML ) THEN
                        !write(*,*) 'ISSETG: ', ISSETG(IMLSW)
                        IF ( ISSETG(IMLSW).GT.0 ) THEN
                           NREF = NREF + 1
                           ISEAOLDR(NREF) = NSEA2
                        END IF
                     END IF
                  END IF
               END IF
            END IF
            QTREE_TEMP%QICELL(IQ,ISUB) = NSEATEM
         END DO  ! Loop over existing quadtree cells
         NSEA = NSEA2
         !write(*,*) 'QA_BMLSTRUC (2A) ILEV,NSEA,NQUAD, NREF: ',ILEV,NSEA,NQUAD,NREF
         !write(*,*) 'QA_BMLSTRUC ISEAOLDR:'
         !do isea=1,NREF
         !   write(*,*) ISEAOLDR(isea)
         !end do
         IF ( ILEV.LE.LVLMAX .AND. NREF.GT.0 ) THEN
!
! Refine all cells that can be refined by one level
            IMODE = 2
! Note: this call to QA_REFINE won't use SPARE approach, 
            CALL QA_REFINE ( IMODE, QTREE_TEMP, NREF, IQREF, ISREF,   &
                             ISEAOLDR, ISEANEWR, ierr=IERS, ndse=IUN, &
                             mapml=MAPML )
            IF ( IERS.GT.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BMLSTRUC'
                  WRITE (IUN,*) ' WHILE CALLING QA_REFINE (2ND PASS)'
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
            !write(*,*) 'QA_BMLSTRUC : QA_REFINE CALLED'
         END IF
         NSEA = QTREE_TEMP%NCELL
         NSEA_DEF = QTREE_TEMP%NCELL_DEF
         NQUAD = QTREE_TEMP%NQUAD
         !write(*,*) 'QA_BMLSTRUC (2B) ILEV,NSEA,NQUAD: ',ILEV,NSEA,NQUAD
!
      END DO
!
! Now we have a quadtree structure fully refined to the finest level 
! allowed by bathymetry data. If MSEA=0, return with the corresponding
! NSEA and NQUAD values, without computing arrays
!
      IF ( PRESENT(QTREE) ) THEN
         MSEA = SIZE(QTREE%INDQUAD,1)
         MQUAD = SIZE(QTREE%QLEVEL,1)
      ELSE
         MSEA = 0
         MQUAD = 0
      END IF
      IF ( MSEA.LE.0 .OR. MQUAD.LE.0 ) RETURN
!
! Otherwise, we may need to coarsen to make NSEA <= MSEA.
! Loop through quads, selecting those with the largest "CFL factor"
! (and containing 4 wet cells) for coarsening
!
      DO ILEV = LVLMAX,0,-1
         IF (ILEV.GT.0) THEN
            NCRS = NINT( (NSEA - MSEA)/(3.*ILEV) )
            NCRS = MAX( NCRS, NQUAD-MQUAD )
         ELSE
            NCRS = NSEA - MSEA
         END IF
         IF (NCRS.LE.0) EXIT
         CALL QA_QMEAN ( QTREE_TEMP, array=CFL, amean=CFLMEAN,         &
                         n_ave=N_AVE, iq_ave=IQ_AVE, exc_mean=EX_TYPE, &
                         nminmax=NCRS, nwetmin=4, iquadmax=IQUADMAX,   &
                         iquadmin=IQUADMIN, exc_minmax=EX_TYPE )
         IMODE = 1
! Note: this call to QA_COARSEN won't use SPARE approach
         CALL QA_COARSEN ( IMODE, QTREE_TEMP, NCRS, IQUADMAX,         &
                           ISEAOLDC, ISEANEWC, NSHIFT, mapml=MAPML )
         NSEA = QTREE_TEMP%NCELL
         NQUAD = QTREE_TEMP%NQUAD
         !write(*,*) 'QA_BMLSTRUC (3A) ILEV,NSEA,NQUAD: ',ILEV,NSEA,NQUAD
!
! Recompute depths at changed cells
         DO ISEA=1,NSEA
            IF ( QTREE_TEMP%CELL_TYPE(ISEA).EQ.UNDEF_TYPE ) CYCLE
            IQ = QTREE_TEMP%INDQUAD(ISEA)
            IML = QTREE_TEMP%INDML(ISEA)
            ISG = ISSETG(IML)
            IF (ISG.GT.0) THEN
               ZT = ZBML(ISG)
            ELSE
               ZT = 1.
            END IF
!
            ZLEVT(ISEA) = ZT
            IBMLT(ISEA) = ISG
!           Compute (nondimensional) wavenumber from the nondimensional
!           frequency omega*sqrt(depth/g): with a first estimate,
            RTDEP = SQRT(ABS(ZT))
            SIGND = SIGFAC*RTDEP
            IF (SIGND .LT.1) THEN
               KH = SIGND
            ELSE
               KH = SIGND*SIGND
            END IF
!
!           Then with iteration
            CALL QA_DISPREL(SIGND,KH,CGND)
            CG = CGND*RTDEP*RTG
            CFL(ISEA) = CG*2.**QTREE_TEMP%QLEVEL(IQ)
            QTREE_TEMP%CELL_TYPE(ISEA) = VALID_TYPE
         END DO
      END DO
!
! Recompute cell neighbours
      CALL QA_FINDNBR ( QTREE_TEMP, ierr=IERS, ndse=IUN )
      IF ( IERS.GT.0 ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE (IUN,*) ' ERROR IN QA_BMLSTRUC'
            WRITE (IUN,*) ' WHILE CALLING QA_FINDNBR'
         END IF
         IF ( PRESENT(IERR) ) IERR = IERS
         RETURN
      END IF
      
      IF( PRESENT(QTREE) ) THEN
         MSEA = SIZE(QTREE%INDQUAD,1)
         MQUAD = SIZE(QTREE%QPARENT,1)
         !
         MSEA = MIN(NMAXB, MSEA)
         MQUAD = MIN(NMAXB,MQUAD)
         !
         QTREE%NCELL = NSEA
         QTREE%NCELL_DEF = NSEA_DEF
         QTREE%NQUAD = NQUAD
         QTREE%NX0 = NX0
         QTREE%NY0 = NY0
         QTREE%LVLREF = LVLREF
         QTREE%LVLMAX = LVLMAX
         QTREE%LVLHI = LVLMAX
         QTREE%UNDEF_TYPE = UNDEF_TYPE
         QTREE%KEEP_REF = .TRUE.
         QTREE%DYNAMIC = QTREE_TEMP%DYNAMIC
!
! Transfer seapoint arrays to optional output versions
         QTREE%INDQUAD(1:MSEA) = QTREE_TEMP%INDQUAD(1:MSEA)
         QTREE%INDSUB(1:MSEA) = QTREE_TEMP%INDSUB(1:MSEA)
         QTREE%XYVAL(1:MSEA,:) = QTREE_TEMP%XYVAL(1:MSEA,:)
         QTREE%NGBR(1:MSEA,:) = QTREE_TEMP%NGBR(1:MSEA,:)
         QTREE%INDLVL(1:MSEA) = QTREE_TEMP%INDLVL(1:MSEA)
         QTREE%CELL_TYPE(1:MSEA) = QTREE_TEMP%CELL_TYPE(1:MSEA)
!
! Transfer quad arrays to optional output versions
         QTREE%QPARENT(1:MQUAD) = QTREE_TEMP%QPARENT(1:MQUAD)
         QTREE%QLEVEL(1:MQUAD) = QTREE_TEMP%QLEVEL(1:MQUAD)
         QTREE%QICELL(1:MQUAD,:) = QTREE_TEMP%QICELL(1:MQUAD,:)
         QTREE%QNBR(1:MQUAD,:) = QTREE_TEMP%QNBR(1:MQUAD,:)
         QTREE%QCHILD(1:MQUAD,:) = QTREE_TEMP%QCHILD(1:MQUAD,:)
      END IF
      IF( PRESENT(IBML) ) THEN
         MSEA = SIZE(IBML,1)
         MSEA = MIN(NMAXB, MSEA)
         IBML(1:MSEA) = IBMLT(1:MSEA)
      END IF
      IF( PRESENT(ZLEV) ) THEN
         MSEA = SIZE(ZLEV,1)
         MSEA = MIN(NMAXB, MSEA)
         ZLEV(1:MSEA) = ZLEVT(1:MSEA)
      END IF
      RETURN
      END SUBROUTINE QA_BMLSTRUC
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_BRSTRUC ( NSEAT, NX, NY, ZBIN, MAPSTA, SIG1,      &
                              GLOBAL, NSEA, NQUAD, IERR, NDSE, QTREE, &
                              DEP )
!/
!       Richard Gorman, NIWA
!         May, 2008:        Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!      Create a threaded quadtree structure based on a rectangular grid,
!      with a limited number of cells produced by coarsening where 
!      a "CFL factor" = CG/cell width is greatest
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ---------------------------------------------------------------
!       NSEAT   Int.   I   Target number of cells
!       NX, NY  Int.   I   Size of input grid , at the finest level
!       ZBIN    R.A.   I   Elevation of the seabed on fine grid 
!                          (m above datum)
!       MAPSTA  I.A.   I   Flag for sea/boundary/excluded points 
!                          (on rectangular grid)
!       SIG1    Real   I   Radian frequency of lowest bin
!       GLOBAL  Log.   I   Flag for global grid (i.e. cyclic in X)
!       NSEA    Int.   O   Actual number of cells
!       NQUAD   Int.   O   Actual number of quads
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!       QTREE   QA_TREE O* Quadtree structure, of which the following 
!                          components are affected:
!         NCELL   Int.       Number of cells (including undefined cells)
!         NCELL_DEF Int.     Number of active cells
!         LVLMAX  Int.       Maximum refinement level allowed
!         LVLHI   Int.       Maximum refinement level actually present
!         NQUAD   Int.       Number of quads
!         INDLVL  I.A.       Level of each cell
!         INDML   I.A.     
!         NGBR    I.A.       Neighbours (up to 8) of each cell
!         XYVAL   R.A.       Coordinates of cells (w.r.t.  
!                            reference grid)
!         CELL_TYPE  I.A.    Flag for sea/boundary/excluded points
!         INDQUAD I.A.       Index of the quad containing each input 
!                            cell
!         INDSUB  I.A.       Index of the sub-quad containing each sea
!                            point:
!                              0 = none (home cell is at level 0) 
!                              1,2,3,4 = SW, SE, NW, NE corner of quad
!         QICELL  I.A.       Sea-point indices for the cells of the quad
!         QLEVEL  I.A.       Level of the cells within each quad
!         QPARENT I.A.       Quad index for the parent of each quad
!         QNBR    I.A.       Quad indices for the neighbours of each quad
!         QCHILD  I.A.       Quad indices for the children of each quad
!       DEP     R.A.   O*  Elevation of the seabed at each sea
!                          point (m above datum)
!     ---------------------------------------------------------------
!                       * = OPTIONAL
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     QA_ALLOC   Subr. qa_utils Allocate a quadtree structure
!     QA_RECT    Subr. qa_utils Set up a quadtree representation (level 0)
!     QA_REFINE  Subr. qa_utils Refine selected cells of a quadtree grid
!     QA_COARSEN Subr. qa_utils Coarsen selected quads of a quadtree grid
!     QA_QMEAN   Subr. qa_utils Compute the mean of a variable over 
!                               the cells of each quad
!     QA_FINDNBR Subr. qa_utils Compute the seapoint neighbour indices
!     QA_DISPREL Subr. qa_utils Compute CG using the dispersion relation
!     ----------------------------------------------------------------
!
!  5. Called by :
!
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
      INTEGER, INTENT(IN)   :: NSEAT, NX, NY
      REAL, INTENT(IN)      :: ZBIN(NX,NY), SIG1
      INTEGER, INTENT(IN)   :: MAPSTA(NY,NX)
      LOGICAL, INTENT(IN)   :: GLOBAL
      INTEGER, INTENT(OUT)  :: NSEA
      INTEGER, INTENT(OUT)  :: NQUAD
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
      TYPE(QA_TREE), INTENT(INOUT), OPTIONAL  :: QTREE
      REAL,    INTENT(OUT), OPTIONAL  :: DEP(:)     
!
! Local variables.
!
!  Temporary versions of the output quadtree arrays, declared for
!  full number of seapoints on the finest grid
      TYPE(QA_TREE)        :: QTREE_TEMP
!
      INTEGER              :: MSEA, MQUAD
      INTEGER              :: NSEA2, NQUAD2, NSEA_DEF
      REAL                 :: CFL(NX*NY), DEPT(NX*NY)
      INTEGER              :: ISREF(NX*NY), ISEAOLDR(NX*NY),          &
                              ISEANEWR(NX*NY,4)
      INTEGER              :: ISEAOLDC(NX*NY,4), ISEANEWC(NX*NY)
      INTEGER              :: IQREF(NX*NY)
      INTEGER              :: NX0, NY0, MAXLEVEL, LVLREF, ILEV
      INTEGER              :: IQ, ISEA, IMODE, NREF, ISUB,            &
                              NCRS, IX, IY, NSEATEM, NSHIFT
      REAL                 :: DELX0, DELY0, DPT, RX, RY, WTSUM, WT
      INTEGER              :: IXP, IYP, NSUM
      INTEGER              :: IQUADMAX(NX*NY), IQUADMIN(NX*NY)
      REAL                 :: CFLMEAN(NX*NY)
      REAL                 :: RTG, SIGFAC, RTDEP, KH, SIGND, CGND, CG
      INTEGER              :: VALID_TYPE, UNDEF_TYPE
      INTEGER              :: N_AVE, IQ_AVE(1), EX_TYPE(1)
      INTEGER              :: IUN, IERS
!
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE      
      IF ( PRESENT(IERR) ) IERR = 0
!
!  Constants:
      RTG = 3.1321
      SIGFAC = SIG1/RTG
!
!  Initialise scalars:
      NSEA = 0
      NQUAD2 = 0
      IMODE = 0
      NREF = 0
      CGND = 0.
!
      VALID_TYPE = 1
      UNDEF_TYPE = 0
!
      N_AVE = 0
      IQ_AVE = 0
      EX_TYPE = UNDEF_TYPE
!
!  Initialise arrays:
!
      DEPT = 0.
      CFL = 0.
      IQUADMAX = 0
      IQUADMIN = 0
      IQREF = 0
      ISREF = 0
      ISEAOLDC = 0
      ISEANEWC = 0
      ISEAOLDR = 0
      ISEANEWR = 0
      CFLMEAN = 0.
!
!  Allocate temporary quadtree arrays
!
      MSEA = NX*NY
      CALL QA_ALLOC ( QTREE_TEMP, MSEA, MSEA, undef_type=UNDEF_TYPE )
!
!      NSPARE = -1
!      ISSPARE = 0
!      PARTNO = NPART
!
!  Compute how many levels are allowed
      NX0 = NX
      NY0 = NY
      MAXLEVEL = 0
      DO WHILE( MOD(NX0,2).EQ.0 .AND. MOD(NY0,2).EQ.0 )
         NX0 = NX0/2
         NY0 = NY0/2
         MAXLEVEL = MAXLEVEL + 1
      END DO
      LVLREF = MAXLEVEL
      DELX0 = 2.**MAXLEVEL
      DELY0 = 2.**MAXLEVEL
!
! Loop over increasingly fine levels
      DO ILEV=0,MAXLEVEL
         IF (ILEV.EQ.0) THEN
!
! Set up a quadtree structure for the level zero grid
            CALL QA_RECT ( NX0, NY0, MAXLEVEL, GLOBAL, VALID_TYPE,    &
                       UNDEF_TYPE, QTREE_TEMP, ierr=IERS, ndse=IUN )
            IF ( IERS.GT.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BRSTRUC'
                  WRITE (IUN,*) ' WHILE CALLING QA_RECT WITH:'
                  WRITE (IUN,*) '          NX0 = ', NX0
                  WRITE (IUN,*) '          NY0 = ', NY0
                  WRITE (IUN,*) '     MAXLEVEL = ', MAXLEVEL
                  WRITE (IUN,*) '       GLOBAL = ', GLOBAL
                  WRITE (IUN,*) '   VALID_TYPE = ', VALID_TYPE
                  WRITE (IUN,*) '   UNDEF_TYPE = ', UNDEF_TYPE
                  WRITE (IUN,*) '  QTREE_ALLOC = ', MSEA, MSEA
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
         ELSE
!
! Note: this call to QA_REFINE won't use SPARE approach, 
            CALL QA_REFINE ( IMODE, QTREE_TEMP, NREF, IQREF, ISREF,    &
                             ISEAOLDR, ISEANEWR, ierr=IERS, ndse=IUN )
            IF ( IERS.NE.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BRSTRUC'
                  WRITE (IUN,*) ' WHILE CALLING QA_REFINE'
                  WRITE (IUN,*) '    ILEV = ', ILEV
                  WRITE (IUN,*) '    IMODE = ', IMODE
                  WRITE (IUN,*) '    NREF = ', NREF
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
!               
            NSEA = QTREE_TEMP%NCELL
            NSEA_DEF = QTREE_TEMP%NCELL_DEF
            NQUAD = QTREE_TEMP%NQUAD
         END IF
!
      END DO
!
! Loop through cells, only keeping those with positive depth
      NSEA2 = 0
      DO ISEA=1,NSEA
         IQ = QTREE_TEMP%INDQUAD(ISEA)
         ISUB = QTREE_TEMP%INDSUB(ISEA)
         IX = NINT( QTREE_TEMP%XYVAL(ISEA,1) )
         IY = NINT( QTREE_TEMP%XYVAL(ISEA,2) )
         IF ( IX.LT.1 .OR. IX.GT.NX .OR. IY.LT.1 .OR. IY.GT.NY ) THEN
            DPT = 1.
            IX = 1
            IY = 1            
         ELSE
            DPT = ZBIN(IX,IY)
         END IF
         IF (DPT.GE.0. .OR. MAPSTA(IY,IX).EQ.0) THEN
            NSEATEM = 0
         ELSE
            NSEA2 = NSEA2 + 1
            NSEATEM = NSEA2
            QTREE_TEMP%XYVAL(NSEA2,1) = QTREE_TEMP%XYVAL(ISEA,1)
            QTREE_TEMP%XYVAL(NSEA2,2) = QTREE_TEMP%XYVAL(ISEA,2)
            DEPT(NSEA2) = DPT
            QTREE_TEMP%INDQUAD(NSEA2) = IQ
            QTREE_TEMP%INDSUB(NSEA2) = ISUB
            QTREE_TEMP%CELL_TYPE(NSEA2) = MAPSTA(IY,IX)
!           Compute (nondimensional) wavenumber from the nondimensional
!           frequency omega*sqrt(depth/g): with a first estimate,
            RTDEP = SQRT(ABS(DPT))
            SIGND = SIGFAC*RTDEP
            IF (SIGND .LT.1) THEN
               KH = SIGND
            ELSE
               KH = SIGND*SIGND
            END IF
!
!           Then with iteration
            CALL QA_DISPREL( SIGND, KH, CGND, ierr=IERS, ndse=IUN )
            IF ( IERS.GT.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE (IUN,*) ' ERROR IN QA_BRSTRUC'
                  WRITE (IUN,*) ' WHILE CALLING QA_DISPREL'
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
            CG = CGND*RTDEP*RTG
            CFL(NSEA2) = CG*2.**QTREE_TEMP%QLEVEL(IQ)
         END IF
         QTREE_TEMP%QICELL(IQ,ISUB) = NSEATEM
      END DO
      NSEA = NSEA2
!
! Search for cells to coarsen, by averaging the CFL factor over quads,
! and selecting quads with the smallest quad-averaged values
!
! Loop through quads, selecting those with no wet cells for coarsening
      DO ILEV = MAXLEVEL,0,-1
         IF (ILEV.GT.0) THEN
            NCRS = NINT( (NSEA - NSEAT)/(3.*ILEV) )
         ELSE
            NCRS = NSEA - NSEAT
         END IF
         IF (NCRS.LE.0) EXIT
         CALL QA_QMEAN ( QTREE_TEMP, array=CFL, amean=CFLMEAN,        &
                         n_ave=N_AVE, iq_ave=IQ_AVE,                  &
                         iquadmax=IQUADMAX, iquadmin=IQUADMIN,        &
                         exc_minmax=EX_TYPE )
         IMODE = 1
! Note: this call to QA_COARSEN won't use SPARE approach
         CALL QA_COARSEN ( IMODE, QTREE_TEMP, NCRS, IQUADMAX,         &
                           ISEAOLDC, ISEANEWC, NSHIFT )
         NSEA = QTREE_TEMP%NCELL
         NQUAD = QTREE_TEMP%NQUAD
!
! Recompute depths at changed  cells
         DO ISEA=1,NSEA
            IQ = QTREE_TEMP%INDQUAD(ISEA)
            IF (QTREE_TEMP%QLEVEL(IQ) .EQ.MAXLEVEL) THEN
               IX = NINT( QTREE_TEMP%XYVAL(ISEA,1) )
               IY = NINT( QTREE_TEMP%XYVAL(ISEA,2) )
               DPT = ZBIN(IX,IY)
            ELSE
               IX = INT( QTREE_TEMP%XYVAL(ISEA,1) )
               RX = QTREE_TEMP%XYVAL(ISEA,1) - IX
               IY = INT( QTREE_TEMP%XYVAL(ISEA,2) )
               RY = QTREE_TEMP%XYVAL(ISEA,2) - IY
               DPT = 0.
               WTSUM = 0.
               NSUM = 0
               IX = MAX(IX,1)
               IXP = MIN(IX+1,NX)
               IY = MAX(IY,1)
               IYP = MIN(IY+1,NY)
               IF (ZBIN(IX,IY).LT.0.) THEN
                  NSUM = NSUM+1
                  WT = (1.-RX)*(1.-RY)
                  WTSUM = WTSUM + WT
                  DPT = DPT + WT*ZBIN(IX,IY)
               END IF
               IF (ZBIN(IXP,IY).LT.0.) THEN
                  NSUM = NSUM+1
                  WT = RX*(1.-RY)
                  WTSUM = WTSUM + WT
                  DPT = DPT + WT*ZBIN(IXP,IY)
               END IF
               IF (ZBIN(IX,IYP).LT.0.) THEN
                  NSUM = NSUM+1
                  WT = (1.-RX)*RY
                  WTSUM = WTSUM + WT
                  DPT = DPT + WT*ZBIN(IX,IYP)
               END IF
               IF (ZBIN(IXP,IYP).LT.0.) THEN
                  NSUM = NSUM+1
                  WT = RX*RY
                  WTSUM = WTSUM + WT
                  DPT = DPT + WT*ZBIN(IXP,IYP)
               END IF
               IF (NSUM.GT.0) THEN
                  DPT = DPT/WTSUM
               END IF
            END IF
!
            DEPT(ISEA) = DPT
!           Compute (nondimensional) wavenumber from the nondimensional
!           frequency omega*sqrt(depth/g): with a first estimate,
            RTDEP = SQRT(ABS(DPT))
            SIGND = SIGFAC*RTDEP
            IF (SIGND .LT.1) THEN
               KH = SIGND
            ELSE
               KH = SIGND*SIGND
            END IF
!
!           Then with iteration
            CALL QA_DISPREL(SIGND,KH,CGND)
            CG = CGND*RTDEP*RTG
            CFL(ISEA) = CG*2.**QTREE_TEMP%QLEVEL(IQ)
            IX = NINT( QTREE_TEMP%XYVAL(ISEA,1) )
            IY = NINT( QTREE_TEMP%XYVAL(ISEA,2) )
            QTREE_TEMP%CELL_TYPE(ISEA) = MAPSTA(IY,IX)
         END DO
      END DO
!
! Recompute cell neighbours
      CALL QA_FINDNBR ( QTREE_TEMP, ierr=IERS, ndse=IUN )
      IF ( IERS.GT.0 ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE (IUN,*) ' ERROR IN QA_BRSTRUC'
            WRITE (IUN,*) ' WHILE CALLING QA_FINDNBR'
         END IF
         IF ( PRESENT(IERR) ) IERR = IERS
         RETURN
      END IF
!
! Output quadtree:
!
      IF( PRESENT(QTREE) ) THEN
         MSEA = SIZE(QTREE%INDQUAD,1)
         MQUAD = SIZE(QTREE%QPARENT,1)
         !
         MSEA = MIN(NX*NY, MSEA)
         MQUAD = MIN(NX*NY,MQUAD)
         !
         QTREE%NCELL = NSEA
         QTREE%NCELL_DEF = NSEA_DEF
         QTREE%NQUAD = NQUAD
         QTREE%NX0 = NX0
         QTREE%LVLREF = LVLREF
         QTREE%LVLMAX = MAXLEVEL
         QTREE%LVLHI = MAXLEVEL
         QTREE%UNDEF_TYPE = UNDEF_TYPE
         QTREE%KEEP_REF= .TRUE.
         QTREE%DYNAMIC = QTREE_TEMP%DYNAMIC
!
! Transfer seapoint arrays to optional output versions
         QTREE%INDQUAD(1:MSEA) = QTREE_TEMP%INDQUAD(1:MSEA)
         QTREE%INDSUB(1:MSEA) = QTREE_TEMP%INDSUB(1:MSEA)
         QTREE%XYVAL(1:MSEA,:) = QTREE_TEMP%XYVAL(1:MSEA,:)
         QTREE%NGBR(1:MSEA,:) = QTREE_TEMP%NGBR(1:MSEA,:)
         QTREE%INDLVL(1:MSEA) = QTREE_TEMP%INDLVL(1:MSEA)
         QTREE%CELL_TYPE(1:MSEA) = QTREE_TEMP%CELL_TYPE(1:MSEA)
!
! Transfer quad arrays to optional output versions
         QTREE%QPARENT(1:MQUAD) = QTREE_TEMP%QPARENT(1:MQUAD)
         QTREE%QLEVEL(1:MQUAD) = QTREE_TEMP%QLEVEL(1:MQUAD)
         QTREE%QICELL(1:MQUAD,:) = QTREE_TEMP%QICELL(1:MQUAD,:)
         QTREE%QNBR(1:MQUAD,:) = QTREE_TEMP%QNBR(1:MQUAD,:)
         QTREE%QCHILD(1:MQUAD,:) = QTREE_TEMP%QCHILD(1:MQUAD,:)
      END IF
      IF( PRESENT(DEP) ) THEN
         MSEA = SIZE(DEP,1)
         MSEA = MIN(NX*NY, MSEA)
         DEP(1:MSEA) = DEPT(1:MSEA)
      END IF
!      
      RETURN
      END SUBROUTINE QA_BRSTRUC
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_COARSEN ( IMODE, QTREE, NCRS, IQCRS, ISEAOLD,     &
                              ISEANEW, NSHIFT, MAPML, QSPARE )
!/
!/
!       Richard Gorman, NIWA
!         April, 2008:      Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!
!/
!  1. Purpose :
!
!      Coarsen selected cells of a threaded quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMODE   Int.   I   Input option: =0 to do all quads
!                                        =1 to use IQCRS
!                                        =2 to use ISEAOLD
!       QTREE   QA_TREE I/O  Quadtree structure, with the following components
!                             affected by this subroutine:
!          NCELL   Int.  I/O  Number of cells
!          NCELL_DEF Int. I/O Number of active cells
!          NQUAD   Int.  I/O  Number of quads
!          LVLREF  Int.   I   Level of reference grid
!          LVLMAX  Int.   I   Maximum refinement level allowed
!          LVLHI   Int.   O   Maximum refinement level present
!          NX0     Int.   I   X dimension of level-zero base grid
!          NY0     Int.   I   Y dimension of level-zero base grid
!          UNDEF_TYPE Int. I  value of flag for unused cell indices
!          INDQUAD I.A.  I/O  Index of the quad containing each input cell
!          INDSUB  I.A.  I/O  Index of the sub-quad containing each cell:
!                               0 = none (home cell is at level zero)
!                               1,2,3,4 = SW, SE, NW, NE corner of quad
!          INDLVL  I.A.  I/O  Level of each cell
!          CELL_TYPE I.A. Flag for wet/dry/boundary/etc. cell
!          INDML   I.A.  I/O  Multilevel index of cells
!          XYVAL   R.A.  I/O  Coordinates of cells (w.r.t. 
!                                reference rectangular grid)
!          QICELL  I.A.  I/O  Sea-point indices for the cells of the quad
!          QLEVEL  I.A.  I/O  Level of the cells within each quad
!          QPARENT I.A.  I/O  Quad index for the parent of each quad
!          QNBR    I.A.  I/O  Quad indices for the neighbours of each quad
!          QCHILD  I.A.  I/O  Quad indices for the children of each quad
!       NCRS    Int.  I/O  Number of cells to coarsen
!       IQCRS   I.A.  I/O  Quad numbers of cells to coarsen
!       ISEAOLD I.A.  I/O  Old sea indices of cells that have been replaced
!       ISEANEW I.A.   O   New sea indices of cells that have been coarsened
!       NSHIFT  Int.   O   Number of cells renumbered (without coarsening)
!       MAPML   I.A.   I*  Array of cell flags on multilevel grid 
!       QSPARE  QA_SPARE I/O* Arrays for treating spare cells in quadtree 
!                             structure, with the following components used:
!          NPART   Int.   I   Number of partitions
!          PARTNO  I.A.   I   Partition number for each cell
!          JSVAL   I.A.   I   Local cell index for each global index
!          NSPARE  I.A.  I/O  Number of spare local cell indices
!          JSSPARE I.A.  I/O  List of spare local cell indices
!     ----------------------------------------------------------------
!                        * = optional
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     ML_PARENT  Func. qa_utils Multilevel index of the parent of a given cell
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_BMLSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC   Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMODE
      TYPE(QA_TREE), INTENT(INOUT)  :: QTREE
      INTEGER, INTENT(INOUT)  :: NCRS
      INTEGER, INTENT(INOUT)  :: IQCRS(:), ISEAOLD(:,:)
      INTEGER, INTENT(OUT)    :: ISEANEW(:)
      INTEGER, INTENT(OUT)    :: NSHIFT
      INTEGER, OPTIONAL, INTENT(IN)     :: MAPML(:)
      TYPE(QA_SPARE), OPTIONAL, INTENT(INOUT)  :: QSPARE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER  :: MQUAD
      INTEGER  :: LVLHI
      INTEGER  :: NPART
!
      INTEGER           :: NSEA1, NQUAD1
      INTEGER           :: NSEA2, NQUAD2, NSEA_DEF
      REAL              :: DELX0, DELY0
      INTEGER         :: ISEA, QNC, ICRS, IQ, ISUB, QISEA0, II,       &
                         IQPAR, NDELQ, QISEA1, NLOOP, IQNEW, IQC,     &
                         QSHIFT, IDELQ, NEXTDEL, ISTMP, PV, NSP,      &
                         JSEA, NINDML, IML
      INTEGER, ALLOCATABLE   :: IQDEL(:), NEWQUAD(:)
      REAL, ALLOCATABLE      :: QDELS(:)
      INTEGER, ALLOCATABLE   :: INDSORT(:)
      REAL            :: DELX, DELY
      LOGICAL         :: DOSPARE, DOPART, NOTCRS
!/
!/ ------------------------------------------------------------------- /
!/
      NSEA1 = QTREE%NCELL
      NQUAD1 = QTREE%NQUAD
!
      DOSPARE = PRESENT(QSPARE) 
      IF ( DOSPARE ) THEN
         NPART = QSPARE%NPART
      ELSE
         NPART = 1
      END IF
      DOPART = NPART.GT.1
!
      MQUAD = SIZE(QTREE%QICELL,1)
      DELX0 = 2.**QTREE%LVLREF
      DELY0 = DELX0
      !
      IF ( PRESENT(MAPML) ) THEN
         NINDML = SIZE(MAPML)
      ELSE
         NINDML = 0
      END IF
!
      ISEANEW = 0
!
      NDELQ = 0
      IF (IMODE.EQ.0) THEN
         NLOOP = NQUAD1
      ELSE
         NLOOP = NCRS
      END IF
!
      ALLOCATE ( IQDEL(NLOOP), NEWQUAD(MQUAD) )
      IQDEL = 0
      NEWQUAD = 0
!
      NCRS = 0
      DO ICRS=1, NLOOP
!
         IF (IMODE.EQ.0) THEN
           IQ = ICRS
         ELSEIF (IMODE.EQ.1) THEN
           IQ = IQCRS(ICRS)
         ELSE
           IQ = QTREE%INDQUAD(ISEAOLD(ICRS,1))
         END IF
         IF (IQ.LE.0) CYCLE
         IQPAR = QTREE%QPARENT(IQ)
!  If this is a top-level quad which is a single cell, 
!  don't try to coarsen:
         IF (IQPAR.EQ.0 .AND. QTREE%QICELL(IQ,0).GT.0) CYCLE
!         
         NCRS = NCRS + 1
         IQCRS(NCRS) = IQ
         ! If cell indices are not retained for quad centres:
         !QISEA0 = 0
         ! If cell indices MAY be retained for quad centres:
         QISEA0 = QTREE%QICELL(IQ,0)
         DELX = -DELX0*2.**(-QTREE%QLEVEL(IQ)-1)
         DELY = DELY0*2.**(-QTREE%QLEVEL(IQ)-1)
         IF (IQPAR.NE.0) THEN
!  This is not a top-level quad: flag it for removal
            NDELQ = NDELQ + 1  
            IQDEL(NDELQ) = IQ  
         END IF
         DO II=1,4
            DELX = -DELX
            IF (II.EQ.3) DELY = -DELY
            QISEA1 = QTREE%QICELL(IQ,II)
            ISEAOLD(NCRS,II) = QISEA1
            IF(QISEA1.GT.0) THEN
               IF(QISEA0.EQ.0) THEN
                  QISEA0 = QISEA1
                  QTREE%XYVAL(QISEA0,1) = QTREE%XYVAL(QISEA0,1) + DELX
                  QTREE%XYVAL(QISEA0,2) = QTREE%XYVAL(QISEA0,2) + DELY
                  IML = ML_PARENT(QTREE%INDML(QISEA1),QTREE%NX0,QTREE%NY0)
                  QTREE%INDML(QISEA0) = IML
                  ! If available, use the MAPML array to give the cell type
                  ! of the new, coarsened cell: otherwise it will be inherited
                  ! from the child who's index it took over.
                  IF ( NINDML.GE.IML .AND. PRESENT(MAPML) ) THEN
                     QTREE%CELL_TYPE(QISEA0) = MAPML(IML)
                  END IF
                  IF (IQPAR.EQ.0) THEN
!  This is a top-level quad: don't remove it, but depopulate it.
!  Find the cell index of the first cell that is a cell
!  Set all children and cell indices to nil:
                     QTREE%QICELL(IQ,0) = QISEA0
                     ISEANEW(NCRS) = QISEA0
                     QTREE%INDQUAD(QISEA0) = IQ
                     QTREE%INDSUB(QISEA0) = 0
                     QTREE%INDLVL(QISEA0) = 0                 
                  END IF
               ELSE
                  QTREE%INDQUAD(QISEA1) = 0
                  QTREE%INDSUB(QISEA1) = 0
                  QTREE%INDLVL(QISEA1) = 0
                  QTREE%INDML(QISEA1) = 0
                  QTREE%CELL_TYPE(QISEA1) = QTREE%UNDEF_TYPE
               END IF
            END IF
            QTREE%QICELL(IQ,II) = 0
            QTREE%QCHILD(IQ,II) = 0   ! WAS COMMENTED OUT
         END DO
         IF (IQPAR.NE.0) THEN
!  This is not a top-level quad: flag it for removal
            ISEANEW(NCRS) = QISEA0
!  Transfer this cell index to the corresponding cell of the parent
!  Remove this quad as a child from the parent quad:
            DO II=1,4
               QNC = QTREE%QCHILD(IQPAR,II)
               IF (QNC.EQ.IQ) THEN
                  QTREE%QICELL(IQPAR,II) = QISEA0
                  QTREE%QCHILD(IQPAR,II) = 0   ! WAS COMMENTED OUT
                  IF (QISEA0.GT.0) THEN
                     QTREE%INDQUAD(QISEA0) = IQPAR
                     QTREE%INDSUB(QISEA0) = II
                     QTREE%INDLVL(QISEA0) = QTREE%QLEVEL(IQPAR)
                  END IF
                  EXIT
               END IF
            END DO
         END IF
      END DO
! 
!   Find mapping from old to new quad indices.
!  
      QSHIFT = 0
      IF ( NDELQ.GT.0 ) THEN
         ALLOCATE (QDELS(NDELQ))
         ALLOCATE (INDSORT(NDELQ))
         QDELS = FLOAT(IQDEL(1:NDELQ))
         DO II=1,NDELQ
            INDSORT(II) = II
         END DO
         CALL QSORT(QDELS,INDSORT,NDELQ)
         IQDEL(1:NDELQ) = INT(QDELS)
         IDELQ = 1
         NEXTDEL = IQDEL(1)
      ELSE
         IDELQ = 0
         NEXTDEL = 0
      END IF
      DO IQ=1,NQUAD1
         IF ( IQ.EQ.NEXTDEL ) THEN
            !write(*,*) 'QA_COARSEN: deleting IQ = ',IQ
            NEWQUAD(IQ) = 0
            QSHIFT = QSHIFT + 1
            IF ( IDELQ.LT.NDELQ ) THEN
               IDELQ = IDELQ + 1
               NEXTDEL = IQDEL(IDELQ)
            END IF
         ELSE
            NEWQUAD(IQ) = IQ - QSHIFT
         END IF
      END DO
! 
!  Loop through removing cells, and remapping to new cell
!  indices.
      NSEA2 = 0
      NSEA_DEF = 0
      PV = 1
      NSP = 0
      LVLHI = 0
      NSHIFT = 0
      DO ISEA=1,NSEA1
         IF ( QTREE%CELL_TYPE(ISEA).EQ.QTREE%UNDEF_TYPE ) THEN
            IQ = 0
         ELSE
            IQ = QTREE%INDQUAD(ISEA)
         END IF
         IF ( DOPART ) PV = QSPARE%PARTNO(ISEA)
         IF ( DOSPARE ) NSP = QSPARE%NSPARE(PV)
         IF ( IQ.NE.0 ) THEN
            NSEA2 = NSEA2 + 1
            NSEA_DEF = NSEA_DEF + 1
            QTREE%INDQUAD(NSEA2) = NEWQUAD(IQ)
            IF(NSEA2.NE.ISEA) THEN
               ISUB = QTREE%INDSUB(ISEA)

               QTREE%QICELL(IQ,ISUB) = NSEA2
               QTREE%INDSUB(NSEA2) = ISUB
               QTREE%INDLVL(NSEA2) = QTREE%INDLVL(ISEA)
               QTREE%XYVAL(NSEA2,1) = QTREE%XYVAL(ISEA,1)
               QTREE%XYVAL(NSEA2,2) = QTREE%XYVAL(ISEA,2)
               QTREE%INDML(NSEA2) = QTREE%INDML(ISEA)
               QTREE%CELL_TYPE(NSEA2) = QTREE%CELL_TYPE(ISEA)
               NOTCRS = .TRUE.
               DO ICRS=1,NCRS
                  IF (ISEANEW(ICRS).EQ.ISEA) THEN
                     ISEANEW(ICRS) = NSEA2
                     NOTCRS = .FALSE.
                     EXIT
                  END IF
               END DO
               IF ( NOTCRS) THEN
                  NSHIFT = NSHIFT + 1
                  ISEAOLD(NCRS+NSHIFT,1) = ISEA
                  ISEANEW(NCRS+NSHIFT) = NSEA2
               END IF
            END IF
            LVLHI = MAX(LVLHI,QTREE%INDLVL(NSEA2))
         ELSE IF ( DOSPARE .AND. NSP.GE.0 ) THEN
            NSEA2 = NSEA2 + 1
            NSP = NSP + 1
            IF (DOPART) THEN
               JSEA = QSPARE%JSVAL(ISEA)
            ELSE
               JSEA = ISEA
            END IF
            QSPARE%JSSPARE(NSP,PV) = JSEA
            QSPARE%NSPARE(PV) = NSP
         END IF
      END DO

      DO ISEA=NSEA2+1,NSEA1
         QTREE%INDQUAD(ISEA) = 0
         QTREE%INDSUB(ISEA) = 0
         QTREE%INDLVL(ISEA) = 0
         QTREE%XYVAL(ISEA,1) = 0.
         QTREE%XYVAL(ISEA,2) = 0.
         QTREE%INDML(ISEA) = 0
         QTREE%CELL_TYPE(ISEA) = QTREE%UNDEF_TYPE
      END DO
!
!  Reverse the spare cells lists, so the lowest numbers will be used first
!
      IF ( DOSPARE ) THEN
         DO PV=1,NPART
            NSP = QSPARE%NSPARE(PV)
            IF (NSP.GT.1) THEN
               DO II=1,NSP/2
                  ISTMP = QSPARE%JSSPARE(II,PV)
                  QSPARE%JSSPARE(II,PV) = QSPARE%JSSPARE(NSP-II+1,PV)
                  QSPARE%JSSPARE(NSP-II+1,PV) = ISTMP
               END DO
            END IF
         END DO
      END IF
!
!  Loop through removing unused quads, and renumbering
!
      NQUAD2 = NQUAD1 - NDELQ
      DO IQ = 1,NQUAD1
         IQNEW = NEWQUAD(IQ)
         IF ( IQNEW.GT.0 ) THEN
            QTREE%QLEVEL(IQNEW) = QTREE%QLEVEL(IQ)
            IQPAR = QTREE%QPARENT(IQ)
            IF ( IQPAR.GT.0 ) THEN
                QTREE%QPARENT(IQNEW) = NEWQUAD(IQPAR)
            ELSE
                QTREE%QPARENT(IQNEW) = 0
            END IF
            QTREE%QICELL(IQNEW,0) = QTREE%QICELL(IQ,0) 
            DO II=1,4
                QTREE%QICELL(IQNEW,II) = QTREE%QICELL(IQ,II) 
                IQC = QTREE%QCHILD(IQ,II) 
                IF ( IQC.GT.0 ) THEN
                   QTREE%QCHILD(IQNEW,II) = NEWQUAD(IQC)
                ELSE
                   QTREE%QCHILD(IQNEW,II) = 0
                END IF
                IQC = QTREE%QNBR(IQ,II) 
                IF ( IQC.GT.0 ) THEN
                   QTREE%QNBR(IQNEW,II) = NEWQUAD(IQC)
                ELSE
                   QTREE%QNBR(IQNEW,II) = 0
                END IF
            END DO
         END IF
      END DO
      DO IQ = NQUAD2+1,NQUAD1
         QTREE%QLEVEL(IQ) = 0
         QTREE%QPARENT(IQ) = 0
         QTREE%QICELL(IQ,0) = 0
         DO II=1,4
            QTREE%QICELL(IQ,II) = 0
            QTREE%QCHILD(IQ,II) = 0
         END DO
      END DO
!
      QTREE%NCELL = NSEA2
      QTREE%NQUAD = NQUAD2
      QTREE%LVLHI = LVLHI
      QTREE%NCELL_DEF = NSEA_DEF
!
      RETURN
      END SUBROUTINE QA_COARSEN
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_CPQT( QTREE1, QTREE2 )
!/
!       Richard Gorman, NIWA
!         May, 2014:    Origination.
!
!  1. Purpose :
!
!     Copy contents of one quadtree structure to another
!
!  2. Method :
!
!
!  3. Parameters :
!
!     ----------------------------------------------------------------
!
!       QTREE1 QA_TREE  I   Input quadtree structure
!       QTREE2 QA_TREE I/O  Output quadtree structure
!
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!
!  Subroutine parameters:
!      
      TYPE(QA_TREE), INTENT(IN)        :: QTREE1
      TYPE(QA_TREE), INTENT(INOUT)     :: QTREE2
!
! Local variables:
!
      INTEGER               :: MSEA, MQUAD, MWT
      !
      QTREE2%NCELL      = QTREE1%NCELL
      QTREE2%NCELL_DEF  = QTREE1%NCELL_DEF
      QTREE2%NQUAD      = QTREE1%NQUAD
      QTREE2%NX0        = QTREE1%NX0
      QTREE2%NY0        = QTREE1%NY0
      QTREE2%LVLREF     = QTREE1%LVLREF
      QTREE2%LVLMAX     = QTREE1%LVLMAX
      QTREE2%LVLHI      = QTREE1%LVLMAX
      QTREE2%UNDEF_TYPE = QTREE1%UNDEF_TYPE
      QTREE2%KEEP_REF   = QTREE1%KEEP_REF
      QTREE2%DYNAMIC    = QTREE1%DYNAMIC
      QTREE2%IWTORDER   = QTREE1%IWTORDER
!
! Transfer quad arrays to output versions
!
      MQUAD = SIZE(QTREE1%QPARENT,1)
      MQUAD = MIN(MQUAD,SIZE(QTREE2%QPARENT,1))
!      
      QTREE2%QPARENT = 0
      QTREE2%QLEVEL = 0
      QTREE2%QICELL = 0
      QTREE2%QNBR = 0
      QTREE2%QCHILD = 0
!      
      QTREE2%QPARENT(1:MQUAD)  = QTREE1%QPARENT(1:MQUAD)
      QTREE2%QLEVEL(1:MQUAD)   = QTREE1%QLEVEL(1:MQUAD)
      QTREE2%QICELL(1:MQUAD,:) = QTREE1%QICELL(1:MQUAD,:)
      QTREE2%QNBR(1:MQUAD,:)   = QTREE1%QNBR(1:MQUAD,:)
      QTREE2%QCHILD(1:MQUAD,:) = QTREE1%QCHILD(1:MQUAD,:)
!
! Transfer seapoint arrays to output versions
!      
      MSEA = SIZE(QTREE1%INDQUAD,1)
      MSEA = MIN(MSEA,SIZE(QTREE2%INDQUAD,1))
!      
      QTREE2%INDQUAD = 0
      QTREE2%INDSUB = 0
      QTREE2%XYVAL = 0.
      QTREE2%NGBR = 0
      QTREE2%INDLVL = 0
      QTREE2%CELL_TYPE = QTREE1%UNDEF_TYPE
      QTREE2%INDML     = 0
      !
      QTREE2%INDQUAD(1:MSEA)   = QTREE1%INDQUAD(1:MSEA)
      QTREE2%INDSUB(1:MSEA)    = QTREE1%INDSUB(1:MSEA)
      QTREE2%XYVAL(1:MSEA,:)   = QTREE1%XYVAL(1:MSEA,:)
      QTREE2%NGBR(1:MSEA,:)    = QTREE1%NGBR(1:MSEA,:)
      QTREE2%INDLVL(1:MSEA)    = QTREE1%INDLVL(1:MSEA)
      QTREE2%CELL_TYPE(1:MSEA) = QTREE1%CELL_TYPE(1:MSEA)
      QTREE2%INDML(1:MSEA)     = QTREE1%INDML(1:MSEA)
      !
      IF ( QTREE1%IWTORDER .GE. 1 ) THEN
         MSEA = SIZE(QTREE1%NCASE,1)
         MSEA = MIN(MSEA,SIZE(QTREE2%NCASE,1))
         QTREE2%NCASE = 0
         QTREE2%NCASE(1:MSEA)  = QTREE1%NCASE(1:MSEA)
      END IF
      !
      IF ( QTREE1%IWTORDER .GE. 2 ) THEN
         MSEA = SIZE(QTREE1%INDWT,1)
         MSEA = MIN(MSEA,SIZE(QTREE2%INDWT,1))
         MWT = SIZE(QTREE1%GNBR,2)
         MWT = MIN(MSEA,SIZE(QTREE2%GNBR,2))
         QTREE2%INDWT = 0
         QTREE2%NGNBR = 0
         QTREE2%GNBR  = 0
         QTREE2%INDWT(1:MSEA)      = QTREE1%INDWT(1:MSEA)
         QTREE2%NGNBR(1:MSEA)      = QTREE1%NGNBR(1:MSEA)
         QTREE2%GNBR(1:MSEA,1:MWT) = QTREE1%GNBR(1:MSEA,1:MWT)
      END IF
!
      RETURN
!
      END SUBROUTINE QA_CPQT
!
!/ ------------------------------------------------------------------- /
!
      
      SUBROUTINE QA_D2WTS( QTREE, QWTS1, QWTS2, IERR, NDSE,           &
                           NCALC, ICALC )
!/
!       Richard Gorman, NIWA
!         October, 2010:    Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!     Compute weights for derivatives at the centres of the
!     active cells, using up to second order neighbours 
!
!  2. Method :
!
!
!  3. Parameters :
!
!     ----------------------------------------------------------------
!
!     QTREE  QA_TREE   I   Quadtree structure, including
!         NSEA     Int.    Actual number of cells
!         IWTORDER Int.    Max. order of interpolation weights
!         NGBR     I.A.    [MSEA x 0:NNGBR] array of neighbours of each cell
!         NCASE    I.A.    Index (for each cell) into the precomputed tables of 
!			   possible configurations of 1st order neighbours
!         INDWT    I.A.   O   Index, for each cell, into the lookup  table of 2nd
!                             order weights
!         GNBR     I.A.   O   Grand neighbours of each cell
!         NGNBR    I.A.   O   Number of grand neighbours of each cell
!     QWTS1  QA_WEIGHTS1 I  Structure of first order weights, including:
!         IQUADDIR I.A     [MCASE x 4] array coding for type of neighbour(s) in
!                            [W,E,S,N] direction for each possible configuration case
!                            = 1 for a neighbour of equal level
!                            = 2 for two neighbours of higher level
!                            = 3 for a neighbour of lower level, displaced in
!                                -ve direction (S for E,W nbrs, W for S,N nbrs)
!                            = 4 for a neighbour of lower level, displaced in 
!                                +ve direction (N for E,W nbrs, E for S,N nbrs)
!                            = 5 for one neighbour of higher level, displaced in 
!                               -ve direction (S for E,W nbrs, W for S,N nbrs)
!                            = 6 for one neighbour of higher level, displaced in 
!                                +ve direction (N for E,W nbrs, E for S,N nbrs)
!                            = 7 for no neighbour
!         LAMBDA0  R.A     Natural neighbour coordinates relative 
!                              to cell centre
!     QWTS2  QA_WEIGHTS2 I/O  Structure of second order weights, including:
!         NCSAV    Int.    Actual number of configurations of 1st and 2nd order 
!                           neighbours that have been saved in the look up table 
!         IC2SAV   I.A.    Array of configurations of 1st and 2nd order neighbours
!                           for which weights have been saved in the lookup table.
!         DERWTSAV R.A.    Look up table of weights for derivatives at cell centres:
!                           DERWTSAV(ISAV,NB,IDER) = weight for "type IDER"
!                           derivative at the cell centre, from nearby cell NB, i.e.
!                           for a variable A, it's various derivatives are
!                             D(ISEA) =     sum     DERWTSAV(ISAV,NB,IDER)*A(JSEA(NB))
!                                       {NB=0:MGNBR} 
!                           where ISAV = INDWT(ISEA), JSEA(NB) = {ISEA          for NB=0
!                                                                {GNBR(ISEA,NB) for NB>0
!                           and "D" is d/dX for IDER=1,
!                                      d/dY          2,
!                                     d2/dx2         3,
!                                     d2/dY2         4,
!                                     d2/dXdY        5,
!                                     d3/dX3         6,
!                                     d3/dY3         7,
!                                     d3/dX2dY       8,
!                                     d3/dXdY2       9
!
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!       NCALC    Int.   I*  Number of cells to process (default = all)
!       ICALC    Int.   I*  Array of cell indices to process
!
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     QA_DERIVWTS Subr. qa_utils Compute weights for higher order derivatives
!                               on an extended stencil 
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
      INTEGER, PARAMETER    :: NNGBR=8, MGNBR=36
      TYPE(QA_TREE), INTENT(INOUT)     :: QTREE
      TYPE(QA_WEIGHTS1), INTENT(IN)    :: QWTS1
      TYPE(QA_WEIGHTS2), INTENT(INOUT) :: QWTS2
      INTEGER, OPTIONAL, INTENT(IN)    :: NCALC
      INTEGER, OPTIONAL, INTENT(IN)    :: ICALC(:)
      INTEGER, OPTIONAL, INTENT(OUT)   :: IERR
      INTEGER, OPTIONAL, INTENT(IN)    :: NDSE
!
! Local variables:
!
      INTEGER               :: ISEA, I, J, NGN, IGN, ICASE
      INTEGER               :: ISUB, IDIR, ICIND, NB, NTYPE, KCIND,   &
                               KSEA, KSUB, KDIR, QRET, KCASE, ICS
      INTEGER               :: ICASE2(0:NNGBR)
      INTEGER               :: GNEIGH(MGNBR)
      INTEGER               :: NGNEIGH
      REAL                  :: LAM(MGNBR), X(MGNBR), Y(MGNBR),        &
                               X0, Y0, XYSCALE
      REAL                  :: WDX(0:MGNBR), WDY(0:MGNBR),            &
                               WDXX(0:MGNBR), WDXY(0:MGNBR),          &
                               WDYY(0:MGNBR), WDXXX(0:MGNBR),         &
                               WDXXY(0:MGNBR), WDXYY(0:MGNBR),        &
                               WDYYY(0:MGNBR)
      REAL                  :: CSIZE(7), XCASE(8,7), YCASE(8,7)
      REAL                  :: DEFLAMBDA, DET
      LOGICAL               :: NOMATCH, NEWCASE
      INTEGER               :: NSEEK, ISEEK
      INTEGER               :: IUN, IERS
!      
      INTEGER               :: MCSAV
      INTEGER               :: MDSAV
      INTEGER               :: NCSAV
!
    
      DATA  DEFLAMBDA / 0.25 /
      DATA  (CSIZE(J),J=1,7) /                                        &
              1.0, 0.5, 2.0, 2.0, 0.5, 0.5, 9. /

      DATA  ((XCASE(I,J),J=1,7),I=1,8) /                              &
             -1., -0.75, -1.25, -1.25, -0.75,   9.,  9.,              &
              1.,  0.75,  1.25,  1.25,  0.75,   9.,  9.,              &
              0., -0.25, -0.5,   0.5,  -0.25,   9.,  9.,              &
              0., -0.25, -0.5,   0.5,  -0.25,   9.,  9.,              &
              9., -0.75,   9.,    9.,     9., -0.75, 9.,              &
              9.,  0.75,   9.,    9.,     9.,  0.75, 9.,              &
              9.,  0.25,   9.,    9.,     9.,  0.25, 9.,              &
              9.,  0.25,   9.,    9.,     9.,  0.25, 9. /
      DATA  ((YCASE(I,J),J=1,7),I=1,8) /                              &
              0., -0.25, -0.5,   0.5,  -0.25,   9.,  9.,              &
              0., -0.25, -0.5,   0.5,  -0.25,   9.,  9.,              &
             -1., -0.75, -1.25, -1.25, -0.75,   9.,  9.,              &
              1.,  0.75,  1.25,  1.25,  0.75,   9.,  9.,              &
              9.,  0.25,   9.,    9.,     9.,  0.25, 9.,              &
              9.,  0.25,   9.,    9.,     9.,  0.25, 9.,              &
              9., -0.75,   9.,    9.,     9., -0.75, 9.,              &
              9.,  0.75,   9.,    9.,     9.,  0.75, 9. /
!
!
      NCSAV = QWTS2%NCSAV
      MCSAV = SIZE(QWTS2%IC2SAV,1)
      MDSAV = SIZE(QWTS2%DERWTSAV,3)
!
!   Determine which cells to process:
!   If ICALC and NCALC are specified, process the first NCALC indices
!   in the array ICALC.
!   If NCALC is provided but not ICALC, process the first NCALC cells
!   If NCALC is not provided, or is zero, process all cells
!
      IF ( PRESENT(NCALC) ) THEN
         IF ( NCALC.LE.0 ) THEN
            NSEEK = QTREE%NCELL
         ELSE
            NSEEK = MIN(NCALC, QTREE%NCELL)
         END IF
      ELSE
         NSEEK = QTREE%NCELL
      END IF
      IF ( PRESENT(ICALC) ) THEN
         NSEEK = MIN(SIZE(ICALC),NSEEK)
      END IF
      
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0
      IF ( QTREE%IWTORDER.LT.2 ) THEN
         IF ( IUN.GT.0 )  THEN
            WRITE(IUN,*) 'ERROR IN QA_D2WTS: IWTORDER MUST BE >=2'
            IF ( PRESENT(IERR) ) IERR = 1
            RETURN 
         END IF
      END IF
!
      QTREE%GNBR = 0     
!
!     Loop over cells
!
      SEALOOP: DO ISEEK=1,NSEEK
         IF ( PRESENT(ICALC) ) THEN
            ISEA = ICALC(ISEEK)
         ELSE
            ISEA = ISEEK
         END IF
         !
         IF ( QTREE%CELL_TYPE(ISEA).EQ.QTREE%UNDEF_TYPE ) CYCLE
         GNEIGH = 0
         NGNEIGH = 0
         ICASE2 = 0
         !
         !  Loop over the 2 possible neighbours per side
         ICASE = QTREE%NCASE(ISEA)
         ICASE2(0) = ICASE
         DO ICIND=1,8
            NB = QTREE%NGBR(ISEA,ICIND)
            IF ( NB.GT.0 ) THEN
               ICASE2(ICIND) = QTREE%NCASE(NB)
            END IF
         END DO
         !
         ! See if this configuration is in the table of already-computed weights:
         ! 
         NEWCASE = .TRUE.
         DO ICS=1,NCSAV
            IF ( ALL(ICASE2.EQ.QWTS2%IC2SAV(ICS,:)) ) THEN
               QTREE%INDWT(ISEA) = ICS
               NEWCASE = .FALSE.
               EXIT
            END IF
         END DO
	 !
	 ! If this is a new configuration, make a new entry in the lookup tables,
	 ! 
         IF (NEWCASE) THEN
            IF ( NCSAV.LT.MCSAV ) THEN
               NCSAV = NCSAV + 1
               QWTS2%NCSAV = NCSAV
               QTREE%INDWT(ISEA) = NCSAV
               QWTS2%IC2SAV(NCSAV,:) = ICASE2
            ELSE
               IF ( IUN.GT.0 )  THEN
                  WRITE(IUN,*) 'ERROR IN QA_D2WTS:'
                  WRITE(IUN,*) 'INSUFFICIENT SPACE IN WEIGHTS TABLE'
                  WRITE(IUN,*) 'Reallocate arrays with size >',MCSAV
               END IF
               IF ( PRESENT(IERR) ) IERR = 1
               RETURN 
            END IF
         END IF
	 !
         ! Get relative X, Y and Natural Neighbour coordinates of 
         ! neighbours
         ! Loop over 1st order neighbours
         NGN = 0
         !  Loop over the 2 possible neighbours per side
         DO ISUB=0,1
            !  Loop over the four sides [W E S N]:
            DO IDIR=1,4
               ICIND = 4*ISUB+IDIR
               KSEA = QTREE%NGBR(ISEA,ICIND)
               IF ( KSEA.GT.0 ) THEN
                  NGN = NGN+1
                  GNEIGH(NGN) = KSEA
                  QTREE%GNBR(ISEA,NGN) = KSEA
                  KCASE = QTREE%NCASE(KSEA)
                  IF ( NEWCASE ) THEN
                     NTYPE = QWTS1%IQUADDIR(ICASE,IDIR)
                     X0 = XCASE(ICIND,NTYPE)
                     Y0 = YCASE(ICIND,NTYPE)
                     X(NGN) = X0
                     Y(NGN) = Y0
                     IF ( ICASE.GT.0 ) THEN
                        LAM(NGN) = QWTS1%LAMBDA0(ICASE,ICIND)
                     ELSE
                        LAM(NGN) = DEFLAMBDA
                     END IF
                     XYSCALE = CSIZE(NTYPE)
                  END IF
                  !  Loop over neighbours of KSEA 
                  !  Loop over the 2 possible neighbours per side
                  DO KSUB=0,1
                     !  Loop over the four sides [W E S N]:
                     DO KDIR=1,4
                        KCIND = 4*KSUB+KDIR
                        NB = QTREE%NGBR(KSEA,KCIND)
                        ! Check that this 2nd order neighbour exists, and isn't the
                        ! central cell:
                        IF ( NB.GT.0 .AND. NB.NE.ISEA ) THEN
                           !  Check that it hasn't already been counted              
                           NOMATCH = .TRUE.
                           DO IGN=1,NGN
                              IF (GNEIGH(IGN).EQ.NB) THEN
                                 NOMATCH = .FALSE.
                                 EXIT
                              END IF
                           END DO
                           IF ( NOMATCH ) THEN
                              ! New cell: add it to ISEA's grand-neighbour list
                              NGN = NGN + 1
                              IF ( NGN.GT.MGNBR ) THEN
                                 IF ( IUN.GT.0 )  THEN
                                    WRITE(IUN,*) 'ERROR IN QA_D2WTS:'
                                    WRITE(IUN,*) 'NUMBER OF 2ND ',    &
                                       'ORDER NEIGHBOURS > ALLOCATION'
                                    WRITE(IUN,*) 'NGN:', NGN
                                    WRITE(IUN,*) 'MGNBR:',MGNBR
                                    WRITE(IUN,*) 'ISEA:',ISEA
                                    WRITE(IUN,*) 'NBR:',              &
                                       (QTREE%NGBR(ISEA,I),I=1,NNGBR)
                                    WRITE(IUN,*) 'GNEIGH, X, Y:'
                                    DO I=1,MGNBR
                                       WRITE(IUN,*) GNEIGH(I),X(I),Y(I)
                                    END DO
                                    WRITE(IUN,*) 'NEW GNEIGH:', NB
                                 END IF
                                 IF ( PRESENT(IERR) ) IERR = 1
                                 RETURN 
                              END IF   
                              GNEIGH(NGN) = NB
                              QTREE%GNBR(ISEA,NGN) = NB
                              IF ( NEWCASE ) THEN
                                 NTYPE = QWTS1%IQUADDIR(KCASE,KDIR)
                                 X(NGN) = X0 + XYSCALE*XCASE(KCIND,NTYPE)
                                 Y(NGN) = Y0 + XYSCALE*YCASE(KCIND,NTYPE)
                                 LAM(NGN) = DEFLAMBDA
                              END IF
                           END IF
                        END IF
                     END DO
                  END DO   !  Loop over neighbours of KSEA 
               END IF
            END DO
        END DO  ! Loop over neighbours (KSEA) of ISEA
        !
        QTREE%NGNBR(ISEA) = NGN
        IF ( NEWCASE ) THEN
           QRET = 3
           CALL QA_DERIVWTS( MGNBR, NGN, X, Y, LAM, QRET, WDX, WDY,   &
                WDXX, WDYY, WDXY, WDXXX, WDYYY, WDXXY, WDXYY, DET,    &
                ierr=IERS, ndse=IUN )
           IF ( IERS.GT.0 )  THEN
              IF ( IUN.GT.0 )  THEN
                 WRITE(IUN,*) 'ERROR IN QA_D2WTS CALLING QA_DERIVWTS'
                 WRITE(IUN,*) 'ISEA:', ISEA
                 WRITE(IUN,*) 'NBR:', (QTREE%NGBR(ISEA,I),I=1,NNGBR)
                 WRITE(IUN,*) 'NGN:', NGN
                 WRITE(IUN,*) 'GNEIGH, X, Y:'
                 DO I=1,NGN
                    WRITE(IUN,*) GNEIGH(I),X(I),Y(I)
                 END DO
              END IF
              IF ( PRESENT(IERR) ) IERR = IERS
              RETURN 
           END IF
           !
    	   !  Load the computed weights into the lookup table
           !
           DO IGN=0,NGN
              IF (QRET.GE.1) THEN
                 IF (MDSAV.GE.1) QWTS2%DERWTSAV(NCSAV,IGN,1) = WDX(IGN)
                 IF (MDSAV.GE.2) QWTS2%DERWTSAV(NCSAV,IGN,2) = WDY(IGN)
              END IF
              IF (QRET.GE.2) THEN
                 IF (MDSAV.GE.3) QWTS2%DERWTSAV(NCSAV,IGN,3) = WDXX(IGN)
                 IF (MDSAV.GE.4) QWTS2%DERWTSAV(NCSAV,IGN,4) = WDYY(IGN)
                 IF (MDSAV.GE.5) QWTS2%DERWTSAV(NCSAV,IGN,5) = WDXY(IGN)
              END IF
              IF (QRET.GE.3) THEN
                 IF (MDSAV.GE.6) QWTS2%DERWTSAV(NCSAV,IGN,6) = WDXXX(IGN)
                 IF (MDSAV.GE.7) QWTS2%DERWTSAV(NCSAV,IGN,7) = WDYYY(IGN)
                 IF (MDSAV.GE.8) QWTS2%DERWTSAV(NCSAV,IGN,8) = WDXXY(IGN)
                 IF (MDSAV.GE.9) QWTS2%DERWTSAV(NCSAV,IGN,9) = WDXYY(IGN)
              END IF
           END DO
        END IF
        !write(*,*) 'ARRAYS transferred'
!  End of loop over cells        
      END DO SEALOOP
!
      END SUBROUTINE QA_D2WTS
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_DERIVC( ICELL, ITBL, NSOL, ISOL, NSTEN, ISTEN,    &
                            VAL, WEIGHT, DERIV, DELXIN, DELYIN, IERR, &
                            NDSE )
!/
!       Richard Gorman, NIWA
!         14-Oct-2010:      Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!/
!  1. Purpose :
!
!     Compute up to 3rd order derivatives at a cell centre using weighted 
!     values from  neighbouring cells.
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ICELL    Int.   I   Target cell
!       MSTEN    Int.   I   Max. number of stencil cells
!       NSOL     Int.   I   Number of derivatives to be calculated (up to 9)
!       ISOL     I.A.   I   Array of codes for derivatives to be calculated:
!                            1 = d/dX, 2=d/dY, 3=d2/dX2, 4=d2/dY2, 5= d2/dXdY 
!                            6 = d3/dX3, 7=d3/dY3, 8=d3/dX2dY, 9=d3/dXdY2
!       NSTEN    I.A.   I   Number of stencil cells around each cell
!       ISTEN    I.A.   I   indices of stencil cells
!       VAL      R.A.   I   Data values at all cells
!       WEIGHT   R.A    I   Array of weights 
!                            WEIGHT(ITBL,NB,ISOL(K)) = weight to compute 
!                            DERIV(ISOL(K)) at cell ISEA from values at cell 
!                            ISTEN(NB), or ISEA if NB=0.
!       DELXIN   Real   I   1/deltaX for the target cell
!       DELYIN   Real   I   1/deltaY for the target cell
!       DERIV    R.A    O   Computed derivatives at cell centres
!                           DERIV(K) = derivative defined by ISOL(K)
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       None:
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_INTERP  Subr. qa_utils Interpolate to a target location 
!                               using weighted values from 
!                               neighbouring cells.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: ICELL, ITBL, NSOL
      INTEGER, INTENT(IN)     :: ISOL(NSOL)
      INTEGER, INTENT(IN)     :: NSTEN, ISTEN(:)     
      REAL, INTENT(IN)        :: VAL(:)
      REAL, INTENT(IN)        :: WEIGHT(1:,0:,1:)
      REAL, INTENT(OUT)       :: DERIV(NSOL)
      REAL, OPTIONAL, INTENT(IN)      :: DELXIN, DELYIN
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER             :: NGB, INBR, IS, IDER
      INTEGER             :: MTBL, MDER, MSTEN, MCELL, MNGU
      INTEGER             :: IUN
      REAL                :: DELXINV, DELYINV
!
      DELXINV = 1.
      DELYINV = 1.
      IF ( PRESENT(DELXIN) ) DELXINV = DELXIN
      IF ( PRESENT(DELYIN) ) DELYINV = DELYIN
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE      
      IF ( PRESENT(IERR) ) IERR = 0
!      
      MCELL = SIZE(VAL,1)
      MTBL = SIZE(WEIGHT,1)
      !MNGL = LBOUND(WEIGHT,2)
      MNGU = UBOUND(WEIGHT,2)
      MDER = SIZE(WEIGHT,3)
      MSTEN = SIZE(ISTEN,1)
      IF ( NSTEN.GT.MNGU )  THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_DERIVC:'
            WRITE(IUN,*) 'NEIGHBOUR INDEX OUT OF BOUNDS OF WEIGHTS TABLE'
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
      DERIV = 0.
      IF ( ICELL.LE.0 .OR. ICELL.GT.MCELL )  THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_DERIVC: '
            WRITE(IUN,*) 'CELL INDEX OUT OF BOUNDS OF VAL ARRAY '
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
      IF ( ITBL.LE.0 .OR. ITBL.GT.MTBL )  THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_DERIVC:'
            WRITE(IUN,*) 'TABLE INDEX OUT OF BOUNDS OF WEIGHTS TABLE'
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
      DO IS=1,NSOL
         IDER = ISOL(IS)
         IF ( IDER.LE.0 .OR. IDER.GT.MDER )  THEN
            IF ( IUN.GT.0 ) THEN
               WRITE(IUN,*) 'ERROR IN QA_DERIVC:'
               WRITE(IUN,*) 'DERIVATIVE INDEX OUT OF BOUNDS OF WEIGHTS TABLE'
            END IF
            IF ( PRESENT(IERR) ) IERR = 1
            RETURN
         END IF
         DERIV(IS) = VAL(ICELL)*WEIGHT(ITBL,0,IDER)
         DO INBR=1,NSTEN
            NGB = ISTEN(INBR)
            IF ( NGB.GT.0 .AND. NGB.LE.MCELL ) THEN
               DERIV(IS) = DERIV(IS) + VAL(NGB)*WEIGHT(ITBL,INBR,IDER)
            END IF
         END DO
         IF (IDER.EQ.1) THEN
            DERIV(IS) = DERIV(IS)*DELXINV
         ELSEIF (IDER.EQ.2) THEN
            DERIV(IS) = DERIV(IS)*DELYINV
         ELSEIF (IDER.EQ.3) THEN
            DERIV(IS) = DERIV(IS)*DELXINV*DELXINV
         ELSEIF (IDER.EQ.4) THEN
            DERIV(IS) = DERIV(IS)*DELYINV*DELYINV
         ELSEIF (IDER.EQ.5) THEN
            DERIV(IS) = DERIV(IS)*DELXINV*DELYINV
         ELSEIF (IDER.EQ.6) THEN
            DERIV(IS) = DERIV(IS)*DELXINV**3
         ELSEIF (IDER.EQ.7) THEN
            DERIV(IS) = DERIV(IS)*DELYINV**3
         ELSEIF (IDER.EQ.8) THEN
            DERIV(IS) = DERIV(IS)*DELXINV*DELXINV*DELYINV
         ELSEIF (IDER.EQ.9) THEN
            DERIV(IS) = DERIV(IS)*DELXINV*DELYINV*DELYINV
         END IF
      END DO
!/
!/ End of QA_DERIVC ----------------------------------------------------- /
!/
      END SUBROUTINE QA_DERIVC
!
!/ ------------------------------------------------------------------- /
!
      RECURSIVE SUBROUTINE QA_DERIVWTS( NM, N, X, Y, LAMBDA, QMAX,    &
                                        WDX, WDY, WDXX, WDYY, WDXY,   &
                                        WDXXX, WDYYY, WDXXY, WDXYY,   &
                                        DET, IERR, NDSE )
!/
!       Richard Gorman, NIWA
!         07-Oct-2010:      Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!/
!  1. Purpose :
!
!     Compute weights for higher order derivatives on an extended stencil 
!     (i.e. beyond 1st order neighbours)
!
!  2. Method :
!
!     We use an extension of the method of Sibson which estimates the gradient 
!     of a function at a cell by least squares fit of the gradient plane through 
!     its natural neighbours.
!     In this extended version, a wider set of cells can be used, fitted by a
!     higher level polynomial, i.e.
!       Zfit(x,y) =  Z(0) + Sum{  Sum{  a(m,n) * x^m * y^n,  with n=q-m
!                          q=1,Q  m=0,q
!
!     We wish to compute the factors a(m,n), i.e. various derivatives up to order Q 
!     at the origin, in terms of weighted sums of the values known at a set of 
!     surrounding points, with coordinates X(k), Y(k), k=1,N, i.e.
!         a(m,n) = Sum{ w(m,n,k) * Z(k)
!                 k=0,N
!     with k=0 being the origin (x,y) = (0,0).
!
!     The least squares fit minimises
!       eta =  Sum{ s(k) *( Zfit(X(k),Y(k)) - Z(k) )^2
!             k=1,N
!     where Z(k) is the data value at the kth cell.
!     Following Sibson, we take the weights s(k) to be 
!           s(k) = lambda(k)/(X(k)^2 + Y(k)^2), where lambda(k) is the Natural
!     Neighbour coordinate, where known (i.e. for first order neighbours),
!     and a default constant otherwise: one option is to use 0.25, the lambda 
!     value for a rectangular grid
!     Taking derivatives with respect to each parameter a(m,n) and setting to 
!     zero results in a linear system
!               
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NM      Int.   I   Maximum number of cells
!       N       Int.   I   Actual number of cells
!       X,Y     R.A.   I   coordinates of cells, relative to the target cell
!                          at (0,0), normalised by the target cell size
!       LAMBDA  R.A.   I   Natural neighbour coordinates of the cells
!       QMAX    Int.   I   Highest order derivatives to solve for
!       QRET    Int.   IO  Highest order derivatives returned (up to 3)
!       WDX     R.A.   O   Weights for d/dX at target cell
!       WDY     R.A.   O   Weights for d/dY at target cell
!       WDXX    R.A.   O   Weights for d2/dX2 at target cell
!       WDXY    R.A.   O   Weights for d2/dXdY at target cell
!       WDYY    R.A.   O   Weights for d2/dY2 at target cell
!       WDXXX   R.A.   O   Weights for d3/dX3 at target cell
!       WDXXY   R.A.   O   Weights for d3/dX2dY at target cell
!       WDXYY   R.A.   O   Weights for d3/dXdY2 at target cell
!       WDYYY   R.A.   O   Weights for d3/dY3 at target cell
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     Name        Type  Module   Description
!     ---------------------------------------------------------------
!     QA_DERIVWTS Subr. qa_utils Compute weights for higher order derivatives
!                                on an extended stencil 
!     QA_INCR     Subr. qa_utils Increment a multiple-digit counter in 
!                                arbitrary base
!     QA_LINGAUSS Subr. qa_utils Solve a linear sytem by Gaussian elimination
!     svdcmp      Subr. SVD      Singular value decomposition
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     Name        Type  Module   Description
!     ---------------------------------------------------------------
!     QA_DERIVWTS Subr. qa_utils Compute weights for higher order derivatives
!                                on an extended stencil 
!     QA_D2WTS    Subr. qa_utils Compute weights for higher derivatives
!                                using second-order neighbours
!     QA_DERWTS   Subr. qa_utils Creates arrays of precomputed 
!                                interpolation weights for all
!                                possible arrangements of quadtree 
!                                neighbours
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NM, N
      REAL, INTENT(IN)        :: X(NM), Y(NM), LAMBDA(NM)
      INTEGER, INTENT(INOUT)  :: QMAX
      REAL, INTENT(OUT)       :: WDX(0:NM), WDY(0:NM)
      REAL, INTENT(OUT)       :: WDXX(0:NM), WDYY(0:NM), WDXY(0:NM)
      REAL, INTENT(OUT)       :: WDXXX(0:NM), WDYYY(0:NM),            &
                                 WDXXY(0:NM), WDXYY(0:NM)
      REAL, INTENT(OUT)       :: DET
      INTEGER, INTENT(OUT), OPTIONAL    :: IERR
      INTEGER, INTENT(IN), OPTIONAL     :: NDSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL, PARAMETER      ::  DETMIN = 1.E-6
      REAL                 ::  S(N)
      REAL, ALLOCATABLE    ::  SS(:,:), A(:,:), B(:,:), WT(:,:)
      REAL                 ::  LBDA, WSUM, WTEMP, WMIN
      INTEGER              ::  I, J, K, IND1, IND2, IST, KP, P, INCSTAT
      INTEGER              ::  M, NEQ, PX(9), PY(9), ITEST, NNULL
      INTEGER              ::  CHILD(9,9)
      LOGICAL, ALLOCATABLE ::  ISNULL(:), SETNULL(:)
      INTEGER, ALLOCATABLE :: IND(:)
      INTEGER              ::  LWORK, IERS
      REAL, ALLOCATABLE    :: WORK(:), AP(:,:)
      REAL, ALLOCATABLE    :: BCHK(:,:)
      REAL, ALLOCATABLE    :: U(:,:), W(:), V(:,:), BVEC(:), XS(:)
      REAL, ALLOCATABLE    :: UP(:,:), WP(:), VP(:,:), VF(:,:),       &
                              TMPP(:), XP(:)
      REAL                 :: XKP(NM), YKP(NM), XPOW, YPOW
      REAL                 :: DENOM
      INTEGER              :: IUN
      INTEGER              :: ISOLVE, I2, J2, PX1, PX2, PY1, PY2
      LOGICAL              :: RECALC
      DATA PX/1,0,2,0,1,3,0,2,1/
      DATA PY/0,1,0,2,1,0,3,1,2/
      DATA ((CHILD(I,J),J=1,9),I=1,9) / 0, 0, 1, 0, 1, 2, 0, 2, 2,    &
                                        0, 0, 0, 1, 1, 0, 2, 2, 2,    &
                                       -1, 0, 0, 0, 0, 1, 0, 1, 0,    &
                                        0,-1, 0, 0, 0, 0, 1, 0, 1,    &
                                       -1,-1, 0, 0, 0, 0, 0, 1, 1,    &
                                       -2, 0,-1, 0, 0, 0, 0, 0, 0,    &
                                        0,-2, 0,-1, 0, 0, 0, 0, 0,    &
                                       -2,-2,-1, 0,-1, 0, 0, 0, 0,    &
                                       -2,-2, 0,-1,-1, 0, 0, 0, 0 /
!      INTERFACE
!         SUBROUTINE svdcmp(a,w,v)
!            REAL, DIMENSION(:,:), INTENT(INOUT) :: a
!            REAL, DIMENSION(:), INTENT(OUT) :: w
!            REAL, DIMENSION(:,:), INTENT(OUT) :: v
!         END SUBROUTINE svdcmp
!         SUBROUTINE svbksb(u,w,v,b,x)
!            REAL, DIMENSION(:,:), INTENT(IN) :: u,v
!            REAL, DIMENSION(:), INTENT(IN) :: w,b
!            REAL, DIMENSION(:), INTENT(OUT) :: x
!         END SUBROUTINE svbksb
!      END INTERFACE
!/
!/ ------------------------------------------------------------------- /
!/
! ITEST = 1 for test output:
      ITEST = 0
! Solution method: 1 for SVD, 0 for Gaussian elimination    
      ISOLVE = 1
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0
!
! 1. Compute the size of the linear system. 
!    If this is too big for the number of data points,
!    reduce the order and call again recursively:
!
      M = 0
      DO I=1,QMAX
        M = M + I + 1
        !QFULL = I
        !IF ( M .GE. N ) EXIT
      END DO
      !M = MFULL
      M = MIN(M,9)
      !NEQ = MIN(M,N)
      NEQ = M
      !QMAX = 3
      IF(ITEST.GT.0) write(*,*) 'NM, N, NEQ, M, QMAX: ',              &
                                 NM, N, NEQ, M, QMAX
!
!   2. Default outputs:
!
      WDX = 0.
      WDY = 0.
      WDXX = 0.
      WDYY = 0.
      WDXY = 0.
      WDXXX = 0.
      WDYYY = 0.
      WDXXY = 0.
      WDXYY = 0.
!
!   3. Compute weights for the least squares fitting procedure.
!
!
      DO I = 1,N
!       3.1 Lambda may only be specified for the immediate neighbours:
!           apply a default for others:
        IF ( LAMBDA(I).LE.0. ) THEN
          LBDA = 0.25
        ELSE
          LBDA = LAMBDA(I)
        END IF
        DENOM = X(I)*X(I) + Y(I)*Y(I)
        IF ( DENOM.LT.TINY(DENOM) ) THEN
           IF ( IUN.GT.0 ) THEN
              WRITE(IUN,*) 'ERROR IN QA_DERIVWTS: (X,Y) = 0.'
           END IF
           IF ( PRESENT(IERR) ) IERR = 1
           RETURN
        END IF
        !S(I) = LBDA/SQRT( DENOM )
        S(I) = LBDA/DENOM
      END DO
      IF(ITEST.GT.0) write(*,*) 'S: ', (S(I),I=1,N)
!
!   4. Compute moments of X, Y:
!
      ALLOCATE( SS(0:2*QMAX,0:2*QMAX) )
      XKP = 1.
      DO I=0,2*QMAX
        IF (I.GT.0 ) THEN
          DO K=1,N
            XKP(K) = XKP(K)*X(K)
          END DO
        END IF
        YKP = 1.
        DO J=0,2*QMAX
          WSUM = 0.
          DO K=1,N
            IF( J.GT.0 ) THEN
               YKP(K) = YKP(K)*Y(K)
            END IF
            WSUM = WSUM + S(K)*XKP(K)*YKP(K)
          END DO
          SS(I,J) = WSUM         
        END DO
      END DO
      IF(ITEST.GT.0) write(*,*) 'SS computed'
!
!   5. Set up linear system
!
      ALLOCATE( A(M,M), B(M,N+1), WT(M,N+1) )
      A = 0.
      B = 0.
      IF(ITEST.GT.0) write(*,*) 'A, B, WT allocated'
!
!     N equations in M unknowns:
!
      DO IND1=1,NEQ
          PX1 = PX(IND1)
          PY1 = PY(IND1)
          DO IND2=1,M
            PX2 = PX(IND2)
            PY2 = PY(IND2)
            A(IND1,IND2) = SS(PX1+PX2,PY1+PY2)
          END DO
          WSUM = 0.
          DO K=1,N
            IF (PX1.EQ.0) THEN
              XPOW = 1.
            ELSE
              XPOW = X(K)**PX1
            END IF
            IF (PY1.EQ.0) THEN
              YPOW = 1.
            ELSE
              YPOW = Y(K)**PY1
            END IF
            WTEMP = S(K)*XPOW*YPOW
            B(IND1,K+1) = WTEMP
            WSUM = WSUM - WTEMP
          END DO
          B(IND1,1) = WSUM
      END DO
!
!   6. Solve the linear system.
!
      RECALC = .FALSE.
      IF ( ISOLVE.EQ.1 ) THEN
!
!        Option 1: use SVD
!
         IST = 0
         ALLOCATE( U(M,M), W(M), V(M,M), BVEC(M), XS(M) )
         ALLOCATE( ISNULL(M), SETNULL(M) )
         IF(ITEST.GT.0) write(*,*) 'U, W, V, BVEC, XS allocated'
         !  
         ! Singular Value Decomposition  
         !   A = U*W*VT
         ! option 1: Numerical Recipes
         !U = A               ! Copy a into u if you don't want it to be destroyed.
         !CALL svdcmp(U,W,V)  !  SVD the square matrix a.
         !                   ! V contains V, not VT
         !
         ! option 2: LAPACK
         LWORK = 5*M
         ALLOCATE ( WORK(LWORK) )
         CALL SGESVD( 'A', 'A', M, M, A, M, W, U, M, V, M, WORK,      &
                      LWORK, IERS )
         IF (IERS.NE.0 ) THEN
            IF ( IUN.GT.0 ) THEN
               WRITE(IUN,*) 'ERROR IN QA_DERIVWTS CALLING SGESVD(1)'
               WRITE(IUN,*) 'IRET_SVD = ', IERS
               WRITE(IUN,*) 'M = ', M
               WRITE(IUN,*) 'A:'
               DO I=1,m
                 WRITE(IUN,*) (A(i,k),k=1,m)
               END DO
            END IF
            IF ( PRESENT(IERR) ) IERR = IERS
            RETURN
         END IF
         !                   
         ! V contains VT, not V
         V = TRANSPOSE(V)
         ! Determinant:
         DET = 1.
         DO IND2=1,M
            DET = DET*W(IND2)
         END DO
         IF(ITEST.GT.0) THEN
            write(ITEST,*) 'A:'
            do i=1,m
               write(ITEST,*) (A(i,k),k=1,m)
            end do
            write(ITEST,*) 'B:'
            do i=1,m
               write(ITEST,*) (B(i,k),k=1,n+1)
            end do
            write(ITEST,*) 'SVD done'
            write(ITEST,*) 'U:'
            do i=1,m
               write(ITEST,*) (U(i,k),k=1,m)
            end do
            write(ITEST,*) 'W:'
            write(ITEST,*) (W(k),k=1,m)
            write(ITEST,*) 'V:'
            do i=1,m
               write(ITEST,*) (V(i,k),k=1,m)
            end do
         END IF
         !  Set the threshold for singular values allowed to be nonzero.
         IF(ITEST.GT.0) write(*,*) 'SVD called'
         WMIN = 1.E-6*MAXVAL(W) 
         !  Now we can backsubstitute to get the "minimum norm" (approximate) solution
         DO I=1,N+1
            BVEC = B(:,I)
            WHERE (W .GT. WMIN)
               XS = MATMUL(BVEC,U)/W          ! Calculate diag(1/wj )U T B,
            ELSEWHERE
               XS = 0.0                    ! but replace 1/wj by zero if wj is small.
            END WHERE
            WT(:,I) = MATMUL(V,XS)
         !    CALL svbksb(U,W,V,BVEC,XS) 
         !    WT(:,I) = XS
         END DO
         ! Identify and count null elements of W
         NNULL = 0
         ISNULL = .FALSE.
         DO  I=1,M
            IF ( W(I).LT.WMIN ) THEN
               NNULL = NNULL + 1
               ISNULL(I) = .TRUE.
            END IF
         END DO
         ! For each null W value, we need to choose a derivative term 
         ! to set to zero. This won't necessarily be the term corresponding
         ! directly to the null W value, rather we choose its highest order 
         ! derivative
         IF ( NNULL.GT.0 ) THEN
            SETNULL = ISNULL
            DO I=1,M
               IF ( SETNULL(I) ) THEN
                  DO J=M,1,-1
                     IF ( CHILD(I,J).GT.0 .AND. .NOT.SETNULL(J) ) THEN
                        SETNULL(I) = .FALSE.
                        SETNULL(J) = .TRUE.
                        EXIT
                     END IF 
                  END DO
               END IF 
            END DO
            !
            ALLOCATE ( VF(M,NNULL), VP(NNULL,NNULL), UP(NNULL,NNULL), WP(NNULL) )
            ALLOCATE ( TMPP(NNULL), XP(NNULL) )
            J2 = 0
            DO J=1,M
               IF ( ISNULL(J) ) THEN
                  J2 = J2 + 1
                  I2 = 0
                  DO I=1,M
                     VF(I,J2) = V(I,J)
                     IF ( SETNULL(I) ) THEN
                        I2 = I2 + 1
                        UP(I2,J2) = V(I,J)
                     END IF
                  END DO
               END IF
            END DO
            !
            ! option 1: Numerical Recipes
            !CALL svdcmp(UP,WP,VP)  !  SVD the square matrix UP.
            !
            ! option 2: LAPACK          
            ALLOCATE ( AP(NNULL,NNULL) )
            AP = UP
            CALL SGESVD( 'A', 'A', NNULL, NNULL, AP, NNULL, WP, UP,   &
                         NNULL, VP, NNULL, WORK, LWORK, IERS )
            IF (IERS.NE.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE(IUN,*) 'ERROR IN QA_DERIVWTS CALLING SGESVD(2)'
                  WRITE(IUN,*) 'IRET_SVD = ', IERS
                  WRITE(IUN,*) 'NNULL = ', NNULL
                  WRITE(IUN,*) 'AP:'
                  DO I=1,NNULL
                     WRITE(IUN,*) (AP(i,k),k=1,NNULL)
                  END DO
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            ENDIF
            !
            ! Adjust that solution by the combination of null vectors that
            ! zeroes all the selected components
            DO J=1,N+1
               I2 = 0
               DO I=1,M
                  IF ( ISNULL(I) ) THEN
                     I2 = I2 + 1
                     XP(I2) = WT(I,J)
                  END IF
               END DO
               WHERE (WP .GT. WMIN)
                  TMPP=MATMUL(XP,UP)/WP          ! Calculate diag(1/WP )UP T WT,
               ELSEWHERE
                  TMPP=0.0                    ! but replace 1/WPj by zero if WPj is small.
               END WHERE
               WT(:,J) = WT(:,J) - MATMUL(VF,MATMUL(VP,TMPP))
            END DO
         END IF
         IF(ITEST.GT.0) THEN
            write(ITEST,*) 'WT:'
            do i=1,m
               write(ITEST,*) (WT(i,k),k=1,n+1)
            end do
         END IF
      ELSE
!
!        Option 2: use Gaussian
!
         RECALC = .FALSE.
         P = M - NEQ
!
!        To complete the system we need P = M-N more equations,
!        that set P of the variables to zero
!
         IF ( P.GT.0 ) THEN
            ALLOCATE(IND(P))
            IND = 0
         END IF
         INCSTAT = 0
         DO WHILE (INCSTAT.EQ.0)
            IF ( P.EQ.0 ) THEN
               INCSTAT = 1
            ELSE
            ! Put ones in the array to set the selected indices to zero
               DO K=1,P
                  A(NEQ+K,NEQ+K-IND(K)) = 1.
               END DO
            END IF
            CALL QA_LINGAUSS(M, N+1, A, B, WT, DET, IST)
            IF(ITEST.GT.0) THEN
               write(ITEST,*) 'A:'
               do i=1,m
                  write(ITEST,*) (A(i,k),k=1,m)
               end do
               write(ITEST,*) 'B:'
               do i=1,m
                  write(ITEST,*) (B(i,k),k=1,n+1)
               end do
               write(ITEST,*) 'QA_LINGAUSS done, IST, DET =', IST, DET
               write(ITEST,*) 'WT:'
               do i=1,m
                  write(ITEST,*) (WT(i,k),k=1,n+1)
               end do
            END IF
            RECALC = IST.NE.0 .OR. ABS(DET).LE.DETMIN
            IF ( P.GT.0 ) THEN
               ! If this system has a satisfactory solution, exit
               IF ( .NOT.RECALC ) EXIT
               ! Otherwise, reset the selected indices and try the next set:
               DO K=1,P
                  A(NEQ+K,NEQ+K-IND(K)) = 0.
               END DO
               CALL QA_INCR( IND, N, INCSTAT )
            END IF 
         END DO
      END IF
!
!   6b. Check the solution.
!
      
      IF(ITEST.GT.0) THEN
         ALLOCATE( BCHK(M,N+1) )
         WRITE(ITEST,*) 'BCHK:'
         DO I=1,M
           DO K=1,N+1
             BCHK(I,K) = 0.
             DO J=1,M
               BCHK(I,K) = BCHK(I,K) + A(I,J)*WT(J,K)
             END DO
           END DO
           WRITE (ITEST,*) (BCHK(I,K), K=1,N+1)
         END DO
      END IF
!
!   7. If this failed, retry (recursively) with lower order
!
      IF ( RECALC ) THEN
        QMAX = QMAX - 1
        IF ( QMAX.GT.0 ) THEN
          CALL QA_DERIVWTS( NM, N, X, Y, LAMBDA, QMAX, WDX, WDY,      &
                         WDXX, WDXY, WDYY, WDXXX, WDXXY, WDXYY,       &
                         WDYYY, DET, ierr=IERS, ndse=IUN )
          IF (IERS.NE.0 ) THEN
            IF ( IUN.GT.0 ) THEN
               WRITE(IUN,*) 'ERROR IN QA_DERIVWTS CALLING ITSELF'
               WRITE(IUN,*) 'IERR = ', IERS
            END IF
            IF ( PRESENT(IERR) ) IERR = IERS
          END IF
        END IF
        RETURN
      END IF
!   
!   8. Transfer results to derivative weight arrays:
!
      DO K=0,N
        KP = K + 1
        DO I=1,M
           IF(I.EQ.1) THEN
              WDX(K) = WT(I,KP)
           ELSE IF(I.EQ.2) THEN
              WDY(K) = WT(I,KP)
           ELSE IF(I.EQ.3) THEN
              WDXX(K) = 2.*WT(I,KP)
           ELSE IF(I.EQ.4) THEN
              WDYY(K) = 2.*WT(I,KP)
           ELSE IF(I.EQ.5) THEN
              WDXY(K) = WT(I,KP)
           ELSE IF(I.EQ.6) THEN
              WDXXX(K) = 6.*WT(I,KP)
           ELSE IF(I.EQ.7) THEN
              WDYYY(K) = 6.*WT(I,KP)
           ELSE IF(I.EQ.8) THEN
              WDXXY(K) = 2.*WT(I,KP)
           ELSE IF(I.EQ.9) THEN
              WDXYY(K) = 2.*WT(I,KP)
           END IF
        END DO
      END DO
      !write(*,*) 'Results transferred'
      RETURN
!
!/
!/ End of QA_DERIVWTS ------------------------------------------------ /
!/
      END SUBROUTINE QA_DERIVWTS
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_DERWTS( QWTS )
!/
!       Richard Gorman, NIWA
!         November, 2009:   Origination.
!         Jan-April 2014:   Rewrite using QWTS structure
!
!  1. Purpose :
!
!      Creates arrays of precomputed interpolation weights for all
!      possible arrangements of quadtree neighbours
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     QWTS   QA_WEIGHTS1 I/O  Structure containing the following fields
!       MCASE    Int.   I   Maximum number of configurations of neighbours
!                           around each cell in a quadtree structure
!       NNGBR    Int.   I   Maximum number of neighbouring cells (8 for quadtree)
!       IQUADDIR I.A    O   Code for type of neighbour(s) in
!                           each direction, for each possible case IC
!                           IQUADDIR(IC,P)
!                            = 1 for a neighbour of equal level
!                   		 = 2 for two neighbours of higher level
!                   		 = 3 for a neighbour of lower level,
!                                displaced in -ve direction
!                   		     (S for E,W nbrs, W for S,N nbrs)
!                   		 = 4 for a neighbour of lower level,
!                                displaced in +ve direction
!                   		     (N for E,W nbrs, E for S,N nbrs)
!                   		 = 5 for one neighbour of higher level,
!                                displaced in -ve direction
!                   		     (S for E,W nbrs, W for S,N nbrs)
!                   		 = 6 for one neighbour of higher level,
!                                displaced in +ve direction
!                   		     (N for E,W nbrs, E for S,N nbrs)
!                   		 = 7 for no neighbour
!                          in the [W,E,S,N] direction for P = [1,2,3,4]
!       MAXORDER I.A   O   Maximum order of derivative computed for each case.
!       DERWTC   R.A   O   Weights for derivatives at cell centre:
!                          DERWTC(IC,K,1) = wt for ngbr K in d/dX for case IC,
!                          DERWTC(IC,K,2) = wt for ngbr K in d/dY for case IC,
!                          DERWTC(IC,K,3) = wt for ngbr K in d2/dX2 for case IC,
!                          DERWTC(IC,K,4) = wt for ngbr K in d2/dY2 for case IC,
!                          DERWTC(IC,K,5) = wt for ngbr K in d2/dXdY for case IC
!       XYREL    R.A   O   X- and Y- coordinates of target and neighbour 
!                            cell centres (relative)
!       DETERM   R.A   O   Determinant of matrix used to derive weights,
!                            for each case.
!       LAMBDA0  R.A.  O   Natural neighbour coordinates relative 
!                            to cell centre
!
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     QA_DERIVWTS Subr. qa_utils Compute weights for higher order derivatives
!                               on an extended stencil 
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_WEIGHTS1), INTENT(INOUT) :: QWTS
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER      :: MCASE=2401, NNGBR=8
      INTEGER   :: ICASE, IW, IE, IS, IN, NCELL, QMAX, II, NCASE0
      REAL      :: LAM(NNGBR), INFVAL, DET
      REAL      :: X(NNGBR), Y(NNGBR), XC(NNGBR), YC(NNGBR)
      REAL      :: WDX(0:NNGBR), WDY(0:NNGBR) 
      REAL      :: WDXX(0:NNGBR), WDYY(0:NNGBR), WDXY(0:NNGBR) 
      REAL      :: WDXXX(0:NNGBR), WDYYY(0:NNGBR), WDXXY(0:NNGBR),    &
                   WDXYY(0:NNGBR)
      LOGICAL   :: ISCELL(NNGBR)
 
      NCASE0 = 7
      INFVAL = -9.
      ICASE = 0
      ! Default:
      LAM = 0.25
      QWTS%DERWTC = 0.
! Loop over possible configurations of western neighbours:
      DO IW=1,NCASE0
         ! Possible configurations of western neighbours:
         X(1) = INFVAL
         Y(1) = INFVAL
         X(5) = INFVAL
         Y(5) = INFVAL
         ISCELL(1) = .FALSE.
         ISCELL(5) = .FALSE.
         IF ( IW.EQ.1 ) THEN
         !  1: one cell of equal level 
            X(1) = -1.
            Y(1) = 0.
            ISCELL(1) = .TRUE.
         ELSE IF ( IW.EQ.2 ) THEN
         !  2: two cells of higher level
            X(1) = -0.75
            Y(1) = -0.25
            X(5) = -0.75
            Y(5) = 0.25
            ISCELL(1) = .TRUE.
            ISCELL(5) = .TRUE.
         ELSE IF ( IW.EQ.3 ) THEN
         ! 3: one cell of lower level, offset sMAPSTPouthward
            X(1) = -1.5
            Y(1) = -0.5
            ISCELL(1) = .TRUE.
         ELSE IF ( IW.EQ.4 ) THEN
         ! 4: one cell of lower level, offset northward
            X(1) = -1.5
            Y(1) = 0.5
            ISCELL(1) = .TRUE.
         ELSE IF ( IW.EQ.5 ) THEN
         ! 5: one cell of higher level, offset southward
            X(1) = -0.75
            Y(1) = -0.25
            ISCELL(1) = .TRUE.
         ELSE IF ( IW.EQ.6 ) THEN
         ! 6: one cell of higher level, offset northward
            X(5) = -0.75
            Y(5) = 0.25
            ISCELL(5) = .TRUE.
         END IF
         DO IE=1,NCASE0
            X(2) = INFVAL
            Y(2) = INFVAL
            X(6) = INFVAL
            Y(6) = INFVAL
            ISCELL(2) = .FALSE.
            ISCELL(6) = .FALSE.
            ! Possible eastern neighbours:
            IF ( IE.EQ.1 ) THEN
               X(2) = 1.
               Y(2) = 0.
               ISCELL(2) = .TRUE.
            ELSE IF ( IE.EQ.2 ) THEN
               X(2) = 0.75
               Y(2) = -0.25
               X(6) = 0.75
               Y(6) = 0.25
               ISCELL(2) = .TRUE.
               ISCELL(6) = .TRUE.
            ELSE IF ( IE.EQ.3 ) THEN
               X(2) = 1.5
               Y(2) = -0.5
               ISCELL(2) = .TRUE.
            ELSE IF ( IE.EQ.4 ) THEN
               X(2) = 1.5
               Y(2) = 0.5
               ISCELL(2) = .TRUE.
            ELSE IF ( IE.EQ.5 ) THEN
               X(2) = 0.75
               Y(2) = -0.25
               ISCELL(2) = .TRUE.
            ELSE IF ( IE.EQ.6 ) THEN
               X(6) = 0.75
               Y(6) = 0.25
               ISCELL(6) = .TRUE.
            END IF
            DO IS=1,NCASE0
               X(3) = INFVAL
               X(7) = INFVAL
               Y(3) = INFVAL
               Y(7) = INFVAL
               ISCELL(3) = .FALSE.
               ISCELL(7) = .FALSE.
               ! Possible southern neighbours:
               IF ( IS.EQ.1 ) THEN
                  X(3) = 0.
                  Y(3) = -1.
                  ISCELL(3) = .TRUE.
               ELSE IF ( IS.EQ.2 ) THEN
                  X(3) = -0.25
                  Y(3) = -0.75
                  X(7) =  0.25
                  Y(7) = -0.75
                  ISCELL(3) = .TRUE.
                  ISCELL(7) = .TRUE.
               ELSE IF ( IS.EQ.3 ) THEN
                  X(3) = -0.5
                  Y(3) = -1.5
                  ISCELL(3) = .TRUE.
               ELSE IF ( IS.EQ.4 ) THEN
                  X(3) = 0.5
                  Y(3) = -1.5
                  ISCELL(3) = .TRUE.
               ELSE IF ( IS.EQ.5 ) THEN
                  X(3) = -0.25
                  Y(3) = -0.75
                  ISCELL(3) = .TRUE.
               ELSE IF ( IS.EQ.6 ) THEN
                  X(7) = 0.25
                  Y(7) = -0.75;
                  ISCELL(7) = .TRUE.
               END IF
               DO IN=1,NCASE0
                  X(4) = INFVAL
                  X(8) = INFVAL
                  Y(4) = INFVAL
                  Y(8) = INFVAL
                  ISCELL(4) = .FALSE.
                  ISCELL(8) = .FALSE.
                  ! Possible northern neighbours:
                  IF ( IN.EQ.1 ) THEN
                     X(4) = 0.
                     Y(4) = 1.
                     ISCELL(4) = .TRUE.
                  ELSE IF ( IN.EQ.2 ) THEN
                     X(4) = -0.25
                     Y(4) = 0.75
                     X(8) = 0.25
                     Y(8) = 0.75
                     ISCELL(4) = .TRUE.
                     ISCELL(8) = .TRUE.
                  ELSE IF ( IN.EQ.3 ) THEN
                     X(4) = -0.5
                     Y(4) = 1.5
                     ISCELL(4) = .TRUE.
                  ELSE IF ( IN.EQ.4 ) THEN
                     X(4) = 0.5
                     Y(4) = 1.5
                     ISCELL(4) = .TRUE.
                  ELSE IF ( IN.EQ.5 ) THEN
                     X(4) = -0.25
                     Y(4) = 0.75
                     ISCELL(4) = .TRUE.
                  ELSE IF ( IN.EQ.6 ) THEN
                     X(8) = 0.25
                     Y(8) = 0.75
                     ISCELL(8) = .TRUE.
                  END IF
                    ICASE = ICASE + 1
                    QWTS%IQUADDIR(ICASE,1) = IW
                    QWTS%IQUADDIR(ICASE,2) = IE
                    QWTS%IQUADDIR(ICASE,3) = IS
                    QWTS%IQUADDIR(ICASE,4) = IN
                    NCELL = 0
                    DO II=1,8
                       QWTS%XYREL(ICASE,II,1) = X(II)
                       QWTS%XYREL(ICASE,II,2) = Y(II)
                       IF ( ISCELL(II) ) THEN
                          NCELL = NCELL + 1
                          XC(NCELL) = X(II)
                          YC(NCELL) = Y(II)
                       END IF
                    END DO
                    QMAX = 2
                    CALL QA_DERIVWTS( NNGBR, NCELL, XC, YC, LAM,      &
                                   QMAX, WDX, WDY, WDXX, WDYY, WDXY,  &
                                   WDXXX, WDYYY, WDXXY, WDXYY, DET )
                    QWTS%MAXORDER(ICASE) = QMAX
                    QWTS%DETERM(ICASE) = DET
                    IF ( QMAX.GE.1 ) THEN
                       QWTS%DERWTC(ICASE,0,1) = WDX(0)
                       QWTS%DERWTC(ICASE,0,2) = WDY(0)
                    END IF
                    IF ( QMAX.GE.2 ) THEN
                       QWTS%DERWTC(ICASE,0,3) = WDXX(0)
                       QWTS%DERWTC(ICASE,0,4) = WDYY(0)
                       QWTS%DERWTC(ICASE,0,5) = WDXY(0)
                    END IF
                    NCELL = 0
                    DO II=1,8
                       IF ( ISCELL(II) ) THEN
                          NCELL = NCELL + 1
                          IF ( QMAX.GE.1 ) THEN
                             QWTS%DERWTC(ICASE,II,1) = WDX(NCELL)
                             QWTS%DERWTC(ICASE,II,2) = WDY(NCELL)
                          END IF
                          IF ( QMAX.GE.2 ) THEN
                             QWTS%DERWTC(ICASE,II,3) = WDXX(NCELL)
                             QWTS%DERWTC(ICASE,II,4) = WDYY(NCELL)
                             QWTS%DERWTC(ICASE,II,5) = WDXY(NCELL)
                          END IF
                       END IF
                    END DO
               END DO
            END DO
         END DO
      END DO
!/
!/ ------------------------------------------------------------------- /
!/
      END SUBROUTINE QA_DERWTS
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_DISPREL( SIGND, KH, CGND, IERR, NDSE )
!/
!         Richard Gorman, NIWA, April, 2008
!
!  1. Purpose :
!
!      Calculates wave number from dispersion relation
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       KH      Real    I    depth multiplied by wave number
!       SIGND   Real    I/O  radian frequency times sqrt(depth/g)
!       CGND    Real    O    group velocity /sqrt(g*depth)
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_BMLSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC   Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
!
      IMPLICIT NONE
      REAL, INTENT(IN)    :: SIGND
      REAL, INTENT(INOUT) :: KH
      REAL, INTENT(OUT)   :: CGND
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!
!  Local variables:
!      
      REAL      ::  TOL, F, FDASH, A1, BB, MAXKH, MAXSIGND
      REAL      ::  KH0
      INTEGER   ::  NITER, ITER
      INTEGER   ::  IUN
!
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE      
      IF ( PRESENT(IERR) ) IERR = 0
      NITER = 20
      TOL = 0.00001
      MAXKH = 100.
      MAXSIGND = 2.
      KH0 = KH
      IF (SIGND.GT.MAXSIGND) THEN
         KH = SIGND*SIGND
         CGND = 0.5*SIGND/SQRT(KH)
         RETURN
      END IF
      DO ITER = 1,NITER
         F=KH*TANH(KH)-SIGND**2
         FDASH=KH/COSH(KH)**2+TANH(KH)
         A1=KH-F/FDASH
         BB=A1-KH
         IF (ABS(BB).LT.TOL) THEN
            CGND = 0.5*(1+2.*KH/SINH(2.*KH))*TANH(KH)*SIGND/SQRT(KH)
            RETURN
         END IF
         KH=A1
         IF (KH.GT.MAXKH) THEN
            KH = SIGND*SIGND
            CGND = 0.5*SIGND/SQRT(KH)
            RETURN
         END IF
      END DO
      IF ( IUN.GT.0 ) THEN
         WRITE(IUN,*) ' QA_DISPREL did not converge'
         WRITE(IUN,*) ' SIGND, KH0, KH = ', SIGND, KH0, KH
      END IF
      IF ( PRESENT(IERR) ) IERR = 1
      RETURN
      END SUBROUTINE QA_DISPREL
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_FILLSPARE( NSEA, QSPARE, ISFROM, ISTO, NMOVEG,    &
                               NSEA2, NMOVLOC, JFROM, JTO )
!/
!       Richard Gorman, NIWA
!         21-Jun-2012:      Origination.
!         Jan-April 2014:   Rewrite using QSPARE structure
!/
!  1. Purpose :
!
!     Renumber cells to eliminate spare cells 
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NSEA    Int.     I   Initial number of cells
!       QSPARE  QA_SPARE I/O structure with components:
!          NPART   Int.  Number of partitions
!          PARTNO  I.A.  Partition number for each cell
!          JSVAL   I.A.  Local cell index for each global index
!          NCLOC   I.A.  Actual number of cells in each partition
!          ISVAL   I.A.  Global cell index for each local index & partition
!          NSPARE  I.A.  Number of spare local cell indices
!          JSSPARE I.A.  List of spare local cell indices
!       ISFROM  I.A.   O   (global) Sea indices that have been removed
!       ISTO    I.A.   O   (global) Sea indices that replace them
!       NMOVEG  Int.   O   (Global) Number of cells that have been swapped
!       NSEA2   Int.   O   Final number of cells
!       NMOVLOC I.A.   O*  Number of cells that have been swapped
!       JFROM   I.A.   O*  (local) Sea indices that have been removed
!       JTO     I.A.   O*  (local) Sea indices that replace them
!     ----------------------------------------------------------------
!                       * = optional
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       None.
!
!  5. Called by :
!
!       WQADGR   Quadtree grid adaptation routine.
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NSEA
      TYPE(QA_SPARE), INTENT(INOUT)  :: QSPARE
      INTEGER, INTENT(OUT)    :: ISFROM(:), ISTO(:)
      INTEGER, INTENT(OUT)    :: NMOVEG, NSEA2
      INTEGER, OPTIONAL, INTENT(OUT)    :: NMOVLOC(:)
      INTEGER, OPTIONAL, INTENT(OUT)    :: JFROM(:,:), JTO(:,:)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, ALLOCATABLE  :: NMOVE(:)
      INTEGER, ALLOCATABLE  :: JSFROM(:,:), JSTO(:,:)
      INTEGER, ALLOCATABLE  :: ISSPARE(:)
      INTEGER     ::  IH, IL, IM, ISPH, ISPL, ISEA1, ISEA2, PV,       &
                      JSEA1, JSEA2, JSPH, JSPL, NSPAREG, MSPG, MSP,   &
                      MPART, NPART
      LOGICAL    :: DOPART
!
      NPART = QSPARE%NPART

      MSPG = SIZE(ISFROM,1)
      MSP = SIZE(QSPARE%JSSPARE,1)
      MPART = SIZE(QSPARE%JSSPARE,2)
      MPART = MAX(MPART,SIZE(QSPARE%NSPARE,1))
      DOPART =  NPART.GT.1
      IF ( .NOT. DOPART ) MPART = 1
      ALLOCATE ( NMOVE(MPART), JSFROM(MSP,MPART),                     &
                 JSTO(MSP,NPART), ISSPARE(MSPG) )
  !        
      NMOVE = 0
      NSPAREG = 0
      DO PV = 1,NPART
         IH = 1
         IL = QSPARE%NSPARE(PV)
         IMLOOP: DO IM=1,QSPARE%NSPARE(PV)
            ! highest spare seapoint index:
            JSPH = QSPARE%JSSPARE(IH,PV)
            ! lowest spare seapoint index:
            JSPL = QSPARE%JSSPARE(IL,PV)
            DO WHILE (QSPARE%NCLOC(PV).EQ.JSPH)
               QSPARE%NCLOC(PV) = QSPARE%NCLOC(PV) - 1
               IH = IH + 1
               IF (IH.GT.IL) THEN
                  EXIT IMLOOP
               END IF
               JSPH = QSPARE%JSSPARE(IH,PV)
            END DO
            JSEA1 = QSPARE%NCLOC(PV)
            ISEA1 = QSPARE%ISVAL(JSEA1,PV)
            ISEA2 = QSPARE%ISVAL(JSPL,PV)
            JSFROM(IM,PV) = JSEA1
            JSTO(IM,PV) = JSPL
            QSPARE%ISVAL(JSPL,PV) = ISEA1
            QSPARE%JSVAL(ISEA1) = JSPL
            NSPAREG = NSPAREG + 1
            NMOVE(PV) = NMOVE(PV) + 1
            QSPARE%NCLOC(PV) = QSPARE%NCLOC(PV) - 1
            ISSPARE(NSPAREG) = ISEA2 
            IL = IL - 1
         END DO IMLOOP
      END DO 
      NMOVEG = 0
      NSEA2 = NSEA
      IH = 1
      IL = NSPAREG
      IMGLOOP: DO IM=1,NSPAREG
         ! highest spare seapoint index:
         ISPH = ISSPARE(IH)
         ! lowest spare seapoint index:
         ISPL = ISSPARE(IL)
         DO WHILE (NSEA2.EQ.ISPH)
            NSEA2 = NSEA2 - 1
            IH = IH + 1
            IF (IH.GT.IL) THEN
               EXIT IMGLOOP
            END IF
            ISPH = ISSPARE(IH)
         END DO
         ISFROM(IM) = NSEA2
         ISTO(IM) = ISPL
         JSEA2 = QSPARE%JSVAL(NSEA2)
         PV = QSPARE%PARTNO(NSEA2)
         QSPARE%PARTNO(ISPL) = PV
         QSPARE%JSVAL(ISPL) = JSEA2
         QSPARE%ISVAL(JSEA2,PV) = ISPL
         NMOVEG = NMOVEG + 1
         NSEA2 = NSEA2 - 1
         IL = IL - 1
      END DO IMGLOOP
      IF (PRESENT(NMOVLOC)) NMOVLOC = NMOVE
      IF (PRESENT(JFROM)) JFROM = JSFROM
      IF (PRESENT(JTO)) JTO = JSTO

!/ ------------------------------------------------------------------- /
      END SUBROUTINE QA_FILLSPARE
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_FINDNBR ( QTREE, IERR, NDSE, NCALC, ICALC )
!/
!       Richard Gorman, NIWA
!         April, 2008:      Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!      Determine the neighbouring cells for a set of leaf points
!      given a threaded quadtree structure
!
!      TO DO/CHECK: "inverse neighbour"
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE   QA_TREE I/O  Quadtree structure, including the following
!                             components used in this subroutine:
!          NCELL   Int.  I   Number of cells
!          NQUAD   Int.  I   Number of quads
!          UNDEF_TYPE Int. I Value of CELL_TYPE representing undefined cells
!          INDQUAD I.A.  I   Index of the quad containing each input cell
!          INDSUB  I.A.  I   Index of the sub-quad containing each cell:
!                              0 = none (home cell is at level zero)
!                              1,2,3,4 = SW, SE, NW, NE corner of quad
!          QICELL  I.A.  I   Sea-point indices for the cells of the quad
!          QNBR    I.A.  I   Quad indices for the neighbours of each quad
!          QPARENT I.A.  I   Quad index for the parent of each quad
!          QCHILD  I.A.  I   Quad indices for the children of each quad
!          CELL_TYPE I.A. I  Flag for wet/dry/boundary/etc. cell
!          NGBR    I.A.  O   Sea-point indices for the neighbours of input cells
!                               NGBR(I,1:4) = index of primary (i.e. of equal 
!                               or lower level) W,E,S,N neighbours,
!                               respectively, if any, of cell I
!                               NGBR(I,5:8) = index of secondary (i.e. of
!                               higher level) W,E,S,N neighbours, respectively, 
!                               if any, of cell I
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!       NCALC    Int.   I*  Number of cells to process (default = all)
!       ICALC    Int.   I*  Array of cell indices to process
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_BMLSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC   Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(INOUT)    :: QTREE
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
      INTEGER, OPTIONAL, INTENT(IN)   :: NCALC
      INTEGER, OPTIONAL, INTENT(IN)   :: ICALC(:)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: ISEA, QNHOME, QSUBHOME, IDIR, QNC, LEVREL,   &
                         QISEA1, QISEA2, QPAR, II, QP, QP2, QQ, JJ
      INTEGER         :: IUN
      INTEGER         :: ISEEK, NSEEK
      INTEGER         :: QN(4), QS(4), QS1(4), QS2(4)
      INTEGER         :: IREV(4), IREVL(4,4)
      INTEGER           :: MSEA, MQUAD
      DATA QS1/ 2, 1, 3, 1 /
      DATA QS2/ 4, 3, 4, 2 /
      DATA IREV/ 2, 1, 4, 3 /    ! Opposite direction 
      DATA IREVL/ 2, 1, 4, 3,                                         &
                  2, 1, 8, 3,                                         &
                  6, 1, 4, 3,                                         &
                  2, 5, 4, 7 /   
!/
!/ ------------------------------------------------------------------- /
!
      IF ( PRESENT(IERR) ) IERR = 0
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
!/
      MQUAD = SIZE(QTREE%QICELL,1)
      MSEA = SIZE(QTREE%INDQUAD,1)
!
!   Determine which cells to process:
!   If ICALC and NCALC are specified, process the first NCALC indices
!   in the array ICALC.
!   If NCALC is provided but not ICALC, process the first NCALC cells
!   If NCALC is not provided, or is zero, process all cells
!
      IF ( PRESENT(NCALC) ) THEN
         IF ( NCALC.LE.0 ) THEN
            NSEEK = QTREE%NCELL
         ELSE
            NSEEK = MIN(NCALC, QTREE%NCELL)
         END IF
      ELSE
         NSEEK = QTREE%NCELL
      END IF
      IF ( PRESENT(ICALC) ) THEN
         NSEEK = MIN(SIZE(ICALC),NSEEK)
      END IF
!
! Loop over requested cells
!
      DO ISEEK=1, NSEEK
         IF ( PRESENT(ICALC) ) THEN
            ISEA = ICALC(ISEEK)
         ELSE
            ISEA = ISEEK
         END IF
!
!  Initialise neighbours of the home cell ISEA to nul:
!
         QTREE%NGBR(ISEA,:) = 0
!
!  No further processing if ISEA is an undefined cell:
!
         IF ( QTREE%CELL_TYPE(ISEA).EQ.QTREE%UNDEF_TYPE ) CYCLE
!
!  Find the quad number and subquad index for the home cell
         QNHOME = QTREE%INDQUAD(ISEA)
         QSUBHOME = QTREE%INDSUB(ISEA)
         IF (QNHOME.LT.1 .OR. QNHOME.GT.QTREE%NQUAD) THEN
            CYCLE
         END IF
         IF (QSUBHOME.LT.0 .OR. QSUBHOME.GT.4) THEN
            CYCLE
         END IF
!  Find the quad number QN and subquad index QS for the four equal-level
!  neighbours (W(1), E(2), S(3), N(4)) of the home cell, if they exist
         IF (QSUBHOME.EQ.0) THEN
!     0. Cell is a central leaf (at level zero, or retained)
            QN(1) = QTREE%QNBR(QNHOME,1)
            QS(1) = 0
            QN(2) = QTREE%QNBR(QNHOME,2)
            QS(2) = 0
            QN(3) = QTREE%QNBR(QNHOME,3)
            QS(3) = 0
            QN(4) = QTREE%QNBR(QNHOME,4)
            QS(4) = 0
            LEVREL = 0
         ELSE IF (QSUBHOME.EQ.1) THEN
!  ISEA is in the SW quadrant of QNHOME:
!  ... W leads to the SE quadrant of the W neighbour
            QN(1) = QTREE%QNBR(QNHOME,1)
            QS(1) = 2
!  ... E leads to the SE quadrant of the same quad
            QN(2) = QNHOME
            QS(2) = 2
!  ... S leads to the NW quadrant of the S neighbour
            QN(3) = QTREE%QNBR(QNHOME,3)
            QS(3) = 3
!  ... N leads to the NW quadrant of the same quad
            QN(4) = QNHOME
            QS(4) = 3
            LEVREL = -1
         ELSE IF(QSUBHOME.EQ.2) THEN
!  ISEA is in the SE quadrant of QNHOME:
!  ... W leads to the SW quadrant of the same quad
            QN(1) = QNHOME
            QS(1) = 1
!  ... E leads to the SW quadrant of the E neighbour
            QN(2) = QTREE%QNBR(QNHOME,2)
            QS(2) = 1
!  ... S leads to the NE quadrant of the S neighbour
            QN(3) = QTREE%QNBR(QNHOME,3)
            QS(3) = 4
!  ... N leads to the NE quadrant of the same quad
            QN(4) = QNHOME
            QS(4) = 4
            LEVREL = -1
         ELSE IF(QSUBHOME.EQ.3) THEN
!  ISEA is in the NW quadrant of QNHOME:
!  ... W leads to the NE quadrant of the W neighbour
            QN(1) = QTREE%QNBR(QNHOME,1)
            QS(1) = 4
!  ... E leads to the NE quadrant of the same quad
            QN(2) = QNHOME
            QS(2) = 4
!  ... S leads to the SW quadrant of the same quad
            QN(3) = QNHOME
            QS(3) = 1
!  ... N leads to the SW quadrant of the N neighbour
            QN(4) = QTREE%QNBR(QNHOME,4)
            QS(4) = 1
            LEVREL = -1
         ELSE IF(QSUBHOME.EQ.4) THEN
!  ISEA is in the NE quadrant of QNHOME:
!  ... W leads to the NW quadrant of the same quad
            QN(1) = QNHOME
            QS(1) = 3
!  ... E leads to the NW quadrant of the E neighbour
            QN(2) = QTREE%QNBR(QNHOME,2)
            QS(2) = 3
!  ... S leads to the SE quadrant of the same quad
            QN(3) = QNHOME
            QS(3) = 2
!  ... N leads to the SE quadrant of the N neighbour
            QN(4) = QTREE%QNBR(QNHOME,4)
            QS(4) = 2
            LEVREL = -1
         END IF
         !WRITE(*,*) 'ISEA, QNHOME, QSUBHOME, QN, QS '
         !WRITE(*,'(11I4)') ISEA, QNHOME, QSUBHOME, QN, QS 
!
!  Loop over W, E, S, N directions:INDSUB
         DO IDIR=1,4
            IF (QN(IDIR) .LT. 0 .OR. QN(IDIR) .GT. QTREE%NQUAD) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                  WRITE(IUN,*) 'INVALID QUAD NUMBER'
                  WRITE(IUN,*) 'ISEA, IDIR, QN(IDIR): ',              &
                                ISEA, IDIR, QN(IDIR)
               END IF
               IF ( PRESENT(IERR) ) IERR = 1
               RETURN
            ELSEIF (QN(IDIR) .EQ. 0) THEN
!  No equal-level quad in that direction: the neighbour could be a 
!  lower-level cell. Identify the parent (QP) of the home quad:
               QP = QTREE%QPARENT(QNHOME)
               IF (QP.EQ.0) THEN
               ! No parent: QNHOME must be a level-zero quad, with no
               ! neighbour because of a domain boundary
                  QTREE%NGBR(ISEA,IDIR) = 0
                  CYCLE
               END IF
!  Also identify the parent quad's neighbour in the target direction
!  (QP2):               
               QP2 = QTREE%QNBR(QP,IDIR)
!  Identify the neighbouring lower-level cell as QTREE%QICELL(QQ,JJ):
!  i.e. the JJth cell of quad QQ. QQ could be either the parent QP
!  or its neighbour QP2.
!  Identifying QQ and JJ requires us to find out which child quadrant 
!  of QP is QNHOME
!   
               IF (QNHOME.EQ.QTREE%QCHILD(QP,1)) THEN
               !  QNHOME is in the SW quadrant of QP:
                  IF (IDIR.EQ.1) THEN
               !  ... W leads to the SE quadrant of the neighbour
                     QQ = QP2
                     JJ = 2
                  ELSEIF (IDIR.EQ.2) THEN
               !  ... E leads to the SE quadrant of the same quad
                     QQ = QP
                     JJ = 2
                  ELSEIF (IDIR.EQ.3) THEN
               !  ... S leads to the NW quadrant of the neighbour
                     QQ = QP2
                     JJ = 3
                  ELSEIF (IDIR.EQ.4) THEN
               !  ... N leads to the NW quadrant of the same quad
                     QQ = QP
                     JJ = 3
                  END IF
               ELSEIF (QNHOME.EQ.QTREE%QCHILD(QP,2)) THEN
               !  QNHOME is in the SE quadrant of QP:
                  IF (IDIR.EQ.1) THEN
               !  ... W leads to the SW quadrant of the same quad
                     QQ = QP
                     JJ = 1
                  ELSEIF (IDIR.EQ.2) THEN
               !  ... E leads to the SW quadrant of the neighbour
                     QQ = QP2
                     JJ = 1
                  ELSEIF (IDIR.EQ.3) THEN
               !  ... S leads to the NE quadrant of the neighbour
                     QQ = QP2
                     JJ = 4
                  ELSEIF (IDIR.EQ.4) THEN
               !  ... N leads to the NE quadrant of the same quad
                     QQ = QP
                     JJ = 4
                  END IF
               ELSEIF (QNHOME.EQ.QTREE%QCHILD(QP,3)) THEN
               !  QNHOME is in the NW quadrant of QP:
                  IF (IDIR.EQ.1) THEN
               !  ... W leads to the NE quadrant of the neighbour
                     QQ = QP2
                     JJ = 4
                  ELSEIF (IDIR.EQ.2) THEN
               !  ... E leads to the NE quadrant of the same quad
                     QQ = QP
                     JJ = 4
                  ELSEIF (IDIR.EQ.3) THEN
               !  ... S leads to the SW quadrant of the same quad
                     QQ = QP
                     JJ = 1
                  ELSEIF (IDIR.EQ.4) THEN
               !  ... N leads to the SW quadrant of the neighbour
                     QQ = QP2
                     JJ = 1
                  END IF
               ELSEIF (QNHOME.EQ.QTREE%QCHILD(QP,4)) THEN
               !  QNHOME is in the NE quadrant of QP:
                  IF (IDIR.EQ.1) THEN
               !  ... W leads to the NW quadrant of the same quad
                     QQ = QP
                     JJ = 3
                  ELSEIF (IDIR.EQ.2) THEN
               !  ... E leads to the NW quadrant of the neighbourINDSUB
                     QQ = QP2
                     JJ = 3
                  ELSEIF (IDIR.EQ.3) THEN
               !  ... S leads to the SE quadrant of the same quad
                     QQ = QP
                     JJ = 2
                  ELSEIF (IDIR.EQ.4) THEN
               !  ... N leads to the SE quadrant of the neighbour
                     QQ = QP2
                     JJ = 2
                  END IF
               ELSE
               !  QNHOME is not listed as a child quad of QP: ERROR!
                  IF ( IUN.GT.0 ) THEN
                     WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                     WRITE(IUN,*) 'QHOME IS NOT A CHILD OF ITS PARENT'
                     WRITE(IUN,*) 'ISEA, QNHOME, QP, QCHILD(QP,:): '
                     WRITE(IUN,*) ISEA, QNHOME, QP, QTREE%QCHILD(QP,:)
                  END IF
                  IF ( PRESENT(IERR) ) IERR = 1
                  RETURN
               END IF
               IF (QQ.EQ.0) THEN
!  No target (neighbouring) quad: domain limit, so no neighbouring cells
!  in this direction
                  QTREE%NGBR(ISEA,IDIR) = 0
                  CYCLE
               END IF
               QISEA1 = QTREE%QICELL(QQ,JJ)
               QTREE%NGBR(ISEA,IDIR) = QISEA1
!
!  Inverse neighbour relationship:
               !IF ( QISEA1.GT.0 ) THEN
               !   II = IREVL(IDIR,QSUBHOME)
               !   QTREE%NGBR(QISEA1,II) = ISEA
               !END IF
            ELSE
!  Check if there is a cell neighbour at equal level
               QISEA1 = QTREE%QICELL(QN(IDIR),QS(IDIR))
               IF (QISEA1.LT.0 .OR. QISEA1.GT.MSEA) THEN
                  IF ( IUN.GT.0 ) THEN
                     WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                     WRITE(IUN,*) 'INVALID CELL'
                     WRITE(IUN,*) 'ISEA, IDIR, QN(IDIR), ',           &
                             'QS(IDIR), QICELL(QN(IDIR),QS(IDIR)): ', &
                              ISEA, IDIR, QN(IDIR), QS(IDIR), QISEA1
                  END IF
                  IF ( PRESENT(IERR) ) IERR = 1
                  RETURN
               ELSEIF (QSUBHOME.EQ.0) THEN
!  Home quad is a leaf cell, at level zero
                  IF (QISEA1 .GT. 0) THEN
!  Target quad is also a leaf cell, so is a neighbour at equal order
                     QTREE%NGBR(ISEA,IDIR) = QISEA1
!  Inverse neighbour:
                     !QTREE%NGBR(QISEA1,IREV(IDIR)) = ISEA
                  ELSE
!  Otherwise, target quad has cells: choose the adjacent two
                     QISEA1 = QTREE%QICELL(QN(IDIR),QS1(IDIR))
                     IF (QISEA1.LT.0 .OR. QISEA1.GT.MSEA) THEN
                        IF ( IUN.GT.0 ) THEN
                           WRITE(IUN,*) 'ERROR IN QA_FINDNBR: '
                           WRITE(IUN,*) 'INVALID CELL'
                           WRITE(IUN,*) 'ISEA, IDIR, QN, QS1, ',      &
                                        'QICELL(QN,QS1): ', ISEA,     &
                                     IDIR, QN(IDIR), QS1(IDIR), QISEA1
                        END IF
                        IF ( PRESENT(IERR) ) IERR = 1
                        RETURN
                     END IF
                     QISEA2 = QTREE%QICELL(QN(IDIR),QS2(IDIR))
                     IF (QISEA2.LT.0 .OR. QISEA2.GT.MSEA) THEN
                        IF ( IUN.GT.0 ) THEN
                           WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                           WRITE(IUN,*) 'INVALID CELL'
                           WRITE(IUN,*) 'ISEA, IDIR, QN, QS12, ',     &
                                        'QICELL(QN,QS2): ', ISEA,     &
                                     IDIR, QN(IDIR), QS2(IDIR), QISEA2
                        END IF
                        IF ( PRESENT(IERR) ) IERR = 1
                        RETURN
                     END IF
                     IF(QISEA1.GT.0 .AND. QISEA2.GT.0) THEN
                        QTREE%NGBR(ISEA,IDIR) = QISEA1
                        QTREE%NGBR(ISEA,IDIR+4) = QISEA2
!  Inverse neighbours:
                        !QTREE%NGBR(QISEA1,IREV(IDIR)) = ISEA
                        !QTREE%NGBR(QISEA2,IREV(IDIR)) = ISEA
                     END IF
                  END IF
               ELSEIF (QISEA1.GT.0) THEN
!  There is a neighbour of equal order
                  QTREE%NGBR(ISEA,IDIR) = QISEA1
!  Inverse neighbour:
                  !QTREE%NGBR(QISEA1,IREV(IDIR)) = ISEA
               ELSE
!  There is no neighbouring cell of equal order
!  Check if there is a child quad instead of a child cell
                  QNC = QTREE%QCHILD(QN(IDIR),QS(IDIR))
                  IF (QNC.LT.0 .OR. QNC.GT.QTREE%NQUAD) THEN
                     IF ( IUN.GT.0 ) THEN
                        WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                        WRITE(IUN,*) 'INVALID CHILD'
                        WRITE(IUN,*) 'ISEA, IDIR, QN(IDIR),',         &
                            ' QS(IDIR),QCHILD(QN(IDIR),QS(IDIR)): ',  &
                              ISEA, IDIR, QN(IDIR), QS(IDIR), QNC
                     END IF
                     IF ( PRESENT(IERR) ) IERR = 1
                     RETURN
                  ELSEIF (QNC.NE.0) THEN
!  Target cell has children:
                     QISEA1 = QTREE%QICELL(QNC,0)
                     IF (QISEA1.LT.0 .OR. QISEA1.GT.MSEA) THEN
                        IF ( IUN.GT.0 ) THEN
                           WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                           WRITE(IUN,*) 'INVALID CELL'
                           WRITE(IUN,*) 'ISEA, IDIR, QNC, ',          &
                                        'QICELL(QNC,0): ',ISEA, IDIR,  &
                                        QNC, QISEA1
                        END IF
                        IF ( PRESENT(IERR) ) IERR = 1
                        RETURN
                     ELSEIF (QISEA1 .GT. 0) THEN
!  Target cell has a child which is a leaf cell of equal order
                        QTREE%NGBR(ISEA,IDIR) = QISEA1
!  Inverse neighbour:
                        !QTREE%NGBR(QISEA1,IREV(IDIR)) = ISEA
                        CYCLE
                     ELSE
!  Target quad has children: choose the adjacent two
                        QISEA1 = QTREE%QICELL(QNC,QS1(IDIR))
                        IF (QISEA1.LT.0 .OR. QISEA1.GT.MSEA) THEN
                           IF ( IUN.GT.0 ) THEN
                              WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                              WRITE(IUN,*) 'INVALID CELL'
                              WRITE(IUN,*) 'ISEA, IDIR, QNC, QS1, ',  &
                                       'QICELL(QNC,QS1): ',ISEA,      &
                                        IDIR, QNC, QS1(IDIR), QISEA1
                           END IF
                           IF ( PRESENT(IERR) ) IERR = 1
                           RETURN
                        END IF
                        QISEA2 = QTREE%QICELL(QNC,QS2(IDIR))
                        IF (QISEA2.LT.0 .OR. QISEA2.GT.MSEA) THEN
                           IF ( IUN.GT.0 ) THEN
                              WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                              WRITE(IUN,*) 'INVALID CELL'
                              WRITE(IUN,*) 'ISEA, IDIR, QNC, QS2, ',  &
                                       'QICELL(QNC,QS2): ',ISEA,      &
                                        IDIR, QNC, QS2(IDIR), QISEA2
                           END IF
                           IF ( PRESENT(IERR) ) IERR = 1
                           RETURN
                        END IF
                        IF(QISEA1.GT.0 .AND. QISEA2.GT.0) THEN
                           QTREE%NGBR(ISEA,IDIR) = QISEA1
                           QTREE%NGBR(ISEA,IDIR+4) = QISEA2
!  Inverse neighbours:
                           !QTREE%NGBR(QISEA1,IREV(IDIR)) = ISEA
                           !QTREE%NGBR(QISEA2,IREV(IDIR)) = ISEA
                           CYCLE
                        END IF
                     END IF
                  ELSE
!
! No corresponding child quad. Check if the neighbouring quad has other
! child cells/quads at the same level as the home cell, in which case what
! would be an equal neighbour is absent because it is dry, etc. In that 
! case, return no neighbours in this direction
                     IF ( ANY(QTREE%QICELL(QN(IDIR),1:4).GT.0) ) CYCLE
                     IF ( ANY(QTREE%QCHILD(QN(IDIR),1:4).GT.0) ) CYCLE
                  END IF
!
! No neighbours at equal or higher level: look for one at lower level
                  QPAR = QTREE%QPARENT(QN(IDIR))
                  IF (QPAR.LT.0 .OR. QPAR.GT.QTREE%NQUAD) THEN
                     IF ( IUN.GT.0 ) THEN
                        WRITE(IUN,*) 'ERROR IN QA_FINDNBR:'
                        WRITE(IUN,*) 'INVALID PARENT QUAD'
                        WRITE(IUN,*) 'ISEA, IDIR, QN(IDIR), ',        &
                                'QPARENT(QN(IDIR)): ',                &
                              ISEA, IDIR, QN(IDIR), QPAR
                     END IF
                     IF ( PRESENT(IERR) ) IERR = 1
                     RETURN
                  ELSEIF (QPAR.GT.0) THEN
                     DO II=1,4
                        IF(QTREE%QCHILD(QPAR,II).EQ.QN(IDIR)) THEN
                           QISEA1 = QTREE%QICELL(QPAR,II)
                           QTREE%NGBR(ISEA,IDIR) = QISEA1
!
!  Inverse neighbour relationship:
                           !IF ( QISEA1.GT.0 )                          &
                           !   QTREE%NGBR(QISEA1,IREVL(IDIR,QSUBHOME))  &
                           !       = ISEA
                           EXIT
                        END IF
                     END DO
                  ELSE
                  !  No parent: neighbour may be a level-0 leaf cell
                     QISEA1 = QTREE%QICELL(QN(IDIR),0)
                     IF(QISEA1.NE.0) THEN
                        QTREE%NGBR(ISEA,IDIR) = QISEA1
!
!  Inverse neighbour relationship:
                        !QTREE%NGBR(QISEA1,IREVL(IDIR,QSUBHOME)) = ISEA
                     END IF
                  END IF
               END IF
            END IF
         END DO
      END DO
      RETURN
!/
!/ End of FINDNBR ----------------------------------------------------- /
!/
      END SUBROUTINE QA_FINDNBR
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_INCR( IND, BASE, ISTAT )
!/
!/    28-Sep-2012 : Origination. Richard Gorman, NIWA
!/
!  1. Purpose :
!
!     Increment a multiple-digit counter in arbitrary base so that
!       IND(1) + IND(2)*BASE + ... + IND(P)*BASE**P
!     inreases by 1 until the condition
!       BASE>IND(1)>=IND(2)>=...>=IND(P)>=0
!     is (again) met
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IND      I.A.   I   Vector of indices, with values 0,...,BASE-1
!                           Equivalent to digits of an integer
!       BASE     Int.   I   Base for arithmetic
!       ISTAT    Int.   O   Return status, = 1 if the maximuminteger 
!                           representable
!
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_DERIVWTS Subr. qa_utils Compute weights for higher order
!                                derivatives on an extended stencil 
!     ---------------------------------------------------------------- 
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(INOUT)  :: IND(:)
      INTEGER, INTENT(IN)     :: BASE
      INTEGER, INTENT(OUT)    :: ISTAT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER     :: P, K, J, ITARG, N
      P = SIZE(IND)
      ISTAT = 0
      DO K=1,P
         IND(K) = IND(K) + 1
         IF ( IND(K) .LT. BASE ) RETURN
         IF ( K.EQ.P ) THEN
            ISTAT = 1
            RETURN
         END IF
         DO J=1,K
            DO N=1,P-J
               ITARG = IND(J+N)+1
               IF ( ITARG.LT.BASE) EXIT
            END DO
            IND(J) = ITARG
         END DO
      END DO    
      END SUBROUTINE QA_INCR
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_INTERP( QTREE, XTARG, YTARG, VAL, VAL_OUT,        &
                            NGORDER, WEIGHT, IERR, NDSE, SCALEFAC,    &
                            MSOL )
!/
!       Richard Gorman, NIWA
!         14-Oct-2010:      Origination.
!         Jan-April 2014:   Rewrite using QTREE structure
!/
!  1. Purpose :
!
!     Interpolate to a target location using weighted values from 
!     neighbouring cells. Up to 3rd order derivatives may be used
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!
!       QTREE   QA_TREE I  Quadtree structure, of which the following 
!                          components are affected:
!          LVLREF  Int.  Refinement level of reference grid
!          NCASE   I.A.  Index identifying the relative configuration of 
!                        neighbouring cells (out of 2401 possibilities)
!          INDWT   I.A.  Index, for each cell, into the lookup table of
!                        (second order) weights
!          GNBR    I.A.  Grand neighbours of each cell
!          NGNBR   I.A.  Number of grand neighbours of each cell
!       XTARG    Real   I   Target X coordinate
!       YTARG    Real   I   Target Y coordinate
!       VAL      R.A.   I   Data values at all cells
!       VAL_OUT  Real   O   Interpolated value at the target location
!       NGORDER  Int.   I/O Order of neighbours to use in interpolation
!                             = 0 to do nearest-neighbour interpolation,
!                                 i.e. take values from the cell in which 
!                                 XTARG,YTARG lies
!                             = 1 to use first order (i.e. immediate) ngbrs
!                             = 2 to use second order neighbours
!       WEIGHT   R.A    I*  Array of weights 
!                            WEIGHT(ITBL,NB,ISOL(K)) = weight to compute 
!                            DERIV(ISOL(K)) at cell ISEA from values at cell 
!                            ISTEN(NB), or ISEA if NB=0.
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!       SCALEFAC R.A.   I*   Scaling factors for extrapolation
!       MSOL     Int.   I/O*   Number of derivatives to be calculated 
!                           (up to 5 for NGORDER=1, up to 9 for NGORDER=2)
!
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     QA_XY2CELL Subr. qa_utils Identify the quadtree cell containing a
!                               given point
!     QA_DERIVC  Subr. qa_utils Compute up to 3rd order derivatives using
!                               values from neighbouring cells.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(IN)        :: QTREE
      REAL, INTENT(IN)                 :: XTARG, YTARG
      REAL, INTENT(IN)                 :: VAL(:)
      REAL, INTENT(OUT)                :: VAL_OUT
      INTEGER, INTENT(INOUT)           :: NGORDER
      REAL, OPTIONAL, INTENT(IN)       :: WEIGHT(1:,0:,1:)
      INTEGER, OPTIONAL, INTENT(OUT)   :: IERR
      INTEGER, OPTIONAL, INTENT(IN)    :: NDSE
      REAL, OPTIONAL, INTENT(IN)       :: SCALEFAC(2,2)
      INTEGER, OPTIONAL, INTENT(INOUT) :: MSOL
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: MSTEN
      INTEGER                 :: LVLREF
      REAL                    :: DELXIN, DELYIN, DELX, DELY
      REAL                    :: XCELL, YCELL, XFAC, YFAC
      REAL, ALLOCATABLE       :: DERIV(:)
      INTEGER, ALLOCATABLE    :: ISOL(:)
      INTEGER                 :: IS, IDER, ISEA, ITBL, NSTEN, LVCELL, &
                                 NSOL, MWT1, MWT2L, MWT2U, MWT3,      &
                                 IQUAD, ISUB
      INTEGER                 :: ISTAT(2)
      INTEGER, ALLOCATABLE    :: ISTEN(:)
      REAL                    :: SFAC(2,2)
      INTEGER                 :: IUN, IERS
!
!   Error handling:
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0
      IF ( QTREE%IWTORDER.LT.NGORDER ) THEN
         IF ( IUN.GT.0 )  THEN
            WRITE(IUN,*) 'ERROR IN QA_INTERP: IWTORDER < NGORDER'
            WRITE(IUN,*) 'IWTORDER: ', QTREE%IWTORDER
            WRITE(IUN,*) 'NGORDER: ', NGORDER
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN 
      END IF
!
      LVLREF = QTREE%LVLREF
      
      IF ( NGORDER.GT.0 .AND. PRESENT(WEIGHT) ) THEN
         MWT1 = SIZE(WEIGHT,1)
         MWT2L = LBOUND(WEIGHT,2)
         MWT2U = UBOUND(WEIGHT,2)
         MWT3 = SIZE(WEIGHT,3)
         IF ( NGORDER.GT.2 ) NGORDER = 2
         IF ( NGORDER.EQ.1 ) THEN
            ! Interpolation using first order neighbours:
            ! can calculate linear (MSOL=2) and 2nd order (MSOL=5)
            ! derivatives: 
            IF ( PRESENT(MSOL) ) THEN
               NSOL = MIN(MSOL,5)
            ELSE
               ! default to linear:
               !NSOL = 2
               ! default to second order derivatives:
               NSOL = 5
            END IF
            ! Check on size of weight array:
            IF ( MWT1.LT.2401 .OR. MWT2L.GT.0 .OR. MWT2U.LT.8 ) THEN
               IF (IUN.GT.0 ) THEN
                  WRITE(IUN,*) ' ERROR IN QA_INTERP:'
                  WRITE(IUN,*) ' 1ST ORDER WEIGHT TABLE IS TOO SMALL'
                  WRITE(IUN,*) '   SIZE(WEIGHT,1) = ', MWT1
                  WRITE(IUN,*) ' LBOUND(WEIGHT,2) = ', MWT2L
                  WRITE(IUN,*) ' UBOUND(WEIGHT,2) = ', MWT2U
                  WRITE(IUN,*) '   SIZE(WEIGHT,3) = ', MWT3
                  WRITE(IUN,*) ' SHOULD BE (2401,0:8,5)'
               END IF
               IF ( PRESENT(IERR) ) IERR = 1
               RETURN
            END IF
            MSTEN = SIZE(QTREE%NGBR,2)
         ELSEIF ( NGORDER.EQ.2 ) THEN
            IF ( PRESENT(MSOL) ) THEN
               NSOL = MIN(MSOL,9)
            ELSE
               ! default to linear:
               !NSOL = 2
               ! default to second order derivatives:
               !NSOL = 5
               ! default to third order derivatives:
               NSOL = 9
            END IF
            ! Check on size of weight array:
            IF ( MWT2L.GT.0 .OR. MWT2U.LT.36 ) THEN
               IF (IUN.GT.0 ) THEN
                  WRITE(IUN,*) ' ERROR IN QA_INTERP:'
                  WRITE(IUN,*) ' 2ND ORDER WEIGHT TABLE IS TOO SMALL'
                  WRITE(IUN,*) '   SIZE(WEIGHT,1) = ', MWT1
                  WRITE(IUN,*) ' LBOUND(WEIGHT,2) = ', MWT2L
                  WRITE(IUN,*) ' UBOUND(WEIGHT,2) = ', MWT2U
                  WRITE(IUN,*) '   SIZE(WEIGHT,3) = ', MWT3
                  WRITE(IUN,*) ' SHOULD BE (????,0:36,9)'
               END IF
               IF ( PRESENT(IERR) ) IERR = 1
               RETURN
            END IF
            MSTEN = SIZE(QTREE%GNBR,2)
         END IF
         ! can't solve for more derivatives than we have weights for:
         NSOL = MIN(NSOL,MWT3)
         
         ALLOCATE ( ISTEN(MSTEN) )
         ALLOCATE ( DERIV(NSOL) )
         ALLOCATE ( ISOL(NSOL) )
      ELSE
         NSOL = 0
      END IF
      IF ( PRESENT(MSOL) ) MSOL = NSOL
      IF ( PRESENT(SCALEFAC) ) THEN
         SFAC = SCALEFAC
      ELSE
         ! Defaults: 
         !   no extrapolation, but set to zero:
         SFAC = 0.
         !   extrapolation in all directions by applying last cell:
         !SFAC = 1.
         !   extrapolation for spectra grid:
         !     in theta (X) direction, apply last cell 
         !     (shouldn't be needed if cyclic):
         !SFAC(1,1) = 1.
         !SFAC(1,2) = 1.
         !     in frequency (Y) direction, set to zero below lowest bin
         !     power law above highest bin:
         !SFAC(2,1) = 0.
         !SFAC(2,2) = FAC**(-PWR)
         !
      END IF
!
!  Find the cell (if any) containing (XTARG, YTARG):
!
      CALL QA_XY2CELL ( QTREE, XTARG, YTARG, ISEA, XCELL, YCELL,      &
                        LVCELL, IQUAD, ISUB, ISTAT )
      IF ( ISTAT(1).LT.0 ) THEN
      !  (XTARG,YTARG) is in an internal part of the domain that is not
      !  in a cell (e.g. dry)
          VAL_OUT = 0.
          RETURN
      END IF
      ! 
      ! Offsets from the cell centre:
      DELX = XTARG - XCELL
      DELY = YTARG - YCELL
      ! 
      ! If (XTARG, YTARG) is outside the domain, ISEA is the nearest cell
      ! on the domain boundary. In that case, interpolate in the 'along-boundary'
      ! direction, but extrapolate by a power law in the direction normal
      ! to the boundary: SCALEFAC(i,j) is the scaling factor per unit increase 
      ! in X (i=1) or Y (i=2) beyond the lower (j=1) or upper (j=2) bound 
      !
      XFAC = 1.
      YFAC = 1.
      IF ( ISTAT(1).EQ.1 ) THEN
         XFAC = SFAC(1,1)**DELX 
         DELX = 0.
      ELSEIF ( ISTAT(1).EQ.2 ) THEN
         XFAC = SFAC(1,2)**DELX 
         DELX = 0.
      END IF
      IF ( ISTAT(2).EQ.1 ) THEN
         YFAC = SFAC(2,1)**DELY 
         DELY = 0.
      ELSEIF ( ISTAT(1).EQ.2 ) THEN
         YFAC = SFAC(2,2)**DELY 
         DELY = 0.
      END IF
      !
      ! Value at the cell centre:
      !
      VAL_OUT = VAL(ISEA)
      !  
      !  If the necessary weights and indices are provided to estimate derivatives, 
      !  improve this  estimate with a Taylor series
      !
      IF ( NSOL.GT.0 ) THEN
         DELXIN = 2.**(LVCELL-LVLREF)
         DELYIN = DELXIN
         !
         DO IS=1,NSOL
            ISOL (IS) = IS
         END DO
         IF ( NGORDER.EQ.1 ) THEN
         !
         ! First-order neighbours, linear interpolation:
            NSTEN = MSTEN
            ISTEN = QTREE%NGBR(ISEA,:)
            ITBL = QTREE%NCASE(ISEA)
         ELSE
         !
         ! Second-order neighbours:
            NSTEN = QTREE%NGNBR(ISEA)
            ISTEN = QTREE%GNBR(ISEA,:)
            ITBL = QTREE%INDWT(ISEA)
         END IF
         !
         ! Compute derivatives in the centre of the containing cell:
         CALL QA_DERIVC( ISEA, ITBL, NSOL, ISOL, NSTEN, ISTEN, VAL,   &
                         WEIGHT, DERIV, delxin=DELXIN, delyin=DELYIN, &
                         ierr=IERS, ndse=IUN )
         IF ( IERS.GT.0 ) THEN
            IF ( IUN.GT.0 ) THEN
               WRITE(IUN,*) 'ERROR IN QA_INTERP CALLING QA_DERIVC'
            END IF
            IF ( PRESENT(IERR) ) IERR = IERS
            RETURN
         END IF
         DO IS=1,NSOL
            IDER = ISOL(IS)
            IF (IDER.EQ.1) THEN
               VAL_OUT = VAL_OUT + DERIV(IS)*DELX
            ELSEIF (IDER.EQ.2) THEN
               VAL_OUT = VAL_OUT + DERIV(IS)*DELY
            ELSEIF (IDER.EQ.3) THEN
               VAL_OUT = VAL_OUT + 0.5*DERIV(IS)*DELX*DELX
            ELSEIF (IDER.EQ.4) THEN
               VAL_OUT = VAL_OUT + 0.5*DERIV(IS)*DELY*DELY
            ELSEIF (IDER.EQ.5) THEN
               VAL_OUT = VAL_OUT + DERIV(IS)*DELX*DELY
            ELSEIF (IDER.EQ.6) THEN
               VAL_OUT = VAL_OUT + (DERIV(IS)*DELX**3)/6.
            ELSEIF (IDER.EQ.7) THEN
               VAL_OUT = VAL_OUT + (DERIV(IS)*DELY**3)/6.
            ELSEIF (IDER.EQ.8) THEN
               VAL_OUT = VAL_OUT + 0.5*DERIV(IS)*DELX*DELX*DELY
            ELSEIF (IDER.EQ.9) THEN
               VAL_OUT = VAL_OUT + 0.5*DERIV(IS)*DELX*DELY*DELY
            END IF
         END DO
      END IF
      !
      !  Apply power-law extrapolation if required:
      !
      IF ( ISTAT(1).GT.0 ) VAL_OUT = XFAC*VAL_OUT
      IF ( ISTAT(2).GT.0 ) VAL_OUT = YFAC*VAL_OUT
!/
!/ End of QA_INTERP ----------------------------------------------------- /
!/
      END SUBROUTINE QA_INTERP
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_IOQT( IUN, QTREE, TASK, IERR, NDSE, NREC, NSIZE )
!/
!       Richard Gorman, NIWA
!         19-April-2014:      Origination.
!         28-April-2014:      Extension to direct access
!/
!  1. Purpose :
!
!      IO of quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!
!       IUN      Int.    I   Unit number
!       QTREE   QA_TREE I/O  Quadtree structure
!       TASK     Int.    I  
!                   TASK = -1 or -10 to READ header, sequential binary
!                        = +1 or +10 to WRITE   ""  ""
!                        = -1 or -11 to READ quad & cell arrays, sequential binary
!                        = +1 or +11 to WRITE   ""  ""
!                        = -2 or -20 to READ header, direct access binary
!                        = +2 to +20 to WRITE   ""  ""
!                        = -2 or -21 to READ quad & cell arrays, direct access binary
!                        = +2 to +21 to WRITE   ""  ""
!                        = -3 or -30 to READ quad header, formatted
!                        = +3 or +30 to WRITE   ""  ""
!                        = -3 or -31 to READ quad array data, formatted
!                        = +3 or +31 to WRITE   ""  ""
!                        = -4 or -40 to READ cell header, formatted
!                        = +4 or +40 to WRITE   ""  ""
!                        = -4 or -41 to READ cell array data, formatted
!                        = +4 or +41 to WRITE   ""  ""
!                        = -5 or -50 to READ weights header, formatted
!                        = +5 or +50 to WRITE   ""  ""
!                        = -5 or -51 to READ weights array data, formatted
!                        = +5 or +51 to WRITE   ""  ""
!                        = -6 or -60 to READ header, stream access binary
!                        = +6 to +60 to WRITE   ""  ""
!                        = -6 or -61 to READ quad & cell arrays, stream access binary
!                        = +6 to +61 to WRITE   ""  ""
!       IERR     Int.    O   Return flag = 1 for error, else 0
!       NDSE     Int.    I*  Unit number for error output (if >0)
!       NREC     Int.   I/O* Last record number for direct access 
!       NSIZE    Int.    I*  Number of integers or reals per record
!
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)              :: IUN
      TYPE(QA_TREE), INTENT(INOUT)     :: QTREE
      INTEGER, INTENT(IN)              :: TASK
      INTEGER, INTENT(OUT)             :: IERR
      INTEGER, OPTIONAL, INTENT(IN)    :: NDSE
      INTEGER, OPTIONAL, INTENT(INOUT) :: NREC
      INTEGER, OPTIONAL, INTENT(IN)    :: NSIZE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER    ::  LRB=4
      INTEGER    ::  NDE, MQ, MC, MC2, MC3, MC4, NQ, NC, NR, IC, IR,  &
                     IPART, NPART, NHD, NCOL, IHD, IR2, ATASK
      LOGICAL    ::  DOWTS, DOQA, DOCA, STREAM
      INTEGER    ::  LRECL, RPOS
      REAL(KIND=LRB), ALLOCATABLE :: WRITEBUFF(:)
!
!  Default for optional parameters
!
      NDE = 6
      IF ( PRESENT(NDSE) ) NDE = NDSE
      NR = -1
      IF ( PRESENT(NREC) ) NR = NREC
!
      IF ( PRESENT(NSIZE) ) THEN
         LRECL = LRB*NSIZE
         ALLOCATE ( WRITEBUFF(NSIZE) )
      END IF
!
      IERR = 0
!
!  For the time being, hardwire NOT to include weights tables
!  in binary I/O:
!
      ATASK = ABS(TASK)
!  Tasks which do I/O on weights:
      DOWTS = ATASK.EQ.5 .OR. ATASK.EQ.50 .OR. ATASK.EQ.51  
!  Tasks which do I/O on quad arrays:
      DOQA =  ATASK.EQ.1 .OR. ATASK.EQ.2 .OR. ATASK.EQ.3 .OR.         &
              ATASK.EQ.11 .OR. ATASK.EQ.21 .OR. ATASK.EQ.31
!  Tasks which do I/O on cell arrays:
      DOCA =  ATASK.EQ.1 .OR. ATASK.EQ.2 .OR. ATASK.EQ.4 .OR.         &
              ATASK.EQ.11 .OR. ATASK.EQ.21 .OR. ATASK.EQ.41 
!  Tasks which do stream I/O:
      STREAM =  ATASK.EQ.6 .OR. ATASK.EQ.60 .OR. ATASK.EQ.61
!
      MQ = SIZE(QTREE%QICELL,1)
      MQ = MIN(MQ,SIZE(QTREE%QPARENT,1))
      MQ = MIN(MQ,SIZE(QTREE%QCHILD,1))
      MQ = MIN(MQ,SIZE(QTREE%QLEVEL,1))
!
      MC = SIZE(QTREE%INDQUAD,1)
      MC = MIN(MC,SIZE(QTREE%INDSUB,1))
      MC = MIN(MC,SIZE(QTREE%INDLVL,1))
      MC = MIN(MC,SIZE(QTREE%NGBR,1))
      MC = MIN(MC,SIZE(QTREE%CELL_TYPE,1))
      MC = MIN(MC,SIZE(QTREE%INDML,1))
      MC = MIN(MC,SIZE(QTREE%XYVAL,1))
!
      IF ( DOWTS ) THEN
         MC2 = SIZE(QTREE%NCASE,1)
!
         MC3 = SIZE(QTREE%NGNBR,1)
         MC3 = MIN(MC3,SIZE(QTREE%GNBR,1))
         MC3 = MIN(MC3,SIZE(QTREE%INDWT,1))
         MC4 = MIN(36,SIZE(QTREE%GNBR,2))
      END IF
!
!    1. Scalars
!
      IF ( TASK.EQ.-1 .OR. TASK.EQ.-10 ) THEN 
            ! Sequential read
            NR = NR + 1
            READ(IUN,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL,           &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
      ELSEIF ( TASK.EQ.1 .OR. TASK.EQ.10 ) THEN
            ! Sequential write
            NR = NR + 1
            WRITE(IUN,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL,          &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
      ELSEIF ( TASK.EQ.-2 .OR. TASK.EQ.-20 ) THEN
            ! Direct access read
            NR = NR + 1
            READ(IUN,REC=NR,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL,    &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
      ELSEIF ( TASK.EQ.2 .OR. TASK.EQ.20 ) THEN
            ! Direct access write
            NR = NR + 1
            WRITE(IUN,REC=NR,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL,   &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
      ELSEIF ( TASK.EQ.-6 .OR. TASK.EQ.-60 ) THEN
            ! Stream access read
            NR = NR + 1
            RPOS  = 1_8 + LRECL*(NR-1_8)
            READ(IUN,POS=RPOS,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL,  &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
      ELSEIF ( TASK.EQ.6 .OR. TASK.EQ.60 ) THEN
            ! Stream access write
            NR = NR + 1
            RPOS  = 1_8 + LRECL*(NR-1_8)
            WRITE(IUN,POS=RPOS) WRITEBUFF
            WRITE(IUN,POS=RPOS,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL, &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
      ELSEIF ( TASK.EQ. -3 .OR. TASK.EQ. -4 .OR. TASK.EQ. -5 .OR.     &
               TASK.EQ.-30 .OR. TASK.EQ.-40 .OR. TASK.EQ.-50) THEN
            ! Formatted read
         READ(IUN,*,IOSTAT=IERR) NHD
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NHD'
         READ(IUN,*,IOSTAT=IERR) NCOL
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NCOL'
         READ(IUN,*,IOSTAT=IERR) QTREE%NQUAD
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NQUAD'
         READ(IUN,*,IOSTAT=IERR) QTREE%NCELL
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NCELL'
         READ(IUN,*,IOSTAT=IERR) QTREE%NCELL_DEF
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NCELL_DEF'
         READ(IUN,*,IOSTAT=IERR) QTREE%LVLREF
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'LVLREF'
         READ(IUN,*,IOSTAT=IERR) QTREE%LVLMAX
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NVLMAX'
         READ(IUN,*,IOSTAT=IERR) QTREE%LVLHI
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'LVLHI'
         READ(IUN,*,IOSTAT=IERR) QTREE%NX0
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NX0'
         READ(IUN,*,IOSTAT=IERR) QTREE%NY0
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NY0'
         READ(IUN,*,IOSTAT=IERR) QTREE%UNDEF_TYPE
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'UNDEF_TYPE'
         READ(IUN,*,IOSTAT=IERR) QTREE%KEEP_REF
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'KEEP_REF'
         READ(IUN,*,IOSTAT=IERR) QTREE%DYNAMIC
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'DYNAMIC'
         READ(IUN,*,IOSTAT=IERR) QTREE%IWTORDER
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'IWTORDER'
         DO IHD=1,NHD-14
            READ(IUN,*,IOSTAT=IERR)
            IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'Extra Header'
         END DO
      ELSEIF ( TASK.EQ.3  .OR. TASK.EQ.4  .OR. TASK.EQ.5 .OR.         &
               TASK.EQ.30 .OR. TASK.EQ.40 .OR. TASK.EQ.50) THEN
            ! Formatted write
         IF ( TASK.EQ.3 .OR. TASK.EQ.30 ) THEN
            !  quad variables only
            NCOL = 16
            NHD = NCOL + 18
         ELSEIF ( TASK.EQ.4 .OR. TASK.EQ.40 ) THEN
            !  cell variables only
            NCOL = 16
            NHD = NCOL + 18
         ELSEIF ( TASK.EQ.5 .OR. TASK.EQ.50 ) THEN
            !  cell weights variables only
            IF ( QTREE%IWTORDER.EQ.1 ) THEN
               NCOL = 2
               NHD = NCOL + 18
            ELSE IF ( QTREE%IWTORDER.GE.2 ) THEN
               NCOL = MC4 + 2
               NHD = 23
            END IF
         END IF 
         WRITE(IUN,1000,IOSTAT=IERR) NHD, NCOL
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NHD, NCOL'
         WRITE(IUN,1001,IOSTAT=IERR) QTREE%NQUAD, QTREE%NCELL,        &
               QTREE%NCELL_DEF, QTREE%LVLREF, QTREE%LVLMAX,           &
               QTREE%LVLHI, QTREE%NX0, QTREE%NY0, QTREE%UNDEF_TYPE,   &
               QTREE%KEEP_REF, QTREE%DYNAMIC, QTREE%IWTORDER
         IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'NQUAD-IWTORDER'
         IF ( TASK.EQ.3 .OR. TASK.EQ.30 ) THEN
            WRITE(IUN,FMT=1003,IOSTAT=IERR) 
            IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'Extra Header'
         ELSEIF ( TASK.EQ.4 .OR. TASK.EQ.40 ) THEN
            WRITE(IUN,FMT=1004,IOSTAT=IERR) 
            IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'Extra Header'
         ELSEIF ( TASK.EQ.5 .OR. TASK.EQ.50 ) THEN
            IF ( QTREE%IWTORDER.EQ.1 ) THEN
               WRITE(IUN,FMT=1005,IOSTAT=IERR) 
               IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'Extra Header'
            ELSE IF ( QTREE%IWTORDER.GE.2 ) THEN
               WRITE(IUN,FMT=1006,IOSTAT=IERR) 
               IF (IERR.GT.0 .AND. NDE.GT.0 ) WRITE(NDE,*) 'Extra Header'
            END IF
         END IF
      END IF
      IF ( IERR.NE.0 ) THEN
         IF ( NDE.GT.0 )                                              &
            WRITE(NDE,*) 'QA_IOQT: ERROR READING/WRITING SCALARS'
         RETURN
      END IF
!
! Bail if no weight only header data required
!
      IF ( ABS(TASK).EQ.10 .OR. ABS(TASK).EQ.20 .OR.                   &
           ABS(TASK).EQ.50 ) THEN
         IF ( PRESENT(NREC) ) NREC = NR
         RETURN
      END IF
!
!    2. Quad arrays
!
      NQ = QTREE%NQUAD
!
! Allocation check (only if quad array I/O will be called)
!
      IF ( NQ.GT.MQ .AND. DOQA ) THEN
         IF ( NDE.GT.0 ) THEN
            WRITE(NDE,*) 'QA_IOQT: INSUFFICIENT QUAD ALLOCATION'
            WRITE(NDE,*) 'ALLOCATED (MQ): ', MQ
            WRITE(NDE,*) 'REQUIRED (NQ): ', NQ
         END IF
         IERR = 1
         RETURN
      END IF
      IF ( TASK.EQ.-1 .OR. TASK.EQ.-11 ) THEN 
         ! Sequential read
         NR = NR + 1
         READ(IUN,IOSTAT=IERR) QTREE%QICELL(1:NQ,:),                   &
                           QTREE%QLEVEL(1:NQ), QTREE%QPARENT(1:NQ),    &
                           QTREE%QNBR(1:NQ,:), QTREE%QCHILD(1:NQ,:)
      ELSEIF ( TASK.EQ.1 .OR. TASK.EQ.11 ) THEN
         ! Sequential write
         NR = NR + 1
         WRITE(IUN,IOSTAT=IERR) QTREE%QICELL(1:NQ,:),                  &
                           QTREE%QLEVEL(1:NQ), QTREE%QPARENT(1:NQ),    &
                           QTREE%QNBR(1:NQ,:), QTREE%QCHILD(1:NQ,:)
      ELSEIF ( TASK.EQ.-2 .OR. TASK.EQ.-21 .OR.                        &
               TASK.EQ.-6 .OR. TASK.EQ.-61 ) THEN
         ! Direct access or stream access read
         NPART  = 1 + (NQ-1)/NSIZE
         DO IC=0,4
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QICELL(IR,IC),   &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QICELL(IR,IC),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
           END DO
         END DO
         DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
                RPOS  = 1_8 + LRECL*(NR-1_8)
                READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QLEVEL(IR),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
                READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QLEVEL(IR),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
         END DO
         DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
                RPOS  = 1_8 + LRECL*(NR-1_8)
                READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QPARENT(IR),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
                READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QPARENT(IR),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
         END DO
         DO IC=1,4
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QNBR(IR,IC),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QNBR(IR,IC),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
           END DO
         END DO
         DO IC=1,4
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QCHILD(IR,IC),   &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QCHILD(IR,IC),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
            END DO
         END DO
      ELSEIF ( TASK.EQ.2 .OR. TASK.EQ.21 .OR.                          &
               TASK.EQ.6 .OR. TASK.EQ.61 ) THEN
            ! Direct access or stream access write
         NPART  = 1 + (NQ-1)/NSIZE
         DO IC=0,4
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QICELL(IR,IC),  &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QICELL(IR,IC),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
           END DO
         END DO
         DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
                 RPOS  = 1_8 + LRECL*(NR-1_8)
                WRITE(IUN,POS=RPOS) WRITEBUFF
                WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QLEVEL(IR),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
                WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QLEVEL(IR),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
         END DO
         DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
                RPOS  = 1_8 + LRECL*(NR-1_8)
                WRITE(IUN,POS=RPOS) WRITEBUFF
                WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QPARENT(IR),   &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
                WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QPARENT(IR),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
         END DO
         DO IC=1,4
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QNBR(IR,IC),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QNBR(IR,IC),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
             IF ( IERR.NE.0 ) THEN
                IF ( NDE.GT.0 )                                        &
                   WRITE(NDE,*) 'QA_IOQT ERROR WRITING QNBR ',         &
                                     'IPART, NR = ', IPART, NR
                RETURN
             END IF
           END DO
         END DO
         DO IC=1,4
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%QCHILD(IR,IC),  &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%QCHILD(IR,IC),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NQ,IPART*NSIZE) )
             END IF
             IF ( IERR.NE.0 ) THEN
                IF ( NDE.GT.0 )                                        &
                   WRITE(NDE,*) 'QA_IOQT ERROR WRITING QCHILD ',       &
                                     'IPART, NR = ', IPART, NR
                RETURN
             END IF
           END DO
         END DO
      ELSEIF ( TASK.EQ.-3 .OR. TASK.EQ.-31 ) THEN
            ! Formatted read, quad variables
         DO IR=1,NQ
            READ(IUN,2003,IOSTAT=IERR) IR2,                           &
                           (QTREE%QICELL(IR,IC),IC=0,4),              &
                           QTREE%QLEVEL(IR), QTREE%QPARENT(IR),       &
                           (QTREE%QNBR(IR,IC),IC=1,4),                &
                           (QTREE%QCHILD(IR,IC),IC=1,4)
            IF ( IERR.NE.0 ) THEN
               IF ( NDE.GT.0 )                                        &
                  WRITE(NDE,*) 'QA_IOQT: ERROR READING QUAD INDEX ',IR
               RETURN
            END IF
         END DO
      ELSEIF ( TASK.EQ.3 .OR. TASK.EQ.31 ) THEN
            ! Formatted write, quad variables
         DO IR=1,NQ
            WRITE(IUN,2003,IOSTAT=IERR) IR,                           &
                           (QTREE%QICELL(IR,IC),IC=0,4),              &
                           QTREE%QLEVEL(IR), QTREE%QPARENT(IR),       &
                           (QTREE%QNBR(IR,IC),IC=1,4),                &
                           (QTREE%QCHILD(IR,IC),IC=1,4)
            IF ( IERR.NE.0 ) THEN
               IF ( NDE.GT.0 )                                        &
                  WRITE(NDE,*) 'QA_IOQT: ERROR WRITING QUAD INDEX ',IR
               RETURN
            END IF
         END DO
      END IF
      IF ( IERR.NE.0 ) THEN
         IF ( NDE.GT.0 )                                              &
            WRITE(NDE,*) 'QA_IOQT: I/O ERROR FOR QUAD ARRAYS'
         RETURN
      END IF
!
!    3. Cell arrays
!
      NC = QTREE%NCELL
      IF ( NC.GT.MC .AND. DOCA ) THEN
         IF ( NDE.GT.0 ) THEN
            WRITE(NDE,*) 'QA_IOQT: INSUFFICIENT CELL ALLOCATION'
            WRITE(NDE,*) 'ALLOCATED (MC): ', MC
            WRITE(NDE,*) 'REQUIRED (NC): ', NC
         END IF
         IERR = 1
         RETURN
      END IF
      IF ( TASK.EQ.-1 .OR. TASK.EQ.-11 ) THEN
            ! Sequential read
            NR = NR + 1
            READ(IUN,IOSTAT=IERR) QTREE%INDQUAD(1:NC),                 &
                         QTREE%INDSUB(1:NC), QTREE%INDLVL(1:NC),       &
                         QTREE%NGBR(1:NC,:), QTREE%CELL_TYPE(1:NC),    &
                         QTREE%INDML(1:NC), QTREE%XYVAL(1:NC,:)
      ELSEIF ( TASK.EQ.1 .OR. TASK.EQ.11 ) THEN
            ! Sequential write
            NR = NR + 1
            WRITE(IUN,IOSTAT=IERR) QTREE%INDQUAD(1:NC),                &
                         QTREE%INDSUB(1:NC), QTREE%INDLVL(1:NC),       &
                         QTREE%NGBR(1:NC,:), QTREE%CELL_TYPE(1:NC),    &
                         QTREE%INDML(1:NC), QTREE%XYVAL(1:NC,:)
      ELSEIF ( TASK.EQ.-2 .OR. TASK.EQ.-21 .OR.                        &
               TASK.EQ.-6 .OR. TASK.EQ.-61 ) THEN
            ! Direct access or stream access read
         NPART  = 1 + (NC-1)/NSIZE
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDQUAD(IR),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDQUAD(IR),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDSUB(IR),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDSUB(IR),        &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDLVL(IR),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDLVL(IR),        &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IC=1,8
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%NGBR(IR,IC),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%NGBR(IR,IC),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             END IF
           END DO
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%CELL_TYPE(IR),   &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%CELL_TYPE(IR),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDML(IR),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDML(IR),         &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IC=1,2
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               READ(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%XYVAL(IR,IC),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             ELSE
               READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%XYVAL(IR,IC),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             END IF
           END DO
         END DO
      ELSEIF ( TASK.EQ.2 .OR. TASK.EQ.21 .OR.                          &
               TASK.EQ.6 .OR. TASK.EQ.61 ) THEN
         ! Direct access write
         NPART  = 1 + (NC-1)/NSIZE
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDQUAD(IR),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDQUAD(IR),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDSUB(IR),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDSUB(IR),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDLVL(IR),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDLVL(IR),       &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IC=1,8
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%NGBR(IR,IC),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%NGBR(IR,IC),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             END IF
           END DO
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%CELL_TYPE(IR),  &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%CELL_TYPE(IR),    &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IPART=1,NPART
            NR = NR + 1
            IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%INDML(IR),      &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDML(IR),        &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
            END IF
         END DO
         DO IC=1,2
           DO IPART=1,NPART
             NR = NR + 1
             IF ( STREAM ) THEN
               RPOS  = 1_8 + LRECL*(NR-1_8)
               WRITE(IUN,POS=RPOS) WRITEBUFF
               WRITE(IUN,POS=RPOS,IOSTAT=IERR) ( QTREE%XYVAL(IR,IC),   &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             ELSE
               WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%XYVAL(IR,IC),     &
                           IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
             END IF
           END DO
         END DO
      ELSEIF ( TASK.EQ.-4 .OR. TASK.EQ.-41 ) THEN
            ! Formatted read, cell variables
         DO IR=1,NC
            READ(IUN,2004,IOSTAT=IERR) IR2, QTREE%INDQUAD(IR),         &
                         QTREE%INDSUB(IR), QTREE%INDLVL(IR),           &
                         (QTREE%NGBR(IR,IC),IC=1,8),                   &
                         QTREE%CELL_TYPE(IR), QTREE%INDML(IR),         &
                         (QTREE%XYVAL(IR,IC),IC=1,2)
            IF ( IERR.NE.0 ) THEN
               IF ( NDE.GT.0 )                                         &
                  WRITE(NDE,*) 'QA_IOQT: ERROR READING CELL INDEX ',IR
               RETURN
            END IF
         END DO
      ELSEIF ( TASK.EQ.4 .OR. TASK.EQ.41 ) THEN
            ! Formatted write, cell variables
         DO IR=1,NC
            WRITE(IUN,2004,IOSTAT=IERR) IR, QTREE%INDQUAD(IR),         &
                         QTREE%INDSUB(IR), QTREE%INDLVL(IR),           &
                         (QTREE%NGBR(IR,IC),IC=1,8),                   &
                         QTREE%CELL_TYPE(IR), QTREE%INDML(IR),         &
                         (QTREE%XYVAL(IR,IC),IC=1,2)
            IF ( IERR.NE.0 ) THEN
               IF ( NDE.GT.0 )                                         &
                  WRITE(NDE,*) 'QA_IOQT: ERROR WRITING CELL INDEX ',IR
               RETURN
            END IF
         END DO
      END IF
      IF ( IERR.NE.0 ) THEN
         IF ( NDE.GT.0 )                                              &
            WRITE(NDE,*) 'QA_IOQT: I/O ERROR FOR CELL ARRAYS'
         RETURN
      END IF
!
! Bail if no weight arrays required
!
      IF ( .NOT.DOWTS ) THEN
         IF ( PRESENT(NREC) ) NREC = NR
         RETURN
      END IF
!
!    4. Indices for first order weights:
!
      IF ( QTREE%IWTORDER.GE.1 .AND. MC2.GT.0 ) THEN
         IF ( NC.GT.MC2 ) THEN
            IF ( NDE.GT.0 ) THEN
               WRITE(NDE,*) 'QA_IOQT: CELL ALLOCATION ERROR'
               WRITE(NDE,*) 'ALLOCATED (MC2): ', MC2
               WRITE(NDE,*) 'REQUIRED (NC): ', NC
            END IF
            IERR = 2
            RETURN
         END IF
         !IF ( TASK.EQ.-1 ) THEN
         !   ! Sequential read
         !   NR = NR + 1
         !   READ(IUN,IOSTAT=IERR) QTREE%NCASE(1:NC)
         !ELSEIF ( TASK.EQ.1 ) THEN
         !   ! Sequential write
         !   NR = NR + 1
         !   WRITE(IUN,IOSTAT=IERR) QTREE%NCASE(1:NC)
         !ELSEIF ( TASK.EQ.-2 ) THEN
         !   ! Direct access read
         !   NPART  = 1 + (NC-1)/NSIZE
         !   DO IPART=1,NPART
         !      NR = NR + 1
         !      READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%NCASE(IR),        &
         !               IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !   END DO
         !ELSEIF ( TASK.EQ.2 ) THEN
         !   ! Direct access write
         !   NPART  = 1 + (NC-1)/NSIZE
         !   DO IPART=1,NPART
         !      NR = NR + 1
         !      WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%NCASE(IR),       &
         !               IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !   END DO
         !ELSEIF ( TASK.EQ.-5 .AND. QTREE%IWTORDER.EQ.1 ) THEN
         IF ( (TASK.EQ.-5 .OR. TASK.EQ.-51) .AND.                     &
               QTREE%IWTORDER.EQ.1 ) THEN
            ! Formatted read
            DO IR=1,NC
               READ(IUN,2005,IOSTAT=IERR) IR2, QTREE%NCASE(IR)
            END DO
         ELSEIF ( (TASK.EQ.5 .OR. TASK.EQ.51) .AND.                   &
                  QTREE%IWTORDER.EQ.1 ) THEN
            ! Formatted write
            DO IR=1,NC
               WRITE(IUN,2005,IOSTAT=IERR) IR, QTREE%NCASE(IR)
            END DO
         END IF
         IF ( IERR.NE.0 ) THEN
            IF ( NDE.GT.0 )                                          &
               WRITE(NDE,*) 'QA_IOQT: I/O ERROR FOR NCASE ARRAY'
            RETURN
         END IF
      END IF
!
!    5. Tables for higher order weights:
!
      IF ( QTREE%IWTORDER.GE.2 .AND. MC3.GT.0 ) THEN
         IF ( NC.GT.MC3 ) THEN
            IF ( NDE.GT. 0 ) THEN
               WRITE(NDE,*) 'QA_IOQT: CELL ALLOCATION ERROR'
               WRITE(NDE,*) 'ALLOCATED (MC3): ', MC3
               WRITE(NDE,*) 'REQUIRED (NC): ', NC
               IERR = 3
               RETURN
            END IF
         END IF
         !IF ( TASK.EQ.-1 ) THEN
         !      ! Sequential read
         !      NR = NR + 1
         !      READ(IUN,IOSTAT=IERR) QTREE%INDWT(1:NC),               &
         !                   QTREE%NGNBR(1:NC), QTREE%GNBR(1:NC,:)
         !ELSEIF ( TASK.EQ.1 ) THEN
         !      ! Sequential write
         !      NR = NR + 1
         !      WRITE(IUN,IOSTAT=IERR) QTREE%INDWT(1:NC),              &
         !                   QTREE%NGNBR(1:NC), QTREE%GNBR(1:NC,:)
         !ELSEIF ( TASK.EQ.-2 ) THEN
         !      ! Direct access read
         !      NPART  = 1 + (NC-1)/NSIZE
         !      DO IPART=1,NPART
         !         NR = NR + 1
         !         READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDWT(IR),     &
         !                  IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !      END DO
         !      DO IPART=1,NPART
         !         NR = NR + 1
         !         READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%NGNBR(IR),     &
         !                  IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !      END DO
         !      DO IC=1,MC4
         !        DO IPART=1,NPART
         !          NR = NR + 1
         !          READ(IUN,REC=NR,IOSTAT=IERR) ( QTREE%GNBR(IR,IC),  &
         !                  IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !        END DO
         !      END DO
         !ELSEIF ( TASK.EQ.2 ) THEN
         !      ! Direct access write
         !      NPART  = 1 + (NC-1)/NSIZE
         !      DO IPART=1,NPART
         !         NR = NR + 1
         !         WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%INDWT(IR),    &
         !                  IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !      END DO
         !      DO IPART=1,NPART
         !         NR = NR + 1
         !         WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%NGNBR(IR),    &
         !                  IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !      END DO
         !      DO IC=1,MC4
         !        DO IPART=1,NPART
         !          NR = NR + 1
         !          WRITE(IUN,REC=NR,IOSTAT=IERR) ( QTREE%GNBR(IR,IC), &
         !                  IR=1+(IPART-1)*NSIZE, MIN(NC,IPART*NSIZE) )
         !        END DO
         !      END DO
         !ELSEIF ( TASK.EQ.-5 .AND. QTREE%IWTORDER.GE.2 ) THEN
         IF ( (TASK.EQ.-5 .OR. TASK.EQ.-51) .AND.                     &
               QTREE%IWTORDER.GE.2 ) THEN
            ! Formatted read
            DO IR=1,NC
               READ(IUN,2006,IOSTAT=IERR) IR2, QTREE%NCASE(IR),       &
                            QTREE%INDWT(IR), QTREE%NGNBR(IR),         &
                            (QTREE%GNBR(IR,IC),IC=1,MC4)
            END DO
         ELSEIF ( (TASK.EQ.5 .OR. TASK.EQ.51) .AND.                   &
                   QTREE%IWTORDER.GE.2 ) THEN
            ! Formatted write
            DO IR=1,NC
               WRITE(IUN,2006,IOSTAT=IERR) IR, QTREE%NCASE(IR),       &
                            QTREE%INDWT(IR), QTREE%NGNBR(IR),         &
                            (QTREE%GNBR(IR,IC),IC=1,MC4)
            END DO
         END IF
         IF ( IERR.NE.0 ) THEN
            IF ( NDE.GT. 0 )                                          &
               WRITE(NDE,*) 'QA_IOQT: I/O ERROR FOR 2ND ORDER ARRAYS'
            RETURN
         END IF
      END IF
!
! 6. Return value of NREC
!
      IF ( PRESENT(NREC) ) NREC = NR
      RETURN
!
! 7. Format statements
!
 1000 FORMAT(1X,I4,'  No. of lines in the header'/                    &
             1X,I4,'  No. of columns of data')
!
 1001 FORMAT(1X,I8,'  NQUAD      No. of quads'/                       &
             1X,I8,'  NCELL      Max. No. of cells'/                  &
             1X,I8,'  NCELL_DEF  No. of valid cells'/                 &
             1X,I8,'  LVLREF     Reference level'/                    &
             1X,I8,'  LVLMAX     Max. refinement level allowed'/      &
             1X,I8,'  LVLHI      Max. refinement level reached'/      &
             1X,I8,'  NX0        No. of X cells in level-0 grid'/     &
             1X,I8,'  NY0        No. of Y cells in level-0 grid'/     &
             1X,I8,'  UNDEF_TYPE Marker for non-valid cells'/         &
             1X,L8,'  KEEP_REF   Retain refined cells?'/              &
             1X,L8,'  DYNAMIC    Nonstationary grid?'/                &
             1X,I8,'  IWTORDER   Order of weight tables included')
!
 1003 FORMAT('!---------------------------------------------------'/  &
             '!  Quadtree quad data'/                                 &
             '! Columns:'/                                            &
             '!  1. iq  = quad index'/                                &
             '!  2. QICELL(iq,0) = central cell index'/               &
             '!  3. QICELL(iq,1) = SW cell index'/                    &
             '!  4. QICELL(iq,2) = SE cell index'/                    &
             '!  5. QICELL(iq,3) = NW cell index'/                    &
             '!  6. QICELL(iq,4) = NE cell index'/                    &
             '!  7. QLEVEL(iq)   = level of child cells'/             &
             '!  8. QPARENT(iq)  = index of parent quad'/             &
             '!  9. QNBR(iq,1)   = W neighbour quad index'/           &
             '! 10. QNBR(iq,2)   = E neighbour quad index'/           &
             '! 11. QNBR(iq,3)   = S neighbour quad index'/           &
             '! 12. QNBR(iq,4)   = N neighbour quad index'/           &
             '! 13. QCHILD(iq,1) = SW child quad index'/              &
             '! 14. QCHILD(iq,2) = SE child quad index'/              &
             '! 15. QCHILD(iq,3) = NW child quad index'/              &
             '! 16. QCHILD(iq,4) = NE child quad index'/              &
             '!---------------------------------------------------')
!
 1004 FORMAT('!---------------------------------------------------'/  &
             '!  Quadtree cell data'/                                 &
             '! Columns:'/                                            &
             '!  1. ic  = cell index'/                                &
             '!  2. INDQUAD(ic)   = quad index'/                      &
             '!  3. INDSUB(ic)    = sub-quad index(1-4=SW,SE,NW,NE)'/ &
             '!  4. INDLVL(ic)    = refinement level'/                &
             '!  5. NGBR(ic,1)    = W primary neighbour index'/       &
             '!  6. NGBR(ic,2)    = E primary neighbour index'/       &
             '!  7. NGBR(ic,3)    = S primary neighbour index'/       &
             '!  8. NGBR(ic,4)    = N primary neighbour index'/       &
             '!  9. NGBR(ic,5)    = W secondary neighbour index'/     &
             '! 10. NGBR(ic,6)    = E secondary neighbour index'/     &
             '! 11. NGBR(ic,7)    = S secondary neighbour index'/     &
             '! 12. NGBR(ic,8)    = N secondary neighbour index'/     &
             '! 13. CELL_TYPE(ic) = cell status'/                     &
             '! 14. INDML(ic)     = multilevel grid index'/           &
             '! 15. XYVAL(ic,1)   = X-coordinate in reference grid'/  &
             '! 16. XYVAL(ic,2)   = Y-coordinate in reference grid'/  &
             '!---------------------------------------------------')
!
 1005 FORMAT('!---------------------------------------------------'/  &
             '!  Quadtree weights data ( 1st order)'/                 &
             '! Columns:'/                                            &
             '!  1. ic  = cell index'/                                &
             '!  2. NCASE(ic)   = index into 1st-order weights table'/&
             '!---------------------------------------------------')
!
 1006 FORMAT('!---------------------------------------------------'/  &
             '!  Quadtree weights data (1st and 2nd order)'/          &
             '! Columns:'/                                            &
             '!  1. ic  = cell index'/                                &
             '!  2. NCASE(ic)  = index into 1st-order weights table'/ &
             '!  3. INDWT(ic)  = index into 2nd-order weights table'/ &
             '!  4. NGNBR(ic)  = number of 2nd-order neighbours'/     &
             '! 5+. GNBR(ic,:) = 2nd-order neighbour indices'/        &
             '!---------------------------------------------------')
!
 2003 FORMAT(1X,I8,5I8,1X,I4,1X,I8,1X,4I8,1X,4I8)
!
 2004 FORMAT(1X,I8,1X,I8,1X,I4,1X,I4,1X,8I8,1X,I4,1X,I8,1X,2F12.6)
!
 2005 FORMAT(2(1X,I8))
!
 2006 FORMAT(4(1X,I8),1X,36(1X,I8))
!/
!/ End of QA_IOQT  ----------------------------------------------------- /
!/
      END SUBROUTINE QA_IOQT
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_LINGAUSS(N, NSOL, A, B, X, DET, STATUS)
!/
!/    07-Oct-2010 : Origination. Richard Gorman, NIWA
!/
!  1. Purpose :
!
!     Solve a system of linear equations A*X = B 
!
!  2. Method :
!
!     Gaussian elimination with pivoting. Code adapted from PL/I algorithm
!     on pp 76-77 of 
!         Pizer, Stephen M. (1975), "Numerical Computing and Mathematical
!         Analysis", SRA Inc., Chigago, 529 pp.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       N       Int.   I   Number of equations, unknowns.
!       NSOL    Int.   I   Number of cases to solve
!       A       R.A.   I   Equations matrix
!       B       R.A.   I   Matrix of RHS vectors
!       X       R.A.   O   Matrix of solution vectors
!       DET     Real   O   Determinant
!       STATUS  Int.   O   Error flag
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_DERIVWTS Subr. qa_utils Compute weights for higher order derivatives
!                                on an extended stencil 
!     ---------------------------------------------------------------- 
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: N, NSOL
      REAL, INTENT(INOUT)     :: A(N,N), B(N,NSOL)
      REAL, INTENT(OUT)       :: X(N,NSOL)
      REAL, INTENT(OUT)       :: DET
      INTEGER, INTENT(OUT)    :: STATUS
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER           ::  P(N)
      REAL              ::  SCAL(N)
      INTEGER           ::  JMAX, I, J, K, PTEMP
      REAL              ::  SCALVAL, PMAX, AMULT, TOL
!/
!/ ------------------------------------------------------------------- /
!/
!
!   1. Initialisation
!
      TOL = 1.E-9
      DET = 0.
      DO I = 1, N
        P(I) = I
        SCAL(I) = 0.
        DO J=1,N
          SCAL(I) = SCAL(I) + ABS(A(I,J))
        END DO
        IF ( ABS(SCAL(I)).LT.TOL ) THEN
          STATUS = I
          RETURN
        END IF
      END DO
      !write(*,*) 'A, B:'
      !do i=1,n
      !  write(*,*) (A(i,k),k=1,n),'  ',(B(i,k),k=1,Nsol)
      !end do
!
!   2. Triangularisation with pivoting
!
      DO I = 1, N-1
        PMAX = 0.
        JMAX = 1
        ! Find the pivot:
        DO J=I,N
          SCALVAL = ABS(A(P(J),I)/SCAL(P(J)))
          IF ( SCALVAL.GT.PMAX ) THEN
            PMAX = SCALVAL
            JMAX = J
          END IF
        END DO
        ! Switch rows:
        PTEMP = P(JMAX)
        P(JMAX) = P(I)
        P(I) = PTEMP
        IF ( ABS(A(P(I),I)).LT.TOL ) THEN
          STATUS = -2
          RETURN
        END IF
        ! Below ith row:
        DO J=I+1,N
          ! Multiplier for jth row:
          AMULT = A(P(J),I)/A(P(I),I)
          ! Compute nonzero elements of A, B for jth row:
         !   A(P(J),I) = 0.
          DO K=I+1,N
            A(P(J),K) = A(P(J),K) - AMULT*A(P(I),K)
          END DO
          DO K=1,NSOL
            B(P(J),K) = B(P(J),K) - AMULT*B(P(I),K)
          END DO        
        END DO
        !write(*,*) 'A, B, step:',I
        !do j=1,n
        !  write(*,*) (A(j,k),k=1,n),'  ',(B(j,k),k=1,Nsol)
        !end do
      END DO
!
!   3. Determinant
!
      DET = 1.
      DO I=1,N
         DET = DET*A(P(I),I)
      END DO
!
!   4. Back substitution
!
      DO K = 1,NSOL
        DO I = N,1,-1
          X(I,K) = B(P(I),K)
          DO J=I+1,N
            X(I,K) = X(I,K) - A(P(I),J)*X(J,K)
          END DO
          X(I,K) = X(I,K)/A(P(I),I)
        END DO
      END DO
        !write(*,*) 'A, X, final:'
        !do j=1,n
        !  write(*,*) (A(j,k),k=1,n),'  ',(X(j,k),k=1,Nsol)
        !end do
      STATUS = 0
!/
!/ End of QA_LINGAUSS --------------------------------------------------- /
!/
      END SUBROUTINE QA_LINGAUSS
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_MLG2QT( LVLREL, LVLSG, LVLREF, NX0, NY0, X0, Y0,  &
                            ARRIN, QTREE, MAPML, ARRQ, ATYPE, ARANGE, &
                            MAP_TYPES, IERR, NDSE )
!/
!       Richard Gorman, NIWA
!         
!         April 2014:   Origination, derived from QA_MLGRID, QA_BMLSTRUC.
!         Nov   2017:   Add ATYPE = 5.
!
!  1. Purpose :
!
!      Precompute array values, e.g. depth and transmission arrays, 
!      for cells in part of a quadtree bathymetry grid corresponding
!      to a given subgrid for which high resolution data is provided.
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       LVLREL  Int.   I   Number of refinement levels below the input 
!                          subgrid to compute
!       LVLSG   Int.   I   Refinement level of the subgrid (relative 
!                          to reference level)
!       LVLREF  Int.   I   Reference refinement level
!       NX0,NY0 Int.   I   Size of level-zero grid
!       X0      Real   I   West wall coordinate of the input subgrid
!       Y0      Real   I   South wall coordinate of the input subgrid
!       ARRIN   R.A.   I   [NX x NY x NARR] Array of NARR variables defined
!                          on the [NX x NY] fine grid 
!       QTREE   QA_TREE IO Quadtree structure.
!       MAPML   I.A.   IO  Flags for valid/invalid cells on the full multilevel 
!                          grid
!       ARRQ    R.A.   IO  [MCELL x NARR] Array of NARR variables defined
!                          on the quadtree grid 
!       ATYPE   I.A.   I*  [NARR x 1] array of flags of how each of the NARR 
!                          variables should be averaged over 4 cells in going 
!                          from a finer to a coarser grid.
!                             ATYPE = 1: simple averaging over wet cells [default]
!                             ATYPE = 2: direction averaging over wet cells
!                             ATYPE = 3: x-transmission averaging
!                             ATYPE = 4: y-transmission averaging
!                             ATYPE = 5: retain value from the first wet cell
!       ARANGE  R.A.   I*  [NARR x 2] array of min. & max valid values for
!                          the NARR variables in ARRIN. [default = -Inf to Inf].
!       MAP_TYPES I.A. I*  Values of MAPML for cells that: 
!                             haven't been determined (MAP_TYPES(1))
!                             are default active sea points (MAP_TYPES(2))
!                             should not be adapted (MAP_TYPES(3:end))
!       IERR     Int.  O*  Return flag = 1 for error, else 0
!       NDSE     Int.  I*  Unit number for error output (if >0)
!
!     ----------------------------------------------------------------
!                        * optional
!
!  4. Subroutines used :
!
!      Name             Type  Module   Description
!     ---------------------------------------------------------------
!     MULTILEVEL_INDEX  Func. qa_utils Multilevel index of a given cell
!     QA_XY2CELL        Subr. qa_utils Identify quadtree cell for a given (X,Y)
!     QA_REFINE         Subr. qa_utils refine a quadtree grid
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
      INTEGER, INTENT(IN)            :: LVLREL
      INTEGER, INTENT(IN)            :: LVLSG, LVLREF, NX0, NY0
      REAL, INTENT(IN)               :: X0, Y0
      REAL,    INTENT(IN)            :: ARRIN(:,:,:)
      TYPE(QA_TREE), INTENT(INOUT)   :: QTREE
      INTEGER, INTENT(INOUT)         :: MAPML(:)
      REAL,    INTENT(OUT)           :: ARRQ(:,:)
      INTEGER, OPTIONAL, INTENT(IN)  :: ATYPE(:)
      REAL, OPTIONAL,    INTENT(IN)  :: ARANGE(:,:)
      INTEGER, OPTIONAL, INTENT(IN)  :: MAP_TYPES(:)
      INTEGER, OPTIONAL, INTENT(OUT) :: IERR
      INTEGER, OPTIONAL, INTENT(IN)  :: NDSE
!
! Local variables.
      INTEGER              :: MINDBG, NX, NY, NARR, MARR
      INTEGER              :: NXSUB, NYSUB, LVL, II, JJ, INDBG,       &
                              NXP, NYP
      INTEGER              :: INDBGMAX, IJS(4), IS, IBG, NWET, IARR,  &
                              IQ, ISUB
      INTEGER              :: INDBG1, IJT, LEVEL, NB, MC, IMLG, ITER
      REAL                 :: ZVAL
      REAL                 :: X, Y, XNP, YNP, DELXY, ZDIFF,           &
                              XCELL, YCELL
      INTEGER              :: ICELL, LVCELL, ISTAT(2), NREF, IREF,    &
                              IMODE
      INTEGER              :: IERS, NT_NOADAPT, NMT
      INTEGER              :: VALID_TYPE, UNDEF_TYPE, UNSET_TYPE
      LOGICAL              :: NOREF, SETWET
      INTEGER              :: IUN
      REAL, ALLOCATABLE    :: ARR_RANGE(:,:)
      INTEGER, ALLOCATABLE :: ARRTYPE(:)
      INTEGER, ALLOCATABLE :: NO_ADAPT_TYPE(:)
      REAL, ALLOCATABLE    :: ARRVAL(:,:), ARRMEAN(:), ARRMAX(:)
      INTEGER, ALLOCATABLE :: INDML(:)
      LOGICAL, ALLOCATABLE :: ISWET(:)
      REAL, ALLOCATABLE    :: XML(:), YML(:), ARRML(:,:)
      INTEGER, ALLOCATABLE :: ISEAOLD(:), ISEANEW(:,:)
      INTEGER, ALLOCATABLE :: IQREF(:), ISREF(:)
      !    integer ::       NWETALL
          integer ::      NINSUB, NSET, NFINE, NRPT, NCRSNB, NWETALL
!
! Array sizes
!      
      NX = SIZE(ARRIN,1)
      NY = SIZE(ARRIN,2)
      NARR = SIZE(ARRIN,3)
      NARR = MIN(NARR,SIZE(ARRQ,2))
!
! Defaults for optional parameters:
!
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0
      !
      ALLOCATE ( ARR_RANGE(NARR,2) )
      ARR_RANGE(:,1) = -HUGE(1.)
      ARR_RANGE(:,2) = HUGE(1.)
      IF ( PRESENT (ARANGE) ) THEN
         MARR = MIN(NARR,SIZE(ARANGE,1))
         MC = MIN(2,SIZE(ARANGE,2))
         ARR_RANGE(1:MARR,1:MC) = ARANGE(1:MARR,1:MC)
      END IF
      !
      ALLOCATE ( ARRTYPE(NARR) )
      ARRTYPE = 1
      IF ( PRESENT (ATYPE) ) THEN
         MARR = MIN(NARR,SIZE(ATYPE,1))
         ARRTYPE(1:MARR) = ATYPE(1:MARR)
      END IF
!
!     MAPML value signifying: 
      UNDEF_TYPE = QTREE%UNDEF_TYPE  ! permanently excluded cells
      VALID_TYPE = 1 - UNDEF_TYPE    ! default for non-excluded cells
      UNSET_TYPE = -9                ! yet to be determined
      NT_NOADAPT = 0
      IF ( PRESENT(MAP_TYPES) ) THEN
         NMT = SIZE(MAP_TYPES,1) 
         IF ( NMT.GE.1 ) THEN
            UNSET_TYPE = MAP_TYPES(1)
            IF ( NMT.GE.2 ) THEN
               VALID_TYPE = MAP_TYPES(2)
               IF ( NMT.GE.3 ) THEN
                  NT_NOADAPT = NMT - 2
                  ALLOCATE ( NO_ADAPT_TYPE(NT_NOADAPT) )
                  NO_ADAPT_TYPE(:) = MAP_TYPES(3:NMT)
               END IF
            END IF
         END IF
      END IF
            
      !
      ALLOCATE ( ARRMEAN(NARR), ARRMAX(NARR) )
      ALLOCATE ( ARRVAL(4,NARR) )
!
!    Number of grid cells covered by this subgrid at the coarsest
!    level to be computed:
      XNP = FLOAT(NX)*2.**(-LVLREL)
      NXP = NINT(XNP)
      YNP = FLOAT(NY)*2.**(-LVLREL)
      NYP = NINT(YNP)
      IF ( ABS(XNP-NXP).GT.1.E-9 .OR. ABS(YNP-NYP).GT.1.E-9 ) THEN

         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_MLG2QT: '
            WRITE(IUN,*) 'MISMATCH IN NXP, NYP'
            WRITE(IUN,*) 'NX, NY = ', NX, NY
            WRITE(IUN,*) 'LVLREL = ', LVLREL
            WRITE(IUN,*) 'NXP, NYP = ', NXP, NYP
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
!
!  No. of cells in the "local" multilevel grid covering this subgrid, 
!  up to the reference level
      MINDBG = NXP*NYP*NINT( (4.**(LVLREL + 1) - 1)/3. )
      ALLOCATE ( INDML(MINDBG), ISWET(MINDBG) )
      ALLOCATE ( XML(MINDBG), YML(MINDBG) )
      ALLOCATE ( ARRML(MINDBG,NARR) )
      ALLOCATE ( ISEAOLD(MINDBG), ISEANEW(MINDBG,4))
      ALLOCATE ( IQREF(MINDBG), ISREF(MINDBG) )
!      
      ISWET = .TRUE.
      !write(*,*) 'QA_MLG2QT: LVLREL, LVLSG, LVLREF: ', LVLREL, LVLSG, LVLREF
      !write(*,*) 'QA_MLG2QT: ARR_RANGE: '
      !do iarr=1,NARR
      !   write(*,*) ARR_RANGE(IARR,:)
      !end do
      
!
! Loop over level relative to the input subgrid, from finest to coarsest:
!
      NXSUB = NX
      NYSUB = NY
      NWETALL = 0
      DO LVL=0,-LVLREL,-1
!
! Absolute refinement level:
         LEVEL = LVL + LVLSG + LVLREF
!
!       Highest (local) multi-level bathy grid index for this level:
         INDBGMAX = NXP*NYP*NINT( (4.**(LVL + LVLREL +1) - 1)/3. ) 
!
!       Highest (local) multi-level bathy grid index for the next-lowest level:
         INDBG = NXP*NYP*NINT( (4.**(LVL + LVLREL) - 1)/3. ) 
!
         DELXY = 2.**(-LVL-LVLSG)
         DO JJ=1,NYSUB
            Y = Y0 + (FLOAT(JJ) - 0.5)*DELXY
            DO II=1,NXSUB
               X = X0 + (FLOAT(II) - 0.5)*DELXY
               ! Index into local multilevel grid:
               INDBG = INDBG + 1
               ! Index into global multilevel grid:
               IMLG = MULTILEVEL_INDEX ( NX0, NY0, X, Y, LEVEL,       &
                                         LVLREF )
               INDML(INDBG) = IMLG
               ! Skip points outside the reference grid
               IF ( IMLG.EQ.0 ) CYCLE
               XML(INDBG) = X
               YML(INDBG) = Y
               IJS = 0
               IF (LVL.EQ.0) THEN
                  ! For cells on the input regular grid, use MAPML as 
                  ! provided as long as it has been preset   
                  SETWET = MAPML(IMLG).EQ.UNSET_TYPE
                  IF ( SETWET ) THEN
                     ISWET(INDBG) = .TRUE.
                     MAPML(IMLG) = VALID_TYPE
                  ELSE
                     ISWET(INDBG) = MAPML(IMLG).NE.UNDEF_TYPE
                  END IF
                  DO IARR=1,NARR
                     ZVAL = ARRIN(II,JJ,IARR)
                     ARRML(INDBG,IARR) = ZVAL
                     IF ( SETWET .AND. ( ZVAL.LT.ARR_RANGE(IARR,1)    &
                         .OR. ZVAL.GT.ARR_RANGE(IARR,2) ) ) THEN
                        ISWET(INDBG) = .FALSE.
                        MAPML(IMLG) = UNDEF_TYPE
                     END IF
                  END DO
               ELSE
!
!             Four higher-level cells within this cell:
!               (1: SW, 2: SE, 3: NW, 4: NE)
                  IBG = INDBGMAX + 4*NXSUB*(JJ-1) + 2*(II-1) + 1
                  IJS(1) = IBG
                  IJS(2) = IBG+1
                  IBG = IBG + 2*NXSUB
                  IJS(3) = IBG
                  IJS(4) = IBG+1
!
!             Average data over wet cells:
                  NWET = 0
                  ARRMEAN = 0.
                  ARRMAX = -HUGE(1.)
                  DO IS=1,4
                     IJT = IJS(IS)
                     IF ( ISWET(IJT) ) NWET = NWET + 1
                     DO IARR=1,NARR
                        ZVAL = ARRML(IJT,IARR)
                        ARRMAX(IARR) = MAX(ARRMAX(IARR), ZVAL)
                        !ARRMIN(IARR) = MIN(ARRMIN(IARR), ZVAL)
                        IF ( ISWET(IJT) ) THEN
                           ! difference from running mean
                           ZDIFF = ZVAL - ARRMEAN(IARR)
                           ! for direction, put this in range [-180 180]
                           IF ( ARRTYPE(IARR).EQ.2 ) ZDIFF =          &
                                       MOD( ZDIFF+180,.360) - 180.
                           ! Accumulate the running mean
                           ARRMEAN(IARR) = ARRMEAN(IARR) + ZDIFF/NWET
                           ! Store individual value
                           ARRVAL(IS,IARR) = ZVAL
                           ! Type 5 (e.g. for interpolation weights and 
                           ! indices) move on now thay we have the first
                           ! wet value
                           IF ( ARRTYPE(IARR).EQ.5 ) CYCLE
                        ELSE
                           ARRVAL(IS,IARR) = 0.
                        END IF  ! subcell wet/dry
                     END DO   ! loop over data arrays
                  END DO  ! loop over 4 subcells
                  IF ( NWET.GT.0 ) THEN
                  ! If there are ANY wet higher-level cells, consider this wet:
                     ISWET(INDBG) = .TRUE.
                     ! If the cell's status is unset, or previously marked as
                     ! excluded, change that to the default wet value
                     IF ( MAPML(IMLG).EQ.UNDEF_TYPE .OR.              &
                          MAPML(IMLG).EQ.UNSET_TYPE )                 &
                                  MAPML(IMLG) = VALID_TYPE
                     DO IARR=1,NARR
                        IF ( ARRTYPE(IARR).EQ.3 ) THEN
                        !  Compute X Transmission factors, taking 
                        !  products of successive cells in the downwave 
                        !  direction and averaging in the cross-wave direction:
                           ARRML(INDBG,IARR) =                        &
                              0.5*( ARRVAL(1,IARR)*ARRVAL(2,IARR) +   &
                                    ARRVAL(3,IARR)*ARRVAL(4,IARR) )
                        ELSEIF ( ARRTYPE(IARR).EQ.4 ) THEN
                        !  Compute Y Transmission factors, " "
                           ARRML(INDBG,IARR) =                        &
                              0.5*( ARRVAL(1,IARR)*ARRVAL(3,IARR) +   &
                                    ARRVAL(2,IARR)*ARRVAL(4,IARR) )
                        ELSE
                        !  Standard or direction scalar: take mean over valid cells
                        !  This will also work for ATYPE=5, where only the first
                        !  valid cell was used
                           ARRML(INDBG,IARR) = ARRMEAN(IARR)
                        END IF
                     END DO
                  ELSE
                  ! If there are no wet higher-level cells:
                     ISWET(INDBG) = .FALSE.
                     MAPML(IMLG) = UNDEF_TYPE
                     DO IARR=1,NARR
                        IF ( ARRTYPE(IARR).EQ.1 ) THEN
                        !  use the maximum for standard scalars, e.g. "mean" 
                        !  bed elevation:
                           ARRML(INDBG,IARR) = ARRMAX(IARR)
                        ELSE
                           ARRML(INDBG,IARR) = 0.
                        END IF
                     END DO  ! number of data arrays
                  END IF  ! coarse cell wet/dry
               END IF   ! LVL.EQ.0
               !write(*,*) 'LVL,LEVEL,II,JJ,INDBG,IMLG,X,Y,ISWET,IJS,ARRML(1),MAPML: ',  &
               !         LVL,LEVEL,II,JJ,INDBG,IMLG,X,Y,ISWET(INDBG), &
               !            IJS, ARRML(INDBG,1), MAPML(IMLG)
               IF ( ISWET(INDBG) ) NWETALL = NWETALL + 1
            END DO  ! X grid index
         END DO  ! Y grid index
         write(*,*) 'QA_MLG2QT: LVL, LEVEL, NWETALL: ',LVL, LEVEL, NWETALL
         IF ( MOD(NXSUB,2).NE.0 ) EXIT
         IF ( MOD(NYSUB,2).NE.0 ) EXIT
         NXSUB = NXSUB/2
         NYSUB = NYSUB/2
      END DO
!
! Loop over level relative to the input subgrid, this time from coarsest
! to finest:
!
      DO LVL=-LVLREL,0,+1
!
! Absolute refinement level:
         LEVEL = LVL + LVLSG + LVLREF
!
!       Highest multi-level bathy grid index for this level:
         INDBGMAX = NXP*NYP*NINT( (4.**(LVL + LVLREL +1) - 1)/3. ) 
!
!       Lowest multi-level bathy grid index for this level:
         INDBG1 = NXP*NYP*NINT( (4.**(LVL + LVLREL) - 1)/3. ) + 1
!
! iterations: to assign data to quadtree cells created
! by refinement in a previous iteration. May need enough to refine
! all the way from level 0.
!
         DO ITER=1,LEVEL+1
!
!       Loop over all cells at this level
            NREF = 0
            ISEAOLD = 0
            IQREF = 0
            ISREF = 0
            NINSUB = 0
            NSET = 0
            NFINE = 0
            NRPT = 0
            NCRSNB = 0
            NWETALL = 0
            DO INDBG = INDBG1, INDBGMAX
               IMLG = INDML(INDBG)
               ! Skip points outside the reference grid
               IF ( IMLG.EQ.0 ) CYCLE
               X = XML(INDBG)
               Y = YML(INDBG)
               CALL QA_XY2CELL ( QTREE, X, Y, ICELL, XCELL, YCELL,    &
                                 LVCELL, IQ, ISUB, ISTAT )
               !write(*,*) 'INDBG,X,Y,ISWET,ICELL,XCELL,YCELL,LVCELL,IQ,ISUB,ISTAT:',  &
               !            INDBG,X,Y,ISWET(INDBG),ICELL,XCELL,YCELL,LVCELL,IQ,ISUB,ISTAT
               ! If (X,Y) is outside the subgrid, go to the next point
               IF ( .NOT.ALL(ISTAT.LE.0) ) CYCLE
               ! Otherwise, we have found a quad (IQUAD) containing (X,Y):
               ! and perhaps a cell (ICELL) containing (X,Y):
               NINSUB = NINSUB + 1
               IF ( ICELL.GT.0 .AND. LVCELL.EQ.LEVEL ) THEN
                   ! There is an exact matching cell in the quadtree:
                   ! apply the data to that
                  DO IARR=1,NARR
                     ARRQ(ICELL,IARR) = ARRML(INDBG,IARR)
                  END DO
                  QTREE%CELL_TYPE(ICELL) = MAPML(IMLG)
                  NSET = NSET + 1
                  !write(*,*) 'SET, ARRQ(1) = ',ARRQ(ICELL,1) 
                  IF ( .NOT.ISWET(INDBG) ) CYCLE
               ELSEIF ( IQ.GT.0 .AND. LVCELL.LT.LEVEL ) THEN
                   ! ICELL is a lower level cell containing (X,Y),
                   ! so needs to be refined, as long as it is wet, doesn't 
                   ! have any coarser neighbours, and hasn't already 
                   ! been selected
                  IF ( .NOT.ISWET(INDBG) ) CYCLE
                  NFINE = NFINE+1
                  NOREF = .FALSE.
                  ! Check that this isn't the central cell of a quad that
                  ! has already been refined
                  IF ( ISUB.EQ.0 ) THEN
                     !IQ = QTREE%INDQUAD(ICELL)
                     DO IS=1,4
                        IF ( QTREE%QCHILD(IQ,IS).GT.0 .OR.            &
                             QTREE%QICELL(IQ,IS).GT.0 ) THEN
                           !write(*,*) 'central cell of already refined quad'
                           NOREF = .TRUE.
                           EXIT
                        END IF
                     END DO
                     IF ( NOREF ) CYCLE
                  END IF
                  ! Check that this cell doesn't have any coarser neighbours,
                  ! and isn't of a type that shouldn't be refined:
                  IF ( ICELL.GT.0 ) THEN
                     DO IS=1,8
                        NB = QTREE%NGBR(ICELL,IS)
                        IF ( NB.GT.0 ) THEN
                           IF ( QTREE%INDLVL(NB).LT.LVCELL .AND.      &
                              QTREE%CELL_TYPE(NB).NE.UNDEF_TYPE ) THEN
                              !write(*,*) 'cell has coarser neighbours'
                              NOREF = .TRUE.
                              NCRSNB = NCRSNB + 1
                              EXIT
                           END IF
                        END IF
                     END DO
                     IF ( NOREF ) CYCLE
                     IF ( NT_NOADAPT.GT.0 ) THEN
                        IF ( ANY(NO_ADAPT_TYPE.EQ.                    &
                                 QTREE%CELL_TYPE(ICELL)) ) CYCLE
                     END IF
                  END IF
                  ! Check that this cell isn't already on the list of cells 
                  ! to be refined:
                  DO IREF=1,NREF
                     IF ( IQREF(IREF).EQ.IQ .AND.                     &
                          ISREF(IREF).EQ.ISUB ) THEN
                        !write(*,*) 'already on refinement list'
                        !NOREF = .TRUE.
                        NRPT = NRPT + 1
                        !EXIT
                     END IF
                  END DO
                  IF ( NOREF ) CYCLE
                  NREF = NREF + 1
                  IQREF(NREF) = IQ
                  ISREF(NREF) = ISUB
                  !write(*,*) 'REFINE'
               END IF
            END DO  ! Loop over cells at this level
            !
            write(*,*) 'QA_MLG2QT: INDBG1, INDBGMAX, ITER, NREF = ',  &
                                   INDBG1, INDBGMAX,ITER,NREF
            write(*,*) 'QA_MLG2QT: NWETALL, NINSUB, NSET, NFINE, NRPT, NCRSNB = ',  &
                        NWETALL, NINSUB, NSET, NFINE, NRPT, NCRSNB
            !
            ! If there is nothing to refine, no further iterations are required
            IF ( NREF.EQ.0 ) EXIT
            ! Refine the selected cells
            IMODE = 1
            CALL QA_REFINE ( IMODE, QTREE, NREF, IQREF, ISREF,         &
                             ISEAOLD, ISEANEW, ierr=IERS, ndse=IUN,    &
                             mapml=MAPML )
            IF ( IERS.NE.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE(IUN,*) 'ERROR IN QA_MLG2QT CALLING QA_REFINE'
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
            !do iref=1,NREF
            !   write(*,*) 'IREF, ISEAOLD, ISEANEW: ',IREF, ISEAOLD(IREF), ISEANEW(IREF,:)
            !end do
            !
            ! Recompute neighbours
            CALL QA_FINDNBR ( QTREE, ierr=IERS, ndse=IUN )
            IF ( IERS.NE.0 ) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE(IUN,*) 'ERROR IN QA_MLG2QT CALLING QA_FINDNBR'
               END IF
               IF ( PRESENT(IERR) ) IERR = IERS
               RETURN
            END IF
            !call qa_ioqt( 6, QTREE, 4, IERS, ndse=IUN )
         END DO  ! iteration loop
      END DO  ! Loop over refinement levels
!
      RETURN
      END SUBROUTINE QA_MLG2QT
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_MLGRID( MINDBG, NX, NY, LVLREL, LVLSG, LVLREF,    &
                            NX0, NY0, X0, Y0, INDBG0, ZBIN,           &
                            TRNX, TRNY, ZLIM, INDML, ISWET, XML, YML, &
                            ZBML, TRNXML, TRNYML, IERR, NDSE )
!/
!       Richard Gorman, NIWA
!         Oct, 2009:        Origination.
!         Jan-April 2014:   Rename
!
!  1. Purpose :
!
!      Precompute depth and transmission arrays, for all possible cells of
!      a multilevel bathymetry grid.
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       MINDBG  Int.   I   Max. number of points in the bathymetry 
!                          multi-grid
!       NX,NY   Int.   I   Size of input (sub)grid
!       LVLREL  Int.   I   Number of refinement levels below the input 
!                          subgrid to compute
!       LVLSG   Int.   I   Refinement level of the subgrid (relative 
!                          to reference level)
!       LVLREF  Int.   I   Reference refinement level
!       NX0,NY0 Int.   I   Size of level-zero grid
!       X0      Real   I   West wall coordinate of the input subgrid
!       Y0      Real   I   South wall coordinate of the input subgrid
!       INDBG0  Int.   IO  Present number of bathymetry grid points at 
!                          all levels
!       ZBIN    R.A.   I   Elevation of the seabed on fine grid 
!                          (m above datum)
!       TRNX,Y  R.A.   I   Transmission factors
!       ZLIM    Real   I   Limiting bottom depth, used to define land.
!       INDML   I.A.   IO  Indices in the full multilevel grid of 
!                          bathymetry grid cells
!       XML     R.A.   IO  X coordinates of (multilevel) bathymetry
!                          grid cells
!       YML     R.A.   IO  Y coordinates of (multilevel) bathymetry
!                          grid cells
!       ZBML    R.A.   IO  Elevation of the seabed on multilevel cells
!       TRNXML  R.A.   IO  X Transmission factors on multilevel cells
!       TRNYML  R.A.   IO  Y Transmission factors on multilevel cells
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name             Type  Module   Description
!     ---------------------------------------------------------------
!     MULTILEVEL_INDEX  Func. qa_utils Multilevel index of a given cell
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
      INTEGER, INTENT(IN)    :: MINDBG, NX, NY, LVLREL
      INTEGER, INTENT(IN)    :: LVLSG, LVLREF, NX0, NY0
      REAL, INTENT(IN)       :: X0, Y0
      INTEGER, INTENT(INOUT) :: INDBG0
      REAL,    INTENT(IN)    :: ZBIN(NX,NY), TRNX(NY,NX), TRNY(NY,NX)
      REAL,    INTENT(IN)    :: ZLIM
      INTEGER, INTENT(INOUT) :: INDML(MINDBG)
      LOGICAL, INTENT(INOUT) :: ISWET(MINDBG)
      REAL,    INTENT(INOUT) :: XML(MINDBG), YML(MINDBG), ZBML(MINDBG), &
                                TRNXML(MINDBG), TRNYML(MINDBG)
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!
! Local variables.
      INTEGER           :: NXSUB, NYSUB, LVL, II, JJ, INDBG, NXP, NYP
      INTEGER           :: INDBGMAX, IJS(4), IS, IBG, NWET
      INTEGER           :: INDBG1, IJT, LEVEL
      REAL              :: ZVAL, ZSUM, ZMAX, TXVAL(4), TYVAL(4)
      REAL              :: X, Y, XNP, YNP, DELXY
      INTEGER           :: IUN
      
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0

      XNP = FLOAT(NX)*2.**(-LVLREL)
      NXP = NINT(XNP)
      YNP = FLOAT(NY)*2.**(-LVLREL)
      NYP = NINT(YNP)
      IF ( ABS(XNP-NXP).GT.1.E-9 .OR. ABS(YNP-NYP).GT.1.E-9 ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_MLGRID: '
            WRITE(IUN,*) 'MISMATCH IN NXP, NYP'
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
      NXSUB = NX
      NYSUB = NY
      INDBG1 = INDBG0
      INDBG0 = INDBG0 + NXP*NYP*NINT( (4.**(LVLREL + 1) - 1)/3. )
      IF (INDBG0.GT.MINDBG) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_MLGRID: '
            WRITE(IUN,*) 'INDBG0 OUT OF RANGE: '
            WRITE(IUN,*) '       INDBG0: ', INDBG0
            WRITE(IUN,*) '       MINDBG: ', MINDBG
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
!
! Loop over level relative to the input subgrid:
!
      DO LVL=0,-LVLREL,-1
!
! Absolute refinement level:
         LEVEL = LVL + LVLSG + LVLREF
!
!       Highest multi-level bathy grid index for this level:
         INDBGMAX = NXP*NYP*NINT( (4.**(LVL + LVLREL +1) - 1)/3. )    &
                    + INDBG1
!
!       Highest multi-level bathy grid index for the next-lowest level:
         INDBG = NXP*NYP*NINT( (4.**(LVL + LVLREL) - 1)/3. ) + INDBG1
!
         DELXY = 2.**(-LVL-LVLSG)
         DO JJ=1,NYSUB
            Y = Y0 + (FLOAT(JJ) - 0.5)*DELXY
            DO II=1,NXSUB
               X = X0 + (FLOAT(II) - 0.5)*DELXY
               INDBG = INDBG + 1
               INDML(INDBG) = MULTILEVEL_INDEX ( NX0, NY0, X, Y,      &
                                                 LEVEL, LVLREF )
               XML(INDBG) = X
               YML(INDBG) = Y
               IF (LVL.EQ.0) THEN
                  ZBML(INDBG) = ZBIN(II,JJ)
                  TRNXML(INDBG) = TRNX(JJ,II)
                  TRNYML(INDBG) = TRNY(JJ,II)
                  ISWET(INDBG) = ZBIN(II,JJ).LE.ZLIM
               ELSE
!
!             Four higher-level cells within this cell:
!               (1: SW, 2: SE, 3: NW, 4: NE)
                  IBG = INDBGMAX + 2*NXSUB*(JJ-1) + 2*(II-1) + 1
                  IJS(1) = IBG
                  IJS(2) = IBG+1
                  IBG = IBG + 2*NXSUB*(JJ-1)
                  IJS(3) = IBG
                  IJS(4) = IBG+1
!
!             Average bed elevation over wet cells:
                  NWET = 0
                  ZSUM = 0.
                  ZMAX = 0.
                  DO IS=1,4
                     IJT = IJS(IS)
                     ZVAL = ZBML(IJT)
                     ZMAX = MAX(ZMAX, ZVAL)
                     TXVAL(IS) = TRNXML(IJT)
                     TYVAL(IS) = TRNYML(IJT)
                     IF ( ISWET(IJT) ) THEN
                        NWET = NWET + 1
                        ZSUM = ZSUM + ZVAL
                     ELSE
                        TXVAL(IS) = 0.
                        TYVAL(IS) = 0.
                     END IF
                  END DO
                  IF ( NWET.GT.0 ) THEN
                     ZBML(INDBG) = ZSUM/NWET
                     ISWET(INDBG) = .TRUE.
!
!            Compute transmission factors, taking products of successive cells in
!            the downwave direction and averaging in the cross-wave direction:
                     TRNXML(INDBG) = 0.5*( TXVAL(1)*TXVAL(2) +        &
                                           TXVAL(3)*TXVAL(4) )
                     TRNYML(INDBG) = 0.5*( TYVAL(1)*TYVAL(3) +        &
                                           TYVAL(2)*TYVAL(4) )
                  ELSE
!
!            If there are no wet cells, use the maximum for "mean" bed elevation:
                     ZBML(INDBG) = ZMAX
                     ISWET(INDBG) = .FALSE.
                     TRNXML(INDBG) = 0.
                     TRNYML(INDBG) = 0.
                  END IF
               END IF
            END DO
         END DO
         IF ( MOD(NXSUB,2).NE.0 ) EXIT
         IF ( MOD(NYSUB,2).NE.0 ) EXIT
         NXSUB = NXSUB/2
         NYSUB = NYSUB/2
      END DO
!
      RETURN
      END SUBROUTINE QA_MLGRID
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_NBDIST( QTREE, NCALC, ICALC, IERR, NDSE )
!/
!         Richard Gorman, NIWA
!           November, 2009: origination
!           Jan-April 2014: rewrite using QTREE structure
!
!  1. Purpose :
!
!     Compute indices into tables of values of weights for values and  
!     gradients at the centres and walls of the active cells, from values
!     at the cell centre (K=0) and its quadtree neighbours (K=1,...,8).
!     The tables have been precomputed for all possible configurations,
!     by QA_DERWTS
!
!  2. Method :
!
!
!  3. Parameters :
!
!     ----------------------------------------------------------------
!
!      QTREE  QA_TREE   I/O  Quadtree structure, including components
!                            affected by this routine:
!          NCELL   Int.  I  Number of cells
!          INDLVL  I.A.  I  Level of each cell
!          NGBR    I.A.  I  Sea-point indices for the neighbours of input cells
!                           NGBR(I,1:4) = index of primary (i.e. of equal or 
!                           lower level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!                           NGBR(I,5:8) = index of secondary (i.e. of higher 
!                           level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!          NCASE   I.A.  O  Index identifying the relative configuration of 
!                           neighbouring cells (out of 2401 possibilities)
!          XYVAL   R.A.  I  Coordinates of cells (w.r.t. 
!                             reference rectangular grid)
!       NCALC    Int.   I*  Number of cells to process (default = all)
!       ICALC    I.A.   I*  Array of cell indices to process
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
      TYPE(QA_TREE), INTENT(INOUT)  :: QTREE
      INTEGER, OPTIONAL, INTENT(IN) :: NCALC
      INTEGER, OPTIONAL, INTENT(IN) :: ICALC(:)
      INTEGER, OPTIONAL, INTENT(OUT) :: IERR
      INTEGER, OPTIONAL, INTENT(IN) :: NDSE
!
! Local variables:
!
      INTEGER               :: NSEEK, ISEEK, IUN
      INTEGER               :: ITYPE(4)
      INTEGER               :: ISEA, IDIR, INBR1, INBR2, LVL0,        &
                               LVLN, ILEV, IDISP, ICASE
      REAL                  :: XVAL0, YVAL0, XVALN, YVALN
!
      IF ( PRESENT(IERR) ) IERR = 0
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
!
      IF ( QTREE%IWTORDER.LT.1 ) THEN
         IF ( IUN.GT.0 )  THEN
            WRITE(IUN,*) 'ERROR IN QA_NBDIST: IWTORDER MUST BE >=1'
            WRITE(*,*) 'ERROR IN QA_NBDIST: IWTORDER MUST BE >=1'
            IF ( PRESENT(IERR) ) IERR = 1
            RETURN 
         END IF
      END IF
!
!   Determine which cells to process:
!   If ICALC and NCALC are specified, process the first NCALC indices
!   in the array ICALC.
!   If NCALC is provided but not ICALC, process the first NCALC cells
!   If NCALC is not provided, or is zero, process all cells
!
      IF ( PRESENT(NCALC) ) THEN
         IF ( NCALC.LE.0 ) THEN
            NSEEK = QTREE%NCELL
         ELSE
            NSEEK = MIN(NCALC, QTREE%NCELL)
         END IF
      ELSE
         NSEEK = QTREE%NCELL
      END IF
      IF ( PRESENT(ICALC) ) THEN
         NSEEK = MIN(SIZE(ICALC),NSEEK)
      END IF
!
!     Loop over cells
!
      DO ISEEK=1,NSEEK
        IF ( PRESENT(ICALC) )  THEN
          ISEA = ICALC(ISEEK)
        ELSE
          ISEA = ISEEK
        END IF
        IF ( QTREE%CELL_TYPE(ISEA).EQ.QTREE%UNDEF_TYPE ) CYCLE
!
!     Loop over directions [W, E, S, N]:
        DO IDIR=1,4
!
!     Evaluate the distribution of neighbours in this direction
          INBR1 = QTREE%NGBR(ISEA,IDIR)
          INBR2 = QTREE%NGBR(ISEA,IDIR+4)
          LVL0 = QTREE%INDLVL(ISEA)
          XVAL0 = QTREE%XYVAL(ISEA,1)
          YVAL0 = QTREE%XYVAL(ISEA,2)
!
!         Treat invalid neighbours as if non-existent:
          IF (INBR1 .GT. 0) THEN
            IF ( QTREE%CELL_TYPE(INBR1).EQ.QTREE%UNDEF_TYPE ) INBR1 = 0
          END IF
          IF (INBR2 .GT. 0) THEN
            IF ( QTREE%CELL_TYPE(INBR2).EQ.QTREE%UNDEF_TYPE ) INBR2 = 0
          END IF
          IF (INBR1 .GT. 0) THEN
            IF (INBR2 .GT. 0 ) THEN
!
!           Two neighbours
              ITYPE(IDIR) = 2
            ELSE
!
!           One neighbour:
              LVLN = QTREE%INDLVL(INBR1)
              IF (LVLN.EQ.LVL0) THEN
!
!           ... at same level:
                ITYPE(IDIR) = 1
              ELSE
!
!           ... at different level:
                ILEV = 3
                IF (LVLN.GT.LVL0) ILEV = 5
!
!           ... offset in + or - direction?:
                XVALN = QTREE%XYVAL(INBR1,1)
                YVALN = QTREE%XYVAL(INBR1,2)
                IDISP = 0
                IF (IDIR.LE.2) THEN
                  IF (YVALN.GT.YVAL0) IDISP = 1
                ELSE
                  IF (XVALN.GT.XVAL0) IDISP = 1
                END IF
                ITYPE(IDIR) = ILEV+IDISP
              END IF
            END IF
          ELSEIF (INBR2 .GT. 0 ) THEN
!
!           One neighbour:
            LVLN = QTREE%INDLVL(INBR2)
            IF (LVLN.EQ.LVL0) THEN
!
!           ... at same level:
              ITYPE(IDIR) = 1
            ELSE
!
!           ... at different level:
              ILEV = 3
              IF (LVLN.GT.LVL0) ILEV = 5
!
!           ... offset in + or - direction?:
              XVALN = QTREE%XYVAL(INBR2,1)
              YVALN = QTREE%XYVAL(INBR2,2)
              IDISP = 0
              IF (IDIR.LE.2) THEN
                IF (YVALN.GT.YVAL0) IDISP = 1
              ELSE
                IF (XVALN.GT.XVAL0) IDISP = 1
              END IF
              ITYPE(IDIR) = ILEV+IDISP
            END IF
          ELSE
!
!    No neighbours:
            ITYPE(IDIR) = 7
          END IF
        END DO
!
        ICASE = ITYPE(1)
        DO IDIR=2,4
           ICASE = 7*(ICASE-1) + ITYPE(IDIR)
        END DO
        QTREE%NCASE(ISEA) = ICASE
        !WRITE(*,*) 'isea, icase: ', isea, icase
      END DO
!
      END SUBROUTINE QA_NBDIST
!/
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_PRCELL( IUNS, QTREE, DEP )
!/
!         Richard Gorman, NIWA 
!          May, 2008:        Origination
!          Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose:
!
!     Print out the cell variables for a quadtree grid
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IUNS    Int.   I   Logical unit number of output file 
!                          (opened previously)
!       QTREE   QA_TREE I  Quadtree structure, including the following
!                             components used in this subroutine:
!         NSEA    Int.   I   Actual number of cells
!         XYVAL   R.A.   I   Coordinates of cells (w.r.t. fine 
!                          rectangular grid)
!         NGBR    I.A.   I   Sea-point indices (up to 8) for the 
!                          neighbours of input cells
!         INDLVL  I.A.   I   Level of each cell
!         MAPSTP  Int.   I   Flag for sea/boundary/excluded points
!       DEP     R.A.   I   Elevation of the seabed at each cell 
!                          (m above datum)
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!
      IMPLICIT NONE
      TYPE(QA_TREE), INTENT(IN)  :: QTREE
      INTEGER, INTENT(IN)  :: IUNS
      REAL, INTENT(IN)     :: DEP(:)
!
!  local variables
!
      INTEGER   :: ISEA, II
!
! Print the cell structure
      WRITE(IUNS,*) '   3   No. of lines in the header'
      WRITE(IUNS,*) '  14   No. of columns of data'
      WRITE(IUNS,*) ' ISEA    X       Y       DEPTH      LEVEL  MAPSTP | NGBR'
      DO ISEA=1,QTREE%NCELL
         IF ( QTREE%CELL_TYPE(ISEA).EQ.QTREE%UNDEF_TYPE ) CYCLE
         WRITE(IUNS,2100) ISEA, QTREE%XYVAL(ISEA,1),                  &
                          QTREE%XYVAL(ISEA,2), DEP(ISEA),             &
                          QTREE%INDLVL(ISEA), QTREE%CELL_TYPE(ISEA),  &
                          (QTREE%NGBR(ISEA,II),II=1,8)
      END DO
      CLOSE(IUNS)
      RETURN
 2100 FORMAT(I6,1X,2F8.2,F10.2,I6,1X,I6,1X,8I6)
!
      END SUBROUTINE QA_PRCELL
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_Q2NOKEEP ( QTREE1, NCELL2, QTREE2, IERR, NDSE )
!/
!         Richard Gorman, NIWA
!          April, 2014: Origination
!
!  1. Purpose :
!
!      Given a quadtree structure QTREE1, which retains refined cells,
!      find the number of cells in the corresponding structure (QTREE2)
!      that doesn't retain refined cells. Optionally, return QTREE2.
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE1  QA_TREE  I  Quadtree structure, that retains refined cells
!       NCELL2   Int.    O  Number of cells in the revised structure
!       QTREE2  QA_TREE  O* Quadtree structure, revised to not retain refined
!                            cells
!       IERR     Int.    O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*   Unit number for error output (if >0)
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(IN)            :: QTREE1
      INTEGER, INTENT(OUT)                 :: NCELL2
      TYPE(QA_TREE), OPTIONAL, INTENT(OUT) :: QTREE2
      INTEGER, OPTIONAL, INTENT(OUT)       :: IERR
      INTEGER, OPTIONAL, INTENT(IN)        :: NDSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: ICELL, ICELL2, ICELL3, NCELL, MCELL
      INTEGER         :: IUN
      INTEGER         :: II, NCKILL, IQ, ISHIFT, NCDEF, MN, MG, MCG,  &
                         MCN, MQUAD, NQUAD
      LOGICAL, ALLOCATABLE :: KILLCELL(:)
      INTEGER, ALLOCATABLE :: NEWCELL(:)
!/
!/ ------------------------------------------------------------------- /
!
      IF ( PRESENT(IERR) ) IERR = 0
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
!
      NCELL = QTREE1%NCELL
      NQUAD = QTREE1%NQUAD
!
      ALLOCATE( KILLCELL(NCELL) )
      KILLCELL = .FALSE.
      NCKILL = 0
!
!     Loop through quads to identify retained central cell indices
!
      DO IQ=1, NQUAD
         ICELL = QTREE1%QICELL(IQ,0)
         IF ( ICELL.GT.0 ) THEN
            DO II=1,4
               IF ( QTREE1%QICELL(IQ,II).GT.0 .OR.                    &
                    QTREE1%QCHILD(IQ,II).GT.0 ) THEN
                  KILLCELL(ICELL) = .TRUE.
                  NCKILL = NCKILL + 1
                  EXIT
               END IF
            END DO
         END IF
      END DO
!
!     Number of cells after these are removed:
!
      NCELL2 = NCELL - NCKILL
!
!     Return if no replacement QTREE required
!
      IF ( .NOT.PRESENT(QTREE2) ) RETURN
!
!     Check allocations
!
      MCELL = SIZE(QTREE2%INDQUAD,1)
      IF ( MCELL.LT.NCELL2 ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 )                                              &
            WRITE(IUN,*) 'QA_Q2NOKEEP: insufficient cell allocation'
            WRITE(IUN,*) 'ALLOCATED: ', MCELL
            WRITE(IUN,*) ' REQUIRED: ', NCELL2
         RETURN
      END IF
      MQUAD = SIZE(QTREE2%QLEVEL,1)
      IF ( MQUAD.LT.NQUAD ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 )                                              &
            WRITE(IUN,*) 'QA_Q2NOKEEP: insufficient quad allocation'
            WRITE(IUN,*) 'ALLOCATED: ', MQUAD
            WRITE(IUN,*) ' REQUIRED: ', NQUAD
         RETURN
      END IF
!
!     Loop through cells to determine cell index mappings
!
!
      ALLOCATE( NEWCELL(NCELL) )
      NEWCELL = 0
      ISHIFT = 0
      DO ICELL = 1, NCELL
         IF ( KILLCELL(ICELL) ) THEN
            ISHIFT = ISHIFT + 1
         ELSE
            NEWCELL(ICELL) = ICELL - ISHIFT
         END IF
      END DO
!
!     Set up new QTREE
!
      QTREE2%NQUAD = NQUAD
      QTREE2%NCELL = NCELL2
      QTREE2%LVLREF = QTREE1%LVLREF
      QTREE2%LVLMAX = QTREE1%LVLMAX
      QTREE2%LVLHI = QTREE1%LVLHI
      QTREE2%NX0 = QTREE1%NX0
      QTREE2%NY0 = QTREE1%NY0
      QTREE2%UNDEF_TYPE = QTREE1%UNDEF_TYPE
      QTREE2%KEEP_REF = .FALSE.
      QTREE2%DYNAMIC = QTREE1%DYNAMIC
      !QTREE2%IWTORDER = QTREE1%IWTORDER

      QTREE2%QICELL = 0
      QTREE2%QLEVEL = 0
      QTREE2%QPARENT = 0
      QTREE2%QNBR = 0
      QTREE2%QCHILD = 0
      DO IQ=1,NQUAD
         DO II=0,4
            ICELL = QTREE1%QICELL(IQ,II)
            IF ( ICELL.GT.0 ) THEN
               ICELL2 = NEWCELL(ICELL)
               QTREE2%QICELL(IQ,II) = ICELL2
            END IF
         END DO
         QTREE2%QLEVEL(IQ) = QTREE1%QLEVEL(IQ)
         QTREE2%QPARENT(IQ) = QTREE1%QPARENT(IQ)
         QTREE2%QNBR(IQ,:) = QTREE1%QNBR(IQ,:)
         QTREE2%QCHILD(IQ,:) = QTREE1%QCHILD(IQ,:)
      END DO
      
      MCN = SIZE(QTREE2%NCASE,1)
      MCG = SIZE(QTREE2%GNBR,1)
      MN = SIZE(QTREE1%NGBR,2)
      MN = MIN(MN,SIZE(QTREE2%NGBR,2))
      MG = SIZE(QTREE1%GNBR,2)
      MG = MIN(MG,SIZE(QTREE1%GNBR,2))
      QTREE2%INDQUAD = 0
      QTREE2%INDSUB = 0
      QTREE2%INDLVL = 0
      QTREE2%NGBR = 0
      QTREE2%CELL_TYPE = 0
      QTREE2%INDML = 0
      QTREE2%XYVAL = 0.
      IF ( MCN.GT.0 .AND. QTREE2%IWTORDER.GE.1 ) QTREE2%NCASE = 0
      IF ( MCG.GT.0 .AND. QTREE2%IWTORDER.GE.2 ) THEN
         QTREE2%INDWT = 0
         QTREE2%GNBR = 0
         QTREE2%NGNBR = 0
      END IF
      !
      NCDEF = 0
      DO ICELL=1,NCELL
         ICELL2 = NEWCELL(ICELL)
         IF ( ICELL2.EQ.0 ) CYCLE
         IF ( QTREE1%CELL_TYPE(ICELL).NE.QTREE1%UNDEF_TYPE )          &
            NCDEF = NCDEF + 1
         QTREE2%INDQUAD(ICELL2) = QTREE1%INDQUAD(ICELL)
         QTREE2%INDSUB(ICELL2) = QTREE1%INDSUB(ICELL)
         QTREE2%INDLVL(ICELL2) = QTREE1%INDLVL(ICELL)
         QTREE2%CELL_TYPE(ICELL2) = QTREE1%CELL_TYPE(ICELL)
         QTREE2%INDML(ICELL2) = QTREE1%INDML(ICELL)
         QTREE2%XYVAL(ICELL2,:) = QTREE1%XYVAL(ICELL,:)
         DO II=1,MN
            ICELL3 = QTREE1%NGBR(ICELL,II)
            IF ( ICELL3.GT.0 ) QTREE2%NGBR(ICELL2,II) = NEWCELL(ICELL3)
         END DO
         IF ( QTREE1%IWTORDER.GE.1 .AND. QTREE2%IWTORDER.GE.1 .AND.   &
              ICELL2.LE.MCN ) THEN
            QTREE2%NCASE(ICELL2) = QTREE1%NCASE(ICELL)
         END IF
         IF ( QTREE1%IWTORDER.GE.2 .AND. QTREE2%IWTORDER.GE.2 .AND.   &
              ICELL2.LE.MCG ) THEN
            QTREE2%INDWT(ICELL2) = QTREE1%INDWT(ICELL)
            QTREE2%NGNBR(ICELL2) = QTREE1%NGNBR(ICELL)
            DO II=1,MG
               ICELL3 = QTREE1%GNBR(ICELL,II)
               IF ( ICELL3.GT.0 ) QTREE2%GNBR(ICELL2,II) =            &
                                  NEWCELL(ICELL3)
            END DO
         END IF
      END DO
      QTREE2%NCELL_DEF = NCDEF


      RETURN
      END SUBROUTINE QA_Q2NOKEEP
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_Q2QMAP ( QTREE1, QTREE2, GET_TYPE, IND12,         &
                             EXACT12, NCALC, ICALC, IERR, NDSE,       &
                             IND21, EXACT21 )
!/
!         Richard Gorman, NIWA
!          February, 2014: Origination
!
!  1. Purpose :
!
!      Given two quadtree structures QTREE1 and QTREE2, identify cells in 
!      QTREE2 matching selected cells in QTREE1 (and vice versa) 
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE1  QA_TREE I/O  Quadtree structure, including the following
!                             components used in this subroutine:
!          NCELL     Int.  Number of cells
!          LVLREF    Int.  Refinement level of reference grid
!          NX0       Int.  X dimension of level-zero base grid
!          NY0       Int.  Y dimension of level-zero base grid
!          INDLVL    I.A.  Level of each cell
!          CELL_TYPE I.A.  Type of cell
!       QTREE2  QA_TREE  I   Quadtree structure, including the same 
!                             components as above used in this subroutine

!       GET_TYPE Log.    I   Flag to derive CELL_TYPE in QTREE1 from CELL_TYPE in QTREE2
!       IND12    I.A.   I/O  Cell indices in QTREE2 corresponding to cells in QTREE1
!       EXACT12  L.A.   I/O  Flag whether IND12 is an exact (equal level) match
!       NCALC    Int.    I*  Number of cells in QTREE1 to process (default = all)
!       ICALC    Int.    I*  Array of cell indices in QTREE1 to process
!       IERR     Int.    O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*   Unit number for error output (if >0)
!       IND21    I.A.   I/O* Cell indices in QTREE2 corresponding to cells in QTREE1
!       EXACT21  L.A.   I/O* Flag whether IND12 is an exact (equal level) match
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(INOUT)      :: QTREE1
      TYPE(QA_TREE), INTENT(IN)         :: QTREE2
      LOGICAL, INTENT(IN)               :: GET_TYPE
      INTEGER, INTENT(INOUT)            :: IND12(:)
      LOGICAL, OPTIONAL, INTENT(INOUT)  :: EXACT12(:)
      INTEGER, OPTIONAL, INTENT(IN)     :: NCALC
      INTEGER, OPTIONAL, INTENT(IN)     :: ICALC(:)
      INTEGER, OPTIONAL, INTENT(OUT)    :: IERR
      INTEGER, OPTIONAL, INTENT(IN)     :: NDSE
      INTEGER, OPTIONAL, INTENT(INOUT)  :: IND21(:)
      LOGICAL, OPTIONAL, INTENT(INOUT)  :: EXACT21(:)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: ICELL1, ICELL2, NCELL
      INTEGER         :: IQUAD, ISUB
      INTEGER         :: IUN
      INTEGER         :: ISEEK, NSEEK
      INTEGER         :: LVL1, LVL2
      INTEGER         :: ISTAT(2)
      REAL            :: X1, Y1, X2, Y2
!/
!/ ------------------------------------------------------------------- /
!
      IF ( PRESENT(IERR) ) IERR = 0
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
!
!   Check that the two quadtrees match at level 0:
!
      IF ( QTREE1%LVLREF .NE. QTREE2%LVLREF .OR.                      &
           QTREE1%NX0 .NE. QTREE2%NX0 .OR.                            &
           QTREE1%NY0 .NE. QTREE2%NY0 ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 )                                              &
            WRITE(IUN,*) 'QA_Q2QMAP: QTREE1 AND QTREE2 inconsistent'
         RETURN
      END IF
!
      NCELL = QTREE1%NCELL
!
!   Determine which cells to process:
!   If ICALC and NCALC are specified, process the first NCALC indices
!   in the array ICALC.
!   If NCALC is provided but not ICALC, process the first NCALC cells
!   If NCALC is not provided, or is zero, process all cells
!
      IF ( PRESENT(NCALC) ) THEN
         IF ( NCALC.LE.0 ) THEN
            NSEEK = NCELL
         ELSE
            NSEEK = MIN(NCALC, NCELL)
         END IF
      ELSE
         NSEEK = NCELL
      END IF
      IF ( PRESENT(ICALC) ) THEN
         NSEEK = MIN(SIZE(ICALC),NSEEK)
      END IF
!
      DO ISEEK=1, NSEEK
         IF ( PRESENT(ICALC) ) THEN
            ICELL1 = ICALC(ISEEK)
         ELSE
            ICELL1 = ISEEK
         END IF
         IF ( QTREE1%CELL_TYPE(ICELL1).EQ.QTREE1%UNDEF_TYPE ) CYCLE
         X1 = QTREE1%XYVAL(ICELL1,1)
         Y1 = QTREE1%XYVAL(ICELL1,2)
         LVL1 = QTREE1%INDLVL(ICELL1)
         CALL QA_XY2CELL ( QTREE2, X1, Y1, ICELL2, X2, Y2, LVL2,      &
                           IQUAD, ISUB, ISTAT )
         IF ( ISTAT(1).EQ.0 ) THEN
            IND12(ICELL1) = ICELL2
            IF ( PRESENT(EXACT12) ) EXACT12(ICELL1) = LVL2.EQ.LVL1
            IF ( GET_TYPE ) QTREE1%CELL_TYPE(ICELL1) =                &
                            QTREE2%CELL_TYPE(ICELL2) 
            IF ( PRESENT(IND21) ) IND21(ICELL2) = ICELL1 
            IF ( PRESENT(EXACT21) ) EXACT21(ICELL2) = LVL2.EQ.LVL1 
         ELSE
            IND12(ICELL1) = 0
            IF ( PRESENT(EXACT12) ) EXACT12(ICELL1) = .FALSE.
         END IF        
      END DO
      RETURN
      END SUBROUTINE QA_Q2QMAP
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_QMEAN ( QTREE, ARRAY, AMEAN, N_AVE, IQ_AVE,       &
                            A_UNDEF, EXC_MEAN, NMINMAX, NWETMIN,      &
                            IQUADMAX, IQUADMIN, EXC_MINMAX, ICMEAN,   &
                            WCMEAN )
!/
!         Richard Gorman, NIWA
!          April, 2008:      Origination
!          Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!      Determine the mean values of a variable over cells of each quad,
!      and/or the indices and weights needed to compute this
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE   QA_TREE I  Quadtree structure, including the following
!                             components used in this subroutine:
!          NQUAD      Int. I   Number of quads
!          UNDEF_TYPE Int. I   value of flag for unused cell indices
!          CELL_TYPE  I.A. I   Flag for wet/dry/boundary/etc. cell
!          QICELL     I.A. I   Sea-point indices for the cells of the quad
!       ARRAY    R.A.   I   Array values at cells
!       AMEAN    R.A.   O   Mean value of ARRAY over the cells 
!					in each quad
!       N_AVE    Int.   I*  Number of quad indices for which to compute averages
!                           [default if absent or zero: all quads]
!       IQ_AVE   I.A.   I*  Array of quad indices for which to compute averages
!       A_UNDEF  Real   I*  Value to assign AMEAN if undefined.
!       EXC_MEAN I.A.   I*  Array of values of cell_types to be excluded from 
!                           computing the means. [default: QTREE%UNDEF_TYPE]
!       NMINMAX  Int.  I/O* Number of values to compare with 
!                           for max/min
!       NWETMIN  Int.   I*  Min. number of 'valid' cells for quad to be
!				            included in min/max search
!       IQUADMAX I.A.   O*  Quad indices sorted in descending order by AMEAN
!       IQUADMIN I.A.   O*  Quad indices sorted in ascending order by AMEAN
!       EXC_MINMAX I.A. I*  Array of values of cell_types to be be considered 'invalid' 
!                           for the min/max search. [default: QTREE%UNDEF_TYPE]
!       ICMEAN   I.A.   O*  Array of (up to four) cell indices to
!                           averaged over for each output index
!       WCMEAN   R.A.   O*  Array of corresponding weights
!     ----------------------------------------------------------------
!              * = OPTIONAL
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_BMLSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC   Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(IN)        :: QTREE
      REAL, INTENT(IN), OPTIONAL       :: ARRAY(:)
      REAL, INTENT(OUT), OPTIONAL      :: AMEAN(:)
      INTEGER, INTENT(IN),OPTIONAL     :: N_AVE
      INTEGER, INTENT(IN),OPTIONAL     :: IQ_AVE(:)
      REAL, INTENT(IN), OPTIONAL       :: A_UNDEF
      INTEGER, INTENT(IN),OPTIONAL     :: EXC_MEAN(:)
      INTEGER, INTENT(INOUT),OPTIONAL  :: NMINMAX
      INTEGER, INTENT(IN),OPTIONAL     :: NWETMIN
      INTEGER, INTENT(OUT),OPTIONAL    :: IQUADMAX(:), IQUADMIN(:)
      INTEGER, INTENT(IN),OPTIONAL     :: EXC_MINMAX(:)
      INTEGER, INTENT(OUT), OPTIONAL   :: ICMEAN(:,:)
      REAL, INTENT(OUT), OPTIONAL      :: WCMEAN(:,:)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: ISEA, NCMP, II, ICMP, IQ, NSUM, NNOTBND
      INTEGER         :: IMAXSAVE, IMINSAVE, IMAXTEMP, IMINTEMP
      INTEGER         :: NSEEK, NWMIN, NEXC1, NEXC2, ITEST, NTEST
      REAL            :: ASUM, UNDEF
      LOGICAL         :: FOUNDB, FOUNDS, FINDMINMAX
      INTEGER         :: MSEA, MQUAD, NQUAD
      INTEGER         :: UNDEF_TYPE
      INTEGER, ALLOCATABLE :: MEAN_EXCLUDE_TYPE(:)
      INTEGER, ALLOCATABLE :: MINMAX_EXCLUDE_TYPE(:)
!/
!/ ------------------------------------------------------------------- /
!/
!
      MSEA = SIZE(ARRAY,1)
      MQUAD = SIZE(QTREE%QICELL,1)
!
      NQUAD = QTREE%NQUAD
      UNDEF_TYPE = QTREE%UNDEF_TYPE
      !
      IF ( PRESENT(N_AVE) ) THEN
         IF ( N_AVE.LE.0 ) THEN
            NTEST = NQUAD
         ELSE
            NTEST = MIN(N_AVE,NQUAD)
         END IF
      ELSE
         NTEST = NQUAD
      END IF
      IF ( PRESENT(IQ_AVE) ) THEN
         NTEST = MIN(NTEST,SIZE(IQ_AVE,1))
      END IF 
      IF ( PRESENT(A_UNDEF) ) THEN
         UNDEF = A_UNDEF
      ELSE
         UNDEF = 0.
      END IF 
            
      IF ( PRESENT(IQUADMAX) ) IQUADMAX = 0
      IF ( PRESENT(IQUADMIN) ) IQUADMIN = 0
      IF ( PRESENT(EXC_MEAN) ) THEN
         NEXC1 = SIZE(EXC_MEAN)
         ALLOCATE( MEAN_EXCLUDE_TYPE(NEXC1) )
         MEAN_EXCLUDE_TYPE = EXC_MEAN
      ELSE
         NEXC1 = 1
         ALLOCATE( MEAN_EXCLUDE_TYPE(NEXC1) )
         MEAN_EXCLUDE_TYPE(NEXC1) = UNDEF_TYPE
      END IF
      IF ( PRESENT(EXC_MINMAX) ) THEN
         NEXC2 = SIZE(EXC_MINMAX)
         ALLOCATE( MINMAX_EXCLUDE_TYPE(NEXC2) )
         MINMAX_EXCLUDE_TYPE = EXC_MINMAX
      ELSE
         NEXC2 = 1
         ALLOCATE( MINMAX_EXCLUDE_TYPE(NEXC2) )
         MINMAX_EXCLUDE_TYPE(NEXC2) = UNDEF_TYPE
      END IF
      !
      IMAXSAVE = 0
      IMINSAVE = 0
      NCMP = 0
      NSEEK = 0
      NWMIN = 1
      FINDMINMAX = .FALSE.
      IF ( PRESENT(NWETMIN) ) NWMIN = NWETMIN
      IF ( PRESENT(ARRAY) .AND. PRESENT(AMEAN) .AND.                  &
           PRESENT(NMINMAX) .AND.                                     &
           (PRESENT(IQUADMAX).OR.PRESENT(IQUADMIN)) ) THEN
         NSEEK = NMINMAX
         NSEEK = MIN(NQUAD,NSEEK)
         FINDMINMAX = .TRUE.
         IF ( PRESENT(IQUADMAX) ) NSEEK = MIN(NSEEK,SIZE(IQUADMAX,1))
         IF ( PRESENT(IQUADMIN) ) NSEEK = MIN(NSEEK,SIZE(IQUADMIN,1))
      END IF
      IF ( PRESENT(NMINMAX) ) NMINMAX = NSEEK
      IF ( PRESENT(ICMEAN) ) ICMEAN = 0
      IF ( PRESENT(WCMEAN) ) WCMEAN = 0.
      !DO IQ=1, NQUAD
      DO ITEST=1, NTEST
         IF ( PRESENT(IQ_AVE) ) THEN
            IQ = IQ_AVE(ITEST)
         ELSE
            IQ = ITEST
         END IF
!  Average the array at the home cell
         ASUM = 0.
         NSUM = 0
         NNOTBND = 0
         DO II=0,4
            ISEA = QTREE%QICELL(IQ,II)
            IF ( ISEA.LE.0 .OR. ISEA.GT.MSEA ) CYCLE
            IF ( ANY(QTREE%CELL_TYPE(ISEA).EQ.MEAN_EXCLUDE_TYPE) ) CYCLE
            IF ( PRESENT(ARRAY) ) ASUM = ASUM + ARRAY(ISEA)
            NSUM = NSUM + 1
            IF ( PRESENT(ICMEAN) ) ICMEAN(ITEST,NSUM) = ISEA
            IF ( PRESENT(WCMEAN) ) WCMEAN(ITEST,NSUM) = 1.
            IF (.NOT.ANY(QTREE%CELL_TYPE(ISEA).EQ.                    &
                         MINMAX_EXCLUDE_TYPE))  NNOTBND = NNOTBND + 1
            IF (II.EQ.0) EXIT
         END DO
         IF (NSUM.EQ.0) THEN
!
!  Can't compute/compare mean where there are no cells in the quad
            IF ( PRESENT(AMEAN) ) AMEAN(ITEST) = UNDEF
         ELSE
            IF ( PRESENT(AMEAN) ) AMEAN(ITEST) = ASUM/NSUM
            IF ( PRESENT(WCMEAN) ) WCMEAN(ITEST,:) = WCMEAN(ITEST,:)/NSUM
            !WRITE(*,*) 'IQ, AMEAN =', IQ, AMEAN(ITEST)
            IF (NSEEK.GT.0) THEN
               IF ( (PRESENT(IQUADMAX).OR.PRESENT(IQUADMIN)) .AND.    &
                    NNOTBND.GE.NWMIN ) THEN
!
!  Compare this value with the current lists of highest/lowest values
                  FOUNDB = .FALSE.
                  FOUNDS = .FALSE.
                  DO ICMP=1,NCMP
                     IF (PRESENT(IQUADMAX)) THEN
                        IMAXTEMP = IQUADMAX(ICMP)
                        IF(FOUNDB) THEN
                           IQUADMAX(ICMP) = IMAXSAVE
                        ELSEIF (AMEAN(ITEST).GT.AMEAN(IMAXTEMP)) THEN
                           IQUADMAX(ICMP) = IQ
                           FOUNDB = .TRUE.
                        END IF
                        IMAXSAVE = IMAXTEMP
                     END IF
                     !
                     IF (PRESENT(IQUADMIN)) THEN
                        IMINTEMP = IQUADMIN(ICMP)
                        IF(FOUNDS) THEN
                           IQUADMIN(ICMP) = IMINSAVE
                        ELSEIF (AMEAN(ITEST).LT.AMEAN(IMINTEMP)) THEN
                           IQUADMIN(ICMP) = IQ
                           FOUNDS = .TRUE.
                        END IF
                        IMINSAVE = IMINTEMP
                     END IF
                  END DO
                  IF (NCMP.LT.NSEEK) THEN
                     NCMP = NCMP + 1
                     IF (PRESENT(IQUADMAX)) THEN
                        IF (FOUNDB) THEN
                           IQUADMAX(NCMP) = IMAXSAVE
                        ELSE
                           IQUADMAX(NCMP) = IQ
                        END IF
                     END IF
                     IF (PRESENT(IQUADMIN)) THEN
                        IF (FOUNDS) THEN
                           IQUADMIN(NCMP) = IMINSAVE
                        ELSE
                           IQUADMIN(NCMP) = IQ
                        END IF
                     END IF
                  END IF
               END IF
            END IF
         END IF
      END DO
      RETURN
!/
!/ End of QA_QMEAN ----------------------------------------------------- /
!
      END SUBROUTINE QA_QMEAN
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_QT2SMC ( QTREE, LBC, IJKCel, IJKUFc, IJKVFc,       &
                       NLvUFc, NLvVFc, NLvUFc2, NLvVFc2, IERR, NDSE )
!/
!       Richard Gorman, NIWA
!         October, 2017:      Origination.
!
!  1. Purpose :
!
!      Create SMC arrays from a threaded quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE   QA_TREE I/O  Quadtree structure, including the following
!                             components used in this subroutine:
!          LVLREF     Int.  Refinement level of reference grid
!          LVLMAX     Int.  Maximum refinement level allowed
!          UNDEF_TYPE Int.  Value of flag for unused cell indices
!          INDLVL     I.A.  Level of each cell
!          NGBR       I.A.  Sea-point indices for the neighbours of input cells
!                           NGBR(I,1:4) = index of primary (i.e. of equal or 
!                           lower level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!                           NGBR(I,5:8) = index of secondary (i.e. of higher 
!                           level) W,E,S,N neighbours, respectively, 
!                           if any, of cell I
!          CELL_TYPE  I.A.  Flag for wet/dry/boundary/etc. cell
!          XYVAL      R.A.  Coordinates of cells (w.r.t. 
!                                reference rectangular grid)
!          NCELL      Int.  Number of cells
!          UNDEF_TYPE Int.  Value of CELL_TYPE representing undefined cells
!       LBC      Int.   I   Lower bound for SMC cell indices (to allow for
!                           indices <= 0 representing extra dummy cells)
!       IJKCel   I.A.   O   SMC array containing cell locations and sizes
!       IJKUFc   I.A.   O   SMC array containing U face locations, sizes, and 
!                           indices (into IJKCel) of neighbouring cells 
!                           The array is sorted from finest to coarsest
!       IJKVFc   I.A.   O   SMC array containing V face locations, sizes, and 
!                           indices (into IJKCel) of neighbouring cells 
!                           The array is sorted from finest to coarsest
!       NLvUFc   I.A.   O*  SMC array containing the number of U faces at 
!                           or above each refinement level, ie. there are 
!                           NLvUFc(k) U faces with a face refinement level >=
!                           L = highest refinement level of adjoining cells
!                             = Max. refinement level - k + 1
!       NLvVFc   I.A.   O*  as above for V faces
!       NLvUFc2  I.A.   O*  NLvUFc2(k)  = No. of U faces with both adjoining
!                           cells at level >= L
!                           with L = Max. refinement level - k +1
!       NLvVFc2  I.A.   O"  as above for V faces
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(IN)       :: QTREE
      INTEGER, INTENT(IN)             :: LBC
      INTEGER, INTENT(OUT)            :: IJKCel(:,LBC:)
      INTEGER, INTENT(OUT)            :: IJKUFc(:,:)
      INTEGER, INTENT(OUT)            :: IJKVFc(:,:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: NLvUFc(:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: NLvVFc(:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: NLvUFc2(:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: NLvVFc2(:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: NCELL
      INTEGER         :: ICELL, IDIR
      INTEGER         :: LVLMAX, LVLREF, LVLHOME, LVLNBR
      INTEGER         :: IUN
      INTEGER         :: ID_REF
      INTEGER         :: IC_HOME, JC_HOME, ID_HOME, JD_HOME
      INTEGER         :: IC_FACE, JC_FACE, ID_FACE, JD_FACE
      INTEGER         :: IC_NBR, JC_NBR, ID_NBR, JD_NBR
      INTEGER         :: IOPP, INBR1, INBR2, INBR, IDIRSAV
      INTEGER         :: IP2, IP1, I00, IM1
      INTEGER         :: IFW, IFE, IFS, IFN
      INTEGER         :: IFW1, IFE1, IFS1, IFN1
      INTEGER         :: IFW2, IFE2, IFS2, IFN2
      INTEGER, PARAMETER  :: IUDRY=-9, IVDRY=-9
      INTEGER         :: NUFACE, NVFACE
      INTEGER         :: MUFACE1, MVFACE1, MUFACE2, MVFACE2, LI, NN
      INTEGER         :: IDONE, IND
      INTEGER, ALLOCATABLE  :: IFACE(:,:) 
      REAL            :: OFFSET
      INTEGER, ALLOCATABLE  :: NLISTU1(:), NLISTU2(:)
      INTEGER, ALLOCATABLE  :: NLISTV1(:), NLISTV2(:)
      INTEGER, ALLOCATABLE  :: LISTU1(:,:), LISTU2(:,:)
      INTEGER, ALLOCATABLE  :: LISTV1(:,:), LISTV2(:,:)
      INTEGER, ALLOCATABLE  :: ITEMP(:,:)
!/
!/ ------------------------------------------------------------------- /
!
      IF ( PRESENT(IERR) ) IERR = 0
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
!
      NCELL = QTREE%NCELL
!
! Check the output array allocations
!
      IF ( UBOUND(IJKCel,2).LT.NCELL ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 ) WRITE(IUN,*)                                  &
              'QA_QT2SMC: ERROR IN IJKCel ALLOCATION(2)'
         RETURN
      END IF
      IF ( UBOUND(IJKCel,1).LT.4 ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 ) WRITE(IUN,*)                                  &
              'QA_QT2SMC: ERROR IN IJKCel ALLOCATION(1)'
         RETURN
      END IF
      MUFACE1 = UBOUND(IJKUFc,1)
      MUFACE2 = UBOUND(IJKUFc,2)
      IF ( MUFACE1.LT.7 ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 ) WRITE(IUN,*)                                  &
              'QA_QT2SMC: ERROR IN IJKUFc ALLOCATION(1)'
         RETURN
      END IF
      MVFACE1 = UBOUND(IJKVFc,1)
      MVFACE2 = UBOUND(IJKVFc,2)
      IF ( MVFACE1.LT.7 ) THEN
         IF ( PRESENT(IERR) ) IERR = 1
         IF ( IUN.GT.0 ) WRITE(IUN,*)                                  &
              'QA_QT2SMC: ERROR IN IJKVFc ALLOCATION(1)'
         RETURN
      END IF
!
! Initialise outputs
      IJKCel = 0
      IJKUFc = 0
      IJKVFc = 0
!
! Initialise counters and working array of face indices
      NUFACE = 0
      NVFACE = 0
      ALLOCATE ( IFACE(NCELL,8) )
      IFACE = 0
!
      LVLMAX = QTREE%LVLMAX
      LVLREF = QTREE%LVLREF
!
! Initialise sorting lists
      ALLOCATE ( NLISTU1(LVLMAX+1), NLISTU2(LVLMAX+1),                 &
                 NLISTV1(LVLMAX+1), NLISTV2(LVLMAX+1) ) 
      ALLOCATE ( LISTU1(MUFACE2,LVLMAX+1) )
      ALLOCATE ( LISTU2(MUFACE2,LVLMAX+1) )
      ALLOCATE ( LISTV1(MVFACE2,LVLMAX+1) )
      ALLOCATE ( LISTV2(MVFACE2,LVLMAX+1) )
      NLISTU1 = 0
      NLISTU2 = 0
      NLISTV1 = 0
      NLISTV2 = 0
!
! First loop over cells, to populate:
!  1. the cell array IJKCel, containing cell locations and sizes
!  2. Parts of the face arrays IJKUFc and IJKVFc containing
!     the face locations, sizes and first order neighbours. 
!     Default values only will be given for the second order neighbours
!  3. a temporary array IFACE giving the indices of faces (up to two 
!     in each direction) of each cell
!
!     Cell locations are for the cell's SW corner, in SMC integer
!     coordinates, in units of the  size of a cell at the highest
!     refinement level
!
!     Size of reference cell relative to finest cell:
      ID_REF = NINT(2.**(LVLMAX - LVLREF))
      DO ICELL=1,NCELL
!
!  No further processing if ICELL is an undefined cell:
!
         IF ( QTREE%CELL_TYPE(ICELL).EQ.QTREE%UNDEF_TYPE ) CYCLE
!
!  Properties of this cell:
         ! refinement level:
         LVLHOME = QTREE%INDLVL(ICELL)
         ! cell size relative to the finest:
         ID_HOME = NINT(2.**(LVLMAX -LVLHOME))
         JD_HOME = ID_HOME
         ! SMC coordinates (location of the SW corner in multiples
         ! of the cell size at highest refinement level):
         OFFSET = 0.5*(1. + 2.**(LVLREF-LVLHOME))
         IC_HOME = NINT((QTREE%XYVAL(ICELL,1) - OFFSET)*ID_REF)
         JC_HOME = NINT((QTREE%XYVAL(ICELL,2) - OFFSET)*ID_REF)
         IJKCel(1,ICELL) = IC_HOME
         IJKCel(2,ICELL) = JC_HOME
         IJKCel(3,ICELL) = ID_HOME
         IJKCel(4,ICELL) = JD_HOME
         !write(*,*) 'QA_QT2SMC: ICELL, LVLHOME, XYVAL = ',             &
         !    icell, lvlhome, QTREE%XYVAL(ICELL,:)

!  Loop over W, E directions: 
         DO IDIR=1,2
            IF (IDIR.EQ.1) THEN
              IP2 = 4
              IP1 = 5
              I00 = 6
              IM1 = 7
              IOPP = 2
              IC_FACE = IC_HOME
            ELSE
              IM1 = 4
              I00 = 5
              IP1 = 6
              IP2 = 7
              IOPP = 1
              IC_FACE = IC_HOME + ID_HOME
            END IF
            INBR1 = QTREE%NGBR(ICELL,IDIR)
            INBR2 = QTREE%NGBR(ICELL,IDIR+4)
            IF ( INBR1.EQ.0 .AND. INBR2.EQ.0 ) THEN
               ! No neighbouring wet cells in this direction
               ! Create a new (external boundary) face
               NUFACE = NUFACE + 1
               IF ( MUFACE2.LT.NUFACE ) THEN
                  IF ( PRESENT(IERR) ) IERR = 1
                  IF ( IUN.GT.0 ) WRITE(IUN,*)                         &
                     'QA_QT2SMC: ERROR IN IJKUFc ALLOCATION(1)'
                  RETURN
               END IF
               IFACE(ICELL,IDIR) = NUFACE
               IJKUFc(1,NUFACE) = IC_FACE
               IJKUFc(2,NUFACE) = JC_HOME 
               IJKUFc(3,NUFACE) = JD_HOME
               ! Cells each side of the face:
               IJKUFc(IP1,NUFACE) = IUDRY
               IJKUFc(I00,NUFACE) = ICELL
               ! Cells one further step from the face:
               IJKUFc(IP2,NUFACE) = IUDRY
               IJKUFc(IM1,NUFACE) = ICELL  ! to be overwritten
               !
               ! Add this face to the list of U faces adjoining (only)
               ! cells of level LVLHOME
               LI = LVLMAX - LVLHOME + 1
               NLISTU1(LI) = NLISTU1(LI) +1
               IF ( MUFACE2.LT.NLISTU1(LI) ) THEN
                  IF ( PRESENT(IERR) ) IERR = 1
                  IF ( IUN.GT.0 ) WRITE(IUN,*)                         &
                     'QA_QT2SMC: ERROR IN LISTU1 ALLOCATION(1)',       &
                        'LI, NLISTU1, MUFACE2: ',                      &
                         LI, NLISTU1(LI), MUFACE2
                  RETURN
               END IF
               LISTU1(NLISTU1(LI),LI) = NUFACE
            ELSE IF ( INBR1.GT.0 .AND. INBR2.GT.0 ) THEN
               ! Two neighbouring wet cells in this direction
               IF ( INBR1.GT.ICELL ) THEN
                  ! This face hasn't been created yet
                  NUFACE = NUFACE + 1
                  IF ( MUFACE2.LT.NUFACE ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN IJKUFc ALLOCATION(2)',    &
                        'NUFACE, MUFACE2: ',                           &
                         NUFACE, MUFACE2
                     RETURN
                  END IF
                  IFACE(ICELL,IDIR) = NUFACE
                  !
                  ! Populate the face array
                  IJKUFc(1,NUFACE) = IC_FACE
                  IJKUFc(2,NUFACE) = JC_HOME 
                  IJKUFc(3,NUFACE) = JD_HOME/2
                  ! Cells each side of the face:
                  IJKUFc(IP1,NUFACE) = INBR1
                  IJKUFc(I00,NUFACE) = ICELL
                  ! Cells one further step from the face:
                  IJKUFc(IP2,NUFACE) = INBR1 ! to be overwritten
                  IJKUFc(IM1,NUFACE) = ICELL  ! to be overwritten
                  !
                  ! Add this face to the list of U faces adjoining
                  ! both cells of level LVLHOME+1 and LVLHOME
                  LI = LVLMAX - LVLHOME
                  NLISTU2(LI) = NLISTU2(LI) +1
                  IF ( MUFACE2.LT.NLISTU2(LI) ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN LISTU2 ALLOCATION(1)',    &
                        'LI, NLISTU2, MUFACE2: ',                      &
                         LI, NLISTU2(LI), MUFACE2
                     RETURN
                  END IF
                  LISTU2(NLISTU2(LI),LI) = NUFACE
               ELSE
                  ! This face was previously created from the neighbour's side:
                  IF ( QTREE%NGBR(INBR1,IOPP) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR1,IOPP)
                  ELSE IF ( QTREE%NGBR(INBR1,IOPP+4) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR1,IOPP+4)
                  END IF
               END IF
               IF ( INBR2.GT.ICELL ) THEN
                  ! This face hasn't been created yet
                  NUFACE = NUFACE + 1
                  IF ( MUFACE2.LT.NUFACE ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN IJKUFc ALLOCATION(2)',    &
                        'NUFACE, MUFACE2: ',                           &
                         NUFACE, MUFACE2
                     RETURN
                  END IF
                  IFACE(ICELL,IDIR+4) = NUFACE
                  !
                  ! Populate the face array
                  IJKUFc(1,NUFACE) = IC_FACE
                  IJKUFc(2,NUFACE) = JC_HOME + JD_HOME/2
                  IJKUFc(3,NUFACE) = JD_HOME/2
                  ! Cells each side of the face:
                  IJKUFc(IP1,NUFACE) = INBR2
                  IJKUFc(I00,NUFACE) = ICELL
                  ! Cells one further step from the face:
                  IJKUFc(IP2,NUFACE) = INBR2 ! to be overwritten
                  IJKUFc(IM1,NUFACE) = ICELL  ! to be overwritten
                  !
                  ! Add this face to the list of U faces adjoining
                  ! both cells of level LVLHOME+1 and LVLHOME
                  LI = LVLMAX - LVLHOME
                  NLISTU2(LI) = NLISTU2(LI) +1
                  IF ( MUFACE2.LT.NLISTU2(LI) ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN LISTU2 ALLOCATION(2)',    &
                        'LI, NLISTU2, MUFACE2: ',                      &
                         LI, NLISTU2(LI), MUFACE2
                     RETURN
                  END IF
                  LISTU2(NLISTU2(LI),LI) = NUFACE
               ELSE
                  ! This face was previously created from the neighbour's side:
                  IF ( QTREE%NGBR(INBR2,IOPP) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR2,IOPP)
                  ELSE IF ( QTREE%NGBR(INBR2,IOPP+4) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR2,IOPP+4)
                  END IF
               END IF
            ELSE
               ! One neighbouring wet cell in this direction
               IF ( INBR2.GT.0 ) THEN
                  INBR = INBR2
                  IDIRSAV = IDIR + 4
               ELSE
                  INBR = INBR1
                  IDIRSAV = IDIR
               END IF
               IF ( INBR.GT.ICELL ) THEN
                  ! This face hasn't been created yet
                  NUFACE = NUFACE + 1
                  IF ( MUFACE2.LT.NUFACE ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN IJKUFc ALLOCATION(2)',    &
                           'NUFACE, MVFACE2: ',                        &
                            NUFACE, MUFACE2
                     RETURN
                  END IF
                  IFACE(ICELL,IDIRSAV) = NUFACE
                  ! refinement level of the neighbour:
                  LVLNBR = QTREE%INDLVL(INBR)
                  ! cell size of the neighbour:
                  JD_NBR = NINT(2.**(LVLMAX -LVLHOME))
                  ! CFS coordinates of the neighbour:
                  JC_NBR = NINT((QTREE%XYVAL(INBR,2) - 1.)*JD_NBR)
                  ! wall length is the smaller of home and neighbour cell sizes
                  IF ( LVLNBR.GT.LVLHOME ) THEN
                     ! neighbour is half the size of the home cell:
                     ! use neighbour's coordinate and face length:
                     JC_FACE = JC_NBR
                     JD_FACE = JD_NBR
                     !
                     ! Add this face to the list of U faces adjoining
                     ! both cells of level LVLHOME+1 and LVLHOME
                     LI = LVLMAX - LVLHOME
                     NLISTU2(LI) = NLISTU2(LI) + 1
                     IF ( MUFACE2.LT.NLISTU2(LI) ) THEN
                        IF ( PRESENT(IERR) ) IERR = 1
                        IF ( IUN.GT.0 ) WRITE(IUN,*)                   &
                           'QA_QT2SMC: ERROR IN LISTU2 ALLOCATION(3)', &
                           'LI, NLISTU2, MUFACE2: ',                   &
                            LI, NLISTU2(LI), MUFACE2
                        RETURN
                     END IF
                     LISTU2(NLISTU2(LI),LI) = NUFACE
                  ELSE IF ( LVLNBR.EQ.LVLHOME ) THEN
                     ! neighbour is the same size as the home cell
                     ! use home cell's coordinate and face length:
                     JC_FACE = JC_HOME
                     JD_FACE = JD_HOME
                     !
                     ! Add this face to the list of U faces adjoining (only)
                     ! cells of level LVLHOME
                     LI = LVLMAX - LVLHOME + 1
                     NLISTU1(LI) = NLISTU1(LI) + 1
                     IF ( MUFACE2.LT.NLISTU1(LI) ) THEN
                        IF ( PRESENT(IERR) ) IERR = 1
                        IF ( IUN.GT.0 ) WRITE(IUN,*)                   &
                           'QA_QT2SMC: ERROR IN LISTU1 ALLOCATION(2)', &
                           'LI, NLISTU1, MVFACE2: ',                   &
                            LI, NLISTU1(LI), MUFACE2
                        RETURN
                     END IF
                     LISTU1(NLISTU1(LI),LI) = NUFACE
                  ELSE
                     ! neighbour is twice the size of the home cell
                     ! use home cell's coordinate and face length:
                     JC_FACE = JC_HOME
                     JD_FACE = JD_HOME
                     !
                     ! Add this face to the list of U faces adjoining
                     ! both cells of level LVLHOME and LVLHOME-1
                     LI = LVLMAX - LVLHOME + 1
                     NLISTU2(LI) = NLISTU2(LI) + 1
                     IF ( MUFACE2.LT.NLISTU2(LI) ) THEN
                        IF ( PRESENT(IERR) ) IERR = 1
                        IF ( IUN.GT.0 ) WRITE(IUN,*)                   &
                           'QA_QT2SMC: ERROR IN LISTU2 ALLOCATION(4)', &
                           'LI, NLISTU2, MVFACE2: ',                   &
                            LI, NLISTU2(LI), MUFACE2
                        RETURN
                     END IF
                     LISTU2(NLISTU2(LI),LI) = NUFACE
                  ENDIF
                  IJKUFc(1,NUFACE) = IC_FACE
                  IJKUFc(2,NUFACE) = JC_FACE
                  IJKUFc(3,NUFACE) = JD_FACE
                  ! Cells each side of the face:
                  IJKUFc(IP1,NUFACE) = INBR
                  IJKUFc(I00,NUFACE) = ICELL
                  ! Cells one further step from the face:
                  IJKUFc(IP2,NUFACE) = INBR ! to be overwritten
                  IJKUFc(IM1,NUFACE) = ICELL ! to be overwritten
               ELSE
                  ! This face was previously created from the neighbour's side:
                  IF ( QTREE%NGBR(INBR,IOPP) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIRSAV) = IFACE(INBR,IOPP)
                  ELSE IF ( QTREE%NGBR(INBR,IOPP+4) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIRSAV) = IFACE(INBR,IOPP+4)
                  END IF
                  !IF ( IFACE(INBR,IOPP) .GT. 0 ) THEN
                  !   IFACE(ICELL,IDIR) = IFACE(INBR,IOPP)
                  !ELSE
                  !   IFACE(ICELL,IDIR) = IFACE(INBR,IOPP+4)
                  !END IF
               END IF
            END IF
         END DO  ! End of loop over W, E directions
!
!  Loop over S, N directions: 
         DO IDIR=3,4
            IF (IDIR.EQ.3) THEN
              IP2 = 4
              IP1 = 5
              I00 = 6
              IM1 = 7
              IOPP = 4
              JC_FACE = JC_HOME
            ELSE
              IM1 = 4
              I00 = 5
              IP1 = 6
              IP2 = 7
              IOPP = 3
              JC_FACE = JC_HOME + JD_HOME
            END IF
            INBR1 = QTREE%NGBR(ICELL,IDIR)
            INBR2 = QTREE%NGBR(ICELL,IDIR+4)
            IF ( INBR1.EQ.0 .AND. INBR2.EQ.0 ) THEN
               ! No neighbouring wet cells in this direction
               ! Create a new (external boundary) face
               NVFACE = NVFACE + 1
               IF ( MVFACE2.LT.NVFACE ) THEN
                 IF ( PRESENT(IERR) ) IERR = 1
                  IF ( IUN.GT.0 ) WRITE(IUN,*)                         &
                     'QA_QT2SMC: ERROR IN IJKVFc ALLOCATION(2)',       &
                           'NVFACE, MVFACE2: ', NVFACE, MVFACE2
                  RETURN
               END IF
               IFACE(ICELL,IDIR) = NVFACE
               IJKVFc(1,NVFACE) = IC_HOME
               IJKVFc(2,NVFACE) = JC_FACE 
               IJKVFc(3,NVFACE) = ID_HOME
               ! Cells each side of the face:
               IJKVFc(IP1,NVFACE) = IVDRY
               IJKVFc(I00,NVFACE) = ICELL
               ! Cells one further step from the face:
               IJKVFc(IP2,NVFACE) = IVDRY
               IJKVFc(IM1,NVFACE) = ICELL  ! to be overwritten
               !
               ! Add this face to the list of V faces adjoining (only)
               ! cells of level LVLHOME
               LI = LVLMAX - LVLHOME + 1
               NLISTV1(LI) = NLISTV1(LI) + 1
               IF ( MVFACE2.LT.NLISTV1(LI) ) THEN
                  IF ( PRESENT(IERR) ) IERR = 1
                  IF ( IUN.GT.0 ) WRITE(IUN,*)                         &
                     'QA_QT2SMC: ERROR IN LISTV1 ALLOCATION(1)',       &
                           'LI, NLISTV1, MVFACE2: ',                   &
                            LI, NLISTV1(LI), MVFACE2
                  RETURN
               END IF
               LISTV1(NLISTV1(LI),LI) = NVFACE
            ELSE IF ( INBR1.GT.0 .AND. INBR2.GT.0 ) THEN
               ! Two neighbouring wet cells in this direction
               IF ( INBR1.GT.ICELL ) THEN
                  ! This face hasn't been created yet
                  NVFACE = NVFACE + 1
                  IF ( MVFACE2.LT.NVFACE ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN IJKVFc ALLOCATION(2)',    &
                           'NVFACE, MVFACE2: ', NVFACE, MVFACE2
                     RETURN
                  END IF
                  IFACE(ICELL,IDIR) = NVFACE
                  !
                  ! Populate the face array
                  IJKVFc(1,NVFACE) = IC_HOME
                  IJKVFc(2,NVFACE) = JC_FACE 
                  IJKVFc(3,NVFACE) = ID_HOME/2
                  ! Cells each side of the face:
                  IJKVFc(IP1,NVFACE) = INBR1
                  IJKVFc(I00,NVFACE) = ICELL
                  ! Cells one further step from the face:
                  IJKVFc(IP2,NVFACE) = INBR1 ! to be overwritten
                  IJKVFc(IM1,NVFACE) = ICELL  ! to be overwritten
                  !
                  ! Add this face to the list of V faces adjoining
                  ! both cells of level LVLHOME+1 and LVLHOME
                  LI = LVLMAX - LVLHOME
                  NLISTV2(LI) = NLISTV2(LI) + 1
                  IF ( MVFACE2.LT.NLISTV2(LI) ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN LISTV2 ALLOCATION(1)',    &
                           'LI, NLISTV2, MVFACE2: ',                   &
                            LI, NLISTV2(LI), MVFACE2
                     RETURN
                  END IF
                  LISTV2(NLISTV2(LI),LI) = NVFACE
               ELSE
                  ! This face was previously created from the neighbour's side:
                  IF ( QTREE%NGBR(INBR1,IOPP) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR1,IOPP)
                  ELSE IF ( QTREE%NGBR(INBR1,IOPP+4) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR1,IOPP+4)
                  END IF
               END IF
               IF ( INBR2.GT.ICELL ) THEN
                  ! This face hasn't been created yet
                  NVFACE = NVFACE + 1
                  IF ( MVFACE2.LT.NVFACE ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN IJKVFc ALLOCATION(2)'
                     RETURN
                  END IF
                  IFACE(ICELL,IDIR+4) = NVFACE
                  !
                  ! Populate the face array
                  IJKVFc(1,NVFACE) = IC_HOME + ID_HOME/2
                  IJKVFc(2,NVFACE) = JC_FACE 
                  IJKVFc(3,NVFACE) = ID_HOME/2
                  ! Cells each side of the face:
                  IJKVFc(IP1,NVFACE) = INBR2
                  IJKVFc(I00,NVFACE) = ICELL
                  ! Cells one further step from the face:
                  IJKVFc(IP2,NVFACE) = INBR2 ! to be overwritten
                  IJKVFc(IM1,NVFACE) = ICELL  ! to be overwritten
                  !
                  ! Add this face to the list of V faces adjoining
                  ! both cells of level LVLHOME+1 and LVLHOME
                  LI = LVLMAX - LVLHOME
                  NLISTV2(LI) = NLISTV2(LI) + 1
                  IF ( MVFACE2.LT.NLISTV2(LI) ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN LISTV2 ALLOCATION(2)',    &
                           'LI, NLISTV2, MVFACE2: ',                   &
                            LI, NLISTV2(LI), MVFACE2
                     RETURN
                  END IF
                  LISTV2(NLISTV2(LI),LI) = NVFACE
               ELSE
                  ! This face was previously created from the neighbour's side:
                  IF ( QTREE%NGBR(INBR2,IOPP) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR2,IOPP)
                  ELSE IF ( QTREE%NGBR(INBR2,IOPP+4) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIR) = IFACE(INBR2,IOPP+4)
                  END IF
               END IF
            ELSE
               ! One neighbouring wet cell in this direction
               IF ( INBR2.GT.0 ) THEN
                  INBR = INBR2
                  IDIRSAV = IDIR + 4
               ELSE
                  INBR = INBR1
                  IDIRSAV = IDIR
               END IF
               IF ( INBR.GT.ICELL ) THEN
                  ! This face hasn't been created yet
                  NVFACE = NVFACE + 1
                  IF ( MVFACE2.LT.NVFACE ) THEN
                     IF ( PRESENT(IERR) ) IERR = 1
                     IF ( IUN.GT.0 ) WRITE(IUN,*)                      &
                        'QA_QT2SMC: ERROR IN IJKVFc ALLOCATION(2)'
                     RETURN
                  END IF
                  IFACE(ICELL,IDIRSAV) = NVFACE
                  ! refinement level of the neighbour:
                  LVLNBR = QTREE%INDLVL(INBR)
                  ! cell size of the neighbour:
                  ID_NBR = NINT(2.**(LVLMAX -LVLHOME))
                  ! CFS coordinates of the neighbour:
                  IC_NBR = NINT((QTREE%XYVAL(INBR,2) - 1.)*ID_NBR)
                  ! wall length is the smaller of home and neighbour cell sizes
                  IF ( LVLNBR.GT.LVLHOME ) THEN
                     ! neighbour is half the size of the home cell:
                     ! use neighbour's coordinate and face length:
                     IC_FACE = IC_NBR
                     ID_FACE = ID_NBR
                     !
                     ! Add this face to the list of V faces adjoining
                     ! both cells of level LVLHOME+1 and LVLHOME
                     LI = LVLMAX - LVLHOME
                     NLISTV2(LI) = NLISTV2(LI) + 1
                     IF ( MVFACE2.LT.NLISTV1(LI) ) THEN
                        IF ( PRESENT(IERR) ) IERR = 1
                        IF ( IUN.GT.0 ) WRITE(IUN,*)                   &
                           'QA_QT2SMC: ERROR IN LISTV2 ALLOCATION(3)', &
                           'LI, NLISTV2, MVFACE2: ',                   &
                            LI, NLISTV2(LI), MVFACE2
                        RETURN
                     END IF
                     LISTV2(NLISTV2(LI),LI) = NVFACE
                  ELSE IF ( LVLNBR.EQ.LVLHOME ) THEN
                     ! neighbour is the same size as the home cell
                     ! use home cell's coordinate and face length:
                     IC_FACE = IC_HOME
                     ID_FACE = ID_HOME
                     !
                     ! Add this face to the list of V faces adjoining (only)
                     ! cells of level LVLHOME
                     LI = LVLMAX - LVLHOME + 1
                     NLISTV1(LI) = NLISTV1(LI) + 1
                     IF ( MVFACE2.LT.NLISTV1(LI) ) THEN
                        IF ( PRESENT(IERR) ) IERR = 1
                        IF ( IUN.GT.0 ) WRITE(IUN,*)                   &
                           'QA_QT2SMC: ERROR IN LISTV1 ALLOCATION(2)', &
                           'LI, NLISTV1, MVFACE2: ',                   &
                            LI, NLISTV1(LI), MVFACE2
                        RETURN
                     END IF
                     LISTV1(NLISTV1(LI),LI) = NVFACE
                  ELSE
                     ! neighbour is twice the size of the home cell
                     ! use home cell's coordinate and face length:
                     IC_FACE = IC_HOME
                     ID_FACE = ID_HOME
                     !
                     ! Add this face to the list of V faces adjoining
                     ! both cells of level LVLHOME and LVLHOME-1
                     LI = LVLMAX - LVLHOME + 1
                     NLISTV2(LI) = NLISTV2(LI) + 1
                     IF ( MVFACE2.LT.NLISTV2(LI) ) THEN
                        IF ( PRESENT(IERR) ) IERR = 1
                        IF ( IUN.GT.0 ) WRITE(IUN,*)                   &
                           'QA_QT2SMC: ERROR IN LISTV2 ALLOCATION(4)', &
                           'LI, NLISTV2, MVFACE2: ',                   &
                            LI, NLISTV2(LI), MVFACE2
                        RETURN
                     END IF
                     LISTV2(NLISTV2(LI),LI) = NVFACE
                  ENDIF
                  IJKVFc(1,NVFACE) = IC_FACE
                  IJKVFc(2,NVFACE) = JC_FACE
                  IJKVFc(3,NVFACE) = ID_FACE
                  ! Cells each side of the face:
                  IJKVFc(IP1,NVFACE) = INBR
                  IJKVFc(I00,NVFACE) = ICELL
                  ! Cells one further step from the face:
                  IJKVFc(IP2,NVFACE) = INBR ! to be overwritten
                  IJKVFc(IM1,NVFACE) = ICELL ! to be overwritten
               ELSE
                  ! This face was previously created from the neighbour's side:
                  IF ( QTREE%NGBR(INBR,IOPP) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIRSAV) = IFACE(INBR,IOPP)
                  ELSE IF ( QTREE%NGBR(INBR,IOPP+4) .EQ. ICELL ) THEN
                     IFACE(ICELL,IDIRSAV) = IFACE(INBR,IOPP+4)
                  END IF
                  !IF ( IFACE(INBR,IOPP) .GT. 0 ) THEN
                  !   IFACE(ICELL,IDIR) = IFACE(INBR,IOPP)
                  !ELSE
                  !   IFACE(ICELL,IDIR) = IFACE(INBR,IOPP+4)
                  !END IF
               END IF
            END IF
         END DO  ! End of loop over directions
      END DO  ! End of first loop over  sea points
!
! Second loop over cells, to get second-order neighbouring cells
! for each face
!
      DO ICELL=1, NCELL
!
!  No further processing if ICELL is an undefined cell:
!
         IF ( QTREE%CELL_TYPE(ICELL).EQ.QTREE%UNDEF_TYPE ) CYCLE
         !
         ! West face(s)
         ! If there is only one, store its index as IFW, otherwise IFW=0
         IFW1 = IFACE(ICELL,1)
         IFW2 = IFACE(ICELL,5)
         IF ( IFW1.GT.0 .AND. IFW2.EQ. 0 ) THEN
            IFW = IFW1
         ELSEIF ( IFW1.EQ.0 .AND. IFW2.GT. 0 ) THEN
            IFW = IFW2
         ELSE
            IFW = 0
         END IF
         ! East face(s)
         ! If there is only one, store its index as IFE, otherwise IFE=0
         IFE1 = IFACE(ICELL,2)
         IFE2 = IFACE(ICELL,6)
         IF ( IFE1.GT.0 .AND. IFE2.EQ. 0 ) THEN
            IFE = IFE1
         ELSEIF ( IFE1.EQ.0 .AND. IFE2.GT. 0 ) THEN
            IFE = IFE2
         ELSE
            IFE = 0
         END IF
         !
         ! If there is only one east face, the cell on its east side
         ! is the 2nd-order east cell for the western face(s)
         IF ( IFE.GT.0 ) THEN
            IF ( IFW1.GT.0 ) IJKUFc(7,IFW1) = IJKUFc(6,IFE)
            IF ( IFW2.GT.0 ) IJKUFc(7,IFW2) = IJKUFc(6,IFE)
         END IF
         !
         ! If there is only one west face, the cell on its west side
         ! is the 2nd-order west cell for the eastern face(s)
         IF ( IFW.GT.0 ) THEN
            IF ( IFE1.GT.0 ) IJKUFc(4,IFE1) = IJKUFc(5,IFW)
            IF ( IFE2.GT.0 ) IJKUFc(4,IFE2) = IJKUFc(5,IFW)
         END IF
         !
         ! South face(s)
         IFS1 = IFACE(ICELL,3)
         IFS2 = IFACE(ICELL,7)
         IF ( IFS1.GT.0 .AND. IFS2.EQ. 0 ) THEN
            IFS = IFS1
         ELSEIF ( IFS1.EQ.0 .AND. IFS2.GT. 0 ) THEN
            IFS = IFS2
         ELSE
            IFS = 0
         END IF
         ! North face(s)
         IFN1 = IFACE(ICELL,4)
         IFN2 = IFACE(ICELL,8)
         IF ( IFN1.GT.0 .AND. IFN2.EQ. 0 ) THEN
            IFN = IFN1
         ELSEIF ( IFN1.EQ.0 .AND. IFN2.GT. 0 ) THEN
            IFN = IFN2
         ELSE
            IFN = 0
         END IF
         !
         ! If there is only one north face, the cell on its north side
         ! is the 2nd-order north cell for the southern face(s)
         IF ( IFN.GT.0 ) THEN
            IF ( IFS1.GT.0 ) IJKVFc(7,IFS1) = IJKVFc(6,IFN)
            IF ( IFS2.GT.0 ) IJKVFc(7,IFS2) = IJKVFc(6,IFN)
         END IF
         !
         ! If there is only one south face, the cell on its south side
         ! is the 2nd-order south cell for the northern face(s)
         IF ( IFS.GT.0 ) THEN
            IF ( IFN1.GT.0 ) IJKVFc(4,IFN1) = IJKVFc(5,IFS)
            IF ( IFN2.GT.0 ) IJKVFc(4,IFN2) = IJKVFc(5,IFS)
         END IF
      END DO  ! End of loop over sea points
!
!     Sort the face arrays in descending order of:
!       1. Level of the face = higher of the adjoining cell levels, then
!       2. Lower of the adjoining cell levels 
!
      ALLOCATE ( ITEMP(MAX(MUFACE1,MVFACE1),MAX(MUFACE2,MVFACE2)) )
      IDONE = 0
      DO LI=1,LVLMAX+1
         DO NN = 1, NLISTU1(LI)
            IND = LISTU1(NN,LI)
            ITEMP(1:MUFACE1,IDONE+NN) = IJKUFc(1:MUFACE1,IND) 
         END DO
         IDONE = IDONE + NLISTU1(LI)
         DO NN = 1, NLISTU2(LI)
            IND = LISTU2(NN,LI)
            ITEMP(1:MUFACE1,IDONE+NN) = IJKUFc(1:MUFACE1,IND) 
         END DO
         IDONE = IDONE + NLISTU2(LI)
      END DO
      IJKUFc(1:MUFACE1,1:NUFACE) = ITEMP(1:MUFACE1,1:NUFACE)
!
      IDONE = 0
      DO LI=1,LVLMAX+1
         DO NN = 1, NLISTV1(LI)
            IND = LISTV1(NN,LI)
            ITEMP(1:MVFACE1,IDONE+NN) = IJKVFc(1:MVFACE1,IND) 
         END DO
         IDONE = IDONE + NLISTV1(LI)
         DO NN = 1, NLISTV2(LI)
            IND = LISTV2(NN,LI)
            ITEMP(1:MVFACE1,IDONE+NN) = IJKVFc(1:MVFACE1,IND) 
         END DO
         IDONE = IDONE + NLISTV2(LI)
      END DO
      IJKVFc(1:MVFACE1,1:NVFACE) = ITEMP(1:MVFACE1,1:NVFACE)
!
!     Optional output arguments:
      IF ( PRESENT(NLvUFc) ) THEN
         NLvUFC = 0
         DO LI = 1, MIN(LVLMAX+1,UBOUND(NLvUFc,1))
            IF ( LI.EQ.1 ) THEN
               NLvUFc(LI) = NLISTU1(LI) + NLISTU2(LI) 
            ELSE
               NLvUFc(LI) = NLvUFc(LI-1) + NLISTU1(LI) + NLISTU2(LI) 
            END IF
         END DO
      END IF
      IF ( PRESENT(NLvVFc) ) THEN
         NLvVFC = 0
         DO LI = 1, MIN(LVLMAX+1,UBOUND(NLvVFc,1))
            IF ( LI.EQ.1 ) THEN
               NLvVFc(LI) = NLISTV1(LI) + NLISTV2(LI) 
            ELSE
               NLvVFc(LI) = NLvVFc(LI-1) + NLISTV1(LI) + NLISTV2(LI) 
            END IF
         END DO
      END IF
      IF ( PRESENT(NLvUFc2) ) THEN
         NLvUFC2 = 0
         DO LI = 1, MIN(LVLMAX+1,UBOUND(NLvUFc2,1))
            IF ( LI.EQ.1 ) THEN
               NLvUFc2(LI) = 1
            ELSEIF ( LI.EQ.2 ) THEN
               NLvUFc2(LI) = 1 + NLISTU1(LI-1)
            ELSE
               NLvUFc2(LI) = NLvUFc2(LI-1) + NLISTU1(LI-1)             &
                                           + NLISTU2(LI-2)
            END IF
         END DO
      END IF
      IF ( PRESENT(NLvVFc2) ) THEN
         NLvVFC2 = 0
         DO LI = 1, MIN(LVLMAX+1,UBOUND(NLvVFc2,1))
            IF ( LI.EQ.1 ) THEN
               NLvVFc2(LI) = 1
            ELSEIF ( LI.EQ.2 ) THEN
               NLvVFc2(LI) = 1 + NLISTV1(LI-1)
            ELSE
               NLvVFc2(LI) = NLvVFc2(LI-1) + NLISTV1(LI-1)             &
                                           + NLISTV2(LI-2)
            END IF
         END DO
      END IF
!
      RETURN
!/
!/ End of QT2SMC ----------------------------------------------------- /
!/
      END SUBROUTINE QA_QT2SMC
!/
!
      SUBROUTINE QA_QTCHECK ( QTREE, IERR, NDSE )
!/
!       Richard Gorman, NIWA
!         October, 2017:      Origination.
!
!  1. Purpose :
!
!      Consistency checks on a threaded quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE   QA_TREE I/O  Quadtree structure, including the following
!                             components used in this subroutine:
!           NCELL
!           n
!       IERR     Int.   O   Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(IN)       :: QTREE
      INTEGER, INTENT(OUT)            :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: NCELL
      INTEGER         :: ICELL, IDIR, IR, INBREV, INBREV1, INBREV2
      INTEGER         :: LVLHOME, LVLNBR, LVLMAX
      INTEGER         :: IUN
      INTEGER         :: INBR1, INBR2, INBR
      INTEGER         :: LVL2, IQ, IS, ICELL2

      INTEGER         :: IREV(4)
!/
!/ ------------------------------------------------------------------- /
!
      IERR = 0
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IREV(1) = 2
      IREV(2) = 1
      IREV(3) = 4
      IREV(4) = 3
!
      NCELL = QTREE%NCELL
      LVLMAX = QTREE%LVLMAX
!
      DO ICELL=1,NCELL
!
!  No further processing if ICELL is an undefined cell:
!
         !IF ( QTREE%CELL_TYPE(ICELL).EQ.QTREE%UNDEF_TYPE ) CYCLE
!
!  Properties of this cell:
         ! refinement level:
         LVLHOME = QTREE%INDLVL(ICELL)
         IF ( LVLHOME .GT. LVLMAX ) THEN
            IERR = 1
            IF (IUN.GT.0) WRITE (IUN,1000) ICELL, LVLHOME, LVLMAX
         END IF
!
!  Neighbour checks:
!  Loop over W, E, S, N directions: 
         DO IDIR=1,4
            INBR1 = QTREE%NGBR(ICELL,IDIR)
            INBR2 = QTREE%NGBR(ICELL,IDIR+4)
            IF ( INBR1.GT.0 .AND. INBR2.GT.0 ) THEN
               ! home cell has two neighbours in this direction:
               ! they should both have the home cell as the sole
               ! neighbour in the reverse direction
               ! and should both be at refinement level +1 higher
               ! First check the primary neighbour:
               LVLNBR = QTREE%INDLVL(INBR1)
               IF ( LVLNBR-LVLHOME .NE. 1 ) THEN
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1002)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                       INBR2, LVLNBR, LVLHOME+1
               END IF
               !
               IR = IREV(IDIR)
               INBREV = QTREE%NGBR(INBR1,IR)
               IF ( INBREV .NE. ICELL ) THEN
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1001)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                       INBR1, IR, INBREV, ICELL
               END IF
               IR = IREV(IDIR) + 4
               INBREV = QTREE%NGBR(INBR1,IR)
               !
               IF ( INBREV .NE. 0 ) THEN
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1001)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                       INBR1, IR, INBREV, 0
               END IF
               ! Now check the secondary neighbour:
               LVLNBR = QTREE%INDLVL(INBR2)
               IF ( LVLNBR-LVLHOME .NE. 1 ) THEN
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1002)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                       INBR2, LVLNBR, LVLHOME+1
               END IF
               !
               IR = IREV(IDIR)
               INBREV = QTREE%NGBR(INBR2,IR)
               IF ( INBREV .NE. ICELL ) THEN
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1001)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                       INBR2, IR, INBREV, ICELL
               END IF
               IR = IREV(IDIR) + 4
               INBREV = QTREE%NGBR(INBR2,IR)
               IF ( INBREV .NE. 0 ) THEN
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1001)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                       INBR1, IR, INBREV, 0
               END IF
            ELSE IF ( ( INBR1.GT.0 .AND. INBR2.EQ.0 ) .OR.             &
                      ( INBR1.EQ.0 .AND. INBR2.GT.0 ) ) THEN
               ! home cell has one neighbour in this direction:
               INBR = MAX ( INBR1, INBR2 )
               LVLNBR = QTREE%INDLVL(INBR)
               !
               IF ( LVLNBR-LVLHOME .EQ. 1 .OR. LVLNBR.EQ.LVLHOME ) THEN
                  ! It is at refinement level +1 higher: 
                  ! OK (e.g. because the would-be other neighbour is dry) 
                  ! or the same refinement level
                  ! It should have the home cell as the sole
                  ! (primary) neighbour in the reverse direction
                  IR = IREV(IDIR)
                  INBREV = QTREE%NGBR(INBR,IR)
                  IF ( INBREV .NE. ICELL ) THEN
                     IERR = 1
                     IF (IUN.GT.0) WRITE (IUN,1001)                       &
                          ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                          INBR, IR, INBREV, ICELL
                  END IF
                  IR = IREV(IDIR) + 4
                  INBREV = QTREE%NGBR(INBR,IR)
                  !
                  IF ( INBREV .NE. 0 ) THEN
                     IERR = 1
                     IF (IUN.GT.0) WRITE (IUN,1001)                       &
                          ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                          INBR, IR, INBREV, 0
                  END IF
               ELSE IF ( LVLNBR-LVLHOME .EQ. -1 ) THEN
                  ! It is at refinement level 1 lower: 
                  ! It should have the home cell as one of its
                  ! neighbours in the reverse direction
                  IR = IREV(IDIR)
                  INBREV1 = QTREE%NGBR(INBR,IR)
                  INBREV2 = QTREE%NGBR(INBR,IR+4)
                  IF ( INBREV1 .NE. ICELL .AND. INBREV2.NE.ICELL ) THEN
                     IERR = 1
                     IF (IUN.GT.0) WRITE (IUN,1003)                    &
                          ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,  &
                          INBR, IR, INBREV1, IR+4, INBREV2, ICELL
                  END IF
               ELSE
                  ! Illegal difference  in refinement levels between
                  !
                  IERR = 1
                  IF (IUN.GT.0) WRITE (IUN,1004)                       &
                       ICELL, LVLHOME, IDIR, INBR1, IDIR+4, INBR2,     &
                          INBR, LVLNBR
               END IF
            END IF
         END DO  ! End of loop over W, E, S, N directions
!
!  Check quad & subquad indices
         IQ = QTREE%INDQUAD(ICELL)
         IS = QTREE%INDSUB(ICELL)
         ICELL2 = QTREE%QICELL(IQ,IS)
         IF ( ICELL2 .NE. ICELL ) THEN
            IERR = 1
            IF (IUN.GT.0) WRITE (IUN,1005) ICELL, IQ, IS,              &
                                           IQ, IS, ICELL2
         END IF
         LVL2 = QTREE%QLEVEL(IQ)
         IF ( LVL2 .NE. LVLHOME .AND. IS.GT.0 ) THEN
            IERR = 1
            IF (IUN.GT.0) WRITE (IUN,1006) ICELL, LVLHOME, IQ,         &
                                           IQ, LVL2
         END IF
      END DO  ! End of loop
!
      RETURN
!
! Formats
!
 1000 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS LEVEL = ',I2/                         &
               '        BUT LVLMAX = ',I6 )
 1001 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS LEVEL = ',I2,' NGBR(',I1,') = ',I6/   &
               '                        AND NGBR(',I1,') = ',I6/       &
               '        BUT CELL ',I6,' HAS NGBR(',I1,') = ',I6/       &
               '                          SHOULD BE = ',I6 )
 1002 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS LEVEL = ',I2,' NGBR(',I1,') = ',I6/   &
               '                        AND NGBR(',I1,') = ',I6/       &
               '          BUT CELL ',I6,' HAS LEVEL = ',I6/            &
               '                          SHOULD BE = ',I6 )
 1003 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS LEVEL = ',I2,' NGBR(',I1,') = ',I6/   &
               '                        AND NGBR(',I1,') = ',I6/       &
               '        BUT CELL ',I6,' HAS NGBR(',I1,') = ',I6/       &
               '                        AND NGBR(',I1,') = ',I6/       &
               '             ONE OF THESE SHOULD BE = ',I6 )
 1004 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS LEVEL = ',I2,' NGBR(',I1,') = ',I6/   &
               '                        AND NGBR(',I1,') = ',I6/       &
               '        BUT CELL ',I6,' HAS LEVEL = ',I2/              &
               '             DIFFERENCE SHOULD BE = -1, 0 OR +1' )
 1005 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS INDQUAD = ',I6,' INDSUB = ',I2/       &
               '          BUT QICELL(',I6,',',I2,') = ',I6)
 1006 FORMAT ( ' QA_QTCHECK: ERROR: '/                                 &
               ' CELL ',I6,' HAS LEVEL = ',I2,' INDQUAD = 'I6/         &
               '        BUT QLEVEL(',I6,') HAS LEVEL = ',I2)
!/
!/ End of QTCHECK ----------------------------------------------------- /
!/
      END SUBROUTINE QA_QTCHECK
!/
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_READWT( IUN, NDSO, NDSE, QWEIGHT )
!/
!         Richard Gorman, NIWA
!          November, 2009:   Origination
!          Jan-April 2014:   Rewrite using QWEIGHT structure
!
!  1. Purpose :
!
!      Reads in arrays of precomputed interpolation weights for all
!      possible arrangements of quadtree neighbours
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IUN     Int.   I   Unit number for weights file
!       NDSO    Int.   I   Unit number for standard output
!       NDSE    Int.   I   Unit number for standard error
!       QWEIGHT QA_WEIGHT O  Structure containing weights, with the 
!                            following components:
!         MCASE        Int.  Allocated number of possible cases (2401)
!         NNGBR        Int.  Maximum number of neighbouring cells (8)
!         WTSREAD      L.A.  Flags for whether arrays have been read fully:
!                           WTSREAD(1) = .TRUE. if centre gradient weights read
!                           WTSREAD(2) = .TRUE. if wall value weights read
!                           WTSREAD(3) = .TRUE. if wall gradient weights read
!                           WTSREAD(4) = .TRUE. if centre 2nd deriv. weights read
!         NCASEQUAD    Int.  Number of precomputed cell configurations
!         IQUADDIR     I.A   Code for type of neighbour(s)
!         XGRADWTC     R.A   Weights for X-gradient at cell centre
!         YGRADWTC     R.A   Weights for Y-gradient at cell centre
!         XGRADWTW     R.A   Weights for X-gradient at west wall
!         XGRADWTE     R.A   Weights for X-gradient at east wall
!         YGRADWTS     R.A   Weights for Y-gradient at south wall
!         YGRADWTN     R.A   Weights for Y-gradient at north wall
!         LAMBDA0      R.A   Natural neighbour coordinates relative 
!                                 to cell centre
!         VALWTW       R.A   Weights for value at west wall
!         VALWTE       R.A   Weights for value at east wall
!         VALWTS       R.A   Weights for value at south wall
!         VALWTN       R.A   Weights for value at north wall
!         XCOORDQUAD   R.A   X-coordinates of target and neighbour 
!                              cell centres
!         YCOORDQUAD   R.A   Y-coordinates of target and neighbour
!                              cell centres
!         XXDER2WTC    R.A   Weights for d2/dXdX at cell centre
!         XYDER2WTC    R.A   Weights for d2/dXdY at cell centre
!         YYDER2WTC    R.A   Weights for d2/dYdY at cell centre
!
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)          :: IUN
      INTEGER, INTENT(IN)          :: NDSO, NDSE
      TYPE(QA_WEIGHT), INTENT(OUT) :: QWEIGHT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, PARAMETER         :: MCASE=2401, NNGBR=8
      INTEGER                    :: NCASEQUAD
      INTEGER   :: NREAD(16), IDIR(4)
      INTEGER   :: ILIN, NHDR, NARRAY, NCOL, IARR, ICIN, ITYP,         &
                   IRD, ICASE, I, IK
      REAL      :: VALS(0:NNGBR)
!/
!/ ------------------------------------------------------------------- /
!/
!
      QWEIGHT%WTSREAD = .FALSE.
!
      ILIN = 0
!
!  Read the first header line: contains the total number of header lines
      READ(IUN,*,END=901,ERR=902) NHDR
      ILIN = ILIN+1
!
!  Read the second header line: contains the number of cases
      READ(IUN,*,END=901,ERR=902) NCASEQUAD
      ILIN = ILIN+1
      IF (NCASEQUAD.GT.MCASE) THEN
         WRITE(NDSO,*) 'QA_READWT: Insufficient Case Nos. Allocated: '
         WRITE(NDSO,*) '         Allocated: ', MCASE
         WRITE(NDSO,*) '           In file: ', NCASEQUAD
      END IF
!
!  Read the third header line: contains the number of input arrays
      READ(IUN,*,END=901,ERR=902) NARRAY
      ILIN = ILIN+1
!
!  Read the fourth header line: contains the number of columns
      READ(IUN,*,END=901,ERR=902) NCOL
      ILIN = ILIN+1
      IF (NCOL.LT.NNGBR+7) THEN
         WRITE(NDSO,*) 'QA_READWT: WARNING: Unexpected NCOL: '
         WRITE(NDSO,*) '              Read: ', NCOL
         WRITE(NDSO,*) '               Max: ', NNGBR+7
      END IF
!
!  Read the rest of the header lines:
      DO IRD=1,NHDR-ILIN
        READ(IUN,*,END=901,ERR=902)
      END DO
!
!  Initialise weight arrays to zero:
      !ICASEQUAD = 0
      QWEIGHT%IQUADDIR = 0
      QWEIGHT%XGRADWTC = 0.
      QWEIGHT%YGRADWTC = 0.
      QWEIGHT%XGRADWTW = 0.
      QWEIGHT%XGRADWTE = 0.
      QWEIGHT%YGRADWTS = 0.
      QWEIGHT%YGRADWTN = 0.
      QWEIGHT%LAMBDA0 = 0.
      QWEIGHT%VALWTW = 0.
      QWEIGHT%VALWTE = 0.
      QWEIGHT%VALWTS = 0.
      QWEIGHT%VALWTN = 0.
      QWEIGHT%XCOORDQUAD = 0.
      QWEIGHT%YCOORDQUAD = 0.
      QWEIGHT%XXDER2WTC = 0.
      QWEIGHT%XYDER2WTC = 0.
      QWEIGHT%YYDER2WTC = 0.
!
!     Initialise counter for number of cases read into each array:
      NREAD = 0
!
!  Loop over expected number of input arrays:
      DO IARR=1,NARRAY
!
!  Loop over expected number of cases:
        DO ICASE=1,NCASEQUAD
           READ(IUN,*,END=901,ERR=902) ICIN, (IDIR(I),I=1,4),       &
                                       ITYP, (VALS(IK),IK=0,NNGBR)
           IF (ICIN.LT.1 .OR. ICIN.GT.NCASEQUAD) THEN
             WRITE(NDSO,*) 'QA_READWT: WARNING: Out-of-range Case No.: '
             WRITE(NDSO,*) '              Read: ', ICIN
             WRITE(NDSO,*) '               Max: ', NCASEQUAD
             CYCLE
           ELSEIF (ICIN.NE.ICASE) THEN
             WRITE(NDSO,*) 'QA_READWT: WARNING: Unexpected Case No.: '
             WRITE(NDSO,*) '              Read: ', ICIN
             WRITE(NDSO,*) '          Expected: ', ICASE
           END IF
           IF (ITYP.GE.1 .AND. ITYP.LE.16) NREAD(ITYP) = NREAD(ITYP)+1
           IF (ICIN.GT.MCASE) CYCLE
           !ICASEQUAD(ICIN) = ICIN
           IF (ICASE.EQ.1) THEN
             QWEIGHT%IQUADDIR(ICIN,:) = IDIR
           END IF
           IF (ITYP.EQ.1) THEN
             ! Weights for X gradient at cell centre
             QWEIGHT%XGRADWTC(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.2) THEN
             ! Weights for Y gradient at cell centre
             QWEIGHT%YGRADWTC(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.3) THEN
             ! Lambda0
             QWEIGHT%LAMBDA0(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.4) THEN
             ! Weights for value at western cell wall
             QWEIGHT%VALWTW(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.5) THEN
             ! Weights for value at eastern cell wall
             QWEIGHT%VALWTE(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.6) THEN
             ! Weights for value at southern cell wall
             QWEIGHT%VALWTS(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.7) THEN
             ! Weights for value at northern cell wall
             QWEIGHT%VALWTN(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.8) THEN
             ! Weights for X gradient at western cell wall
             QWEIGHT%XGRADWTW(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.9) THEN
             ! Weights for X gradient at eastern cell wall
             QWEIGHT%XGRADWTE(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.10) THEN
             ! Weights for Y gradient at southern cell wall
             QWEIGHT%YGRADWTS(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.11) THEN
             ! Weights for Y gradient at northern cell wall
             QWEIGHT%YGRADWTN(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.12) THEN
             ! X coordinates of neighbours
             QWEIGHT%XCOORDQUAD(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.13) THEN
             ! Y coordinates of neighbours
             QWEIGHT%YCOORDQUAD(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.14) THEN
             ! cell centre dXdX 2nd derivative weights
             QWEIGHT%XXDER2WTC(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.15) THEN
             ! cell centre dXdY 2nd derivative weights
             QWEIGHT%XYDER2WTC(ICIN,:) = VALS(:)
           ELSEIF (ITYP.EQ.16) THEN
             ! cell centre dYdY 2nd derivative weights
             QWEIGHT%YYDER2WTC(ICIN,:) = VALS(:)
           ELSE
             WRITE(NDSO,*) 'QA_READWT: WARNING: Unknown array type:',  &
                           ITYP
           END IF
        END DO
      END DO
      CLOSE(IUN)
!
! Check if the expected number of lines have been read into each array,
!  and set flags accordingly
      QWEIGHT%WTSREAD(1) = NREAD(1).EQ.NCASEQUAD .AND.                &
                           NREAD(2).EQ.NCASEQUAD
      QWEIGHT%WTSREAD(2) = NREAD(4).EQ.NCASEQUAD .AND.                &
                           NREAD(5).EQ.NCASEQUAD .AND.                &
                           NREAD(6).EQ.NCASEQUAD .AND.                &
                           NREAD(7).EQ.NCASEQUAD
      QWEIGHT%WTSREAD(3) = NREAD(8).EQ.NCASEQUAD .AND.                &
                           NREAD(9).EQ.NCASEQUAD .AND.                &
                           NREAD(10).EQ.NCASEQUAD .AND.               &
                           NREAD(11).EQ.NCASEQUAD
      QWEIGHT%WTSREAD(4) = NREAD(14).EQ.NCASEQUAD .AND.               &
                           NREAD(15).EQ.NCASEQUAD .AND.               &
                           NREAD(16).EQ.NCASEQUAD
      QWEIGHT%NCASEQUAD = NCASEQUAD
      QWEIGHT%MCASE = MCASE
      QWEIGHT%NNGBR = NNGBR
      RETURN
!
! Error traps
 901  CONTINUE
      WRITE(NDSE,*) 'QA_READWT: Unexpected END OF FILE encountered '
      RETURN
 902  CONTINUE
      WRITE(NDSE,*) 'QA_READWT: ERROR reading file '
      RETURN
!
      END SUBROUTINE QA_READWT
!/ ------------------------------------------------------------------- /
!/ ------------------------------------------------------------------- /
      SUBROUTINE QA_RECT ( NX0, NY0, LVLREF, GLOBAL, DEFAULT_TYPE,      &
                           UNDEF_TYPE, QTREE, IERR, NDSE )
!/
!         Richard Gorman, NIWA
!          April, 2008:      Origination
!          Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!      Set up a threaded quadtree structure (all at level 0) on a rectangular grid
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NX0,NY0      Int.     I  Size of rectangular grid at level 0
!       LVLREF       Int.     I  Reference level
!       GLOBAL       Log.     I  Flag for global grid (i.e. cyclic in X)
!       DEFAULT_TYPE Int.     I  Value of flag for created cells
!       UNDEF_TYPE   Int.     I  Value of flag for invalid cells
!       QTREE      QA_TREE   I/O Quadtree structure, with the following 
!                                components affected:
!          NCELL   Int.  Number of cells
!          NQUAD   Int.  Number of quads
!          LVLREF  Int.  Level of reference grid
!          LVLMAX  Int.  Maximum refinement level allowed
!          LVLHI   Int.  Maximum refinement level actually present
!          NX0     Int.  X dimension of level-zero base grid
!          NY0     Int.  Y dimension of level-zero base grid
!          INDQUAD I.A.  Index of the quad containing each input cell
!          INDSUB  I.A.  Index of the sub-quad containing each cell:
!                               0 = none (home cell is at level zero)
!                               1,2,3,4 = SW, SE, NW, NE corner of quad
!          INDLVL  I.A.  Level of each cell
!          INDML   I.A.  Multilevel index of cells
!          XYVAL   R.A.  Coordinates of cells (w.r.t. 
!                          reference rectangular grid)
!          QICELL  I.A.  Sea-point indices for the cells of the quad
!          QLEVEL  I.A.  Level of the cells within each quad
!          QPARENT I.A.  Quad index for the parent of each quad
!          QNBR    I.A.  Quad indices for the neighbours of each quad
!          QCHILD  I.A.  Quad indices for the children of each quad
!       IERR         Int.   O*  Return flag = 1 for error, else 0
!       NDSE         Int.   I*  Unit number for error output (if >0)
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_BMLSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC   Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NX0, NY0, LVLREF
      LOGICAL, INTENT(IN)     :: GLOBAL
      INTEGER, INTENT(IN)     :: DEFAULT_TYPE, UNDEF_TYPE
      TYPE(QA_TREE), INTENT(INOUT)    :: QTREE
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER           :: IX, IY, IQ
      REAL              :: DELX0, DELY0, X0, Y0
      INTEGER           :: MSEA, MQUAD
      INTEGER           :: NSEA, NQUAD
      INTEGER           :: IUN
      LOGICAL           :: MAKE_CELLS
      !
      NQUAD = NX0*NY0
      !
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0
      ! Array allocation checks:
      MQUAD = SIZE(QTREE%QICELL,1)
      MQUAD = MIN(MQUAD,SIZE(QTREE%QLEVEL,1))
      MQUAD = MIN(MQUAD,SIZE(QTREE%QPARENT,1))
      MQUAD = MIN(MQUAD,SIZE(QTREE%QCHILD,1))
      MQUAD = MIN(MQUAD,SIZE(QTREE%QNBR,1))
      !
      IF ( NQUAD.GT.MQUAD ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_RECT: INSUFFICIENT MEMORY ALLOCATED'
            WRITE(IUN,*) '                  FOR A QTREE QUAD ARRAY:'
            WRITE(IUN,*) ' ALLOCATED:', MQUAD
            WRITE(IUN,*) '  REQUIRED:', NQUAD
            write(IUN,*) 'Quadtree quad allocations:'
            write(IUN,*) 'QICELL:  ', SIZE(QTREE%QICELL,1)
            write(IUN,*) 'QLEVEL:  ', SIZE(QTREE%QLEVEL,1)
            write(IUN,*) 'QPARENT: ', SIZE(QTREE%QPARENT,1)
            write(IUN,*) 'QCHILD:  ', SIZE(QTREE%QCHILD,1)
            write(IUN,*) 'QNBR:    ', SIZE(QTREE%QNBR,1)
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
      !
      MSEA = SIZE(QTREE%INDQUAD,1)
      MSEA = MIN(MSEA,SIZE(QTREE%INDSUB,1))
      MSEA = MIN(MSEA,SIZE(QTREE%INDLVL,1))
      MSEA = MIN(MSEA,SIZE(QTREE%INDML,1))
      MSEA = MIN(MSEA,SIZE(QTREE%CELL_TYPE,1))
      MSEA = MIN(MSEA,SIZE(QTREE%XYVAL,1))
      !
      IF ( NQUAD.GT.MSEA ) THEN
         IF ( IUN.GT.0 ) THEN
            WRITE(IUN,*) 'ERROR IN QA_RECT: INSUFFICIENT MEMORY ALLOCATED'
            WRITE(IUN,*) '                  FOR A QTREE CELL ARRAY:'
            WRITE(IUN,*) ' ALLOCATED:', MSEA
            WRITE(IUN,*) '  REQUIRED:', NQUAD
            write(IUN,*) 'Quadtree cell allocations:'
            write(IUN,*) 'INDQUAD:   ', SIZE(QTREE%INDQUAD,1)
            write(IUN,*) 'INDSUB:    ', SIZE(QTREE%INDSUB,1)
            write(IUN,*) 'INDLVL:    ', SIZE(QTREE%INDLVL,1)
            write(IUN,*) 'INDML:     ', SIZE(QTREE%INDML,1)
            write(IUN,*) 'CELL_TYPE: ', SIZE(QTREE%CELL_TYPE,1)
            write(IUN,*) 'XYVAL:     ', SIZE(QTREE%XYVAL,1)
         END IF
         IF ( PRESENT(IERR) ) IERR = 1
         RETURN
      END IF
      !
      DELX0 = 2.**LVLREF
      DELY0 = 2.**LVLREF
      X0 = (1 + DELX0)/2.
      Y0 = (1 + DELY0)/2.
      !
!
!  Initialise
      QTREE%INDQUAD = 0
      QTREE%INDSUB = 0
      QTREE%QICELL = 0
      QTREE%QLEVEL = 0
      QTREE%QPARENT = 0
      QTREE%QNBR = 0
      QTREE%QCHILD = 0
      QTREE%XYVAL = 0.
      QTREE%INDLVL = 0
      QTREE%INDML = 0
      QTREE%NGBR = 0
      QTREE%CELL_TYPE = UNDEF_TYPE
!
      MAKE_CELLS = ( DEFAULT_TYPE .NE. UNDEF_TYPE )
!
      IQ = 0
      NSEA = 0
      DO IY=1,NY0
         !WRITE(*,*) 'IY: ',IY
         DO IX=1,NX0
            IQ = IQ + 1
            IF ( MAKE_CELLS ) NSEA = IQ
            !WRITE(*,*) IX, IY, NSEA
            QTREE%QICELL(IQ,0) = NSEA
            QTREE%QLEVEL(IQ) = 1
            IF (IX.EQ.1) THEN
               IF (GLOBAL) QTREE%QNBR(IQ,1) = IQ+NX0-1
            ELSE
               QTREE%QNBR(IQ,1) = IQ-1
            END IF
            IF (IX.EQ.NX0) THEN
               IF (GLOBAL) QTREE%QNBR(IQ,2) = IQ-NX0+1
            ELSE
               QTREE%QNBR(IQ,2) = IQ+1
            END IF
            IF (IY.GT.1) THEN
               QTREE%QNBR(IQ,3) = IQ-NX0
            END IF
            IF (IY.LT.NY0) THEN
               QTREE%QNBR(IQ,4) = IQ+NX0
            END IF
            !
            IF ( MAKE_CELLS ) THEN
               QTREE%INDQUAD(NSEA) = NSEA
               QTREE%XYVAL(NSEA,1) = X0 + DELX0*(IX-1)
               QTREE%XYVAL(NSEA,2) = Y0 + DELY0*(IY-1)
               QTREE%INDML(NSEA) = NSEA
               QTREE%NGBR(NSEA,1:4) = QTREE%QNBR(NSEA,1:4)
               QTREE%CELL_TYPE(NSEA) = DEFAULT_TYPE
            END IF
         END DO
      END DO
!
      QTREE%NX0 = NX0
      QTREE%NY0 = NY0
      QTREE%LVLREF = LVLREF
      IF ( QTREE%LVLMAX .LT. LVLREF ) QTREE%LVLMAX = LVLREF
      QTREE%LVLHI = LVLREF
      QTREE%NCELL = NSEA
      QTREE%NCELL_DEF = NSEA
      QTREE%NQUAD = NQUAD
      QTREE%UNDEF_TYPE = UNDEF_TYPE
!
      RETURN
      END SUBROUTINE QA_RECT
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_REFINE ( IMODE, QTREE, NREF, IQREF, ISREF,        &
                             ISEAOLD, ISEANEW, IERR, NDSE, MAPML,     &
                             QSPARE )
!/
!         Richard Gorman, NIWA
!          April, 2008:      Origination
!          Jan-April 2014:   Rewrite using QTREE, QSPARE structure
!
!  1. Purpose :
!
!      Refine selected cells of a threaded quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMODE   Int.   I   Input option: =0 all cells
!                                        =1 to use IQREF, ISREF,
!                                        =2 to use ISEAOLD
!       QTREE   QA_TREE I/O  Quadtree structure, with the following
!                             components affected by this subroutine:
!          NCELL   Int.  I/O  Number of cells
!          NQUAD   Int.  I/O  Number of quads
!          LVLREF  Int.   I   Level of reference grid
!          LVLMAX  Int.   I   Maximum refinement level allowed
!          LVLHI   Int.  I/O  Maximum refinement level actually present
!          NX0     Int.   I   X dimension of level-zero base grid
!          NY0     Int.   I   Y dimension of level-zero base grid
!          UNDEF_TYPE Int. value of flag for unused cell indices
!          KEEP_REF Log.  I   retain cell indices for refined cells
!          INDQUAD I.A.  I/O  Index of the quad containing each input cell
!          INDSUB  I.A.  I/O  Index of the sub-quad containing each cell:
!                               0 = none (home cell is at level zero)
!                               1,2,3,4 = SW, SE, NW, NE corner of quad
!          INDLVL  I.A.  I/O  Level of each cell
!          CELL_TYPE I.A. Flag for wet/dry/boundary/etc. cell
!          INDML   I.A.  I/O  Multilevel index of cells
!          XYVAL   R.A.  I/O  Coordinates of cells (w.r.t. 
!                             reference rectangular grid)
!          QICELL  I.A.  I/O  Sea-point indices for the cells of the quad
!          QLEVEL  I.A.  I/O  Level of the cells within each quad
!          QPARENT I.A.  I/O  Quad index for the parent of each quad
!          QNBR    I.A.  I/O  Quad indices for the neighbours of each quad
!          QCHILD  I.A.  I/O  Quad indices for the children of each quad
!       NREF    Int.  I/O  Number of cells to refine
!       IQREF   I.A.  I/O  Quad numbers of cells to refine
!       ISREF   I.A.  I/O  Sub-Quad indices of cells to refine
!       ISEAOLD I.A.  I/O  Old sea indices of cells that have been refined
!       NSEA2   Int.   O   Final number of cells
!       NQUAD2  Int.   O   Final number of quads
!       ISEANEW I.A.   O   New sea indices of cells that have been refined
!       IERR     Int.   O*  Return flag = 1 for error, else 0
!       NDSE     Int.   I*  Unit number for error output (if >0)
!       MAPML   I.A.   I*  Array of cell flags on multilevel grid 
!       QSPARE  QA_SPARE I/O* Arrays for treating spare cells in quadtree structure,
!                             with the following components used:
!          NPART   Int.   I   Number of partitions
!          PARTNO  I.A.  I/O  Partition number for each cell
!          JSVAL   I.A.  I/O  Local cell index for each global index
!          NCLOC   I.A.  I/O  Actual number of cells in each partition
!          ISVAL   I.A.  I/O  Global cell index for each local index & partition
!          NSPARE  I.A.  I/O  Number of spare local cell indices
!          JSSPARE I.A.  I/O  List of spare local cell indices
!     ----------------------------------------------------------------
!                        * = optional
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     ML_CHILD   Func. qa_utils Multilevel index of the SW child of a given cell
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_ADGR     Subr. QA_UTILS  General grid adaptation routine 
!     QA_BMLSTRUC Subr. QA_UTILS  Create a threaded quadtree structure 
!                                  from a multilevel bathymetry grid
!     QA_BRSTRUC  Subr. QA_UTILS  Create a threaded quadtree structure
!                                  from a rectangular grid
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMODE
      TYPE(QA_TREE), INTENT(INOUT)     :: QTREE
      INTEGER, INTENT(INOUT)  :: NREF
      INTEGER, INTENT(INOUT)  :: IQREF(:), ISREF(:), ISEAOLD(:)
      INTEGER, INTENT(OUT)    :: ISEANEW(:,:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
      INTEGER, OPTIONAL, INTENT(IN)     :: MAPML(:)
      TYPE(QA_SPARE), OPTIONAL, INTENT(INOUT)  :: QSPARE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER           :: MSEA, MQUAD
      INTEGER           :: NPART
      INTEGER           :: NSEA1, NQUAD1
      INTEGER           :: NSEA2, NQUAD2
      REAL              :: DELX0, DELY0
      INTEGER           :: IDIR, QNC, IREF, IQ, ISUB, QISEA0, IQMOD,  &
                           II, NSTEPS, IML, ISC, IMLP, ISN, JSN,      &
                           PV, NSP, NINDML, IQPAR, MAPNEW, NSTEPS2
      INTEGER           :: QN(4), QS(4)
      INTEGER           :: IDREV(4)
      INTEGER           :: IUN
      INTEGER           :: IRFAC
      INTEGER, ALLOCATABLE  :: IQM(:)
      LOGICAL           :: DOSPARE, DOPART
      REAL              :: DELX, DELY, X0, Y0
      DATA IDREV/ 2, 1, 4, 3 /
!/
!/ ------------------------------------------------------------------- /
!/
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE
      IF ( PRESENT(IERR) ) IERR = 0
      NSEA1 = QTREE%NCELL
      NQUAD1 = QTREE%NQUAD
!
      MQUAD = SIZE(QTREE%QICELL,1)
      MSEA = SIZE(QTREE%INDQUAD,1)
      ALLOCATE ( IQM(MSEA) )
      DELX0 = 2.**QTREE%LVLREF
      DELY0 = DELX0
      !WRITE(*,*) 'QA_REFINE:  IMODE, MSEA, MQUAD: ', IMODE, MSEA, MQUAD
      !WRITE(*,*) 'QA_REFINE:  NSEA1, NQUAD1, NREF: ', NSEA1, NQUAD1, NREF
      !WRITE(*,*) 'QA_REFINE:  DELX0, DELY0: ', DELX0, DELY0
      IF ( PRESENT(MAPML) ) THEN
         NINDML = SIZE(MAPML)
      ELSE
         NINDML = 0
      END IF
!
      DOSPARE = PRESENT(QSPARE)
      IF ( DOSPARE ) THEN 
         NPART = QSPARE%NPART
      ELSE
         NPART = 1
      END IF
      DOPART = NPART.GT.1
!
      ISEANEW = 0               
!
      NSEA2 = NSEA1
      NQUAD2 = NQUAD1
      IQ = 1
      ISUB = 0
      IF (IMODE.EQ.0) THEN
         NSTEPS = 4*NQUAD1
      ELSE
         NSTEPS = NREF
      END IF
      ! Cells created per refinement:
      IRFAC = 3
      IF (QTREE%KEEP_REF) IRFAC = 4
      IF (IRFAC*NSTEPS.GT.MSEA-NSEA1) THEN
         NSTEPS = (MSEA-NSEA1)/IRFAC
         NREF = NSTEPS
         IF ( IUN.GT. 0 ) THEN
            WRITE(IUN,*) 'QA_REFINE: '
            WRITE(IUN,*) 'REDUCING NUMBER OF REFINEMENTS TO '
            WRITE(IUN,*) 'FIT CELL ARRAY SIZE. NREF = ', NREF
         END IF
         IF ( PRESENT(IERR) ) IERR = -1
      END IF
! First loop, to create IQREF, ISREF arrays, if IMODE .NE. 1
      IF ( IMODE.NE.1 ) THEN
        NSTEPS2 = NSTEPS
        DO IREF=1, NSTEPS
!
         !WRITE(*,*) 'IREF: ',IREF
          IF (IMODE.EQ.0) THEN
          ! Refine all cells (that don't have children)
            ISUB = ISUB+1
            IF (ISUB.EQ.5) THEN
              ISUB = 1
              IQ = IQ + 1
            END IF
            IF (IQ.GT.NQUAD1) THEN
              NSTEPS2 = IREF - 1
              NREF = NSTEPS2
              EXIT
            END IF
            IQREF(IREF) = IQ
            ISREF(IREF) = ISUB
          ELSEIF (IMODE.GT.1 ) THEN
            IQ = QTREE%INDQUAD(ISEAOLD(IREF))
            ISUB = QTREE%INDSUB(ISEAOLD(IREF))
            IQREF(IREF) = IQ
            ISREF(IREF) = ISUB
          END IF
          IF (IQ.GT.MQUAD) THEN
            NSTEPS2 = IREF - 1
            NREF = NSTEPS2
            IF ( IUN.GT. 0 ) THEN
               WRITE(IUN,*) 'QA_REFINE: '
               WRITE(IUN,*) 'REDUCING NUMBER OF REFINEMENTS TO '
               WRITE(IUN,*) 'FIT QUAD ARRAY SIZE. NREF = ', NREF
               WRITE(IUN,*) 'ALLOCATED: ', MQUAD
            END IF
            IF ( PRESENT(IERR) ) IERR = -1
            EXIT
          END IF
        END DO
      END IF
      NSTEPS = NSTEPS2
      DO IREF=1, NSTEPS
         IQ = IQREF(IREF)
         ISUB = ISREF(IREF)
         !write(*,*) 'QA_REFINE: IREF,IQ,SUB: ',IREF,IQ, ISUB
         !
         ! Don't refine if this quad is at maximum level:
         !
         !IF ( QTREE%QLEVEL(IQ).GE.QTREE%LVLMAX ) CYCLE
         !
         IF ( IQ.EQ.0 ) CYCLE
         QISEA0 = QTREE%QICELL(IQ,0)
         IQPAR = QTREE%QPARENT(IQ)
         !IF (IQPAR.EQ.0 .AND. QISEA0.NE.0) THEN
         IF (IQPAR.EQ.0 .AND. ISUB.EQ.0) THEN
!  This is a top-level "empty" quad: don't add a new quad,
!  but populate this one
!  Don't refine if this is already a split cell:
            IF (QTREE%QICELL(IQ,ISUB).EQ.0) THEN
               IQM(IREF) = 0
               CYCLE
            END IF
            IQMOD = IQ
            ISUB = 4
            !QTREE%QLEVEL(IQMOD) = QTREE%QLEVEL(IQ)+1
            IQM(IREF) = IQMOD
            ISREF(IREF) = 0
         ELSE
!  Don't try to add a new quad if there is already a child:
            IF ( QTREE%QCHILD(IQ,ISUB).NE.0) THEN
               IQM(IREF) = 0
!  Don't refine if this is already a split cell:
               IF (QTREE%QICELL(IQ,ISUB).EQ.0) THEN
                  CYCLE
               END IF
               IQMOD = QTREE%QCHILD(IQ,ISUB)
            ELSE
               IF ( NQUAD2.GT.MQUAD ) THEN
!
!  At the limit of allocated quads: exit the loop
                  IF ( IUN.GT.0 ) THEN
                     WRITE(IUN,*) 'QA_REFINE: '
                     WRITE(IUN,*) 'MAX NUMBER OF QUADS EXCEEDED'
                     WRITE(IUN,*) 'ALLOCATED: ', MQUAD
                  END IF
                  IF ( PRESENT(IERR) ) IERR = 1
                  NREF = IREF-1
                  NSTEPS = NREF
                  EXIT
               END IF
               NQUAD2 = NQUAD2 + 1
               IQMOD = NQUAD2
               IQM(IREF) = IQMOD
!  Add a child in the parent quad:
               QTREE%QCHILD(IQ,ISUB) = IQMOD
               DO II=1,4
                  QTREE%QCHILD(IQMOD,II) = 0
               END DO
               QTREE%QLEVEL(IQMOD) = QTREE%QLEVEL(IQ)+1
               QTREE%QPARENT(IQMOD) = IQ
            END IF
!  Save the cell index to pass down, and set to zero here
            QISEA0 = QTREE%QICELL(IQ,ISUB)
            QTREE%QICELL(IQ,ISUB) = 0
!
!  Fill in properties of the new quad:
            !WRITE(*,*) 'IQMOD: ',IQMOD
            QTREE%QICELL(IQMOD,0) = QTREE%QICELL(IQ,0)
         END IF
!
!        Properties of the old cell:
         IF ( QISEA0.GT.MSEA ) THEN
            IF ( IUN.GT.0 ) THEN
               WRITE(IUN,*) 'QA_REFINE:'
               WRITE(IUN,*) 'SEA POINT INDEX OUT OF RANGE'
               WRITE(IUN,*) 'ALLOCATED: ', MSEA
               WRITE(IUN,*) 'QISEA0 = ',QISEA0
               WRITE(IUN,*) 'IQ, ISUB = ', IQ,ISUB
               WRITE(IUN,*) 'QTREE%QICELL(IQ,0) = ', QTREE%QICELL(IQ,0)
               WRITE(IUN,*) 'QCHILD(IQ,ISUB) = ', QTREE%QCHILD(IQ,ISUB)
               WRITE(IUN,*) 'IQMOD = ', IQMOD
               WRITE(IUN,*) 'IREF = ', IREF
            END IF
            IF ( PRESENT(IERR) ) IERR = 1
            RETURN
         ELSE IF ( QISEA0.EQ.0 ) THEN
            ! Empty level-0 quad
            X0 = MODULO(IQMOD-1,QTREE%NX0) + 1
            Y0 = (IQMOD - X0)/QTREE%NX0 + 1
            IMLP = IQMOD
            MAPNEW = QTREE%UNDEF_TYPE
            IML = ML_CHILD(IMLP, QTREE%NX0, QTREE%NY0)
         ELSE
            X0 = QTREE%XYVAL(QISEA0,1)
            Y0 = QTREE%XYVAL(QISEA0,2)
            IMLP = QTREE%INDML(QISEA0)
            MAPNEW = QTREE%CELL_TYPE(QISEA0)
            IML = ML_CHILD(IMLP, QTREE%NX0, QTREE%NY0)
         END IF
!        Assign cell indices to the cells of the (new) quad
         IF ( QTREE%KEEP_REF .OR. QISEA0.EQ.0 ) THEN
            QTREE%QICELL(IQMOD,0) = QISEA0
            ISC = 0
            DELX = DELX0*2.**(-QTREE%QLEVEL(IQMOD)-1)
            DELY = -DELY0*2.**(-QTREE%QLEVEL(IQMOD)-1)
            IML = IML - 1
         ELSE
            QTREE%QICELL(IQMOD,0) = 0
!        First cell takes over index of the parent : SW cell if wet
            ISC = 1
            DELX = -DELX0*2.**(-QTREE%QLEVEL(IQMOD)-1)
            DELY = -DELY0*2.**(-QTREE%QLEVEL(IQMOD)-1)
!
!        If a MAPML array identifying wet/dry cells is present, use it
!        to only create a seapoint for this child if it is wet. If not, 
!        look for the first wet child:
!
            IF ( NINDML.GT.1 .AND. PRESENT(MAPML) ) THEN
               DO II=1,4
                  IF (IML.LE.NINDML) THEN
                     IF ( MAPML(IML).EQ.QTREE%UNDEF_TYPE ) THEN
                        DELX = -DELX
                        IF (II.EQ.2) THEN
                           DELY=-DELY
                           IML = IML + QTREE%NX0*2**QTREE%QLEVEL(IQMOD) - 1
                        ELSE
                           IML = IML + 1
                        END IF
                        CYCLE
                     END IF
                  END IF
                  MAPNEW = MAPML(IML)
                  ISC = II
                  EXIT
               END DO
            END IF
            QTREE%QICELL(IQMOD,ISC) = QISEA0
            QTREE%INDQUAD(QISEA0) = IQMOD
            QTREE%INDSUB(QISEA0) = ISC
            QTREE%INDLVL(QISEA0) = QTREE%QLEVEL(IQMOD)
            QTREE%XYVAL(QISEA0,1) = X0 + DELX
            QTREE%XYVAL(QISEA0,2) = Y0 + DELY
            QTREE%INDML(QISEA0) = IML
            QTREE%CELL_TYPE(QISEA0) = MAPNEW
            IF (IMODE.GE.1) ISEANEW(IREF,ISC) = QISEA0
         END IF
         IF (IMODE.EQ.1) ISEAOLD(IREF) = QISEA0
!  Other cells take new indices
         PV = 1
         IF ( DOPART ) PV = QSPARE%PARTNO(QISEA0)
         NSP = 0
         DO II=ISC+1,4
            DELX = -DELX
            IF (II.EQ.3) THEN
              DELY=-DELY
              IML = IML + QTREE%NX0*2**QTREE%QLEVEL(IQMOD) - 1
            ELSE
              IML = IML + 1
            END IF
!
!           If a MAPML array identifying wet/dry cells is present, use it
!           to only create a seapoint for this child if it is wet:
!
            IF ( NINDML.GT.0 .AND. PRESENT(MAPML) ) THEN
               IF (IML.LE.NINDML) THEN
                  IF ( MAPML(IML).EQ.QTREE%UNDEF_TYPE ) CYCLE
                  MAPNEW = MAPML(IML)
               END IF
            END IF
            IF ( DOSPARE ) NSP = QSPARE%NSPARE(PV)
            IF ( DOSPARE .AND. NSP.GT.0) THEN
               ! spare cell indices available: use one
               JSN = QSPARE%JSSPARE(NSP,PV)
               IF ( DOPART ) THEN
                  ISN = QSPARE%ISVAL(JSN,PV)
               ELSE
                  ISN = JSN
               END IF
               QSPARE%JSSPARE(NSP,PV) = 0
               QSPARE%NSPARE(PV) = NSP - 1
            ELSE
               ! no spares: append to the list
               NSEA2 = NSEA2+1
               ISN = NSEA2
               IF ( DOPART ) THEN
                  QSPARE%PARTNO(ISN) = PV
                  QSPARE%NCLOC(PV) = QSPARE%NCLOC(PV)+1
                  JSN = QSPARE%NCLOC(PV)
                  QSPARE%ISVAL(JSN,PV) = ISN
                  QSPARE%JSVAL(ISN) = JSN
               END IF
            END IF
            QTREE%NCELL_DEF = QTREE%NCELL_DEF + 1
            QTREE%QICELL(IQMOD,II) = ISN
            !write(*,*) 'QA_REFINE QICELL(IQMOD,II) = ISN',IQMOD,II,ISN
            IF (ISN.GT.MSEA) THEN
               IF ( IUN.GT.0 ) THEN
                  WRITE(IUN,*) 'QA_REFINE'
                  WRITE(IUN,*) 'EXCEEDING MAX. NO. OF CELLS'
                  WRITE(IUN,*) 'ALLOCATED: ', MSEA
                  WRITE(IUN,*) 'INITIAL NSEA: ', NSEA1
                  WRITE(IUN,*) 'INITIAL NQUAD: ', NQUAD1
                  WRITE(IUN,*) 'ATTEMPTING TO CREATE: ', ISN
                  WRITE(IUN,*) 'IREF: ', IREF
                  WRITE(IUN,*) 'NSTEPS: ', NSTEPS
               END IF
               IF ( PRESENT(IERR) ) IERR = 1
               RETURN
            ELSE
               QTREE%INDQUAD(ISN) = IQMOD
               QTREE%INDSUB(ISN) = II
               QTREE%INDLVL(ISN) = QTREE%QLEVEL(IQMOD)
               QTREE%LVLHI = MAX(QTREE%LVLHI,QTREE%INDLVL(ISN))
               QTREE%INDML(ISN) = IML
               QTREE%CELL_TYPE(ISN) = MAPNEW
               QTREE%XYVAL(ISN,1) = X0 + DELX
               QTREE%XYVAL(ISN,2) = Y0 + DELY
            END IF
            IF (IMODE.GE.1) ISEANEW(IREF,II) = ISN
         END DO
      END DO
      !WRITE(*,*) 'NSTEPS: ',NSTEPS
!
!  Second loop to compute neighbours
      DO IREF=1, NSTEPS
!
         !WRITE(*,*) 'IREF: ',IREF
         IQ = IQREF(IREF)
         ISUB = ISREF(IREF)
         IQMOD = IQM(IREF)
         !WRITE(*,*) 'IQ, ISUB, IQMOD: ',IQ, ISUB, IQMOD
         IF (IQ.GT.0 .AND. ISUB.GT.0 .AND. IQMOD.GT.0) THEN
!  Find the quad number and subquad index for the four equal-level
!  neighbours (W(1), E(2), S(3), N(4)) of the home cell, if they exist
!            IF (ISUB.EQ.0) THEN
!!     0. Cell is a leaf at level zero
!               QN(1) = QTREE%QNBR(IQ,1)
!               QS(1) = 0
!               QN(2) = QTREE%QNBR(IQ,2)
!               QS(2) = 0
!               QN(3) = QTREE%QNBR(IQ,3)
!               QS(3) = 0
!               QN(4) = QTREE%QNBR(IQ,4)
!               QS(4) = 0
!            ELSE IF (ISUB.EQ.1) THEN
            IF (ISUB.EQ.1) THEN
!  Cell is in the SW quadrant of the quad:
!  ... W leads to the SE quadrant of the W neighbour
               QN(1) = QTREE%QNBR(IQ,1)
               QS(1) = 2
!  ... E leads to the SE quadrant of the same quad
               QN(2) = IQ
               QS(2) = 2
!  ... S leads to the NW quadrant of the S neighbour
               QN(3) = QTREE%QNBR(IQ,3)
               QS(3) = 3
!  ... N leads to the NW quadrant of the same quad
               QN(4) = IQ
               QS(4) = 3
            ELSE IF(ISUB.EQ.2) THEN
!  Cell is in the SE quadrant of the quad:
!  ... W leads to the SW quadrant of the same quad
               QN(1) = IQ
               QS(1) = 1
!  ... E leads to the SW quadrant of the E neighbour
               QN(2) = QTREE%QNBR(IQ,2)
               QS(2) = 1
!  ... S leads to the NE quadrant of the S neighbour
               QN(3) = QTREE%QNBR(IQ,3)
               QS(3) = 4
!  ... N leads to the NE quadrant of the same quad
               QN(4) = IQ
               QS(4) = 4
            ELSE IF(ISUB.EQ.3) THEN
!  Cell is in the NW quadrant of the quad:
!  ... W leads to the NE quadrant of the W neighbour
               QN(1) = QTREE%QNBR(IQ,1)
               QS(1) = 4
!  ... E leads to the NE quadrant of the same quad
               QN(2) = IQ
               QS(2) = 4
!  ... S leads to the SW quadrant of the same quad
               QN(3) = IQ
               QS(3) = 1
!  ... N leads to the SW quadrant of the N neighbour
               QN(4) = QTREE%QNBR(IQ,4)
               QS(4) = 1
            ELSE IF(ISUB.EQ.4) THEN
!  Cell is in the NE quadrant of the quad:
!  ... W leads to the NW quadrant of the same quad
               QN(1) = IQ
               QS(1) = 3
!  ... E leads to the NW quadrant of the E neighbour
               QN(2) = QTREE%QNBR(IQ,2)
               QS(2) = 3
!  ... S leads to the SE quadrant of the same quad
               QN(3) = IQ
               QS(3) = 2
!  ... N leads to the SE quadrant of the N neighbour
               QN(4) = QTREE%QNBR(IQ,4)
               QS(4) = 2
            END IF
            !WRITE(*,*) 'QN: ', QN
            !WRITE(*,*) 'QS: ', QS
!
!  Loop over W, E, S, N directions, targetting the quad in that direction:
            DO IDIR=1,4
               IF (QN(IDIR) .EQ. 0) THEN
!  No target quad (e.g. domain limit), so no neighbour in this direction
                  QTREE%QNBR(IQMOD,IDIR) = 0
!               ELSEIF (QTREE%QICELL(QN(IDIR),0) .GT. 0) THEN
!  Target quad is a leaf cell, so is a neighbour at lower order
!                  QTREE%QNBR(IQMOD,IDIR) = QN(IDIR)
!                    inverse neighbour:                     
                     !QTREE%QNBR(QN(IDIR),IDREV(IDIR)) = IQMOD
               ELSE
!  See if neighbouring quad has a child at the target subquad index
                  !WRITE(*,*) 'IDIR,QN, QS: ',IDIR,QN(IDIR), QS(IDIR)
                  QNC = QTREE%QCHILD(QN(IDIR),QS(IDIR))
                  IF (QNC.EQ.0) THEN
!                No: neighbour is the parent quad
                     !QTREE%QNBR(IQMOD,IDIR) = QN(IDIR)
!                    inverse neighbour:                     
                     !QTREE%QNBR(QN(IDIR),IDREV(IDIR)) = IQMOD
!                No: no equal-level neighbour 
                     QTREE%QNBR(IQMOD,IDIR) = 0
                  ELSE
!                Yes: neighbour is the child quad
                     QTREE%QNBR(IQMOD,IDIR) = QNC
!                    inverse neighbour:                     
                     QTREE%QNBR(QNC,IDREV(IDIR)) = IQMOD
                  END IF
               END IF
            END DO
         END IF
      END DO
!
      QTREE%NCELL = NSEA2
      QTREE%NQUAD = NQUAD2

      RETURN
      END SUBROUTINE QA_REFINE
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_REFVAL( A, NREF, NNGBR, ISO, ISN, NGBRV, DXWT,    &
                            DYWT )
!/
!/       Richard Gorman, NIWA
!/         16-Aug-2012 : Origination.
!          Jan-April 2014:   Rename
!
!  1. Purpose :
!
!     Assign values of a field in refined cells using estimated gradients 
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NREF    Int.   I   Number of cells to refine
!       NNGBR   Int.   I   Max. Number of neighbours per cell
!       ISO     I.A.   I   List of (old) cells to be refined 
!       ISN     I.A.   I   List of (new) cells resulting from refinement 
!       NGBRV   I.A.   I   Array of neighbours of old cells
!       DXWT    R.A.   I   Array of weights for computing x-gradients at old cells
!       DYWT    R.A.   I   Array of weights for computing y-gradients at old cells
!       A       R.A.  I/O  Array that is refined
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       None.
!
!  5. Called by :
!
!       WQADGR   Quadtree grid adaptation routine.
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(INOUT)  ::  A(:)
      INTEGER, INTENT(IN)  ::  NREF, NNGBR
      INTEGER, INTENT(IN)  ::  ISO(NREF), ISN(NREF,4), NGBRV(NREF,NNGBR)
      REAL, INTENT(IN)     ::  DXWT(NREF,0:NNGBR), DYWT(NREF,0:NNGBR)
!
!  Local parameters
!
      REAL       :: V0(NREF), XGRAD(NREF), YGRAD(NREF)
      INTEGER    :: K, II, ISEA1, IREF, INB, MSEA
!
      MSEA = SIZE(A,1)      
!
!  First loop to compute central values and gradients of A at each of the (old)
!  cells
!
      DO IREF=1,NREF
         ISEA1 = ISO(IREF)
         IF ( ISEA1.LE. 0 .OR. ISEA1.GT.MSEA ) CYCLE
         V0(IREF) = A(ISEA1)
         XGRAD(IREF) = DXWT(IREF,0)*A(ISEA1)
         YGRAD(IREF) = DYWT(IREF,0)*A(ISEA1)
         DO K=1,NNGBR
            INB = NGBRV(IREF,K)
            IF (INB.GT.0 .AND. INB.LE.MSEA) THEN
               XGRAD(IREF) = XGRAD(IREF) + DXWT(IREF,K)*A(INB)
               YGRAD(IREF) = YGRAD(IREF) + DYWT(IREF,K)*A(INB)
            END IF
         END DO
      END DO
!
!  Second loop to use central values and gradients to estimate A at the 4 child
!  cells of each refined cell
!
      DO IREF=1,NREF
         ISEA1 = ISO(IREF)
         IF ( ISEA1.LE.0 .OR. ISEA1.GT.MSEA ) CYCLE
         II = ISN(IREF,1)
         IF (II.GT.0) A(II) = V0(IREF) - .25*XGRAD(IREF) - .25*YGRAD(IREF)
         II = ISN(IREF,2)
         IF (II.GT.0) A(II) = V0(IREF) + .25*XGRAD(IREF) - .25*YGRAD(IREF)
         II = ISN(IREF,3)
         IF (II.GT.0) A(II) = V0(IREF) - .25*XGRAD(IREF) + .25*YGRAD(IREF)
         II = ISN(IREF,4)
         IF (II.GT.0) A(II) = V0(IREF) + .25*XGRAD(IREF) + .25*YGRAD(IREF)
      END DO
      RETURN
      END SUBROUTINE QA_REFVAL
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_REMUNDEF( QTREE, NEWCELL, IERR, NDSE )
!/
!       Richard Gorman, NIWA
!         May 2014:   Origination.
!
!  1. Purpose :
!
!      Remove undefined points from a threaded quadtree structure 
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ---------------------------------------------------------------
!       QTREE   QA_TREE I/O  Quadtree structure
!       NEWCELL  I.A.    O*  Mapping from old to new cell indices
!       IERR     Int.    O*  Return flag = 1 for error, else 0
!       NDSE     Int.    I*  Unit number for error output (if >0)
!     ---------------------------------------------------------------
!                         * optional
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  5. Called by :
!
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ Parameter list
      TYPE(QA_TREE), INTENT(INOUT)   :: QTREE
      INTEGER, OPTIONAL, INTENT(OUT) :: NEWCELL(:)
      INTEGER, OPTIONAL, INTENT(OUT)  :: IERR
      INTEGER, OPTIONAL, INTENT(IN)   :: NDSE
!
! Local variables.
      INTEGER              :: NCELL1, NCELL2, IUN, IC, IC2, NR1, NR2, &
                              NC2, LVLHI, II, IQ, IQ2, NCH, NQUAD1,   &
                              NQUAD2, NB
      INTEGER, ALLOCATABLE :: NEW(:), NEWQ(:)
!
!
      IUN = 6
      IF ( PRESENT(NDSE) ) IUN = NDSE      
      IF ( PRESENT(IERR) ) IERR = 0
      NR1 = SIZE(QTREE%NCASE,1)
      NR2 = SIZE(QTREE%GNBR,1)
      NC2 = SIZE(QTREE%GNBR,2)
      NCELL1 = QTREE%NCELL
      NQUAD1 = QTREE%NQUAD
      !
      IF ( PRESENT(NEWCELL) ) THEN
         IF ( SIZE(NEWCELL,1).LT.NCELL1 ) THEN
            IF ( IUN.GT.0 ) THEN
               WRITE(IUN,*) 'QA_REMUNDEF: NEWCELL ALLOCATION ERROR'
               WRITE(IUN,*) 'ALLOCATED: ', SIZE(NEWCELL,1)
               WRITE(IUN,*) '   NEEDED: ', NCELL1
            END IF
            IF ( PRESENT(IERR) ) IERR = -1
            RETURN
         END IF
      END IF
      ALLOCATE ( NEW(NCELL1) )
      NEW = 0
!
! Loop through cells, to identify those which are undefined, and
! create a renumbering map.
!
      IC2 = 0
      DO IC=1,NCELL1
         IF ( QTREE%CELL_TYPE(IC).EQ.QTREE%UNDEF_TYPE ) CYCLE
         IC2 = IC2 + 1
         NEW(IC) = IC2
      END DO
      NCELL2 = IC2
!
      IF ( PRESENT(NEWCELL) ) NEWCELL(1:NCELL1) = NEW(1:NCELL1)
!
!  No undefined cells: return with no further work to do
!
      IF (NCELL2.EQ.NCELL1 ) RETURN
!
!  Loop through cells, renumbering cell indices:
!
      LVLHI = 0
      DO IC=1,NCELL1
         IC2 = NEW(IC)
         ! For an old cell being removed, move on:
         IF ( IC2.EQ.0 ) CYCLE
         LVLHI = MAX(LVLHI,QTREE%INDLVL(IC))
         ! even if IC=IC2, the neighbour indices may have changed:
         DO II=1,8
            NB = QTREE%NGBR(IC,II)
            IF ( NB.GT.0 ) THEN
               QTREE%NGBR(IC2,II) = NEW(NB)
            ELSE
               QTREE%NGBR(IC2,II) = 0
            END IF
         END DO
         IF ( QTREE%IWTORDER.GE.2 .AND. IC.GE.NR2 ) THEN
            DO II=1,NC2
               NB = QTREE%GNBR(IC,II)
               IF ( NB.GT.0 ) THEN
                  QTREE%GNBR(IC2,II) = NEW(NB)
               ELSE
                  QTREE%GNBR(IC2,II) = 0
               END IF
            END DO
         END IF
         ! "single cell" data doesn't need to be moved for unchanged indices:
         IF ( IC2.EQ.IC ) CYCLE
         QTREE%XYVAL(IC2,:) = QTREE%XYVAL(IC,:)
         QTREE%INDQUAD(IC2) = QTREE%INDQUAD(IC)
         QTREE%INDSUB(IC2) = QTREE%INDSUB(IC)
         QTREE%INDML(IC2) = QTREE%INDML(IC)
         QTREE%INDLVL(IC2) = QTREE%INDLVL(IC)
         QTREE%CELL_TYPE(IC2) = QTREE%CELL_TYPE(IC)
         IF ( QTREE%IWTORDER.GE.1 .AND. IC.GE.NR1 )                   &
              QTREE%NCASE(IC2) = QTREE%NCASE(IC)
         IF ( QTREE%IWTORDER.GE.2 .AND. IC.GE.NR2 ) THEN
            QTREE%INDWT(IC2) = QTREE%INDWT(IC)
            QTREE%NGNBR(IC2) = QTREE%NGNBR(IC)
         END IF
         IF ( IC.GT.NCELL2 ) THEN
            ! Set null values for cell indices that won't be overwritten:
            QTREE%XYVAL(IC,:) = 0.
            QTREE%INDQUAD(IC) = 0
            QTREE%INDSUB(IC) = 0
            QTREE%INDML(IC) = 0
            QTREE%INDLVL(IC) = 0
            QTREE%CELL_TYPE(IC) = QTREE%UNDEF_TYPE
            QTREE%NGBR(IC,:) = 0
            IF ( QTREE%IWTORDER.GE.1 .AND. IC.GE.NR1 )                &
                 QTREE%NCASE(IC) = 0
            IF ( QTREE%IWTORDER.GE.2 .AND. IC.GE.NR2 ) THEN
               QTREE%INDWT(IC) = 0
               QTREE%NGNBR(IC) = 0
               QTREE%GNBR(IC,:) = 0
            END IF
         END IF
      END DO
      !
      QTREE%NCELL = NCELL2
      QTREE%NCELL_DEF = NCELL2
      QTREE%LVLHI = LVLHI
      !
      ! Loop through quads, updating cell indices in QICELL,
      ! and also locating quads (level>0) with neither child quads 
      ! nor leaf cells for removal.
      !
      ALLOCATE ( NEWQ(NQUAD1) )
      NEWQ = 0
      !
      IQ2 = 0
      DO IQ=1,NQUAD1
         NCH = 0
         DO II=0,4
            IC = QTREE%QICELL(IQ,II)
            IF ( IC.GT.0 ) THEN
               QTREE%QICELL(IQ,II) = NEW(IC)
               NCH = NCH + 1
            END IF
         END DO
         IF ( NCH.EQ.0 ) THEN
            DO II=1,4
               IC = QTREE%QCHILD(IQ,II)
               IF ( IC.GT.0 ) NCH = NCH + 1
            END DO
         END IF 
         IF ( NCH.EQ.0 .AND. QTREE%QLEVEL(IQ).GT.1 ) CYCLE
         IQ2 = IQ2 + 1
         NEWQ(IQ) = IQ2
      END DO
      NQUAD2 = IQ2
!
!  Loop through quads, renumbering quad indices
!
      DO IQ=1,NQUAD1
         IQ2 = NEWQ(IQ)
         IF ( IQ2.EQ.0 ) CYCLE
         QTREE%QICELL(IQ2,:) = QTREE%QICELL(IQ,:)
         QTREE%QLEVEL(IQ2) = QTREE%QLEVEL(IQ)
         IC = QTREE%QPARENT(IQ)
         IF ( IC.GT.0 ) THEN
            QTREE%QPARENT(IQ2) = NEWQ(IC)
         ELSE
            QTREE%QPARENT(IQ2) = 0
         END IF
         DO II=1,4
           IC = QTREE%QCHILD(IQ,II)
           IF ( IC.GT.0 ) THEN
              QTREE%QCHILD(IQ2,II) = NEWQ(IC)
           ELSE
              QTREE%QCHILD(IQ2,II) = 0
           END IF
         END DO
         DO II=1,4
           IC = QTREE%QNBR(IQ,II)
           IF ( IC.GT.0 ) THEN
              QTREE%QNBR(IQ2,II) = NEWQ(IC)
           ELSE
              QTREE%QNBR(IQ2,II) = 0
           END IF
         END DO
         IF ( IQ.GT.NQUAD2 ) THEN
            ! Set null values for quad indices that won't be overwritten:
            QTREE%QICELL(IQ,:) = 0
            QTREE%QLEVEL(IQ) = 0
            QTREE%QPARENT(IQ) = 0
            QTREE%QCHILD(IQ,:) = 0
            QTREE%QNBR(IQ,:) = 0
         END IF
      END DO
      !
      QTREE%NQUAD = NQUAD2
      !
      RETURN
      END SUBROUTINE QA_REMUNDEF
!
!/ ------------------------------------------------------------------- /
!
      SUBROUTINE QA_TSORDER( NTSUB, LVLMAX, LVLREF, LVLSTEP, IOPT )
!/
!/       Richard Gorman, NIWA
!/         23-Sep-2010:      Origination.
!          Jan-April 2014:   Rename
!          Nov 2017:         New method for UNO
!/
!  1. Purpose :
!
!     Order of time substeps in which different quadtree levels are
!     propagated.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NTSUB   Int.   I   Number of time substeps.
!       LVLMAX  Int.   I   Maximum quadtree level
!       LVLREF  Int.   I   Reference quadtree level
!       LVLSTEP I.A.   O   Sequence of levels
!       IOPT    Int.   I*  Option: =0 for "delayed" ordering,
!                                     where a step at level L is taken immediately
!                                     after 2 steps at level L+1
!                                  =1 for "centred" ordering
!                                     where a step at level L is taken mid way
!                                     between 2 steps at level L+1
!                                  default = "centred" (old method)
!     ----------------------------------------------------------------
!
!     Local variables.
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       None
!
!  5. Called by :
!
!       W3XYPX   Quadtree spatial propagation routine.
!
!  6. Error messages :
!
!       None.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NTSUB, LVLMAX, LVLREF
      INTEGER, INTENT(OUT)    :: LVLSTEP(:)
      INTEGER, OPTIONAL, INTENT(IN)     :: IOPT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER           ::  IT2, LVL, LVLDIF, ITSUB, MTSUB
      INTEGER           ::  LS, IOFF, K, KK
!/
!/ ------------------------------------------------------------------- /
!/
      MTSUB = SIZE(LVLSTEP,1)
      IF ( PRESENT (IOPT) ) THEN
!
!  New algorithm for either "delayed" (IOPT=0) or "centred" (IOPT=1) ordering
        IOFF = MAX(IOPT,0)
        IOFF = MIN(IOPT,1)
        LS = 1
        LVLSTEP(1) = LVLMAX
        DO LVL = LVLMAX-1,LVLREF,-1
          KK = (2-IOFF)*LS+1
          IF (KK.LE.MTSUB) LVLSTEP(KK) = LVL
          DO K=1,LS
            KK = LS + K + IOFF
            IF (KK.LE.MTSUB) LVLSTEP(KK) = LVLSTEP(K)
          END DO
          LS = 2*LS+1
          IF ( LS.GE.MTSUB ) EXIT
        END DO
      ELSE
!
!  Original algorithm for "centred" ordering
        DO ITSUB = 1, MIN ( NTSUB, MTSUB )
          IT2 = ITSUB - 1
          DO LVL = LVLREF, LVLMAX
            LVLDIF = LVLMAX - LVL
            IF ( BTEST(ITSUB,LVLDIF) .NEQV. BTEST(IT2,LVLDIF) ) THEN
              LVLSTEP(ITSUB) = LVL
              EXIT
            END IF
          END DO
        END DO
      END IF
!/
!/ End of QA_TSORDER ----------------------------------------------------- /
!/
      END SUBROUTINE QA_TSORDER
!/
!/ ------------------------------------------------------------------- /
!/
      SUBROUTINE QA_XY2CELL ( QTREE, XTARG, YTARG, ICELL, XCELL,      &
                              YCELL, LVCELL, IQUAD, ISUB, ISTAT )
!/
!/       Richard Gorman, NIWA
!/         March, 2013:      Origination.
!          Jan-April 2014:   Rewrite using QTREE structure
!
!  1. Purpose :
!
!      Given the (reference grid) coordinates of a point, locate the cell 
!      (if any) in a threaded quadtree structure containing it.
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       QTREE   QA_TREE I  Quadtree structure, including components:
!          NX0, NY0 Int.   Size of level-0 grid
!          LVLREF  Int.    Reference refinement level
!          LVLMAX  Int.    Max. number of refinement levels
!          QISEA   I.A.    Sea-point index for the centre of the quad (if it
!                          has no children), and/or cells
!          QCHILD  I.A.    Quad indices for the children of each quad
!       XTARG   Real   I   target X coordinate (reference grid coordinates)
!       YTARG   Real   I   target Y coordinate (reference grid coordinates)
!       ICELL   Int.   O   Index of the cell containing the target coordinates
!       XCELL   Real   O   cell X coordinate (reference grid coordinates)
!       YCELL   Real   O   cell Y coordinate (reference grid coordinates)
!       LVCELL  Int.   O   Level of the cell containing the target coordinates
!       ISTAT   I.A.   O   Return flags: 
!		 ISTAT(1:2) = -1  no cell found (the target point is in a 
!                                 part of the domain not covered by a cell)
!                ISTAT(1:2) =  0  cell found containing (XTARG, YTARG)
!                ISTAT(1) =  1(2)  (XTARG, YTARG) is outside the level-0 grid
!                                 in the -ve (+ve) x direction
!                                 (W, E) wall, respectively
!                ISTAT(2) =  1(2)  (XTARG, YTARG) is outside the level-0 grid
!                                 in the -ve (+ve) y direction
!                                 returned ISEA is for the nearest point on the
!                                 (S, N) wall, respectively
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     QA_MLG2QT  Subr. qa_utils Derive a quadtree structure from multiple 
!                               regular grids
!     QA_INTERP  Subr. qa_utils Interpolate to a target location 
!                               using weighted values from 
!                               neighbouring cells.
!     QA_Q2QMAP  Subr. qa_utils Mapping from one quadtree structure to
!                               another
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      TYPE(QA_TREE), INTENT(IN) :: QTREE
      REAL, INTENT(IN)        :: XTARG, YTARG
      REAL, INTENT(OUT)       :: XCELL, YCELL
      INTEGER, INTENT(OUT)    :: ICELL, LVCELL, IQUAD, ISUB, ISTAT(2)

!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER          :: IX, IY, IQ, IS, IC, ILVL
      !INTEGER          :: ICMIN, LVMIN
      REAL             :: EPS, X, Y, DELXY, DELXY0, EPS2, D2, XC, YC
      !REAL             :: XMIN, YMIN, D2MIN
!/
!/ ------------------------------------------------------------------- /
!/
!
      ISTAT = 0
      EPS = 0.999*2.**(QTREE%LVLREF - QTREE%LVLMAX - 1)
      EPS2 = EPS**2
      DELXY0 = 2.**QTREE%LVLREF
      !D2MIN = (2.*DELXY0)**2
      !ICMIN = 0
      !XMIN = 0.
      !YMIN = 0.
      ICELL = 0
      IQUAD = 0
      ISUB = 0
!
!  Identify the level-0 quad to start searching in
      X = XTARG
      Y = YTARG
      IX = NINT((X-0.5)/DELXY0 + 0.5)
      IY = NINT((Y-0.5)/DELXY0 + 0.5)
      IF ( IX.LT.1 ) THEN
         !  Target point is outside the domain to the west
         !  Flag that we are searching for a point just inside the W boundary
         X = 0.5 + 0.1*EPS
         IX = 1
         ISTAT(1) = 1
      ELSEIF ( IX.GT.QTREE%NX0 ) THEN
         !  Target point is outside the domain to the east
         !  Flag that we are searching for a point just inside the E boundary
         X = QTREE%NX0 + 0.5 - 0.1*EPS
         IX = QTREE%NX0
         ISTAT(1) = 2
      END IF
      IF ( IY.LT.1 ) THEN
         !  Target point is outside the domain to the south
         !  Flag that we are searching for a point just inside the S boundary
         Y = 0.5 + 0.1*EPS
         IY = 1
         ISTAT(2) = 1
      ELSEIF ( IY.GT.QTREE%NY0 ) THEN
         !  Target point is outside the domain to the north
         !  Flag that we are searching for a point just inside the N boundary
         Y = QTREE%NY0 + 0.5 - 0.1*EPS
         IY = QTREE%NY0
         ISTAT(2) = 2
      END IF
      IQ = IX + QTREE%NX0*(IY-1)
      XC = (IX-1)*DELXY0 + 0.5*(1. + DELXY0)
      YC = (IY-1)*DELXY0 + 0.5*(1. + DELXY0)
      IQUAD = IQ
      ISUB = 0
      LVCELL = 0
!
      DO ILVL=1, QTREE%LVLMAX
         IC = QTREE%QICELL(IQ,0)
         IF ( IC.GT.0 ) THEN
            ! The (level-0) quad is a leaf cell: see how close this is
            ICELL = IC
            LVCELL = ILVL-1
            XCELL = XC
            YCELL = YC
            D2 = (X-XC)**2 + (Y-YC)**2 
            ! As close as possible: seek no further
            IF ( D2.LT.EPS2 ) EXIT
            ! Otherwise, see if it is closest yet:
            !IF ( D2.LE.1.001*D2MIN ) THEN
            ! No: take the higher level refinement, even if it is no closer
               !ICMIN = ICELL
               !LVMIN = LVCELL
               !D2MIN = D2
               !XMIN = XCELL
               !YMIN = YCELL
            !END IF
         END IF
         DELXY = DELXY0*2.**(-ILVL-1)
         !
         ! decide which sub-quad to search, based on displacement of
         ! the target point from the quad centre:
         IF ( X.LT.XC .AND. Y.LT.YC ) THEN
            IS = 1
            XC = XC - DELXY
            YC = YC - DELXY
         ELSEIF ( X.GE.XC .AND. Y.LT.YC ) THEN
            IS = 2
            XC = XC + DELXY
            YC = YC - DELXY
         ELSEIF ( X.LT.XC .AND. Y.GE.YC ) THEN
            IS = 3
            XC = XC - DELXY
            YC = YC + DELXY
         ELSEIF ( X.GE.XC .AND. Y.GE.YC ) THEN
            IS = 4
            XC = XC + DELXY
            YC = YC + DELXY
         END IF
         !
         IC = QTREE%QICELL(IQ,IS)
         IF ( IC.GT.0 ) THEN
            ! The relevant sub-quad contains a leaf cell: see how close this is
            ICELL = IC
            ISUB = IS
            XCELL = XC
            YCELL = YC
            LVCELL = ILVL
            D2 = (X-XC)**2 + (Y-YC)**2 
            ! As close as possible: seek no further
            IF ( D2.LT.EPS2 ) EXIT
            ! Otherwise, see if it is closest yet:
            !IF ( D2.LE.1.001*D2MIN ) THEN
            ! No: take the higher level refinement, even if it is no closer
               !ICMIN = ICELL
               !LVMIN = LVCELL
               !D2MIN = D2
               !XMIN = XC
               !YMIN = YC
            !END IF
         END IF
         !
         IQ = QTREE%QCHILD(IQ,IS)
         IF ( IQ.EQ.0 ) THEN
            !  The quad has no child quad in the relevant quadrant
            !IF ( ICMIN.EQ.0 ) THEN
               ! No candidate cells found.
               !  Flag as a failure and exit (perhaps later we can look for
               !  the nearest point in the mesh)
               !ISTAT = -1
               !ICELL = 0
            !ELSE
               ! Return the closest candidate cell:
               !ICELL = ICMIN
               !LVCELL = LVMIN
               !XCELL = XMIN
               !YCELL = YMIN
               !IQUAD = IQMIN
               !ISUB = ISMIN
            !END IF
            IF (ICELL.EQ.0) ISTAT = -1
            EXIT
         ELSE
            IQUAD = IQ
         END IF
      END DO
      RETURN
!/
!/ End of QA_XY2CELL ----------------------------------------------------- /
!/
      END SUBROUTINE QA_XY2CELL
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER FUNCTION ML_CHILD ( IML, NX0, NY0 )
!/
!         Richard Gorman, NIWA, September, 2010
!
!  1. Purpose :
!
!      Return an index into multilevel arrays for the first (SW) child
!      of a cell in a quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IML        Int.   multilevel index of the input cell
!       NX0,NY0    Int.   Size of level-0 rectangular grid
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IML, NX0, NY0
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: IX, JY, IXY, LEVEL, MLIMIN, NXY, NX, NY,     &
                         MLIMAX, IXC, JYC, IXYC
!  Corresponding multilevel index:
      LEVEL = 0
      MLIMIN = 0
      NXY = NX0*NY0
      NX = NX0
      NY = NY0
      MLIMAX = NXY
      DO WHILE ( MLIMAX.LT.IML )
        MLIMIN = MLIMAX
        LEVEL = LEVEL + 1
        NX = 2*NX
        NY = 2*NY
        NXY = NXY*4
        MLIMAX = MLIMAX + NXY
      END DO
!
!  Indices in the grid at the same level as this cell:
      IXY = IML - MLIMIN
      JY = (IXY-1)/NX + 1
      IX = IXY - NX*(JY-1)
!
!  Indices in the grid at the child level:
      IXC = (IX-1)*2 + 1
      JYC = (JY-1)*2 + 1
!
!  Corresponding 1D cell index:
      IXYC = (JYC - 1)*NX*2 + IXC
!
!  Corresponding multilevel index:
      ML_CHILD = MLIMAX + IXYC
!
      RETURN
!/
!/ ------------------------------------------------------------------- /
      END FUNCTION ML_CHILD
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER FUNCTION ML_INDEX ( XI, YJ, LEVEL, LVLMAX, LVLREF,      &
                                  NSUBGRID, NSGL, XYRANGE, NXSGL,     &
                                  MLOFFSET, INDSGL )
!/
!         Richard Gorman, NIWA, April, 2012
!
!  1. Purpose :
!
!      Return an index into multilevel arrays for a cell in a quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       XI,YJ      Real   Cell coordinates relative to the highest level grid
!       LEVEL      Int.   Level of the cell
!       LVLMAX     Int.   Maximum level
!       LVLREF     Int.   Reference level
!       NSUBGRID   Int.   Number of subgrids under the reference grid
!       NSGL       I.A.   Number of subgrids that include each level
!       XYRANGE    R.A.   Extent of each subgrid [Xmin Xmax Ymin Ymax]
!                           (in reference grid coordinates)
!       NXSGL      I.A.   Number of cells in the X direction for each level
!                         of each subgrid 
!       MLOFFSET   I.A.   Offset for the multilevel index for each level of 
!                         each subgrid 
!       INDSGL     I.A.   List of subgrids having a layer at a given level
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: XI, YJ
      INTEGER, INTENT(IN)     :: LEVEL, LVLMAX, LVLREF, NSUBGRID
      INTEGER, INTENT(IN)     :: NSGL(0:LVLMAX)
      REAL, INTENT(IN)        :: XYRANGE(0:NSUBGRID,4)
      INTEGER, INTENT(IN)     :: NXSGL(0:NSUBGRID,0:LVLMAX)
      INTEGER, INTENT(IN)     :: MLOFFSET(0:NSUBGRID,0:LVLMAX)
      INTEGER, INTENT(IN)     :: INDSGL(NSUBGRID+1,0:LVLMAX)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: IX, JY, IXY, INDSG, K, ISG
      REAL            :: DELXYINV
!
!  Identify which subgrid this point is in by looping through all the 
!  subgrids which include the target level, and checking if X, Y is in
!  the spatial range:
!
      INDSG = -1
      DO K=1,NSGL(LEVEL)
         ISG = INDSGL(K,LEVEL)
         IF ( XI.GE.XYRANGE(ISG,1) .AND. XI.LE.XYRANGE(ISG,2) .AND.   &
              YJ.GE.XYRANGE(ISG,3) .AND. YJ.LE.XYRANGE(ISG,4) ) THEN
            INDSG = ISG
            EXIT
         END IF
      END DO
!
! If none found, return with ML_INDEX = 0      
      IF ( INDSG.LT.0 ) THEN
         ML_INDEX = 0
         RETURN
      END IF
!
!  Grid spacing (relative to the reference level) in the rectangular
!  grid at this level of the subgrid      
      DELXYINV = 2.**(LEVEL-LVLREF)
!
!  Indices in the rectangular grid at the same level as this cell:
      IX = NINT( (XI - XYRANGE(INDSG,1))*DELXYINV + 0.5 )
      JY = NINT( (YJ - XYRANGE(INDSG,3))*DELXYINV + 0.5 )
!
!  Corresponding 1D cell index:
      IXY = (JY - 1)*NXSGL(INDSG,LEVEL) + IX
!
!  Corresponding multilevel index:
      ML_INDEX = MLOFFSET(INDSG,LEVEL) + IXY
!
      RETURN
!/
!/ ------------------------------------------------------------------- /
      END FUNCTION ML_INDEX
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER FUNCTION ML_PARENT ( IML, NX0, NY0 )
!/
!         Richard Gorman, NIWA, September, 2010
!
!  1. Purpose :
!
!      Return an index into multilevel arrays for the parent of a cell
!      in a quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IML        Int.   multilevel index of the input cell
!       NX0,NY0    Int.   Size of level-0 rectangular grid
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IML, NX0, NY0
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: IX, JY, IXY, LEVEL, MLIMIN, NXY, NX, NY,     &
                         MLIMINP, MLIMAX, NXP, NYP, IXP, JYP, IXYP
!  Corresponding multilevel index:
      LEVEL = 0
      MLIMIN = 0
      NXY = NX0*NY0
      NX = NX0
      NY = NY0
      MLIMAX = NXY
      DO WHILE ( MLIMAX.LT.IML )
        NXP = NX
        NYP = NY
        MLIMINP = MLIMIN
        MLIMIN = MLIMAX
        LEVEL = LEVEL + 1
        NX = 2*NX
        NY = 2*NY
        NXY = NXY*4
        MLIMAX = MLIMAX + NXY
      END DO
      IF (LEVEL.EQ.0) THEN
        ML_PARENT = 0
        RETURN
      END IF
!
!  Indices in the grid at the same level as this cell:
      IXY = IML - MLIMIN
      JY = (IXY-1)/NX + 1
      IX = IXY - NX*(JY-1)
!
!  Indices in the grid at the parent level:
      IXP = (IX-1)/2 + 1
      JYP = (JY-1)/2 + 1
!
!  Corresponding 1D cell index:
      IXYP = (JYP - 1)*NXP + IXP
!
!  Corresponding multilevel index:
      ML_PARENT = MLIMINP + IXYP
!
      RETURN
!/
!/ ------------------------------------------------------------------- /
      END FUNCTION ML_PARENT
!
!/ ------------------------------------------------------------------- /
!
      INTEGER FUNCTION MULTILEVEL_INDEX ( NX0, NY0, XI, YJ, LEVEL, LVLREF )
!/
!         Richard Gorman, NIWA, October, 2009
!
!  1. Purpose :
!
!      Return an index into multilevel arrays for a cell in a quadtree structure
!
!  2. Method :
!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NX0,NY0    Int.   Size of level-0 rectangular grid
!       XI,YJ      Real   Cell coordinates relative to the reference grid
!       LEVEL      Int.   Level of the cell
!       LVLREF     Int.   Reference level, at which the grid spacing is 1
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NX0, NY0, LEVEL, LVLREF
      REAL, INTENT(IN)        :: XI, YJ
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER         :: IX, JY, IXY, NX, NY
      REAL            :: STEPINV
!
!  Inverse grid step at the same level as this cell:
      STEPINV = 2.**(LEVEL - LVLREF)
!
!  Indices in the grid at the same level as this cell:
      IX = NINT( (XI - 0.5)*STEPINV + 0.5 )
      JY = NINT( (YJ - 0.5)*STEPINV + 0.5 )
!
!  Grid size at this level
      NX = NX0*2**LEVEL
      NY = NY0*2**LEVEL
      IF ( IX.LT.1 .OR. IX.GT.NX .OR. JY.LT.1 .OR. JY.GT.NY ) THEN
         MULTILEVEL_INDEX = 0
      ELSE
!
!  Corresponding 1D cell index:
         IXY = (JY - 1)*NX0*2**LEVEL + IX
!
!  Corresponding multilevel index:
         MULTILEVEL_INDEX = NX0*NY0*NINT( (4.**LEVEL - 1)/3. ) + IXY
      END IF
!
      RETURN
!/
!/ ------------------------------------------------------------------- /
      END FUNCTION MULTILEVEL_INDEX
!/ ------------------------------------------------------------------- /
!
recursive subroutine QSort(A,ind,na)
!/
!         Richard Gorman, NIWA, February, 2014
!
!  1. Purpose :
!
!      Sort an array A into ascending order, and also return a ranking
!      index into the original array, i.e.
!              Asort = A
!              call Qsort(Asort,Ind,size(A,1))
!      returns Asort in ascending order, and Asort(i) = A(ind(i))
!
!  2. Method :
!
!      Quicksort algorithm
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A          R.A    I/O   Array to be sorted
!       Ind        I.A    I/O   Indices into original order of A
!       na         Int.    I    Length of A and Ind
!
!  4. Subroutines used :
!
!     None
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ---------------------------------------------------------------- 
!     ----------------------------------------------------------------
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
!/ ------------------------------------------------------------------- /
!/
IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
! DUMMY ARGUMENTS
integer, intent(in) :: nA
real, intent(in out) :: A(nA)
integer, intent(in out) :: Ind(nA)
 
! LOCAL VARIABLES
integer :: left, right
real :: random
real :: pivot
real :: temp
integer :: itemp
integer :: marker
!
! Check that ind(1) is valid: if not, initialise it
if ( ind(1).lt.1 ) then
    do left=1,nA
        Ind(left) = left
    end do  
end if 
! For a length-1 array, there is nothing to do 
if (nA > 1) then

    call random_number(random)
    pivot = A(int(random*real(nA-1))+1)   ! random pivot (not best performance, 
                                          ! but avoids worst-case)
    left = 0
    right = nA + 1
    !
    ! assign elements into left and right partitions,
    ! according to whether A is <= or > pivot, respecively
    !
    ! Work through the whole array from both ends simultaneously
    do while (left < right)
        ! 
        right = right - 1
        do while (A(right) > pivot)
            right = right - 1
        end do
        left = left + 1
        do while (A(left) < pivot)
            left = left + 1
        end do
        ! Searching from the right found an element <= pivot, and
        ! searching from the left we have reached an element >= pivot.
        ! If the two searches haven't met, swap the elements into the correct partitions,
        ! and carry on.
        ! Otherwise this partitioning is complete and we can exit the loop. 
        if (left < right) then
            ! 
            temp = A(left)
            A(left) = A(right)
            A(right) = temp
            itemp = Ind(left)
            Ind(left) = Ind(right)
            Ind(right) = itemp
        end if
    end do
    !
    ! Mark the start of the right partiion
    if (left == right) then
        marker = left + 1
    else
        marker = left
    end if
    !
    ! Apply the same process to each of these left and right partitions:
    call QSort(A(:marker-1),Ind(:marker-1),marker-1)
    call QSort(A(marker:),ind(marker:),nA-marker+1)
 
end if
 
end subroutine QSort
 
!/
!/ End of module QA_UTILS ----------------------------------------------------- /
!/
      END MODULE QA_UTILS
