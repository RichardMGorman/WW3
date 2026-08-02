#include "w3macros.h"
!/ ------------------------------------------------------------------- /
  MODULE W3ADGRMD
  !/ ------------------------------------------------------------------- /
  !/
  !/ Required modules
  !/
  USE QA_UTILS
  !/
  !/ Specify default accessibility
  !/
  PUBLIC
  !/
CONTAINS
  !/
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3ADGR ( DIAG_TYPE, ADAPT_VARS, NDSE )
  !/
  !/    Richard Gorman, NIWA, May 2014
  !/
  !  1. Purpose :
  !
  !     Adapt a quadtree grid
  !
  !  2. Method :
  !
  !
  !  3. Parameters :
  !
  !     Parameter list
  !     ----------------------------------------------------------------
  !       DIAG_TYPE  Int.  I   Code for diagnostic variable selected:
  !				1 for depth-dependent CFL factor
  !				2 for 2nd derivative of wind speed
  !				12 for the product of types 1 & 2
  !       ADAPT_VARS L.A.  I   Flags for sets of arrays to recompute for 
  !                            adapted grid
  !                              ADAPT_VARS(1)  Adapt the quadtree grid
  !                              ADAPT_VARS(2)  Recompute geographic grid arrays
  !                              ADAPT_VARS(3)  Recompute WN and CG
  !                              ADAPT_VARS(4)  Recompute spectra
  !       NDSE       Int.  I   Unit number for error output
  !     ----------------------------------------------------------------
  !
  !     Local data
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      QA_ADGR    Subr. QA_UTILS  Generic grid adaptation.
  !      QA_ADVAR   Subr. QA_UTILS  Revision of a data array after adaptation.
  !      EXTCDE     Subr. W3SERVMD  Exit on error.
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      WSTRT     Prog. W3_STRT  Wave model initial conditions.
  !     ----------------------------------------------------------------
  !
  !  6. Error messages :
  !
  !
  !  7. Remarks :
  !
  !
  !  8. Structure :
  !
  !     -------------------------------------------
  !     -------------------------------------------
  !
  !  9. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !       !/T     Test output.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  USE CONSTANTS
  USE QA_UTILS, ONLY: QA_ADGR, QA_ADVAR
  USE W3GDATMD, ONLY: NSEA, ZB, IQGW, IQGB, WTS1_QA,           &
                      DIAGVAR_QA, DVMAX_QA, DVTOLFAC_QA,       &
                      NCTARGET_QA,  NMOD_QA, NFROM_QA,         &
                      IFROM_QA, WFROM_QA, ITO_QA, NAUX_QA,     &
                      IAUX_QA, EXAUX_QA, DATB_QA, NSPEC,       &
                      QTREE, WTS1_QA, NK, DMIN, MAPML_QA,      &
                      NCMXQ
#ifdef W3_UNO                      
  USE W3GDATMD, ONLY: NUFc, NVFc, NRLv, IJKCel, IJKUFc, IJKVFc,    &
                      NLvUFc, NLvVFc   
#endif
  USE W3WDATMD, ONLY: VA
  USE W3ADATMD, ONLY: U10, CG, WN
  USE W3SERVMD, ONLY: EXTCDE
  !USE W3DISPMD, ONLY: WAVNU1
  USE W3ODATMD, ONLY: IAPROC, NAPERR
#ifdef W3_S
  USE W3SERVMD, ONLY: STRACE
#endif
  !
  IMPLICIT NONE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  INTEGER, INTENT(IN)     :: DIAG_TYPE
  LOGICAL, INTENT(IN)     :: ADAPT_VARS(:)
  INTEGER, INTENT(IN)     :: NDSE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  INTEGER                 :: ISEA, IMOD, ISPEC, IK, NSEA1
  REAL                    :: DPT
  INTEGER                 :: IERR, IUN, ICASE, NSOL
  REAL, ALLOCATABLE       :: FIELD(:), CFLFAC(:), D2WIND(:)
  INTEGER, ALLOCATABLE    :: ISOL(:)
  REAL, ALLOCATABLE       :: DERVAL(:)
  !/
  !/ ------------------------------------------------------------------- /
  !/
#ifdef W3_s
  CALL STRACE (IENT, 'W3ADGR')
#endif
  !
  ! 0. Initialisation:
  !
  !  Constants:
  IUN = 0
  IF ( IAPROC .EQ. NAPERR ) IUN = NDSE
  !
  IF ( SIZE(ADAPT_VARS,1).LT.1 ) then
     IF ( IUN.GT.0 ) THEN
        WRITE (IUN,*) ' ERROR IN W3ADGRMD'
        WRITE (IUN,*) ' ADAPT_VARS IS EMPTY'
     END IF
     CALL EXTCDE ( 1 )
  END IF
  IF ( ADAPT_VARS(1) ) THEN
     IF ( DIAG_TYPE.NE.1 .AND. DIAG_TYPE.NE.2 .AND.               &
          DIAG_TYPE.NE.12 .AND. DIAG_TYPE.NE.13 ) THEN
        IF ( IUN.GT.0 ) THEN
           WRITE (IUN,*) ' ERROR IN W3ADGR'
           WRITE (IUN,*) ' INVALID DIAG_TYPE = ', DIAG_TYPE
           WRITE (IUN,*) ' VALID OPTIONS ARE: '
           WRITE (IUN,*) '     1: depth-dependent CFL factor'
           WRITE (IUN,*) '     2: 2nd derivative of wind speed'
           WRITE (IUN,*) '     3: wind speed'
           WRITE (IUN,*) '    12: product of types 1 & 2'
           WRITE (IUN,*) '    13: product of types 1 & 3'
        END IF
        CALL EXTCDE ( 1 )
     END IF
  !
  ! 1. Adapt the quadtree grid      
  !
  !        1.2. Compute the selected diagnostic variable
  !
     DIAGVAR_QA = 0.
     NSEA1 = QTREE(IQGW)%NCELL
  !         
     IF ( DIAG_TYPE.EQ.1 .OR. DIAG_TYPE.EQ.12 .OR.                &
          DIAG_TYPE.EQ.13 ) THEN
  !             Use "CFL factor" as a diagnostic variable
        ALLOCATE ( CFLFAC(NSEA1) )
        CFLFAC = 0.
  !             Loop through cells
        DO ISEA=1,NSEA1
           IF ( QTREE(IQGW)%CELL_TYPE(ISEA).EQ.                   &
                QTREE(IQGW)%UNDEF_TYPE ) CYCLE
           DPT = MAX( DMIN, -ZB(ISEA) )
           IF ( DPT.LE.0. ) CYCLE
           !CALL WAVNU1 ( SIG(1), DPT, WN1, CG1 )
           CFLFAC(ISEA) = (2.**(-QTREE(IQGW)%INDLVL(ISEA)))       &
                          /SQRT(DPT)
           DIAGVAR_QA(ISEA) = CFLFAC(ISEA)
        END DO
     END IF
  !
     IF ( DIAG_TYPE.EQ.2 .OR. DIAG_TYPE.EQ.12 ) THEN
  !             Use 2nd deriv of wind speed as a diagnostic variable
        ALLOCATE ( D2WIND(NSEA1) )
        D2WIND = 0.
        !NSOL = 5
        !ISOL(1) = 1
        !ISOL(2) = 2
        !ISOL(3) = 3
        !ISOL(4) = 4
        !ISOL(5) = 5
        ALLOCATE ( ISOL(2) )
        ALLOCATE ( DERVAL(2) )
        NSOL = 2
        ISOL(1) = 3
        ISOL(2) = 4
  !             Loop through cells
        DO ISEA=1,NSEA1
           IF ( QTREE(IQGW)%CELL_TYPE(ISEA).EQ.                   &
                QTREE(IQGW)%UNDEF_TYPE ) CYCLE
           ICASE = QTREE(IQGW)%NCASE(ISEA)
           CALL QA_DERIVC( ISEA, ICASE, NSOL, ISOL,             &
                       WTS1_QA%NNGBR, QTREE(IQGW)%NGBR(ISEA,:), &
                       U10, WTS1_QA%DERWTC, DERVAL, ierr=IERR,  &
                       ndse=IUN )
           IF ( IERR.GT.0 ) THEN
              IF ( IUN.GT.0 ) THEN
                 WRITE (IUN,*) ' ERROR IN W3ADGR'
                 WRITE (IUN,*) ' WHILE CALLING QA_DERIVC'
                 WRITE (IUN,*) '    RETURNED IERR = ', IERR
              END IF
              CALL EXTCDE ( 2 )
           END IF
           D2WIND(ISEA) = ABS(DERVAL(1)) + ABS(DERVAL(2))
           DIAGVAR_QA(ISEA) = D2WIND(ISEA)
        END DO
     END IF
     !
     IF ( DIAG_TYPE.EQ.3 ) DIAGVAR_QA(1:NSEA1) = U10(1:NSEA1)
     IF ( DIAG_TYPE.EQ.12 ) DIAGVAR_QA(1:NSEA1) =                 &
                                CFLFAC(1:NSEA1)*D2WIND(1:NSEA1)
     IF ( DIAG_TYPE.EQ.13 ) DIAGVAR_QA(1:NSEA1) =                 &
                                CFLFAC(1:NSEA1)*U10(1:NSEA1)
  !
  !        1.3. Adapt the grid
  !
     IF ( NAUX_QA.GT.1 ) THEN
        CALL QA_ADGR( QTREE(IQGW), WTS1_QA, DIAGVAR_QA, DVMAX_QA, &
                    DVTOLFAC_QA, NCTARGET_QA, NMOD_QA, NFROM_QA,  &
                    IFROM_QA, WFROM_QA, ITO_QA,                   &
                    !qspare = QSPARE,                              &
                    qtree_aux=QTREE(2:NAUX_QA),                   &
                    ind_aux=IAUX_QA(:,2:NAUX_QA),                 &
                    exact_aux=EXAUX_QA(:,2:NAUX_QA),              & 
                    no_adapt_type=(/0,-1/), ndse=IUN, ierr=IERR )
  !                        no_adapt_type=(/0,-1/), ndse=IUN, ierr=IERR,  &
  !                        mapml=MAPML_QA )
     ELSE
        CALL QA_ADGR( QTREE(IQGW), WTS1_QA, DIAGVAR_QA, DVMAX_QA, &
                   DVTOLFAC_QA, NCTARGET_QA, NMOD_QA, NFROM_QA,   &
                   IFROM_QA, WFROM_QA, ITO_QA,                    &
                   !qspare = QSPARE,                               &
                   no_adapt_type=(/0,-1/), ndse=IUN, ierr=IERR )
  !                       no_adapt_type=(/0,-1/), ndse=IUN, ierr=IERR,   &
  !                       mapml=MAPML_QA )
     END IF
     IF ( IERR.GT.0 ) THEN
        IF ( IUN.GT.0 ) THEN
           WRITE (IUN,*) ' ERROR IN W3ADGR'
           WRITE (IUN,*) ' WHILE CALLING QA_ADGR'
           WRITE (IUN,*) '    RETURNED IERR = ', IERR
        END IF
        CALL EXTCDE ( 3 )
     END IF
     NSEA = QTREE(IQGW)%NCELL
  !
#ifdef W3_UNO  
  !  Recompute maps for propagation
     CALL QA_QT2SMC ( QTREE(IQGW), -9, IJKCel, IJKUFc, IJKVFc,     &
                      NLvUFc, NLvVFc, ierr=IERR, ndse=IUN )
     IF ( IERR.GT.0 ) THEN
        IF ( IUN.GT.0 ) THEN
           WRITE (IUN,*) ' ERROR IN W3ADGR'
           WRITE (IUN,*) ' WHILE CALLING QA_QT2SMC'
           WRITE (IUN,*) '    RETURNED IERR = ', IERR
        END IF
        CALL EXTCDE ( 3 )
     END IF
#endif  
  END IF
  !
  IF ( SIZE(ADAPT_VARS,1).LT.2 ) RETURN
  IF ( ADAPT_VARS(2) .AND. NAUX_QA.GE.2 ) THEN
  !
  ! 2. Recompute Geographic grid arrays on the adapted grid
  !
     DO IMOD=1,NMOD_QA
        ISEA = ITO_QA(IMOD)
        CALL W3ADG1 ( ISEA )
     END DO
  END IF
  !
  IF ( SIZE(ADAPT_VARS,1).LT.3 ) RETURN
  ALLOCATE ( FIELD(NCMXQ(IQGW)) )
  IF ( ADAPT_VARS(3) ) THEN
  !
  ! 3. Recompute CG and WN on the adapted grid, by interpolation
  !    These will be redone (for shallow water) by W3ULEV, but
  !    that needs an old value for each cell, which this can provide
  !
     DO IK=0,NK+1
        FIELD(1:NSEA1) = CG(IK,1:NSEA1)
        CALL QA_ADVAR ( NMOD_QA, NFROM_QA, IFROM_QA, WFROM_QA,    &
                        ITO_QA, FIELD )
        CG(IK,1:NSEA) = FIELD(1:NSEA)
        FIELD(1:NSEA1) = WN(IK,1:NSEA1)
        CALL QA_ADVAR ( NMOD_QA, NFROM_QA, IFROM_QA, WFROM_QA,    &
                        ITO_QA, FIELD )
        WN(IK,1:NSEA) = FIELD(1:NSEA)
     END DO
  END IF
  !
  IF ( SIZE(ADAPT_VARS,1).LT.4 ) RETURN
  IF ( ADAPT_VARS(4) ) THEN
  !
  ! 4. Recompute spectra on the adapted grid
  !
#ifdef W3_MPI  
  IF ( IUN.GT.0 ) THEN
     WRITE (IUN,*) ' ERROR IN W3ADGR'
     WRITE (IUN,*) ' MPI CALLS TO ADAPT SPECTRA NOT YET ENABLED'
  END IF
  CALL EXTCDE ( 1 )
#endif
     DO ISPEC=1, NSPEC
        FIELD(1:NSEA1) = VA(ISPEC,1:NSEA1)
        CALL QA_ADVAR ( NMOD_QA, NFROM_QA, IFROM_QA, WFROM_QA,    &
                        ITO_QA, FIELD )
        VA(ISPEC,1:NSEA) = FIELD(1:NSEA)
     END DO
  END IF
  !      
  RETURN
  !/
  !/ End of W3ADGR ----------------------------------------------------- /
  !/
  END SUBROUTINE W3ADGR
  !
  !/ ------------------------------------------------------------------- /
  !
  SUBROUTINE W3ADG1 ( ISEA )
  !/
  !/    Richard Gorman, NIWA, May 2014
  !/
  !  1. Purpose :
  !
  !     Recompute auxiliary spatial data for one cell of a quadtree grid
  !
  !  2. Method :
  !
  !
  !  3. Parameters :
  !
  !     Parameter list
  !     ----------------------------------------------------------------
  !       ISEA       Int.  I   Index of cell to be recalculated
  !     ----------------------------------------------------------------
  !
  !     Local data
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3ADGR     Prog. W3ADGRMD  Adapting the wave quadtree
  !      W3INIT     Prog. W3INITMD  Wave model initialisation
  !     ----------------------------------------------------------------
  !
  !  6. Error messages :
  !
  !
  !  7. Remarks :
  !
  !   The following arrays are not recomputed, as they are not used in the 
  !   quadtree case: 
  !               
  !     DXDP(MY,MX), DXDQ(MY,MX), DYDP(MY,MX), DYDQ(MY,MX), 
  !     DPDX(MY,MX), DPDY(MY,MX), DQDX(MY,MX), DQDY(MY,MX), 
  !     GSQRT(MY,MX), HPFAC(MY,MX), HQFAC(MY,MX)
  !
  !    Used for /BT4, presently unsupported for quadtrees:
  !      SED_D50(0:MSEA), SED_PSIC(0:MSEA)
  !
  !    Used for /SMC, presently unsupported for quadtrees:
  !      NLvCel(0:MRLv), NLvUFc(0:MRLv), NLvVFc(0:MRLv), 
  !      IJKCel(5, -9:MCel), IJKUFc(7,MUFc), IJKVFc(8,MVFc),
  !      CTRNX(-9:MCel), CTRNY(-9:MCel), CLATF(MVFc) 
  !
  !    Used for /ARC, presently unsupported for quadtrees:
  !      ICLBAC(MBAC), ANGARC(MARC), SPCBAC(MSPEC,MBAC) )
  !
  !    Used for /REF1, presently unsupported for quadtrees:
  !      REFLC(4,0:NSEA)), REFLD(6,0:NSEA)
  !
  !  8. Structure :
  !
  !     -------------------------------------------
  !     -------------------------------------------
  !
  !  9. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !       !/T     Test output.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  USE CONSTANTS
  USE QA_UTILS
  USE W3GDATMD, ONLY  :  TRFLAG, MAPSTA, MAPST2, MAPFS,           &
                         MAPSF, SX, SY, X0, Y0, ZB, CLATS,        &
                         CLATIS, CTHG0S, TRNX, TRNY, XGRD, YGRD,  &
                         FLAGLL, IQGB, IQGW, QTREE, NAUX_QA,      &
                         IAUX_QA, EXAUX_QA, DATB_QA
#ifdef W3_S
  USE W3SERVMD, ONLY: STRACE
#endif
  !
  IMPLICIT NONE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  INTEGER, INTENT(IN)     :: ISEA
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  INTEGER                 :: IX, ISEA2, MAPT
  REAL                    :: X, Y
  !/
  !/ ------------------------------------------------------------------- /
  !/
#ifdef W3_S
  CALL STRACE (IENT, 'W3ADG1')
#endif
  !
  !
  ! 2. Recompute Geographic grid arrays on the adapted grid
  !
  IX = ISEA
  MAPT = QTREE(IQGW)%CELL_TYPE(ISEA)
  MAPSTA(1,IX) = MOD(MAPT+2,8) - 2
  MAPST2(1,IX) = (MAPT-MAPSTA(1,IX)) / 8
  MAPFS(1,IX)  = ISEA
  MAPSF(ISEA,1) = IX
  MAPSF(ISEA,2) = 1
  MAPSF(ISEA,3) = IX
  !
  ! Defaults:
  !
  ZB(ISEA) = 0.
  TRNX(1,IX) = 1.
  TRNY(1,IX) = 1.
  CLATS(ISEA)  = 1.
  CLATIS(ISEA) = 1.
  CTHG0S(ISEA) = 0.
  XGRD(1,IX) = 0.
  YGRD(1,IX) = 0.
  !
  ! Bail for undefined point
  !
  IF ( MAPT.EQ.QTREE(IQGW)%UNDEF_TYPE ) RETURN
  !
  X = X0 + (QTREE(IQGW)%XYVAL(ISEA,1)-1.)*SX
  Y = Y0 + (QTREE(IQGW)%XYVAL(ISEA,2)-1.)*SY
  XGRD(1,IX) = X
  YGRD(1,IX) = Y
  IF ( FLAGLL ) THEN
      CLATS(ISEA)    = COS(Y*DERA)
      CLATIS(ISEA)   = 1. / CLATS(ISEA)
      CTHG0S(ISEA)   = - TAN(DERA*Y) / RADIUS
    END IF
  IF ( IQGB.GT.0 .AND. IQGB.LE.NAUX_QA ) THEN
     ISEA2 = IAUX_QA(ISEA,IQGB)
     IF ( ISEA2.GT.0 ) THEN
         ZB(ISEA) = DATB_QA(ISEA2,1)
         IF ( TRFLAG.NE.0 ) THEN
             TRNX(1,IX) = DATB_QA(ISEA2,2)
             TRNY(1,IX) = DATB_QA(ISEA2,3)
           END IF
       END IF
    END IF
  RETURN
  !/
  !/ End of W3ADG1 ----------------------------------------------------- /
  !/
  END SUBROUTINE W3ADG1
  !
  !/ ------------------------------------------------------------------- /
  !
  SUBROUTINE W3QTGR ( A, DADX, DADY, NDSE )
  !/
  !/    Richard Gorman, NIWA, May 2014
  !/
  !  1. Purpose :
  !
  !     Compute X- and Y- gradients for one variable on a quadtree grid
  !
  !  2. Method :
  !
  !
  !  3. Parameters :
  !
  !     Parameter list
  !     ----------------------------------------------------------------
  !       A          R.A.  I   Input array on quadtree grid
  !       DADX       R.A.  O   dA/dx on quadtree grid
  !       DADY       R.A.  O   dA/dy on quadtree grid
  !     ----------------------------------------------------------------
  !
  !     Local data
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !     QA_DERIVC  Subr. qa_utils   computes (normalised) derivatives on
  !                                 a quadtree grid
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3WAVE     Prog. W3WAVEMD  Wave model
  !     ----------------------------------------------------------------
  !
  !  6. Error messages :
  !
  !
  !  7. Remarks :
  !
  !
  !  8. Structure :
  !
  !     -------------------------------------------
  !     -------------------------------------------
  !
  !  9. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  USE QA_UTILS
  USE W3GDATMD, ONLY  :  NSEA, IQGW, QTREE, WTS1_QA
  USE W3SERVMD, ONLY: EXTCDE
#ifdef W3_S
  USE W3SERVMD, ONLY: STRACE
#endif
  !
  IMPLICIT NONE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  REAL, INTENT(IN)     :: A(:)
  REAL, INTENT(OUT)    :: DADX(:,:)
  REAL, INTENT(OUT)    :: DADY(:,:)
  INTEGER, INTENT(IN)  :: NDSE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  INTEGER                 :: ISOL(2), NSOL
  INTEGER                 :: ISEA, ICASE, IERR
  REAL                    :: SCFAC
  REAL                    :: DERVAL(2)
  !/
  !/ ------------------------------------------------------------------- /
  !/
#ifdef W3_S
  CALL STRACE (IENT, 'W3QTGR')
#endif
  !
  !
  NSOL = 2
  ISOL(1) = 1
  ISOL(2) = 2
  !
  IF ( SIZE(A,1).LT.NSEA .OR. SIZE(DADX,2).LT.NSEA .OR.           &
       SIZE(DADY,2).LT.NSEA ) THEN
     IF ( NDSE.GT.0 ) THEN
        WRITE (NDSE,*) ' ERROR IN W3QTGR'
        WRITE (NDSE,*) ' INSUFFICIENT ARRAY ALLOCATION'
     END IF
     CALL EXTCDE ( 1 )
  END IF
  !             Loop through cells
  DO ISEA=1,NSEA
     DADX(1,ISEA) = 0.
     DADY(1,ISEA) = 0.
     IF ( QTREE(IQGW)%CELL_TYPE(ISEA).EQ.                         &
                QTREE(IQGW)%UNDEF_TYPE ) CYCLE
     ICASE = QTREE(IQGW)%NCASE(ISEA)
     SCFAC = 2.**(QTREE(IQGW)%INDLVL(ISEA) - QTREE(IQGW)%LVLREF)
     CALL QA_DERIVC( ISEA, ICASE, NSOL, ISOL, WTS1_QA%NNGBR,      &
                     QTREE(IQGW)%NGBR(ISEA,:), A, WTS1_QA%DERWTC, &
                     DERVAL, ierr=IERR, ndse=NDSE )
     IF ( IERR.GT.0 ) THEN
        IF ( NDSE.GT.0 ) THEN
           WRITE (NDSE,*) ' ERROR IN W3QTGR'
           WRITE (NDSE,*) ' WHILE CALLING QA_DERIVC'
        END IF
        CALL EXTCDE ( 3 )
     END IF
     DADX(1,ISEA) = DERVAL(1)*SCFAC
     DADY(1,ISEA) = DERVAL(2)*SCFAC
  END DO
  RETURN
  !/
  !/ End of W3QTGR ----------------------------------------------------- /
  !/
  END SUBROUTINE W3QTGR

  !/
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3FLQO ( J, IUN, NCELL, NQUAD, IERR )
  !/
  !/    Richard Gorman, NIWA, May 2014
  !/
  !  1. Purpose :
  !
  !     Allocate quadtree structures for different inputs (levels, currents,
  !     wind, ice) and read their initial states
  !
  !  2. Method :
  !
  !
  !  3. Parameters :
  !
  !     Parameter list
  !     ----------------------------------------------------------------
  !       J          Int.  I   Code for variable selected:
  !                               1,2,3,4 for level, current, wind, ice
  !       IUN        Int.  I   Unit number for input file
  !       NCELL      Int.  I   Number of cells for allocation
  !       NQUAD      Int.  I   Number of quads for allocation
  !       IERR       Int.  O   Return status
  !     ----------------------------------------------------------------
  !
  !     Local data
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3QALL    Subr. w3gdatmd  Allocation of quadtree arrays
  !      QA_IOQT   Subr. QA_UTILS  IO of quadtree arrays
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      WSHELL    Prog. WW3_SHEL  Wave model shell
  !     ----------------------------------------------------------------
  !
  !  6. Error messages :
  !
  !
  !  7. Remarks :
  !
  !
  !  8. Structure :
  !
  !     -------------------------------------------
  !     -------------------------------------------
  !
  !  9. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !       !/T     Test output.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  USE QA_UTILS, ONLY  :  QA_IOQT, QA_CPQT
  USE W3GDATMD, ONLY  :  W3QALL
  USE W3GDATMD, ONLY  :  IQGA0, IQGC0, IQGL0, IQGI0,              &
                         IQGAN, IQGCN, IQGLN, IQGIN,              &
                         QTREE, NAUX_QA,NCMXQ, NQMXQ
#ifdef W3_S
  USE W3SERVMD, ONLY: STRACE
#endif
  !
  IMPLICIT NONE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  INTEGER, INTENT(IN)     :: J, IUN, NCELL, NQUAD
  INTEGER, INTENT(OUT)    :: IERR
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  !/
  !/ ------------------------------------------------------------------- /
  !/
#ifdef W3_S  
  CALL STRACE (IENT, 'W3FLQO')
#endif  
  !
  ! 
  !
  IF ( J.EQ.1 ) THEN
     ! sea levels
     NAUX_QA = NAUX_QA + 1
     IQGL0 = NAUX_QA
     NCMXQ(IQGL0) = NCELL
     NQMXQ(IQGL0) = NQUAD
     IQGLN = IQGL0
     ! allocate
     CALL W3QALL ( IQGL0, NCELL, NQUAD, 0, 0 )
     ! read the (initial) quadtree grid data
     CALL QA_IOQT( IUN, QTREE(IQGL0), -1, IERR, ndse=0 )
     IF ( IERR .NE. 0 ) RETURN
  ELSE IF ( J.EQ.2 ) THEN
     ! currents
     NAUX_QA = NAUX_QA + 1
     IQGC0 = NAUX_QA
     NCMXQ(IQGC0) = NCELL
     NQMXQ(IQGC0) = NQUAD
     IQGCN = IQGC0
     ! allocate
     CALL W3QALL ( IQGC0, NCELL, NQUAD, 0, 0 )
     ! read the (initial) quadtree grid data
     CALL QA_IOQT( IUN, QTREE(IQGC0), -1, IERR, ndse=0 )
     IF ( IERR .NE. 0 ) RETURN
     ! If currents are on a time-varying adaptive grid, set up a 
     ! second quadtree structure to allow for time interpolation 
     ! between the different grids
     IF ( QTREE(IQGC0)%DYNAMIC ) THEN
        NAUX_QA = NAUX_QA + 1
        IQGCN = NAUX_QA
        CALL W3QALL ( IQGCN, NCELL, NQUAD, 0, 0 )
        CALL QA_CPQT( QTREE(IQGC0), QTREE(IQGCN) )
        !QTREE(IQGCN)%DYNAMIC = .TRUE.
     END IF
  ELSE IF ( J.EQ.3 ) THEN
     ! winds
     NAUX_QA = NAUX_QA + 1
     IQGA0 = NAUX_QA
     NCMXQ(IQGA0) = NCELL
     NQMXQ(IQGA0) = NQUAD
     IQGAN = IQGA0
     ! allocate
     CALL W3QALL ( IQGA0, NCELL, NQUAD, 0, 0 )
     CALL QA_IOQT( IUN, QTREE(IQGA0), -1, IERR, ndse=0 )
     ! If winds are on a time-varying adaptive grid, set up a 
     ! second quadtree structure to allow for time interpolation 
     ! between the different grids
     IF ( QTREE(IQGA0)%DYNAMIC ) THEN
        NAUX_QA = NAUX_QA + 1
        IQGAN = NAUX_QA
        CALL W3QALL ( IQGAN, NCELL, NQUAD, 0, 0 )
        CALL QA_CPQT( QTREE(IQGA0), QTREE(IQGAN) )
        !QTREE(IQGAN)%DYNAMIC = .TRUE.
     END IF
#ifdef W3_T     
     OPEN(91,FILE='Wind0_quad.txt')
     CALL QA_IOQT(91,QTREE(IQGA0),3,ierr=IERR)
     CLOSE(91)         
     OPEN(91,FILE='Wind0_cell.txt')
     CALL QA_IOQT(91,QTREE(IQGA0),4,ierr=IERR)
     CLOSE(91)         
     OPEN(91,FILE='WindN_quad.txt')
     CALL QA_IOQT(91,QTREE(IQGAN),3,ierr=IERR)
     CLOSE(91)         
     OPEN(91,FILE='WindN_cell.txt')
     CALL QA_IOQT(91,QTREE(IQGAN),4,ierr=IERR)
     CLOSE(91)         
#endif  
  ELSE IF ( J.EQ.4 ) THEN
     NAUX_QA = NAUX_QA + 1
     IQGI0 = NAUX_QA
     NCMXQ(IQGI0) = NCELL
     NQMXQ(IQGI0) = NQUAD
     IQGIN = IQGI0
     ! allocate
     CALL W3QALL ( IQGI0, NCELL, NQUAD, 0, 0 )
     CALL QA_IOQT( IUN, QTREE(IQGI0), -1, IERR, ndse=0 )
     IF ( IERR .NE. 0 ) RETURN
  END IF
  RETURN
  !/
  !/ End of W3FLQO ----------------------------------------------------- /
  !/
  END SUBROUTINE W3FLQO
  !/
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3STQT ( NSTGT, NDSEN, IERR, NDST )
  !/
  !/    Richard Gorman, NIWA, May 2014
  !/
  !  1. Purpose :
  !
  !     Set up an initial wave quadtree
  !
  !  2. Method :
  !
  !
  !  3. Parameters :
  !
  !     Parameter list
  !     ----------------------------------------------------------------
  !       NSTGT      Int.  I   Target number of cells
  !       NDSEN      Int.  I   Unit number for error output
  !       IERR       Int.  O   Return status
  !       NDST       Int.  I   Unit number for test output
  !     ----------------------------------------------------------------
  !
  !     Local data
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3QALL    Subr. w3gdatmd  Allocation of quadtree arrays
  !      QA_IOQT   Subr. QA_UTILS  IO of quadtree arrays
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3STRT    Prog. WW3_STRT  Wave model initialization
  !     ----------------------------------------------------------------
  !
  !  6. Error messages :
  !
  !
  !  7. Remarks :
  !
  !
  !  8. Structure :
  !
  !     -------------------------------------------
  !     -------------------------------------------
  !
  !  9. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !       !/T     Test output.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  USE QA_UTILS, ONLY: QA_DERWTS, QA_Q2NOKEEP, QA_FINDNBR,         &
                      QA_NBDIST, QA_Q2QMAP, QA_IOQT
  USE W3GDATMD, ONLY  :  W3DIMQ, W3QALL
  USE W3GDATMD, ONLY  :  NSEA, NQUAD, IQGW, IQGB, QTREE, WTS1_QA, &
                         NAUX_QA, IAUX_QA, EXAUX_QA, DVTOLFAC_QA, &
                         NCTARGET_QA, DVMAX_QA, MAPML_QA
#ifdef W3_S
  USE W3SERVMD, ONLY: STRACE
#endif
  !
  IMPLICIT NONE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  INTEGER, INTENT(IN)     :: NSTGT
  INTEGER, INTENT(IN)     :: NDSEN, NDST
  INTEGER, INTENT(OUT)    :: IERR
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  INTEGER  :: ISEA, NINDML, IS, ITER, NITER, NCT_FINAL, NCT_INIT
  LOGICAL  :: ADAPT_VARS(2) = .TRUE.
  !/
  !/ ------------------------------------------------------------------- /
  !/
#ifdef W3_S  
  CALL STRACE (IENT, 'W3STQT')
#endif
  !
  !  Precompute weights for derivative operators using only first-order neighbours:
  ! 
  CALL QA_DERWTS( WTS1_QA )
  !        Allocate quadtree arrays. Note, calling with the same number of cells
  !        means W3DIMQ will not deallocate grid arrays from the model definition
  !        file that would then need to be reallocated for the wave quadtree:          
  NAUX_QA = 2
  NINDML = MAXVAL(QTREE(IQGB)%INDML)
  CALL W3DIMQ  ( 1, NSEA, NSEA, NAUX_QA, NINDML, NDSEN, NDST )
  !        Allocate a wave quadtree with the same size as the bathy quadtree
  CALL W3QALL ( IQGW, NSEA, NQUAD, 0, 1 )
  !                            iwtorder=1, qspare=QSPARE,                &
  !                            mpart=NAPROC, mcloc=NSEAL )
  !        Derive an initial wave quadtree from the bathy. quadtree, but
  !        not retaining refined cell indices
  CALL  QA_Q2NOKEEP ( QTREE(IQGB), NSEA, QTREE(IQGW), ierr=IERR,  &
                          ndse=NDSEN )
  IF ( IERR.NE.0 ) THEN
     IF ( NDSEN.GT.0 ) THEN
        WRITE(NDSEN,*) 'ERROR IN W3STQT CALLING QA_2NOKEEP'
        WRITE(NDSEN,*) '    RETURNED IERR = ', IERR
     END IF
     RETURN
  END IF
  QTREE(IQGW)%IWTORDER = 1
  QTREE(IQGW)%DYNAMIC = .TRUE.
  !        Recompute neighbours and 1st order weights
  CALL  QA_FINDNBR ( QTREE(IQGW), ierr=IERR, ndse=NDSEN )
  IF ( IERR.NE.0 ) THEN
     IF ( NDSEN.GT.0 ) THEN
        WRITE(NDSEN,*) 'ERROR IN W3STQT CALLING QA_FINDNBR'
        WRITE(NDSEN,*) '    RETURNED IERR = ', IERR
     END IF
     RETURN
  END IF
  !
  CALL  QA_NBDIST ( QTREE(IQGW), ierr=IERR, ndse=NDSEN )
  IF ( IERR.NE.0 ) THEN
     IF ( NDSEN.GT.0 ) THEN
        WRITE(NDSEN,*) 'ERROR IN W3STQT CALLING QA_NBDIST'
        WRITE(NDSEN,*) '    RETURNED IERR = ', IERR
     END IF
     RETURN
  END IF
  !
  !     Compute status map on multilevel grid
  MAPML_QA = QTREE(IQGB)%UNDEF_TYPE
  DO ISEA=1,QTREE(IQGB)%NCELL
     IS = QTREE(IQGB)%INDML(ISEA)
     MAPML_QA(IS) = QTREE(IQGB)%CELL_TYPE(ISEA)
  END DO
  !
  !        Trivial mapping between wave quadtree and itself:
  DO ISEA=1,NSEA
     IAUX_QA(ISEA,IQGW) = ISEA
  END DO
  EXAUX_QA(:,IQGW) = .TRUE.      
  !
  !        Compute the mapping between the wave and bathy quadtrees          
  CALL QA_Q2QMAP( QTREE(IQGW), QTREE(IQGB), .TRUE.,               &
                  IAUX_QA(:,IQGB), EXAUX_QA(:,IQGB), ierr=IERR,   &
                  ndse=NDSEN )
  IF ( IERR.NE.0 ) THEN
     IF ( NDSEN.GT.0 ) THEN
        WRITE(NDSEN,*) 'ERROR IN W3STQT CALLING QA_Q2QMAP'
        WRITE(NDSEN,*) '    RETURNED IERR = ', IERR
     END IF
     RETURN
  END IF
  !
  !        Re-compute spatial arrays on the new wave quadtree:
  DO ISEA=1,NSEA
     CALL W3ADG1 ( ISEA )
  END DO
  !
  !      Adapt the wave quadtree down to the required limits
  !      Only want to coarsen, so ensure upper tolerance is high:
  DVTOLFAC_QA = 1.E9
  DVMAX_QA = 0.
  NCT_INIT = QTREE(IQGW)%NCELL_DEF
  NCT_FINAL = NCT_INIT
  IF ( NSTGT.GT.0 ) NCT_FINAL = MIN(NCT_FINAL, NSTGT) 
  IF ( NCT_FINAL.LT.NCT_INIT ) THEN
     NITER = QTREE(IQGW)%LVLHI
     DO ITER = 1, NITER
        NCTARGET_QA = NCT_FINAL + (NITER-ITER)*(NCT_INIT-NCT_FINAL)/NITER
        CALL W3ADGR ( 1, ADAPT_VARS, NDSEN )
        IF ( QTREE(IQGW)%NCELL_DEF.LE.NCT_FINAL ) EXIT
     END DO
  END IF
  NCTARGET_QA = NCT_FINAL
  NSEA = QTREE(IQGW)%NCELL
  NQUAD = QTREE(IQGW)%NQUAD
  !      
  RETURN
  !/
  !/ End of W3STQT ----------------------------------------------------- /
  !/
  END SUBROUTINE W3STQT
  !/
  !/ ------------------------------------------------------------------- /
  SUBROUTINE W3BCRESET ( NDST, NDSE )
  !/
  !/    Richard Gorman, July 2026
  !/
  !  1. Purpose :
  !
  !     Recompute interpolation points and indices for boundary outputs on
  !     a quadtree grid
  !
  !  2. Method :
  !
  !
  !  3. Parameters :
  !
  !     Parameter list
  !     ----------------------------------------------------------------
  !       NDST       Int.  I   Unit number for test output
  !       NDSE       Int.  I   Unit number for error output
  !     ----------------------------------------------------------------
  !
  !     Local data
  !     ----------------------------------------------------------------
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines used :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      QA_INTERP1WT Subr. QA_UTILS  Compute interpolation weights
  !     ----------------------------------------------------------------
  !
  !  5. Called by :
  !
  !      Name      Type  Module   Description
  !     ----------------------------------------------------------------
  !      W3WAVE    Prog. W3WAVE   Wave model
  !     ----------------------------------------------------------------
  !
  !  6. Error messages :
  !
  !
  !  7. Remarks :
  !
  !
  !  8. Structure :
  !
  !     -------------------------------------------
  !     -------------------------------------------
  !
  !  9. Switches :
  !
  !       !/S     Enable subroutine tracing.
  !       !/T     Test output.
  !
  ! 10. Source code :
  !
  !/ ------------------------------------------------------------------- /
  USE QA_UTILS, ONLY:  QA_INTERP1WT
  USE W3GDATMD, ONLY:  IQGW, QTREE, WTS1_QA, SX, X0, SY, Y0, MAPSTA, MAPFS
  USE W3ODATMD, ONLY:  NWTBO, NFBPO, NBO, NBO2, XBPO, YBPO, RDBPO,  &
                       IPBPO, ISBPO
#ifdef W3_MPI
  USE W3GDATMD, ONLY:  NSPEC
  USE W3WDATMD, ONLY:  VA
  USE W3ADATMD, ONLY:  MPI_COMM_WAVE
  USE W3ODATMD, ONLY:  IT0BPT, NRQBP, NRQBP2, IRQBP1, IRQBP2, NAPBPT, IAPROC, ABPOS
  USE W3PARALL, ONLY:  INIT_GET_JSEA_ISPROC
  USE MPI_F08
#endif
#ifdef W3_S
  USE W3SERVMD, ONLY: STRACE
#endif
  !
  IMPLICIT NONE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  INTEGER, INTENT(IN)     :: NDST
  INTEGER, INTENT(IN)     :: NDSE
  !/
  !/ ------------------------------------------------------------------- /
  !/ Local parameters
  !/
  INTEGER  :: IFIL, IPT, J, IST, IERR
  INTEGER, ALLOCATABLE :: IXR(:), IYR(:), ISEAI(:)
  REAL     :: XO, YO, RDTOT
  REAL, ALLOCATABLE :: RD(:)
  LOGICAL  :: FLNEW
#ifdef W3_MPI
  INTEGER  :: ISPROC, IH, IT, IT0, ISEA, JSEA, ITARG, IROOT
#endif
  !/
  !/ ------------------------------------------------------------------- /
  !/
#ifdef W3_S  
  CALL STRACE (IENT, 'W3BCRESET')
#endif
  !
  ALLOCATE ( IXR(NWTBO), IYR(NWTBO), ISEAI(NWTBO), RD(NWTBO) )
#ifdef W3_T
  WRITE (NDST,9090)
#endif
  NBO2 = 0
  DO IFIL=1,NFBPO
    DO IPT=NBO(IFIL-1)+1,NBO(IFIL)
      XO = XBPO(IPT)
      YO = YBPO(IPT)
      !
      ! Assign weights for interpolation on the wave quadtree: 
      ! This uses up to 9 cells/weights
      !
      CALL QA_INTERP1WT( QTREE(IQGW), 1.+(XO-X0)/SX, 1.+(YO-Y0)/SY, &
                               WTS1_QA%DERWTC, IXR, RD, IERR, NDSE )
      IYR = 1
#ifdef W3_T
      WRITE (NDST,9097) XO, YO,                       &
               (IXR(J), IYR(J), RD(J), J=1,NWTBO)
#endif
      !
      ! ... Interpolation factors
      !
      RDTOT = 0.
      DO J=1, NWTBO
        RDBPO(IPT,J) = 0.
        ISEAI(J) = 0
        IF ( IYR(J).GT.0 .AND. IXR(J).GT.0 ) THEN
          IF ( MAPSTA(IYR(J),IXR(J)).GT.0 .AND.               &
                       RD(J).GT.0.05 ) THEN
            RDBPO(IPT,J) = RD(J)
          END IF
          ISEAI(J) = MAPFS(IYR(J),IXR(J))
        END IF
        RDTOT = RDTOT + RDBPO(IPT,J)
      END DO
      !
      DO J=1, NWTBO
        RDBPO(IPT,J) = RDBPO(IPT,J) / RDTOT
      END DO
              !
#ifdef W3_T
      WRITE (NDST,9098) RDTOT, (RDBPO(IPT,J),J=1,NWTBO)
#endif
      !
      ! ... Determine sea and interpolation point counters
      !
      DO J=1, NWTBO
        IF ( ISEAI(J).EQ.0 .OR. RDBPO(IPT,J).EQ. 0. ) THEN
          IPBPO(IPT,J) = 0
        ELSE
          FLNEW   = .TRUE.
          DO IST=NBO2(IFIL-1)+1, NBO2(IFIL)
            IF ( ISEAI(J) .EQ. ISBPO(IST) ) THEN
              FLNEW  = .FALSE.
              IPBPO(IPT,J) = IST - NBO2(IFIL-1)
            END IF
          END DO
          IF ( FLNEW ) THEN
            NBO2(IFIL)        = NBO2(IFIL) + 1
            IPBPO(IPT,J)      = NBO2(IFIL) - NBO2(IFIL-1)
            ISBPO(NBO2(IFIL)) = ISEAI(J)
          END IF
        END IF
      END DO
      !
#ifdef W3_T
      WRITE (NDST,9099) ISEAI, (IPBPO(IPT,J),J=1,NWTBO)
#endif
    END DO
  END DO
  !
#ifdef W3_MPI  
  !
  ! 3.  Set-up for W3IOBC ( SENDs ) ------------------------------------ /
  !
  NRQBP  = 0
  NRQBP2 = 0
  IH     = 0
  IT     = IT0BPT
  IROOT  = NAPBPT - 1
  !
#ifdef W3_MPIT
  WRITE (NDST,9030) 'MPI_SEND_INIT'
#endif
  !
  DO J=1, NFBPO
    DO IPT=NBO2(J-1)+1, NBO2(J)
      !
      IT     = IT + 1
      !
      ! 3.b Residence processor of point
      !
      ISEA   = ISBPO(IPT)
      CALL INIT_GET_JSEA_ISPROC(ISEA, JSEA, ISPROC)
      !
      ! 3.c If stored locally, send data
      !
      IF ( IAPROC .EQ. ISPROC ) THEN
        IH     = IH + 1
        CALL MPI_SEND_INIT (VA(1,JSEA),NSPEC,MPI_REAL, IROOT, IT, MPI_COMM_WAVE, &
             IRQBP1(IH), IERR)
#ifdef W3_MPIT
        WRITE (NDST,9031) IH, IPT, J, IROOT, IT, IRQBP1(IH), IERR
#endif
      END IF
      !
    END DO
  END DO
  !
  ! ... End of loops 4.a
  !
  NRQBP  = IH
  !
#ifdef W3_MPIT
  WRITE (NDST,9032)
  WRITE (NDST,9033) NRQBP
#endif
  !
  ! 3.d Set-up for W3IOBC ( RECVs ) ------------------------------------ /
  !
  IF ( IAPROC .EQ. NAPBPT ) THEN
    !
    IH     = 0
    IT     = IT0BPT
    !
    ! 3.e Loops over files and points
    !
#ifdef W3_MPIT
    WRITE (NDST,9030) 'MPI_RECV_INIT'
#endif
    !
    DO J=1, NFBPO
      DO IPT=NBO2(J-1)+1, NBO2(J)
        !
        ! 3.f Residence processor of point
        !
        ISEA   = ISBPO(IPT)
        CALL INIT_GET_JSEA_ISPROC(ISEA, JSEA, ISPROC)
        !
        ! 3.g Receive in correct array
        !
        IH     = IH + 1
        IT     = IT + 1
        ITARG  = ISPROC - 1
        CALL MPI_RECV_INIT (ABPOS(1,IH),NSPEC,MPI_REAL, ITARG, IT, MPI_COMM_WAVE, &
             IRQBP2(IH), IERR)
#ifdef W3_MPIT
        WRITE (NDST,9031) IH, IPT, J, ITARG, IT, IRQBP2(IH), IERR
#endif
        !
      END DO
    END DO
    !
    NRQBP2 = IH
    !
    ! ... End of loops 4.e
    !
#ifdef W3_MPIT
    WRITE (NDST,9032)
    WRITE (NDST,9033) NRQBP2
#endif
    !
  END IF
  !
#endif    
  !
#ifdef W3_T
9090 FORMAT ( ' TEST W3BCRESET : OUTPUT BOUND. POINT DATA LINE SEG.')
9097 FORMAT ( '             ',2F8.2,9(2I4,F7.2))
9098 FORMAT ( '                      ',F7.2,2X,9F7.2)
9099 FORMAT ( '                            ',9I7/                   &
         '                            ',9I7)
#endif
#ifdef W3_MPIT
9030 FORMAT (/' TEST W3BCRESET : ',A,' CALLS FOR W3IOBC'/           &
         ' +------+------+---+------+------+--------------+'/       &
         ' |  IH  | IPT  | F | TARG |  TAG |   handle err |'/       &
         ' +------+------+---+------+------+--------------+')
9031 FORMAT ( ' |',2(I5,' |'),I2,' |',2(I5,' |'),I9,I4,' |')
9032 FORMAT (                                                       &
         ' +------+------+---+------+------+--------------+')
9033 FORMAT ( ' TEST W3BCRESET : NRQBC :',I10)
#endif
  !/
  !/ End of W3BCRESET -------------------------------------------------- /
  !/
  END SUBROUTINE W3BCRESET
  !/
  !/ End of module W3ADGRMD -------------------------------------------- /
  !/
  END MODULE W3ADGRMD

