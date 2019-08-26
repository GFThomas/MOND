Subroutine densitycenter(NDIM,OM,OPOS,NMAX,FRAC,MTHRES,POSC,&
                              RDENS,AVDENS,NSTAT)
!     DENSITY CENTER, DENSITY RADIUS AND DENSITY-WEIGHTED CENTRAL AVERAGE DENSITY
!     AFTER CASERTANO & HUT (1985 ApJ, 298, 80-94)
!
!     INPUT:
!     NDIM        : Number of dimensions. Note: Uses volume scaling only!
!                   To use ndim<3 leave one or two components void (untested!)
!     OM(NMAX)    : Original (complete) mass array
!     OPOS(3,NMAX): Original (complete) positions array
!     NMAX        : Total number of points
!     FRAC        : Fraction of this used (randomly) for calculation
!     MTHRES      : Bodies beyond this mass are skipped (sinks, black holes...)
!
!     OUTPUT:
!     POSC(3)     : Position of density center
!     RDENS       : Density radius
!     AVDENS      : Density-weighted average density
  USE rnd_array_module
  IMPLICIT NONE
  INTEGER :: I,J,K,NDIM,NMAX,NUSED,NNEIBS=6
  INTEGER :: ORD(NMAX),counter,NSTAT
  LOGICAL isMTHRES
  REAL(KIND=8) :: pi=3.141592653589793238D0
  REAL(KIND=8) :: pi43=4.188790204786390985D0
  REAL(KIND=8),DIMENSION(3,NMAX) :: pos
  REAL(KIND=8),DIMENSION(3,NMAX) :: opos
  REAL(KIND=8),DIMENSION(NMAX) :: M,OM,RSQ,RNNEIBS,DENSITIES,XXX
  REAL(KIND=8),DIMENSION(3) :: POSC,zahler
  REAL(KIND=8) :: FRAC,MTHRES,MNEIBS,VOL,nenner,rdens,avdens

  IF (FRAC<=0.d0) THEN
     WRITE(*,*) "DENSITYCENTER: FRAC must be >0 and <=1!"
     STOP
  ENDIF

!  CALL get_rnd_array(0,NMAX,xxx)
  DO i=1,nmax
     CALL wrap_one_rnd(xxx(i))
  ENDDO

  counter=0
  DO I=1,NMAX
     IF ((xxx(i).le.FRAC).and.(OM(i).le.MTHRES)) THEN
        counter=counter+1.d0
        m(counter)=OM(i)
        DO j=1,3
           POS(j,counter)=OPOS(j,i)
        END DO
!        print *,'pos(counter) =',pos(1:3,counter)
     END IF
  END DO
  NUSED=counter


  DO I=1,NUSED
     DO J=1,NUSED
        IF (J.NE.I) THEN
           RSQ(J)=(POS(1,I)-POS(1,J))**2+(POS(2,I)-POS(2,J))**2+    &
                (POS(3,I)-POS(3,J))**2
        END IF
     END DO
     CALL RQSORT (Nused,RSQ,ORD)
     MNEIBS=0.d0
     DO K=1,NNEIBS-1!skip ref star and outermost neighbour
!        MNEIBS=MNEIBS+M(K)!assume equal masses here?
        MNEIBS=MNEIBS+M(ORD(K))!IT 2015-07-09
     END DO
     RNNEIBS(I)=sqrt(RSQ(ORD(NNEIBS)))
     IF (ndim==3) THEN
!        VOL=pi43*sqrt(RSQ(ORD(NNEIBS)))**3
        VOL=pi43*RNNEIBS(I)**3
     ELSE IF (ndim==2) THEN
        VOL=pi*RSQ(ORD(NNEIBS))
     ELSE IF (ndim==1) THEN
!        VOL=sqrt(RSQ(ORD(NNEIBS)))
        VOL=RNNEIBS(I)
     ELSE
        WRITE(*,*) "Bad Dimension (must be 1...3)"
     ENDIF
     DENSITIES(I)=MNEIBS/VOL!=\rho_{nneibs}^{(I)} in Casertano & Hut's notation
     IF (i==nstat*int(i/MAX(nstat,1)).OR.i==nused) THEN
        !            write(*,*) "Did star No.",i,' of',nmax
        write(*,'(F6.2,A)') 100.*float(i)/nused,"% finished"
     ENDIF
  END DO

!  nenner=0.d0!this might be wrong
  DO k=1,3
     zahler(K)=0.d0
     nenner=0.d0!this seems to be better (IT, 2015-07-07)
     DO I=1,NUSED
        zahler(K)=zahler(K)+ POS(K,I)*DENSITIES(I)
        nenner   =nenner   + DENSITIES(I)
     END DO
     POSC(k)=zahler(K)/nenner
  END DO
      
  zahler=0.d0! recycle zahler
  DO I=1,NUSED
     zahler(1)=zahler(1)+sqrt((POS(1,I)-POSC(1))**2+ &
          (POS(2,I)-POSC(2))**2+                     &
          (POS(3,I)-POSC(3))**2)*DENSITIES(I)
     zahler(2)=zahler(2)+DENSITIES(I)**2
  END DO
      
  RDENS =zahler(1)/nenner!re-use nenner here
  AVDENS=zahler(2)/nenner
!     write(*,*) "density center", POSC(1), POSC(2), POSC(3)  
 
END Subroutine densitycenter
