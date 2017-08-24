      subroutine roraf(iyr, imon, idy, flonum, n, da, iyrst, iyren,
     $     xnewmin, idiff, k, te, ta, qp, qa, qb, qc, c,
     $     delq, rech, year, mon, day, npeaks, itbase, ierr)
C          BY AL RUTLEDGE, USGS           2002 VERSION
C          MODIFIED To run in R BY DAVE LORENZ September, 2015
C  THIS PROGRAM ESTIMATES GROUND-WATER RECHARGE FROM A DAILY-VALUES
C  RECORD OF STREAMFLOW, USING THE "RORABAUGH METHOD."  THIS AND
C  OTHER PROGRAMS ARE DOCUMENTED IN USGS WRIR 98-4148. THIS 2006
C  VERSION OF THE PROGRAM runs as a subroutine callable from R
C
C  SEVERAL DECLARATION STATEMENTS PERTAIN TO ARRAY SIZES:
C          MAXIMUM NUMBER OF YEARS = 120
C          MAXIMUM NUMBER OF DAYS = 44000.
C               (ALSO NOTE INITIALIZING STATEMENTS)
C          RETRIEVE DAILY VALUES AFTER THE YEAR 1890.
C               (ALSO NOTE LINES BETWEEN LABELS 330 AND 340)
C
C ---- Input arguments
      integer iyr(*), imon(*), idy(*)
      real*8 flonum(*)
      integer n
      real*8 da
      integer iyrst, iyren
      real*8 xnewmin
      integer idiff
      real*8 k
C ---- Ouput arguments
      INTEGER TE(*)
      real*8 TA(*), QP(*), QA(*), QB(*), QC(*), C(*), 
     $     DELQ(*), RECH(*)
      integer year(*), mon(*), day(*), npeaks, itbase, ierr

      COMMON/BIG/Q(120,12,31)
      COMMON/BIG/EST(120,12,31)
      COMMON/BIG/Q1D(44000)
      COMMON/BIG/B1D(44000,3)
      COMMON/BIG/IYR1D(44000)
      COMMON/BIG/IMO1D(44000)
      COMMON/BIG/IDA1D(44000)
      COMMON/BIG/EST1D(44000)
      COMMON/BIG/ALLGW(44000)
      COMMON/BIG/TP(6000)
      COMMON/BIG/TS(6000)
      COMMON/BIG/TBC(6000)
      real*8 Q, Q1D, TP, TBC, SUM
      INTEGER IYR1D, IMO1D, IDA1D, TS
      INTEGER iyearst, iyearen, IBEFORE
      CHARACTER*1 EST
      CHARACTER*1 EST1D
      CHARACTER*1 ALLGW
 
      IBEFORE=iyr(1) - 1
 
C
C------------------------- INITIALIZE VARIABLES : -------------------
C
      IERR = 0
      DO 10 IYEAR=1,120
      DO 10 IMONTH=1,12
      DO 10 IDAY=1,31
   10 Q(IYEAR,IMONTH,IDAY)=-999.0D0
      DO 11 I=1,44000
      ALLGW(I)= ' '
   11 Q1D(I)= -999.0D0
C
C
C -------- PROCESS THE DAILY-VALUES FILE OF STREAMFLOW: ------------------
C
      IFRSTYR=iyr(1)
      do 35 i=1,n
         IYEAR = IYR(i)-IBEFORE
         Q(IYEAR,imon(i),idy(i)) = FLONUM(i)
 35   CONTINUE

      ILSTYR= IYR(n)
c
c --- flag nonexistent dates with flow=-9999 -----
c
      DO 38 IYEAR=IFRSTYR, ILSTYR
         DO 37 IMONTH=1,12
            DO 36 IDAY=1,31
               IF((IMONTH.EQ.2).AND.(IDAY.GT.29)) THEN
                  Q(IYEAR-IBEFORE,IMONTH,IDAY)= -9999.0D0
               END IF
               IF((IMONTH.EQ.2).AND.(IDAY.EQ.29)) THEN
                  IDIV=INT((IYEAR)/4.0)
                  XDIV=(IYEAR)/4.0
                  DIFFER=ABS(IDIV-XDIV)
                  IF(DIFFER.GT.0.1) THEN
                     Q(IYEAR-IBEFORE,IMONTH,IDAY)= -9999.0D0
                  END IF
               END IF
               IF(IDAY.EQ.31) THEN
                  IF((IMONTH.EQ.4).OR.(IMONTH.EQ.6).OR.(IMONTH.EQ.9)
     $                 .OR.(IMONTH.EQ.11)) THEN
                     Q(IYEAR-IBEFORE,IMONTH,IDAY)= -9999.0D0
                  END IF
               END IF
 36         CONTINUE
 37      CONTINUE
 38   CONTINUE
      IFRSTYR= IFRSTYR-IBEFORE
      ILSTYR= ILSTYR-IBEFORE
C
C   ---------------   SELECT TIME PERIOD OF INTEREST:  ------------------
C
      IIMAX=0
      DO 120 IIYR=IYRST,IYREN
           IDIV=INT(IIYR/4.0D0)
           XDIV=IIYR/4.0D0
           DIFFER=ABS(IDIV-XDIV)
           IF(DIFFER.LT.0.1D0) THEN
               IIMAX=IIMAX+366
             ELSE
               IIMAX=IIMAX+365
            END IF
  120 CONTINUE
      IYEARST= IYRST-IBEFORE
      IYEAREN= IYREN-IBEFORE
C
C ---- ASSIGN VALUES TO 1-DIMENSIONAL ARRAYS OF DISCHARGE AND DATE: ----
C
      ICOUNT= 0
      IBREAK= 0
      ITESTX= 0
      DO 180 IYEAR= IYEARST, IYEAREN
         DO 179 IMONTH= 1,12
            DO 178 IDAY= 1, 31
               SFLOW= Q(IYEAR,IMONTH,IDAY)
               IF(SFLOW .EQ. -9999.D0) GO TO 178
               IF(SFLOW .EQ. -99D0 .OR. SFLOW .EQ. -999D0) THEN
                  ITEST=0
               ELSE
                  ITEST=1
               ENDIF
               IF(ITEST.EQ.1.AND.ITESTX.EQ.0) THEN
                  IBREAK= IBREAK+1
                  IF(IBREAK.GT.1) THEN
                     IERR = 2
                     return
                  ENDIF
               ENDIF
               ITESTX= ITEST
               IF(SFLOW .EQ. -99D0 .OR. SFLOW .EQ .-999D0) GO TO 178
               
               ICOUNT= ICOUNT + 1
               Q1D(ICOUNT)= SFLOW
               IYR1D(ICOUNT)= IYEAR
               IMO1D(ICOUNT)= IMONTH
               IDA1D(ICOUNT)= IDAY
               EST1D(ICOUNT)= EST(IYEAR,IMONTH,IDAY)
 178        CONTINUE
 179     CONTINUE
 180  CONTINUE
C
      IZERO=1
      DO 185 I=1,ICOUNT
      IF(Q1D(I).EQ.0.0) IZERO=0
  185 ALLGW(I)= ' '
C
C ----------- Modify 0 flows
C
      IF (IZERO.EQ.0) THEN
         DO 187 I=1,ICOUNT
            IF(Q1D(I).EQ.0.00D0) THEN
               Q1D(I)= XNEWMIN
            END IF
 187     CONTINUE
      END IF

C
C ------------ Process DA
C
      IF (DA.LT.1.0) IERR = -2
      IF (DA.GT.500.0D0) IERR = -3

C
C  ---  DETERMINE THE MINIMUM NUMBER OF DAYS OF ANTECEDENT RECESSION  --
C  ---  TO INDICATE THAT STREAMFLOW IS TO BE CONSIDERED GROUND-WATER  --
C  ---  DISCHARGE.  OBTAIN FROM THE EQUATION DA**0.2, ROUNDED UPWARD  --
C
      ITBASE=0
  210 ITBASE= ITBASE + 1
      IF(ITBASE .GT. 10) THEN
         IERR = 3
         return
      END IF
      IF (ITBASE.GT.DA**0.2D0) THEN
         GO TO 220
      ELSE
         GO TO 210
      END IF
 220  CONTINUE

C
C ----- SPECIFY THE MAXIMUM ALLOWABLE NUMBER OF DAYS THAT CAN BE USED
C ----- AFTER A PEAK, TO DETERMINE THE GROUND-WATER DISCHARGE AFTER
C ----- THE PEAK: 
C
      IRECMAX= INT(0.2144D0*K)

C ---- ALLOW THE USER TO OVERRIDE THE DEFAULT METHOD OF EXECUTION.  -----
C ---- NOTE: THE ALLOWANCE FOR NON-DEFAULT EXECUTION, IS DIFFERENT  -----
C ---- FROM THE NON-DEFAULT EXECUTION IN PREVIOUS VERSIONS OF RORA ------

      IF (IDIFF .LT. 0 .OR. IDIFF .GT. 3) THEN
         IERR = 1
         RETURN
      ELSE
         ITBASE= ITBASE+IDIFF
      ENDIF
          
      IF(IRECMAX.LT.ITBASE) THEN
         IERR = -1
         IRECMAX= ITBASE
      END IF
C
C-----------------ON DATES PRECEEDED BY A RECESSION PERIOD,-------------
C-------------------------SET VARIABLE ALLGW='*'  ----------------------
C
      DO 270 I= ITBASE+1, ICOUNT
         INDICAT=1
         IBACK= 0
  260    IBACK= IBACK + 1
            IF(Q1D(I-IBACK).LT.Q1D(I-IBACK+1)) INDICAT=0
            IF(IBACK.LT.ITBASE) GO TO 260
            IF(INDICAT.EQ.1) THEN
                   ALLGW(I)= '*'
               ELSE
                   ALLGW(I)= ' '
             END IF
  270 CONTINUE
C
C
C -- THE FOLLOWING LINES UP TO LABEL 800 DETERMINE THE LOCATION (IN TIME) --
C --------------------  OF PEAKS AND RECESSION PERIODS:  -------------------
C
C
C ------------------- LOCATE END OF FIRST RECESSION: ----------------
C
      I= 0
 280  I= I+1
      IF(I.GT.ICOUNT) THEN
         IERR = 4
         return
      END IF
      IF (ALLGW(I).NE.'*') GO TO 280
 290  I=I+1
      IF (I.GT.ICOUNT) THEN
         IERR = 5
         return
      END IF
      IF (ALLGW(I).NE.' ') GO TO 290
      I=I-1
      TA(1)= I
      QA(1)= Q1D(I)
      IPEAK= 0
C
C -------------    FIND STREAMFLOW AND DAY OF THE PEAK:     -----------
C
 330  CONTINUE
      IPEAK= IPEAK+1
      IF(IPEAK .GT. 6000) THEN
         IERR = 6
         return
      END IF
      QP(IPEAK)= Q1D(I)
      TP(IPEAK)= I
      ILOOK= I + 1
 340  CONTINUE
      IF(ILOOK.GE.ICOUNT) GO TO 800
      IF(ALLGW(ILOOK).EQ.'*') GO TO 350
      IF (Q1D(ILOOK).GE.QP(IPEAK)) THEN
         QP(IPEAK)= Q1D(ILOOK)
         TP(IPEAK)= ILOOK
      END IF
      IF(ALLGW(ILOOK).EQ.' ') THEN
         ILOOK= ILOOK+1
         GO TO 340
      END IF
 350  CONTINUE
      TBC(IPEAK)= TP(IPEAK) + (0.2144*K)
      
C     
C     -----   FIND FIRST AND LAST DAYS OF RECESSION FOLLOWING THE PEAK: ----
C     
      I= ILOOK
      TS(IPEAK)= I
 370  CONTINUE
      I= I+1
      IF(I.GT.ICOUNT) GO TO 380
      IF(ALLGW(I).NE.' ') THEN
         GO TO 370
      END IF
  380 CONTINUE
      I= I-1
      TE(IPEAK)= I
      IF (TE(IPEAK)-TP(IPEAK).GT.IRECMAX) THEN
         TE(IPEAK)= TP(IPEAK) + IRECMAX
      END IF
      IF (TE(IPEAK).LT.TS(IPEAK)) THEN
         TE(IPEAK)= TS(IPEAK)
      END IF
      
      I=I+1
      IF(I.LT.ICOUNT) GO TO 330

C
C ----  GO HERE WHEN ALL RECESSION PERIODS AND PEAKS HAVE BEEN LOCATED: ----
C
 800  CONTINUE

      NPEAKS= IPEAK - 1
      
C
C --- EXTRAPOLATE STREAMFLOW, DETERMINE RECESSION CURVE DISPLACEMENT, AND --
C ----------------- CALCULATE RECHARGE FOR FIRST PEAK:  --------------------
C
      QB(1)= QA(1)*10**(-1*(TBC(1)-TA(1))/K)
      SUM= 0.0D0
      DO 810 I= TS(1), TE(1)
         DQ= Q1D(I) - (QA(1)*10D0**(-1.0D0*(I-TA(1))/K))
 810  SUM= SUM + (DQ*SQRT(I-TP(1)))
      NRECS= TE(1)-TS(1)+1
      C(1)= SUM/NRECS
      DELQ(1)= C(1)/(SQRT(TBC(1)-TP(1)))
      RECH(1)= 2.0D0*0.0371901D0*DELQ(1)*K / (2.3025851D0*DA)     
      QC(1)= QB(1) + DELQ(1)
      year(1) = IBEFORE+IYR1D(INT(TP(1)))
      mon(1) = IMO1D(INT(TP(1)))
      day(1) = IDA1D(INT(TP(1)))
C
C -- EXTRAPOLATE STREAMFLOW, DETERMINE RECESSION CURVE DISPLACEMENT, AND --
C --------------- CALCULATE RECHARGE FOR ALL OTHER PEAKS: -----------------
C
      DO 840 IPEAK=2, NPEAKS
         QA(IPEAK)= QC(IPEAK-1)
         TA(IPEAK)= TBC(IPEAK-1)
         QB(IPEAK)= QA(IPEAK)*10**(-1*(TBC(IPEAK)-TA(IPEAK))/K)
         
         SUM= 0.0D0
         DO 820 I= TS(IPEAK), TE(IPEAK)
            IF(I.GT.TA(IPEAK)) THEN
               BASELIN= QA(IPEAK)*10**(-1*(I-TA(IPEAK))/K)
            ELSE
               BASELIN= C(IPEAK-1)/(SQRT(I-TP(IPEAK-1))) +
     $              QA(IPEAK-1)*10D0**(-1.0D0*(I-TA(IPEAK-1))/K)
            END IF
            DQ= Q1D(I) - BASELIN
 820     SUM= SUM + (DQ*SQRT(I-TP(IPEAK)))
         NRECS= TE(IPEAK) - TS(IPEAK) + 1
         C(IPEAK)= SUM/NRECS
         DELQ(IPEAK)= C(IPEAK) / (SQRT(TBC(IPEAK)-TP(IPEAK)))
         RECH(IPEAK)= 2.0D0*0.0371901D0*DELQ(IPEAK)*K / (2.3025851D0*DA)
         QC(IPEAK)= QB(IPEAK) + DELQ(IPEAK)
         year(ipeak) = IBEFORE+IYR1D(INT(TP(IPEAK)))
         mon(ipeak) = IMO1D(INT(TP(IPEAK)))
         day(ipeak) = IDA1D(INT(TP(IPEAK)))
 840  CONTINUE
C
C ----- Adjust TE and TA to be relative to TP
C
      do 850 i=1, npeaks
         TE(i) = TE(i) - int(tp(i))
         TA(i) = TA(i) - tp(i)
 850  continue
      return
      END
