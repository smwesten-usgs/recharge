subroutine roraf(iyr, imon, idy, flonum, n, da, iyrst, iyren,                 &
    xnewmin, idIFf, k, te, ta, qp, qa, qb, qc, c,                         &
    delq, rech, year, mon, day, npeaks, itbase, ierr)

  !          BY AL RUTLEDGE, USGS           2002 VERSION
  !          MODIFIED To run in R BY DAVE LORENZ September, 2015
  !  THIS PROGRAM ESTIMATES GROUND-WATER RECHARGE FROM A DAILY-VALUES
  !  RECORD OF STREAMFLOW, USING THE "RORABAUGH METHOD."  THIS AND
  !  OTHER PROGRAMS ARE doCUMENTED IN USGS WRIR 98-4148. THIS 2006
  !  VERSION OF THE PROGRAM runs as a SUBROUTINE callable from R
  !
  !  SEVERAL DECLARATION STATEMENTS PERTAIN TO ARRAY SIZES:
  !          MAXIMUM NUMBER OF YEARS = 120
  !          MAXIMUM NUMBER OF DAYS = 44000.
  !               (ALSO NOTE INITIALIZING STATEMENTS)
  !          RETRIEVE DAILY VALUES AFTER THE YEAR 1890.
  !               (ALSO NOTE LINES BETWEEN LABELS 330 AND 340)
  !
  ! ---- Input arguments
  integer iyr(*), imon(*), idy(*)
  real*8 flonum(*)
  integer n
  real*8 da
  integer iyrst, iyren
  real*8 xnewmin
  integer idIFf
  real*8 k
  ! ---- Ouput arguments
  integer TE(*)
  real*8 TA(*), QP(*), QA(*), QB(*), QC(*), C(*),                           &
      DELQ(*), RECH(*)
  integer year(*), mon(*), day(*), npeaks, itbase, ierr

  common/BIG/Q(120,12,31)
  common/BIG/EST(120,12,31)
  common/BIG/Q1D(44000)
  common/BIG/B1D(44000,3)
  common/BIG/IYR1D(44000)
  common/BIG/IMO1D(44000)
  common/BIG/IDA1D(44000)
  common/BIG/EST1D(44000)
  common/BIG/ALLGW(44000)
  common/BIG/TP(6000)
  common/BIG/TS(6000)
  common/BIG/TBC(6000)
  real*8 Q, Q1D, TP, TBC, SUM
  integer IYR1D, IMO1D, IDA1D, TS
  integer iyearst, iyearen, IBEFORE
  character*1 EST
  character*1 EST1D
  character*1 ALLGW

  IBEFORE=iyr(1) - 1

  !
  !------------------------- INITIALIZE VARIABLES : -------------------
  !
  IERR = 0
  Q = -999.0D0
  ALLGW = ' '
  Q1D = -999.0D0
  !
  !
  ! -------- PROCESS THE DAILY-VALUES FILE OF STREAMFLOW: ------------------
  !
  IFRSTYR=iyr(1)
  do i=1,n
    IYEAR = IYR(i)-IBEFORE
    Q(IYEAR,imon(i),idy(i)) = FLONUM(i)
  end do

  ILSTYR= IYR(n)
  !
  ! --- flag nonexistent dates WITH flow=-9999 -----
  !
  do IYEAR=IFRSTYR, ILSTYR
    do IMONTH=1,12
      do IDAY=1,31
        if((IMONTH.eq.2).and.(IDAY.gt.29)) then
          Q(IYEAR-IBEFORE,IMONTH,IDAY)= -9999.0D0
        end if
        if((IMONTH.eq.2).and.(IDAY.eq.29)) then
          IDIV=int((IYEAR)/4.0)
          XDIV=(IYEAR)/4.0
          DIFFER=abs(IDIV-XDIV)
          if(DIFFER.gt.0.1) then
            Q(IYEAR-IBEFORE,IMONTH,IDAY)= -9999.0D0
          end if
        end if
        if(IDAY.eq.31) then
          if(  ( IMONTH == 4 )           &
           .or.( IMONTH == 6 )           &
           .or.( IMONTH == 9 )           &
           .or.( IMONTH == 11 ) ) then

            Q(IYEAR-IBEFORE,IMONTH,IDAY)= -9999.0D0

          end if
        end if
      end do
    end do
  end do
  IFRSTYR= IFRSTYR-IBEFORE
  ILSTYR= ILSTYR-IBEFORE
  !
  !   ---------------   SELECT TIME PERIOD OF INTEREST:  ------------------
  !
  IIMAX=0
  do IIYR=IYRST,IYREN
    IDIV=int(IIYR/4.0D0)
    XDIV=IIYR/4.0D0
    DIFFER=abs(IDIV-XDIV)
    if(DIFFER.lt.0.1D0) then
      IIMAX=IIMAX+366
    else
      IIMAX=IIMAX+365
    end if
  end do
  IYEARST= IYRST-IBEFORE
  IYEAREN= IYREN-IBEFORE
  !
  ! ---- ASSIGN VALUES TO 1-DIMENSIONAL ARRAYS OF DISCHARGE AND DATE: ----
  !
  ICOUNT= 0
  IBREAK= 0
  ITESTX= 0
  do IYEAR= IYEARST, IYEAREN
    do IMONTH= 1,12
      do IDAY= 1, 31
        SFLOW= Q(IYEAR,IMONTH,IDAY)
        if(SFLOW .eq. -9999.D0) cycle
        if(SFLOW .eq. -99D0 .or. SFLOW .eq. -999D0) then
          ITEST=0
        else
          ITEST=1
        endif
        if(ITEST.eq.1.and.ITESTX.eq.0) then
          IBREAK= IBREAK+1
          if(IBREAK.gt.1) then
            IERR = 2
            return
          endif
        endif
        ITESTX= ITEST
        if ( ( SFLOW == -99D0 ) .or. ( SFLOW == -999D0 ) ) cycle

        ICOUNT= ICOUNT + 1
        Q1D(ICOUNT)= SFLOW
        IYR1D(ICOUNT)= IYEAR
        IMO1D(ICOUNT)= IMONTH
        IDA1D(ICOUNT)= IDAY
        EST1D(ICOUNT)= EST(IYEAR,IMONTH,IDAY)

      end do
    end do
  end do

  IZERO=1
  do 185 I=1,ICOUNT
    if(Q1D(I).eq.0.0) IZERO=0
185 ALLGW(I)= ' '
    !
    ! ----------- ModIFy 0 flows
    !
    if (IZERO.eq.0) then
      do 187 I=1,ICOUNT
        if(Q1D(I).eq.0.00D0) then
          Q1D(I)= XNEWMIN
        end if
187     continue
      end if

      !
      ! ------------ Process DA
      !
      if (DA.lt.1.0) IERR = -2
      if (DA.gt.500.0D0) IERR = -3

      !
      !  ---  DETERMINE THE MINIMUM NUMBER OF DAYS OF ANTECEDENT RECESSION  --
      !  ---  TO INDICATE THAT STREAMFLOW IS TO BE CONSIDERED GROUND-WATER  --
      !  ---  DISCHARGE.  OBTAIN FROM THE EQUATION DA**0.2, ROUNDED UPWARD  --
      !
      ITBASE=0
210   ITBASE= ITBASE + 1
      if(ITBASE .gt. 10) then
        IERR = 3
        return
      end if
      if (ITBASE.gt.DA**0.2D0) then
        GO TO 220
      else
        GO TO 210
      end if
220   continue

      !
      ! ----- SPECIFY THE MAXIMUM ALLOWABLE NUMBER OF DAYS THAT CAN BE USED
      ! ----- AFTER A PEAK, TO DETERMINE THE GROUND-WATER DISCHARGE AFTER
      ! ----- THE PEAK:
      !
      IRECMAX= int(0.2144D0*K)

      ! ---- ALLOW THE USER TO OVERRIDE THE DEFAULT METHOD OF EXECUTION.  -----
      ! ---- NOTE: THE ALLOWANCE FOR NON-DEFAULT EXECUTION, IS DIFFERENT  -----
      ! ---- FROM THE NON-DEFAULT EXECUTION IN PREVIOUS VERSIONS OF RORA ------

      if (IDIFF .lt. 0 .or. IDIFF .gt. 3) then
        IERR = 1
        return
      else
        ITBASE= ITBASE+IDIFF
      endif

      if(IRECMAX.lt.ITBASE) then
        IERR = -1
        IRECMAX= ITBASE
      end if
      !
      !-----------------ON DATES PRECEEDED BY A RECESSION PERIOD,-------------
      !-------------------------SET VARIABLE ALLGW='*'  ----------------------
      !
      do 270 I= ITBASE+1, ICOUNT
        INDICAT=1
        IBACK= 0
260     IBACK= IBACK + 1
        if(Q1D(I-IBACK).lt.Q1D(I-IBACK+1)) INDICAT=0
        if(IBACK.lt.ITBASE) GO TO 260
        if(INDICAT.eq.1) then
          ALLGW(I)= '*'
        else
          ALLGW(I)= ' '
        end if
270     continue
        !
        !
        ! -- THE FOLLOWING LINES UP TO LABEL 800 DETERMINE THE LOCATION (IN TIME) --
        ! --------------------  OF PEAKS AND RECESSION PERIODS:  -------------------
        !
        !
        ! ------------------- LOCATE end OF FIRST RECESSION: ----------------
        !
        I= 0
280     I= I+1
        if(I.gt.ICOUNT) then
          IERR = 4
          return
        end if
        if (ALLGW(I).ne.'*') GO TO 280
290     I=I+1
        if (I.gt.ICOUNT) then
          IERR = 5
          return
        end if
        if (ALLGW(I).ne.' ') GO TO 290
        I=I-1
        TA(1)= I
        QA(1)= Q1D(I)
        IPEAK= 0
        !
        ! -------------    FIND STREAMFLOW AND DAY OF THE PEAK:     -----------
        !
330     continue
        IPEAK= IPEAK+1
        if(IPEAK .gt. 6000) then
          IERR = 6
          return
        end if
        QP(IPEAK)= Q1D(I)
        TP(IPEAK)= I
        ILOOK= I + 1
340     continue
        if(ILOOK.ge.ICOUNT) GO TO 800
        if(ALLGW(ILOOK).eq.'*') GO TO 350
        if (Q1D(ILOOK).ge.QP(IPEAK)) then
          QP(IPEAK)= Q1D(ILOOK)
          TP(IPEAK)= ILOOK
        end if
        if(ALLGW(ILOOK).eq.' ') then
          ILOOK= ILOOK+1
          GO TO 340
        end if
350     continue
        TBC(IPEAK)= TP(IPEAK) + (0.2144*K)

        !
        !     -----   FIND FIRST AND LAST DAYS OF RECESSION FOLLOWING THE PEAK: ----
        !
        I= ILOOK
        TS(IPEAK)= I
370     continue
        I= I+1
        if(I.gt.ICOUNT) GO TO 380
        if(ALLGW(I).ne.' ') then
          GO TO 370
        end if
380     continue
        I= I-1
        TE(IPEAK)= I
        if (TE(IPEAK)-TP(IPEAK).gt.IRECMAX) then
          TE(IPEAK)= TP(IPEAK) + IRECMAX
        end if
        if (TE(IPEAK).lt.TS(IPEAK)) then
          TE(IPEAK)= TS(IPEAK)
        end if

        I=I+1
        if(I.lt.ICOUNT) GO TO 330

        !
        ! ----  GO HERE WHEN ALL RECESSION PERIODS AND PEAKS HAVE BEEN LOCATED: ----
        !
800     continue

        NPEAKS= IPEAK - 1

        !
        ! --- EXTRAPOLATE STREAMFLOW, DETERMINE RECESSION CURVE DISPLACEMENT, AND --
        ! ----------------- CALCULATE RECHARGE FOR FIRST PEAK:  --------------------
        !
        QB(1)= QA(1)*10**(-1*(TBC(1)-TA(1))/K)
        SUM= 0.0D0
        do 810 I= TS(1), TE(1)
          DQ= Q1D(I) - (QA(1)*10D0**(-1.0D0*(I-TA(1))/K))
810       SUM= SUM + (DQ*sqrt(I-TP(1)))
          NRECS= TE(1)-TS(1)+1
          C(1)= SUM/NRECS
          DELQ(1)= C(1)/(sqrt(TBC(1)-TP(1)))
          RECH(1)= 2.0D0*0.0371901D0*DELQ(1)*K / (2.3025851D0*DA)
          QC(1)= QB(1) + DELQ(1)
          year(1) = IBEFORE+IYR1D(int(TP(1)))
          mon(1) = IMO1D(int(TP(1)))
          day(1) = IDA1D(int(TP(1)))
          !
          ! -- EXTRAPOLATE STREAMFLOW, DETERMINE RECESSION CURVE DISPLACEMENT, AND --
          ! --------------- CALCULATE RECHARGE FOR ALL OTHER PEAKS: -----------------
          !
          do 840 IPEAK=2, NPEAKS
            QA(IPEAK)= QC(IPEAK-1)
            TA(IPEAK)= TBC(IPEAK-1)
            QB(IPEAK)= QA(IPEAK)*10**(-1*(TBC(IPEAK)-TA(IPEAK))/K)

            SUM= 0.0D0
            do 820 I= TS(IPEAK), TE(IPEAK)
              if(I.gt.TA(IPEAK)) then
                BASELIN= QA(IPEAK)*10**(-1*(I-TA(IPEAK))/K)
              else
                BASELIN= C(IPEAK-1)/(sqrt(I-TP(IPEAK-1)))                      &
                         + QA(IPEAK-1)*10D0**(-1.0D0*(I-TA(IPEAK-1))/K)
              end if
              DQ= Q1D(I) - BASELIN
820           SUM= SUM + (DQ*sqrt(I-TP(IPEAK)))
              NRECS= TE(IPEAK) - TS(IPEAK) + 1
              C(IPEAK)= SUM/NRECS
              DELQ(IPEAK)= C(IPEAK) / (sqrt(TBC(IPEAK)-TP(IPEAK)))
              RECH(IPEAK)= 2.0D0*0.0371901D0*DELQ(IPEAK)*K / (2.3025851D0*DA)
              QC(IPEAK)= QB(IPEAK) + DELQ(IPEAK)
              year(ipeak) = IBEFORE+IYR1D(int(TP(IPEAK)))
              mon(ipeak) = IMO1D(int(TP(IPEAK)))
              day(ipeak) = IDA1D(int(TP(IPEAK)))
840           continue
              !
              ! ----- Adjust TE and TA to be relative to TP
              !
              do i=1, npeaks
                TE(i) = TE(i) - int(tp(i))
                TA(i) = TA(i) - tp(i)
              end do

              return
              
end subroutine roraf
