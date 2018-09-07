      SUBROUTINE RADIATORS(T)                                           
      IMPLICIT NONE
      COMMON /SET/   E0SET,EPSET,THSET 
      REAL  E0SET,EPSET,THSET                                    
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                    
      COMMON /CON/   PI,PM,PP,EM,AL  
      REAL PI,PM,PP,EM,AL                                    
      COMMON /OMG/   DELTA,R 
      REAL   DELTA,R                                             
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA     
      COMMON /TRL/TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG 
      REAL  TTARG,TWALL,TBEAM,TSPEC
      INTEGER NSEG,JTARG,ITARG,nprnt          
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM,thwin
      REAL DELP1,DEL01,DEL02,DELP2,DELP,TBEFOR,TAFTER,TAAL,TBAL,
     >  TEQUI,DEL0,T,x_end,ztarg,cryo_cm,wall_cm,cell_diam,
     >  cell_wall
      REAL ECIR,ECOR,ENTEC,TARGLEN,tcm
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      save                                                                        
C If Internal corrections only, no target material only equiv. radiator:
                                                                        
      TBEFOR = T                                                        
      TAFTER = 0.                                                       
      TBAL   = 0.                                                       
      TAAL   = 0.                                                       
      ATBX0  = 0.                                                       
      ATAX0  = 0.                                                       
      BTBI   = 0.                                                       
      BTAI   = 0.                                                       
      DEL0   = 0.                                                       
      DELP   = 0.                                                       
      IF (JTARG.EQ.1) GO TO 1                                           
                                                                        
C Radiator quantities. Ionization spectrum should not have equivalent   
C radiator. Correct factors AX0 are used for the aluminium and target   
C parts:                                                                

C     IF (TARGET.EQ.'V_CYL') CALL V_CYL   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E_087') CALL E_087   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E_089') CALL E_089   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E_140') CALL E_140   (T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'E140XH1') CALL E140XH1 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'E140XH2') CALL E140XH2 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'E140XD1') CALL E140XD1 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'E140XD2') CALL E140XD2 (T,THR,TAFTER,TBAL,TAAL)    
C     IF (TARGET.EQ.'OTHER') CALL USERTARG(T,THR,TAFTER,TBAL,TAAL)      
C     IF (TARGET.EQ.'SOLID') TAFTER = (TTARG-T)/COS(THR)                
C     TBAL  = TBAL+TBEAM                                                
C     TAAL  = TAAL+TSPEC                                                
C
c default for a simple target
      TBEFOR=t                        ! target
      TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th)
      TBAL= TBEAM                     ! material before the target
      TAAL= TSPEC                     ! material after the target

      IF(INDEX(TARGET,'E154').GT.0) THEN                    
       CALL RADLEN154(TARGET,T,TH,TAFTER,TBAL,TAAL,TBEFOR)  
      ELSEIF(INDEX(TARGET,'E142').GT.0) THEN                                
       CALL RADLEN42(TARGET,T,THR,TAFTER,TBAL,TAAL)              
      ELSEIF(INDEX(TARGET,'E143').GT.0) THEN                            
       CALL RADLEN43(TARGET,T,THR,TAFTER,TBAL,TAAL)              
      ELSEIF(INDEX(TARGET,'E149').GT.0) THEN                            
       CALL RADLEN49(TARGET,T,THR,TAFTER,TBAL,TAAL)              
c     ELSEIF(ITARG.LE.7) THEN                                           
c      CALL RADLENGT2(ITARG,T,THR,TAFTER,TBAL,TAAL)                     
      ENDIF
      IF(INDEX(TARGET,'EG1b').GT.0) THEN  ! Eg1btargets
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th)
        TBAL= 0.0024               ! Three 71 um Al windows (incl. banjo)
! get angle of exit window assuming radius of curvature is
! 3.5 times larger than distance to window (i.e. 35 inches
! versus 10 inches from target to window)
        thwin = thr / 3.5
        TAAL= 0.0014 / cos(thr)  + ! Banjo exit, 4K and 100K shields
     >    0.0031 / cos(thr + thwin) + ! 280 um exit widnow bowed inward
     >    0.0016 +                                    ! DC1
     >   (172. - th * 1.3 - 23.1 / cos(thr)) / 30420.  ! air
! extra 50 um of aluminum in 4K shield beyond 29 deg
        if(th.gt.29.) TAAL = TAAL + 0.00056
! extra 50 um of aluminum in 100K shield beyond 21 deg
        if(th.gt.21.) TAAL = TAAL + 0.00056
! 12 layers of superinsulation:
        if(th.gt.19.) TAAL = TAAL + 0.00029

        nprnt= nprnt+1
        if(nprnt.le.100) write(*,111) tbal,tbefor,tafter,taal
 111    format(1x,'EG1b rad len=',4f10.5)
      endif

      IF(INDEX(TARGET,'EG4').GT.0) THEN  ! Eg4 targets
        call eg4radlen(th,tbal,taal)
        nprnt= nprnt+1
        if(nprnt.le.100) write(*,122) tbal,tbefor,tafter,taal
 122    format(1x,'EG4 rad len=',4f10.5)
      endif

! New for E99-118
      IF(INDEX(TARGET,'E99S').GT.0) THEN  ! solid targets
        TBEFOR=t/cos(0.354)             ! target/cos(20.3 deg)
        TAFTER=(ttarg-t)/cos(thr-0.354) ! target(cos(th-20.3 deg)
        TBAL= TBEAM                     ! material before the target
        TAAL= TSPEC                     ! material after the target
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,112) tbal,tbefor,tafter,taal
 112    format(1x,'E99S rad len=',4f10.5)
      endif
      IF(INDEX(TARGET,'E99L').GT.0) THEN ! liquid in LD2 or LH2
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th)
        TBAL= TBEAM                     ! material before the target
        x_end = 4.095 * (ttarg-t)/ttarg * tan(thr)   ! position at endcap
        if(x_end.lt.2.0) then           ! went through endcap
          TAAL= TSPEC                     ! material after the target
        else
          taal=twall/sin(thr)+0.00606     ! went through side wall
        endif
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,113) tbal,tbefor,tafter,taal
 113    format(1x,'E99L rad len=',4f10.5)
      endif
      IF(INDEX(TARGET,'E99D').GT.0) THEN  ! Dummy aluminum 
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th): default
        x_end = 4.09 * tan(thr)         ! position at 2nd endcap from first
        if(x_end.gt.2.38) then             ! misses 2nd endcap +/- 2.38 cm   
          if(abs(t-ttarg/2.).lt.0.00001) then ! middle
            tafter = ttarg/2.0/2.0/cos(thr)  ! average
          else
            if(t.lt.ttarg/2.-0.00001) tafter=(ttarg/2.-t)/cos(thr)
            if(t.gt.ttarg/2.+0.00001) tafter=(ttarg-t)/cos(thr)
          endif
        endif
        TBAL= TBEAM                     ! material before the target
        TAAL= TSPEC                     ! material after the target
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,114) t/ttarg,tbal,tbefor,tafter,taal
 114    format(1x,'E99D rad len=',5f10.5)
      endif
      IF(INDEX(TARGET,'E99E').GT.0) THEN  ! for LH2, LD2 endcaps
        TBEFOR=t                        ! target
        TAFTER=(ttarg-t)/cos(thr)       ! target/cos(th): default
        x_end = 4.09 * tan(thr)         ! position at 2nd endcap from first
        if(x_end.gt.2.0) then           ! misses  2nd endcap
          if(abs(t-ttarg/2.).lt.0.0001) then ! middle
            tafter = ttarg/2.0/2.0/cos(thr)  ! average
          else
            if(t.lt.ttarg/2.-0.0001) tafter=(ttarg/2.-t)/cos(thr)
            if(t.gt.ttarg/2.+0.0001) tafter=(ttarg-t)/cos(thr)
          endif
        endif
        TBAL= TBEAM                      ! material before the target
        TAAL= TSPEC                      ! material after the target
! add liquid and walls to tbal
        if(t.lt.ttarg/2.-0.0001) then
          if(x_end.lt.2.0) then
            taal = taal + 0.0052/cos(thr)
          else
            taal = taal + (2.0/sin(thr)/890.+twall/sin(thr))
          endif
        endif
        if(t.gt.ttarg/2.+0.0001) then
          tbal = tbal + 0.0052
        endif
        if(abs(t-ttarg/2.).lt.0.0001) then ! middle
          tbal = tbal + 0.0052/2.0
          if(x_end.lt.2.0) then
            taal = taal + 0.0052/cos(thr)/2.0
          else
            taal = taal + (2.0/sin(thr)/890.+twall/sin(thr))/2.
          endif
        endif
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,115) t/ttarg,tbal,tbefor,tafter,taal
 115    format(1x,'E99E rad len=',5f10.5)
      endif
! For Tuna Can targets, make sure to use correct input values
! for ttarg, tbeam, and tspec. twall is not used.
      IF(INDEX(TARGET,'TUNA').GT.0) THEN ! liquid in LD2 or LH2
        TBEFOR=t                        ! target
c note: if change cell length, make sure ttarg matches for
c       actual density and material of target
        cell_diam = 3.96
        cell_wall = 0.0125
! z in cm relative to center 
        ztarg = cell_diam * (t / ttarg - 0.5) 
        call tunatarg(cell_diam, cell_wall, 
     >    ztarg, thr, cryo_cm, wall_cm)
        TAFTER = ttarg * cryo_cm / cell_diam
        TBAL= TBEAM                     ! material before the target
        TAAL= TSPEC + wall_cm / 8.9     ! material after the target
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,123) tbal,tbefor,tafter,taal,cryo_cm
 123    format(1x,'TUNA rad len=',5f10.5)
      endif

      IF(INDEX(TARGET,'CRYO17').GT.0) THEN     ! target geometry for
                                               ! 10 cm LH2 and LD2
        ECIR=1.315*2.54                        ! endcap inner radius (cm
        ECOR=1.320*2.54                        ! endcap outer radius (cm
        TARGLEN=3.942*2.54                     ! target length (cm
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm
        
        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! liquid target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen ! liquid target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'10OLD17').GT.0) THEN    ! target geometry for
                                               ! old 10 cm LH2 and LD2
        ECIR=0.795*2.54                        ! endcap inner radius (cm
        ECOR=0.800*2.54                        ! endcap outer radius (cm
        TARGLEN=10                             ! target length (cm
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! liquid target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen ! liquid target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
      ENDIF

      IF(INDEX(TARGET,'TRITIUM').GT.0) THEN     ! target geometry for
                                               ! 25 cm tritium cell
        ECIR=0.635                             ! endcap inner radius (cm)
        ECOR=0.635+0.0324                      ! endcap outer radius (cm)
        TARGLEN=25.0                           ! target length (cm)
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm)
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm
        
        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! target gas
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes through end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen !  target gas
c          if((TAFTER*targlen/ttarg)*sin(thr).gt.0.5) then
c             ECOR=0.635+0.04279 ! thicker part of exit cap
c          endif
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
        nprnt= nprnt+1
        if(nprnt.le.6) write(*,125) tcm,tbal,tbefor,tafter,taal
 125    format(1x,'TRITIUM rad len=',5f10.5)
      ENDIF


      IF(INDEX(TARGET,'4OLD17').GT.0) THEN     ! target geometry for
                                               ! old 4 cm LH2 and LD2
        ECIR=0.795*2.54                        ! endcap inner radius (cm
        ECOR=0.800*2.54                        ! endcap outer radius (cm
        TARGLEN=4                              ! target length (cm
        ENTEC=TARGLEN-ECIR                     ! entrance to end cap (cm
        TBAL=TBEAM                             ! material bf the target
        TBEFOR=T                               ! distance traveled bf
                                               ! vertex
        tcm=t*targlen/ttarg                    ! t in cm

        if((tcm+ECIR/tan(thr)).lt.ENTEC) then  ! e goes through sidewall
          TAFTER=ECIR*ttarg/sin(thr)/targlen   ! liquid target
          TAAL=TSPEC+TWALL/sin(thr)            ! material af the target 
                                               ! & side wall
        else
          TAFTER=                              ! e goes throught end cap
     >     (sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    +(targlen-ECIR-tcm)*cos(thr))*ttarg/targlen ! liquid target
          TAAL=TSPEC                           ! material af the target
     >    +(sqrt(ECOR**2-((targlen-ECIR-tcm)*sin(thr))**2)
     >    -sqrt(ECIR**2-((targlen-ECIR-tcm)*sin(thr))**2))
     >    *TWALL/(ECOR-ECIR)                   ! & end cap
        endif
      ENDIF

      ATBX0 = AX0*TBEFOR+AX0A*TBAL                                      
      ATAX0 = AX0*TAFTER+AX0A*TAAL                                      
      BTBI  = B*TBEFOR+BA*TBAL                                          
      BTAI  = B*TAFTER+BA*TAAL                                          
C Redefine Initial and Final energies subtracTIng the most probable     
C ionization loss and calculate kinematic dependent quantities:         
                                                                        
      DEL01 = 0.                                                        
      DEL02 = 0.                                                        
      DELP1 = 0.                                                        
      DELP2 = 0.                                                        
      IF (TBEFOR.GT.0.) DEL01 = AX0*TBEFOR                              
     >        *(LOG(3.E9*AX0*TBEFOR*E0SET**2/EM**2/iZ**2)-0.5772)       
      IF (TBAL.GT.0.) DEL02 = AX0A*TBAL                                     
     >        *(LOG(3.E9*AX0A*TBAL*E0SET**2/EM**2/13.**2)-0.5772)       
      IF (TAFTER.GT.0.) DELP1 = AX0*TAFTER                              
     >        *(LOG(3.E9*AX0*TAFTER*EPSET**2/EM**2/iZ**2)-0.5772)       
      IF (TAAL.GT.0.) DELP2 = AX0A*TAAL                                 
     >        *(LOG(3.E9*AX0A*TAAL*EPSET**2/EM**2/13.**2)-0.5772)       
      DEL0  = DEL01+DEL02                                               
      DELP  = DELP1+DELP2                                               
                                                                        
ccc pyb TURNED OFF energy loss smearing!
      DEL0 = 0.
      DELP = 0.

1     E0 = E0SET-DEL0                                                   
      EP = EPSET+DELP                                                   
      Q2 = 4.*E0*EP*SIN(THR/2.)**2                                      
      R  = (PM+E0*(1-COS(THR)))/(PM-EP*(1-COS(THR)))                    
C Equivalent radiator for internal bremstrahlung correction             
      TEQUI = AL/PI*(LOG(Q2/EM**2)-1.)/B                                
                                                                        
      TB  = TBEFOR+TBAL+TEQUI                                           
      TA  = TAFTER+TAAL+TEQUI                                           
      BTB = B*(TBEFOR+TEQUI)+BA*TBAL                                    
      BTA = B*(TAFTER+TEQUI)+BA*TAAL                                    
c      write(6,'(''rad'',10f7.3)') b,tequi,tb,ta,btb,bta
      RETURN                                                            
      END                                                               
                                                                      
C=======================================================================
                                                                        
      SUBROUTINE RADLEN42 (TARGET,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E142  radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the in cap of the target cell.             
!  T1G and T2G are same as T1 and T2 but in grams instead of r.l.       
!  Writes diagnostic and error messages to unit IU.                     
!                                                                       
!  Calculations are carried out in target coordinate system where the   
!  x axis is along the axis of the target and y axis perpendicular to   
!  x axis with positive values towards 8 GeV spectrometer.  The origin o
!  this coordinate system is at the in cap of the target cell.          
!  The out cap of the target cell is assumed to be ellipsoidal given by 
!  (B*(X-XC))**2 + (A*Y)**2 - (A*B)**2 = 0                              
!                                                                       
!  Changed r.l. of LH2 to 61.28 (1988 value) from 63.                   
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL LHE3/0./,LWALL/0./,LENDCAP/0./,XCM,TTARG,TWALL,TBEAM,TSPEC   
      REAL X,TR,TBAL,TAAL,XMAX/0./                                      
      REAL TAFTER         ! target material after scattering.           
      REAL RLENDCAP/0./   ! radiation length in end cap after scatter   
      REAL ECRL_F/0./, ECRL_R/0./ !rad length of end cap, front and rear
      INTEGER ITARG,NSEG,JTARG
      CHARACTER*7    TARGET                                             
      REAL      DAT(4,11)                                               
                                                                        
      REAL RADIUS/.825/                        
      LOGICAL FIRST/.TRUE./,FIRST1/.TRUE./                              
      LOGICAL INHE,INEC                                                 
      DATA    DAT/                                                      
! TARGPARAM: radiation lengths of various types of Al from P. 52        
! Info on wire arrayS: Jerry Davis says 5 mil diamter Al wire 25 mil spa
!              and there are two arrays.                                
! Info on Hymen: John Mark says it was 1 mil Al.                        
! Info on In-cap: from Hunter.  .01 to .012 inches total for both ends  
! Info on He3: from Hunter 2.32E20 atoms/cc.                            
!              30 cm long                                               
! Info on cell wall. Value of 4.4 mils is from p. 24                    
!              Actual value varies with x and is hardwired into program.
! Infor on endcap comes Hunter.                                         
! Info on insulation comes from J. Mark (10 layers of 1/4 mil mylar)    
!              it pretty much lines up with endcap of long targets,     
!              so for short targets need to take this into account.     
! Info on 1.6 gev window comes from p. 24 (sample)                      
! Thickness of 8 GeV exit window is 10.7 mils when not crinkled (Eislele
!                                                                       
!       Thickness     density      X0      Z    Name       Materl       
!         (cm)        (g/cm3)    (g/cm2)                                
     >    0.00400,    2.7,       24.0111, 13., 
     >    0.00254,    2.7,       24.0111, 13.,     
     >    0.01400,    2.52,      27.0000, 14.,       
     >    30.0280,    0.00116,   71.0000,  2.,         
     >    0.014,      2.52,      27.,     14.,       
     >    0.014,      2.52,      27.,     14.,       
     >    0.00635,    1.39,      39.95,   5.2,       
     >    0.00762,    2.68,      23.6311, 13.,     
     >    0.03048,    2.68,      23.6311, 13.,     
     >    16.0000,    0.001205,  36.975,  7.3,        
     >    0.02540,    2.68,      23.6311, 13./     
                                                                        
                                                                        
      IF(FIRST) THEN                                                    
       FIRST =.FALSE.                                                   
       DAT(2,4) = TTARG * DAT(3,4) / DAT(1,4) !Density = X0 *rl/len     
       DAT(1,5) = TWALL * DAT(3,5) / DAT(2,5) !len(cm)=  X0 *rl/den     
       XMAX = DAT(1,4) - DAT(1,3) - DAT(1,6)  !targ len minus end caps  
       ECRL_F = DAT(1,3) * DAT(2,3)/DAT(3,3)  !radiation length front ec
       ECRL_R = DAT(1,6) * DAT(2,6)/DAT(3,6)  !radiation length rear ec 
       WRITE( 6,'('' ****USING E142 TARGET MODEL '')')                  
      ENDIF                                                             
      INHE =.FALSE.                                                     
      IF(INDEX(TARGET,'E142F').GT.0) THEN  ! Gas in target goes thru wal
        INHE=.TRUE.                                                     
! Material after the scatter.                                           
! target model is a cylinder with flat end caps.                        
! Length in Helium                                                      
        LHE3 = RADIUS/SIN(TR)    ! in cm                                
        LWALL = DAT(1,5)/SIN(TR)                                        
        LENDCAP = 0.                                                    
      ELSEIF(INDEX(TARGET,'E142R').GT.0) THEN  ! exits thru endcap      
        INHE=.TRUE.                                                     
        XCM=X*DAT(3,4)/DAT(2,4)         ! changed to unit cm            
        LHE3 = (XMAX-XCM)/COS(TR)                                       
        LWALL =0.                                                       
        LENDCAP = DAT(1,6)/COS(TR)                                      
      ENDIF                                                             
      IF(INHE) THEN                                                     
! Radiation  lengths before entering the gas.                           
       TBAL = TBEAM + ! material before the target not including end ca 
     >        ECRL_F                       ! front end cap              
       TAAL= TSPEC +            !Non He3 material after scatter         
     >      LWALL  * DAT(2,5)/DAT(3,5) +                                
     >      LENDCAP* DAT(2,6)/DAT(3,6)                                  
       TAFTER = LHE3   * DAT(2,4)/DAT(3,4)  ! He3 after scatter         
      ENDIF                                                             
                                                                        
!End caps                                                               
      INEC=.FALSE.                                                      
      IF(INDEX(TARGET,'E142ECF').GT.0) THEN  !Front End Cap radiation co
        INEC=.TRUE.                                                     
        LHE3 = RADIUS/SIN(TR)                                           
        LWALL= DAT(1,5)/SIN(TR)                                         
        RLENDCAP = (ECRL_F -X)/COS(TR)    !front end cap after scatter  
        TBAL = TBEAM      !material before target                       
      ELSEIF(INDEX(TARGET,'E142ECR').GT.0) THEN  !rear end cap          
        INEC=.TRUE.                                                     
        LHE3 = 0.                                                       
        LWALL =0.                                                       
        RLENDCAP = (ECRL_F + ECRL_R -X)/COS(TR) !rear ec after scatter  
        TBAL = TBEAM +                                                  
     >        ECRL_F                       ! front end cap              
      ENDIF                                                             
                                                                        
      IF(INEC) THEN                                                     
       TAAL= TSPEC +            !material after target                  
     >      LWALL  * DAT(2,5)/DAT(3,5) +                                
     >      RLENDCAP                                                    
       TAFTER = LHE3   * DAT(2,4)/DAT(3,4)  ! He3 after scatter         
      ENDIF                                                             
                                                                        
      IF((.NOT.INEC) .AND. (.NOT.INHE))THEN !no target model            
       WRITE(6, '(''TARGET='',A8,'' NOT ALLOWED'')') TARGET             
       WRITE(10,'(''TARGET='',A8,'' NOT ALLOWED'')') TARGET             
       STOP                                                             
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
C====================================================================== 
                                                                        
                                                                        
                                                                        
      SUBROUTINE RADLEN43 (TARGET,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E143  radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the beginning of the target material.      
!  Steve Rock 4/9/93                                                    
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL TTARG ! total of amonia (nitrogen and H2 and He3 coolant     
      REAL TWALL,TBEAM,TSPEC                                            
      REAL X,TR,TBAL,TAAL                                         
      REAL TAFTER         ! target material after scattering.           
      INTEGER ITARG,NSEG,JTARG
      CHARACTER*7    TARGET 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM                                             
      LOGICAL FIRST/.TRUE./
      CHARACTER*(*) STRING
      REAL TD                                                                        
                                                                        
      IF(FIRST) THEN                                                    
       FIRST =.FALSE.                                                   
       WRITE(6,'('' ENTERED RADLEN43'')')                               
      ENDIF                                                             
                                                                        
! Radiation  lengths before entering the gas.                           
       TBAL = TBEAM   ! material before the target amonia               
       TAAL= TSPEC    ! material after the target amonia.               
       TAFTER = TTARG -X   ! Amonia after scattering                    
       TD = 57.296*TR
! 10 deg spectrometer fit from test_targ2 6/2/98 with 302 cm of air in mags
       IF(TD.GE.8.3) THEN  
        IF((IA.EQ.14).OR.(IA.EQ.1)) THEN   ! N or p
         TAAL = -1.3171E-01 + 5.1404E-02*TD -5.7571E-03*TD**2 + 
     >     2.1474E-04*TD**3
        ELSEIF((IA.EQ.6).OR.(IA.EQ.2)) THEN  ! Li or d
         TAAL =  -1.3200E-01 + 5.1496E-02*TD -5.7681E-03*TD**2 + 
     >     2.1506E-04*TD**3
        ELSEIF((IA.EQ.9).or.(IA.EQ.12)) THEN   ! Long Be and long C
         TAAL =  -1.2472E-01  +4.9967E-02*TD -5.6091E-03*TD**2 +
     >     2.0969E-04*TD**3

        ENDIF
       ENDIF                                                                 
                                                                        
      RETURN                                                            
     
      ENTRY RADLEN43_INIT(STRING)
         STRING='** 10 deg T-Spec from fit in Targ_test2'
      RETURN

      END                                                               
                                                                        
                                                                        
C====================================================================== 
                                                                        
                                                                        
                                                                        
      SUBROUTINE RADLEN49 (TARGET,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E149  radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the beginning of the target material.      
!  Steve Rock 4/26/93                                                   
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL TTARG ! total of amonia (nitrogen and H2 and He3 coolant     
      REAL TWALL,TBEAM,TSPEC                                            
      REAL X,TR,TBAL,TAAL
      REAL TAFTER         ! target material after scattering.           
      INTEGER ITARG,NSEG,JTARG                                                
      CHARACTER*7    TARGET                                             
!*** changed 4/8/03 pyb
      REAL RADIUS /3./                                                  
      LOGICAL FIRST/.TRUE./                                             
                                                                        
                                                                        
      IF(FIRST) THEN                                                    
       FIRST =.FALSE.                                                   
       WRITE(6,'('' ENTERED RADLEN49'')')                               
      ENDIF                                                             
                                                                        
! Radiation  lengths before entering the gas.                           
       TBAL = TBEAM   ! material before the target amonia               
       TAAL= TSPEC    ! material after the target amonia.               
       TAFTER = RADIUS/SIN(TR)/757. ! radiation lengths leaving target  
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                       
C=====================================================================  
                                                                        
      SUBROUTINE RADLEN154 (TARGET,X,THETAD,TAFTER,TBAL,TAAL,TBEFOR)
!_______________________________________________________________________
!  Calculates E154  radiation length (units of rl) for all materials   
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TR (rad) and scattering taking place at X r.l along the axis   
!  of the target relative to the in cap of the target cell.             
!  T1G and T2G are same as T1 and T2 but in grams instead of r.l.                                             
      IMPLICIT NONE   
      CHARACTER*7    TARGET
      REAL X,THETAD,TBAL,TAAL !TH is scattering angle in degrees
      REAL TAFTER         ! target material after scattering.
      REAL TBEFOR         ! target material before scattering.  
      REAL TTOT,TLEN
      INTEGER ISEG154
      COMMON/E154/ISEG154,TTOT,TLEN    

C For Picard: from linda
       REAL*8 TB2(4)/0.00054D0, 0.00060D0, 0.00065D0, 0.00082D0/,
     >        TA2(4)/0.14631D0, 0.19431D0, 0.04434D0, 0.00067D0/,
     >        TLEN2(4)/3.D0, 5.D0, 2.D0, 20.D0/,
     >        TB5(4)/0.00062D0, 0.00075D0, 0.00081D0, 0.00090D0/,
     >        TA5(4)/0.07416D0, 0.09616D0, 0.04416D0, 0.00059D0/,
     >        TLEN5(4)/14.D0, 4.D0, 2.D0, 10.D0/

!use rho_ec =2.52, X0(rad length of glass)=28 gms/cm2 (guess from Wallet Cards
!end cap front = 69um =.00062 r.l.
!end cap rear  =61.6um .00055 r.l.
      REAL T_ECF/.00062/
      REAL T_ECR/.00055/

!assume all material is in end caps
      TAFTER=0.
      TBEFOR=0.
                            
         IF(THETAD.LT.3.5) THEN
           IF(INDEX(TARGET,'ECF').NE.0) THEN  ! front end cap
            TBAL= .5*TB2(1)
            TAAL = .5*TB2(1) +TA2(1)
            TLEN=1.
            TTOT=1.
           ELSEIF(INDEX(TARGET,'ECR').NE.0) THEN   !rear end cap
            TBAL =TB2(4) +.00014 +T_ECR/2.  !add last  10cm of 3He
            TAAL =TA2(4) -.00014 -T_ECR/2.
            TLEN=1.
            TTOT=1.
           ELSE
            TBAL = TB2(ISEG154)
            TAAL = TA2(ISEG154) 
            TLEN = TLEN2(ISEG154)
            TTOT =  TLEN2(1) +TLEN2(2) +TLEN2(3) +TLEN2(4) 
           ENDIF
         ELSEIF(THETAD.GT.4.0.AND.THETAD.LT.7.0) THEN
          IF(INDEX(TARGET,'ECF').NE.0) THEN  ! front end cap
            TBAL= .5*TB5(1) - .00010     !remove 7 cm of 3He
            TAAL = .5*TB5(1) +TA5(1) + .00010
            TLEN=1.
            TTOT=1.
           ELSEIF(INDEX(TARGET,'ECR').NE.0) THEN  !rear end cap
            TBAL =TB5(4) +.00007 +T_ECR/2.        !add last  5cm of 3He
            TAAL =TA5(4) -.00007 -T_ECR/2.
            TLEN=1.
            TTOT=1.
           ELSE
            TBAL = TB5(ISEG154)
            TAAL = TA5(ISEG154) 
            TLEN = TLEN5(ISEG154)
            TTOT =  TLEN5(1) +TLEN5(2) +TLEN5(3) +TLEN5(4) 
           ENDIF
         ENDIF

         RETURN
         END


                                                 
C=======================================================================
      SUBROUTINE RADLENGT2 (ITARG,X,TR,TAFTER,TBAL,TAAL)                
!_______________________________________________________________________
!  Calculates E140X radiation length (units of rl) for all materials    
!  before (TBAL) and after scattering (TAFTER,TAAL), and scattering     
!  angle TH (deg) and scattering taking place at X r.l along the axis   
!  of the target relative to the in cap of the target cell.             
!  T1G and T2G are same as T1 and T2 but in grams instead of r.l.       
!  Writes diagnostic and error messages to unit IU.                     
!                                                                       
!  Calculations are carried out in target coordinate system where the   
!  x axis is along the axis of the target and y axis perpendicular to   
!  x axis with positive values towards 8 GeV spectrometer.  The origin o
!  this coordinate system is at the in cap of the target cell.          
!  The out cap of the target cell is assumed to be ellipsoidal given by 
!  (B*(X-XC))**2 + (A*Y)**2 - (A*B)**2 = 0                              
!                                                                       
!  Changed r.l. of LH2 to 61.28 (1988 value) from 63.                   
!______________________________________________________________________ 
      COMMON /TRL/   TTARG,TWALL,TBEAM,TSPEC,KTARG,NSEG,JTARG!add 4/16/9
      REAL TTARG, TWALL, TBEAM, TSPEC                                   
      INTEGER KTARG !Dummy replacing ITARG which is in arguement        
      INTEGER JTARG ! 1 = internal rad only, 2= extern + intern         
      DIMENSION DAT(4,11)                                               
      DATA    DAT/                                                      
! TARGPARAM: radiation lengths of various types of Al from P. 52        
! Info on wire arrayS: Jerry Davis says 5 mil diamter Al wire 25 mil spa
!              and there are two arrays.                                
! Info on Hymen: John Mark says it was 1 mil Al.                        
! Info on In-cap is from p. 32 (sample in book, 3 mils thick)           
! Info on LH2: thickness is actually inner radius from 34, corrected    
!              for wall thickness and 0.4% shrinkage.                   
!              X0 value is from new 1988 particle                       
!              data book. (Used to be 63.3 in old books).               
! Values are overwritten if D2 is used                                  
! Info on cell wall. Value of 4.4 mils is from p. 24                    
!              Actual value varies with x and is hardwired into program.
! Infor on endcap comes from p. 24. Actuall thickness varies with y     
!              as is hardwired into program                             
! Info on insulation comes from J. Mark (10 layers of 1/4 mil mylar)    
!              it pretty much lines up with endcap of long targets,     
!              so for short targets need to take this into account.     
! Info on 1.6 gev window comes from p. 24 (sample)                      
! Thickness of 8 GeV exit window is 10.7 mils when not crinkled (Eislele
!                                                                       
!       Thickness     density      X0      Z    Name       Materl       
!         (cm)        (g/cm3)    (g/cm2)                                
     >    0.00400,    2.7,       24.0111, 13., ! Wire array  Al pure    
     >    0.00254,    2.7,       24.0111, 13., ! Hymen       Al pure    
     >    0.00762,    2.68,      23.6311, 13., ! In_cap      Al 5052    
     >    3.19600,    0.0707,    61.2800,  1., ! radius      LH2 or LD2 
     >    0.01143,    2.72,      23.6396, 13., ! Cell wall   Al 3004    
     >    0.01143,    2.72,      23.6396, 13., ! End_cap     Al 3004    
     >    0.00635,    1.39,      39.95,   5.2, ! Insulat     Mylar      
     >    0.00762,    2.68,      23.6311, 13., ! 1.6 window  Al 5052    
     >    0.03048,    2.68,      23.6311, 13., ! 8GEV windo  Al 5052    
     >    16.0000,    0.001205,  36.975,  7.3, ! Air /Vac gap VAC       
     >    0.02540,    2.68,      23.6311, 13./ ! quad_window Al 5052    
                                                                        
! Enter density and radiation length of liquid target                   
        IF(ITARG.EQ.1.OR.ITARG.EQ.2) THEN                               
          DAT(2,4)=0.0707                                               
          DAT(3,4)=61.28                                                
          IF(ITARG.EQ.1)  DAT(1,4)=15.836                               
          IF(ITARG.EQ.2)  DAT(1,4)=4.0470                               
        ELSE                                                            
          DAT(2,4)=0.1693                                               
          DAT(3,4)=122.60                                               
          IF(ITARG.EQ.4)  DAT(1,4)=15.763                               
          IF(ITARG.EQ.3)  DAT(1,4)=4.047                                
        ENDIF                                                           
        IF(ITARG.EQ.1) XMAX=15.818                                      
        IF(ITARG.EQ.2) XMAX=4.0292                                      
        IF(ITARG.EQ.3) XMAX=4.0292                                      
        IF(ITARG.EQ.4) XMAX=15.745                                      
                                                                        
                                                                        
        IF(ITARG.EQ.5.OR.ITARG.EQ.6) THEN                               
          TB=0.00045           ! wire array                             
     >      +0.00029           ! entrance window                        
            IF(ITARG.EQ.5) THEN ! 15 cm target                          
             IF(X.LE.0.0108)                                            
     >        TAFTER=(0.0108-X)/COS(TR)                                 
     >          +0.00011/COS(TR)                                        
             IF(X.GT.0.0108)                                            
     >        TAFTER=(0.0216-X)/COS(TR)                                 
     >        +0.00011+0.00011/COS(TR)                                  
            ELSE                                                        
              IF(TR.GT.0.671) THEN ! Doesn't hit downstream dummy       
               IF(X.LE.0.0108)                                          
     >          TAFTER=(0.0108-X)/COS(TR)                               
     >         +0.00011/COS(TR)                                         
               IF(X.GT.0.0108)                                          
     >          TAFTER=(0.0216-X)/COS(TR)                               
     >         +0.00011+0.00011/COS(TR)                                 
              ELSE                ! Does hit downstream dummy           
                TAFTER=(0.0216-X)/COS(TR) ! Target                      
     >         +0.00011+0.00011/COS(TR)                                 
              ENDIF                                                     
            ENDIF                                                       
            TA=  +0.00346      ! Chamber window                         
     >           +0.00053      ! air                                    
     >           +0.00288      ! spect window                           
          TAAL =TA                                                      
          TBAL =TB                                                      
          RETURN                                                        
          ENDIF                                                         
                                                                        
        IF(ITARG.EQ.7) THEN                                             
          TB=0.00045           ! wire array                             
     >      +0.00029           ! entrance window                        
          TAFTER=(TTARG-X)*COS(0.79)                                    
     >     /COS(ABS(0.79-TR))                                           
          TA=0.00022/SIN(TR)  !Mylar                                    
            TA=TA+0.00346      ! Chamber window                         
     >           +0.00053      ! air                                    
     >           +0.00288      ! spect window                           
          TAAL =TA                                                      
          TBAL =TB                                                      
          RETURN                                                        
          ENDIF                                                         
                                                                        
! Do liquild targets calculation                                        
                                                                        
      T1 = 0.0                                ! Clear sums              
      DO I = 1, 3                                                       
        TX = DAT(1,I)*DAT(2,I)/DAT(3,I)        ! Add wire array,        
        T1 = T1 + TX                           ! window, and in-cao     
      END DO                                                            
      TBAL=T1                                                           
      X=X*DAT(3,4)/DAT(2,4)         ! changed to unit cm                
! Find amount after scattering                                          
      TAAL=0.                                                           
! find horizontal shift for beam spot DY                                
      DY= 0.0                                                           
        DO I = 4,11                ! Do H2, wall, endcap, and mylar     
          IF(I.NE.8) THEN          ! exit window, air, and entr window  
            XEFF=DAT(1,I)          ! calculate geometry for first       
            IF(I.LE.7) CALL TARGET2(X,XMAX,TR,DY,I,ITARG,XEFF)          
            TX = XEFF*DAT(2,I)/DAT(3,I)                                 
            IF(I.EQ.4) TAFTER=TX                                        
            IF(I.GE.5) TAAL=TAAL+TX                                     
          ENDIF                                                         
        END DO                                                          
      X=X*DAT(2,4)/DAT(3,4)     ! changed back to r.l                   
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
      SUBROUTINE TARGET2 (X,XMAX,TR,DY,I,ITARG,XEFF)                    
!______________________________________________________________________ 
!  Calculates a thickness of material traversed by a scattered electron 
!  in a target contained in a target cell whose geometry is considered  
!  here.  Experiment = NE11. Only does 8 GeV side.                      
!  Input variables:                                                     
!  X    - x coordinate of scattering center,                            
!  XMAX - maximum x coordinate of LH2 = reference point (rear endcap)   
!  TR   - scattering angle alpha in rad.,                               
!  DY   - lateral beam offset in cm.  (+ = towards spectrometer)        
!  I    - indexes Ith element in the target cell.                       
!  ITARG - 1 or 4 for long cell, 2  or 3 for short cell                 
!  XEFF - thickness trversed in cm (0. if not in path)                  
!       - as input, it is nominal perpendicular thickness               
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      REAL*4  A /1.300/               !Ellipse x semiaxis (cm)          
      REAL*4  B /3.20000/             !Ellipse y semiaxis (cm)          
      REAL*4  PI/3.14159265/                                            
      REAL X,XMAX,TR,DY,XNOM,XEFF,XC,X1,X2,Y1,Y2,XWALL      
      INTEGER ITARG,I                                                   
      LOGICAL HIT                                                       
                                                                        
      XNOM=XEFF                       ! Save perpendicular thickness    
      XEFF=0.                         ! Set to 0 in case don't go thru  
      XC=XMAX-A                       ! x of ellipse center             
      IF(X .LT. 0.0 .OR. X .GT. XMAX) RETURN  ! Should never happen     
      CALL ENDCAP2(X,DY,XC,A,B,TR,X1,Y1,HIT) ! where his endcap or wal  
      IF(I.EQ.4) XEFF=SQRT((X1-X)**2+(Y1-DY)**2) ! Hydrogen             
      IF(HIT) THEN                                                      
        IF(I.EQ.6) THEN                ! Hit the endcap                 
          IF(Y1.LT.2.60) XNOM=0.0112   ! Put in y-dependence of endcap  
          IF(Y1.GE.2.60.AND.Y1.LT.2.9) XNOM=0.0178  ! based on Peter's m
          IF(Y1.GE.2.90) XNOM=0.0254   ! of an actual endcap.           
          CALL ENDCAP2(X,DY,XC,(A+XNOM),(B+XNOM),TR,X2,Y2,HIT)          
          IF(.NOT.HIT) THEN           ! Didn't hit other side of cap!   
            WRITE(6,'(1X,''ERROR IN PASSING ENDCAP'')')                 
          ELSE                                                          
            XEFF=SQRT((X2-X1)**2+(Y2-Y1)**2)                            
          ENDIF                                                         
        ENDIF                         ! If hit endcap, never hits wall b
        IF(I.EQ.7.AND.(ITARG.EQ.2.OR.ITARG.EQ.3)) THEN ! short target, c
          IF(X+(3.2-DY)*COS(TR)/SIN(TR).LT.9.5) THEN ! thru insultation 
            XEFF=XNOM/SIN(TR)  ! which ends at endcap of long target    
          ENDIF                ! (9.5 cm from begining of short target) 
        ENDIF                  ! where radius is 3.2 cm TARGPARAM:      
      ELSE                            ! Do wall and inslut.             
        IF(I.EQ.5) THEN        ! Put in rough x dependence of thickness 
          XWALL=XMAX-1.3-X     ! TARGPARAM: 1.3 cm is where wall begins 
          XWALL=XWALL-(3.2-DY)*COS(TR)/SIN(TR) ! Go from middle to edg  
          IF(XWALL.LT.0.50) XNOM=0.0254   ! based on Peter's meas.      
          IF(XWALL.GE.0.50.AND.XWALL.LT.1.00) XNOM=0.0238               
          IF(XWALL.GE.1.00.AND.XWALL.LT.1.50) XNOM=0.0152               
          IF(XWALL.GE.1.50.AND.XWALL.LT.2.50) XNOM=0.0127               
          IF(XWALL.GT.2.50) XNOM=0.0114                                 
        ENDIF                                                           
        IF(I.EQ.5.OR.I.EQ.7) XEFF=XNOM/SIN(TR)  ! Goes at angle TR      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
      SUBROUTINE  ENDCAP2 (X0,Y0,XC,A,B,TR,X,Y,HIT)                     
!______________________________________________________________________ 
!  Calculates an intersection point between a straight line and an      
!  ellipse (end cap of a target cell).  The line is given by:           
!  y = Y0 + (x - X0)*tan(TR),  and the ellipse by:                      
!  (B*(x - XC))**2 + (A*y)**2 = (A*B)**2.                               
!  The center of the target is taken as the y axis.                     
!  It set HIT .TRUE. if the solution exists, .FALSE. otherwise.         
!  If doesn't hit endcap, then give x,y of where hits target wall,      
!  assumed to be at radius B                                            
!  Output variables:                                                    
!  X,Y  - are the coordinates of the intersection.                      
!                                                                       
!  Assumes all x dimensions measured from in-cap of target, and         
!  positive Y is towards the 8 gev spectrometer (where electron goes).  
!  Note that the scattered electron can exit ONLY through that portion o
!  the ellipse for which XC<=x<=XC+A and y=>DY                          
!                                                                       
!  ZMS                  Jan. 1985...                                    
!  AFS                  2 May 85  simplified routine.                   
!  pyb  checked routine, added implicit none, provision for 90 deg      
!       and correct treatment of Y0                                     
!______________________________________________________________________ 
      IMPLICIT NONE                                                     
      REAL X0,Y0,XC,A,B,TR,X,Y                                          
      REAL TA,G,H,AA,BB,CC,DD                                           
      LOGICAL HIT                                                       
      REAL*4  PI/ 3.14159265/                                           
                                                                        
      HIT = .FALSE.                                                     
      X = 0.                                                            
      Y = 0.                                                            
                                                                        
! If TR very close to 90 degress, do special solution                   
      IF(ABS(TR-PI/2.).LT.0.017) THEN                                   
        IF(X0.GT.XC) THEN      ! Is it beyond center of endcap          
          HIT=.TRUE.                                                    
          X=X0                                                          
          Y=B*SQRT(1.0-(X0-XC)**2/A**2)                                 
        ELSE                   ! Must be upstream of endcap             
          X=X0                                                          
          Y=B                                                           
        ENDIF                                                           
        RETURN                                                          
      ENDIF                                                             
                                                                        
! Solve quadratic equation to see if hits endcap                        
      TA = TAN (TR)                                                     
      G  = A*TA                                                         
      H  = (X0 - XC)*TA - Y0                                            
      AA = 1.0 + G*G/(B*B)                                              
      BB = 2.0 * H                                                      
      CC = H*H - G*G                                                    
      DD = BB*BB - 4.0*AA*CC                                            
      IF(DD .GE. 0.0) THEN            ! See if radical postive          
        Y = 0.5*(-BB + SQRT (DD))/AA                                    
        X = X0 + (Y-Y0)/TA                                              
      END IF                                                            
      IF(X .GE. XC) THEN   ! Is solution in endcap region?              
        HIT = .TRUE.                                                    
      ELSE                 ! Otherwise, make it hit the wall            
        X=X0+(B-Y0)/TA                                                  
        Y=B                                                             
      ENDIF                                                             
      RETURN                                                            
      END                                                               
C==============================================================         

      subroutine eg4radlen(theta,TBEFOR,TAFTER)
! materials before, after target excluding anything
! between Banjo windows. From Piotr Konczykowski, Mon, 10 Sep 2007 
! theta in degrees, tbefor and tafter are in r.l.
      implicit none

      real He,theta,th
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
      real TBEFOR,TAFTER

      th=abs(theta)

      TBEFOR=0.00167415733

         p1=9.037982e-03
         p2=1.624980e-08
         p3=5.874948e-03
         p4=9.681330e-04
         p5=1.040954e-02
         p6=3.256229e-04
         p7=1.276336e-02
         p8=1.231782e-04
         p9=1.362711e-02
         p10=1.022152e-04
         p11=-3.142974e-06
         p12=7.551468e-08

      if((th.ge.0.).and.(th.lt.3.25)) then
         TAFTER=p1+p2*th
      elseif((th.ge.3.25).and.(th.lt.7.)) then  
         TAFTER=p3+p4*th
      elseif((th.ge.7.).and.(th.lt.11.9)) then  
         TAFTER=p5+p6*th
      elseif((th.ge.11.9).and.(th.lt.16.)) then  
         TAFTER=p7+p8*th
      elseif((th.ge.16.)) then  
         TAFTER=p9+p10*th+p11*th*th+p12*th*th*th
      endif   

      end
