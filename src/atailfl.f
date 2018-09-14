!-----------------------------------------------------------------------

      REAL*8 FUNCTION ATAILFL_139(ES,EP,SNSQ,TB,TA,ZP,ZN,A,IFLG)  
      IMPLICIT NONE
      REAL*8 ES,EP,SNSQ,TB,TA,ZP,ZN,A
      INTEGER IFLG
      REAL*8 FUNXL_139,elastcl_139,PHIS,PHIP,BTR,PEAKS,PEAKP,COEFF,
     >       VAC_TERM_S,VAC_TERM_P,FACP,FACS,PART1,PART2,PART3,ABREM,
     >       ATAIL,COST,QUADMO,TM,Q2,XI,OMP,OMS,Q2S,Q2P,VS,VP,ETA,DLZ,
     >       B, VERTEX_TERM_S, VERTEX_TERM_P,TOT
      INTEGER NLVL
      REAL*4 BREMS_139,DELVAC
      real*8 DDILOG
      EXTERNAL FUNXL_139
      COMMON/TAL1/ESX,EPX,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1  
      REAL*8  ESX,EPX,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1  
      COMMON/TAL2/ZPX,ZNX,AX,IFLGX
      REAL*8 ZPX,ZNX,AX
      INTEGER IFLGX
      REAL*8 ETAIL,QTAIL,QE_SMEAR_PEAK,QE_DELTA_EXACT,QE_DELTA_PEAK, 
     >       TARGET_RADIATION,QE_SOFT                                           
      COMMON/WHICHT/ ETAIL, QTAIL, QE_SMEAR_PEAK, QE_DELTA_EXACT,               
     >               QE_DELTA_PEAK,TARGET_RADIATION,QE_SOFT                     
      REAL*8 EMSQ/.26113D-6/,PI/3.1415926535D0/                           
      REAL*8 ALPHA/7.2973515D-3/,TWO_ALPH_PI/.00464564D0/                                            
      common/testing/ prttst,usegrd
      logical prttst,usegrd

C------------------------------------------------------------------------
C     CALCULATES THE RADIATIVE TAIL OF ELASTIC OR QUASI-ELASTIC         
C     SCATTERING USING THE EXACT EXPRESSION (TSAI)                      
C                                                                       
C     GIVEN ES = INCIDENT ENERGY IN GEV                                 
C           EP = SCATTERED ENERGY IN GEV                                
C           SNSQ OF THE SCATTERING ANGLE                                
C           TB = RADIATOR "BEFORE" IN RADIATION LENGTHS                 
C           TA = RADIATOR "AFTER"  IN RADIATION LENGTHS                 
C           ZP = ATOMIC NUMBER (NO. OF PROTONS)                         
C           ZN = AVERAGE NUMBER OF NEUTRONS                             
C           A= ATOMIC WEIGHT (CARBON-12 SCALE) H=1.00797                
C           IFLG = 0 FOR ELASTIC, 1 FOR QUASI-ELASTIC                   
C                                                                       
C     ANSWER IS IN PB/(GEV-SR)                                          
! Version of ATAILFL, originally from E139/E140 code, 
! adopted for use in EXTERNAL 7/25/96
!------------------------------------------------------------------------------
      ATAILFL_139=0.D0                                                      
      IF (EP.LE.0.D0) RETURN                                            
C                                                                       
C     START FILLING THE COMMON BLOCK FOR THE INTEGRAND                  
      EPX=EP                                                            
      ESX=ES                                                            
      ZPX=ZP                                                            
      ZNX=ZN                                                            
      AX=A                                                              
      IFLGX=IFLG                                                        
C                                                                       
C     CALCULATE TARGET DEPENDENT PARAMETERS                             
C                                                                       
C *** NOTE FACTORS OF 1440 & 183 HAVE BEEN CHANGED TO 1194 & 184.15     
C     TO AGREE WITH TSAI, REV. MOD. PHYS. 46, 4 (1975), ON P. 828       
C                                              25 JUNE 83 AFS           
C                                                                       
      DLZ=DLOG(184.15D0/ZP**(1.D0/3.D0))                                
      ETA=DLOG(1194.D0/ZP**(2.D0/3.D0))/DLZ                             
      B=4.D0/3.D0*(1.D0+(ZP+1.D0)/(9.D0*(ZP+ETA)*DLZ))                  
      XI=.11D0*(TB+TA)/((ZP+ETA)*DLZ)                                   
      TM=.93828D0                                                      
      IF (IFLG.EQ.0) TM=A*TM/1.00797D0                                  
      TMSQ=TM**2                                                        
C                                                                       
C     CALCULATE KINEMATICS                                              
      Q2=4.D0*ES*EP*SNSQ                                                
      OMS=ES-EP/(1.D0-2.D0*EP*SNSQ/TM)                                  
      OMP=ES/(1.D0+2.D0*ES*SNSQ/TM)-EP                                  
      IF (OMP.LE.0.) RETURN                                             
      Q2S=4.D0*(ES-OMS)*EP*SNSQ                                         
      Q2P=4.D0*ES*(EP+OMP)*SNSQ                                         
      VS=OMS/ES                                                         
      VP=OMP/(EP+OMP)                                                   
C                                                                       
C     CALCULATE QED TERMS AND MULT. PHOTON NORMALIZATION                
      TOT=TB+TA                                                         
C Second order term in B*TOT added 5/28/87 5:07 PM- Steve Rock          
      FAC1=(1.D0+0.5772D0*B*TOT- .66*(B*TOT)**2)                        
     >   -ALPHA/(2.D0*PI)*DLOG(ES/EP)**2                                
Cxxx dropped this term temporarily to see if can get code to work
     *+ALPHA/PI*(PI**2/6.D0-DDILOG(1.D0-SNSQ))                           
C     FACS=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2S/EMSQ))  
C     FACP=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2P/EMSQ))  
                                                                        
                                                                        
                                                                        
C*******8/19/87**  Steve Rock ***************************************   
C Use the Bardin Calculation of Vertex terms instead of only electron   
C  and Tsai's Vertex term.                                              
C     FAC2=2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(-QSQ/EMSQ))      
                                                                        
C  Vacuum terms of Bardin including 3 leptons and quarks                
      VAC_TERM_S= 2. * DELVAC(SNGL(Q2S)) !DELVAC was implicitly real*8 before
      VAC_TERM_P= 2. * DELVAC(SNGL(Q2P))                                
                                                                        
C  Vertex term of Tsai (eq. 2.11)                                       
      VERTEX_TERM_S = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG(Q2S/EMSQ) )      
      VERTEX_TERM_P = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG(Q2P/EMSQ) )      
C7/29/96      FAC2 = VAC_TERM + VERTEX_TERM                                     
C     FACS=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2S/EMSQ))  
C     FACP=FAC1+2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(Q2P/EMSQ))  
      FACS = FAC1 + VAC_TERM_S + VERTEX_TERM_S                          
      FACP = FAC1 + VAC_TERM_P + VERTEX_TERM_P                          
                                                                        
C*******************************************************************    
                                                                        
                                                                        
C     EQUIVALENT RADIATOR FOR INTERNAL                                  
      BTR=ALPHA/PI*(DLOG(Q2/EMSQ)-1.D0)                                 
C                                                                       
C     BREMSSTRAHLUNG SPECTRUM                                           
C$$   PHIS=1.D0-VS+0.75D0*VS**2                                         
C$$   PHIP=1.D0-VP+0.75D0*VP**2                                         
                                                                        
C Installed 6/4/87 to get better value of bremstrung spectrum using     
C  Tsai's ReV Mod Physics Article eq (3.83) Steve                       
C     This is the Old-E139.E140 version of BREMS which is different than
c     the one traditionally in External.  SER-7/25/96
      PHIS= BREMS_139(VS,ZP)                                                
      PHIP= BREMS_139(VP,ZP)                                                
C                                                                       
C     S-PEAK                                                            
      PEAKS=(TM+2.D0*(ES-OMS)*SNSQ)/(TM-2.D0*EP*SNSQ)*FACS*             
     * ELASTCL_139(ES-OMS,SNSQ,ZP,ZN,A,IFLG)*(B*TB/OMS*PHIS                 
     * +XI/(2.D0*DMAX1(OMS,10.D0*XI)**2))                               
C                                                                       
C     P-PEAK                                                            
      PEAKP=FACP*                                                       
     * ELASTCL_139(ES,SNSQ,ZP,ZN,A,IFLG)*(B*TA/OMP*PHIP                     
     * +XI/(2.D0*DMAX1(OMP,10.D0*XI)**2))                               
C                                                                       
C     CALCULATE QUANTITIES INDEPENDENT OF PHOTON ANGLE                  
      SVEC=DSQRT(ES**2-EMSQ)                                            
      PVEC=DSQRT(EP**2-EMSQ)                                            
      COST=1.D0-2.D0*SNSQ                                               
      SP=ES*EP-SVEC*PVEC*COST                                           
      USQ=2.D0*EMSQ+TMSQ+2.D0*TM*(ES-EP)-2.D0*SP                        
      U0=ES+TM-EP                                                       
      UVEC=DSQRT(U0**2-USQ)                                             
      COSP=(SVEC*COST-PVEC)/UVEC                                        
      COSS=(SVEC-PVEC*COST)/UVEC                                        
C                                                                       
C     INTEGRATE OVER COS(THETAK) IN THREE PIECES                        
C     MAKING CERTAIN THAT THE S AND P PEAKS ARE FOUND                   
                                                                        
      PART1=QUADMO(FUNXL_139,-1.D0,COSP,1.D-4,NLVL)                         
      PART2=QUADMO(FUNXL_139,COSP,COSS,1.D-4,NLVL)                          
      PART3=QUADMO(FUNXL_139,COSS,1.D0,1.D-4,NLVL)                          
                                                                        
C$$   PART1=DCADRE(FUNXL,-1.D0,COSP,0.D0,1.D-4,ERROR,IER)               
C$$   IF((IER.NE.0).AND.(IER.NE.65).AND.(IER.NE.66))                    
C$$  > WRITE(6,'('' ** DCADRE INTEGRATION:PART1 in ATAILFL; IER='',I4)')
C$$  >  IER                                                             
C$$   PART2=DCADRE(FUNXL,COSP,COSS,0.D0,1.D-4,ERROR,IER)                
C$$   IF((IER.NE.0).AND.(IER.NE.65).AND.(IER.NE.66))                    
C$$  > WRITE(6,'('' ** DCADRE INTEGRATION:PART2 in ATAILFL; IER='',I4)')
C$$  >  IER                                                             
C$$   PART3=DCADRE(FUNXL,COSS,1.D0,0.D0,1.D-4,ERROR,IER)                
C$$   IF((IER.NE.0).AND.(IER.NE.65).AND.(IER.NE.66))                    
C$$  > WRITE(6,'('' ** DCADRE INTEGRATION:PART3 in ATAILFL; IER='',I4)')
C$$  >  IER                                                             
                                                                        
                                                                        
      COEFF=389.44D6*ALPHA**3*TM/PI*EP/ES                               
      ATAIL=COEFF*(PART1+PART2+PART3)                                   
C                                                                       
C     PUT IT ALL TOGETHER                                               
      ABREM=PEAKS+PEAKP                                                 
      QE_SOFT=VS**(B*TB+BTR)*VP**(B*TA+BTR)                             
                                                                        
C     TARGET_RADIATION=TARCOR(ES,EP,SNSQ,ZP,A,IFLG)                     
      TARGET_RADIATION= 1.                                              
                                                                        
      ATAILFL_139=(ATAIL*TARGET_RADIATION+ABREM)*QE_SOFT 
      if(prttst) then
        write(*,111) tb,ta,COEFF,PART1,PART2,PART3
 111    format(/1x,' tb,ta,co,p1,p2,p3,at=',2f7.4,5e10.2)
        write(*,112) iflg,atail,ATAILFL_139,abrem,qe_soft,peaks,peakp
 112    format(1x,' iflg,at,ab,qe,ps,pb=',i2,6e11.3)
      endif
C                                                                       
C DEBUG:                                                                
C                                                                       
c      WRITE(*,1000) COEFF,PART1,PART2,PART3,ATAIL,ABREM,QE_SOFT  
c1000  FORMAT(' COEFF,PART1,PART2,PART3,ATAIL,ABREM,QE_SOFT:',    
c     >       / 7(1PE14.5) )                                             
c      WRITE(*,2000) ATAILFL_139                                             
2000  FORMAT(' ATAILFL =',1PE14.5)                                      
C                                                                       
      RETURN                                                            
      END                                                               
!----------------------------------------------------------------------------



      REAL*4 FUNCTION ATAILFL1(T)
!-----------------------------------------------------------------------------------
! Calls the E139/E140 code function ATAILFL_139 to calculate the 
!  "exact" nuclear tail.
!   SER 7/24/96
!----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL T
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL  E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      REAL           ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
! TBEFOR, TBAL are the target material and windows before scattering point
! TAFTER, TAAL are the target material and windows after scattering point
! NOTE: TB and TA include the equivilent radiator: This is added on in ATAILFL
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL*8 TB1,TA1, AtWt,ATAILFL_139
      common/testing/ prttst,usegrd
      logical prttst,usegrd
      logical smrelasp,usearen
      real smrdEp
      common/smrproton/ smrelasp,smrdEp,usearen
                                                                        
      ATAILFL1 = 0.     
c removed 3/25/08
c      if(iz.lt.1) return
cxx      if(ia.eq.1.and.smrelasp) return

      AtWt= avgM/.93828 * 1.00797  ! atomic weight so that H=1.00797

      CALL RADIATORS(T)
      TB1 = TBEFOR +TBAL
      TA1 = TAFTER +TAAL
      ATAILFL1 =ATAILFL_139(DBLE(E0),DBLE(EP),DBLE(SIN2),TB1,TA1,
c 3/25/08 modifiled to ensure z>0
     >   DBLE(max(1,iZ)),DBLE(avgN),AtWt,0) !last arguement=0 for elastic
      if(prttst) write(*,111) t,tb1,ta1,atailfl1
 111  format(1x,'t,tb1,ta1,atailfl1=',3f8.5,e11.2) 
      RETURN
      END  
!==============================================================



      REAL*4 FUNCTION ATAILFL_QE(T)
!-----------------------------------------------------------------------------------
! Calls the E139/E140 code function ATAILFL_139 to calculate the 
!  quasi-elastic tail integrating over COSK.
!   SER 7/31/96
!----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL T
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL  E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      REAL           ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 

! TBEFOR, TBAL are the target material and windows before scattering point
! TAFTER, TAAL are the target material and windows after scattering point
! NOTE: TB and TA include the equivilent radiator: This is added on in ATAILFL
      REAL*8 TB1,TA1,ATAILFL_139

      CALL RADIATORS(T)
      TB1 = TBEFOR +TBAL
      TA1 = TAFTER +TAAL
      ATAILFL_QE = ATAILFL_139(DBLE(E0),DBLE(EP),DBLE(SIN2),TB1,TA1,
     >  DBLE(iZ),DBLE(avgN),DBLE(AvgA),1) !last arguement=1 for Quasi Elastic
      RETURN
      END  
