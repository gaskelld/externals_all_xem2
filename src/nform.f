      SUBROUTINE NFORM(IG,QQG,GEP,GEN,GMP,GMN)                             

!-----------------------------------------------------------------------
C CALCULATE NUCLEON FORM FACTORS                                        
! Modified by Steve Rock on 6/21/96 adding Peters's Phys Rev C. fit
!  and putting IG in arguements
!
C IG =1 - DIPOLE WITH FORM FACTOR SCALING AND GEN=0.0                
C     2 - IJL FIVE PARAMETER MODEL FIT                               
C     3 - GEP AND GMP FROM IJL, GMN=DIPOLE, GEN=GALSTER              
C            WE CALL THIS THE "BEST FIT" NUCLEON FORM FACTORS           
C     4 - BEST FIT EXCEPT GEN = 0.0                                  
C     5 - BLATNIK + ZOVKO VDM FIT                                    
C     6 - JANNSENS 1966 STANDARD FIT                                 
C     7 - DIPOLE + F1N = 0.0                                         
C     8 - GES = 0.0,  GMS = 1.0                                      
C     9 - GES = 1.0,  GMS = 0.0                                      
C    10 - HOHLER1 - PROTON AND NEUTRON FROM FIT 8.2                  
C    11 - HOHLER2 - PROTON FROM FIT 5.3, NEUTRON FROM 8.2            
C    12 - GARI + KRUMPELMANN, Z PHYS. A322,689(1985)                 
C    13 - KORNER + KURODA, PHYS. REV. D16,2165(1977)  
C    14 - GARI AND KRUMPELMANN WITH NE11 FIT (PETER')                  
C     15 - Peter Bosted's fit from SLAC-PUB-6651 (NE11 data + other) in
c     Phys Rev C
C    16 - Radyushkin,  Acta Physica Polonica B15 403,
c     (1984) 
c    
C QQG = INPUT Q SQUARED (GEV**2)  
!---------------------------------------------------------------------------


      IMPLICIT NONE
      INTEGER IG
      REAL*8 QQG,GEP,GEN,GMP,GMN,QQ,TAU
      REAL*8 GT,T1,T2,ALPH,TOP,BOT,RHO,F1S,F1V,F2S,F2V
      REAL*8 RS,RV,F1E,F2E,F1M,F2M
      REAL*8 F1,F2,F3,GES,GMS,GEV,GMV
      REAL*8 F1RHO,F2RHO,F1P,F2P
      REAL*8 QQP,C1,C2,C3,C4,F2VK,F2SK,F1N,F2N
      REAL*8 Q,Q3,Q4
      REAL*8 FRAC1,FRAC2,GD
      INTEGER I,IN

C IJL PARAMETERS FOR 5 PARAMETER DIPOLE FIT (IN GEV UNITS)  
C PHYS LETT. 43B, 191(1973)              
      REAL*8   GAM, BR, BW, BF, AF
      REAL*8   RMN2, RMW2, RMF2, RMR2, GAMR, PI   
      REAL*8   RMPI/ .139/,  RMPI2/ .019321  /                         
                                                                        
C PARAMETERS FOR BLATNIK AND ZOVKO VDM FIT                              
C ACTA PHYSICA AUSTRIACA 39, 62(1974)                                   
C VECTOR MESON MASSES SQUARED (GEV UNITS)  
      REAL*8     TRO, TROP, TROPP, TFI, TOM, TOMP     
                                                                        
C FITTED PARAMETERS        
      REAL*8     RMUS, RMUV, BS, BV         
                                                                        
      REAL*8     RMRHO2, CRHO, RKRHO, RKS, RMOMG2, COMG, RKOMG, RKV     
      REAL*8     RLAM12, RLAM22, RLAMQ2
                                                                        
C PARAMETERS FOR KORNER AND KURODA                                      
C VECTOR MESON MASSES SQUARED USING REGGE PARAMETER ALPHA=1.    
      REAL*8     VRH1, VRH2, VRH3, VOM1, VOM2, VOM3          
      REAL*8     RMUP/ 2.792782/ ,RMUN/ -1.913148 /                      
      common/testing/ prttst,usegrd
      logical prttst,usegrd

! Parameters for Radyuskin: From Table on p 414 of paper.               
      REAL R_QSQ(12)   / 1.,   2.,   3.,   4.,   5.,   6.,              
     >                   8.,  10.,  12.,  15.,  20.,  30./              
      REAL R_GMP_U_D(12)/ .91, 1.01, 1.05, 1.05, 1.04, 1.02,            
     >                   0.97, 0.91, 0.86, 0.78, 0.67, 0.53/            
      REAL R_GEP_D(12)  /1.00, 1.13, 1.16, 1.15, 1.11, 1.06,            
     >                   0.95, 0.86, 0.77, 0.67, 0.54, 0.38/            
      REAL R_GMN_U_D(12)/0.82, 0.80, 0.79, 0.77, 0.76, 0.74,            
     >                   0.70, 0.65, 0.61, 0.56, 0.49, 0.38/            
      REAL R_GEN_D(12)  /-.13, -.12, -.10, -.06, -.03, 0.00,            
     >                   0.05, 0.08, 0.11, 0.13, 0.14, 0.14/    
      DATA     RMN2, RMW2, RMF2, RMR2, GAMR, PI                         
     *                / 0.8817, .6146, 1.0384, 0.5852 , .112 , 3.14159/
      DATA     GAM, BR, BW, BF, AF  /0.25, 0.672, 1.102, 0.112, -0.052/
      DATA     TRO, TROP, TROPP, TFI, TOM, TOMP                         
     *                       / 0.585, 1.30,  2.10,  1.039, 0.614, 1.40 /
      DATA     RMUS, RMUV, BS, BV  / -0.060, 1.853, -0.91, -1.10 /      
      DATA     RMRHO2, CRHO, RKRHO, RKS, RMOMG2, COMG, RKOMG, RKV       
     *         /0.6022, 0.377, 6.62, -0.12, 0.6147, 0.411, 0.163, 3.706/
      DATA     RLAM12, RLAM22, RLAMQ2  /  0.632, 5.153, 0.0841 /        
      DATA       VRH1, VRH2, VRH3, VOM1, VOM2, VOM3                       
     *             / 0.593, 1.593, 2.593, 0.614, 1.614, 2.614 /
C=======================================================================
                                                                        
      QQ  = QQG/.197328**2                                              
      TAU = QQG/(4.*RMN2)                                               
      GO TO (110,120,120,120,150,160,170,180,190,200,200,               
     > 220,230,240,250,260) IG                                             
C DIPOLE                                                                
  110 GEP = 1./(1.+QQG/0.71)**2                                         
      GEN = 0.0                                                         
      GMP = RMUP*GEP                                                    
      GMN = RMUN*GEP                                                    
      RETURN                                                            
                                                                        
C IJL 5 PARAMTER JOB 
  120 GT  = 0.5/(1.+GAM*QQG)**2                                         
      T1  = SQRT(QQG+4.*RMPI2)                                          
      T2  = SQRT(QQG)                                                   
      ALPH= 2.*T1*LOG((T1+T2)/(2.*RMPI))/(T2*PI)                       
      TOP = RMR2+8.*GAMR*RMPI/PI                                        
      BOT = RMR2+QQG+(4.*RMPI2+QQG)*GAMR*ALPH/RMPI                      
      RHO = TOP/BOT                                                     
      F1S = GT*((1.-BW-BF)+BW/(1.+QQG/RMW2)+BF/(1.+QQG/RMF2))           
      F1V = GT*((1.-BR)+BR*RHO)                                         
      F2S = GT*((-0.12-AF)/(1.+QQG/RMW2)+AF/(1.+QQG/RMF2))              
      F2V = GT*(3.706*RHO)                                              
      GEP = F1V+F1S-TAU*(F2V+F2S)                                       
      GEN = F1S-F1V-TAU*(F2S-F2V)                                       
      GMP = F1V+F1S+F2V+F2S                                             
      GMN = F1S-F1V+F2S-F2V                                             
      IF (IG.EQ.2) RETURN                                               
      GD  = 1./(1.+QQG/.71)**2                                          
      GMN = RMUN*GD                                                     
      GEN = -RMUN*TAU*GD/(1.+5.6*TAU)                                   
      IF (IG.EQ.3) RETURN                                               
      GEN = 0.0                                                         
      RETURN                                                            
                                                                        
C BLATNIK AND ZOVKO                                                     
  150 RS  = 1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))                
      RV  = 1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP))              
      F1E = (0.5-TAU*(RMUS+2.*RMN2*BS))*RS                              
      F2E = (0.5-TAU*(RMUV+2.*RMN2*BV))*RV                              
      F1M = (0.5+RMUS-0.5*BS*QQG)*RS                                    
      F2M = (0.5+RMUV-0.5*BV*QQG)*RV                                    
      GEP = F1E+F2E                                                     
      GMP = F1M+F2M                                                     
      GEN = F1E-F2E                                                     
      GMN = F1M-F2M                                                     
      RETURN                                                            
                                                                        
C JANNSSENS                                                             
  160 F1  = 1.+QQ/15.7                                                  
      F2  = 1.+QQ/26.7                                                  
      F3  = 1.+QQ/8.19                                                  
      GES = 0.5  *(2.5 /F1-1.6 /F2+0.10)                                
      GMS = 0.44 *(3.33/F1-2.77/F2+0.44)                                
      GEV = 0.5  *(1.16/F3-0.16)                                        
      GMV = 2.353*(1.11/F3-0.11)                                        
      GEP = GES+GEV                                                     
      GMP = GMS+GMV                                                     
      GEN = GES-GEV                                                     
      GMN = GMS-GMV                                                     
      RETURN                                                            
                                                                        
C DIPOLE + F1N = 0.0                                                    
  170 GEP = 1./(1.+QQG/0.71)**2                                         
      GEN = -RMUN*TAU*GEP                                               
      GMP =  RMUP*GEP                                                   
      GMN =  RMUN*GEP                                                   
      RETURN                                                            
                                                                        
  180 GEP = 0.0                                                         
      GEN = 0.0                                                         
      GMP = 1.0                                                         
      GMN = 0.0                                                         
      RETURN                                                            
                                                                        
  190 GEP = 1.0                                                         
      GEN = 0.0                                                         
      GMP = 0.0                                                         
      GMN = 0.0                                                         
      RETURN                                                            
                                                                        
C HOHLER1 AND HOHLER2                                                   
  200 F1RHO = 0.5*(0.955+0.090/(1.+QQG/0.355)**2)/(1.+QQG/0.536)        
      F2RHO = 0.5*(5.335+0.962/(1.+QQG/0.268))   /(1.+QQG/0.603)        
      F1S   =  0.71/(0.6129+QQG)-0.64/(1.0404+QQG)-0.13/(3.240+QQG)     
      F2S   = -0.11/(0.6129+QQG)+0.13/(1.0404+QQG)-0.02/(3.240+QQG)     
      F1V   = F1RHO+0.05/(1.464+QQG)-0.52/(6.0025+QQG)+0.28/(8.7025+QQG)
      F2V   = F2RHO-1.99/(1.464+QQG)+0.20/(6.0025+QQG)+0.19/(8.7025+QQG)
      GEP = F1V+F1S-TAU*(F2V+F2S)                                       
      GEN = F1S-F1V-TAU*(F2S-F2V)                                       
      GMP = F1V+F1S+F2V+F2S                                             
      GMN = F1S-F1V+F2S-F2V                                             
      IF (IG.EQ.10) RETURN                                              
                                                                        
C HOHLER2 - USE PROTON FIT 5.3                                          
      F1P = F1RHO+0.67/(0.6129+QQG)-0.39/(0.9216+QQG)-0.54/( 2.7556+QQG)
      F2P = F2RHO+0.04/(0.6129+QQG)-1.88/(1.2996+QQG)+0.24/(10.1761+QQG)
      GEP = F1P-TAU*F2P                                                 
      GMP = F1P+F2P                                                     
      RETURN                                                            
                                                                        
C GARI AND KRUMPELMANN                                                  
 220  QQP  = QQG*LOG(((RLAM22+QQG)/RLAMQ2))/LOG(RLAM22/RLAMQ2)          
      C1   = RLAM12/(RLAM12+QQP)                                        
      C2   = RLAM22/(RLAM22+QQP)                                        
      F1   = C1*C2                                                      
      F2   = F1*C2                                                      
      C3   = RMRHO2/(RMRHO2+QQG)                                        
      C4   = RMOMG2/(RMOMG2+QQG)                                        
      F1V  = (C3*CRHO+(1-CRHO))*F1                                      
      F1S  = (C4*COMG+(1-COMG))*F1                                      
      F2VK = (C3*CRHO*RKRHO+(RKV-CRHO*RKRHO))*F2                        
      F2SK = (C4*COMG*RKOMG+(RKS-COMG*RKOMG))*F2                        
      F1P  = 0.5*( F1S+F1V)                                             
      F1N  = 0.5*( F1S-F1V)                                             
      F2P  = 0.5*(F2SK+F2VK)                                            
      F2N  = 0.5*(F2SK-F2VK)                                            
      GEP  = F1P-TAU*F2P                                                
      GMP  = F1P+F2P                                                    
      GEN  = F1N-TAU*F2N                                                
      GMN  = F1N+F2N                                                    
      RETURN                                                            
                                                                        
C KORNER AND KURODA                                                     
  230 F1S = (1/(1+QQG/VOM1))*(1/(1+QQG/VOM2))                           
      F1V = (1/(1+QQG/VRH1))*(1/(1+QQG/VRH2))                           
      F2S = F1S*(1/(1+QQG/VOM3))                                        
      F2V = F1V*(1/(1+QQG/VRH3))                                        
      F1P = 0.5*F1S+0.5*F1V                                             
      F1N = 0.5*F1S-0.5*F1V                                             
      F2P = (RMUP-1)*(-0.0335*F2S+1.0335*F2V)                           
      F2N =    -RMUN*(-0.0335*F2S-1.0335*F2V)                           
      GEP = F1P-TAU*F2P                                                 
      GMP = F1P+F2P                                                     
      GEN = F1N-TAU*F2N                                                 
      GMN = F1N+F2N                                                     
      RETURN                                                            
                                                                        
C GARI AND KRUMPELMANN WITH NE11 FIT (PETER')                           
 240  QQP  = QQG*LOG(((RLAM22+QQG)/RLAMQ2))/LOG(RLAM22/RLAMQ2)          
      C1   = RLAM12/(RLAM12+QQP)                                        
      C2   = RLAM22/(RLAM22+QQP)                                        
      F1   = C1*C2                                                      
      F2   = F1*C2                                                      
      C3   = RMRHO2/(RMRHO2+QQG)                                        
      C4   = RMOMG2/(RMOMG2+QQG)                                        
      F1V  = (C3*CRHO+(1-CRHO))*F1                                      
      F1S  = (C4*COMG+(1-COMG))*F1                                      
      F2VK = (C3*CRHO*RKRHO+(RKV-CRHO*RKRHO))*F2                        
      F2SK = (C4*COMG*RKOMG+(RKS-COMG*RKOMG))*F2                        
      F1P  = 0.5*( F1S+F1V)                                             
      F1N  = 0.5*( F1S-F1V)                                             
      F2P  = 0.5*(F2SK+F2VK)                                            
      F2N  = 0.5*(F2SK-F2VK)                                            
      GMP  = F1P+F2P                                                    
      GEP  = GMP/RMUP                                                   
      GEN  = 0.0                                                        
      GMN  = GMP/RMUP * RMUN                                            
      RETURN                                                            

! Peter Bosted's fit from SLAC-PUB-6651 (NE11 data + other) in Phys Rev C
 250  CONTINUE
      Q = SQRT(QQG)
      Q3= Q*QQG
      Q4 = QQG*QQG
      TAU = QQG/3.52
      GEP = 1./  (1.+0.14*Q +3.01*QQG + 0.02*Q3 +1.20*Q4 +.32*Q**5)
      GMP = RMUP*GEP
      GMN = RMUN/(1.-1.74*Q +9.29*QQG - 7.63*Q3 +4.63*Q4)
      GEN = 1.25* RMUN*TAU/(1.+18.3*TAU)/((1.+QQG/.71)**2)
c     if(prttst) write(8,'(1x,11f7.3)') qqg,q,q3,q4,gep,gmp,gmn,gen
c !!! this line was missing up until 8/06!!!!
      return

 260  CONTINUE
! Radyushkin:                                                           
      GD = 1./(1.+ QQG/.71)**2                                          
      IF(QQG.LT.R_QSQ(1)) QQG=1.                                        
      IF(QQG.GT.R_QSQ(12))QQG=12.                                       
      DO I=1,11                                                         
       IF(QQG.GE.R_QSQ(I) .AND. QQG.LE.R_QSQ(I+1) ) THEN                
        IN = I                                                          
        GO TO 241                                                       
       ENDIF                                                            
      ENDDO                                                             
! Out of range.                                                         
      GMP=0.                                                            
      GMN=0.                                                            
      GEP=0.                                                            
      GEN=0.                                                            
      RETURN                                                      
! Do linear interpolation                                               
  241 FRAC1 = (QQG - R_QSQ(IN) )/(R_QSQ(IN+1) -R_QSQ(IN) )              
      FRAC2 = (R_QSQ(IN+1) -QQG)/(R_QSQ(IN+1) -R_QSQ(IN) )              
      GMP = RMUP*GD* (R_GMP_U_D(IN) * FRAC2 + R_GMP_U_D(IN+1) *FRAC1)   
      GMN = RMUN*GD* (R_GMN_U_D(IN) * FRAC2 + R_GMN_U_D(IN+1) *FRAC1)   
      GEP =      GD* (R_GEP_D  (IN) * FRAC2 + R_GEP_D  (IN+1) *FRAC1)   
      GEN =      GD* (R_GEN_D  (IN) * FRAC2 + R_GEN_D  (IN+1) *FRAC1)
      RETURN        
                                                                
      END                                 
