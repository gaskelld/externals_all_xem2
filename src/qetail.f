C=======================================================================
                                                                        
      REAL*4 FUNCTION QETAIL(T)                                                
      IMPLICIT NONE
      REAL T                                                            
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
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                  
      COMMON /ERR/   REPS1,AEPS1
      REAL  REPS1,AEPS1                                        
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM                          
      COMMON /SIG/   CSTYPE                                             
      CHARACTER*1    CSTYPE
      REAL*4 ANS,ANS1,ANS2,ANS3,ANS4,EPLO,EPCN,EPHI,E0LO,E0CN,E0HI
      INTEGER NLVL
      REAL*4   FUNC1,FUNC2,FUNC3,FUNC4Q,QUADMO_R
      EXTERNAL       FUNC1,FUNC2,FUNC3,FUNC4Q                      
      logical doing_elas
      common/doelas/ doing_elas

      QETAIL = 0.                                                    
      CSTYPE = 'Q'                                                 
      CALL RADIATORS(T)                                           
                                                                  
C Soft photon part; E0-R*DELTA ---> E0 and EP ---> EP+DELTA:       
                                                                     
      ANS1 = FUNC1(T)                                          
                                                                     
C Soft photons along initial electron only; single integration.     
C E0-R*DELTA ---> E0 analytically and EP Integral numerically.        
C lower limit defined by kinematic point and upper limit by x=M_targ.
                                                                      
      EPLO = EP+DELTA                                                 
      EPCN = E0/(1.+E0*(1.-COS(THR))/PM)                             
      EPHI = E0/(1.+E0*(1.-COS(THR))/avgM)                           
      EPCN = MAX(EPLO,EPCN)                                           
ccc pyb 04/07
ccc      EPHI = E0 - DELTA

      ANS  = QUADMO_R(FUNC2,EPLO,EPCN,REPS1,NLVL)
      ANS2 = QUADMO_R(FUNC2,EPCN,EPHI,REPS1,NLVL)                     
      ANS2 = ANS + ANS2                                               
                                                                 
C Soft photons along final electron only; single integration.         
C EP ---> EP+DELTA analytically and E0 Integral numerically.          
C upper limit defined by kinematic point and lower limit by x=M_targ. 
                                                                      
      E0LO = EP/(1.-EP*(1.-COS(THR))/avgM)                            
      E0CN = (EP+PP+PP**2/2./PM)/(1.-EP*(1.-COS(THR))/PM)             
      E0HI = E0-R*DELTA                                               
      E0CN = MIN(E0CN,E0HI)                                           
ccc pyb 04/07
ccc      E0HI = E0 - DELTA


      ANS  = QUADMO_R(FUNC3,E0LO,E0CN,REPS1,NLVL)  !NEW
      ANS3 = QUADMO_R(FUNC3,E0CN,E0HI,REPS1,NLVL)  !NEW
      ANS3 = ANS + ANS3                                               
                                                                      
C Hard photon part; double integration.  Set up the integration over t
C initial energy of the electrons. upper limit defined by delta cutof 
C and lower limit by x=M_target,above the W2=M_proton+M_pion for EP   
C integral.                                                           
                                                              
      E0LO = EPLO/(1.-EPLO*(1.-COS(THR))/avgM)                       
      ANS4 = QUADMO_R(FUNC4Q,E0LO,E0HI,0.01,NLVL)  
                                     
      QETAIL = ANS1 + ANS2 + ANS3 + ANS4                       

      if(doing_elas) then
c        write(6,'(''q'',7f6.3,4e10.3)') t,e0,ep,eplo,ephi,
c     >    e0lo,e0hi,ans1, ans2,
c     >    ans3, ans4
      endif
      RETURN                                                         
      END                                                               
                                               
