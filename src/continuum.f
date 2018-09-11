C=======================================================================
                                                                        
      REAL*4 FUNCTION CONTINUUM(T)                                             
      IMPLICIT NONE
      REAL T                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /CON/ PI,PM,PP,EM,AL 
      REAL PI,PM,PP,EM,AL                                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                              
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                        
      COMMON /ERR/ REPS1,AEPS1  
      REAL  REPS1,AEPS1                                        
      COMMON /SIG/ CSTYPE                                               
      CHARACTER*1  CSTYPE
      REAL E0LO,E0HI,ANS4,ANS3,ANS2,ANS1,EPHI,EPLO,QUADMO_R
      INTEGER NLVL
      REAL FUNC1,FUNC2,FUNC3,FUNC4
      EXTERNAL     FUNC1,FUNC2,FUNC3,FUNC4                              
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM                          
                                                                        
      CSTYPE = 'I'                                                      
      CALL RADIATORS(T)                                                 
C Soft photon part; E0-R*DELTA ---> E0 and EP ---> EP+DELTA:            
                                                                        
      ANS1 = FUNC1(T)                                                   
                                                                        
C Soft photons along initial electron only; single integration.         
C E0-R*DELTA --> E0 analytically and EP+DELTA --> EPMAX(E0) numerically:
                                                                        
c PEB modified to use avgM for IA>1 10/05
      EPLO = EP+DELTA                                                   
      if(IA.le.1) then
        EPHI = (E0-PP-PP**2/2./PM)/(1.+E0*(1.-COS(THR))/PM)               
      else
        EPHI = (E0-PP-PP**2/2./avgM)/(1.+E0*(1.-COS(THR))/avgM)        
      endif
C$$      ANS2 = RGAUSS(FUNC2,EPLO,EPHI,REPS1) 
      if(ephi.le.eplo) then
cc      write(6,'(1x,''error, eplo,ephi='',2f8.3)') eplo,ephi
        ans2=0.
      else
        ANS2 = QUADMO_R(FUNC2,EPLO,EPHI,REPS1,NLVL) 
      endif

C Soft photons along final electron only; single integration.           
C EP --> EP+DELTA analytically and E0MIN(EP) --> E0-R*DELTA numerically:
                                                                        
c PEB modified to use avgM for IA>1 10/05
      if(IA.le.1) then
        E0LO = (EP+PP+PP**2/2./PM)/(1.-EP*(1.-COS(THR))/PM)               
      else
        E0LO = (EP+PP+PP**2/2./avgM)/(1.-EP*(1.-COS(THR))/avgM)          
      endif
      E0HI = E0-R*DELTA                                                 
C$$     ANS3 = RGAUSS(FUNC3,E0LO,E0HI,REPS1)                              
      if(ephi.le.eplo) then
c        write(6,'(1x,''error, e0lo,e0hi='',2f8.3)') e0lo,e0hi
        ans3=0.
      else
        ANS3 = QUADMO_R(FUNC3,E0LO,E0HI,REPS1,NLVL)    
      endif

C Hard photon part; double integration.  Set up the integration over the
C initial energy of the electrons.                                      
                                                                        
c PEB modified to use avgM for IA>1 10/05
      if(IA.le.1) then
        E0LO = (EPLO+PP+PP**2/2./PM)/(1.-EPLO*(1.-COS(THR))/PM)           
      else
        E0LO = (EPLO+PP+PP**2/2./avgM)/(1.-EPLO*(1.-COS(THR))/avgM)   
      endif
C$$      ANS4 = RGAUSS(FUNC4,E0LO,E0HI,REPS1)                              
      if(ephi.le.eplo) then
cc      write(6,'(1x,''error4, e0lo,e0hi='',2f8.3)') e0lo,e0hi
        ans4=0.
      else
        ANS4 = QUADMO_R(FUNC4,E0LO,E0HI,REPS1,NLVL)  
      endif
      CONTINUUM = ANS1+ANS2+ANS3+ANS4                                   
cc      if(CONTINUUM.gt.0.) 
c       write(6,'(1x,''continuum='',5e12.4)') ans1,ans2,ans3,ans4,
c     >  CONTINUUM                                                                       
      RETURN                                                            
      END                                                               
                               
