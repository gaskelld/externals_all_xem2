C=======================================================================
                                                                        
      REAL*4 FUNCTION QEPEAK(T)                                                
                                                                        
C Function to integrate the quasi elastic radiative cross section in    
C the delta function approximation for the peak in the cross section.   
C The only integral is over the E0.                                     
      IMPLICIT NONE
      REAL*4 T                                                                        
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
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM !put in by SER 4/14/93  
      INTEGER IZ,IA 
      REAL AVGN,AVGA,AVGM,AMUM   
      REAL*4 EPP,EPLO,E0HI,E0LO,E0P,FAC,SIGEFF,AI11,AI1,AI2,
     >       ANS1,ANS2,ANS3,SIGMA
      REAL*4 FTSAI,BREMS,QUADMO_R
      INTEGER NLVL
      CHARACTER*1  CSTYPE                                               
      EXTERNAL     FUNCQE                                               
                                                                        
      QEPEAK = 0.                                                       
      IF (iA.LE.1) RETURN                                               
      CSTYPE = 'P'                                                      
      CALL RADIATORS(T)                                                 
                                                                        
C Soft photon part; E0-R*DELTA ---> E0 done analytically along incident 
C electron. EP for elastic peak for incident electron of energy E0:     
                                                                        
      EPP    = PM*E0/(PM+2.*E0*SIN(THR/2.)**2)                          
      CALL CROSS(E0,EPP,TH,SIGMA)   !EPP not used here
      SIGEFF = SIGMA*FTSAI(E0,EPP,THR)                                  
                                                                        
C Analytical evaluation of the integral assuming constant cross section 
C and DELTA/E0 << 1.                                                    
                                                                        
      AI11 = 1.+0.5772*BTBI-0.62*BTBI**2                                
      AI1  = AI11*(R*DELTA/E0)**BTB*(1.-ATBX0/(1.-BTBI)/R/DELTA)        
      AI2  = BREMS(EPP,EPP-EP,TA,2)                                     
      ANS1 = AI1*SIGEFF*AI2                                             
                                                                        
C Hard photon part; E0LO ---> E0-R*DELTA done numerically:              
                                                                        
      EPLO = EP+DELTA                                                   
      E0HI = E0-R*DELTA                                                 
      E0LO = PM*EPLO/(PM-2.*EPLO*SIN(THR/2.)**2)                        
C$$      ANS2 = RGAUSS(FUNCQE,E0LO,E0HI,REPS1)                             
      ANS2 = QUADMO_R(FUNCQE,E0LO,E0HI,REPS1,NLVL)
C Soft photon part; E0P(EP+DELTA) ---> E0P(EP) done analytically along  
C scattered electron. E0 for elastic peak for scattered electron of     
C energy EP:                                                            
                                                                        
      E0P    = PM*EP/(PM-2.*EP*SIN(THR/2.)**2)                          
      CALL CROSS(E0P,EP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0P,EP,THR)                                  
      FAC    = (PM+2.*E0P*SIN(THR/2.)**2)/(PM-2.*EP*SIN(THR/2.)**2)     
                                                                        
C Analytical evaluation of the integral assuming constant cross section 
C and DELTA/EP << 1.                                                    
                                                                        
      AI1  = BREMS(E0,E0-E0P,TB,1)                                      
      AI11 = 1.+0.5772*BTAI-0.62*BTAI**2                                
      AI2  = AI11*(DELTA/EP)**BTA*(1.-ATAX0/(1.-BTAI)/DELTA)            
      ANS3 = AI1*SIGEFF*FAC*AI2                                         
                                                                        
      QEPEAK = ANS1+ANS2+ANS3                                           
                                                                        
      RETURN                                                            
      END                                                               
                                           
