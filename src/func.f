C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC1(T)                                                 
      IMPLICIT NONE
      REAL*4 T                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                      
      COMMON /OMG/ DELTA,R  
      REAL   DELTA,R                                             
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL*4 SIGMA,SIGEFF,AI11,AI1,AI2,FTSAI                                                                             
C The constant part of the radiative cross section                      
                                                                        
      CALL CROSS(E0,EP,TH,SIGMA)                                        
      SIGEFF = SIGMA*FTSAI(E0,EP,THR)                                   
      AI11   = 1. + 0.5772 * BTBI - 0.62*BTBI**2                              
      AI1    = AI11 * (R*DELTA/E0)**BTB * (1.-ATBX0 / 
     >  (1.-BTBI) / R / DELTA)      
      AI11   = 1.+0.5772*BTAI-0.62*BTAI**2                              
      AI2    = AI11*(DELTA/EP)**BTA*(1.-ATAX0/(1.-BTAI)/DELTA)          
      FUNC1  = SIGEFF*AI1*AI2                                           

c      write(6,'(''func1'',2e10.2,9f7.3)') sigma,sigeff,AI11,AI1,R,
c     >   delta,e0,btb,atbx0,BTBI,r

      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC2(EPP)                                               
      IMPLICIT NONE
      REAL*4 EPP
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                                  
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL*4 SIGMA,SIGEFF,AI11,AI1,AI2,FTSAI,BREMS
C The E0P integral is done analytically to obtain AI1                   
                                                                        
      CALL CROSS(E0,EPP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0,EPP,THR)                                  
      AI11   = 1.+0.5772*BTBI-0.62*BTBI**2                              
      AI1    = AI11*(R*DELTA/E0)**BTB*(1.-ATBX0/(1.-BTBI)/R/DELTA)      
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNC2  = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC3(E0P)                                               
      IMPLICIT NONE
      REAL*4 E0P
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2  
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                            
      COMMON /OMG/ DELTA,R     
      REAL   DELTA,R                                         
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL    
      REAL*4 SIGMA,SIGEFF,FTSAI,AI1,AI11,AI2,BREMS

C The EPP integral is done analytically to obtain AI2                   
                                                                        
      CALL CROSS(E0P,EP,TH,SIGMA)                                       
      SIGEFF = SIGMA*FTSAI(E0P,EP,THR)                                  
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
      AI11   = 1.+0.5772*BTAI-0.62*BTAI**2                              
      AI2    = AI11*(DELTA/EP)**BTA*(1.-ATAX0/(1.-BTAI)/DELTA)          
      FUNC3  = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC4(E0P)                                               
      IMPLICIT NONE
      REAL E0P                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                        
      COMMON /CON/ PI,PM,PP,EM,AL 
      REAL PI,PM,PP,EM,AL                                       
      COMMON /OMG/ DELTA,R 
      REAL   DELTA,R                                              
      COMMON /FUNCOM/ E0PP 
      REAL*4 E0PP                                                
      COMMON /ERR/ REPS1,AEPS1 
      REAL  REPS1,AEPS1
      REAL*4 EPLO,EPHI,QUADMO_U
      INTEGER NLVL
      EXTERNAL     FUNC4P                                               
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
                                                                        
C Set up the EPP integration                                            
                                                                        
      E0PP  = E0P                                                       
      EPLO  = EP+DELTA                                                  
c PEB fixed to use avgM if IA>1 10/05
      if(IA.le.1) then
        EPHI  = (E0P-PP-PP**2/2./PM)/(1.+E0P*(1.-COS(THR))/PM)            
      else
        EPHI  = (E0P-PP-PP**2/2./avgm)/(1.+E0P*(1.-COS(THR))/avgM)     
      endif
C$$      FUNC4 = UGAUSS(FUNC4P,EPLO,EPHI,REPS1)                            
      FUNC4 = QUADMO_U(FUNC4P,EPLO,EPHI,REPS1,NLVL)   
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC4Q(E0P)
      IMPLICIT NONE
      REAL E0P
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                    
      COMMON /CON/   PI,PM,PP,EM,AL
      REAL PI,PM,PP,EM,AL                                      
      COMMON /OMG/   DELTA,R
      REAL   DELTA,R                                             
      COMMON /FUNCOM/   E0PP 
      REAL*4 E0PP                                                
      COMMON /ERR/   REPS1,AEPS1 
      REAL  REPS1,AEPS1                                       
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL*4 EPLO,EPCN,EPHI, ANS1,ANS2,QUADMO_U
      INTEGER NLVL
      EXTERNAL     FUNC4P                                               
                                                                        
C Set up the EPP integration from kinematic limit to x=M_target         
C quasi elastic case                                                    
                                                                        
      E0PP = E0P                                                        
      EPLO = EP+DELTA                                                   
      EPCN = (E0P-PP-PP**2/2./PM)/(1.+E0P*(1.-COS(THR))/PM)             
      EPHI = E0P/(1.+E0P*(1.-COS(THR))/avgM)                            
      EPCN = MAX(EPLO,EPCN)                                             
                                                                        
C$$      CALL SIMP2(EPLO,EPCN,100,FUNC4P,ANS1)                             
C$$      CALL SIMP2(EPCN,EPHI,100,FUNC4P,ANS2)                             

      ANS1 = QUADMO_U(FUNC4P,EPLO,EPCN,REPS1,NLVL)  !NEW
      ANS2 = QUADMO_U(FUNC4P,EPCN,EPHI,REPS1,NLVL)  !NEW
      FUNC4Q = ANS1+ANS2                                                
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNC4P(EPP)                                              
      IMPLICIT NONE
      REAL*4 EPP                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                       
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                      
      COMMON /FUNCOM/ E0P                                                  
      REAL*4 E0P
      REAL*4 SIGMA,SIGEFF,AI1,AI2,BREMS,FTSAI                                                                        
C The complete integrand in (C.1)                                       
C Sigma_r is of course approximated by equivalent radiator method       
                                                                        
      CALL CROSS(E0P,EPP,TH,SIGMA)                                      
      SIGEFF = SIGMA*FTSAI(E0P,EPP,THR)                                 
c     if(E0-E0P.le.0.) write(6,'(1x,''ERROR,e0,e0p='',2f8.4)')
c    >  e0,e0p
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
c     if(EPP-EP.le.0.) write(*,'(1x,''ERROR,epp,ep='',2f8.4)')
c    >  ep,epp
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNC4P = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNCQE(E0P)                                              
      IMPLICIT NONE
      REAL E0P                                                                  
      COMMON /KIN/ E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                      
      COMMON /CON/ PI,PM,PP,EM,AL  
      REAL PI,PM,PP,EM,AL                                       
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                    
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL       
      REAL*4 EPP,SIGMA,SIGEFF,AI1,AI2,FTSAI,BREMS
C Integrand after the delta function killed the EPP integral. EPP is now
C constrained as below:                                                 
                                                                        
      EPP    = PM*E0P/(PM+2.*E0P*SIN(THR/2.)**2)                        
      CALL CROSS(E0P,EPP,TH,SIGMA)                                      
      SIGEFF = SIGMA*FTSAI(E0P,EPP,THR)                                 
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNCQE = SIGEFF*AI1*AI2                                           
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      REAL*4 FUNCTION FUNCE(E0P)              
      IMPLICIT NONE
      REAL E0P                                                                        
      COMMON /KIN/   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2 
      REAL   E0,EP,TH,THR,X,Q2,EPS,W2,Y,ANU,SIN2                     
      COMMON /CON/   PI,PM,PP,EM,AL  
      REAL PI,PM,PP,EM,AL                                     
      COMMON /RAD/   ETA,B,AX0,BA,AX0A,                                 
     >               ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL  
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                       
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      REAL*4 EPP,SIGMA,SIGEFF,AI1,AI2,FTSAI,BREMS
                                                                        
C For Elastic electron-Nucleus collision as opposed to e-nucleon.       
C Integrand after the delta function killed the EPP integral            
C EPP is now constrained as below; Peak is at x=target_mass             
                                                                        
      EPP    = avgM*E0P/(avgM+2.*E0P*SIN(THR/2.)**2)                    
      CALL CROSS(E0P,EPP,TH,SIGMA)                                      
      SIGEFF = SIGMA*FTSAI(E0P,EPP,THR)                                 
      AI1    = BREMS(E0,E0-E0P,TB,1)                                    
      AI2    = BREMS(EPP,EPP-EP,TA,2)                                   
      FUNCE  = SIGEFF*AI1*AI2                                           
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
