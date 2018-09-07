C=======================================================================
                                                                        
      SUBROUTINE ELASTIC(E0,EP,TH,SIGMA)                                
                                                                        
C subroutine to calculate electron-nucleus electron scattering cross    
C section. Formulas from ASTRUC of Mo-Tsai program.                     
      IMPLICIT NONE
      REAL E0,EP,TH,SIGMA       
      REAL CHBAR,ALPHA,PM,PI,FDEL,FSHELL,FGAUSS
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828)                                       
      PARAMETER (PI    = 3.1415927)                                     

      COMMON    /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      REAL THR,QSQ,TAU
      REAL W1,W2,FF,CMOTT,CSMOTT,RECOIL
      REAL A_ULMAR,B_ULMAR,w2_old
      REAL*8 GE,GM,GEP,GEN,GMP,GMN
                                                                  
      sigma = 0.
c removed 3/25/08
c      if(iz.lt.1) return

      THR = TH*PI/180.                                                  
      QSQ = 4.*E0*EP*SIN(THR/2.)**2 
      CALL NUC_FORM_FACTOR(QSQ,W1,W2,FF)                        

c THIS SKIPS REST OF CODE!!! WHY???
      IF(1.EQ.1) GO TO 1111
      TAU = QSQ/4./PM**2                                                
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)                                   
                                                                        
      IF (iA.EQ.1) THEN                                                 
           W1 = TAU*GMP**2                                              
           W2 = (GEP**2+W1)/(1.+TAU)                                    
      ELSEIF (iA.EQ.2) THEN  
        IF((IDUT.GE.11).AND.(IDUT.LE.14).AND.(QSQ.LE.3.5))THEN   !Tjon fits
           !Ulmar d2 elastic mode                                       
           CALL DEUT_U1(IDUT-10,QSQ,A_ULMAR,B_ULMAR)                     
           W1 = B_ULMAR/2.  ! sigma = sig_mot(A + B*tan...)
           W2 = A_ULMAR 
        ELSEIF(IDUT.EQ.1) THEN ! Linda Stuart's Model installed 5/30/96
           CALL FFD(DBLE(QSQ),GE,GM)   
           TAU = QSQ/4./avgM**2  
           W1 = TAU*GM**2                                              
           W2 = (GE**2+W1)/(1.+TAU) 
        ELSE   ! old  elastic deuterium from original code   
           FF  = FDEL(QSQ)                                              
           W1  = FF**2*TAU*.6667*(GMN+GMP)**2                           
           W2  = W1+(FF*(avgN*(GEN+TAU*GMN)+GEP+TAU*GMP)/(1.+TAU))**2   
        ENDIF
      ELSEIF (iA.EQ.3) THEN  !3HE  ! added 5/30/96  SER
        CALL FFHE3(DBLE(QSQ),GE,GM)
           TAU = QSQ/4./avgM**2  
           W1 = TAU*GM**2                                              
           W2 = (GE**2+W1)/(1.+TAU)  
           W2_old  = (iz*FSHELL(QSQ) )**2
      ELSEIF (iA.LE.20) THEN
           FF  = FSHELL(QSQ)                                            
           W1  = 0.                                                     
           W2  = (iZ*FF)**2                                             
      ELSE                                                              
           FF  = FGAUSS(QSQ)                                            
           W1  = 0.                                                     
           W2  = (iZ*FF)**2                                             
      ENDIF                                                             
 1111 CONTINUE

      CMOTT  = CHBAR**2*0.001*ALPHA**2/4.                               
      CSMOTT = CMOTT*COS(THR/2.)**2/(E0*SIN(THR/2.)**2)**2              
      RECOIL = avgM/(avgM+E0*(1.-COS(THR)))                             
      SIGMA  = (W2+2.*W1*TAN(THR/2.)**2)*CSMOTT*RECOIL                  
!Change elastic cross section by asymetry if doing asymetry run on prot 
c     IF((                                                              
c    >    ( (   (INDEX(TARGET,'E142')+INDEX(TARGET,'E143') ).GT.0)      
c    >           .AND.(IA.EQ.1)                                         
c    >    )   .OR.                                                      
c    >    ( (INDEX(TARGET,'E149') .GT.0 ).AND.(IA.EQ.2)!fake el d2 asym 
c    >    )                                                             
c    >   )                                                              
c    >      .AND.( INDEX(TARGET,'_P').GT.0)                             
c    >   )                                                              
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >               CSMOTT*RECOIL,SIGMA)                               
      RETURN                                                            
      END                                                               
                                                                       
