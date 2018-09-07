C=======================================================================
                                                                        
      SUBROUTINE QELASTIC(E0,TH,SIGMA,W1,W2)                               
                                                                        
C Subroutine to calculate the quasi elastic cross section using NFORM.  
      IMPLICIT NONE
      REAL E0,TH,SIGMA                             
      REAL CHBAR,ALPHA,PM,PM24,PI
      PARAMETER (CHBAR = 19732.8)                                       
      PARAMETER (ALPHA = 7.29735E-03)                                   
      PARAMETER (PM    = 0.93828, PM24=3.5216)                                       
      PARAMETER (PI    = 3.1415927)                                     
      REAL FF/1./  !phony initialization to please compiler: SER 4/15/93
      COMMON/TARGT/ iZ,iA,avgN,avgA,avgM,amuM                       
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL THR,QSQ,TAU,W1,W2,CMOTT,CSMOTT,RECOIL,Pauli_sup1,PAULI_SUP2
      REAL EP,SINSQ,RAD,THR2
      REAL*8 GEP,GEN,GMP,GMN
      LOGICAL FIRST/.TRUE./

      IF(FIRST) THEN
         CMOTT  = CHBAR**2*0.001*ALPHA**2/4.  
         RAD = PI/180.
         FIRST=.FALSE.
      ENDIF
                                                                        
      SIGMA = 0.                                                        
      IF (iA.LE.1) RETURN                                               
                                                                        
      THR = TH *RAD
      THR2 = THR/2.
      SINSQ =  SIN(THR2)**2
      EP    = PM*E0/(PM+2.*E0*SINSQ)                                      
      QSQ = 4.*E0*EP*SINSQ                                     
      TAU = QSQ/PM24                                               
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)  

! Pauli suppression model
      CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     > PAULI_SUP1,PAULI_SUP2)

      W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)
      W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU) 
      CSMOTT = CMOTT*(1.-SINSQ)/(E0*SINSQ)**2              
      RECOIL = PM/(PM+E0*(1.-COS(THR)))                                 
      SIGMA  = (W2+2.0*TAN(THR2)**2*W1)*CSMOTT*RECOIL                 
! Are we changing quasi-elastic cross section due to Asymmetry?         
c     IF( ( (INDEX(TARGET,'E142')+INDEX(TARGET,'E143') +                
c    >       INDEX(TARGET,'E149') ).GT.0 )                              
c    >     .AND.(INDEX(TARGET,'_P').GT.0) )                             
c    > CALL ASYM_QEL(E0,EP,THR,QSQ,TARGET,GEN,GMN,GEP,GMP,              
c    >  CSMOTT*RECOIL,SIGMA)                                            
                                                                        
      RETURN                                                            
      END                                                      
