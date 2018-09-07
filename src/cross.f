C=======================================================================
                                                                        
      SUBROUTINE CROSS(E0,EP,TH,SIGMA)                                  
                                                                        
C calculates cross section per NUCLEUS, not per nucleon.                
      IMPLICIT NONE 
      REAL E0,EP,TH,SIGMA
      COMMON /CON/ PI,PM,PP,EM,AL   
      REAL PI,PM,PP,EM,AL,W1,W2
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM       
      COMMON /SIG/ CSTYPE                                               
      CHARACTER*1  CSTYPE  
c Coulomb corrections stuff
      COMMON /COULOMB/ doccor
      logical   doccor        

      real deltae_cc, f1cc ! boost and focus factor for coulomb correction.
      real E0_cc,EP_cc,ccor  


cc aji  ----------------------------------------------
cc apply Coulombic effect to bornmodel and radiate it
cc corrected vertex model, sigborn should be smaller: so divide
cc Actually this means that SIGMA_new=SIGMA_nominal/ccor and we 
cc radiate SIGMA_new

      if(doccor) then
         call get_cc_info(iA,iZ,deltae_cc)
         E0_cc  = E0+ deltae_cc
         EP_cc = EP+ deltae_cc                                        
      else
         E0_cc = E0
         EP_cc = EP
      endif
                                                                        
C inelastic tail case         = 'I'                                     
C quasielastic smeared case   = 'Q'                                     
C quasielastic unsmeared case = 'P'                                     
C elastic case                = 'E'                                     
                                                                        
      IF (CSTYPE.EQ.'I') CALL  SECNUCLW(E0_cc,EP_cc,TH,SIGMA)                  
      IF (CSTYPE.EQ.'Q') CALL  QUASIY8(E0_cc,EP_cc,TH,SIGMA)                  
C     IF (CSTYPE.EQ.'Q') CALL  FYXSEC8(E0_cc,EP,TH,SIGMA)                 
      IF (CSTYPE.EQ.'P') CALL  QELASTIC(E0_cc,TH,SIGMA,W1,W2)                  
      IF (CSTYPE.EQ.'E') CALL  ELASTIC(E0_cc,EP_cc,TH,SIGMA)                  
                      
      f1cc=E0_cc/E0
      IF(SIGMA.gt.0.0) THEN
         SIGMA=SIGMA*f1cc*f1cc
      ENDIF
                                                  
      RETURN                                                            
      END                                                               
                                                          
