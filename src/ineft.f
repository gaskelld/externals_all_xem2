!------------------------------------------------------------------------

c       SUBROUTINE INEFT(QQ,W,W1,W2,amuM)                      
       SUBROUTINE INEFT(QQ,W,W1,W2)
C Modified 6feb87 by lww to accept target information passed through    
C common block /targt/.                                                 
                                                                        
C This program takes the old slac structure function model (Atwood,     
C Bodek, et.al.) and outputs values for W1 and W2 at given kinematics. 
! As of 11/3/95 this version is per NEUCLEON   ! Steve Rock

! amuM is atomic number, ie. 1. 2.xxx etc.
! 11/05 changed to FITEMC_N (with neutron excess) instead of FITEMC
! 8/06 changed back to FITEMC and put in neutron excess from
! diff of d and p fits pyb peb
c this doesnt work: put in CONSTNAT n/p ratio for test
      Implicit None 
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM 
      integer iz,ia
      REAL avgN, avgA, avgM, amuM               
      REAL*8 QQ,W,W1,W2,WW,V,VV,OMEGAP,SP,UNIV,BRES,SLACF2,B,fneut
      REAL*8 VW2,X,EMCFAC,UNIVD,BRESD,UNIVP,BRESP,vW2p,vW2d,vW2n
      REAL*8    C(24),CF(11),CD(24),CFD(11)                          
      REAL*8    EF(7) 
      REAL FITEMC
c      REAL FITEMC_N
                                                  
      REAL*8         PM / .93828/,PMPM/.8803/,TPM/1.876512/            
      REAL*8         R /  .18/,ALPHAX/137.0388/,THCONST/0.0174533/
      LOGICAL GOODFIT
      common/testing/prttst,usegrd
      logical prttst,usegrd
      DATA         EF / -0.00136693,-.00510425,-.0375986,-.0946004,     
     +                  -.122435,-.0112751,0.406435/                    
                                                                        
C FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE TERMS)    
                                                                        
      DATA   C(1) / 0.10741163E 01/,  C(2) / 0.75531124E 00/,           
     *       C(3) / 0.33506491E 01/,  C(4) / 0.17447015E 01/,           
     *       C(5) / 0.35102405E 01/,  C(6) / 0.10400040E 01/,           
     *       C(7) / 0.12299128E 01/,  C(8) / 0.10625394E 00/,           
     *       C(9) / 0.48132786E 00/,  C(10)/ 0.15101467E 01/,           
     *       C(11)/ 0.81661975E-01/,  C(12)/ 0.65587179E 00/,           
     *       C(13)/ 0.17176216E 01/,  C(14)/ 0.12551987E 00/,           
     *       C(15)/ 0.74733793E 00/,  C(16)/ 0.19538129E 01/,           
     *       C(17)/ 0.19891522E 00/,  C(18)/-0.17498537E 00/,           
     *       C(19)/ 0.96701919E-02/,  C(20)/-0.35256748E-01/,           
     *       C(21)/ 0.35185207E 01/,  C(22)/-0.59993696E 00/,           
     *       C(23)/ 0.47615828E 01/,  C(24)/ 0.41167589E 00/            
                                                                        
C CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION) OMEGAW FIT   
                                                                        
      DATA CF(1) / 0.25615498E 00/,  CF(2) / 0.21784826E 01/,           
     *     CF(3) / 0.89783738E 00/,  CF(4) /-0.67162450E 01/,           
     *     CF(5) / 0.37557472E 01/,  CF(6) / 0.16421119E 01/,           
     *     CF(7) / 0.37635747E 00/,  CF(8) / 0.93825625E 00/,           
     *     CF(9) / 0.10000000E 01/,  CF(10)/ 0.0           /,           
     *     CF(11)/ 0.50000000E 00/                                      
                                                                        
C FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONANT TERMS)   
                                                                        
      DATA  CD(1) / 0.10521935E 01/, CD(2) / 0.76111537E 00/,           
     *      CD(3) / 0.41469897E 01/, CD(4) / 0.14218146E 01/,           
     *      CD(5) / 0.37119053E 01/, CD(6) / 0.74847487E 00/,           
     *      CD(7) / 0.12399742E 01/, CD(8) / 0.12114898E 00/,           
     *      CD(9) / 0.11497852E-01/, CD(10)/ 0.14772317E 01/,           
     *      CD(11)/ 0.69579815E-02/, CD(12)/ 0.12662466E 00/,           
     *      CD(13)/ 0.15233427E 01/, CD(14)/ 0.84094736E-01/,           
     *      CD(15)/ 0.74733793E 00/, CD(16)/ 0.19538129E 01/,           
     *      CD(17)/ 0.19891522E 00/, CD(18)/-0.24480414E 00/,           
     *      CD(19)/ 0.14502846E-01/, CD(20)/-0.35256748E-01/,           
     *      CD(21)/ 0.35185207E 01/, CD(22)/-0.21261862E 00/,           
     *      CD(23)/ 0.69690531E 01/, CD(24)/ 0.40314293E 00/            
                                                                        
C CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)        
C OMEGAW FIT                                                            
                                                                        
      DATA CFD(1) / 0.47708776E 00/, CFD(2) / 0.21601918E 01/,          
     *     CFD(3) / 0.36273894E 01/, CFD(4) /-0.10470367E 02/,          
     *     CFD(5) / 0.49271691E 01/, CFD(6) / 0.15120763E 01/,          
     *     CFD(7) / 0.35114723E 00/, CFD(8) / 0.93825625E 00/,          
     *     CFD(9) / 0.10000000E 01/, CFD(10)/ 0.0           /,          
     *     CFD(11)/ 0.50000000E 00/                                     
                                                                        
C COMPUTE SOME KINEMATIC QUANTITIES                                     
                                                                        
      WW     = W**2                                                     
      V      = (WW+QQ-PMPM)/2.D0/PM                                     
      VV     = V*V                                                      
      OMEGAP = TPM*V/QQ+PMPM/QQ                                         
                                                                        
C OVERCOME RISK OF UNDERFLOW IN THE EXPONENTIATION                      
      OMEGAP = DMIN1(20.0D0,OMEGAP)                                       
                                                                        
      SP = 1.0-EXP(-7.7*(OMEGAP-1.0))                                   
C pyb *** modified to get proton too when IA ne IZ*2
      IF (IA. ne. 2*IZ) THEN 
C          UNIVERSAL AND RESONANCE FIT FOR HYDROGEN                     
           UNIVP = SLACF2(W,QQ,CF)                                      
           BRESP = B(W,QQ,C)                                            
      endif
      if(amuM.ge.1.5) then
C          UNIVERSAL AND RESONANCE FIT FOR DEUTERIUM                    
           UNIVd = SLACF2(W,QQ,CFD)/SP
           BRESd = B(W,QQ,CD)          
      ENDIF                                                             
                                                                        
C COMPUTE VW2,W2,W1                                                     
                                                                        
      if(amuM.le.1.5) then
        vW2 = UNIVp * BRESp
      else
        VW2p    = UNIVp * BRESp 
        VW2d    = UNIVd * BRESd / 2. ! per nucleon
        vW2n    = 2.*vW2d - vW2p
        fneut = float(IA - 2*IZ) / float(IA)
c       if(fneut.lt.0.) write(6,'(''ERROR, fneut='',f8.3,2i3)') 
c    >    fneut,ia,iz
c *** set fneut back to zero: assume n/p=1!
        fneut=0.
        vW2 = (1. - fneut) * vW2d + fneut * vW2n
c        if(prttst) write(98,'(8f8.3)')  qq,w,vw2p,vw2d,vw2n,
c     >    vw2,vw2n/vw2d,fneut
      endif

      W2     = VW2/V                                                    
      W1     = (1.0D0+VV/QQ)/(V*(1.0D0+R))*VW2                          
!      if(prttst) write(*,'(1x,''univ...='',6f10.4)') sp,univ,bres,
!     >  vw2,w2,w1

      IF (amuM.LE.2.5) RETURN                                               
      X      = QQ/2./PM/V
c Modified 11/05 to include neutron excess pyb peb
c Modified 8/06 to do neutron excess as above, not in EMC
      EMCFAC= FITEMC(REAL(X),REAL(amuM),GOODFIT)
c      EMCFAC= FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GOODFIT)
                                                                        
      W2     = W2*EMCFAC                                                
      W1     = W1*EMCFAC                                                
                                                                        
      RETURN                                                            
      END                                                               
                                   
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION SLACF2(WM,QSQ,CF)                                        
                                                                        
C UNIVERSAL FUNCTION FOR ATWOOD'S FIT                                   

      Implicit none                                      
      REAL*8    WM,QSQ,CF(11)                                               
      REAL*8    PM2/1.876512/, PMSQ/.8803/, PHTR/.61993/
      REAL*8    V,OMEGA,XX,XPX,OMEGAW,ARG

                                                                        
C OMEGAW FIT...NO PHOTO-PRODUCTION COUPLING                             
                                                                        
      V      = (WM**2+QSQ-PMSQ)/PM2                                     
      OMEGA  = 2.*CF(8)*V/QSQ                                           
      XX     = 1./OMEGA                                                 
      XPX    = CF(9)+CF(10)*(XX-CF(11))**2                              
      OMEGAW = (2.D0*CF(8)*V+CF(6))/(QSQ+CF(7))                         
      ARG    = 1.-1./OMEGAW                                             
                                                                        
      SLACF2 = OMEGAW/OMEGA*ARG**3*(CF(1)+CF(2)*ARG+                    
     >         CF(3)*ARG**2+CF(4)*ARG**3+CF(5)*ARG**4)                  
      SLACF2 = SLACF2*XPX                                               
                                                                        
      RETURN                                                            
      END               

C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION B(WM,QSQ,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(24),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG,BRES
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2,ressv(4)
      INTEGER   LSPIN(4),INDEX,J,K                                             
      REAL*8    PMSQ/.8803/, PM2/1.876512/, PM/.93828/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
      common/testing/prttst,usegrd
      logical prttst,usegrd
      DATA      LSPIN/1,2,3,2/                                       
                                                                        
C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(22)+C(23)*(X-C(24))**2                                 
      PIEMSQ = (C(1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(1)) B1 = C(2)                                         
      EB1 = C(3)*(WM-C(1))                                              
      IF (EB1.LE.25.) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(4)) B2 = (1.-C(2))                                    
      EB2 = C(5)*(WSQ-C(4)**2)                                          
      IF (EB2.LE.25.0) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(INDEX)                                            
           IF (I.EQ.1) RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2              
           RMA    = C(INDEX+1)                                          
           IF (I.EQ.3) RMA=RMA*(1.D0+C(20)/(1.D0+C(21)*QSQ))            
           RWD    = C(INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           ressv(i)=res
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
c      if(prttst) write(*,'(1x,''w,q2,res='',6f7.3)') wm,qsq,
c     >  ressv
                                                                   
C FORM VW2/F2                                                           
                                                                        
      B = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
!      if(prttst) write(*,'(1x,''b...'',6f10.5)') b,bbkg,xpx,ressum                                                                  
      RETURN                                                            
      END                                                               
                  
