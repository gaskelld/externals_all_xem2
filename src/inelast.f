C ----------------------------------------------------------------------
                                                                        
      SUBROUTINE INELAST(QQ,WSQnom,W1sm,W2sm)                              
! Choose which inelastic model is to be used for resonance and for DIS  
! Returns structure functions per NUCLEON (11/3/95 SER)
C CHANGED TO USE WSQ RATHER THAN W peb 10/05
c changed to do smearing here now

      Implicit NONE
      REAL  avgN, avgA, avgM, amuM 
      COMMON /TARGT/ iZ, iA, avgN, avgA, avgM, amuM 
      INTEGER IZ,IA, ism, nsmr
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL   
      REAL*8 QQ,WSQ,W,X,W1,W2,F2_E665,R8,EMCFAC,F1C,RC,wsqnom,w1sm,w2sm
      REAL*8 MP2/0.8803/,MP/.93828/,F2,F1,F2nF2p,pi,alpha,sigL,sigT,Z,A
      REAL Q_RES_LIMIT/.2/                                              
      REAL QQ4,X4,F2_4,W1_4,W2_4,FITEMC_N,R4,DR4,NU,F2Allm,FNP_NMC
      REAL rnphi,rnp,rnplo,fact
      LOGICAL GD
! new variables for Fermi smearing over +/- 3 sigma. Sum of 
! fy values is 1.000, so no effect if W1 is constant with W
      REAL*8 XX(15)/-3.000,-2.571,-2.143,-1.714,-1.286,-0.857,
     >              -0.429, 0.000, 0.429, 0.857, 1.286, 1.714, 
     >               2.143, 2.571, 3.000/
      real*8 fy(15)/0.0019, 0.0063, 0.0172, 0.0394, 0.0749, 0.1186, 
     >              0.1562, 0.1712, 0.1562, 0.1186, 0.0749, 0.0394, 
     >              0.0172, 0.0063, 0.0019/
! use w**2 for d-state
      real*8 fydr(15)/0.0192,0.0303,0.0538,0.1091,0.2544,0.6857,2.0061,
     >  3.5775,2.0061,0.6857,0.2544,0.1091,0.0538,0.0303,0.0192/
! *** with NO D-sate!
       real*8 fydn(15)/0.0040,0.0102,0.0277,0.0764,0.2150,0.6410,1.9589,
     >  3.5301,1.9589,0.6410,0.2150,0.0764,0.0277,0.0102,0.0040/
! using 1.5 * (1-cos**2)**2 for dtate
       real*8 fyd(15)/ 0.0094,0.0187,0.0411,0.0970,0.2462,0.6866,2.0207,
     >  3.6003,2.0207,0.6866,0.2462,0.0970,0.0411,0.0187,0.0094/
      logical first/.true./
      REAL*8 DW2DPF,pf,kf,es,dw2des,pz,qv,PM/0.93828/

      W1sm = 0.
      W2sm = 0.
      pi = 3.141593
      alpha = 1./137.036
      nu = (wsqnom - MP**2 + qq) / 2. / MP
      if(nu .le. 0.) return

! If model is 4, use new code for all values of A
      if(INEL_MODEL.EQ.4) THEN 
        Z = IZ
        A = IA
c        call F1F2IN07(Z, A, QQ, WSQnom, F1, F2)
c        call F1F2IN09(Z, A, QQ, WSQnom, F1, F2)
      call gsmearing(Z,A,WSQnom,QQ,F1,F2)
        W1sm = F1 / MP
        W2sm = F2 / nu
c need to make xsection per nucleon here!
c (they are converted back to per nucleus in SECNUCLW
        w1sm = w1sm / A
        w2sm = w2sm / A
        return
      endif

! Apply Fermi smearing to inelatic, except for H2
      if(first) then
        IF(IA .eq. 1) nsmr=1
        if(IA .ge. 2) nsmr=15
! Modifed to use Superscaling from Sick, Donnelly, Maieron,
! nucl-th/0109032
        if(IA.eq.2) kf=0.085
        if(iA.eq.2) Es=0.0022
        if(IA.eq.3) kf=0.180
        if(iA.eq.3) Es=0.010 
cc       if(IA.eq.4) kf=0.210
        if(IA.eq.4) kf=0.180
        if(iA.eq.4) Es=0.015 
cc        if(iA.eq.4) Es=0.012 
        if(IA.gt.4) kf=0.165
        if(iA.gt.4) Es=0.015 
        if(IA.gt.7) kf=0.228
        if(iA.gt.7) Es=0.020 
        if(IA.gt.16) kf=0.230
        if(iA.gt.16) Es=0.025 
        if(IA.gt.25) kf=0.236
        if(iA.gt.25) Es=0.018 
        if(IA.gt.38) kf=0.241
        if(iA.gt.38) Es=0.028 
        if(IA.gt.55) kf=0.241
        if(iA.gt.55) Es=0.023 
        if(IA.gt.60) kf=0.245
        if(iA.gt.60) Es=0.028 
! adjust pf to give right width based on kf
        pf = 0.5 * kf 
        first = .false.
      endif

! assume this is 2 * pf * qv
      qv = sqrt(nu**2 + qq)
      DW2DPF = 2. * qv
      dw2des = 2. * (nu + PM) 
      do ism = 1,nsmr
       if(IA.le.1) wsq = wsqnom
       if(IA.eq.2) wsq = wsqnom + (-0.3 + 0.6 * 
     >        (float(ism)-0.5)/15.) * dw2dpf
       if(IA.gt.2) WSQ = WSqnom + XX(ISM) * PF * DW2DPF - es * dw2des
       if(WSQ.gt.(0.93828 + 0.13957)**2) then
        W = SQRT(WSQ)

        IF(INEL_MODEL.EQ.0) THEN                                          
         CALL INEFT(QQ,W,W1,W2)   !Bodek fit                              
! 8/06 changed to use proton fit for all nuclei too!
! modified by simple n/p model
        ELSEIF(INEL_MODEL.EQ.3) THEN 
c        CALL CHRISTY705(WSQ,QQ,F1C,RC) ! Christy H2 resonance fit
c         CALL CHRISTY31606(WSQ,QQ,F1C,RC) ! Christy H2 resonance fi
         CALL CHRISTY0507(WSQ,QQ,F1C,RC) ! Christy H2 resonance fit
         NU = (WSQ +QQ -MP2)/(2.*MP)
         W1 = F1C / MP
         W2 = W1 /(1.0 + NU*NU/QQ) * (1.0 + RC)
! Simple model for n/p. Constnat in resonance region,
! and DIS at high W. Smoothly join from W=1.7 to 2.3
! Fixed to work for any Z,A 8/06 pyb
         if(IA.gt.1) then
          rnplo = 0.66
          X4 = QQ/(WSQ -MP2 +QQ) 
          qq4 = qq
          rnphi = FNP_NMC(X4,QQ4)
          fact = max(0., min(1., (wsq-1.7)/0.6))
          rnp = rnplo * (1.-fact) + fact * rnphi
          W1 = W1 * (1. + (avgN/avgA) * rnp) / (1. + (avgN/avgA))
          W2 = W2 * (1. + (avgN/avgA) * rnp) / (1. + (avgN/avgA))
         endif


        ELSEIF((INEL_MODEL.EQ.9).OR.(INEL_MODEL.EQ.12)) THEN              
         X = QQ/(WSQ -MP2 +QQ)                                           
         IF(QQ.GT.Q_RES_LIMIT) THEN
! Only H2 model. Calls F2GLOB for DIS.  
! Remember that F2GLOB returns structure functions/nucleon for H2 and D2
! Convert to Real*4 for this call
          QQ4=QQ
          X4=X                                       
!stuart for res, F2_Glob DI
          CALL INELSTU(QQ4, X4, F2_4, W1_4, W2_4, INEL_MODEL) 
          W2=W2_4
          W1=W1_4 
         ELSE  !out of range of inelstu.                                  
          CALL INEFT(QQ,W,W1,W2)   !Bodek fit                             
         ENDIF

!NMC model for DIS, INEFT for resonance
        ELSEIF(INEL_MODEL.EQ.1) THEN  
         X4 = QQ/(WSQ -MP2 +QQ) 
         QQ4=QQ
         CALL NMCDIS(QQ4,X4,W1_4,W2_4,amuM)
         W2=W2_4
         W1=W1_4
        ELSEIF(INEL_MODEL.EQ.2) THEN   ! E665 Fit (includes resonances
         X = QQ/(WSQ -MP2 +QQ) 
         IF(amuM.LT.1.5) THEN
          CALL  F2PGLO(QQ,X,F2_E665,R8) ! Prot E665
          EMCFAC=1.
         ELSE  ! Deuteron or heavier 
          CALL  F2DGLO(QQ,X,F2_E665,R8) ! Deut E665 
!with neutron excess 8/19/98
          EMCFAC = ABS(FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GD)) 
         ENDIF
         CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)
         NU = (W**2 +QQ -MP2)/(2.*MP)
         W2 = F2_E665/NU *EMCFAC
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R4)
       ELSEIF(INEL_MODEL.EQ.6) THEN   ! F2allm
         X4 = QQ/(WSQ -MP2 +QQ) 
         qq4 = qq
         F2_4 = F2ALLM(X4,QQ4)   ! proton fit
         F2 = F2_4
         x = x4
         IF(amuM.GT.1.5) THEN
! Changed to FNP_NMC 7/19/06: seems to work much better!
c          f2nf2p = 0.976+X*(-1.34+X*(1.319+X*(-2.133+X*1.533)))    
          F2 = (FNP_NMC(X4,QQ4) + 1.) / 2.  * F2
         endif
         if(amuM.gt.2.5) then
           EMCFAC = ABS(FITEMC_N(REAL(X),REAL(IA),REAL(IZ),GD))
           F2 = F2 * EMCFAC
         endif
         CALL R1998(REAL(X),REAL(QQ),R4,DR4,GD)
         NU = (W**2 +QQ -MP2)/(2.*MP)
         W2 = F2/NU 
         W1 = (1.0 + NU*NU/QQ)*W2/(1.0 + R4)
        ELSEIF(INEL_MODEL.GE.100) THEN !joining a DIS model with resonance model
         X4 = QQ/(W**2 -MP2 +QQ) 
         QQ4=QQ
         CALL F2_JOIN(QQ4,X4,amuM,iZ,INEL_MODEL,W1_4,W2_4,iA)
         W2=W2_4
         W1=W1_4
        ELSE                                                              
         WRITE(6,'('' NOT GOOD H2 INELASTIC MODEL='',I3)')INEL_MODEL        
         STOP                                                             
        ENDIF                                                             
        if(IA.eq.1) then
          W1sm = W1
          W2sm = W2
        endif
        if(IA.eq.2) then
          W1sm = W1sm + W1 * FYD(ISM)/10.
          W2sm = W2sm + W2 * FYD(ISM)/10.
        else
          W1sm = W1sm + W1 * FY(ISM)
          W2sm = W2sm + W2 * FY(ISM)
c         if(prttst) write(8,'(1x,''insm'',10f8.3)') 
c    >          qsq, wsqp,w1,w1p,fy(ism)
        endif
       endif ! test of wsq
      enddo ! loop on ism
c     CALL INELAST(DBLE(QSQ),DBLE(WSQ),W1P,W2P)!get s.f. /nucleon 
c     if(prttst) write(98,'(1x,''inelas sm'',6f10.3)') 
c    >     qsq, w,w1,w1p,w2,w2p
                                                                        
      RETURN                                                            
      END                                                               
                                                       

                                                                        
                                                               
