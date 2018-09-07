C=======================================================================
                                                                        
      SUBROUTINE SECNUCLW(E0,EP,TH,SIGMA)                               
!--------------------------------------------------------------------
! SIGMA is cross section nb/GeV/str/nucleus.                                       
C Bastardization of Javier's SECNUC, inspired by patches, patches, and  
C patches splattered all over XSECFIT, INEFT, and this routine. Target  
C information is passed through common /targt/.                         
                                                                        
C SECNUCLW calculates the nuclear cross section in the deep inelastic   
C region using a variety of fits.                                  
C Added Fermi smearing Sept. 2005 P. Bosted
C Changed arguement of INELAST to WSQ rather than W
c Changed 8/06 to use better smearing pyb
!--------------------------------------------------------------------
      implicit none
      real e0,ep,th,sigma,sinsq,cossq,tansq,qsq,wsq,w,csmott,f1c,x
      real avgn,avga,avgM,amuM,r,dr,nu,eps,kappa,sigres,flux
      REAL*8 W1,W2,sigt,rc,q2,w1t,w2t
      INTEGER ISM
      logical goodfit
      integer iz,ia,i,iwp
      logical smrelasp, usearen
      real smrdEp,smrdW
      common/smrproton/ smrelasp,smrdEp,usearen
      COMMON      /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                     
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET                                             
      real PM/0.93828/,THCNST/0.01745329/,ALPHA/137.0388/        
      real*8 PI/3.1415926535D0/
      logical usegrd,first/.true./
      common/testing/prttst,usegrd
      logical prttst
      real*8 xval1(50),xvall(50),temp(4)
      integer iq,iw
      real*8 w1sv(0:1500,0:110),w2sv(0:1500,0:110),w1tmp(0:1500)
      real*8 delw,w1a,w1b,delq2,w2a,w2b,q2lo,q2hi,w2tmp(0:1500),psm(7)

c make grid if first time through
      if(usegrd.and.first) then
        do iq=0,110
          do iw=0,1500
c        do iq=0,100
c          do iw=0,1000
            wsq = -10. + 20. * (iw-0.5) / 1000.
cc            q2 = 0.1 * iq
c change to use log(q2) grid
            q2 = 10**(-2. + 3.*float(iq)/100.)
            nu = (wsq - pm**2 + q2) / 2. / pm
            x = q2 / 2. / pm / nu
            W1sv(iw,iq) = 0.
            W2sv(iw,iq) = 0.
            if(nu.gt.0.0 .and. x.lt.2.) then
             CALL INELAST(q2,DBLE(wsq),W1sv(iw,iq),W2sv(iw,iq)) 
c scale by dipole to make smoother
             W1sv(iw,iq) = W1sv(iw,iq) * (1. + q2/0.71)**2
             W2sv(iw,iq) = W2sv(iw,iq) * (1. + q2/0.71)**2
             if(w1sv(iw,iq).lt.0.0.or.w2sv(iw,iq).lt.0.) 
     >        write(6,'(''negative w1,w2'',2i4,2f10.3)') 
     >        iw,iq, W1sv(iw,iq),W2sv(iw,iq)
            endif
          enddo
! Add resolution smearing. 
          if(smrelasp) then
           do iw=0,1000
            w1tmp(iw) = W1sv(iw,iq)
            w2tmp(iw) = W2sv(iw,iq)
            W1sv(iw,iq) = 0.
            W2sv(iw,iq) = 0.
           enddo
! THIS IS HARD-WIRED for typical experiment
           smrdW = 0.030
           psm(1) = exp(-0.060**2 / smrdW**2) / sqrt(3.1416)
           psm(2) = exp(-0.040**2 / smrdW**2) / sqrt(3.1416)
           psm(3) = exp(-0.020**2 / smrdW**2) / sqrt(3.1416)
           psm(5) = exp(-0.020**2 / smrdW**2) / sqrt(3.1416)
           psm(6) = exp(-0.040**2 / smrdW**2) / sqrt(3.1416)
           psm(7) = exp(-0.060**2 / smrdW**2) / sqrt(3.1416)
           psm(4) = 1. - psm(1)- psm(2)- psm(3)- psm(5)- psm(6)- psm(7) 
           do iw=3,997
             do iwp = iw-3,iw+3
               W1sv(iwp,iq) = W1sv(iwp,iq) + w1tmp(iw) * psm(iwp+4-iw)
               W2sv(iwp,iq) = W2sv(iwp,iq) + w2tmp(iw) * psm(iwp+4-iw)
             enddo
           enddo
          endif
c          write(6,'(1x,''making grid...'',i4,f7.3,4f10.6)') iq,q2,
c     >      w1sv(600,iq),w1sv(900,iq),w2sv(600,iq),w2sv(900,iq)
        enddo
        first=.false.
      endif

C Calculate QSQ and WSQ for this kinematic point:                       
                                                                        
      SINSQ  = SIN(TH*3.1415927/180./2.)**2                             
      COSSQ  = 1.0-SINSQ                                                
      TANSQ  = SINSQ/COSSQ                                              
      QSQ    = 4.0*E0*EP*SINSQ                                          
      WSQ    = PM*PM+2*PM*(E0-EP)-QSQ                                   
      NU     = E0 - EP
      sigma = 0.

      CSMOTT = (.001)*(19732./(2.0*ALPHA*E0*SINSQ))**2*COSSQ            

! Moved Fermi smearing into inelast
      if(usegrd) then
        w1 = 0.
        w2 = 0.
c        iq = int(qsq*10.) 
c changed to log10 q2 bins
        iq = 0
c        if(qsq.gt.0.01) iq = (alog10(qsq) + 2.0) * 100. / 3.0
        if(qsq.gt.0.01) iq = (alog10(qsq) + 2.0) * 100. / 3.0
        iw = -1
        if(wsq .gt.-10.0) iw = int((wsq + 10.)*50.)
c        if(iw.ge.0.and.iw.lt.1000 .and.
c     >     iq.ge.0.and.iq.lt.100) then
        if(iw.ge.0.and.iw.lt.1500 .and.
     >     iq.ge.0.and.iq.lt.110) then
          delw = (wsq + 10.)*50. - iw
          if(delw.lt.0.0.or.delw.gt.1.) write(6,'(''error delw2'',
     >      i2,2f10.4)') iq,sqrt(wsq),delw
          W1a = w1sv(iw,iq)   + delw*(w1sv(iw+1,iq)  -w1sv(iw,iq))
          W1b = w1sv(iw,iq+1) + delw*(w1sv(iw+1,iq+1)-w1sv(iw,iq+1))
          W2a = w2sv(iw,iq)   + delw*(w2sv(iw+1,iq)  -w2sv(iw,iq))
          W2b = w2sv(iw,iq+1) + delw*(w2sv(iw+1,iq+1)-w2sv(iw,iq+1))
cc          delq2 = qsq*10. - iq
          q2lo = 10**(-2. + 3.*float(iq)/100.)
          q2hi = 10**(-2. + 3.*float(iq+1)/100.)
          delq2 = (qsq - q2lo) / (q2hi - q2lo)
          if(delq2.lt.0.0.or.delq2.gt.1.) then
c            write(6,'(''error delq2'',
c     >      i2,4f10.4)') iq,q2lo,qsq,q2hi,delq2
            delq2 = 0.
          endif
          w1 = w1a + delq2 * (w1b - w1a)
          w2 = w2a + delq2 * (w2b - w2a)
c divide out dipole
          w1 = w1 / (1. + qsq/0.71)**2
          w2 = w2 / (1. + qsq/0.71)**2
        endif
      else
        CALL INELAST(DBLE(QSQ),DBLE(WSQ),W1,W2)!get s.f. /nucleon 
        if(w1.lt.0.0.or.w2.lt.0.0) then
          write(6,'(1x,''error, e,ep,th,q2,w,w1,w2='',6f8.3)') e0,ep,
     >      th,qsq,w,w1,w2
        endif
      endif

      SIGMA = (W2+2.*TANSQ*W1)*CSMOTT
C             (per nucleon, including emc effect (which includes neutron
C              excess correction), if any)                              
                                                                        
      SIGMA  = avgA*SIGMA  !per nucleus                                               
      if(sigma.lt.0.0.and.wsq.gt.1.0) write(6,'(1x,''error, sigma='',e10.3,2i4,
     >   2f6.3,4e10.3)') sigma,iw,iq,wsq,qsq,w1,w2,tansq,csmott
      sigma = max(0., sigma)                                                                  
! change cross section for polarized beam and target?                   
c     IF( ( (INDEX(TARGET,'E142') +INDEX(TARGET,'E143') +               
c    >       INDEX(TARGET,'E149')).GT.0) )
c    >     .AND.( (INDEX(TARGET,'_P').GT.0).OR.                         
c    >            (INDEX(TARGET,'_A').GT.0) )   )  !inelastic only      
c    > CALL ASYM_INE(E0,EP,TANSQ,QSQ,TARGET,SIGMA)                      


cc      write(6,'(1x,3f7.3,e11.4)') E0,EP,TH,SIGMA
      RETURN                                                            
      END                                                               
                           
