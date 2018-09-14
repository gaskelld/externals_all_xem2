      program EXTERNALS
c thetop                                                                        
C Program to calculate the measured cross section using Tsai's formula. 
C We use equivalent radiator approximation to estimate the internal     
C bremstrahlung contribution. However the integrals over radiation      
C length initial and final energy losses are done in full glory without 
C any approximations except at the edges of kinematics where the        
C divergences creep in.  Code by S. Dasu, modified by LWW.              
c Modified so list of input/output files is always looked
c for in extern.inp pb                                                       
c Modifed to allow multiple files to be looked at for Jan 05
c Modified to have Fermi smeaing of inelastic. see SECNUCLW line 2544
c Modified to have correct limits of integration for inelastic
c for nuclear target: see CONTINUUM and FUNC4. PEB 10/05
C=======================================================================
      IMPLICIT NONE
      COMMON /SET/   E0SET,EPSET,THSET      
      REAL  E0SET,EPSET,THSET   
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
      COMMON /TRL/   TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL  TTARG,TWALL,TBEAM,TSPEC
      INTEGER NSEG,JTARG,ITARG,ipsv,ithsv,iesv,iqsv,iwsv,izz
      COMMON /TARGT/ iZ,iA,avgN,avgA,avgM,amuM
      logical smrelasp,usearen
      real smrdEp
      common/smrproton/ smrelasp,smrdEp,usearen
      INTEGER IZ,IA                                                     
      REAL AVGN,AVGA,AVGM,AMUM 
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL  
      COMMON /TTYPE/ TARGET                                             
      CHARACTER*7    TARGET  
      REAL TTOT,TLEN
      INTEGER ISEG154  
      COMMON/E154/ISEG154,TTOT,TLEN                       
      REAL  SIGMA_INEL(2),DELTA_INEL(2),                       
     >      SIGMA_QELA(2),DELTA_QELA(2),                       
     >      SIGMA_QEPK(2),DELTA_QEPK(2),                       
     >      SIGMA_QE_PK_COS(2),DELTA_QE_PK_COS(2),  
     >      SIGMA_ELPK(2),DELTA_ELPK(2),                       
     >      TOTAL(2),ZERO(18),SIGMA_CORR,Q2SET,SIGMA_BORN,RAT_SMR,
     >      VTRUN,TTRUN,VTSTOP,TTSTOP,VTSTART,TTSTART,X0
      REAL  SIGMA_INEL_154,SIGMA_QELA_154,SIGMA_QE_PK_COS_154,
     > SIGMA_ELPK_154,SIGMA_BORNPLSQE,SIGMA_elastic
      real e0sv(100000),epsv(100000),thsv(100000),xsv(100000)
      real w2sv(100000),sbsv(100000),srsv(100000),srqesv(100000)
      real sbqesv(100000),sbisv(100000),srinsv(100000)
      real srelsv(100000), ddnom
      INTEGER I,II,ISEG_NUM,j,npts
      REAL           CONTINUUM,QETAIL,QEPEAK,EPEAK,ATAILFL1,ATAILFL_QE
      EXTERNAL       CONTINUUM,QETAIL,QEPEAK,EPEAK,ATAILFL1,ATAILFL_QE
      logical usegrd
      common/testing/prttst,usegrd
      logical prttst,doeg1b,good
      common/experiment/ doeg1b
      real*4  xval(37),tmptmp,psf1(5),psf2(5),dr,w2last,ymax,ymin
      real*8 q28,w8,w18,w28,rcsim
C Stuff for Coulomb corrections
      real E0SET_cc, EPSET_cc, f1cc, ccor, deltae_cc
      real SIGMABORN_CC, SIGMACORR_CC, SIGTOT_CC

      COMMON/COULOMB/ doccor
      logical doccor

      COMMON/ioana/xval

      real*8 jane0(23)/1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,
     >   2.3, 2.3, 2.3, 2.3, 3.4, 3.4, 3.4, 3.4, 3.4, 3.4,
     >  4.6, 4.6, 4.6, 4.6/
      real*8 janth(23)/10.8, 13., 16., 19., 22., 28., 45., 55., 70.,
     >   20., 30., 45., 60., 14., 20., 28., 36., 40., 45.,
     >   10.65, 16.0, 20.0, 25.0/
c      logical dojan05/.true./
      logical dojan05/.false./
c      logical dofromin/.true./
      logical dofromin/.false./
c      logical dosimona/.true./
      logical dosimona/.false./
c      logical donucr/.true./
      logical donucr/.false./
c Put "Eg1b" on first line of extern.inpe for speical eg1b in/out format
      integer nfile, ifile, IAjan, itrgsim
      real*8 sigexp(100000),sigexper(100000),fact,Afromin,Zfromin
      real*8 qv,ams,amf2,a,bb,c,disc,costh,ysv(100000),xplt(100000)
      real*8 ysvp(100000)
      real*8 thfromin(6)/18.,22.,26.,32.,40.,50./
      real*8 thsimona(6)/38.,41.,45.,55.,60.,70./
      character*30 fname,fnamep
! used to tell quasiy8 if really smearing elastic peak
      logical doing_elas
      common/doelas/ doing_elas

C set true if you want to calculate the COulomb corrections AND
c incorporate in your radiative correction
      doccor=.true.

c set ture if doing many points. 
c ONLY USE if Jlab kinematics up to 6 gev
      usegrd=.true.
c      usegrd=.false.

c set true if want proton elastic to be smeared and treated as q.e.
c if so, will smear by amound smrdEp
cc code currently modified so only does this for A=1
cc need to change back if want it for A>1 or q.e., inel. smearing
cc      smrelasp = .true.
      smrelasp = .false.
cc see below for smrdEp since depnds on kinematics

c set true if want to use Arenhoevel for q.e., threshold
cc      usearen = .true.
      usearen = .false.

cc      DELTA = 0.010                                                 
! xxx try making much smaller to get A(q2), q.e. at low Q2
      DELTA = 0.002
c      write(6,'(''A b='',e11.3)') b

      REPS1 = 0.0002  !0.001  changed 21 Nov 00
      AEPS1 = 1.E-19                                                    
      PI    = 3.1415927                                                 
      PM    = 0.93828                                                   
      PP    = 0.13957                                                   
      EM    = 0.000511                                                  
      AL    = 7.29735E-03                                               
                                                                        
      call ARENHOVEL_INIT()

c      write(6,'(1x,''got here'')')
c      open(unit=17,file='output/externals/kpp_shms_488.out')
c      open(unit=17,file='output/externals/test.out')
c      open(unit=18,file='louk.out')
c      open(unit=19,file='ratioy.top')
c      open(unit=29,file='ratiow2.top')
c      write(19,'(1x,''set device postscript portrait'')')
c      write(29,'(1x,''set device postscript portrait'')')
c      write(18,1818)
c 1818 format(
c     >  1x,'                                   internal       ',
c     >  1x,'                     internal  + external'/
c     >  1x,'   E    Ep    Th     s_born      s_rad     s_elas',
c     >      '      s_dis      s_rad',
c     >     '     s_elas      s_dis')
CCCC      CALL READIN(LEN_FILENAME,FILENAME)                                                       

      do i=1,4
        do j=1,4
          x=0.1*i
          Q2=j
          CALL R1998(X,Q2,R,DR,GOOD)
c          write(6,'(1x,''x,q2,r,dr'',4f8.3)') x,q2,r,dr
        enddo
      enddo

c      write(6,'(1x,''got here'')')
c      write(6,'(''A b='',e11.3)') b
      CALL READIN_EXT()                                                       

c      write(6,'(1x,''got here'')')
c      write(6,'(''A b='',e11.3)') b
      CALL RADLENGTH(X0)                                                
c      write(6,'(''A b='',e11.3)') b
                                                                        
c read in params for d2model
c           write(6,*)'Reading parameters for Ioana''s model'
           open(unit=77,file
     >          ='input/d2model/parms.dat'
     >          ,status='old')
           do I=1,37
              read(77,*)j
              read(77,*)xval(i)
*              write(*,*)xval(i)
              read(77,*)tmptmp
              read(77,*)tmptmp
              read(77,*)tmptmp
           enddo	
           close(77)

C Constants for bremstrahlung and ionization spectra for target and for 
C aluminum:                                                             
c changed to uze iz=1 if neutron
      izz = max(iz,1)
      ETA  = LOG(1194./izz**0.6667)/LOG(184.15/izz**0.3333)               
      B    = 4./3.*(1.+1./9.*(izz+1.)/(izz+ETA)/LOG(184.15/izz**0.3333))   
      AX0  = 0.000154*izz/amuM*X0                                        
      BA   = 1.3667511                                                  
      AX0A =  .0017535                                                  
c      write(6,'(''b...'',5f10.3)') eta,b,ax0,ba,ax0a

C     TASUM=0.                                                          
C     TBSUM=0.                                                          
C Unless specified, set default number of target segments for integrals:
c                                                                        
cc took this out. nseg=0 means skip external
cc      IF (NSEG.LE.0) NSEG = 4                                           
                                                                        
      WRITE(10,'(                                                       
     >''    E1    TH    Q2 COS_K         SMR     NUC  '',     
     >''        INT-SMR'',/,                                               
     >''    EP     X    W2 QE-PK  INELA  Q-ELA  ELAST'',     
     >''   TOTAL  MODEL'',/)')                                              

c check out pauli supp.
      do ii=1,5,2
        do j=1,10
          e0set=float(ii)
          q2set=0.05*j
          do i=1,5
            call PAULI_SUPPRESSION(i-1,1,2,q2set,e0set,
     >        psf1(i),psf2(i))
          enddo
c          write(6,'(1x,12f5.2)') e0set,q2set,psf1,psf2
        enddo
      enddo

      
c Loop over kinematic points:                                           
      nfile=0
      if(dojan05) nfile = nfile + 1
      if(dofromin) nfile = nfile + 1
      if(dosimona) nfile = nfile + 1
      if(donucr) nfile = nfile + 1
      if(nfile.gt.1) then
        write(6,'(''cant do more than one special case at same time'')')
        stop
      endif
      nfile = 1
      if(dojan05) nfile=23
      if(dofromin) nfile=6
      if(dosimona) nfile=6
      if(donucr) nfile=28
c      if(donucr) nfile=2 ! for test
      if(donucr) open(unit=3,file='e99118/files.txt')
      do ifile=1,nfile
c      do ifile=20,20
      if(dojan05) then 
        close(unit=7)
        open(unit=7,file='jan05.dat')
      endif
      if(dofromin) then 
        close(unit=7)
c old        open(unit=7,file='fromin.newdat')
        open(unit=7,file='fomin.dat')
      endif
      if(dosimona) then 
        close(unit=7)
        open(unit=7,file='simona.dat')
      endif
      if(donucr) then 
        read(3,'(a)') fnamep
        write(fname,'(''e99118/h''a)') fnamep(2:20)
        if(IA.eq.2) write(fname,'(''e99118/d''a)') fnamep(2:20)
        close(unit=7)
        open(unit=7,file=fname)
      endif

      npts=0                                                                     
      w2last=-100.
      if(dofromin) w2last=100.
      DO 99 ii=1,100000                                                   
c           CALL TIME(VTSTART,TTSTART)                                   
        if(dojan05) then
          read(7,*,end=100) IAjan,E0SET,EPSET,THSET,
     >      sigexp(npts+1),sigexper(npts+1) 
          if(IAjan.ne.IA) goto 99
          if(abs(thset - janth(ifile)).gt.0.5) goto 99
          if(abs(e0set - jane0(ifile)).gt.0.5) goto 99
        elseif(dofromin) then
          e0set=5.77
          read(7,*,end=100) EPSET,THSET,Afromin,Zfromin,
     >      sigexp(npts+1),sigexper(npts+1) 
          if(abs(Afromin - float(IA)).gt.0.3) goto 99
          if(abs(Zfromin - float(IZ)).gt.0.3) goto 99
          if(abs(thset - thfromin(ifile)).gt.0.5) goto 99
c          write(6,'(1x,''fromin'',6f10.3)') EPSET,THSET,Afromin,Zfromin,
c     >      sigexp(npts+1),sigexper(npts+1) 
        elseif(dosimona) then
          read(7,'(i2,3f9.4,27x,3f11.4)',end=100) itrgsim,
     >      E0SET,THSET,EPSET,
     >      sigexp(npts+1),sigexper(npts+1),rcsim
          if(abs(thset - thsimona(ifile)).gt.0.5) goto 99
          if(itrgsim.eq.11.and.IA.eq.2) goto 99
          if(itrgsim.eq.13.and.IA.eq.1) goto 99
! for H2, values are born, so need to divide by rc used
          if(itrgsim.eq.11) then
            CALL SECNUCLW(E0SET,EPSET,TH,SIGMA_BORN)                     
            write(99,'(1x,3f6.2,f10.3)') E0SET,EPSET,TH,
     >        sigexp(npts+1)/SIGMA_BORN,
     >        sigexper(npts+1)/SIGMA_BORN
            sigexp(npts+1) = sigexp(npts+1) / rcsim
            sigexper(npts+1) = sigexper(npts+1) / rcsim
          endif
        elseif(donucr) then
          read(7,*,end=100) E0SET,THSET,EPSET,
     >      sigexp(npts+1),sigexper(npts+1) 
        elseif(doeg1b) then
c          read(7,'(1x,2i3,2f6.3,F7.3)') ipsv,ithsv,
c     >      E0SET,EPSET,THSET
c          read(7,'(3i3,3F8.3)',end=100) iesv,iqsv,iwsv,
c     >      E0SET,EPSET,THSET
          read(7,*,end=100) iesv,iqsv,iwsv,
     >      E0SET,EPSET,THSET
          if(iesv.eq.0) goto 100
        else
c           READ(7,'(1X,2F6.3,F7.3)',END=100) E0SET,EPSET,THSET                
c           READ(7,'(1X,F6.3,2F8.3)',END=100) E0SET,EPSET,THSET                
           READ(7,*,END=100) E0SET,EPSET,THSET                
c           write(6,'(1x,3f10.3)') E0SET,EPSET,THSET
           IF (E0SET.LE.0.) GOTO 100                                    
        endif
C          CALL INTERPOL(E0SET,EPSET,THSET)                             
                                                                        
           TH    = THSET                                                
           THR   = TH*PI/180.                                           
           SIN2  = SIN(THR/2.)**2                                       
           Q2SET = 4.*E0SET*EPSET*SIN2                                  
           ANU   = E0SET-EPSET                                          
           X     = Q2SET/2./PM/ANU                                      
           Y     = ANU/E0SET                                            
           EPS   = 1./( 1.+2.*SIN2/(1.-SIN2)*(1.+Q2SET/(2.*PM*X)**2))   
           W2    = PM**2+2.*PM*ANU-Q2SET                                
ccc        if(dojan05.and.(w2 - w2last).lt.0.2 .and.w2.gt.1.3) goto 99
           if(dojan05.and.(w2.lt.1.1.or.w2.gt.4.0)) goto 99
           if(dofromin.and.(w2last-w2).lt.0.1) goto 99
cc           if(iA.lt.1.5.and.w2.le.1.17.and.(.not.smrelasp)) goto 99
c           if(w2.le.0.) goto 99
           if(IA.gt.2.and.x.gt.float(ia).and.(.not.smrelasp)) goto 99
           if(x.gt.3.0) goto 99
           w2last=w2
           q28 = q2set
           w8 = sqrt(max(0.,w2))
           npts = npts+1
           e0sv(npts)=e0set
           epsv(npts)=epset
           thsv(npts)=thset
           xsv(npts)=x
           w2sv(npts)=w2
           costh = cos(thr)
           E0 = E0set
           EP = EPSEt
           QV = SQRT(E0*E0+EP*EP-2.*E0*EP*COSTH) 
           AMS   = avgM-PM 
           AMF2  = PM**2
           A = 4.*QV**2-4.*(E0-EP+avgM)**2 
           BB = -4.*QV*(AMS**2-AMF2-QV**2+(E0-EP+avgM)**2) 
           C = (E0-EP+avgM)**4+(AMS**2-AMF2-QV**2)**2
     >       -2.*(E0-EP+avgM)**2*(AMF2+QV**2+AMS**2)
           DISC  = BB**2-4.*A*C
           ysv(npts)=0.
           IF (A.ne.0..and.DISC.gt.0.) then
             Ysv(npts) = (-BB-SQRT(DISC))/2./A
             Ysvp(npts) = (-BB+SQRT(DISC))/2./A
           endif
! Default resolutions
! if using grid, then resolution smearing in W2 is
! hardwired: see line 3160!
! good for BARONS SOS with .2% and 4 mr
!*** Th resolution for eA elastic: should use A* 0.98 in future
c          smrdEp = sqrt((ep * 0.002)**2 +
c    >      (0.005 * (2. * e0 * ep * sin(thr)) /
c    >               (2. * 0.9383 + q2set / ep))**2) 
! this is for HMS: 0.15% and 1 mr
c          smrdEp = sqrt((ep * 0.0015)**2 +
c    >      (0.001 * (2. * e0 * ep * sin(thr)) /
c    >               (2. * 0.9383 + q2set / ep))**2)
! Make smaller for test
           smrdEp = sqrt((ep * 0.0005)**2 +
     >      (0.0001 * (2. * e0 * ep * sin(thr)) /
     >               (2. * 0.9383 + q2set / ep))**2)
C Resolution in W at W=M for particular experiment
           IF(INDEX(TARGET,'EG1b').GT.0) THEN 
             smrdEp = ep * 0.02 / sqrt(e0) 
           ENDIF
           IF(INDEX(TARGET,'EG4').GT.0) THEN 
             smrdEp = ep * 0.02 / sqrt(e0) 
           ENDIF
C Born cross section:                                    
           prttst=.true.
           CALL SECNUCLW(E0SET,EPSET,TH,SIGMA_BORN)                     
c           write(98,'(''sig born'',e11.4)') sigma_born

           doing_elas = .false.
           CALL QUASIY8(E0SET,EPSET,TH,SIGMA_CORR)                     

! added smeared elastic peak born
           doing_elas = .true.
           sigma_elastic = 0.
           if(smrelasp) CALL QUASIY8(E0SET,EPSET,TH,SIGMA_elastic)
c           write(6,'(''elastic'',3f8.3,e11.3)') 
c     >       E0SET,EPSET,TH,SIGMA_elastic
           doing_elas = .false.

           prttst=.false.
C          CALL  FYXSEC8(E0SET,EPSET,TH,SIGMA_CORR)                     
           write(*,'(1x,''e,ep,th,sigb,sigqe='',3f8.2,2e12.4)')
     >       E0SET,EPSET,TH,SIGMA_BORN,sigma_corr
C Radiation length integral. Loop over internal and internal+external:  
! If eg1b type target, replace index 1 with results with ttarg/2.
           IF(INDEX(TARGET,'EG1b').GT.0 .or.
     >        INDEX(TARGET,'EG4').GT.0) THEN 
             JTARG = 2               
             ttarg = ttarg / 2.0
             SIGMA_INEL(1) =0.
! inelastic
             if(w2.gt.1.17.or.IA.gt.1) CALL SIMP(
     >          0.,TTARG,NSEG,CONTINUUM,SIGMA_INEL(1))             
! quasielastic
             CALL SIMP(0.,TTARG,NSEG,   QETAIL,SIGMA_QELA(1)) 
! elastic
             if(smrelasp)  then
               doing_elas = .true.
               CALL SIMP(0.,TTARG,NSEG,QETAIL,SIGMA_ELPK(1)) 
               doing_elas = .false.
             else
              IF(NUC_METHOD.EQ.0)  THEN
               CALL SIMP(0.,TTARG,NSEG,EPEAK,SIGMA_ELPK(1))
              ELSEIF(NUC_METHOD.EQ.1)  THEN
               CALL SIMP(0.,TTARG,NSEG,ATAILFL1,SIGMA_ELPK(1)) 
               SIGMA_ELPK(1) = 0.001 *SIGMA_ELPK(1) 
              ENDIF
             endif
             SIGMA_INEL(1) = SIGMA_INEL(1)/TTARG                       
             SIGMA_QELA(1) = SIGMA_QELA(1)/TTARG  
             SIGMA_ELPK(1) = SIGMA_ELPK(1)/TTARG
c             write(6,'(''tst 1'',f7.3,3e11.4)') ttarg,
c     >         sigma_elpk(1) , sigma_qela(1) , sigma_inel(1)
             ttarg = ttarg * 2.0
           ELSE
             JTARG = 1                                                    
             SIGMA_INEL(1) = 0.                               
! inelastic
             if(w2.gt.1.17.or.IA.gt.1) SIGMA_INEL(1) = CONTINUUM(0.)   
! quasielastic
             SIGMA_QELA(1) = QETAIL(0.)                                   
cc no longer used
             SIGMA_QEPK(1) = 0.
cc             SIGMA_QEPK(1) = QEPEAK(0.)
! Q-E using peak and cosk: no longer used
             SIGMA_QE_PK_COS(1) = 0.
cc             SIGMA_QE_PK_COS(1) = .001 *ATAILFL_QE(0.)
! elastic
             if(smrelasp)  then
               doing_elas = .true.
               SIGMA_ELPK(1) = QETAIL(0.)                                   
               doing_elas = .false.
             else
              IF(NUC_METHOD.EQ.0)  THEN
                SIGMA_ELPK(1) = EPEAK(0.)  !original external    
              ELSEIF(NUC_METHOD.EQ.1)  THEN
c              prttst=.true.
! E140 Mo-Tsai code (pb->nb)
               SIGMA_ELPK(1) = 0.001*ATAILFL1(0.) 
c              prttst=.false.
              ELSE
c               WRITE(6,'('' *** BAD parameter***  NUC_METHOD='',I2,
c     >           '' SHOULD BE 0 OR 1'')')NUC_METHOD 
                STOP
              ENDIF
             endif
           ENDIF

! now do with full radiation lengths
! (unless nseg=0)
           if(nseg.eq.0) then
             write(17,'(1x,5f7.3,6e11.4,3i3,4e11.4)') 
     >       e0set,epset,thset,xsv(npts),q2set,SIGMA_BORN,
     >       sigma_elpk(1)+sigma_qepk(1)+sigma_inel(1),
     >       sigma_elpk(1),sigma_qepk(1),sigma_inel(1)
             goto 99
           endif
           JTARG = 2               
! ineasltic
           SIGMA_INEL(2) =0.
           if(w2.gt.1.17.or.IA.gt.1) CALL SIMP(
     >         0.,TTARG,NSEG,CONTINUUM,SIGMA_INEL(2))             
! quasielastic
           CALL SIMP(0.,TTARG,NSEG,   QETAIL,SIGMA_QELA(2)) 
! skip this
           SIGMA_QE_PK_COS(2) = 0.
cc           CALL SIMP(0.,TTARG,NSEG,ATAILFL_QE,SIGMA_QE_PK_COS(2))
C PERM.    CALL SIMP(0.,TTARG,NSEG,   QEPEAK,SIGMA_QEPK(2)) 
! elastic
           if(smrelasp)  then
             doing_elas = .true.
             CALL SIMP(0.,TTARG,NSEG,   QETAIL,SIGMA_ELPK(2)) 
             doing_elas = .false.
           else
            IF(NUC_METHOD.EQ.0)  THEN
              CALL SIMP(0.,TTARG,NSEG,EPEAK,SIGMA_ELPK(2))  !original external
c              write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
            ELSEIF(NUC_METHOD.EQ.1)  THEN
              CALL SIMP(0.,TTARG,NSEG,ATAILFL1,SIGMA_ELPK(2)) 
              SIGMA_ELPK(2) = 0.001 *SIGMA_ELPK(2)   ! pb->nb 
c              write(*,*) '###########################################'
            ENDIF
           endif

cc no longer used
           SIGMA_QE_PK_COS(2) = 0.
cc           SIGMA_QE_PK_COS(2) = 0.001 *SIGMA_QE_PK_COS(2)  !pb->nb
           write(*,111) ttarg,SIGMA_INEL(2), SIGMA_QELA(2),
     >       SIGMA_QE_PK_COS(2), SIGMA_ELPK(2)
 111       format(1x,'ttarg,sig=',f8.4,4e11.3)
           SIGMA_INEL(2) = SIGMA_INEL(2)/TTARG                          
           SIGMA_QELA(2) = SIGMA_QELA(2)/TTARG  
cc           SIGMA_QE_PK_COS(2) = SIGMA_QE_PK_COS(2)/TTARG                          
           SIGMA_ELPK(2) = SIGMA_ELPK(2)/TTARG                          
c             write(6,'(''tst 2'',f7.3,3e11.4)') ttarg,
c     >         sigma_elpk(2) , sigma_qela(2) , sigma_inel(2)
! Total inelastic (DIS + Res. + q.e.) born xsection
           SIGMA_BORNPLSQE = SIGMA_BORN+SIGMA_CORR
           DO 98 I=1,2                                                  
!*** changed to divide by SIGMA_BORNPLSQE in all cases (before
! was only for QELA).
            DELTA_INEL(I) = (SIGMA_INEL(I)/ SIGMA_BORNPLSQE-1.)*100.         
            DELTA_QELA(I) = (SIGMA_QELA(I)/ SIGMA_BORNPLSQE)*100.
            DELTA_QEPK(I) = (SIGMA_QEPK(I)/ SIGMA_BORNPLSQE)*100. 
            DELTA_QE_PK_COS(I) = (SIGMA_QE_PK_COS(I)/ SIGMA_BORNPLSQE)*100.
            DELTA_ELPK(I) = (SIGMA_ELPK(I)/ SIGMA_BORNPLSQE)*100.            
            TOTAL(I)      =  DELTA_INEL(I)+DELTA_QELA(I)+DELTA_ELPK(I)  
 98        CONTINUE                                                     
                                                                        
           RAT_SMR = 1.                                                 
           IF(DELTA_QEPK(1).NE.0.0.and.w2.gt.0.95) 
     >       RAT_SMR = DELTA_QELA(1)/DELTA_QEPK(1)
                                                                        
c           CALL TIME(VTSTOP,TTSTOP)                                     
           VTRUN = VTSTOP-VTSTART                                       
           TTRUN = TTSTOP-TTSTART                                       
C          WRITE(11,'(1X,2F8.4)') TBSUM/12.,TASUM/12.                   
C          TASUM=0.                                                     
C          TBSUM=0.                                                     
c            WRITE(10,'(1X,2f7.3,f6.2,5F8.2,F7.3/
c    >          1x,2f7.3,F6.2,5f8.2,1PE10.3,/)')
c    >          E0SET, TH, Q2SET, DELTA_QE_PK_COS(1),  
c    >          DELTA_INEL(1),DELTA_QELA(1), DELTA_ELPK(1), 
c    >          TOTAL(1), RAT_SMR,EPSET, x, W2, 
c    >          DELTA_QE_PK_COS(2), DELTA_INEL(2),DELTA_QELA(2), 
c    >          DELTA_ELPK(2), TOTAL(2), SIGMA_BORNPLSQE  
             ddnom = SIGMA_BORN + SIGMA_CORR + sigma_elastic
             if(ddnom.ne.0.) WRITE(10,'(4f6.2,e10.3,9f6.2)')
     >          E0, EP, TH, w2, ddnom, 
     >          sigma_elastic/ddnom,  
     >          sigma_corr/ddnom,  
     >          sigma_born/ddnom,  
     >          sigma_elpk(2)/ddnom,
     >          sigma_qela(2)/ddnom,
     >          sigma_inel(2)/ddnom
      WRITE(*,'(1X,2f7.3,f6.2,5F8.2,F7.3/1x,2f7.3,F6.2,5f8.2,    
     >           1PE10.3,/)') 
     >    E0SET, TH, Q2SET, DELTA_QE_PK_COS(1),  DELTA_INEL(1),
     >       DELTA_QELA(1), DELTA_ELPK(1), TOTAL(1), RAT_SMR,
     >    EPSET,  x,    W2, DELTA_QE_PK_COS(2), DELTA_INEL(2),
     >       DELTA_QELA(2), DELTA_ELPK(2), TOTAL(2), SIGMA_BORNPLSQE  
      sbsv(npts)= SIGMA_BORNPLSQE
      sbisv(npts) = SIGMA_BORN
      sbqesv(npts) = SIGMA_CORR

      srsv(npts)= sigma_qela(2) + sigma_elpk(2) + sigma_inel(2)
      srqesv(npts) = sigma_qela(2)
      srelsv(npts) = sigma_elpk(2)
      srinsv(npts) = sigma_inel(2)

      ccor=1.0
c --- aji ----------------------------------------------
c now do Coulomb correction - should be set at top
      if(doccor) then 
cc first get the vertex quantities and do a target dependent
cc boost to beam energy and scattered electron energy
cc and do Coulomb correction.  recipe is given in:
cc Aste et al., Eur. Phys. J. A 26, 167 (2005)
cc  note that these deltae_cc's are calculated as follows:
cc deltae_cc = (1.5)*(Z/R)*(hbar*c)*alpha*0.775
cc the 0.775 comes from Aste EPJ A 26 167 (2005)
cc also, all deltae_cc's are computed for Z-1, not Z!/-* 
         
         call get_cc_info(iA,iZ,deltae_cc)
   
         E0SET_cc  = e0sv(npts)+ deltae_cc
         EPSET_cc = epsv(npts)+ deltae_cc
         CALL SECNUCLW(E0SET_cc,EPSET_cc,TH,SIGMABORN_CC)
         CALL QUASIY8(E0SET_cc,EPSET_cc,TH,SIGMACORR_CC)
         SIGTOT_CC = SIGMABORN_CC + SIGMACORR_CC

         f1cc = E0SET_cc/e0sv(npts)
         ccor=1.

         if(SIGTOT_CC.gt.0.0) then ! avoid nan
            ccor = sbsv(npts)/SIGTOT_CC/f1cc/f1cc
         else
            ccor=1.0
         endif
c------------------------------------------------------------------
c aji the born model with ccor, since we  radiated the model with  Coulombic
c effect for this version
c this means, the "Born" cross section in the output file actually has Coulomb corrections
c applied - you need to remove them if you want the true "Born" cross section
         sbsv(npts)= sbsv(npts)/ccor
         sbisv(npts)= sbisv(npts)/ccor
         sbqesv(npts)= sbqesv(npts)/ccor
      endif                     !doccor 


      if(srsv(npts).ne.0.) then 
         write(17,'(1x,5f9.4,8e13.5)') 
     >        e0sv(npts),epsv(npts),thsv(npts),
     >        xsv(npts),q2set,
     >        max(0.,sbsv(npts)),
     >        max(0.,sbisv(npts)),
     >        max(0.,sbqesv(npts)),
     >        max(0.,srsv(npts)),
     >        max(0.,srelsv(npts)),
     >        max(0.,srqesv(npts)),
     >        max(0.,srinsv(npts)),
     >        ccor
c     >        max(0.,sigma_qela(1) + sigma_elpk(1) + sigma_inel(1)),
c     >        max(0.,sigma_elastic)
      endif
c      write(18,'(1x,4f10.3)') e0sv(npts),epsv(npts),thsv(npts),
c     >  1./(1+total(2)/100.)
c      write(18,'(1x,f4.1,f6.2,f6.1,7e11.4)') 
c     >  e0sv(npts),epsv(npts),thsv(npts),
c     >  sbsv(npts),SIGMA_BORNPLSQE*(1+total(1)/100.),
c     >  sigma_elpk(1),sigma_inel(1),srsv(npts),
c     >  srelsv(npts),srinsv(npts)
 99   CONTINUE                                                          
100   CONTINUE                                                          
      if(dojan05.or.dofromin.or.dosimona.or.donucr) then
       ymax=0.
       do i=1,npts
         ymax = max(ymax, sigexp(i))
! if Fromin, plot versus y instead of W2
         write(44,'(1x,''w2,y'',3f8.3)') w2sv(i),ysv(i),
     >     ysvp(i)
         xplt(i) = sqrt(max(0.,w2sv(i)))
         if(dofromin) xplt(i) = ysv(i)
       enddo
       if(dojan05) ymin=ymax/100.
       if(donucr) ymin=ymax/10.
       if(dofromin) ymin=ymax/100000.
       if(dosimona) ymin=ymax/300.
       if(dofromin) thset=thfromin(ifile)
       if(dosimona) thset=thsimona(ifile)
       if(dojan05) thset = janth(ifile)
       if(dojan05) E0set = jane0(ifile)
c       open(unit=16,file='extern.top')
c       write(16,'(1x,''set device postscript portrait'')')
c       write(16,1626) IA,IZ,E0set,THSET,ymin,ymax*1.1
 1626  format(1x,'new frame'/
     >  1x,'title top ',1h','A=',i3,' Z=',i2,' E0=',f5.3,
     >    ' Th=',f4.1,1h'/
     >  1x,'set bar size 0. ; set order x y dy '/
     >  1x,'set scale y log ; set limits y  ',2e12.5/
     >  1x,'set sym 9O size 0.5'/
     >  1x,'title left ',1h','dsigma (nb/sr/gev)',1h')
c       write(19,1926) IA,IZ,E0set,THSET
c       write(29,1926) IA,IZ,E0set,THSET
 1926  format(1x,'new frame'/
     >  1x,'title top ',1h','A=',i3,' Z=',i2,' E0=',f5.3,
     >    ' Th=',f4.1,1h'/
     >  1x,'set bar size 0. ; set order x y dy '/
     >  1x,'set scale y log ; set limits y  0.3 3.0'/
     >  1x,'set sym 9O size 2.0'/
     >  1x,'title left ',1h','data/model',1h')
       if(dofromin) then
c        write(16,1636)
 1636    format(1x,'title bottom ',1h','y',1h')
       else
c         write(16,1637)
 1637    format(1x,'title bottom ',1h','W2 (GeV)',1h')
       endif
c       write(19,1636)
c       write(29,1637)
c       write(16,'(1x,f8.4,2e12.4)') (xplt(i),
c     >   sigexp(i),sigexper(i),i=1,npts)
c       write(16,'(1x,''plot ; set symbol 9O size 1.0'')')
c       write(16,'(1x,''plot ; set symbol 9O size 1.5'')')
c       write(16,'(1x,''plot ; set symbol 9O size 2.0'')')
       fact = 1000. * float(IA)
       if(dofromin) fact=1000.
       if(donucr) fact=1000.
       if(dosimona) fact=1.
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srsv(i)/fact,i=1,npts)
c       write(16,'(1x,''plot ; set symbol 1O size 2.0'')')
       if(dofromin) write(16,'(1x,''join'')')
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srelsv(i)/fact,i=1,npts)
c       write(16,'(1x,''plot ; set symbol 2O size 2.0'')')
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srqesv(i)/fact,i=1,npts)
c       write(16,'(1x,''plot ; set symbol 3O size 2.0'')')
c       write(16,'(1x,f8.4,e12.4)') (xplt(i),srinsv(i)/fact,i=1,npts)
c      write(16,'(1x,''plot ; set symbol 4O size 2.0'')')
       do i=1,npts
         if(srsv(i).ne.0.) then
c           write(19,'(1x,f8.4,2e12.4)') ysv(i),
c     >       max(0.3,min(3.,
c     >       sigexp(i)/srsv(i)*fact)),
c     >       sigexper(i)/srsv(i)*fact
c           write(29,'(1x,f8.4,2e12.4)') w2sv(i),
c     >       max(0.3,min(3.,
c     >       sigexp(i)/srsv(i)*fact)),
c     >       sigexper(i)/srsv(i)*fact
         endif
       enddo
c       write(19,'(1x,''plot ; set symbol 9O size 1.0'')')
c       write(19,'(1x,''plot ; set symbol 9O size 1.5'')')
c       write(19,'(1x,''plot ; set symbol 9O size 2.0'')')
c       write(29,'(1x,''plot ; set symbol 9O size 1.0'')')
c       write(29,'(1x,''plot ; set symbol 9O size 1.5'')')
c       write(29,'(1x,''plot ; set symbol 9O size 2.0'')')
c       write(19,'(1x,''-10. 1. ; 100. 1. ; join dash'')')
c       write(29,'(1x,''-10. 1. ; 100. 1. ; join dash'')')
      else
c       write(16,1616)
c 1616  format(1x,'title bottom ',1h','W**2 (GeV)**2',1h'/
c     >  1x,'title left ',1h','dsigma (nb/sr/gev)',1h')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),sbsv(i),i=1,npts)
c       write(16,'(1x,''join'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srsv(i),i=1,npts)
c       write(16,'(1x,''join'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srelsv(i),i=1,npts)
c       write(16,'(1x,''join dot'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srqesv(i),i=1,npts)
c       write(16,'(1x,''join dash'')')
c       write(16,'(1x,f8.4,e12.4)') (w2sv(i),srinsv(i),i=1,npts)
c       write(16,'(1x,''join dotdash'')')
      endif
      WRITE(10,'(1X,2F6.3,F7.3,1X,F4.3,2F6.2,4F6.2,F6.3,/,25X,6F6.2,    
     >           1PE10.3,/)') (ZERO(i),i=1,18)                          

! loop over files
      enddo 
      CLOSE(10)
                                                                        
      END
                                                               
C=======================================================================
                                                                        
      SUBROUTINE SIMP(A,B,N,FUNC,ANS)                                   
                                                                        
      REAL LOW                                                          
                                                                        
C Subroutine to do numerical integration using Simpson's rule           
C Input A lower limit                                                   
C Input B upper limit                                                   
C Input N number of bins; FUNC is called N+1 times; N HAS TO BE EVEN    
C Input FUNC Function subroutine address                                
C Output ANS estimate of the definite integral                          
                                                                        
      ANS = 0.                                                          
      IF (A.GE.B) RETURN                                                
      IF (N.LE.0) RETURN                                                
      IF (N.EQ.1) THEN                                                  
           ANS = FUNC((A+B)/2.)*(B-A)                                   
           RETURN                                                       
      ENDIF                                                             
      BINW = (B-A)/FLOAT(N)                                             
                                                                        
C FUNCTION AT LOWER LIMIT                                               
                                                                        
      LOW = FUNC(A)                                                     
                                                                        
C FUNCTION AT ODD BINS                                                  
                                                                        
      ODD = 0.                                                          
      DO 100 I = 1,N-1,2                                                
           ODD = ODD+FUNC(A+BINW*FLOAT(I))                              
100   CONTINUE                                                          
                                                                        
C FUNCTION AT EVEN BINS                                                 
                                                                        
      EVEN = 0.                                                         
      DO 200 I = 2,N-2,2                                                
           EVEN = EVEN+FUNC(A+BINW*FLOAT(I))                            
200   CONTINUE                                                          
                                                                        
C FUNCTION AT UPPER LIMIT                                               
                                                                        
      UP = FUNC(B)                                                      
                                                                        
C ANSWER                                                                
                                                                        
      ANS = (LOW+4.*ODD+2.*EVEN+UP)*BINW/3.                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE SIMP2(A,B,N,FUNC,ANS)                                  
                                                                        
      REAL LOW                                                          
      ANS = 0.                                                          
      IF (A.GE.B) RETURN                                                
      IF (N.LE.0) RETURN                                                
      IF (N.EQ.1) THEN                                                  
           ANS = FUNC((A+B)/2.)*(B-A)                                   
           RETURN                                                       
      ENDIF                                                             
      BINW = (B-A)/FLOAT(N)                                             
      LOW = FUNC(A)                                                     
      ODD = 0.                                                          
      DO 100 I = 1,N-1,2                                                
           ODD = ODD+FUNC(A+BINW*FLOAT(I))                              
100   CONTINUE                                                          
      EVEN = 0.                                                         
      DO 200 I = 2,N-2,2                                                
           EVEN = EVEN+FUNC(A+BINW*FLOAT(I))                            
200   CONTINUE                                                          
      UP = FUNC(B)                                                      
      ANS = (LOW+4.*ODD+2.*EVEN+UP)*BINW/3.                             
      RETURN                                                            
      END                                                               
                                                                        
                                                                        
C=======================================================================
                                                                        
      SUBROUTINE INTERPOL(E0,EP,TH)                                     
                                                                        
      CHARACTER*1 STAR                                                  
      DIMENSION ZERO(18)                                                
      DATA ZERO /18*0./                                                 
                                                                        
100   READ(7,'(1X,2F6.3,F7.3,1X,F4.3,2F6.2,F6.3,1X,A1)')                
     >         E0,EP,TH,X,Q2,W2,TASI,STAR                               
                                                                        
      IF (STAR.NE.'*'.AND.STAR.NE.'X'.AND.STAR.NE.'Q') THEN             
           CALL SECNUCLW(E0,EP,TH,SLAC)                                 
           WRITE(10,'(1X,2F6.3,F7.3,1X,F4.3,2F6.2,4F6.2,F6.3,/,25X,     
     >           6F6.2,1PE10.3,/)') E0,EP,TH,(ZERO(i),i=1,14),SLAC      
           GOTO 100                                                     
      ELSE                                                              
           RETURN                                                       
      ENDIF                                                             
      END                                                               
                                                                        
C=======================================================================
                                                                        
                                                                        
C=======================================================================
                                                                        
      subroutine weiz                                                   
                                                                        
C Calculates weizsacker mass formula from Segre iff avgM not supplied by
C READIN. Accuracy = +/-.002 GeV.                                       

      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
                                                                        
      avgA = iA                                                         
      if (iA.le.5) then                                                 
           if (iA.eq.1)             avgM =  .93828                      
           if (iA.eq.2)             avgM = 1.87537                      
ccc pyb to try to make elastic peak work
cc           if (iA.eq.2)             avgM = 2.0
                      
           if (iA.eq.3.and.iZ.eq.1) avgM = 2.8095                       
           if (iA.eq.3.and.iZ.eq.2) avgM = 2.8094                       
           if (iA.eq.4)             avgM = 3.7284                       
           return                                                       
      endif                                                             
                                                                        
      extra = 0.                                                        
      iN    = iA-iZ                                                     
      iN2   = iN/2*2                                                    
      iZ2   = iZ/2*2                                                    
      if (iN2.eq.iN.and.iZ2.eq.iZ) extra = -.012*(iA**(-.5))            
      if (iN2.ne.iN.and.iZ2.ne.iZ) extra =  .012*(iA**(-.5))            
                                                                        
      avgM = iN*.939573 + iZ*(.938280+.000511) - iA*.01567 +            
     +       .01723*(iA**(2./3.)) + .09315*((.5*iA-iZ)**2.)/iA +        
     +       .0006965*(iZ**2.)*(iA**(-1./3.)) + extra                   
                                                                        
      return                                                            
      end                                                               
                                                                        
C ----------------------------------------------------------------------
c                                                                       
c     subroutine time(vt,tt)                                            
c                                                                       
c      call vttime(ivtime,ittime)                                       0
c     vt = float(ivtime)/100.0                                          
c     tt = float(ittime)/100.0                                          
c     return                                                            
c     end                                                               
                                                                        
C ====================================================================  
                                                                        
       SUBROUTINE ASYM_INE(E0,EP,TANSQ,QSQ,TARGET,SIGMA)                
!-------------------------------------------------------------------    
! Modifies polarized cross section because of E142 asymmetry.           
! A1 parameterization from Emlyn on 3/16/93                             
! Crude numbers for beam and target polarization and dilution factor    
!-------------------------------------------------------------------    
      IMPLICIT NONE                                                     
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER iZ,iA                                                     
      REAL avgN,avgA,avgM,amuM                                          
      REAL E0,EP,TANSQ,QSQ,SIGMA,X,EPS,D,NU,A1/0./,APAR,DELTA           
      REAL N_OVER_P,TARG_DILUTE/0./,A1N/0./,A1P/0./                     
      REAL R/.2/,W_D/.058/ !deuteron d state                            
      CHARACTER*7 TARGET                                                
      LOGICAL FIRST/.TRUE. /,GOOD_TARG/.FALSE./                         
                                                                        
      NU = E0 -EP                                                       
      X = QSQ/(2.*.93828 *NU)                                            
      EPS = 1./( 1. +2.*(1.+ NU**2/QSQ)*TANSQ  )                        
      D =( 1. -EP/E0 *EPS )/(1. +EPS*R)                                 
      N_OVER_P = 1.- .8*X     !ratio of cross sections                  
      IF(INDEX(TARGET,'E142').GT.0) THEN!  Neutron A1                   
       A1 = 1.58 * X**1.4 * (1.- X) + X * (2.5476 * X - 1.5476)         
! We are using He3 cross section, so asymmetry is diluted.              
       TARG_DILUTE = 1. + 2./N_OVER_P  ! ratio of He to neutron sigma   
       GOOD_TARG =.TRUE.                                                
       IF(FIRST) WRITE(6,'('' NEUT TARG'')')                            
      ELSE IF(INDEX(TARGET,'E143').GT.0)THEN ! NH3 target, proton asym  
       A1 = 1.025*(X**.12) * (1.-EXP(-2.7*X))  ! EMC parameterization   
       A1P=A1                                                           
       IF(IA.EQ.1) THEN ! proton target (with N in NH3 as radiator)     
        TARG_DILUTE =1.                                                 
       ELSEIF(IA.EQ.2) THEN !deuteron target                            
        A1N = 1.58 * X**1.4 * (1.- X) + X * (2.5476 * X - 1.5476)       
        A1 = (A1P/(1.+N_OVER_P) +A1N *(N_OVER_P)/(1.+N_OVER_P) )*       
     >    (1.- 1.5*W_D)                                                 
        TARG_DILUTE =1.                                                 
       ELSE  ! Perhaps weird case using NH3 cross section??             
        TARG_DILUTE = (10. + 8.*N_OVER_P)/3. ! NH3/H sigma ratio        
       ENDIF                                                            
       GOOD_TARG =.TRUE.                                                
      ELSEIF(INDEX(TARGET,'E149').GT.0) THEN !Parity violation          
c      IF(FIRST) THEN                                                   
c       WRITE(10,'('' Asym is 100 times parity violation'')')           
c       WRITE( 6,'('' Asym is 100 times parity violation'')')           
c      ENDIF                                                            
       GOOD_TARG=.TRUE.                                                 
       D=1.0
       A1 =  0.8E-4 * QSQ/D         !Fake a Sigma = delta sigma         
c      A1 = 100. * 0.80E-4 *QSQ/D   ! 100 Times parity violation        
       TARG_DILUTE =1.                                                  
      ENDIF                                                             
      APAR = D * A1/TARG_DILUTE !change in He or NH3 cross section      
! alter neutron xsection by differences.                                
      DELTA =     (1.+ APAR)/(1.-APAR)                                  
! If model is sigma parallel, then we calculate sigma anti-parallel     
! Measured asymmetry = antiparrallel - parallel                         
c      SIGMA = SIGMA* DELTA                                              
      SIGMA = SIGMA * (1. + apar)                                             
                                                                        
      IF(FIRST.AND.GOOD_TARG) THEN                                      
       FIRST =.FALSE.                                                   
       WRITE(6 ,'('' sig chngd dueto asym; D,A1,DELT,dilut='',          
     >   F5.2,1P,2E8.1,0P,F6.3,A8)')                                    
     >  D,A1,DELTA,TARG_DILUTE,TARGET                                   
       WRITE(10,'('' sig chngd dueto asym; D,A1,DELT,DILUT='',          
     >   f5.2,3f8.4,A8)')                                    
     >  D,A1,DELTA,TARG_DILUTE,TARGET                                   
       ENDIF                                                            
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
!====================================================================== 
                                                                        
      SUBROUTINE ASYM_QEL(E0,EP,THETAR,QSQ,TARGET,GEN,GMN,GEP,GMP,      
     > MOTR,SIGMA)                                                      
!--------------------------------------------------------------         
! Calculate the Difference in cross sections in quasielastic            
! scattering from polarized Helium.to                                   
! From  Blankleider and Woloshyn, PR C29, 538 (1984)                    
!                                  -Steve Rock  4/13/93                 
!---------------------------------------------------------------        
      IMPLICIT NONE                                                     
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
      INTEGER iZ,iA                                                     
      REAL avgN,avgA,avgM,amuM                                          
      REAL E0,EP,THETAR,QSQ,GEN,GMN,GEP,GMP,SIGMA,GE/0./,GM/0./         
      REAL MOTR ! Mott * recoil factor                                  
      REAL ASYM/0./!Fractional change in cross section(phony initialize 
      CHARACTER*7 TARGET                                                
      REAL TARG_DILUTE/0./,TAN2,SIGN,SIGHE,SIGP,ASYM_P,ASYM_N           
      REAL G1/0./,G2/0./,TAU                                    
      REAL BETA /0./  ! angle between beam and Helium spins             
      REAL COSB/1./                                                     
      REAL MP2/.8803/,MP/.93828/                                           
      REAL GEX,GMX,ASYM_F,FDEL                                          
      INTEGER NTIMES/1/,TIMES/0/                                        
                                                                        
! M.J. Aluard,Phys. Rev. Lett 37, 1258(76) (slac experiment)            
! GE is only in numerator                                               
      ASYM_F(GEX,GMX) = TAU*GM *     !statement function                
     > (2.*MP*GEX/E0 +GMX*(2.*TAU*MP/E0 +2.*(1.+TAU)*TAN2))             
     >  /(GEX**2 + TAU*(GMX)**2 *(1. +2.*(1.+TAU)*TAN2))                
                                                                        
! If model is sigma parallel, then we calculate sigma anti-parallel     
! Measured asymmetry = antiparrallel - parallel                         
                                                                        
      TAU = QSQ/(4.*MP2)                                                
      TAN2 = (TAN(THETAR/2.))**2                                        
      IF(( INDEX(TARGET,'E142')+INDEX(TARGET,'E143')).GT.0) THEN        
       IF(INDEX(TARGET,'E142').GT.0) THEN !neutron asymetry             
C       G1 = - GMN *(GEN+ TAU* GMN)/(2.*(1.+TAU) )                      
C       G2 =   GMN * (GMN -GEN)/(4.*(1.+TAU) )                          
        GE = GEN                                                        
        GM = GMN                                                        
        SIGN = (GEN**2 +TAU*GMN**2)/(1.+TAU) +2.*TAU*GMN**2*TAN2        
        SIGHE=  (GEN**2 +2.*GEP**2 +TAU*(GMN**2+2.*GMP**2))/(1.+TAU) +  
     >   2.*TAU*(GMN**2 +2.*GMP**2)*TAN2                                
        TARG_DILUTE= SIGN/SIGHE                                         
        ASYM =ASYM_F(GE,GM)                                             
        SIGMA = SIGMA * (1.+ TARG_DILUTE*ASYM)/(1.-TARG_DILUTE*ASYM)    
       ELSEIF(INDEX(TARGET,'E143').GT.0) THEN !proton asym, NH3 target  
        IF(IA.EQ.1) THEN ! proton                                       
         GE = GEP                                                       
         GM = GMP                                                       
         TARG_DILUTE =1.                                                
         ASYM =ASYM_F(GE,GM)                                            
         SIGMA = SIGMA * (1.+ TARG_DILUTE*ASYM)/(1.-TARG_DILUTE*ASYM)   
        ELSEIF(IA.EQ.2) THEN ! deuterium                                
         SIGN = (GEN**2 +TAU*GMN**2)/(1.+TAU) +2.*TAU*GMN**2*TAN2       
         SIGP =((1.-FDEL(QSQ))*GEP**2 +TAU*GMP**2)/(1.+TAU)             
     >          +2.*TAU*GMP**2*TAN2                                     
         ASYM_P =ASYM_F(GEP,GMP)                                        
         ASYM_N =ASYM_F(GEN,GMN)                                        
         SIGMA =(SIGN *(1.+TARG_DILUTE*ASYM_N)/(1.-TARG_DILUTE*ASYM_N)  
     >        + SIGP *(1.+ TARG_DILUTE*ASYM_P)/(1.-TARG_DILUTE*ASYM_P) )
     >        * MOTR                                                    
        ENDIF                                                           
       ENDIF                                                            
C       ASYM=  TAU*GM *                                                 
C     > (2.*MP*GE/E0 +GM*(2.*TAU*MP/E0 +2.*(1.+TAU)*TAN2))              
C     >  /(GE**2 + TAU*(GM)**2 *(1. +2.*(1.+TAU)*TAN2))                 
      ELSEIF(INDEX(TARGET,'E149').GT.0) THEN  !parity violation         
       IF(TIMES.EQ.0) THEN                                              
        WRITE(10,'('' Asym is 100 times parity violtation=0.6E-4'')')   
        WRITE( 6,'('' Asymm is 100 times parity violation=0.6E-4'')')   
       ENDIF                                                            
C      ASYM = (0.6E-4 *QSQ -1.)/2.!Fake ASYM to make sigma = delta sig  
       ASYM = 100. *0.6E-4*QSQ    ! 100 Times parity violation.         
       TARG_DILUTE=1.                                                   
       SIGMA = SIGMA * (1.+ TARG_DILUTE*ASYM)/(1.-TARG_DILUTE*ASYM)     
      ENDIF                                                             
                                                                        
      IF(TIMES.LT.NTIMES) THEN                                          
        TIMES = TIMES + 1                                               
        WRITE(6 ,'('' Quasi chnged  asym;A,DILUTE,THR,E0,EP,Q2='',      
     >   1P,E8.1,0P,F5.2,F4.1,2F5.1,F6.2)')                             
     >   ASYM,TARG_DILUTE,THETAR*57.296,E0,EP,QSQ                       
        WRITE(10,'('' Sig Quasi chnged dueto asym: M.J.Aluard formula''/
     >  '' A,DILUTE,THR,E0,EP,Q2='',                                    
     >   1P,E8.1,0P,F5.2,F4.1,2F5.1,F6.2)')                             
     >   ASYM,TARG_DILUTE,THETAR*57.296,E0,EP,QSQ                       
       ENDIF                                                            
                                                                        
      RETURN                                                            
      END                                                               

C $MEMBER=QUADMO, DATE=75061922, USER=KCE                               
      REAL*4 FUNCTION QUADMO_R(FUNCT,LOWER,UPPER,EPSLON,NLVL)  
      REAL*4 FUNCT,LOWER,UPPER,EPSLON                                   
      INTEGER NLVL                                                      
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*4  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS


         IF(LOWER.EQ.UPPER) THEN
          QUADMO_R =0.
          RETURN
         ENDIF

         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)
         FM = FUNCT(.5*(L+R)) 
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=ABS(ESTL)+ABS(ESTR)        
      ABAREA=AREA+ABAREA-ABS(EST)     
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((ABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.     
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO_R = ROMBRG /12.0D0                                         
      RETURN                                                            
      END                                                               

C $MEMBER=QUADMO, DATE=75061922, USER=KCE                               
      REAL*4 FUNCTION QUADMO_U(FUNCT,LOWER,UPPER,EPSLON,NLVL)  
      REAL*4 FUNCT,LOWER,UPPER,EPSLON                                   
      INTEGER NLVL                                                      
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*4  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS
   
         IF(LOWER.EQ.UPPER) THEN
          QUADMO_U =0.
          RETURN
         ENDIF

         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)
         FM = FUNCT(.5*(L+R)) 
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=ABS(ESTL)+ABS(ESTR)        
      ABAREA=AREA+ABAREA-ABS(EST)     
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((ABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.     
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO_U = ROMBRG /12.0D0                                         
      RETURN                                                            
      END                                                               




!=========================================================================

      SUBROUTINE QE_PEAK(QSQ,E0,W1,W2)
!------------------------------------------------------
! E0 is a guess of the incident energy used for the vanOrden calculation
!  of Pauli Suppression.  It is calculated from Es- OM
!   where OM is the photon energy.
! Calculate nucleon W1 and W2 from Q2
!--------------------------------------------------------
      IMPLICIT NONE
      REAL*4 QSQ,E0,W1,W2
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      COMMON/TARGT/ iZ,iA,avgN,avgA,avgM,amuM                       
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM
      REAL TAU,PAULI_SUP1,PAULI_SUP2,PM,PM24
      REAL*8 GEP,GEN,GMP,GMN
      PARAMETER (PM    = 0.93828, PM24=3.5216)    



      TAU = QSQ/PM24                                               
      CALL NFORM(IG,DBLE(QSQ),GEP,GEN,GMP,GMN)


! Pauli suppression model
      CALL PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     > PAULI_SUP1,PAULI_SUP2)

      W1 =  Pauli_sup1* TAU * (iZ*GMP**2+avgN*GMN**2)                           
      W2 = (Pauli_sup2*(iZ*GEP**2+avgN*GEN**2) + W1)/(1.+TAU)  

      RETURN
      END


      SUBROUTINE PAULI_SUPPRESSION(PAULI_MODEL,NUC_MODEL,iA,QSQ,E0,
     > PAULI_SUP1,PAULI_SUP2)
!-----------------------------------------------------------------------
! Gets Pauli Suppression factor for quasi-elastic scattering from
! Several different Models.
! The VanOrden Model is very Slow.
! Used by both INTERNAL and EXTERNAL
!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*4 QSQ,E0,PAULI_SUP1,PAULI_SUP2,FDEL
      REAL Q,TAU,FSHELL,FGAUSS,FF_BESSEL,XX
      INTEGER PAULI_MODEL,NUC_MODEL,iA
      REAL FF
      REAL*4 MP24/3.52/   
      LOGICAL OUT_OF_RANGE
      REAL*4 P_FERMI/.25/

      PAULI_SUP1 = 1.
      PAULI_SUP2 = 1.
      if(IA.eq.1) return

      IF(PAULI_MODEL.EQ.0) THEN
        FF = 0.                                                           
        IF (iA.EQ.2.AND.QSQ.LE.8.)  FF = FDEL(QSQ) 
        IF(IA.GT.2) THEN
         OUT_OF_RANGE =.FALSE.
         IF(NUC_MODEL.EQ.1) THEN
          FF = FF_BESSEL(REAL(QSQ),OUT_OF_RANGE) !only for some Nuclei,feor limited Q2
         ENDIF
         IF(OUT_OF_RANGE.OR.(NUC_MODEL.EQ.0)) THEN !use if FF_BESSEL out of range
          IF (iA.LE.20) THEN
            FF  = FSHELL(QSQ)
          ELSE     !ia >20
            FF  = FGAUSS(QSQ) 
          ENDIF
         ENDIF  
        ENDIF                           
        Pauli_sup2 = (1.-FF**2)    !old model from Stein
        PAULI_sup1 =1.
      ELSE
       IF(PAULI_MODEL.EQ.1) THEN !Tsai RMP 46,816(74) eq.B54
           TAU   = QSQ/MP24     
           Q=SQRT(QSQ*TAU+QSQ)
           IF((Q .GT. 2.*P_FERMI).OR.(iA.EQ.1)) THEN
              PAULI_SUP2 =1.0
           ELSE
              PAULI_SUP2=3.0*Q*(1.0-0.08333*(Q/P_FERMI)**2)
           ENDIF
       ELSEIF(PAULI_MODEL.EQ.3) THEN
           CALL  Q_E_VANORDEN(QSQ,E0,PAULI_SUP2)
       ELSEIF(PAULI_MODEL.EQ.2) THEN ! no suppression
           PAULI_SUP2 =1.0   ! No suppression
       ELSE
           PAULI_SUP2 =1.0   ! No suppression
       ENDIF
        PAULI_SUP1= Pauli_sup2
      ENDIF
      IF(PAULI_MODEL.EQ.4) THEN ! Louk's formula for Deuterium
        TAU   = QSQ/MP24     
        Q=SQRT(QSQ*TAU+QSQ)
        xx = q * 1000. / 100.
        Pauli_sup2 = (1.  - 1. / (1. + xx * (0.0060 + xx * (0.3882 + 
     >    xx * 0.2477))))
        PAULI_SUP1= Pauli_sup2
      endif
           
      RETURN
      END


C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Quadratic Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming a quadratic
C       dependence:  Y = AX^2 + BX + C
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERPQ(X,Y,X0,NBINS)
      IMPLICIT NONE
C
      DOUBLE PRECISION X(100),Y(100),X0,DET,A,B,C
      INTEGER NBINS,I1/0/,I2/0/,I3/0/,I
C
C$$      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
c$$     #    WRITE(10,'('' Extrapolating outside range: X='',G10.3)') X0
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
1     DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y(I1)*(X(I2)-X(I3)) - X(I1)*(Y(I2)-Y(I3)) + Y(I2)*X(I3)
     #        - X(I2)*Y(I3) )/DET
      B = ( X(I1)**2*(Y(I2)-Y(I3)) - Y(I1)*(X(I2)**2-X(I3)**2)
     #        + X(I2)**2*Y(I3) - X(I3)**2*Y(I2) )/DET
      C = ( X(I1)**2*(X(I2)*Y(I3)-X(I3)*Y(I2))
     #        - X(I1)*(X(I2)**2*Y(I3)-X(I3)**2*Y(I2))
     #        + Y(I1)*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C
      RINTERPQ = A*X0**2 + B*X0 + C
C
      RETURN
      END
 
 
 
 
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     Exponential Interpolation Routine:
C
C       Calculates Y(X0) given array Y(X) assuming the exponential
C       form:  Y = EXP(AX^2 + BX + C)
C
C       This routine assumes the X values are arranged in
C       ascending order but does not assume a uniform X spacing
C       of the array.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERPEXP(X,Y,X0,NBINS)
      IMPLICIT NONE
C
      DOUBLE PRECISION X(1),Y(1),X0,DET,A,B,C,Y1,Y2,Y3
      INTEGER NBINS,I1/0/,I2/0/,I3/0/,I
C
C$$      IF(X0 .LT. X(1) .OR. X0 .GT. X(NBINS))
C$$     #    WRITE(10,'('' Extrapolating outside range: X='',G10.3)') X0
C
      IF(X0 .LE. X(1)) THEN
         I1 = 1
         I2 = 2
         I3 = 3
      ELSEIF(X0 .GE. X(NBINS-1)) THEN
         I1 = NBINS-2
         I2 = NBINS-1
         I3 = NBINS
      ELSE
         DO I=2,NBINS-1
            IF(X0 .LE. X(I)) THEN
               I1 = I-1
               I2 = I
               I3 = I+1
               GOTO 1
            ENDIF
         ENDDO
      ENDIF
C
C ----------------------------------------------------------------------
C     If all three Y-values are > 0, perform quadratic interpolation
C     on their logarithms; otherwise return zero as the result.
C ----------------------------------------------------------------------
C
1     IF(Y(I1).GT.0.D0 .AND. Y(I2).GT.0.D0 .AND. Y(I3).GT.0.D0) THEN
         Y1 = LOG(Y(I1))
         Y2 = LOG(Y(I2))
         Y3 = LOG(Y(I3))
      ELSE
 
         WRITE(10,'('' RiNTERPEXP:non-positive y-value; set to 0'')')
         RINTERPEXP = 0.D0
         RETURN
      ENDIF
C
      DET = (X(I2)-X(I3))*(X(I1)**2 - X(I1)*(X(I2)+X(I3)) + X(I2)*X(I3))
      A = ( Y1*(X(I2)-X(I3)) - X(I1)*(Y2-Y3) + Y2*X(I3)
     #        - X(I2)*Y3 )/DET
      B = ( X(I1)**2*(Y2-Y3) - Y1*(X(I2)**2-X(I3)**2)
     #        + X(I2)**2*Y3 - X(I3)**2*Y2 )/DET
      C = ( X(I1)**2*(X(I2)*Y3-X(I3)*Y2)
     #        - X(I1)*(X(I2)**2*Y3-X(I3)**2*Y2)
     #        + Y1*(X(I2)**2*X(I3)-X(I3)**2*X(I2)) )/DET
C
      RINTERPEXP = EXP(A*X0**2 + B*X0 + C)
C
      RETURN
      END
 
 
 
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C     2-Dimensional Linear Interpolation Routine:
C
C       Calculates F(X0,Y0) given array F(X,Y) by two successive
C       interpolations, first in X and then in Y.
C
C       Assumes uniform spacing of array in X and Y.
C------------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RINTERP2D(X,Y,F,X0,Y0,NX,NY,DELX,DELY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200),Y(200),F(200,200)
C
      I = DINT((X0-X(1))/DELX) + 1
      J = DINT((Y0-Y(1))/DELY) + 1
      IF((I+1.GT.NX).OR.(I.LT.1).OR.(J+1.GT.NY).OR.(J.LT.1))THEN
        RINTERP2D = 0.D0
        RETURN
      ENDIF
C
      RINTX1 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J)-F(I,J))
     #                  + F(I,J)
      RINTX2 = ((X0-X(I))/(X(I+1)-X(I)))*(F(I+1,J+1)-F(I,J+1))
     #                  + F(I,J+1)
      RINTERP2D = ((Y0-Y(J))/(Y(J+1)-Y(J)))*(RINTX2-RINTX1) + RINTX1
C
      RETURN
      END

                          

                          
!-----------------------------------------------------------------------
      SUBROUTINE Q_E_VANORDEN(Q2_ELAS_GEV,E_GEV,SUPPRESSION)
!---------------------------------------------------------------------------
C   This program compute the quasi-elastic cross section
C   based on the Van Orden calculation using the fermi gas model
!  input energy now in GeV
! It returns the Suppression factor for the quasi-elastic peak.
!-----------------------------------------------------------------------------

      IMPLICIT NONE
      REAL E_GEV,Z,A,KF_GEV,SUPPRESSION
      REAL SUM,AP,FMT,FACTOR,N,MT,TH,SINSQ,FMOTT,WMOTT,TANSQ,EP_ELAS,
     >  EF,RLR,RTR,SRL,RATIO,W,WBAR,QV,Q2BAR,F2,XETA,ZETA,T1,T2,E1,
     >  E,SL,ST,RS,RSSIG,XSECT,CROSS,QT,C,Q2_ELAS,E2,W0,W1,Q2,EINC,KF,
     >  KF3,Z1,A1,QV2,WSQ,MT4PI,MP2,Q2_ELAS_GEV
      INTEGER IOPT,I0,I00,I,IOPT1
      REAL MUP2/7.784/,MUN2/3.65/
      REAL EBAR/1./,WDEL/2./ !$$ are these OK?
      REAL  ALPH/.0072972/, PI/3.1415927/,AMASS/931.5/,MP/939./
      LOGICAL FIRST/.TRUE./
!-----------------------------------------------------------------------------------------
! In this notation, W=nu
!---------------------------------------------------------------------------------------     
      EINC = 1000.*E_GEV
      Q2_ELAS = 1.E6 * Q2_ELAS_GEV

C$$      OPEN(7,FILE='Q_E_VANORDEN2')
      SUM=0.

      IF(IOPT.EQ.1)THEN
 1     EP_ELAS = EINC - Q2_ELAS/(2.*MP)
       IF(EP_ELAS.GT.0) THEN
        SINSQ = Q2_ELAS/(4.*EINC*EP_ELAS)
        TH = 2.*ASIN(SQRT(SINSQ))
       ELSE
        EINC = EINC +5.
        GO TO 1
       ENDIF

       FMOTT=ALPH*ALPH*COS(TH/2.)**2/4./SINSQ/SINSQ*197.3**2*1.E-26
       WMOTT=FMOTT/EINC**2
       TANSQ=TAN(TH/2.)**2
       QT = sqrt(Q2_ELAS**2/(4.*MP2) + Q2_ELAS)
       IF(QT.GT. 2.*KF) THEN
        SUPPRESSION=1.
        RETURN
       ENDIF
       W0=MAX(EINC-(EP_ELAS+2.*KF),2.)
       W1=MAX(EINC-(EP_ELAS-2.*KF),0.)
C$$       WRITE(7,'(/''E_INC,THET,WDEL,FERMI_MOM='',1P9E10.2)')
C$$     >   EINC,THETA,WDEL,KF
C$$       WRITE(7,'(''EBAR,Z,A,Q2_ELAS,W0,W1='',1P10E9.2)')
C$$     >  EBAR,Z,A,Q2_ELAS,w0,w1
C$$       WRITE(7,'(/''        NU          SIGMA(nb)'')')
      ENDIF
      
      RLR=0.
      RTR=0.
      SRL=0.
C      QV=0.
      RATIO=1.
C
      W=W0
      I00=W0/WDEL
      I0=W1/WDEL
      DO 17 I=I00,I0
       WBAR=W-EBAR
       IF(WBAR.LE.0.)GO TO 15
       IF(IOPT.EQ.1)Q2=4.*EINC*(EINC-W)*SINSQ
c       WSQ = W**2
       QV2 = Q2+WSQ
       QV=SQRT(QV2)
       Q2BAR=QV2-WBAR**2
       E1=SQRT(QV2*MP2/Q2BAR+QV2/4.)-WBAR/2.
       IF(EF.LT.E1)GO TO 18     ! do not calculate 
       RATIO=Q2/(Q2+WSQ)
       F2=1./(1.+Q2/855.**2)**4
       XETA=Q2/(4.*MP2)
       ZETA=1.+XETA
!       T1=F2*Q2/2.*(((1.+2.79*XETA)/ZETA+(1.79/ZETA))**2+N/Z*3.65)
!         T1 = 2Mp**2 * DIPOLE *Tau*(MuP2 + N/Z * MuN2) = Tau*(Gmp**2 + Gmn**2 )
!        = 2Mp*Tau*(F1+(Mu-1)F2)**2 where
!          F1= DIPOLE*(1+TAU*Mu)/(1+Tau)    F2= DIPOLE/(1+Tau)
!       T2=2.*MP22*(((1.+2.79*XETA)/ZETA)**2+XETA*((1.79/ZETA)**2+
!     >  N/Z*1.91**2))*F2  !$$$*** I think the neutron term should be divided by ZETA
!         T2=2Mp**2 *DIPOLE* (Gep +Tau*Gmp)/(1+Tau)  + neutron
!           =2Mp**2 * (F1**2 +Tau*(MuP-1)F2**2)
     
! Below is Steve's Redoing
       T1 =F2*Q2/2.*(MUP2 +N/Z*MUN2)
       T2 =2.*MP2*F2*
     >   ( (1.+MUP2*XETA)/ZETA +N/Z* (0.+MUN2*XETA)/ZETA)
       E2=EF-WBAR
       E=E1
       IF(E2.GT.E1)E=E2

       RLR=(.75*Z/(KF3*QV))*(T2*((EF**3-E**3)/3.+W*(EF**2-E**2)/2.
     >  +WSQ*(EF-E)/4.)/MP2-QV2*T1*(EF-E)/Q2)

       RTR=(.75*Z/(KF3*QV))*(2.*T1*(EF-E)+T2*(Q2BAR*(EF**3-E**3)/
     > (3.*QV2)+Q2BAR*WBAR*(EF**2-E**2)/(2.*QV2)-
     > (Q2BAR**2/(4.*QV2)
     > +MP2)*(EF-E))/MP2)
15     CONTINUE
       SL=RLR*MT4PI
       ST=RTR*MT4PI

       RS=RATIO*RATIO*SL+(0.5*RATIO+TANSQ)*ST
       RSSIG=WMOTT*RS
       XSECT=FACTOR*WMOTT*RS
       IF(SL.EQ.0.AND.ST.EQ.0.)GO TO 18
C       SRL=SRL+RLR
       IF(IOPT.EQ.1)THEN
        SUM = SUM + XSECT* WDEL
C$$        WRITE(7,35)W,XSECT
       ELSE
C$$        WRITE(7,36)W,RLR,RTR
       ENDIF
35     FORMAT(1X,1P1E15.2,1X,1P1E15.2)
36     FORMAT(1X,F6.2,2E15.8)
 18    W=W+WDEL
 17   CONTINUE

      F2=1./(1.+Q2_ELAS/855.**2)**4
      XETA=Q2_ELAS/(4.*MP2)
      ZETA=1.+XETA
 
      CROSS =WMOTT* F2* 1.E33 *
     >         ( Z*( (1.+MUP2*XETA)/ZETA +2*TANSQ*MUP2 *XETA)+
     >           N*( (0.+MUN2*XETA)/ZETA +2*TANSQ*MUN2 *XETA))
C$$   WRITE(7,'(''SIGMA: SUM='',1PE9.2,'';  PEAK='',1PE9.2)')SUM,CROSS
      SUPPRESSION = SUM/CROSS
! Tsai
      QT = sqrt(Q2_ELAS**2/(2.*938.28)**2 + Q2_ELAS)
      IF(QT.GT. 2.*KF) THEN
       C=1.
      ELSE
       C = .75*QT/KF *(1.-.083333*(QT/KF)**2)
      ENDIF
C$$      WRITE(7,'('' SUPPRESSION: TSAI='',F5.3,'':   Meziani='',F5.3)')
C$$     >   C,sum/cross
      RETURN


      ENTRY  Q_E_VANORDEN_INIT(Z1,A1,KF_GEV,IOPT1)  
       Z = Z1
       A = A1
       KF= 1000.*KF_GEV
       KF3 = KF**3
       IOPT=IOPT1
       AP=ALPH/PI
       FMT=A*AMASS
       FACTOR= 4.0*PI/FMT *1.E33
       N=A-Z
       MP2 = MP**2
       MT=931.5*(Z+N)
       EF=SQRT(KF**2+MP2)
       MT4PI = MT/4./PI
      RETURN
      END
                          

                        
!---------------------------------------------------------------------
C     $MEMBER=QUADMO, DATE=75061922, USER=KCE
      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL)   
      REAL*8 FUNCT,LOWER,UPPER,EPSLON                                   
         INTEGER NLVL                                                   
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            
     1   FMRX(50), ESTRX(50), EPSX(50)                                  
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS                         
         LEVEL = 0                                                      
         NLVL = 0                                                       
         ABAREA = 0.0                                                   
         L = LOWER                                                      
         R = UPPER                                                      
         FL = FUNCT(L)                                                  
         FM = FUNCT(0.5*(L+R))                                          
         FR = FUNCT(R)                                                  
         EST = 0.0                                                      
         EPS = EPSLON                                                   
  100 LEVEL = LEVEL+1                                                   
      M = 0.5*(L+R)                                                     
      COEF = R-L                                                        
      IF(COEF.NE.0) GO TO 150                                           
         ROMBRG = EST                                                   
         GO TO 300                                                      
  150 FML = FUNCT(0.5*(L+M))                                            
      FMR = FUNCT(0.5*(M+R))                                            
      ESTL = (FL+4.0*FML+FM)*COEF                                       
      ESTR = (FM+4.0*FMR+FR)*COEF                                       
      ESTINT = ESTL+ESTR                                                
      AREA=DABS(ESTL)+DABS(ESTR)                                        
      ABAREA=AREA+ABAREA-DABS(EST)                                      
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     
         NLVL = NLVL+1                                                  
         ROMBRG = ESTINT                                                
         GO TO 300                                                      
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                         
     1         (LEVEL.LT.MINLVL))  GO TO 400                            
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             
  300    LEVEL = LEVEL-1                                                
         I = RETRN(LEVEL)                                               
         VALINT(LEVEL, I) = ROMBRG                                      
         GO TO (500, 600), I                                            
  400    RETRN(LEVEL) = 1                                               
         MX(LEVEL) = M                                                  
         RX(LEVEL) = R                                                  
         FMX(LEVEL) = FM                                                
         FMRX(LEVEL) = FMR                                              
         FRX(LEVEL) = FR                                                
         ESTRX(LEVEL) = ESTR                                            
         EPSX(LEVEL) = EPS                                              
         EPS = EPS/1.4                                                  
         R = M                                                          
         FR = FM                                                        
         FM = FML                                                       
         EST = ESTL                                                     
         GO TO 100                                                      
  500    RETRN(LEVEL) = 2                                               
         L = MX(LEVEL)                                                  
         R = RX(LEVEL)                                                  
         FL = FMX(LEVEL)                                                
         FM = FMRX(LEVEL)                                               
         FR = FRX(LEVEL)                                                
         EST = ESTRX(LEVEL)                                             
         EPS = EPSX(LEVEL)                                              
         GO TO 100                                                      
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          
      IF(LEVEL.GT.1) GO TO 300                                          
      QUADMO = ROMBRG /12.0D0                                           
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION ELASTCL_139(E,SNSQ,ZP,ZN,A,IFLG)
      IMPLICIT NONE
      REAL*8 E,SNSQ,ZP,ZN,A
      INTEGER IFLG            
      REAL*8 TM,RECOIL,EP,SIGM
      REAL QSQ_4,W1,W2,FF,E_4
      REAL*8 ALPHA/.00729735D0/                                           
                                                                       
C     GIVEN INCIDENT ELECTRON ENERGY E IN GEV                           
C           AND SNSQ OF SCATTERING ANGLE                                
C           ZP IS ATOMIC NUMBER (NO. OF PROTONS)                        
C           ZN IS AVERAGE NUMBER OF NEUTRONS                            
C           A IS ATOMIC WEIGHT (CARBON-12 SCALE) H=1.00797              
C           IFLG = 0 FOR ELASTIC, 1 FOR QUASI-ELASTIC                   
C     RETURNS THE CROSS SECTION FOR ELASTIC OR QUASI ELASTIC SCATTERING 
C     IN PB/SR (1.E-36 CM**2/SR)                                        

      TM=.93828D0                                                      
      IF (IFLG.EQ.0) TM=TM*A/1.00797D0                                  
      RECOIL=1.D0+2.D0*E/TM*SNSQ                                        
      EP=E/RECOIL                                                       
      QSQ_4=4.D0*E*EP*SNSQ                                                
      E_4 = E
      SIGM=(19732.D0*ALPHA/(2.D0*E*SNSQ))**2*(1.D0-SNSQ)/RECOIL         
C$$      CALL ASTRUL_139(QSQ,ZP,ZN,A,W1,W2,IFLG)     
      IF(IFLG.EQ.0)CALL NUC_FORM_FACTOR(QSQ_4,W1,W2,FF)
      IF(IFLG.EQ.1) CALL QE_PEAK(QSQ_4,E_4,W1,W2)   

      ELASTCL_139=SIGM*(W2+2.D0*W1*SNSQ/(1.D0-SNSQ)) 
c      write(6,'(''elastcl139'',i2,2f6.3,e11.3,2f8.3)') iflg,e,ep,
c     >  sigm,w2,w1

      RETURN                                                            
      END                                                               



      SUBROUTINE R1990F(X,QQ2,R,ERR)
C
C    Ref: L.W.Whitlow SLAC Report 357 (1990)         and
C         L.W.Whitlow et al.: PL B250 (1990) 193
C
C    for Q2 < 0.35 we extrapolate R as a constant with a rough error of 100 %
C
      REAL A(3)  / .06723, .46714, 1.89794 /,
     >     B(3)  / .06347, .57468, -.35342 /,
     >     C(3)  / .05992, .50885, 2.10807 /
C
      DATA QMAX /64./
      Q2=QQ2
      IF(Q2.LT.0.35) Q2=0.35
      FAC   = 1+12.*(Q2/(1.+Q2))*(.125**2/(X**2+.125**2))
      RLOG  = FAC/LOG(Q2/.04)
      Q2THR = 5.*(1.-X)**5
      RA   = A(1)*RLOG + A(2)/SQRT(SQRT(Q2**4+A(3)**4))
      RB   = B(1)*RLOG + B(2)/Q2 + B(3)/(Q2**2+.3**2)
      RC   = C(1)*RLOG + C(2)/SQRT((Q2-Q2THR)**2+C(3)**2)
      R     = (RA+RB+RC)/3.
      IF (Q2.GE.0.35) THEN
        Q = MIN(Q2,QMAX)
        S = .006+.03*X**2
        AA= MAX(.05,8.33*X-.66)
        XLOW  = .020+ABS(S*LOG(Q/AA))
        XHIGH = .1*MAX(.1,X)**20/(.86**20+MAX(.1,X)**20)
        D1SQUARE=XLOW**2+XHIGH**2
        D2SQUARE= ((RA-R)**2+(RB-R)**2+(RC-R)**2)/2.
        D3    = .023*(1.+.5*R)
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3
        ERR    = SQRT(D1SQUARE+D2SQUARE+D3**2)
      ELSE
        ERR = R
      ENDIF
      RETURN
      END

!---------------------------------------------------------------------

      REAL FUNCTION FITEMC(X,A,GOODFIT)                                         
!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 11/05 PYB PEB modified to use value at x=0.7 if x>0.7, because
!    beyond that assume Fermi motion is giving the rise, and we
!    already are taking that into account with the y-smearing of
!    the inelastic
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*4 ALPHA, C,LN_C,X,A ,X_U
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           

      REAL*8  ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965D-03,                                   
     >  2.18888887D+00,    3.792D-01,                                   
     > -2.46673765D+01,    6.302D+00,                                   
     >  1.45290967D+02,    4.763D+01,                                   
     > -4.97236711D+02,    1.920D+02,                                   
     >  1.01312929D+03,    4.401D+02,                                   
     > -1.20839250D+03,    5.753D+02,                                   
     >  7.75766802D+02,    3.991D+02,                                   
     > -2.05872410D+02,    1.140D+02 /                                  
                                                     
                              
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        
                                                                        
      fitemc = 1.
      if(A .lt. 2.5) return

c     IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
      IF( (X.GT.0.70).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
c       IF(X.GT. 0.88) X_U =.88
       IF(X.GT. 0.70) X_U =.70
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END                                                               

!---------------------------------------------------------------------


!---------------------------------------------------------------------

      SUBROUTINE R1998(X,Q2,R,DR,GOODFIT)                               
                                                                        
!----------------------------------------------------------------       
! X      : Bjorken x                                                    
! Q2     : Q squared in (GeV/c)**2                                      
! R      :                                                              
! DR     : Absolute error on R                                          
! GOODFIT:  = .TRUE. if the X,Q2 are within the range of the fit.       
!-------------------------------------------------------------------    
! Model for R, based on a fit to world R measurements. Fit performed by 
! program RFIT8 in pseudo-gaussian variable: log(1+.5R).  For details   
! see Reference.                                                        
!                                                                       
! Three models are used, each model has three free parameters.  The     
! functional forms of the models are phenomenological and somewhat      
! contrived.  Each model fits the data very well, and the average of    
! the fits is returned.  The standard deviation of the fit values is    
! used to estimate the systematic uncertainty due to model dependence.  
!                                                                       
! Statistical uncertainties due to fluctuations in measured values have 
! have been studied extensively.  A parametrization of the statistical  
! uncertainty of R1990 is presented in FUNCTION DR1990.                 
!                                                                       
!                                                                       
! Each model fits very well.  As each model has its own strong points   
! and drawbacks, R1998 returns the average of the models.  The          
! chisquare for each fit (237 points with 6 parameters) are:  
!          ALL DATA  #PTS=237         |     X<  0.07 #PTS= 28
! FIT  #PARAM CHISQ  CHISQ/DF PROB(%) |  CHISQ  CHISQ/DF   PROB(%)
! R1998         217.4   0.94    73.1        28.7   1.06    37.5
! R1998a   6    219.4   0.95    69.8        28.8   1.07    37.1
! R1998b   6    217.7   0.94    72.6        27.8   1.03    42.2
! R1998c   6    221.9   0.96    65.5        30.9   1.15    27.4                

!                                                                       
! This subroutine returns reasonable values for R for all x and for all 
! Q2 greater than or equal to .3 GeV.                                   
!                                                                       
! The uncertainty in R originates in three sources:                     
!                                                                       
!     D1 = uncertainty in R due to statistical fluctuations of the data 
!          as reflected in the error of the fit. 
!          It is parameterized in FUNCTION DR1990, for details see     
!          Reference.                                                   
!                                                                       
!     D2 = uncertainty in R due to possible model dependence, approxi-  
!          mated by the variance between the models.                    
!                                                                       
!     D3 = uncertainty in R due to possible epsilon dependent errors    
!          in the radiative corrections, taken to be +/- .025.  See     
!          theses (mine or Dasu's) for details.  This is copied from R1990                       
!                                                                       
! and the total error is returned by the program:                       
!                                                                       
!     DR = is the total uncertainty in R, DR = sqrt(D1+2D) 7
!          DR is my best estimate of how well we have measured R.  At   
!          high Q2, where R is small, DR is typically larger than R.  If
!          you have faith in QCD, then, since R1990 = Rqcd at high Q2,  
!          you might wish to assume DR = 0 at very high Q2.             
!                                                                       
! NOTE:    In many applications, for example the extraction of F2 from  
!          measured cross section, you do not want the full error in R  
!          given by DR.  Rather, you will want to use only the D1 and D2
!          contributions, and the D3 contribution from radiative        
!          corrections propogates complexely into F2.  For more informa-
!          tion, see the documentation to dFRC in HELP.DOCUMENT, or     
!          for explicite detail, see Reference.                         
!                                                                       
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!***** modified 10/28/01 pyb to give R at  max(0.5,Q**2)

      IMPLICIT NONE                                                     
      REAL QP2,FAC,RLOG,Q2THR,R_A,R_B,R_C,R, D1,D2,D3,DR,DR1998,X,Q2
      REAL Q2_SAVE

!!** Note, in S. Rock version, this is 0.2
      REAL Q2_LIMIT/.5/
      REAL A(6) /4.8520E-02,  5.4704E-01,  2.0621E+00,
     >          -3.8036E-01,  5.0896E-01, -2.8548E-02/
      REAL B(6) /4.8051E-02,  6.1130E-01, -3.5081E-01, 
     >          -4.6076E-01,  7.1697E-01, -3.1726E-02/
      REAL C(6) /5.7654E-02,  4.6441E-01,  1.8288E+00,
     >           1.2371E+01, -4.3104E+01,  4.1741E+01/


      LOGICAL GOODFIT                                                   
!    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      Q2_SAVE = Q2
      IF(Q2.LT.Q2_LIMIT) THEN   ! added 10/28/01 to match Blok et al
       Q2 = Q2_LIMIT
      ENDIF
                                 
      FAC = 1.+12.*Q2/(Q2+1.)*(.125**2 /(.125**2+X**2))
      RLOG  = FAC/LOG(Q2/.04)!   <--- we use natural logarithms only!      


!Model A
      QP2 = (A(2)/(Q2**4+A(3)**4)**.25) * (1.+A(4)*X +A(5)*X**2) !1.07   .030
      R_A   = A(1)*RLOG +QP2*X**A(6)                                           
!Model B
      QP2 = (1.+B(4)*X+B(5)*X**2)*(B(2)/Q2 +B(3)/(Q2**2+.09)) !1.06   .042
      R_B   =  B(1)*RLOG+QP2*X**B(6)   
!Model C
      Q2thr =C(4)*(X) +c(5)*x**2 +c(6)*X**3  
      QP2 =  C(2)/SQRT((Q2-Q2thr)**2+C(3)**2) 
      R_C   =  C(1)*RLOG+QP2    

      R     = (R_A+R_B+R_C)/3.                                          

      D1    = DR1998(X,Q2)                                              

      D2    = SQRT(((R_A-R)**2+(R_B-R)**2+(R_C-R)**2)/2.)               
      D3    = .023*(1.+.5*R)                                            
              IF (Q2.LT.1.OR.X.LT..1) D3 = 1.5*D3                       


      DR    = SQRT(D1**2+D2**2+D3**2)                                   
                                                                        
      GOODFIT = .TRUE.                                                  
      IF ((X.LT.0.02).OR.(Q2.LT.0.3)) GOODFIT = .FALSE.                                   

! Restore Q2
      Q2 = Q2_SAVE

!*** commented out 10/28/01
!      IF(Q2.LT.Q2_LIMIT) THEN   ! added 11/15 to avoid low Q2 problems caused by RLOG
!       R = Q2/Q2_LIMIT * R
!      ENDIF
      

      RETURN                                                            
      END                                                               
                                                                        
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                                        
      REAL FUNCTION DR1998(X,Q2)                                             
                                                                        
! Parameterizes the uncertainty in R1990 due to the statistical         
! fluctuations in the data.  Values reflect an average of the R-values  
! about a neighborhood of the specific (x,Q2) value.  That neighborhood 
! is of size [+/-.05] in x, and [+/-33%] in Q2.  For details, see       
! Reference.                                                            
!                                                                       
! This subroutine is accurate over all (x,Q2), not only the SLAC deep   
! inelastic range.  Where there is no data, for example in the resonance
! region, it returns a realistic uncertainty, extrapolated from the deep
! inelastic region (suitably enlarged).  We similarly estimate the      
! uncertainty at very large Q2 by extrapolating from the highest Q2     
! measurments.  For extremely large Q2, R is expected to fall to zero,  
! so the uncertainty in R should not continue to grow.  For this reason 
! DR1990 uses the value at 64 GeV for all larger Q2.                    
!                                                                       
! XHIGH accounts for the rapidly diminishing statistical accuracy for   
! x>.8, and does not contribute for smaller x.                          
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      REAL Q2,X                   

      DR1998 = .0078 -.013*X +(.070 -.39*X+.70*X**2)/(1.7+Q2)
      RETURN                                                            
      END                                                               

                                    
!------------------------------------------------------------------------
                              
!------------------------------------------------------------------------
C      DOUBLE PRECISION FUNCTION FUNCAL_139(ES1)                             
C                                                                                                                                                
CC---------------------------------------------------------------------  
CC 5/6/87: Allows OMS to be < DEELIM * XI instead of 10*XI as in         
CC          previous version called FUNCA                                
CC         Called by RDIATL                                              
CC 5/29/87: Calculates equivilent radiator at QSQ(Q2S) of interaction    
CC           instead of previous nominal QSQ (calculated elsewhere)      
CC  Steve Rock                                                           
CC 8/20/87: Use Bardin Vacuum term, use single precission logs  -Steve   
CC--------------------------------------------------------------------                                                                           FUN00130
C      IMPLICIT REAL*8 (A-H,O-Z)                                         
C      REAL BREMS_139                                         
C      COMMON/RAD/ES,EP,SNSQ,R,FAC1,BTA,BTB,ZP,ZN,A,XI,IFLG              
C      COMMON/IFIT/ IFIT                                                 
C      include '/u/ra/ser/radiative/e139code/func.common'
C      COMMON/DEELIM/DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP  
C      REAL*8  DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP  
C      DATA ALPHA/.00729735D0/,PI/3.1415926535D0/,EMSQ   /.26113D-6/     
C      DATA TWO_ALPH_PI/.00464564D0/                                     
C                                                                        
CC     CALCULATES THE INTEGRAND FOR ES1 INTEGRATION IN RADIAT            
C                                                                        
CC    CALCULATE KINEMATICS                                               
C      Q2S   =4. *    (ES1)*    (EP)*    (SNSQ)                          
C      OMS=ES-ES1                                                        
C      VS=OMS/ES                                                         
CC                                                                       
CC     QED TERMS AND MULT. PHOTON NORMALIZATION                          
C      VAC_TERM_L = 2.D0 * DELVAC(SNGL(Q2S))                             
C      DLOG_Q = DLOG(Q2S   /EMSQ   )                                     
C      VERTEX_TERM_L = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG_Q )              
C                                                                        
C      FACS=FAC1 + VAC_TERM_L + VERTEX_TERM_L                            
C                                                                        
CC     BREMSSTRAHLUNG SPECTRUM                                           
CC$$   PHIS=1.D0-VS+0.75D0*VS**2                                         
C                                                                        
CC Installed 6/4/87 to get better value of bremstrung spectrum using     
CC  Tsai's ReV Mod Physics Article eq (3.83) Steve                       
C      PHIS= BREMS_139(VS,ZP)                                                
C                                                                        
CC     PUT IT TOGETHER                                                   
C      CROSS=0.D0                                                        
C       IF(IFIT.EQ.8) THEN                                               
C        CALL XINELS(ES1,EP,SNSQ,CROSS)                                  
C       ELSE                                                             
C        IF (IFLG.EQ.0) CROSS=QUASIL(ES1,EP,SNSQ,ZP,ZN,A)                
C        IF (IFLG.EQ.1) CALL STRUCT(ES1,EP,SNSQ,CROSS,SIGMOT,110+IFIT)   
C       ENDIF                                                            
C                                                                        
CC     Equivalent raditior at scattering Q**2   5/28/87                  
C      BTEQS = (ALPHA/PI) *(DLOG_Q -1.D0)                                
C                                                                        
C      BTAS=B*DTAFTR+BTEQS                                               
C      BTBS=B*DTBEFR+BTEQS                                               
C                                                                        
CC$$   BTAS =BTA                                                         
CC$$   BTBS =BTB                                                         
C                                                                        
C      FUNCAL_139=(OMS/(R*EP))**BTAS*VS**BTBS*(BTBS/OMS*PHIS+                
C     * XI/(2.D0*DMAX1(OMS,DEELIM*XI)**2))*FACS*CROSS                    
CC                                                                       
C      RETURN                                                            
C      END                                                               
CC---------------------------------------------------------------------  
C 5/6/87: Allows OMP to be < DEELIM * XI instead of 10*XI as in         
C          previous version called FUNCB                                
C         Called by RDIATL                                              
C 5/29/83: Calculates equivilent radiator at QSQ(Q2S) of interaction    
C           instead of previous nominal QSQ (calculated elsewhere)      
C  Steve Rock                                                           
C 8/20/87: Use Bardin Vacuum term, use single precission logs  -Steve   
C--------------------------------------------------------------------   
                                                                        
                                                                        
                                                                        
C      DOUBLE PRECISION FUNCTION FUNCBL_139(EP1)                             
C                                                                        
C      IMPLICIT REAL*8 (A-H,O-Z)                                         
C      REAL BREMS_139                                                        
C      COMMON/RAD/ES,EP,SNSQ,R,FAC1,BTA,BTB,ZP,ZN,A,XI,IFLG              
C      COMMON/IFIT/ IFIT                                                 
C      include '/u/ra/ser/radiative/e139code/func.common'
C      COMMON/DEELIM/DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP
C      REAL*8 DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP                                          
C      DATA ALPHA/.00729735D0/,PI/3.1415926535D0/,EMSQ   /.26113D-6/     
C      DATA TWO_ALPH_PI/.00464564D0/                                     
C                                                                        
CC     CALCULATES THE INTEGRAND FOR EP1 INTEGRATION IN RADIAT            
C                                                                        
CC     CALCULATE KINEMATICS                                              
C      Q2P=4.D0 *(ES)*(EP1)*(SNSQ)                                       
C      OMP=EP1-EP                                                        
C      VP=OMP/EP1                                                        
C      DLOG_Q = DLOG(Q2P/EMSQ)                                           
C                                                                        
CC     QED TERMS AND MULT. PHOTON NORMALIZATION                          
C      VAC_TERM_L = 2.D0 * DELVAC(SNGL(Q2P))                             
C      VERTEX_TERM_L = TWO_ALPH_PI *(-1.D0 +.75D0 *DLOG_Q )              
C                                                                        
C      FACP=FAC1 + VAC_TERM_L + VERTEX_TERM_L                            
CC                                                                       
CC     BREMSSTRAHLUNG SPECTRUM                                           
CC$$   PHIP=1.D0-VP+0.75D0*VP**2                                         
C                                                                        
CC Installed 6/4/87 to get better value of bremstrung spectrum using     
CC  Tsai's ReV Mod Physics Article eq (3.83) Steve                       
C      PHIP= BREMS_139(VP,ZP)                                                
CC                                                                       
CC     PUT IT TOGETHER                                                   
C      CROSS=0.D0                                                        
C       IF(IFIT.EQ.8) THEN                                               
C        CALL XINELS(ES,EP1,SNSQ,CROSS)                                  
C       ELSE                                                             
C        IF (IFLG.EQ.0) CROSS=QUASIL(ES,EP1,SNSQ,ZP,ZN,A)                
C        IF (IFLG.EQ.1) CALL STRUCT(ES,EP1,SNSQ,CROSS,SIGMOT,110+IFIT)   
C       ENDIF                                                            
C                                                                        
C                                                                        
CC     Equivalenr raditior at scattering Q**2   5/28/87                  
C      BTEQP = (ALPHA/PI) *(DLOG_Q -1.D0)                                
C                                                                        
C      BTAP=B*DTAFTR+BTEQP                                               
C      BTBP=B*DTBEFR+BTEQP                                               
C                                                                        
CC$$   BTAP =BTA                                                         
CC$$   BTBP = BTB                                                        
C                                                                        
C      FUNCBL_139=VP**BTAP*(R*OMP/ES)**BTBP  *(BTAP/OMP*PHIP+                
C     * XI/(2.D0*DMAX1(OMP,DEELIM*XI)**2))*FACP*CROSS                    
CC                                                                       
C      RETURN                                                            
C      END                            

!------------------------------------------------------------------------
                                   
C $MEMBER=FUNX, DATE=75061922, USER=KCE                                 
      DOUBLE PRECISION FUNCTION FUNXL_139(COSK) 
      IMPLICIT NONE                            

      REAL*8 COSK
      COMMON/TAL1/ES,EP,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1 
      REAL*8  ES,EP,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1 
      COMMON/TAL2/ZP,ZN,AT,IFLG                                         
      REAL*8 ZP,ZN,AT
      INTEGER IFLG
      REAL*8 OM,QSQ,FAC2,VAC_TERM,VERTEX_TERM,FAC,FW1,FW2,A,AP,BP,GNU
      REAL*8 X,Y,TERM1,TERM21,TERM22,TERM23,TERM24,TERM25,TERM26,TERM2
      REAL*8 EMSQ/.26113D-6/,PI/3.1415926535D0/                           
      REAL*8 ALPHA/7.2973515D-3/
      REAL*8 TWO_ALPH_PI/.00464564D0/  
      REAL*4 QSQ_4,W1,W2,FF,ESOM_4,delvac
      common/testing/ prttst,usegrd
      logical prttst,usegrd

                                                                        
C -----------------------------------------------------------------     
C 8/28/87: Single preciession log - Steve Rock                          
C----------------------------------------------------------------- 
      OM=.5D0*(USQ-TMSQ)/(U0-UVEC*COSK)                                 
      QSQ=2.D0*EMSQ-2.D0*SP-2.D0*OM*(ES-EP-UVEC*COSK)  
      QSQ_4 = -QSQ
      ESOM_4 = ES-OM
      IF(IFLG.EQ.0)CALL NUC_FORM_FACTOR(QSQ_4,W1,W2,FF) 
      IF(IFLG.EQ.1) CALL QE_PEAK(QSQ_4,ESOM_4,W1,W2)   
C$$      CALL ASTRUL_139(-QSQ,ZP,ZN,AT,W1_8,W2_8,IFLG)

                                                                        
C*******8/13/87**  Steve Rock ***************************************   
C Use the Bardin Calculation of Vertex terms instead of only electron   
C  and Tsai's Vertex term.                                              
C     FAC2=2.D0*ALPHA/PI*(-14.D0/9.D0+13.D0/12.D0*DLOG(-QSQ/EMSQ))      
C                                                                       
C  Vacuum terms of Bardin including 3 leptons and quarks                
      VAC_TERM = 2. * DELVAC(SNGL(-QSQ))                                 
C  Vertex term of Tsai (eq. 2.11)                                       
      VERTEX_TERM =TWO_ALPH_PI *(-1.+.75 *ALOG(-SNGL(QSQ)/SNGL(EMSQ)) ) 
      FAC2 = VAC_TERM + VERTEX_TERM                                     
C*******************************************************************    
                                                                        
                                                                        
      FAC=FAC1+FAC2                                                     
      FW1=FAC*W1                                                        
      FW2=FAC*W2                                                        
      A=OM*(EP-PVEC*COSP*COSK)                                          
      AP=OM*(ES-SVEC*COSS*COSK)                                         
                                                                        
C      Tsai on page 40 of SLAC-PUB-898 (1971) says that an uncertainty  
C      of zero / zero occurs when a = a' in his equation A.24, the      
C      integrand of which is calculated by this function.  To avoid this
C                                                                       
      FUNXL_139= 0.D0                                                       
      IF (A.EQ.AP) RETURN                                               
C                                                                       
      BP=-OM*PVEC*DSQRT((1.D0-COSP**2)*(1.D0-COSK**2))                  
      GNU=1.D0/(AP-A)                                                   
      X=DSQRT(A**2-BP**2)                                               
      Y=DSQRT(AP**2-BP**2)                                              
      TERM1=(A/X**3+AP/Y**3)*EMSQ*(2.D0*EMSQ+QSQ)+4.D0                  
     * +4.D0*GNU*(1.D0/X-1.D0/Y)*SP*(SP-2.D0*EMSQ)                      
     * +(1.D0/X-1.D0/Y)*(2.D0*SP+2.D0*EMSQ-QSQ)                         
      TERM21=2.D0*ES*(EP+OM)+.5D0*QSQ                                   
      TERM22=2.D0*EP*(ES-OM)+.5D0*QSQ                                   
      TERM23=2.D0*ES*EP-SP+OM*(ES-EP)                                   
      TERM24=2.D0*(ES*EP+ES*OM+EP**2)+.5D0*QSQ-SP-EMSQ                  
      TERM25=2.D0*(ES*EP-EP*OM+ES**2)+.5D0*QSQ-SP-EMSQ                  
      TERM26=EMSQ*(SP-OM**2)+SP*TERM23                                  
      TERM2=-A*EMSQ/X**3*TERM21-AP*EMSQ/Y**3*TERM22                     
     * -2.D0+2.D0*GNU*(1.D0/X-1.D0/Y)*TERM26+1.D0/X*TERM24              
     * -1.D0/Y*TERM25                                                   
      FUNXL_139=(FW2*TERM2+FW1*TERM1)*OM/(QSQ**2*(U0-UVEC*COSK)) 
c      If(FUNXL_139.GT.1.D15) THEN
c       BP= BP*(1.+1.D-20           )
c      ENDIF
      if(prttst) write(44,'(1x,7e12.5,10f7.3)')cosk,FUNXL_139,
     >   term1,fac,fac1,fac2
      RETURN                                                            
      END                                                               
!------------------------------------------------------------------------------------------

C $MEMBER=FUNXI, DATE=75061922, USErqsq=KCE                                
C      DOUBLE PRECISION FUNCTION FUNCOSK(COSK)                           
C      IMPLICIT REAL*8 (A-H,O-Z)                                         
C      REAL*8 MASS_MAX,MASS_MIN/1.0783D0/                                
C      COMMON/TAL1/ES,EP,SVEC,PVEC,SP,USQ,U0,UVEC,COSP,COSS,TMSQ,FAC1    
C      COMMON/TAL2/ZP,ZN,AT,IFLG                                         
C      COMMON/SOFT/ OM,B,SNSQ                                            
C      COMMON/DEELIM/DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP  
C      REAL*8 DEELIM,DEE,CRSLIM,XI4,CRSES,CRSEP                
C      COMMON/COSK/COSK_PASS                                             
C      EXTERNAL FUNMASS                                                  
C      DATA EMSQ/.26113D-6/,PI/3.1415926535D0/                           
C      DATA TM/.93828D0/                                                
C      DATA ALPHA/7.2973515D-3/                                          
C      DATA TWO_ALPH_PI/.00464564D0/                                     
C      COSK_PASS = COSK                                                  
C                                                                        
C      MASS_MAX = DSQRT(USQ -2.D0 *(U0-UVEC*COSK)*DEE )                  
C                                                                        
C      FUNCOSK = QUADMO1(FUNMASS,MASS_MIN,MASS_MAX,3.D-4,NLVL)           
C                                                                        
C      RETURN                                                            
C      END                                                               
