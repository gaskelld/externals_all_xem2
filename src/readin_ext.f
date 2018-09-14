c      SUBROUTINE READIN_EXT(LEN_FILENAME,FILENAME)     
      SUBROUTINE READIN_EXT()     
! used to be Called from cexternal.c
      IMPLICIT NONE
      COMMON /TARGT/  iZ,iA,avgN,avgA,avgM,amuM 
      INTEGER IZ,IA
      REAL avgN,avgA,avgM,amuM                                   
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      COMMON /TTYPE/  TARGET                                            
      CHARACTER*7     TARGET
      COMMON /TRL/    TTARG,TWALL,TBEAM,TSPEC,ITARG,NSEG,JTARG
      REAL TTARG,TWALL,TBEAM,TSPEC
      INTEGER ITARG,NSEG,JTARG
      INTEGER I,J
      CHARACTER*72    COMMENT(6),RAD_STRING/'  '/                           
      CHARACTER*80    EXTERNAL_RUNPLAN,EXTERNAL_TARGET,EXTERNAL_OUT
      character*80    filename
      character        dbase_file*60
      character        base*40
      logical doeg1b
      common/experiment/ doeg1b
      logical smrelasp,usearen
      real smrdEp
      common/smrproton/ smrelasp,smrdEp,usearen


c      OPEN(UNIT=20,FILE=FILENAME(1:LEN_FILENAME))
c      OPEN(UNIT=20,FILE='input/externals/kpp_shms_488.inp')
      write(6,*) 'Enter the input file name (in INP directory)'
      read(5,'(a)') dbase_file
      j=index(dbase_file,'/')
      if (j.eq.0) then          !no path given
         write(filename,'(a)') 'input/INP/'//dbase_file
      else
         write(filename,'(a)') dbase_file
      endif
      i=index(filename,'.')
      j=index(filename,'/')
      if(i.eq.0) then		!add .inp if not included in filename
         i=index(filename,' ')
         if(i+2.le.len(filename)) write(filename(i:),'(''.inp'')')
      endif
      write(6,'(a10,a69)')'filename=',filename
      if (i.gt.1) base=filename(j+5:i-1)
      write(6,*) 'cheesy poofs',filename
      OPEN(UNIT=20,FILE=filename)
      READ(20,'(A)') COMMENT(1)
      WRITE(6,'(A)') COMMENT(1)
      READ(20,'(A)') EXTERNAL_RUNPLAN
      WRITE(6,'(A)') EXTERNAL_RUNPLAN
      READ(20,'(A)') EXTERNAL_TARGET
      WRITE(6,'(A)') EXTERNAL_TARGET
      READ(20,'(A)') EXTERNAL_OUT
      WRITE(6,'(A)') EXTERNAL_OUT
      close(unit=20)

      OPEN(UNIT=10,FILE=EXTERNAL_OUT)
      OPEN(UNIT=7,FILE=EXTERNAL_RUNPLAN)  
      OPEN(UNIT=5,FILE=EXTERNAL_TARGET)     

      do i=1,6
        READ(7,'(a)') COMMENT(i)             
        write(6,'(a)') comment(i)
      enddo
      doeg1b=.false.
      if(comment(1)(1:4).eq.'EG1b' .or.
     >   comment(1)(2:5).eq.'EG1b') then
        doeg1b=.true.
        write(6,'(''using eg1b format'')')
      endif
      READ(5,100) (COMMENT(i),i=4,6),                                   
     +            iZ,iA,                                                
     +            avgA,avgM,                                            
     +            target,                                               
     +            ttarg,twall,tbeam,tspec,                              
     +            NSEG,                                                 
     +    IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL

    
      IF(TARGET.EQ.'E140XH2') ITARG=1                                   
      IF(TARGET.EQ.'E140XH1') ITARG=2                                   
      IF(TARGET.EQ.'E140XD1') ITARG=3                                   
      IF(TARGET.EQ.'E140XD2') ITARG=4                                   
      IF(TARGET.EQ.'E140XAL') ITARG=5                                   
      IF(TARGET.EQ.'E140XAS') ITARG=6                                   
      IF(TARGET.EQ.'E140XBE') ITARG=7                                   
100   FORMAT      (3(A72,/),/,                                          
     +             2(15X,I5,/),/,                                       
     +             2(15X,F14.8,/),/,                                    
     +              (20X,A7,/),/,                                       
     +             4(15X,F14.8,/),                                      
     +             15X,I5,/,/,/,         ! NSEG, skip blank and ing   
     +             15X,I5,/,/,   ! IG,  elastic nucleon form factor model
     +             15X,I5,/,     ! IDUT                                      
     +             15X,I5/    ! INEL_MODEL
     +             15X,I5/     !Pauli suppression     
     +             15X,I5/     ! Nuclear Tail Method: NUC_METHOD
     +             15X,I5)     ! Nuclear Form Factor: NUC_MODEL

      call Q_E_VANORDEN_INIT(REAL(iZ),REAL(iA),.25,1)!Fermi Momentum=.25
      if (avgA.eq.0..or.avgM.eq.0) call weiz                            
      if(IA.eq.1) avgA = 1.
! pyb change avgM if smearing proton elastic
      if(smrelasp) avgM = max(avgM, 2.0)
      avgN = avgA-iZ                                                    
      amuM = avgM/.931501                                               
      IF(INDEX(TARGET,'E143').GT.0)CALL RADLEN43_INIT(RAD_STRING)
      WRITE(10,'(''EXTERNAL/TSAI:  VERSION OF 8/8/96''/
     > /A72/4(A72/)/,A72/A72/)')   
     +COMMENT(1),COMMENT(2),COMMENT(3),COMMENT(4),COMMENT(5),COMMENT(6),
     > RAD_STRING
      WRITE(10,200) iZ,iA,avgA,amuM,avgM,target,ttarg,twall,tbeam,tspec,
     + Nseg,IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
200   FORMAT('TARGET (iZ,iA) = (',I3,',',I3,');   avgA = ',F7.3,        
     +                 '; amuM=',F7.3,'; avgM=',F7.3,/                  
     +                 ' targ=',A7,'; ttarg(rl)=',F7.5,'; twall=',F7.5, 
     +                 '; tbeam=',F7.5,'; tspec=',F7.5/                 
     +       '         Nseg  = ',I3/' MODELS: P elastc(ikk)=',i2,
     >     '; deut elastic=',i2, '; inelastc=',I3, 
     >     '; Pauli Supresion=',i2 /
     >     ' NUCLEAR METHOD=',I2,';   NUCLEAR_MODEL=',I2/)
                                                               
      i = index(base,' ')      	
      filename = 'output/OUT/'//base(1:i-1)//'.out'                             
      open(unit=17,file=filename)
      write(17,*) '***Ebeam    Eprime  Theta     x        Q2'//
     > '      Sig_Born     Sig_Born_In  Sig_Born_QE  Sig_Rad'//
     > '      Sig_Rad_EL   Sig_Rad_QE   Sig_Rad_DIS  C_cor'
c      write(17,*) '***Ebeam    Eprime  Theta, x, Q2  Sig_Born'//
c     > '  Sig_Born_Inel, Sig_Born_QE, Sig_Rad Sig_Rad_EL'//
c     > '  Sig_Rad_QE Sig_Rad_DIS C_cor'
      i = index(base,' ')      	
      filename = 'output/OUT/'//base(1:i-1)//'_int.out'                             
      open(unit=18,file=filename)

         
      RETURN                                                            
      END                                     
