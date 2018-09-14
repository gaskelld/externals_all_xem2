CDJG: Various radcor functions/routines

C=======================================================================
                                                                        
      REAL*4   FUNCTION FTSAI(E0,EP,TH)                                          
      IMPLICIT NONE
      REAL E0,EP,TH
                                                                  
      COMMON /CON/ PI,PM,PP,EM,AL
      REAL PI,PM,PP,EM,AL 
      REAL ALPI,Q2,VAC,VERTEX,ARG,SMALL,CORR,DELVAC,DSPEN
      PARAMETER    (ALPI=2.32281E-03)                                   
                                                                        
      Q2     = 4.*E0*EP*SIN(TH/2.)**2                                   
      VAC    = 2.*DELVAC(Q2)                                            
      VERTEX = 2.*ALPI*(-1.+0.75*ALOG(Q2/EM**2))                        
      ARG    = COS(TH*.5)**2                                            
      SMALL  = ALPI*(PI**2/6.-DSPEN(ARG))                               
      CORR   = -ALPI/2.*LOG(E0/EP)**2                                   
      FTSAI  = 1.+VAC+VERTEX+SMALL+CORR                                 
                                                                        
      RETURN                                                            
      END          

C=======================================================================
                                                                        
      REAL*4 FUNCTION BREMS(E,EPS,T,ID)                                        
      IMPLICIT NONE
      REAL E,EPS,T,tmp
      INTEGER ID
C Bremstrahlung and ionization energy loss probablity function. E is the
C incoming energy of particle passing through T radiation lengths to    
C lose EPS energy. Formula from Tsai; Rev of Mod Physics 46,(1974)815,  
C equation 3.83. Complete screening.                                    
                                                                        
      COMMON /RAD/ ETA,B,AX0,BA,AX0A,                                   
     >             ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL   
      REAL   ETA,B,AX0,BA,AX0A,                                 
     >       ATBX0,ATAX0,TB,TA,BTBI,BTAI,BTB,BTA,
     >               TBEFOR,TBAL,TAFTER,TAAL                     
      COMMON /CON/ PI,PM,PP,EM,AL   
      REAL PI,PM,PP,EM,AL 
      REAL Y,WB,WI                                    
      REAL ATX0/0/, BT/0./  !initialization to make vm compiler happy   
                                                                        
      IF (T.EQ.0.) THEN                                                 
        BREMS = 0.                                                      
        RETURN                                                          
      ENDIF                                                             
                                                                        
      IF (ID.EQ.1) THEN                                                 
        ATX0 = ATBX0                                                    
        BT   = BTBI                                                     
      ELSE IF (ID.EQ.2) THEN                                            
        ATX0 = ATAX0                                                    
        BT   = BTAI                                                     
      ENDIF                                                             
                                                                        
C BREMS is Tsai eq.(3.83) normalized to be 1 at Y=0.  The higher order  
C order terms, with 2 approximations, result in the "-f" term of eq.    
C (3.66):                                                               
      Y     = EPS/E                                                     
      if(y.le.0.) then
c        write(6,'(1x,''ERROR, eps,e,atx0,bt='',4f10.3)')
c     >    eps,e,atx0,bt
        brems=0.
        return
      endif

      tmp = (1-Y+.75*Y**2)                                            
      WB    = B*T/EPS*tmp                                             
                                                                        
C Ionization:                                                           
      WI    = ATX0/MAX(EPS,10.*ATX0)**2*(1.+EPS**2/E/(E-EPS))**2        
      BREMS = (1.+0.5772*BT-0.62*BT**2) * Y**(B*T) * (WB+WI)                
                                                                        
      RETURN                                                            
      END                                          

C=======================================================================
                                                                        
      REAL*4 FUNCTION DSPEN(X)               

      IMPLICIT NONE                              
      REAL*4 X,FSPENS
      REAL  F1 /1.644934/                                       
                                                                        
      IF (X.LT.-1.) THEN                                                
           DSPEN = -.5*ALOG(1.-X)*ALOG(X**2/(1.-X))-F1+FSPENS(1./(1.-X))
      ELSE IF (X.LT.0.) THEN                                            
           DSPEN = -.5*ALOG(1.-X)**2-FSPENS(X/(X-1.))                   
      ELSE IF (X.LE..5) THEN                                            
           DSPEN = FSPENS(X)                                            
      ELSE IF (X.LE.1.) THEN                                            
           DSPEN = F1-ALOG(X)*ALOG(1.-X+1.E-10)-FSPENS(1.-X)            
      ELSE IF (X.LE.2.) THEN                                            
           DSPEN = F1-.5*ALOG(X)*ALOG((X-1.)**2/X)+FSPENS(1.-1./X)      
      ELSE                                                              
           DSPEN = 2*F1-.5*ALOG(X)**2-FSPENS(1./X)                      
      ENDIF                                                             
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
      REAL*4 FUNCTION FSPENS(X)
      IMPLICIT NONE
      REAL*4 X,A,F,AN,TCH,B

      A   = 1.                                                          
      F   = .0                                                          
      AN  = .0                                                          
      TCH = 1.E-16                                                      
1     AN  = AN+1.                                                       
      A   = A*X                                                         
      B   = A/AN**2                                                     
      F   = F+B                                                         
      IF (B.GT.TCH) GO TO 1                                             
      FSPENS = F                                                        
      RETURN                                                            
      END                                                               
                               
C=======================================================================
                                                                        
      REAL*4 FUNCTION DELVAC(T)                                                
                                                                        
C VACUUM POLARIZATION FOR LEPTONS AND HADRONS, includes terms for       
C electrons, muons, taus, and quarks.  The expression for small mass    
C reduces to Tsai's expression (2.10) in the slac pub.                  
      IMPLICIT NONE
      REAL*4 T, AMF2(9), COL(9), RTMF2,RT,ALT,AI34,SUML,A,B,C,SUMH
      INTEGER IC
      REAL*4  EPST / 1.E-3 /                                         
      REAL*4     ALPH_PI /.0023229 /                                    
      DATA       AMF2 / .26112E-6,.011164,3.1827,.0064,.0064,2.25,      
     +                  .09,1024.,20.25 /                               
      DATA       COL  / 3*1.,6*3. /                                     
                                                                        
C For newer treatment of vacuum correction, comment out the following   
C three lines.  These lines are exactly (2.10) of Tsai's Slac Pub '74:  
                                                                        
C     EM     = .000511                                                  
C     DELVAC = ALPH_PI*(-5./9.+ALOG(T/EM**2)/3.)                        
C     RETURN                                                            
                                                                        
      SUML = 0.                                                         
      DO 121 IC=1,3                                                     
           RTMF2 = T/AMF2(IC)                                           
           IF (RTMF2.LT.EPST) THEN                                      
                AI34 = RTMF2*(2./15.-1./70.*RTMF2)                      
           ELSE                                                         
                RT   = SQRT(1.+4./RTMF2)                                
                ALT  = RT*ALOG(RTMF2*(1.+RT)**2/4.)                     
                AI34 = -5./9.+4./3./RTMF2+(1./3.-2./3./RTMF2)*ALT       
           ENDIF                                                        
           SUML = SUML+COL(IC)*AI34                                     
121   CONTINUE                                                          
                                                                        
C HADRONIC VACUUM POLARIZATION TAKEN FROM BURKHARD, TASSO NOTE 192, 1982
                                                                        
      IF (T.LT.1.) THEN                                                 
           A = -1.345E-9                                                
           B = -2.302E-3                                                
           C =  4.091                                                   
      ELSE IF (T.LE.64.) THEN                                           
           A = -1.512E-3                                                
           B = -2.822E-3                                                
           C =  1.218                                                   
      ELSE                                                              
           A = -1.1344E-3                                               
           B = -3.0680E-3                                               
           C =  9.9992E-1                                               
      ENDIF                                                             
                                                                        
      SUMH   = -(A+B*ALOG(1.+C*T))                                      
      DELVAC = SUML*ALPH_PI+SUMH                                        
                                                                        
      RETURN                                                            
      END                                                               
                                                
C=======================================================================
                                                                        
      subroutine radlength(X0)                                          
                                                                        
C calculates the radiation length based on Tsai74, table III.6, except  
C for iZ=1 values taken from Physics Letters 170B--Review of Particle   
C Properties. amuM is average atomic mass in C12-based atomic mass      
C units:                                                                
                                                                        
      common /targt/ iZ,iA,avgN,avgA,avgM,amuM                          
                                                                        
      if (iZ.le.4) then                                                 
           if (iA.eq.1) X0 =  61.28                                     
           if (iA.eq.2) X0 = 122.6                                      
           if (iA.eq.3.and.iZ.eq.1) X0 = 183.7                          
           if (iZ.eq.2) X0 =  94.322*amuM/4.0026                        
           if (iZ.eq.3) X0 =  82.756*amuM/6.9390                        
           if (iZ.eq.4) X0 =  65.190*amuM/9.0122                        
           return                                                       
      endif                                                             
                                                                        
      zz    = (iZ*.00729735)**2                                         
      f     = 1.202*zz - 1.0369*zz**2 + 1.008*zz**3/(1.+zz)             
      denom = (log(184.15/iZ**.3333)-f)*iZ**2 + log(1194/iZ**.6667)*iZ  
                                                                        
      X0 = 716.405*amuM/denom                                           
                                                                        
      return                                                            
      end   
                           
!-----------------------------------------------------------------------
      REAL FUNCTION BREMS_139 (Y,Z)                                              
                                                                        
C=======================================================================
C                                                                       
C  Returns value of Bremstrung probability multiplied by k (photon energ
C     normalized to 1 at Y=0                                            
C                                                                       
C  Formula from Tsai; Rev of Mod Physics 46,(1974)815,  equation 3.83   
C      Complete screening                                               
C                                                                       
C  Y = k/E : fraction of beam energy carried by photon                  
C  Z =     : atomic number of target                                    
C                                                                       
C      6/4/87 Steve Rock                                                
C=======================================================================
      IMPLICIT NONE
                                                                        
      REAL*8 Z,Y                                                        
      REAL LRAD,LRADP,FC,ZZ                                          
      REAL ALPHA/7.2973515E-3/                                          
      LOGICAL FIRST /.TRUE./                                                                  
                                                                        
C TSAI page 486 Table B.2  and eq. 3.67, 3.68                           
      IF(FIRST) THEN                                                                  
       IF(Z.EQ.1)THEN                                                    
        LRAD = 5.31                                                      
        LRADP= 6.144                                                     
       ELSEIF(Z.EQ.2) THEN                                               
        LRAD = 4.79                                                      
        LRADP=5.62                                                       
       ELSEIF(Z.EQ.3) THEN                                               
        LRAD = 4.74                                                      
        LRADP= 5.805                                                     
       ELSEIF(Z.EQ.4) THEN                                               
        LRAD = 4.71                                                      
        LRADP= 5.924                                                     
       ELSE                                                              
        LRAD = ALOG(184.15) - DLOG(Z)/3.                                 
        LRADP= ALOG(1194.) -2.*DLOG(Z)/3.                                
       ENDIF                                                             
                                                                        
C Get coulomb correction from Tsai eq. 3.3                              
       ZZ=  (Z*ALPHA)**2
       FC = 1.202 * ZZ - 1.0369 * ZZ**2 +                  
     >  1.008 *ZZ**3/(1.+ ZZ)
       FIRST=.FALSE.
      ENDIF                                
                                                                        
                                                                        
C  Tsai eq.(3.83) normalized to be 1 at Y=0                             
      BREMS_139 = (1 -Y +.75 *Y**2) +                                       
     >         (1.-Y)/12. *(Z+1.)/(Z *(LRAD -FC) + LRADP)               
                                                                        
      RETURN                                                            
      END                                       
