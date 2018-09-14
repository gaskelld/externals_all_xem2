!---------------------------------------------------------------------

      FUNCTION f2allm(x,q2)             
!-----------------------------------------------------------------
c       allm97, NMC published measured points Q2>0.75 GeV2
c       for values Q<1 use data of E665!
c       parameterization of F2 , according to
c       H.Abramowicz and A.Levy, hep-ph/9712415
c
c       3*10-6 < x  < 0.85, W2>3GeV2
c       0.   < Q2 < 5000 GeV2, dof=0.97
c
! From Peter Bosted, Nov 00 from Antje Burrel.
!------------------------------------------------------------------
      IMPLICIT NONE
      REAL F2ALLM,X,Q2
      REAL SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
      COMMON/ALLM/SP,AP,BP,SR,AR,BR,S,XP,XR,F2P,F2R
C  POMERON
      REAL S11,S12,S13,A11,A12,A13,B11,B12,B13,M12
      PARAMETER (
     , S11   =   0.28067, S12   =   0.22291, S13   =   2.1979,
     , A11   =  -0.0808 , A12   =  -0.44812, A13   =   1.1709,
     , B11   =   0.60243**2, B12   =   1.3754**2, B13   =   1.8439,
     , M12   =  49.457 )
 
C  REGGEON
      REAL S21,S22,S23,A21,A22,A23,B21,B22,B23,M22
      PARAMETER (
     , S21   =   0.80107, S22   =   0.97307, S23   =   3.4942,
     , A21   =   0.58400, A22   =   0.37888, A23   =   2.6063,
     , B21   =   0.10711**2, B22   =   1.9386**2, B23   =   0.49338,
     , M22   =   0.15052 )
C
      REAL M02,LAM2,Q02,ALFA,XMP2
      PARAMETER ( M02=0.31985, LAM2=0.065270, Q02=0.46017 +LAM2 )
      PARAMETER ( ALFA=112.2, XMP2=0.8802)
      REAL W2,W,Z
C                                                                               
      W2=q2*(1./x -1.)+xmp2
      W=sqrt(w2)
C
      IF(Q2.EQ.0.) THEN                                                        
       S=0.
       Z=1.           
C                                                                               
C   POMERON                                                                     
C                                                                               
       XP=1./(1.+(W2-XMP2)/(Q2+M12))                                   
       AP=A11                                                            
       BP=B11                                                               
       SP=S11                                                             
       F2P=SP*XP**AP                                                
C                                                                               
C   REGGEON                                                                     
C                                                                               
       XR=1./(1.+(W2-XMP2)/(Q2+M22))          
       AR=A21                                 
       BR=B21                                 
       SR=S21
       F2R=SR*XR**AR              
C                                                                               
      ELSE                                                                      
       S=LOG(LOG((Q2+Q02)/LAM2)/LOG(Q02/LAM2))                       
       Z=1.-X                                      
C                                                                               
C   POMERON                                                                     
C                                                                               
       XP=1./(1.+(W2-XMP2)/(Q2+M12))                
       AP=A11+(A11-A12)*(1./(1.+S**A13)-1.)         
       BP=B11+B12*S**B13                                                  
       SP=S11+(S11-S12)*(1./(1.+S**S13)-1.)         
       F2P=SP*XP**AP*Z**BP                            
C                                                                               
C   REGGEON                                                                     
C                                                                               
       XR=1./(1.+(W2-XMP2)/(Q2+M22))                                          
       AR=A21+A22*S**A23                                                     
       BR=B21+B22*S**B23                                                
       SR=S21+S22*S**S23                                                     
       F2R=SR*XR**AR*Z**BR                                                   
    
C                                                                               
      ENDIF

c      CIN=ALFA/(Q2+M02)*(1.+4.*XMP2*Q2/(Q2+W2-XMP2)**2)/Z              
c      SIGal=CIN*(F2P+F2R)                                             
c      f2allm=sigal/alfa*(q2**2*(1.-x))/(q2+4.*xmp2*x**2)
      F2ALLM = q2/(q2+m02)*(F2P+F2R)
 

      RETURN                                                                    
      END                                                                       
