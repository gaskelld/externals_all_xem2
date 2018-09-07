C====================================================================== 
      SUBROUTINE FYXSEC8(E0,E1,THT,XSEC)                                
C                                                                       
      DIMENSION x(1),a(4)                                               
      common/sigc /qx,qy,qz,q2,px,py,pz,p2,tpq                          
     >,pxf,pyf,pzf                                                      
      COMMON     /TARGT/ iZ,iA,avgN,avgA,avgM,amuM                      
C initialize some values                                                
      nz = iz                                                           
      n = iA - iZ                                                       
      e0m=e0*1000.                                                      
      e1m=e1*1000.                                                      
!  iron coeff.                                                          
      nf=11                                                             
      nt=4                                                              
                                                                        
       if(iz.eq.26) then                                                
        a(1)=2.8316                                                     
        a(2)=44.624                                                     
        a(3)=0.3785                                                     
        a(4)=12.003                                                     
!   carbon coeff                                                        
       else                                                             
        a(1)=2.8757                                                     
        a(2)=41.922                                                     
        a(3)=0.3380                                                     
        a(4)=13.824                                                     
       endif                                                            
      XSEC = 0.0                                                        
      ep=e1m                                                            
      pz = pper                                                         
C     separation energy                                                 
      if(iz.eq.4) eps = 16.                                             
      if(iz.eq.26) eps=10.                                              
      hbc=0.1973289e+03                                                 
      fscnst=1./137.04                                                  
      rads=3.14159/180.                                                 
      rmp=.93828e+03                                                     
      tpm=2.*rmp                                                        
      tpm1=1./tpm                                                       
      thit=tht                                                          
      z=nz                                                              
      na=n+nz                                                           
      na1=na-1                                                          
      tma=2.*na*rmp                                                     
      tma1=2.*rmp*(na-1)                                                
      th2r=tht*rads/2.                                                  
      rma=na*rmp                                                        
      tn2=(tan(th2r))**2                                                
      thr=tht*rads                                                      
      tn2=(tan(th2r))**2.                                               
      Cost=COS(THT*RADS)                                                
      thtr=tht*rads                                                     
      sint=sin(thtr)                                                    
      cs2=(cos(th2r))**2.                                               
      s2=(sin(th2r))**2.                                                
      om=e0m-e1m                                                        
      ef=e0m-om                                                         
      q42=4.*e0m*ef*s2                                                  
      qx=e0m-ef*cost                                                    
      qy=-ef*sint                                                       
      qz=0.0                                                            
      qv2=q42+om**2                                                     
      qv=sqrt(qv2)                                                      
      q2=qv2                                                            
      qmp=q42/hbc/hbc                                                   
      tau=q42/4./rmp**2.                                                
      thpr=thp*rads                                                     
      csp=cos(thpr)                                                     
      tangle=qx/(-qy)                                                   
      thx=atan(tangle)                                                  
      e0gev = e0m/1000.                                                 
      epgev = e1m/1000.                                                 
      epsgev = eps/1000.                                                
      call yvalue(e0gev,epgev,tht,na,                                   
     >     epsgev,ygev)                                                 
      y = ygev                                                          
      p = ygev*1000.                                                    
      yp2=p                                                             
      px = p*sin(thx)                                                   
      py = -p*cos(thx)                                                  
      pxf=px+qx                                                         
      pyf=py+qy                                                         
      pzf=pz                                                            
      tpq=2.*e0m*px-(2.*px*cost+2.*py*sint)*ef                          
      p2=px*px+py*py+pz*pz                                              
      p=sqrt(p2)                                                        
      call sigcc(e0m,e1m,tht,sigp,sign)                                 
      sigt=z*sigp+n*sign                                                
      if (y .le. 0.0)x(1) = y                                           
      if (y .ge. 0.0)x(1) =-abs(y) ! sym                                
C     sum of two gaussians one located at zero                          
      functn = a(1)*exp(-a(2)*x(1)**2) +                                
     >    a(3)*exp(-a(4)*x(1)**2)                                       
C    ciofi's version                                                    
      top = qv                                                          
      bottom = sqrt(rmp**2 + qv2 + yp2**2 +                             
     > 2.*qv*yp2)                                                       
      domkdcos = top/bottom                                             
      dkcosdom = 1./domkdcos                                            
      dydom = dkcosdom                                                  
      xsec=fy*dydom*sigt*1000000.                                       
C      the om was in MeV                                                
22    return                                                            
      END                                                               
                                                                        
       subroutine yvalue(e0gev,epgev,theta,iatomic,eps,y)               
cccc       Implicit real*8 (a-h,o-z)                                        
!                                                                       
!subroutine to calulate the value of y given e,ep,theta,eps             
!                                                                       
       data rmp/0.93828/                                                 
       data rads/.01745329/                                             
       data thp/0.0/                                                    
       y = 0.0                                                          
       tpm=2.*rmp                                                       
       tpm1=1./tpm                                                      
!                                                                       
       na=iatomic                                                       
       na1=na-1                                                         
       tma=2.*na*rmp                                                    
       tma1=2.*rmp*(na-1)                                               
       th2r=theta*rads/2.                                               
       rma=na*rmp                                                       
       tn2=(tan(th2r))**2                                               
       thr=theta*rads                                                   
       tn2=(tan(th2r))**2.                                              
       cost=cos(theta*rads)                                             
       thtr=theta*rads                                                  
       sint=sin(thtr)                                                   
       cs2=(cos(th2r))**2.                                              
       s2=(sin(th2r))**2.                                               
                                                                        
       om=e0gev-epgev                                                   
       q42=4.*e0gev*epgev*s2                                            
       qx=e0gev-epgev*cost                                              
       qy=-epgev*sint                                                   
       qz=0.0                                                           
       qv2=q42+om**2                                                    
       qv=sqrt(qv2)                                                     
       q2=qv2                                                           
       thpr=thp*rads                                                    
       csp=cos(thpr)                                                    
       tangle=qx/(-qy)                                                  
       thx=atan(tangle)                                                 
                                                                        
       w=om - eps +na*rmp                                               
       wp=w**2+(na1*rmp)**2-rmp*rmp                                     
       c=4.*w*w*(na1*rmp)**2 +2.*wp*qv2-qv2*qv2-wp*wp                   
       b=qv*(4.*wp-4.*qv2)                                              
       a=4.*w*w-4.*qv2                                                  
       rad = b*b - 4.*a*c                                               
       if(rad .lt. 0.0)return                                           
       yp1 = (-b-sqrt(rad))/(2.*a)                                      
       yp2 = (-b+sqrt(rad))/(2.*a)                                      
       p = yp2                                                          
       if(p .ge. 0.0) p = abs(p)                                        
       y = p                                                            
                                                                        
       return                                                           
       end                                                              
                                                                        
C=========================================================              
       subroutine sigcc(e0m,e1m,tht,sigp,sign)                          
C gives the elementary sigm ep and en cross section for                 
C moving nucleon                                                        
      common/sigc /qx,qy,qz,q2,px,py,pz,p2,tpq,                         
     > pxf,pyf,pzf
      COMMON/IKK12/IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      INTEGER IG,IDUT,INEL_MODEL,PAULI_MODEL,NUC_METHOD,NUC_MODEL
      REAL*8 GEP,GEN,GMP,GMN
                                                                        
      pm=938.28                                                          
      cost=cos(tht*3.1416/180.)                                         
      sint=sin(tht*3.1416/180.)                                         
      fscnst = 1./137.04                                                
      hbc=197.3289                                                      
      pm2=pm*pm                                                         
      tpm1=1./2./pm                                                     
      teomc=2.*e0m*(1.-cost)                                            
      tan2=(1.-cost)/(1.+cost)                                          
      sig=5./137.04**2/(1.-cost)/e0m**2/tan2                            
      qm2=teomc*e1m                                                     
      ep=sqrt(pm2+p2)                                                   
      epf=sqrt(pm2+p2+q2+tpq)                                           
      qm2k=q2-(epf-ep)**2                                               
      rq2=qm2/q2                                                        
      v2=-qm2*tpm1**2                                                   
      omv2=1./(1.-v2)                                                   
      pt2=p2-tpq*tpq*0.25/q2                                            
      qmp=qm2/1000000.                                                  
      call nform(IG,DBLE(qmp),gep,gen,gmp,gmn)                                   
      f1p=gep-v2*gmp                                                    
      f1n=gen-v2*gmn                                                    
      cf2p=gmp-gep                                                      
      cf2n=gmn-gen                                                      
      pstcaq=qx*py-qy*px                                                
      ve=((ep+epf)*rq2*0.5-sqrt((rq2+tan2)/q2)*pstcaq)**2               
     >    +TAN2*PZ**2                                                   
      vm=((rq2*0.5+tan2)*qm2k-rq2*0.5*qm2)*0.5                          
c     efpf=pxf*cost+pyf*sint                                            
      sigm=sig*omv2**2/ep/epf !(epf-efpf) in ingo's version             
      sigp=sigm*(ve*(f1p**2+qm2k*(tpm1*cf2p)**2)+vm*(f1p+cf2p)**2)      
      sign=sigm*(ve*(f1n**2+qm2k*(tpm1*cf2n)**2)+vm*(f1n+cf2n)**2)      
      sigp=sigp*hbc*hbc                                                 
      sign=sign*hbc*hbc                                                 
      return                                                            
      end                                             
