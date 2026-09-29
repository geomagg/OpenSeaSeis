c-----------------------------------------------------------------------
c     Rotina copiada sem alteracao de xvspmig2.f (SEISPAK).
c     Usada pelo modulo SPK_XVSPMIG2.
c-----------------------------------------------------------------------
      subroutine xsgimrg2r(z,u,u1,y,y1,shft,nl,nh,nt1,nx,itzr) 
      complex a,b,c,s1,s2,s3,s4,s5,s6
      real u(*),y(*),u1(*),y1(*),z(*),shft(*)
      pi=3.1415927
      rz=1./itzr
      nws=nh-nl+1 
c     delw=2.*pi/nt1
      delw=4.*pi/nt1
      nl1=nl+1
      wf=(nl-1)*delw
      do 200 i=1,nx*itzr
      z(i)=0.
 200  continue
c    ***    initialize frequencies
      do 100 i=1,nx 
      sr1=1.
      sr2=1. 
      sr3=cos(-delw*shft(i)) 
      sr4=cos(-wf*shft(i)) 
      sr5=cos(delw*rz*shft(i))  
      sr6=cos(wf*rz*shft(i)) 
      si1=0.
      si2=0. 
      si3=sin(-delw*shft(i)) 
      si4=sin(-wf*shft(i)) 
      si5=sin(delw*rz*shft(i))  
      si6=sin(wf*rz*shft(i)) 
      do 20 k=1,nws
      ar=sr2
      ai=si2
      br=sr4
      bi=si4
      cr=sr6
      ci=si6
      ll=0
      jr=2*((k-1)*nx+i)-1
      ji=jr+1 
      ujr=u(jr)*u1(jr)+u(ji)*u1(ji)
      uji=u(ji)*u1(jr)-u(jr)*u1(ji)
      yjr=y(jr)*y1(jr)+y(ji)*y1(ji)
      yji=y(ji)*y1(jr)-y(jr)*y1(ji)
      do 10 l=0,itzr-1
      f1=(itzr-l)*rz
      f2=1.-f1
      z(ll+i)=z(ll+i)+f1*(ar*ujr-ai*uji)+f2*(br*yjr-bi*yji)
c     z(ll+i)=z(ll+i)+f1*(real(a)*u(jr)-aimag(a)*u(ji) )
c    *               +f2*(real(b)*y(jr)-aimag(b)*y(ji) )
c     a=a*c
c     b=b*c
      tmp=ar*cr-ai*ci
      ai=ar*ci+ai*cr
      ar=tmp
      tmp=br*cr-bi*ci
      bi=br*ci+bi*cr
      br=tmp
      ll=ll+nx
10    continue
c   ***   update next frequency  
c     s2=s2*s1 
c     s4=s4*s3  
c     s6=s6*s5 
      tmp=sr2*sr1-si2*si1
      si2=sr2*si1+si2*sr1
      sr2=tmp
      tmp=sr4*sr3-si4*si3
      si4=sr4*si3+si4*sr3
      sr4=tmp
      tmp=sr6*sr5-si6*si5
      si6=sr6*si5+sr5*si6
      sr6=tmp
20    continue
100   continue  
      return
      end 
