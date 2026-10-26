c     SEASEIS: analise de coerencia (velocity focussing) de xvspmig1.f
c     (FOR LOR RNI WINDOW VMIN VMAX EPS), para o SPK_XVSPMIG1.
c     SPKV1FOC: trecho final do programa principal de xvspmig1.f (apos o
c       laco 1000): TDAT**2, boxflt, t2d2, transformada de inclinacao,
c       OFFLINE4, VAVG, tau2v2, picks (vtps), OFFLINE5; mesmos indices
c       do vetor s do original (lhdr=2, i1, i2, i3).
c     tau2v2 (xvspmig1.f) e tau2v2m (= tau2v2 de xvspmig2.f, sem zerar
c       os primeiros 200 ms): copiados, so a listagem dos picks removida.
c     boxflt, spl1: fontes da biblioteca SEISPAK; vshift1 como no
c       SPK_DVSMIG3D.
c     t2d2, relmax: reconstruidas do codigo de maquina da biblioteca
c       SEISPAK (binarios do SEISPAK), sem a fonte.
c
      SUBROUTINE SPKV1FOC(TDAT,IFITV,NT,SR,NZ2,DZ,IFOR,IRNI,JNT,NGSN,
     :  VMN,VMX,SEPS,BETA,NNZT0,OTAU,OVC,OVAVG,PICKS,IPROG)
c     IPROG 1: xvspmig1 (tau2v2), 2: xvspmig2 (tau2v2m)
      REAL TDAT(*),IFITV(*),OTAU(*),OVC(*),OVAVG(NT,*),PICKS(201,*)
      REAL, ALLOCATABLE :: S(:),DATA(:),VEL(:)
      CHARACTER*4 VINT
      vint='no  '
      dt=.001*sr
      lhdr=2
      nnzt=max0(nnzt0,nt)
      i1=lhdr+1+nnzt
      i2=i1+nnzt
      i3=i2+nnzt
      allocate(s(max0(nt*jnt,i3+nz2+nnzt)+4000),data(2*nt*jnt+100),
     :         vel(nt+1000))
      iisr=8./sr+.1
      jnt1=(jnt-1)/iisr+1
      do j=1,nz2*jnt*ngsn
      tdat(j)=tdat(j)**2
      enddo
      do i=1,ngsn
      ii=(i-1)*nz2*jnt+1
c   ***   filter semblance data
      call rtran(tdat(ii),data,jnt,nz2)
      call boxflt(data,tdat(ii),s,nz2,jnt,3)
      mjnt=jnt/2+1
      do j=1,jnt
      j1=(j-1)*nz2+ii
      j2=(j-1)*nt+1
      call movlev(tdat(j1),s(i1),nz2)
      dpn=ifor+irni*(i-1)
      call t2d2(s,s(i1),s(i2),s(i3),ifitv,nt,sr,nz2,dz,dpn,beta,vint)
c   ***   apply skewing transform
      if(mjnt-j.gt.0)then
      call clear(s(i1+nt-mjnt+j),mjnt-j)
      call movlev(s(mjnt-j+1),s(i1),nt-mjnt+j)
      end if
      if(mjnt-j.eq.0)call movlev(s,s(i1),nt)
      if(mjnt-j.lt.0)then
      call clear(s(i1),-mjnt+j)
      call movlev(s,s(i1-mjnt+j),nt+mjnt-j)
      end if
      call movlev(s(i1),data(j2),nt)
      end do
      jj1=1+nt*jnt
      call boxflt(data,data(jj1),s,nt,jnt,5)
      call movlev(data(jj1),data,nt*jnt)
c   ***   tau-t data (OFFLINE4); tracos jnt/2 e jnt/2+1 marcados -.05
      do j=1,jnt
      j1=(j-1)*nt+1
      if(mod(j-1,iisr).eq.0) then
      itrc=(j-1)/iisr + 1
      io=((i-1)*jnt1+itrc-1)*nt
      do n=1,nt
      otau(io+n)=data(j1+n-1)
      enddo
      if(j.eq.jnt/2 .or. j.eq.jnt/2+1)then
      do n=1,nt
      otau(io+n)=-.05
      enddo
      endif
      endif
      end do
c   ***   compute average velocity from depths returned in s(i2) from t2d
      do j=1,nt-1
      vel(j+1)=2.*dz*s(i2+j)/(j*dt)
      end do
      vel(1)=vel(2)
      do n=1,nt
      ovavg(n,i)=vel(n)
      enddo
      tx=1.
      if(iprog.eq.1)then
      call tau2v2(data,s,s(i1),s(i2),vel,tx,sr,seps,vmn,vmx,nt,jnt)
      else
      call tau2v2m(data,s,s(i1),s(i2),vel,tx,sr,seps,vmn,vmx,nt,jnt)
      endif
c   ***   SEASEIS: xmax(1)=umax (s(i2)) cai dentro do painel de coerencia
c   ***   ja transposto em s; zera para nao gerar um spike no OFFLINE5
      s(i2)=0.
c   ***   picks (vtps): picks(1)=n, (vel,tempo) em picks(2..)
      npk=min0(100,nint(vel(1)))
      picks(1,i)=npk
      do j=1,npk
      picks(2*j,i)=vel(2*j)
      picks(2*j+1,i)=vel(2*j+1)
      enddo
c   ***   velocity coherence (OFFLINE5)
      do j=1,jnt
      j1=(j-1)*nt+1
      if(mod(j-1,iisr).eq.0) then
      itrc=(j-1)/iisr + 1
      io=((i-1)*jnt1+itrc-1)*nt
      do n=1,nt
      ovc(io+n)=s(j1+n-1)
      enddo
      endif
      end do
      end do
      deallocate(s,data,vel)
      RETURN
      END
c
      subroutine t2d2(u,x,e,s,v,nt,sr,nnz,ddz,dpn,beta,vint)
c   ***   SEASEIS: reconstruida do codigo de maquina de t2d2_ (SEISPAK)
c   ***   traco em profundidade x (nnz amostras de ddz) -> tempo u (nt
c   ***   amostras de sr ms) com as funcoes v(z) de v (pares vz, como
c   ***   cacomt2); e(n) = profundidade (em amostras) no tempo n
      real u(*),x(*),e(*),s(*),v(*)
      character*4 vint
      nxv=v(1)
      nyv=v(2)
      np=v(3)
      mindpn=v(5)
      idpn=dpn-mindpn+beta-1.
      idpn=max0(1,min0(nxv,idpn))
      ii=6+(idpn-1)*(2*np+3)
      ii0=ii+2*np
      rddz=1./ddz
      dt2=.0005*sr
111   continue
      if(v(ii).eq.0.)then
      ii=ii+2
      goto 111
      end if
      v1=v(ii)
      v2=v1
      z=0.
      z1=0.
      z2=v(ii+1)
      t=0.
      eps=.1
      dzr=1./(z2-z1+eps)
      do 101 n=1,nt
102   continue
      if(ii.gt.ii0.or.z.le.z2)goto 103
      if(v(ii).eq.0. .or. z .gt. v(ii+1))then
      ii=ii+2
      goto 102
      end if
      v1=v2
      v2=v(ii)
      z1=z2
      z2=v(ii+1)
      goto 102
103   continue
      vel=v2
      if(vint.eq.'yes' .and. ii.lt.ii0 )
     *vel=((z-z1)*v2+(z2-z)*v1)*dzr
      e(n)=z*rddz
      z=z+dt2*vel
      t=t+sr
101   continue
      call spl1(x,1.,u,e,s,nnz,nt,1.)
      return
      end
c
      subroutine relmax(u,n,m,eps,nmax,xmax)
c   ***   SEASEIS: reconstruida do codigo de maquina de relmax_ (SEISPAK)
c   ***   normaliza u(n,m) pelo maximo interno e acha maximos locais
c   ***   (> eps, vizinhanca +-2) e depois maximos numa janela +-13;
c   ***   xmax(3k-1..3k+1) = valor, i, j ; xmax(1) = maximo; ate 100
      real u(n,*),xmax(*)
      nmax=0
      umax=1.e-6
      do j=4,m-3
      do i=4,n-3
      umax=amax1(umax,u(i,j))
      end do
      end do
      scl=1./umax
      do j=1,m
      do i=1,n
      u(i,j)=u(i,j)*scl
      end do
      end do
      do j=4,m-3
      do i=4,n-3
      a=u(i,j)
      if(a.gt.eps .and. a.ge.u(i-2,j) .and. a.ge.u(i-1,j) .and.
     *   a.ge.u(i+1,j) .and. a.ge.u(i+2,j) .and. a.ge.u(i,j-2) .and.
     *   a.ge.u(i,j-1) .and. a.ge.u(i,j+1) .and. a.ge.u(i,j+2))then
      nmax=nmax+1
      xmax(3*nmax-1)=a
      xmax(3*nmax)=i
      xmax(3*nmax+1)=j
      end if
      end do
      end do
      nn=0
      do k=1,nmax
      i=xmax(3*k)
      j=xmax(3*k+1)
      i1=max0(1,i-13)
      i2=min0(n,i+13)
      j1=max0(1,j-13)
      j2=min0(m,j+13)
      am=0.
      do jj=j1,j2
      do ii=i1,i2
      am=amax1(am,u(ii,jj))
      end do
      end do
      if(u(i,j).ge.am)then
      nn=nn+1
      xmax(3*nn-1)=am
      xmax(3*nn)=i
      xmax(3*nn+1)=j
      if(nn.eq.100)goto 99
      end if
      end do
99    continue
      xmax(1)=umax
      nmax=nn
      return
      end
c
      subroutine tau2v2(u,s,imax,xmax,vel,tx,sr,eps,vmin,vmax,nt,ntr)
      real s(*),u(*),xmax(*),vel(*)
      integer imax(*)
      call rtran(u,s,nt,ntr)
c   set top 200 mls to 0
      ii=ntr*nint(200/sr)
      call clear(s,ii)
      dv=(vmax-vmin)/(ntr-1)
      mntr=ntr/2
      do 40 i=1,nt
      i1=(i-1)*ntr
      do 50 j=1,ntr
      v=vmin+(j-1)*dv
      r2=(vel(i)/v)**2
      r2=vel(i)/v
      tt=(i-1)*sr
c     call phitab(tx,tt,r2,gam)
      gam=.8*tt*(1. - v/vel(i) )/sr
      x=gam+mntr+1
      ix=x
      if(1.lt.ix .and. ix.lt.ntr)then
      f1=x-ix
      f2=1.-f1
      u(i1+j)=f2*s(i1+ix)+f1*s(i1+ix+1)
      else
      u(i1+j)=0.
      end if
50    continue
c   ***   
40    continue
      call relmax(u,ntr,nt,eps,nmax,xmax)
      aamax=xmax(1)
c   ***   set gam=0 location in velocity-time plot to 0.
      do i=1,nt
      i1=(i-1)*ntr
      ii=(vel(i)-vmin)/dv + 1.5
      u(i1+ii)=0.
      end do  
c     SEASEIS (listagem): write(*,*)':::',nmax,' SIGNIFICANT LOCAL COHERENCE MAXIMA FOUND'
c     SEASEIS (listagem): write(*,*)':::'
c     SEASEIS (listagem): write(*,19)':::','VEL','TIME','TAU','POW','DEPTH'
      t2=0.
      velint=vel(nt+1)
      do i=1,nmax
      t1=t2
      t2=(xmax(3*i+1)-1.)*sr
      vel(nt+i)=vmin+(xmax(3*i)-1.)*dv
      if(t2-t1.gt.sr)then
      velint=(t2*vel(nt+i) - t1*vel(nt+i-1) )/(t2-t1)
      else
      velint=0.
      end if
      i1=t2+.5
      i2=vel(nt+i)+.5
      i3=velint+.5
      i3=.0005*vel(nt+i)*t2+.5
      i4=1000.*xmax(3*i-1)
      i6=t2/sr+1.5
      r2=(vel(i6)/vel(nt+i))**2
c     call phitab(tx,t2,r2,gam)
      gam=t2*(1.-1./r2)
c     SEASEIS (listagem): ttx=t2/tx
c     SEASEIS (listagem): gam=gam*fudge1(ttx)
c     SEASEIS (listagem): i5=gam+.5
c     SEASEIS (listagem): write(*,29)':::',i2,i1,i5,i4,i3
c     SEASEIS (listagem): 19    format(a3,5a10)
c     SEASEIS (listagem): 29    format(a3,5i10)
      end do
c   ***   return nmax in vel(1) and nmax vtp pairs in vel(2)->vel(2*nmax+1)
      vel(1)=nmax
      do i=1,nmax
      vel(2*i)=vel(nt+i)
      vel(2*i+1)=(xmax(3*i+1)-1.)*sr
      end do
10    continue
c   ***   transpose velocity file 
      call rtran(u,s,ntr,nt)
      xmax(1)=aamax
      return
      end
      subroutine tau2v2m(u,s,imax,xmax,vel,tx,sr,eps,vmin,vmax,nt,ntr)
      real s(*),u(*),xmax(*),vel(*)
      integer imax(*)
      call rtran(u,s,nt,ntr)
      dv=(vmax-vmin)/(ntr-1)
      mntr=ntr/2
      do 40 i=1,nt
      i1=(i-1)*ntr
      do 50 j=1,ntr
      v=vmin+(j-1)*dv
c     r2=(vel(i)/v)**2
      r2=(vel(i)/v)
      tt=(i-1)*sr
c     call phitab(tx,tt,r2,gam)
c   ***   factor for calibration
c   ***
      factor=.5
      gam=factor*tt*(1.-1./r2)
c   ***
c   ***
c     SEASEIS (listagem): ttx=tt/tx
      gam=gam/sr
      x=gam+mntr+1
      ix=x
      if(1.lt.ix .and. ix.lt.ntr)then
      f1=x-ix
      f2=1.-f1
      u(i1+j)=f2*s(i1+ix)+f1*s(i1+ix+1)
      else
      u(i1+j)=0.
      end if
50    continue
c   ***   
40    continue
      call relmax(u,ntr,nt,eps,nmax,xmax)
      aamax=xmax(1)
c   ***   set gam=0 location in velocity-time plot to 0.
      do i=1,nt
      i1=(i-1)*ntr
      ii=(vel(i)-vmin)/dv + 1.5
      u(i1+ii)=0.
      end do  
c     SEASEIS (listagem): write(*,*)':::',nmax,' SIGNIFICANT LOCAL COHERENCE MAXIMA FOUND'
c     SEASEIS (listagem): write(*,*)':::'
c     SEASEIS (listagem): write(*,19)':::','VEL','TIME','TAU','POW','DEPTH'
      t2=0.
      velint=vel(nt+1)
      do i=1,nmax
      t1=t2
      t2=(xmax(3*i+1)-1.)*sr
      vel(nt+i)=vmin+(xmax(3*i)-1.)*dv
      if(t2-t1.gt.sr)then
      velint=(t2*vel(nt+i) - t1*vel(nt+i-1) )/(t2-t1)
      else
      velint=0.
      end if
      i1=t2+.5
      i2=vel(nt+i)+.5
      i3=velint+.5
      i3=.0005*vel(nt+i)*t2+.5
      i4=1000.*xmax(3*i-1)
      i6=t2/sr+1.5
      r2=(vel(i6)/vel(nt+i))**2
c     call phitab(tx,t2,r2,gam)
      gam=t2*(1.-1./r2)
c     SEASEIS (listagem): ttx=t2/tx
c     SEASEIS (listagem): gam=gam*fudge1(ttx)
c     SEASEIS (listagem): i5=gam+.5
c     SEASEIS (listagem): write(*,29)':::',i2,i1,i5,i4,i3
c     SEASEIS (listagem): 19    format(a3,5a10)
c     SEASEIS (listagem): 29    format(a3,5i10)
      end do
c   ***   return nmax in vel(1) and nmax vtp pairs in vel(2)->vel(2*nmax+1)
      vel(1)=nmax
      do i=1,nmax
      vel(2*i)=vel(nt+i)
      vel(2*i+1)=(xmax(3*i+1)-1.)*sr
      end do
10    continue
c   ***   transpose velocity file 
      call rtran(u,s,ntr,nt)
      xmax(1)=aamax
      return
      end
      subroutine vshift1(y,a,b,s,shft,dw,nl,nh,nx,sgn)
      complex y(*),s(*),a(*),b(*)
      real shft(*)
c   ***   efficient  complex exponential table computation
c   ***   a and b are equivalent.
      w1=sgn*(nl-1)*dw
      dws=sgn*dw
      nws=nh-nl+1
      if(nx.gt.1)then
      do 10 i=1,nx
      s(i)=cexp(cmplx(0.,dws*shft(i)))
10    a(i)=cexp(cmplx(0.,w1*shft(i)))
      do 20 j=1,nws
      ii=(j-1)*nx
      do 20 i=1,nx
      y(ii+i)=a(i)*y(ii+i)
      a(i)=s(i)*a(i)
20    continue
      else
      log2=alog(real(nws)-1.)/alog(2.)
      s(1)=cexp(cmplx(0.,dws*shft(1)))
      s(2)=cexp(cmplx(0.,w1*shft(1)))
      a(1)=cmplx(1.,0.)
      a(2)=s(1)
      i0=1
      do l=1,log2
      i0=2*i0
      s(1)=s(1)*s(1)
      i1=min0(i0,nws-i0)
      do i=1,i1
      b(i0+i)=s(1)*a(i)
      end do
      end do
      do i=1,nws
      y(i)=s(2)*b(i)*y(i)
      end do
      end if
      return
      end
      subroutine boxflt(s1,u,s2,nt,nx,l)
      real u(*),s1(*),s2(*)
c   ***   clear work array
      do j=1,nt
      s2(j)=0.
      end do
      x=1./real(2*l+1)
      do i=1,l
      j1=(i-1)*nt
      do j=1,nt
      s2(j)=s2(j)+s1(j1+j)
      end do
      end do
      do i=1,nx
      ii=(i-1)*nt
      j1=i-l-1
      j2=i+l
      if(j1.lt.1)then
      do j=1,nt
      s2(j)=s2(j)+s1((j2-1)*nt+j)
      u(ii+j)=x*s2(j)
      end do
      end if
      if(j2.gt.nx)then
      do j=1,nt
      s2(j)=s2(j)-s1((j1-1)*nt+j)
      u(ii+j)=x*s2(j)
      end do
      end if
      if(j1.ge.1 .and. j2.le.nx)then
      do j=1,nt
      s2(j)=s2(j)+s1((j2-1)*nt+j)-s1((j1-1)*nt+j)
      u(ii+j)=x*s2(j)
      end do
      end if
      end do
      return
      end
      subroutine spl1(y,hy,f,xf,s,n,m,scf)
c  resample n samples of y whose sample rate is hy to m samples of f
c   at xf(i),i=1,m .  where  0 <= xf(i) <= xf(i+1) <= (n-1)*hy
      real xf(*),f(*),y(*), s(*) 
      data b1/.25/,b2/.2667/,b3/.2679/ 
c   compute spline coefficients  
      r=6./(hy*hy) 
      s(1)=0.  
      s(2)=b1*r*(y(3)-2*y(2)+y(1)) 
      s(3)=b2*(r*(y(4)-2*y(3)+y(2))-s(2))  
      do 10 i=4,n-1  
10    s(i)=b3*(r*(y(i+1)-2*y(i)+y(i-1)) -s(i-1)) 
      s(n)=0.  
      do 20 i=n-2,4,-1 
20    s(i)=s(i)-b3*s(i+1)  
      s(3)=s(3)-b2*s(4)  
      s(2)=s(2)-b1*s(3)  
c   compute polynomial coefficients and interpolate  
      r1=1./hy 
      r2=1./(6.*hy)  
      r3=hy/6. 
      j=1  
      xy=0.  
      a=s(2)*r2  
      b=0. 
      c=(y(2)-y(1))*r1-s(2)*r3 
      do 30 i=1,m  
13    d1=xf(i)-xy 
      d2=hy-d1 
      if(d1.ge.0..and.d2.ge.0.)goto 12 
      j=j+1  
      xy=xy+hy 
      if(j.eq.n)goto 40  
      a=(s(j+1)-s(j))*r2 
      b=.5*s(j)  
      c=(y(j+1)-y(j))*r1-(2.*s(j)+s(j+1))*r3 
      goto 13  
12    continue 
      e=d1*(c+d1*(b+d1*a)) 
      f(i)=scf*(y(j)+e)  
30    continue 
40    continue
      if(i.lt.m)then
      do 25 j=i,m
      f(j)=0.
25    continue
      end if
      return
      end  
