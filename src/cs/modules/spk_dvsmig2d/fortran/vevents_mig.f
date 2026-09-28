c-----------------------------------------------------------------------
c     VEVENTS e rotinas auxiliares, copiadas de dvslib.f (SEISPAK) sem
c     alteracao, exceto onde marcado "SEASEIS:". Leitura da lista VXT
c     (horizontes) para o modulo SPK_DVSMIG2D.
c-----------------------------------------------------------------------
      subroutine vevents(ee,e,xx,x,nx,nny,nt,vint,nwds,ia,ib,ddd) 
      dimension ee(*),e(*),xx(*),x(*),nevents(100),ntodo(100)  
      dimension mode(100)  
      real ia(*),ib(*)
      common/rdd/rdd 
      save
      character*4 topdwn,vint,depth,vtype
      data topdwn/'yes'/,depth/'no'/,vtype/'int'/
      vint='no'
      nm=0 
      ifwb=1 
      nww=nwds 
      if(e(1).ge.10.)topdwn='no' 
      if(e(1).ge.10.)e(1)=e(1)-10. 
      write(*,*)'topdwn =',topdwn
      rdd=1. 
      if(e(1).ge.5.)then
      rdd=1./ddd
      e(1)=e(1)-4.
      end if
      iloc=2 
      if(e(1).eq.3.5 .or. e(1).eq.4.5)vtype='avg'
      if(e(1).eq.3.1 .or. e(1).eq.4.1)vtype='rms'
      if(e(1).le.2.1)depth='yes' 
      if(e(1).eq.2..or.e(1).eq.4.)vint='yes' 
      if(e(2).lt.0.)iloc=iloc+1  
      if(e(2).lt.0.)ill=-e(2)  
      if(e(3).lt.0.)iloc=iloc+1  
      if(e(3).lt.0.)irr=-e(3)  
      i=irr  
      if(ill.gt.0.and.i.eq.0)irr=ill 
      if(ill.gt.0.and.i.eq.0)ill=1 
      if(e(iloc).lt.1.1)iloc=iloc+1  
      igrp=0 
      ee(1)=1. 
200   continue 
      igrp=igrp+1  
      call movlev(e(iloc),ee(2),nww) 
      if(iloc.gt.1 .and. igrp.gt.1 )then
      ee(1)=e(iloc-1)
      write(*,*)'ny = ',ny,' at igrp = ',igrp
      end if
      ny=ee(1)
      nsects=1 
      icheck=icheck+1 
      if(icheck.eq.2)go to 1111  
      write(*,*)' vxt velocity specification edit check' 
      k=1  
      dpnmax=0.  
      dpnmin=1000000.  
      vmax=0.  
      vmin=1000000.  
      x1=0.  
      anx=0. 
      i=1  
      vv=ee(2) 
      x2=ee(3) 
      if(rdd.ne.1.)x2=aint(x2*rdd+1.5) 
      if(rdd.ne.1.)e(iloc+1)=x2  
      if(rdd.ne.1.)ee(3)=x2  
      tt=ee(4) 
      dpnmax=amax1(dpnmax,x2)  
      dpnmin=amin1(dpnmin,x2)  
      vmax=amax1(vmax,vv)  
      vmin=amin1(vmin,vv)  
      i=4  
      j=5  
  80  continue 
      vv=ee(j) 
      x2=ee(j+1) 
      if(rdd.ne.1.)x2=aint(x2*rdd+1.5) 
      if(rdd.ne.1.)e(iloc+j-1)=x2  
      if(rdd.ne.1.)ee(j+1)=x2  
      tt=ee(j+2) 
      anx=amax1(anx,x2)  
      dpnmax=amax1(dpnmax,x2)  
      dpnmin=amin1(dpnmin,x2)  
      vmax=amax1(vmax,vv)  
      vmin=amin1(vmin,vv)  
      if(1.+x1.le.x2)go to 1100  
      write(*,*)'dpns not strictly increasing in horizon ',k 
1100    continue 
      i=i+3  
      j=j+3  
      x1=x2  
      if(ee(j).gt..5)go to 80  
      j=j+1  
      k=k+1  
      x1=0.  
      if(ee(j).gt..5)go to 80  
      k=k-1  
      write(*,*)' number of events found = ',k 
      write(*,*)' minimum dpn found = ',dpnmin 
      write(*,*)' maximum dpn found = ',dpnmax 
      write(*,*)' minimum velocity found = ',vmin  
      write(*,*)' maximum velocity found = ',vmax  
      if(vint.eq.'no')write(*,*)'velocity constant between horizons' 
      if(vint.eq.'yes')write(*,*)'velocity varies between horizons'  
      nx=dpnmax-dpnmin+1.1 
      if(icheck.eq.1)return  
1111   continue  
      nevents(igrp)=1  
      mindpn=dpnmin+.1 
      maxdpn=dpnmax+.1 
      do 99 i=1,nx 
      xx(i)=mindpn+i-1 
99    continue 
      i=2  
1     continue 
      if(ee(i).le.0.)go to 2 
      i=i+3  
      go to 1  
2     continue 
      mode(nevents(igrp))=ee(i)  
      if(ee(i+1).eq.0.)go to 12  
      i=i+1  
      nevents(igrp)=nevents(igrp)+1  
      go to 1  
 12   continue 
      nwords=i+1 
      nsects=1 
      iy1=1  
      iy2=iy1+nwords 
      if(ny.eq.1)go to 5 
3     continue 
      fy=e(iy1)  
      if(1.+ee(iy1).ge.ee(iy2))go to 5 
      iy1=iy1+nwords 
      iy2=iy1+nwords 
      nsects=nsects+1  
      go to 3  
5     continue 
      nwds2=2*nwords 
      iloc=iloc+iy2-1  
      if(ny.gt.1)call trans(e,ee,nwords,nsects)  
      if(ny.eq.1)call movlev(ee,e,nwords)  
      nwm=nwords-1 
      i1=nsects+1  
      i2=1 
      iy=e(1)  
      ntodo(igrp)=e(nsects)-e(1)+1 
      ntd=ntodo(igrp)  
      if(ny.eq.1 .or. ntd.eq.1)goto 21
      l1=ntd+1 
      l2=l1+ntd  
      l3=l2+ntd  
      l4=l3+ntd  
      do 20 i=1,nwm  
      call clear(x,ntd)  
      if(ntd.gt.1)then 
      if(mode(i).eq.-1)then
      call splyn(e,xx(iy),e(i1),x,x(l1),x(l2),x(l3),x(l4),nsects,ntd)
      else
      if(mode(i).eq.0)mode(i)=-2
      call joehint(e,e(i1),xx(iy),x,x(l1),nsects,ntd,mode(i))
      end if
      end if
      call movlev(x,ee(i2),ntd)  
      i1=i1+nsects 
      i2=i2+ntd  
20    continue 
      call trans(x,ee,ntd,nwm) 
21    continue 
      if(ny.eq.1 .or. ntd .eq. 1)call movlev(e(2),x,nwords)  
      kk=1 
      ifw=1  
      l1=ntd*nwords+1  
      l2=l1+500  
      l3=l2+500  
      l4=l3+500  
      l5=l4+2*nx 
      do 40 i=1,ntd  
      iyy=iyy+1  
      write(*,*)'      y = ',iyy 
      nnn=nevents(igrp)  
      do 41 j=1,nnn  
c     if(mode(j).eq.0)write(*,1101)j 
c     if(mode(j).eq.-1)write(*,1102)j  
c     if(mode(j).le.-2)write(*,1103)j  
c     if(mode(j).eq.-3)write(*,1104)j  
1104  format(' interpolation of velocities for event ',i3,' is linear')  
1101  format('  interpolation of event',i3,' by linear interpolation') 
1102  format('  interpolation of event',i3,' by spline interpolation') 
1103  format('  interpolation of event',i3,' by joeh   interpolation') 
c      if(mode(j).eq.0)mode(j)=-2
      call interpv(x(kk),x(l1),x(l2),x(l3),x(l4),xx,x(l5),np,nx,mode(j)) 
      kk=kk+3*np+1 
c     call p(ia,ifw,2*nx,x(l4))  
      call movlev(x(l4),ia(ifw),2*nx)
      ifw=ifw+2*nx 
 41   continue 
      kk=kk+1  
40    continue 
c****    sparse velocity specification has been made dense.  
      ifw=1  
      do 31 i=1,ntd  
      call trnsp(ia,ifw,ib,ifwb,nx,2*nevents(igrp),x,ee,10000,1.)  
      ifw=ifw+2*nx*nevents(igrp) 
      ifwb=ifwb+2*nx*nevents(igrp) 
 31   continue 
      nww=nww-nwords 
      nm=max0(nm,nevents(igrp))  
      if(iyy.lt.ny)go to 200 
c****    trace headers are affixed and velocity function is  
c****    stored in form accepted by velcom in wem. 
      ee(1)=nx 
      ee(2)=ny 
      ee(3)=nm 
      ifw1=1 
c     call p(ia,ifw1,3,ee) 
      call movlev(ee,ia(ifw1),3)
      ifw1=ifw1+3  
      ifw2=1 
      ii=1 
      do 100 n=1,igrp  
      nn1=nx*ntodo(n)  
      do 81 k=1,nn1  
      nwds=nevents(n)  
      iy=(ii-1)/nx+1 
      ix=ii+mindpn-1-(iy-1)*nx 
      ee(1)=iy 
      ee(2)=ix 
c     call g(ib,ifw2,2*nwds,x) 
      call movlev(ib(ifw2),x,2*nwds)
      ifw2=ifw2+2*nwds 
c****    get velocity-depth pairs into increasing depth order  
      iv1=nevents(n)-1 
      if(iv1.eq.0)go to 54 
      do 50 j=1,iv1  
      iv2=nevents(n)-j 
      do 51 jj=1,iv2 
      z0=x(2*jj)
      z1=x(2*jj+2)
      if(z0.eq.0. .or. (z1.lt.z0 .and. topdwn .eq. 'no'))then
      if(z1.gt.0.)then
      tmp=x(2*jj-1)  
      x(2*jj-1)=x(2*jj+1)  
      x(2*jj+1)=tmp  
      tmp=x(2*jj)  
      x(2*jj)=x(2*jj+2)  
      x(2*jj+2)=tmp  
      end if
      end if
53    continue 
51    continue 
50    continue 
54    continue 
c****    unpack velocity-depth pairs and put away. 
      ntoput=2*nm+3  
      jj=3 
      do 60 i=1,nwds 
      ee(jj)=x(2*i-1)  
      ee(jj+1)=x(2*i)  
      jj=jj+2  
60    continue 
61    continue 
      do 62 l=jj,ntoput  
      ee(l)=0. 
62    continue 
      if(depth.eq.'yes')go to 900  
      if(vtype.eq.'rms')then
      mm=3
1800   continue    
      if(ee(mm).eq.0.)then
      mm=mm+2
      goto 1800
      end if
      v1=ee(mm)
      t1=ee(mm+1)
      vel=v1
      do 1801 j=mm,ntoput-2,2
      v2=ee(j)
      t2=ee(j+1)
      if(t2-t1.gt. 20.)then
      vel=(t2*v2**2-t1*v1**2)/(t2-t1)
      vel=sqrt(vel)
      v1=v2
      t1=t2
      end if
      ee(j)=vel
1801   continue
      end if
      if(vtype.eq.'avg')then
      mm=3
800   continue    
      if(ee(mm).eq.0.)then
      mm=mm+2
      goto 800
      end if
      v1=ee(mm)
      t1=ee(mm+1)
      vel=v1
      do 801 j=mm,ntoput-2,2
      v2=ee(j)
      t2=ee(j+1)
      if(t2-t1.gt. 20.)then
      vel=(v2*t2-v1*t1)/(t2-t1)
      v1=v2
      t1=t2
      end if
      ee(j)=vel
801   continue
      end if
      mm=3 
      z=0. 
      t1=0.  
      v1=ee(mm)  
      do 901 l=1,nwds  
      v2=ee(mm)  
      if(v2.lt..01)go to 900 
      t2=ee(mm+1)  
      ddt=.0005*(t2-t1)  
      v=v2 
      if(vint.eq.'yes'.and.abs(v2-v1).gt..1)v=(v2-v1)/alog(v2/v1)  
      z=z+v*ddt  
      ee(mm+1)=z 
      mm=mm+2  
      t1=t2  
      v1=v2  
901   continue 
 900  continue 
      ii=ii+1  
      if(ix.lt.ill.or.ix.gt.irr)go to 85 
      if(ny.eq.1)write(*,*)' x = ',int(ee(2))  
      if(ny.gt.1)write(*,*)' y = ',int(ee(1)),'   x = ',int(ee(2)) 
      nsteps=(nm-1)/12 +1  
      ll1=3  
      ll2=4  
      ll3=min0(ntoput-1,26)  
      do 86 l=1,nsteps 
      print 112,(ee(ll),ll=ll1,ll3,2)  
      print 112,(ee(ll),ll=ll2,ll3,2)  
      write(*,*)' '  
      ll1=ll1+24 
      ll2=ll2+24 
      ll3=ll3+24 
      ll3=min0(ll3,ntoput-1) 
 86   continue 
 85   continue 
112    format(5x,12f10.0)  
c     call p(ia,ifw1,ntoput,ee)  
      call movlev(ee,ia(ifw1),ntoput)
      ifw1=ifw1+ntoput
81    continue 
100    continue  
c     SEASEIS: zera contadores (SAVE) para permitir nova chamada
      icheck=0
      iyy=0
      end
      subroutine interpv(e,v,x,t,vv,xx,tt,npnts,nx,mode) 
      dimension e(*),x(*),v(*),t(*),vv(*)  
      common/rdd/rdd 
      dimension xx(1),tt(1)  
      vnm=4hv(i) 
      xnm=4hx(i) 
      tnm=4ht(i) 
      mindpn=xx(1)+.1  
      maxdpn=mindpn+nx-1 
      v(1)=e(1)  
      i=1  
2     continue 
      x(i)=e(3*i-1)  
      t(i)=e(3*i)  
      v(i+1)=e(3*i+1)  
      if(v(i+1).le.0.)go to 3  
      i=i+1  
      go to 2  
3     continue 
      npnts=i  
      npnts1=npnts-1 
      ixl=x(1)+.1  
      ixlm=ixl-mindpn+1  
      ixr=x(npnts)+.0001 
      nnx=ixr-ixl+1  
      nsteps=(npnts-1)/15 +1 
      il=1 
      ir=min0(15,npnts)  
      do 55 i=1,nsteps 
c     print 56,vnm,(v(l),l=il,ir)  
c     print 56,xnm,(x(l),l=il,ir)  
c     print 56,tnm,(t(l),l=il,ir)  
      il=il+15 
        ir=min0(ir+15,npnts) 
55    continue 
56    format(2x,a4,15f8.0) 
      l1=nx+1  
        l2=l1+npnts  
        l3=l2+npnts  
        l4=l3+npnts  
      ixl=x(1)+.1  
      ixr=x(npnts)+.0001 
      ixlm=ixl-mindpn+1  
      if(mode.le.-2.or.mode.eq.0)goto 17 
      call splyn(x,xx(ixlm),v,vv(ixlm),tt,tt(l1),tt(l2)  
     *  ,tt(l3),npnts,nnx) 
      call splyn(x,xx(ixlm),t,tt,tt(l1),tt(l2),tt(l3)  
     *  ,tt(l4),npnts,nnx) 
      go to 99 
17    continue 
      call joehint(x,v,xx(ixlm),vv(ixlm),tt(l1),npnts,nnx,mode)  
      if(mode.eq.-3)mode=-2  
c     if(mode.eq.0)mode=-2
      call joehint(x,t,xx(ixlm),tt,tt(l1),npnts,nnx,mode)  
99    continue 
      ixl=x(1)+.1 
      ixlm=ixl-mindpn+1 
      do 6 j=1,nnx 
      vv(ixlm+nx+j-1)=tt(j)  
6     continue 
      if(ixl.le.mindpn)go to 70  
      ixl1=ixl-1 
      do 61 i=mindpn,ixl1  
      vv(i-mindpn+1)=0.  
      vv(i-mindpn+nx+1)=0. 
61    continue 
70    continue 
      if(ixr.ge.maxdpn)go to 80  
      ixr1=ixr+1 
      do 71 i=ixr1,maxdpn  
      vv(i-mindpn+1)=0.  
      vv(i-mindpn+nx+1)=0. 
71    continue 
80    continue 
      return 
      end
      subroutine joehint(x,y,xx,yy,dy,mm,nn,mode)  
      dimension x(*),y(*),xx(*),yy(*),dy(*)  
      mm1=mm-1 
c   ***   compute slopes at end points 
      mmm=mm1  
      if(mode.eq.-2)mmm=1  
      do 101 i=1,mmm 
      dy(i)=(y(i+1)-y(i))/(x(i+1)-x(i))  
101   continue 
      dy(mm)=(y(mm)-y(mm-1))/(x(mm)-x(mm-1)) 
      if(mm.lt.3)go to 10  
      if(mode.eq.-3)goto 10  
      if(mode.eq.0)goto 10 
c   ***   compute slopes at interior points  
      do 100 i=2,mm1 
      x0=x(i)  
      xm=x(i-1)-x0 
      xp=x(i+1)-x0 
      fm=y(i-1)  
      f0=y(i)  
      fp=y(i+1)  
      dy(i)=0. 
      pm=fm-f0 
      pp=fp-f0 
      wm=xm**2+pm**2 
      wp=xp**2+pp**2 
      xm2=xm**2/wm 
      xm3=xm*xm2 
      xm4=xm*xm3 
      xm5=xm*xm4 
      xp2=xp**2/wm 
      xp3=xp*xp2 
      xp4=xp*xp3 
      xp5=xp*xp4 
      a1=.25*(xm4-xp4) 
      b1=.333333*(xm3-xp3) 
      g1=.333333*(pm*xm2-pp*xp2) 
      a2=.20*(xm5-xp5) 
      b2=a1  
      g2=.25*(pm*xm3-pp*xp3) 
      den=a2*b1-a1*b2  
      top=a2*g1-a1*g2  
      if(den*1000000..le.top)go to 11  
      dy(i)=top/den  
11    continue 
100   continue 
10    continue 
c   ***   interpolate  
      j=1  
      ii=1 
      do 30 i=1,nn 
3     continue 
      if(ii.ne.i)go to 31  
      delx=x(j+1)-x(j) 
      eps=.000001*delx 
c   ***   compute polynomial coefficients  
      dxr=1./delx**2 
      a1=y(j)  
      a2=dy(j) 
      a3=((y(j+1)-a1)-(x(j+1)-x(j))*a2)*dxr  
      a4=((dy(j+1)-a2)-2.*(x(j+1)-x(j))*a3)*dxr  
      if(mode.eq.-2) goto  31  
      a3=0.  
      a4=0.  
31    continue 
      del1=xx(i)-x(j)  
      del2=x(j+1)-xx(i)  
      if(del1+eps.ge.0..and.del2+eps.ge.0.)go to 2 
      if(del1.lt.0..and.j.eq.1)go to 2 
      if(j.eq.mm1)go to 2  
      j=j+1  
      ii=i 
      go to 3  
2     continue 
      yy(i)=a1+del1*(a2+del1*(a3-del2*a4)) 
30    continue 
4     continue 
      return 
      end
      subroutine splyn(x,xx,y,f,a,b,c,d,n,m) 
      dimension x(*),xx(*),y(*),f(*),a(*),b(*),c(*),d(*) 
c       cubic spline interpolation. input points are (x(i),y(i)),i=1,n.  
c     output points are (xx(i),f(i)),i=1,m. output independent variables 
c     xx(i) must safisfy x(1)@xx(i). if xx(i)>x(n)then f(i)=0. 
      n1=n-1 
      n2=n-2 
      n3=n-3 
      if(n1.gt.1)go to 10  
      den=1./(x(2)-x(1)) 
      do 9 i=1,m 
      cl=(x(2)-xx(i))*den  
      cr=(xx(i)-x(1))*den  
      f(i)=cl*y(1)+cr*y(2) 
9     continue 
      go to 15 
10    continue 
c       compute lhs and rhs of tridiagonal matrix system 
      do 1 i=1,n2  
      h1=x(i+1)-x(i) 
      h2=x(i+2)-x(i+1) 
      a(i)=h1*h1*h2  
      c(i)=h1*h2*h2  
      b(i)=2.*(a(i)+c(i))  
      d(i+1)=6.*(h1*(y(i+2)-y(i+1))-h2*(y(i+1)-y(i)))  
1     continue 
c        solve tridiagonal matrix system 
      c6=1./6. 
      b(1)=1./b(1) 
      c(1)=c(1)*b(1) 
      d(2)=d(2)*b(1) 
      if(n3.eq.1)go to 6 
      if(n3.eq.0)go to 5 
      do 2 i=2,n3  
      b(i)=1./(b(i)-c(i-1)*a(i-1)) 
      c(i)=b(i)*c(i) 
      d(i+1)=b(i)*(d(i+1)-a(i-1)*d(i)) 
2     continue 
6     continue 
      b(n2)=b(n2)-c(n3)*a(n3)  
      d(n2+1)=(d(n2+1)-a(n3)*d(n2))/b(n2)  
      do 3 k=1,n3  
      i=n2-k 
      d(i+1)=d(i+1)-c(i)*d(i+2)  
3     continue 
5     continue 
      d(1)=0.  
      d(n)=0.  
c       spline coeficients are in d  
c       computation of cubic spline polynomial coefficients  
      do 4 i=1,n1  
      hi=x(i+1)-x(i) 
      hib=1./hi  
      b(i)=.5*d(i) 
      a(i)=(d(i+1)-d(i))*c6*hib  
      c(i)=hib*(y(i+1)-y(i))-c6*(hi*(2.*d(i)+d(i+1)))  
4     continue 
      eps=(x(2)-x(1))*.0001  
      j=1  
      do 11 i=1,m  
13     continue  
      del1=xx(i)-x(j)  
      del2=x(j+1)-xx(i)  
      if(del1+eps.ge.0..and.del2+eps.ge.0.)go to 12  
      if(del1.lt.0..and.j.eq.1)go to 12  
      if(j.eq.n1)go to 12  
      j=j+1  
      go to 13 
12     continue  
      f(i)=y(j)+del1*(c(j)+del1*(b(j)+del1*a(j)))  
11     continue  
14     continue  
15    continue 
      return 
      end
      subroutine trans(b,a,nr,nc) 
      dimension a(*),b(*) 
      jj=1
      do 1 i=1,nc 
       ii=i 
      do 1 j=1,nr 
      b(ii)=a(jj) 
      ii=ii+nc
1     jj=jj+1 
      return
      end
      subroutine trnsp(if1,iw1,if2,iw2,nr,nc,a,b,lbuf,sc)  
      dimension a(*),b(*)  
      real if1(*),if2(*)
c     call g(if1,iw1,nr*nc,a)  
c     call trans(b,a,nr,nc)  
c     call p(if2,iw2,nr*nc,b)  
      call trans(if2(iw2),if1(iw1),nr,nc)
      return 
      end
