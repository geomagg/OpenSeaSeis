c-----------------------------------------------------------------------
c     Rotinas copiadas sem alteracao de dvssyn.f (SEISPAK), exceto onde
c     marcado com "SEASEIS:". Usadas pelo modulo SPK_DVSMIG3D (XVEVENTS, leitura da lista VXT).
c-----------------------------------------------------------------------
      subroutine dwnv2(a,b,e,r,kx2,nx,nl,nh,dw,shft,vv,sgn)
      real kx2(*),r(*)
      complex a(*),b(*),e(*)
      nws=nh-nl+1 
      dkx=2.*3.141593/nx
c     ***  compute kx**2
      nx2=nx/2+1
      kx2(1)=0. 
      do 10 i=2,nx2 
10    kx2(i)=((i-1)*dkx)**2 
      vvrw2=0.
      cos2=.1
      cos1=sqrt(cos2)
      xx1=1./cos2
      xx2=cos1/cos2
      do 20 n=1,nws 
      w=(nl+n-2)*dw 
      if(w.eq.0.)then 
      do 30 i=1,nx
30    b(i)=0. 
      else
      vvrw2=vv/w**2 
      t=sgn*w*shft
      do 40 i=1,nx2 
40    r(i)=1.-vvrw2*kx2(i)
      do 50 i=1,nx2 
c     if(r(i).lt.cos2)then
c     e(i)=cexp(cmplx(0.,t*(sqrt(cos2))))
c     else
c     e(i)=cexp(cmplx(0.,t*(sqrt(r(i)))))
c     end if
      if(r(i) .lt. cos2)then
      e(i)=cexp(cmplx(0.,t*cos1))
      else
      e(i)=cexp(cmplx(0.,t*sqrt(r(i))))
      end if
50    continue
      ii=(n-1)*nx 
      iii=nx+ii+2
      b(ii+1)=e(1)*a(ii+1)
      b(ii+nx2)=e(nx2)*a(ii+nx2)
      do 60 i=2,nx/2 
60    b(ii+i)=e(i)*a(ii+i)
      do 61 i=2,nx/2
61    b(iii-i)=e(i)*a(iii-i)
      end if
20    continue
      return
      end
      subroutine xdwnn(b,nx,ny,dx,dy,w,shft,vv,sgn)
*     b has length 2*nx*ny
      dimension b(*)
      if (w .eq. 0.) then
         do 1 i  =  1,2*nx*ny
            b(i)  =  0.
    1    continue
         return
      end if
      dkx = 2.*3.141593/nx
      dky = 2.*3.141593/ny
      dkx = dkx*dy/dx
      t = sgn*w*shft
      vvrw2 = vv/w**2 
         k = 1
         ind1=3
         ind2=2*nx-1
         do 21 i = 2,nx/2+1
            akx=(i-1)*dkx
            xk2 = akx**2
            r = 1.-vvrw2*xk2
c   set taper region
            r2=.04
            r1=.01
            if (r .lt.r1)then 
               b(ind1) = 0.
               b(ind1+1) = 0.
               b(ind2) = 0.
               b(ind2+1) = 0.
            else
            if (r1.le.r.and.r.le.r2)scl= (r-r1)/(r2-r1)
            if (r2.le.r)scl=1.
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind1)-ei*b(ind1+1)
               b(ind1+1) = er*b(ind1+1)+ei*b(ind1)
               b(ind1) = t1
               t1  = er*b(ind2)-ei*b(ind2+1)
               b(ind2+1) = er*b(ind2+1)+ei*b(ind2)
               b(ind2)  = t1
                b(ind1)  =scl*b(ind1)
                b(ind1+1)=scl*b(ind1+1)
                b(ind2+1)=scl*b(ind2+1)
                b(ind2)  =scl*b(ind2)
            end if
            ind1 = ind1 + 2
            ind2 = ind2 - 2
   21    continue
         if(ny .eq. 1)return
         k = ny/2+1
         ind = 1 + (k-1)*2*nx
         aky = min0(k-1,ny-k+1)*dky
         aky2 = aky**2
         do 22 i = 1,nx
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .lt. .01) then 
               b(ind) = 0.
               b(ind+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
            end if
            ind = ind + 2
   22    continue
         do 24 k = 2,ny/2
            kk = ny-k+2
            aky = (k-1)*dky
            aky2 = aky**2
            i = 1
            ind = i + (k-1)*2*nx
            indd = i + (kk-1)*2*nx
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .lt. .01) then 
               b(ind) = 0.
               b(ind+1) = 0.
               b(indd) = 0.
               b(indd+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
               t2 = er*b(indd)-ei*b(indd+1)
               b(indd+1) = er*b(indd+1)+ei*b(indd)
               b(indd) = t2
            end if
            i = nx/2+1
            ind = 2*i-1 + (k-1)*2*nx
            indd = 2*i-1 + (kk-1)*2*nx
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .lt. .01) then 
               b(ind) = 0.
               b(ind+1) = 0.
               b(indd) = 0.
               b(indd+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
               t2 = er*b(indd)-ei*b(indd+1)
               b(indd+1) = er*b(indd+1)+ei*b(indd)
               b(indd) = t2
            end if
            i = 2
            ind = 2*i-1 + (k-1)*2*nx
            indd = 2*i-1 + (kk-1)*2*nx
            inb = 2*(nx-i+2)-1 + (k-1)*2*nx
            inbb = 2*(nx-i+2)-1 + (kk-1)*2*nx
CDIR$ IVDEP
         do 43 i = 2,nx/2
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .lt. .01) then 
               b(ind) = 0.
               b(ind+1) = 0.
               b(indd) = 0.
               b(indd+1) = 0.
               b(inb) = 0.
               b(inb+1) = 0.
               b(inbb) = 0.
               b(inbb+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
               t2 = er*b(indd)-ei*b(indd+1)
               b(indd+1) = er*b(indd+1)+ei*b(indd)
               b(indd) = t2
               t3  = er*b(inb)-ei*b(inb+1)
               b(inb+1) = er*b(inb+1)+ei*b(inb)
               b(inb) = t3
               t3 = er*b(inbb)-ei*b(inbb+1)
               b(inbb+1) = er*b(inbb+1)+ei*b(inbb)
               b(inbb) = t3
            end if
            ind = ind + 2
            indd = indd + 2
            inb = inb - 2
            inbb = inbb - 2
   43    continue
   24    continue
      return
      end
      subroutine xint2(f,x,y,x1,dx,n,m) 
c   interpolate m points (x(i),y(i))  for x = x1 to x1+(n-1)*dx  
c   x1 <= x(1) < x(i) < x(i+1) must be satisfied.  
c   f(x)=0 if x < x1 or x > x1+(n-1)*dx  
c   interpolated results in f  
      double precision d1,d2,a,r1  
      dimension f(*),x(*),y(*) 
      j=1  
      xf=x1-dx 
      a=0. 
      eps=.0001  
      if(x(2)-x(1).gt.eps)a=(y(2)-y(1))/(x(2)-x(1))  
      do 30 i=1,n  
      xf=xf+dx 
13    d1=xf-x(j) 
      d2=x(j+1)-xf 
      if(d1.lt.0.)then 
      f(i)=0.  
      goto 29  
      end if 
      if(d2.lt.0.)then 
      j=j+1  
      if(j.eq.m)goto 40  
      r1=x(j+1)-x(j) 
      a=y(j+1)-y(j)  
      if(r1.gt.eps)then  
      a=a/r1 
      else 
      a=0. 
      end if 
      goto 13  
      else 
      e=d1*a 
      f(i)=y(j)+e  
      end if 
29    continue 
30    continue 
40    do 50 j=i,n  
50    f(j)=0.  
      return 
      end
      subroutine xinterpv(e,v,x,t,vv,xx,tt,npnts,nx,mode) 
      dimension e(*),x(*),v(*),t(*),vv(*)  
      common/rdd/rdd 
      dimension xx(*),tt(*)  
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
      print 56,vnm,(v(l),l=il,ir)  
      print 56,xnm,(x(l),l=il,ir)  
      print 56,tnm,(t(l),l=il,ir)  
      il=il+15 
        ir=min0(ir+15,npnts) 
55    continue 
56    format(2x,a4,15f8.0) 
      l1=nx+1  
        l2=l1+npnts  
        l3=l2+npnts  
        l4=l3+npnts  
      go to 15 
      jj=1 
      do 5 i=1,npnts1  
      xl=x(i)  
      xr=x(i+1)  
      ixl=xl+.1  
      ixlm=ixl-mindpn+1  
      nxpnts=xr-xl+1.1 
      nxpnts1=nxpnts-1 
      den=1./(xr-xl) 
      do 5 j=1,nxpnts1 
      ix=ixlm+j-1  
      cr=(j-1)*den 
      cl=(xr-xl+1.-j)*den  
      vv(ix)=cl*v(i)+cr*(v(i+1)) 
      tt(jj)=cl*t(i)+cr*t(i+1) 
      jj=jj+1  
5     continue 
      vv(ixlm+nxpnts-1)=v(npnts) 
      tt(jj)=t(npnts)  
      go to 99 
15    continue 
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
      subroutine xmultpspi(y,a,b,s,rs,nl,nh,nt1,nx,shft,vv, 
     *   sgn,nvels)
      complex y(*),a(*),b(*),s(*) 
      real shft(*),vv(*),rs(*)
      nws=nh-nl+1 
      dw=2.*3.1415927/nt1 
      dx=1.
      f1=1./nx
      f2=-1.
c   ***   find min and max shifts 
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=-1000000.
      do 5 i=1,nx
      v1=amax1(v1,vv(i))
      v2=amin1(v2,vv(i))
      shft1=amin1(shft1,shft(i))
      shft2=amax1(shft2,shft(i))
5     continue
      vmax=1./sqrt(v2)
      vmin=1./sqrt(v1)
      shftmax=shft2
      shftmin=shft1
      dzx=shft1/vmin
c   ***   fft, downward continue, and inverse fft for shift=vmin
      call fftc(y,s,2,2*nx,nx,nws,-1)
c   ***   copy kx-w field to scratch array 
      if(vmax-vmin .ge. .01*vmax)call movlev(y,a,2*nx*nws)
c   ***   downward continue by phase-shift for shft1,v1
c   ***   if nvels = 1 the v1 = ( .5*(1./vmin+1./vmax) )**2
c   ***                shft1 = .5*(shftmin+shftmax)
      if(nvels .eq. 1)then
      v1=(.5*(1./vmin + 1./vmax) )**2
      shft1=.5*(shftmin+shftmax)
      end if
      do iw=nl,nh
      w=(iw-1)*dw
      ii=(iw-nl)*nx+1
      call xdwnn(y(ii),nx,1,dx,dx,w,shft1,v1,sgn)
      enddo
c   ***   inverse fft data and scale
      call fftc(y,s,2,2*nx,nx,nws,+1)
      do j=1,nx*nws
      y(j)=f1*y(j)
      end do
      if(vmax-vmin .lt. .01*vmax .or. nvels .eq. 1)then
c     call vshift1(y,a,a,b,shft,dw,nl,nh,nx,sgn)
      call sgivshift1(y,shft,dw,nl,nh,nx,sgn)
      return
      end if
c   ***   velocity non-constant in this interval
c   ***   apply multi-velocity pspi
c   ***   make vmax slightly  larger than max of 1./sqrt(vv)
c   ***   make vmin slightly smaller than min of 1./sqrt(vv)
        vmax=1.001*vmax
        vmin=.9999*vmin
      nsteps=(vmax-vmin)/(.2*vmax) +1
      nsteps=nvels-1
      vvinc= (vmax-vmin)/nsteps
      shftinc=(shftmax-shftmin)/nsteps
      call movlev(y,b,2*nx*nws)
      call clear(y,2*nx*nws)
      do istep=1,nsteps
      vv1=vmin+(istep-1)*vvinc
      vv2=vv1+vvinc
      v1 = (1./vv1)**2
      v2 = (1./vv2)**2
      shft1=shftmin+(istep-1)*shftinc
      shft2=shft1+shftinc
c   ***   compute contribution of field downward continued using vv1
      do k=1,nx
      rs(k) = ( vv2-1./sqrt(vv(k)) )/(vv2-vv1)
      rs(k) = amax1(0.,rs(k))
      if(rs(k) .gt. 1.)rs(k)=0.
      end do
c   ***   add contrubution of field downward continued using vv1
c   ***
      do j=1,nws
      i1=(j-1)*nx
      do k=1,nx
      y(i1+k) = y(i1+k) + rs(k)*b(i1+k)
      end do
      end do
c   ***
      call movlev(a,b,2*nx*nws)
      do iw=nl,nh
      w=(iw-1)*dw
      ii=(iw-nl)*nx+1
      call xdwnn(b(ii),nx,1,dx,dx,w,shft2,v2,sgn)
      enddo
c   ***   inverse fft data and scale
      call fftc(b,s,2,2*nx,nx,nws,+1)
      do j=1,nx*nws
      b(j)=f1*b(j)
      end do
c   ***   compute contribution of field downward continued using vv2
c   ***
      do k=1,nx
      rs(k)=(1./sqrt(vv(k))-vv1)/(vv2-vv1)
      rs(k)=amax1(0.,rs(k))
      if(rs(k) .ge. 1.)rs(k)=0.
      end do
c   ***   add contribution of field downward continued using vv2
      do j=1,nws
      i1=(j-1)*nx
      do k=1,nx
      y(i1+k) = y(i1+k) + rs(k)*b(i1+k)
      end do
      end do
c   ***
       end do
c   ***
c   ***   apply thin lens correction
c   ***
c     call vshift1(y,a,a,b,shft,dw,nl,nh,nx,sgn)
      call sgivshift1(y,shft,dw,nl,nh,nx,sgn)
      return
      end
      subroutine xvevents(ee,e,xx,x,nx,nny,nt,vint,nwds,ia,ib,ddd) 
      dimension ee(*),e(*),xx(*),x(*),nevents(100),ntodo(100)  
      dimension mode(100)  
      real ia(*),ib(*)
      common/rdd/rdd 
      save
      character*4 topdwn,vint,depth,rms
      data topdwn/'yes'/,depth/'no'/,rms/'no'/
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
      write(*,*)'rdd  e(1) ',rdd,e(1)
      e(1)=e(1)-4.
      end if
      iloc=2 
      if(e(1).eq.3.5 .or. e(1).eq.4.5)rms='yes'
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
      if(icheck.ge.2)go to 1111 
      k=1 
      dpnmax=0. 
      dpnmin=1000000. 
      vmax=0. 
      vmin=1000000. 
      anx=0.
      x1=-1000.
      j=2 
      jj=iloc
      ii=1
  80  continue
      vv=e(jj)
      x2=e(jj+1)
      tt=e(jj+2)
      ii=ii+1
      if(rdd.ne.1.)x2=aint(x2*rdd+1.5)
      anx=amax1(anx,x2) 
      dpnmax=amax1(dpnmax,x2) 
      dpnmin=amin1(dpnmin,x2) 
      vmax=amax1(vmax,vv) 
      vmin=amin1(vmin,vv) 
      if(1.+x1.le.x2)go to 1100 
      print 1106,k,ii
1106  format('dpns not increasing in horizon ',i3,' triple ',i3,
     *' ejected')
      jj=jj+3
      if(e(jj).lt..5)goto 1109
      goto 80
1100    continue
      ee(j)=vv
      ee(j+1)=x2
      ee(j+2)=tt
      j=j+3 
      jj=jj+3
      x1=x2 
      if(e(jj).gt..5)go to 80 
1109  continue
      ee(j)=e(jj)
      j=j+1 
      jj=jj+1
      k=k+1 
      x1=0. 
      ii=0
      if(e(jj).gt..5)go to 80 
      ee(j)=e(jj)
      k=k-1 
      write(*,*)' number of events found = ',k 
      write(*,*)' minimum dpn found = ',dpnmin 
      write(*,*)' maximum dpn found = ',dpnmax 
      write(*,*)' minimum velocity found = ',vmin  
      write(*,*)' maximum velocity found = ',vmax  
      if(vint.eq.'no')write(*,*)'velocity constant between horizons' 
      if(vint.eq.'yes')write(*,*)'velocity varies between horizons'  
      nx=dpnmax-dpnmin+1.1 
      write(*,*)' j  jj ',j,jj
      ee(j+1)=0
      ee(j+2)=0
      call movlev(ee(2),e(iloc),j+1)
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
      fy=ee(iy1)  
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
      if(mode(j).eq.0)write(*,1101)j 
      if(mode(j).eq.-1)write(*,1102)j  
      if(mode(j).le.-2)write(*,1103)j  
      if(mode(j).eq.-3)write(*,1104)j  
1104  format(' interpolation of velocities for event ',i3,' is linear')  
1101  format('  interpolation of event',i3,' by linear interpolation') 
1102  format('  interpolation of event',i3,' by spline interpolation') 
1103  format('  interpolation of event',i3,' by joeh   interpolation') 
      call xinterpv(x(kk),x(l1),x(l2),x(l3),x(l4),xx,x(l5), 
     *   np,nx,mode(j))
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
c     call trnsp(ia,ifw,ib,ifwb,nx,2*nevents(igrp),x,ee,10000,1.)  
      call trans(ib(ifwb),ia(ifw),nx,2*nevents(igrp)) 
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
      if(rms.eq.'yes')then
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
      subroutine xxconv(y,a,b,s,rs,nl,nh,nt1,nx,shft,vv,dx,dt,dz,sgn) 
      complex y(*),a(*),b(*),s(*) 
      complex wvtab(17,2001)
      real shft(*),vv(*),rs(*)
      integer iwv(5000)
      data kk/0/
      save
      kk=kk+1
      nws=nh-nl+1 
      dw=2.*3.1415927/nt1 
c   ***   set up downward continuation convolution table at first call
c   ***   convolution components for this shft and vv stored in wvtab
      if(kk.eq.1)then
c   ***   wvtab is l by ll
c      call etime(rs,xr)
c      t1=rs(1)
c      t2=rs(2)
      l=13
      nxsave=nx
      nx=128
      ll=2001
c     call c1dfft(a,nx,s,-3,ier)
      ddw=nh*(dw/ll)
      xshft=shft(1)
      xvv=vv(1)
      write(*,*)'shft and vv at table initiation ',xshft,xvv
      do i=1,ll
      do j=1,nx
      a(j)=cmplx(1.,0.)
      end do
      ww=i*ddw
      call dwnv2(a,a,b,rs,rs(nx+1),nx,2,2,ww,xshft,xvv,sgn)
c     call c1dfft(a,nx,s,-1,ier)
      call fftc(a,s,2,2*nx,nx,1,+1)
      do jj=1,nx
      a(jj)=(1./nx)*a(jj)
      end do
      nx0=nx/2
      nx0=min0(nx0,2*l)
      do k=l,nx0-1
      a(l)=a(l)+a(k+1)
      end do
      do j=1,l
      wvtab(j,i)=a(j)
      end do
      end do
c     call etime(rs,xr)
c     t1=rs(1)-t1
c     t2=rs(2)-t2
c     write(*,*)'time to compute table',t1,t2
      nx=nxsave
      return
      end if
c   ***   downward continue by convolution with pre-computed operator
c   ***   load convolution operator for this w/v
      wvmax=nh*dw/sqrt(xvv)
      do j=1,nx
      rs(j)=1./sqrt(vv(j))
      end do
      do 104 m=1,nws
      xx=ll*(m+nl-2)*dw/wvmax
      do j=1,nx
      iwv(j)=.5+xx*rs(j)
      end do
      ii=(m-1)*nx+1
c   ***
c   ***   convolution
c   ***
      l1=l-1
      iii=ii+nx-l
      do  j=1,l1 
      s(j)=y(iii+j)
      end do
      iii=ii-1
      do j=1,l1
      s(nx+l1+j)=y(iii+j) 
      end do
      do j=1,nx
      s(l1+j)=y(iii+j)
      end do
      do  j=1,nx
      y(iii+j)=wvtab(1,iwv(j))*s(l1+j)
      end do
      do i=1,l1
      do j=1,nx
      y(iii+j)=y(iii+j)+wvtab(i+1,iwv(j))*(s(l1-i+j)+s(l1+i+j))
      end do
      end do
c   ***
c   ***   end convolution
c   ***
104   continue
      return
      end
