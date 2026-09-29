c-----------------------------------------------------------------------
c     Rotinas copiadas sem alteracao de dvslib.f (SEISPAK).
c     Usadas pelo modulo SPK_XVSPMIG2.
c-----------------------------------------------------------------------
      subroutine cacomt2(v,ifitvv,s) 
      real s(*),v(*),ifitvv(*)
      character*4 vint,rfc
      common/velneed/nrecs,irec,lev,vinc,ipass,npass,nxv,lag,locv,ispj
      common/migneed/nx,nt,nz,itzr,vint,dt,nnz,ddz,beta,dx,apx,maxscr, 
     *rfc,alph1,gam1,frq,rfstir,itype,cyc,raa,padl,padr,zzy,ltm,ny,dy,  
     *frq1,frq2,alph2,rni,alph3,
     *vmin,vmax,seps,lev1,lev2,ifor,ilor,window,lunin,lunout
      ritzr=1./itzr
      nz=(nnz-1)/itzr
      lltm=ltm/itzr 
      nz1=nz+lltm+2 
      dz=itzr*ddz
      dt2=.5*dt
c     call gh(ifitv,1,5,v)
c     call movlev(ifitv(1),v,5) 
      nxvv=v(1) 
      mindpn=v(5) 
      ntv=v(3)
      dpnv1=mindpn
      vnx=nxvv
      idpn=dpnv1
      ll=ll+1 
      lag=gam1+1
c   ***   velocity location will be offset by oty=beta dpns 
      ispj=zzy
      if(lev.lt.1 .or. lev .gt. ntv)lev=ntv
      if(ispj.eq.0)ispj=1 
c   ***   set velocity location for first record
      vvdpn=rfstir+beta 
      iwa=1 
      nxv=nx+(nrecs-1)*iabs(ispj) 
      if(itype.eq.2)nxv=nx
      do 100 i=1,nxv
      vdpn=vvdpn+i-1. 
      if(cyc.eq.2.)vdpn=vdpn-2. 
      vdpn=amax1(vdpn,dpnv1)
      vdpn=amin1(vdpn,dpnv1+vnx-1.) 
      ix=vdpn-dpnv1+1.
      ivwa=4+(2*ntv+3)*(ix-1) 
c     call gh(ifitv,ivwa,2*ntv+3,v) 
c     call movlev(ifitv(ivwa),v,2*ntv+3)
      ivv=2*ntv+4
c     ii=3
      ii=ivwa+2
111   continue

      if(v(ii).eq.0. )then
      ii=ii+2 
      goto 111
      end if
      v1=v(ii)
      v2=v1 
      z=0.
      z2=v(ii+1)
      eps=.1
      z1=0. 
      dzr=1./(z2-z1+eps)
      do 101 n=1,nz1
      rvbar=0.
      do 104 l=1,itzr 
102   continue
      if(ii .gt. ivwa+2*lev .or. z .le. z2)goto 103
      if(v(ii).eq.0. .or. z .gt. v(ii+1))then 
      ii=ii+2 
      goto 102
      end if
      v1=v2 
      v2=v(ii)
      z1=z2 
      z2=v(ii+1)
      dzr=1./(z2-z1+eps)
103   continue
      vel=v2
      if(vint.eq.'yes'.and.ii.lt.2*lev+2) 
     *vel=((z-z1)*v2+(z2-z)*v1)*dzr 
      if(rfc.eq.'YES')then
      z=z+ddz 
      else
      z=z+.5*vel*dt
      end if
      rvbar=rvbar+1./vel
c   *** 
104   continue
      s(n)=rvbar*ritzr
101   continue
c   ***   filter velocity reciprocals
c     call filtv(s,s(nz1+1),nz1)
c   ***   store velocity reciprocals
      do 201 j=1,nz1
      ifitvv(i+(j-1)*nxv)=s(j) 
201   continue
100   continue
c   ***   apply lateral filter to velocity reciprocals
      do 300 j=1,nz1
      ii=(j-1)*nxv+1
c     call filtv(ifitvv(ii),s,nxv)
      do 300 k=0,nxv-1
      ifitvv(ii+k)=1./ifitvv(ii+k)
300   continue
      return
      end 
      subroutine arriv6(t1,t2,vv,shft,s,s1,s2,t,jltm,lltm,nxx,xk)
      real vv(*),shft(*) 
      real s(*),s1(*),s2(*),t(*),t1(*),t2(*) 
      character*4 vint,rfc
      common/migneed/nx,nt,nz,itzr,vint,dt,nnz,ddz,beta,dx,apx,maxscr, 
     *rfc,alph1,gam1,frq,rfstir,itype,cyc,raa,padl,padr,zzy,ltm,ny,dy,  
     *frq1,frq2,alph2,rni,alph3,
     *vmin,vmax,seps,lev1,lev2,ifor,ilor,window,lunin,lunout
      save
      k=xk
      k=max0(1,k)
      if(k.eq.1)then 
      dz=itzr*ddz  
      xdz=(xk-k)*dz
      vmax=0.  
      vmin=100000. 
      lltm=0 
      jltm=0 
      shft1=0  
      iltm=0 
      do 9 i=1,nxx 
      t1(i)=0. 
9     t2(i)=0. 
      end if 
      do 10 i=1,nxx  
      vmin=amin1(vmin,shft(i)) 
10    vmax=amax1(vmax,shft(i)) 
      if(vmax-vmin.lt..02*vmax)iltm=iltm+1 
      if(iltm.lt.k)jltm=jltm+1 
      if(jltm.eq.1)lltm=iltm 
      if(itype.eq.2)return
      if(jltm.eq.0)then  
c ***   compute arrival times for z=(iltm-1)*dz+xdz and z=iltm*dz+xdz  
      al1=alph1  
      if(cyc.eq.2.)al1=alph1+2 
      ial1=al1 
      c1=((iltm-1)*dz+xdz)**2  
      c2=(iltm*dz+xdz)**2  
      c3=shft(1)/dz  
      do 63 l=1,nxx  
      s(l)=(l-al1)*dx  
63    continue 
c   ***   if raa=1 then start plane wave
      if(raa.eq.1.)then
      do l=1,nxx
      s(l)=0.
      end do
      end if
      do 20 i=1,nxx  
      s(i)=s(i)*s(i) 
      t1(i)=c3*sqrt(c1+s(i)) 
20    t2(i)=c3*sqrt(c2+s(i)) 
      return 
      end if 
      if(jltm.gt.0)then  
c   ***  
c   ***   compute imaging times by downward continuation 
c   ***  
      if(k.eq.1.and.jltm.eq.1)then 
c   ***   compute arrival times for z=0. and z=dz. 
      al1=alph1  
      if(cyc.eq.2.)al1=alph1+2 
      ial1=al1 
      c1=xdz**2
      c2=(dz+xdz)**2  
      c3=shft(ial1)/dz  
      do 633 l=1,nxx  
      s(l)=(l-al1)*dx  
633    continue 
c   ***   if raa=1 then start plane wave
      if(raa.eq.1.)then
      do l=1,nxx
      s(l)=0.
      end do
      end if
      do 634 i=1,nxx  
      s(i)=s(i)*s(i) 
      t1(i)=c3*sqrt(c1+s(i)) 
634    t2(i)=c3*sqrt(c2+s(i)) 
      return 
      end if 
c   ***  
c   ***   downward continue arrival times in samples by ray-trace. 
c   ***  
      do 71 i=1,nxx  
      t1(i)=t2(i)  
71    s2(i)=t2(i)  
c   ***  
      ritzr=.5/itzr  
      do 501 j=1,2*itzr  
      do 502 i=1,nxx 
502   t(i)=s2(i) 
c   ***   compute arrival times for i=1 and i=nxx  
      sn=(t(2)-t(1))*sqrt(vv(1)) 
      r=amax1(1.-sn**2,0.)
      cn=sqrt(r)
      s2(1)=t(1)+cn*ritzr*shft(1)
      sn=(t(nxx)-t(nxx-1))*sqrt(vv(nxx)) 
      r=amax1(1.-sn**2,0.)
      cn=sqrt(r)
      s2(nxx)=t(nxx)+cn*ritzr*shft(nxx)  
      do 500 i=2,nxx-1 
      dt1=t(i)-t(i-1)  
      dt2=t(i+1)-t(i)  
      dt3=.5*(dt1+dt2) 
      if(dt1.ge.0.)delt=dt1  
      if(dt2.le.0.)delt=dt2  
      if(dt1*dt2.lt.0.)delt=dt3  
      sn=sqrt(vv(i))*delt  
      r=1.-sn**2 
      if(r.ge.0.)then
      s2(i)=t(i)+sqrt(r)*ritzr*shft(i)
      else
      if(dt2.le.0.)s2(i)=t(i+1)+1./sqrt(vv(i))
      if(dt1.ge.0.)s2(i)=t(i-1)+1./sqrt(vv(i))
      if(dt2.gt.0. .and. dt1.lt.0.)
     * s2(i)=(1./3.)*(t(i-1)+t(i)+t(i+1))-sqrt(-r)*ritzr*shft(i)
      end if
500   continue 
501   continue 
      do 75 j=1,nxx
      t2(j)=s2(j)
75    continue
      end if
      return
      end 
      subroutine multpspi(y,a,b,s,rs,nl,nh,nt1,nx,shft,vv, 
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
      call dwnn(y(ii),nx,1,dx,dx,w,shft1,v1,sgn)
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
        vmax=1.0001*vmax
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
      call dwnn(b(ii),nx,1,dx,dx,w,shft2,v2,sgn)
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
      subroutine dwnn(b,nx,ny,dx,dy,w,shft,vv,sgn)
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
            if (r .le. 0.) then 
               b(ind1) = 0.
               b(ind1+1) = 0.
               b(ind2) = 0.
               b(ind2+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind1)-ei*b(ind1+1)
               b(ind1+1) = er*b(ind1+1)+ei*b(ind1)
               b(ind1) = t1
               t1  = er*b(ind2)-ei*b(ind2+1)
               b(ind2+1) = er*b(ind2+1)+ei*b(ind2)
               b(ind2) = t1
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
            if (r .le. 0.) then 
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
            if (r .le. 0.) then 
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
            if (r .le. 0.) then 
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
            if (r .le. 0.) then 
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
      subroutine sgivshift1(y,shft,dw,nl,nh,nx,sgn)
      real y(*)
      real shft(*)
c   ***   efficient  complex exponential table computation
      w1=sgn*(nl-1)*dw
      dws=sgn*dw
      nws=nh-nl+1
      nx2=2*nx
      do 10 i=1,nx
c     s1=cexp(cmplx(0.,dws*shft(i)))
c     a=cexp(cmplx(0.,w1*shft(i)))
      sr=cos(dws*shft(i))
      si=sin(dws*shft(i))
      ar=cos( w1*shft(i))
      ai=sin( w1*shft(i))
      ii=2*(i-1)+1
      do 20 j=1,nws
c     y(ii+i)=a*y(ii+i)
c     a=s1*a
      t1=ar*y(ii)-ai*y(ii+1)
      y(ii+1)=ar*y(ii+1)+ai*y(ii)
      y(ii)=t1
      t1=ar*sr-ai*si
      ai=ar*si+ai*sr
      ar=t1
      ii=ii+nx2
20    continue
10    continue
      return
      end
      SUBROUTINE FFTC(DATA,WORK,INC,JUMP,LEN,LOT,IDIR)
C
C ABSTRACT *************************************************************
C
C     COMPLEX TO COMPLEX, PERIODIC, MULTIPLE FFTS, POWERS OF 2,3,5
C     DATA HAS REAL AND IMAGINARY PARTS IN SINGLE ARRAY, WITH EITHER
C     REAL/IMAG PARTS IN EVERY OTHER ROW, OR EVERY OTHER COLUMN
C
C KEYWORDS
C  FFT
C
C PURPOSE
C
C     COMPLEX TO COMPLEX, PERIODIC, MULTIPLE FFTS, POWERS OF 2,3,5
C     DATA HAS REAL AND IMAGINARY PARTS IN SINGLE ARRAY, WITH EITHER
C     REAL/IMAG PARTS IN EVERY OTHER ROW, OR EVERY OTHER COLUMN
C
C     DATA       REAL ARRAY      IN/OUT      DATA ARRAY OF LENGTH
C                                            2*LEN*LOT (REAL)
C     WORK       REAL ARRAY      SCRATCH     WORK ARRAY OF LENGTH
C         LENGTH SHOULD BE FOUND WITH ROUTINE FFTCL
C     INC        INTEGER         IN    INCREMENT BETWEEN SUCCESSIVE
C                                      REAL ELEMENTS OF DATA
C     JUMP       INTEGER         IN    INCREMENT BETWEEN FIRST ELEMENT
C                                      OF SUCCESSIVE VECTORS
C     LEN        INTEGER         IN    LENGTH OF EACH FFT ( PRIME
C                                      FACTORS MUST BE 2, 3, 5 )
C     LOT        INTEGER         IN    NUMBER OF TRANSFORMS
C     IDIR       INTEGER         IN    SIGN IN EXPONENT; NORMALIZATION
C                                      IS 1 INDEPENDENT OF IDIR
C
C  X(K) = SUM(N=0,LEN-1) X(N)EXP(2*PI*I*IDIR*K*N/LEN)
C
C
C
C  NORMAL USE IS EITHER
C     1)  FFT DOWN COLUMNS (ROW 1,3,... ARE REAL PARTS; 2,4,... ARE IMAG)
C         JUMP ROWS * LOT COLUMNS, INC = 2, JUMP = 2*LEN
C  OR
C     2)  FFT ACROSS ROWS (COL 1,3,... ARE REAL PARTS; 2,4,... ARE IMAG)
C         INC ROWS * LEN COLUMNS, JUMP = 2, INC = 2*LOT
C
C  EXAMPLES
C     1)  INPUT IS A REAL MATRIX WITH 2*LEN ROWS (ROW 1 REAL, ROW 2 IMAG)
C         AND LOT COLUMNS  -  FFT DOWN EACH OF LOT COLUMNS
C         CALL FFTC(DATA,WORK,INC=2,JUMP=2*LEN,LEN,LOT,-1)
C     2)  INPUT IS A MATRIX WITH LOT ROWS AND 2*LEN COLUMNS (COL 1 REAL,
C         COL 2 IMAG )
C         CALL FFTC(DATA,WORK,INC=2*LOT,JUMP=2,LEN,LOT,-1)
C
C     CONVERSION FOR CRAY FFTS IS: (REPLACE C999 WITH FFT)
C     C999(DATAR,DATAI,WORK,TRIGS,IFAX,INC=2,JUMP=2*LEN,LEN,LOT,IDIR)
C     FFTC(DATAR,DATAI,WORK,INC=2,JUMP=2*LEN,LEN,LOT,IDIR)
C
C HISTORY
C
C   MP01     ,  09/20/ 89,  M PORTNEY,  INITIAL RELEASE.
C   SS02     ,  29/SEP/89,  SL SMITH ,  GENERIC VERSION
C   MP03     ,  10/NOV/90,  M PORTNEY   MAKE VECTOR LENGTH PARAMETER STMNT,
C                                       USE FFTCL TO FIND WORK ARRAY,
C                                       REMOVE SEPARATE ARRAY FOR TRIG FACTS
C
* END ******************************************************************
*     PORTABLE VERSION
*     LENV MUST MATCH IN FFTC1L,FFTC1,FFTCL,FFTC,VC999
      PARAMETER (lenv = 8)
      DIMENSION DATA(*),WORK(*),IFAX(13)
      IT = 1 + 2*LEN
      CALL VFTFAX(LEN,IFAX,WORK)
      CALL VC999(DATA,WORK(IT),WORK,IFAX,INC,JUMP,LEN,LOT,IDIR)
C     CALL GENERC('FFTC','WORKING')
      RETURN
      END
      SUBROUTINE FACTR ( N, IP, NP)
      INTEGER ISCR(300), IP(6), NP(6), IPP(300), NPP(300)
      NMAX = SQRT(REAL(N))
      DO 1 I = 1,6
            IP(I) = 0     ! PRIME FACTORS
            NP(I) = 0     ! NUMBER OF OCCURRENCES
    1 CONTINUE
      CALL SIEVE ( NMAX, ISCR, IPP, NPR )
      DO 2 I = 1,NPR
            NPP(I) = 0
    2 CONTINUE
      NREDUC = N
      LFACT = 0
    5 CONTINUE
      DO 20 I = 1,NPR
            IMOD = MOD(NREDUC,IPP(I))
            IF (IMOD .EQ. 0.) THEN
                 NPP(I) = NPP(I) + 1     ! NUMBER OF OCCURS
                 NREDUC = NREDUC/IPP(I)          ! NEW NUMBER
                 IF (NREDUC .NE. 1) THEN
                      GO TO 5            ! FIND NEXT FACTOR
                 ELSE
                      GO TO 999          ! DONE
                 END IF
            END IF
            IF (I .EQ. NPR) THEN
                LFACT = NREDUC              ! NO MORE FACTORS
            END IF
   20 CONTINUE
  999 CONTINUE
      J = 0
      DO 40 I = 1,NPR
            IF (NPP(I) .NE. 0) THEN
                  J = J + 1          ! NUMBER OF PRIMES
                  NP(J) = NPP(I)     ! NUMBER OF OCCURS
                  IP(J) = IPP(I)          ! PRIMES
            END IF
   40 CONTINUE
      NFACT = J                         ! NUMBER OF PRIMES
      IF (LFACT .NE. 0) THEN
           NFACT = NFACT + 1
           NP(NFACT) = 1
           IP(NFACT) = LFACT
      END IF
      RETURN
      END
      SUBROUTINE FFTFLY(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
*     COMPLEX F(*),T(*),CS(*)
*     COMPLEX COEF,COEF1,COEF2,COEF3
      real f(*),t(*),cs(*),coef,coef1,coef2,coef3,coef4
      IT=1
      LT0=1
      DO 400 JF=1,NW
         JT=IT
         NDX=LT0
         LT=LT0
         IF(IP.EQ.2)THEN
*           COEF=CS(NDX)
            coefr = cs(2*ndx-1)
            coefi = cs(2*ndx  )
            DO 210 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               DO 211 I=1,LOT
*                  T(I1+I)=F(I2+I)+COEF*F(I3+I)
                   t(2*(i1+i)-1) = f(2*(i2+i)-1) + coefr*f(2*(i3+i)-1)
     :                                           - coefi*f(2*(i3+i)  )
                   t(2*(i1+i)  ) = f(2*(i2+i)  ) + coefi*f(2*(i3+i)-1)
     :                                           + coefr*f(2*(i3+i)  )
  211          CONTINUE
  210       CONTINUE
         ELSEIF(IP.EQ.3)THEN
*           COEF=CS(NDX)
*           COEF1=COEF*COEF
            coefr = cs(2*ndx-1)
            coefi = cs(2*ndx  )
            coef1r = coefr*coefr - coefi*coefi
            coef1i = coefr*coefi + coefi*coefr
               DO 220 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
                  DO 221 I=1,LOT
*                    T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)
                    t(2*(i1+i)-1) = f(2*(i2+i)-1) + coefr*f(2*(i3+i)-1)
     :                                            - coefi*f(2*(i3+i)  )
     :                                            + coef1r*f(2*(i4+i)-1)
     :                                            - coef1i*f(2*(i4+i)  )
                    t(2*(i1+i)  ) = f(2*(i2+i)  ) + coefi*f(2*(i3+i)-1)
     :                                            + coefr*f(2*(i3+i)  )
     :                                            + coef1i*f(2*(i4+i)-1)
     :                                            + coef1r*f(2*(i4+i)  )
  221             CONTINUE
  220          CONTINUE
         ELSEIF(IP.EQ.5)THEN
*           COEF=CS(NDX)
*           COEF1=COEF*COEF
*           COEF2=COEF*COEF1
*           COEF3=COEF*COEF2
            coefr = cs(2*ndx-1)
            coefi = cs(2*ndx  )
            coef1r = coefr*coefr - coefi*coefi
            coef1i = coefr*coefi + coefi*coefr
            coef2r = coefr*coef1r - coefi*coef1i
            coef2i = coefr*coef1i + coefi*coef1r
            coef3r = coefr*coef2r - coefi*coef2i
            coef3i = coefr*coef2i + coefi*coef2r
            DO 530 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
               I5 = (K+LF+LF+LF-1)*INC
               I6 = (K+LF+LF+LF+LF-1)*INC
               DO 531 I=1,LOT
*                 T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)+
*    :                    COEF2*F(I5+I)+COEF3*F(I6+I)
                  t(2*(i1+i)-1) = f(2*(i2+i)-1) + coefr*f(2*(i3+i)-1)
     :                                          - coefi*f(2*(i3+i)  )
     :                                          + coef1r*f(2*(i4+i)-1)
     :                                          - coef1i*f(2*(i4+i)  )
     :                                          + coef2r*f(2*(i5+i)-1)
     :                                          - coef2i*f(2*(i5+i)  )
     :                                          + coef3r*f(2*(i6+i)-1)
     :                                          - coef3i*f(2*(i6+i)  )
                  t(2*(i1+i)  ) = f(2*(i2+i)  ) + coefi*f(2*(i3+i)-1)
     :                                          + coefr*f(2*(i3+i)  )
     :                                          + coef1i*f(2*(i4+i)-1)
     :                                          + coef1r*f(2*(i4+i)  )
     :                                          + coef2i*f(2*(i5+i)-1)
     :                                          + coef2r*f(2*(i5+i)  )
     :                                          + coef3i*f(2*(i6+i)-1)
     :                                          + coef3r*f(2*(i6+i)  )
  531          CONTINUE
  530       CONTINUE
         ENDIF
         IT=IT+LFOLD
         IF(IT.GT.N)IT=1
         LT0=LT0+LF
  400 CONTINUE
      RETURN
      END
      SUBROUTINE SIEVE (N, ISCR, IP, NP)
*
* Abstract *************************************************************
*
*   prime number generator
*
* Keywords
*
*   Fourier transform, math.
*
* Purpose
*
*   Finds prime numbers using sieve of Eratosthenes
*
* Arguments
*
*   N       INPUT         INTEGER           SCALAR
*     Input number
*
*   ISCR    SCRATCH       INTEGER           ARRAY(*)
*     Scratch array, length is N
*
*   IP      OUTPUT        INTEGER           ARRAY(*)
*     Array of primes, increasing order
*
*   NP      OUTPUT        INTEGER           ARRAY(*)
*     Number of primes found
*
* Common
*
*   None
*
* Errors
*
*   None
*
* Notes
*
*   None
*
* History
*
*   MNP      ,  04/12/88,  Mark Portney,  Initial release.
*
* End ******************************************************************
*
      INTEGER ISCR(*),IP(*)
      N = MAX0(N,3)
      DO 10 I = 1,N
            ISCR(I) = 0
   10 CONTINUE
      DO 20 I = 2,N-1
            DO 30 J =2*I,N,I
                  ISCR(J) = ISCR(J) + 1  ! NUMBER OF TIMES J IS FOUND
C                                    ! IN SIEVE
   30       CONTINUE
   20 CONTINUE
      J = 0
      DO 40 I = 2,N
            IF (ISCR(I) .EQ. 0) THEN
                  J = J + 1          ! NUMBER OF PRIMES
                  IP(J) = I          ! PRIMES
            END IF
   40 CONTINUE
      NP = J                         ! NUMBER OF PRIMES
      RETURN
      END
      SUBROUTINE VC64(F,T,CS,IFAX,INC,JUMP,N,LOT,IDIR)
      INTEGER IPRM(6),NPRM(6),IFAX(*)
      COMPLEX F(*),T(*),CS(*)
      SAVE
C     JUMP = 1 FOR EFFICIENCY
C     INC NOT EQ. 1 FOR EFFICIENCY
      DO 100 I=1,6
         IPRM(I)=IFAX(1+I)
         NPRM(I)=IFAX(7+I)
  100 CONTINUE
      JNVFWD=IDIR
      TEST = JNVFWD*AIMAG(CS(2))
      IF (TEST .GT. 0.) THEN
         DO 50 I = 1,N
            CS(I) = CONJG(CS(I))
   50    CONTINUE
      END IF
      IPC=0
      LF=N
      NW=1
      ISWT=1
  101 CONTINUE
         IPC=IPC+1
         IP=IPRM(IPC)
         NP=NPRM(IPC)
         DO 800 JP=1,NP
            LFOLD=LF
            LF=LF/IP
            NWOLD=NW
            NW=NW*IP
            IF(ISWT.GT.0)THEN
               CALL FFTFLY(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
            ELSEIF(ISWT.LT.0)THEN
               CALL FFTFLY(T,F,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
            ENDIF
            ISWT=-ISWT
  800    CONTINUE
      IF (LF.GT.1) GO TO 101
C
      IF (ISWT.LT.0) THEN
            LOT1 = 2*((LOT+1)/2)
            DO 951 I=1,N*LOT1
               F(I)=T(I)
  951       CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE VC999(DATA,WORK,TRIG,IFAX,INC,JUMP,LEN,LOT,IDIR)
*     LENV MUST MATCH IN FFTC1L,FFTC1,FFTCL,FFTC,VC999
      PARAMETER (lenv = 8)
      DIMENSION DATA(*),WORK(*),TRIG(*)
      IW = 1+2*LEN*MIN0(LOT,LENV)
      LDIR = - IDIR
      DO 10 I = 1,LOT,LENV
         NLOT = MIN0(LOT-I+1,LENV)
         DO 20 J = 1,NLOT
            I1 = 1+(I+J-2)*JUMP
            I2 = 1+2*(J-1)
            DO 30 K = 1,LEN
               WORK(I2+(K-1)*NLOT*2) = DATA(I1+(K-1)*INC)
   30       CONTINUE
            DO 31 K = 1,LEN
               WORK(1+I2+(K-1)*NLOT*2) = DATA(1+I1+(K-1)*INC)
   31       CONTINUE
   20    CONTINUE
         CALL VC64(WORK(1),WORK(IW),TRIG,IFAX,NLOT,1,LEN,NLOT,LDIR)
         DO 40 J = 1,NLOT
            I1 = 1+(I+J-2)*JUMP
            I2 = 1+2*(J-1)
            DO 50 K = 1,LEN
               DATA(I1+(K-1)*INC) = WORK(I2+(K-1)*NLOT*2)
   50       CONTINUE
            DO 51 K = 1,LEN
               DATA(1+I1+(K-1)*INC) = WORK(1+I2+(K-1)*NLOT*2)
   51       CONTINUE
   40    CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VFTFAX(N,IFAX,CS)
      INTEGER IFAX(*),IPRM(6),NPRM(6)
      COMPLEX*8 CS(*), COEF
      PI=3.14159265378
      PI2=PI*2.0
      DTH=-PI2/N
      C=1.
      S=0.
      CD=COS(DTH)
      SD=SIN(DTH)
      COEF=CMPLX(CD,SD)
      CS(1)=CMPLX(C,S)
      NSTEP=10
      DO 30 J=2,N,NSTEP
         IEND=J+NSTEP-1
         IF(IEND.GT.N)IEND=N
         DO 20 I=J,IEND-1
            CS(I)=COEF*CS(I-1)
   20    CONTINUE
         TH=-PI2*(IEND-1)/N
         CS(IEND)=CMPLX(COS(TH),SIN(TH))
   30 CONTINUE
      CALL FACTR ( N, IPRM, NPRM)
      DO 100 I=1,6
         IFAX(1+I)=IPRM(I)
         IFAX(7+I)=NPRM(I)
  100 CONTINUE
      RETURN
      END
      subroutine movlev(a,b,n)
      dimension a(*),b(*) 
      do 1 i=1,n
1     b(i)=a(i) 
      return
      end 
      subroutine clear(a,n) 
      dimension a(*)
      do 10 i=1,n 
10     a(i)=0.
      return
      end 
      subroutine rtran(a,b,nr,nc) 
      real a(*),b(*) 
      do 10 k=1,nc
      do 10 i=1,nr
      b((i-1)*nc+k)=a((k-1)*nr+i)
10    continue
      return
      end 
      subroutine fac235(nx,ifax)
c   ***   subroutine fac235 returns in ifax(1) the nearest even integer 
c   ***   >= nx whose only nontrivial factors are 2 3 or 5
c   ***   ifax(2), ifax(3), and ifax(4) contain the number of factors of 2
c   ***   3, and 5 respectively
      integer ifax(*)
      nxx=((nx-1)/2)*2
200   continue
      nxx=nxx+2
      i2=0
      i3=0
      i5=0
      n1=nxx
      n3=n1
300   continue
      n2=n1
      if(mod(n1,2).eq.0)then
      i2=i2+1
      n1=n1/2
      n3=2*n3
      end if
      if(mod(n1,3).eq.0)then
      i3=i3+1
      n1=n1/3
      n3=3*n3
      end if
      if(mod(n1,5).eq.0)then
      i5=i5+1
      n1=n1/5
      n3=5*n3
      end if
c   ***   if n1=1 then finished
      if(n1.eq.1)goto 1000
c   ***   if n1 > 1 but n1=n2 then nxx not factorable by 2,3,5
      if(n1.gt.1 .and. n1.eq.n2)goto 200
c   ***   if n1 > 1 and n1 < n2 then a 2 3 or 5 factor was found
c   ***   make another pass with this nxx
      if(n1.gt.1 .and. n1.lt.n2)goto 300
1000   continue
      ifax(1)=nxx
      ifax(2)=i2
      ifax(3)=i3
      ifax(4)=i5
      return
      end
