c     SEASEIS: analise de velocidade por focalizacao do test3dbcpmig
c     (Texaco/Geonunes). xvelcoh, xz2t, spl2 e xxtau2v copiados sem
c     mudanca de test3dbcpmig.f; boxflt.f original da biblioteca
c     SEISPAK; SPKBFOC substitui o trecho do
c     programa principal que le TDAT/velscr e grava OFFLINE4/5 e VAVG.
c
      SUBROUTINE SPKBFOC(TDAT,VELS,NXOUT,NYOUT,NZ,DZ,NT,DT,JNT,
     :  NGSNX,NGSNY,IFORX,IFORY,IRNIX,IRNIY,VMN,VMX,
     :  OTAU,OVC,OVAVG,ZD,TD,S)
c     TDAT(jnt*ngsns,nz) acumulado; para cada ponto de analise:
c       OTAU (nt x jnt)  desvio de tau (OFFLINE4)
c       OVC  (nt x jnt)  coerencia x velocidade (OFFLINE5)
c       OVAVG(nt)        velocidade media (VAVG)
      REAL TDAT(JNT*NGSNX*NGSNY,*),VELS(NXOUT*NYOUT,*)
      REAL OTAU(NT*JNT,*),OVC(NT*JNT,*),OVAVG(NT,*),ZD(*),TD(*),S(*)
      REAL VEL(20001)
      ngsns=ngsnx*ngsny
      sr=1000.*dt
      twindow=(jnt-1)*sr
      lf=2
      do iy=1,ngsny
      do ix=1,ngsnx
      dpny=(iy-1)*irniy+ifory
      dpnx=(ix-1)*irnix+iforx
      ii=(iy-1)*ngsnx + ix
c   ***   tau-deviation data for this point: zdata(nz,jnt)
      do j=1,jnt
      do k=1,nz
      zd((j-1)*nz+k)=tdat((ii-1)*jnt+j,k)
      enddo
      enddo
      idpnx=max0(1,min0(nint(dpnx),nxout))
      idpny=max0(1,min0(nint(dpny),nyout))
      ivrec=nxout*(idpny-1)+idpnx
      do k=1,nz
      vel(k)=vels(ivrec,k)
      enddo
      ivd=1+nt*jnt
      call xvelcoh(zd,td,td(ivd),vel,s,twindow,vmn,vmx,
     :              nt,dt,nz,dz,lf)
      do j=1,nt
      ovavg(j,ii)=s(j)
      enddo
      do j=1,nt*jnt
      otau(j,ii)=td(j)
      ovc(j,ii)=td(ivd+j-1)
      enddo
      end do
      end do
      RETURN
      END
c
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
       subroutine xvelcoh(zdata,tdata,vdata,vel,s,
     : win,vmin,vmax,nt,dt,nz,dz,lf)
c *** zdata  real  input tau-deviaton data (nz x jnt = win/sr+1)
c *** tdata  real  output tau dev data (nt x jnt)
c *** vdata  real  output velocity coherence (nt x jnt)
c *** vel    real  input interval velocity (nz samples at dz)
c *** s      real  scratch 
c *** win    real  window centered at t=0 in milliseconds
c *** vmin   real  min velocity to process
c *** vmax   real  max velocity to process
c *** nt     real  # time samples for input/output
c *** dt     real  time sample interval in seconds
c *** nz     real  # depth samples for input data
c *** dz     real  depth sample interal in depth units
c *** lf     int   length of box filter is 2*lf+1
       real zdata(*),tdata(*),vdata(*),vel(*),s(*)
c ***   square data
       sr=1000*dt
       jnt=win/sr+1
       do j=1,nz*jnt
       tdata(j)=zdata(j)**2
       end do
c   ***   apply box filter to smooth
      call boxflt(tdata,zdata,s,nz,jnt,lf)
      nzt=max0(nt,nz)
      i1=1
      i2=i1+nzt
      mjnt=jnt/2+1
      do j=1,jnt
      j1=(j-1)*nz+1
      j2=(j-1)*nt+1
c   ***   convert from depth -> time
       call xz2t(zdata(j1),s,vel,s(i2),nz,dz,nt,dt)
c   ***   xz2t returns interval velocity in s(i2)
c   ***           sampled at (j-1)*dt for j=1,nt
c   ***   apply skewing transform
      do n=1,nt
      vdata(j2+n-1)=0.
      end do
      ishift=mjnt-j
      if(ishift.gt.0)then
c   ***  pull trace up by ishift samples for ishift>0
      do n=1+ishift,nt
      vdata(j2+n-1-ishift)=s(n)
      end do
      else
c   ***  push trace down by ishift samples for ishift<0
      do n=1,nt+ishift
      vdata(j2-ishift+n-1)=s(n)
      end do
      end if
      end do
c   ***   compute rms velocity trace from interval velocity trace in s(i2)
c     write(*,*)'interval velocity trace '
c     write(*,1234)(nint(s(i2+n-1)),n=1,nt)
c   ***   compute avg velocity trace from interval velocity trace in s(i2)
      t=0.
      vavg=s(1)
      do n=1,nt
      vavg=(t*vavg+dt*s(i2+n-1))/(t+dt)
      s(n)=vavg
      t=t+dt
      end do
c   ***   write avg velocity

c     write(*,*)'avg velocity in xvelcoh '
c     write(*,1234)(nint(s(n)),n=1,nt)
1234  format(10i8)
c   ***   smooth again by applying boxflt
      j1=i1+nt
      call boxflt(vdata,tdata,s(j1),nt,jnt,lf)
c   ***   scale data so largest value is 1.
c     call scl(tdata,1.,nt*jnt)
      iv=1
      i1=iv+nt
      i2=i1+nt
      call xxtau2v(tdata,vdata,s(iv),vmin,vmax,nt,dt,jnt)
      return
      end
      subroutine xz2t(y,u,v,s,nz,dz,nt,dt)
c   ***   convert from depth to time
c   ***   input  y real nz samples at sample rate dz 
c   ***   ouput  u real nt samples at sample rate dt
c   ***   input  v real  nz interval velocities sampled at dz
c   ***   input  s real  nt+nz
c   ***   output s real nt interval velocities sampled at dt
c   ***   input nz  int # depth samples
c   ***   input dz  real depth sample interval in depth units
c   ***   input nt  int # time samples
c   ***   input dt  real time sample interval in seconds
      real u(*),y(*),s(*),v(*)
      dzr=1./dz
c   ***   set dt0 to one-way sample interval in seconds
      dt0=dt/2
      z=0.
      t=0.
      z1=0.
      z2=dz
      v1=v(1)
      v2=v(2)
      iz=1
      s(1)=z1
      do n=2,nt
      t=t+dt0
      zz=z+v1*dt0
      if(zz.lt.z2)then
      z=zz
      zr=z*dzr
      iz=int(zr+1.)
      else
      dt1=(z2-z)/v1
      dt2=dt0-dt1
      z=z+v1*dt1+v2*dt2
      zr=z*dzr
      iz=iz+1
      z1=z2
      z2=z1+dz
      v1=v2
      v2=v(min0(iz+1,nz))
      end if
      s(n)=z
      end do
c   ***   s contains depths at equal time samples : s(j)=z((j-1)*sr)
c   ***   spline interpolation is used to interpolate y at array s
      call spl2(y,dz,u,s,s(nt+1),nz,nt,1.)
c   ***   depths at equal time samples are converted to interval
c   ***   velocities at equal time samples
      dtr=1./dt0
      z2=s(1)
      s(1)=v(1)
      do j=2,nt
      z1=z2
      z2=s(j)
      s(j)=(z2-z1)*dtr
      end do
      return
      end
      subroutine spl2(y,hy,f,xf,s,n,m,scf)
c  resample n samples of y whose sample rate is hy to m samples of f
c   at xf(i),i=1,m .  where  0 <= xf(i) <= xf(i+1) <= (n-1)*hy
c   spline interpolation using free endpoint condition is used
c  **  y   real   input (n samples at sample interval hy)
c  **  hy  real   input  sample interval for y
c  **  f   real   output (m samples at sample locations xf(i)
c  **  xf  real   input   f(i)=y(xf(i)) i=1,m
c  **  s   real   scratch  max(n,m)
c  **  n   integer  # input samples of y
c  **  m   integer  # output samples for f
c  **  scf real   scale factor to apply 
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
      subroutine xxtau2v(tdata,vdata,vel,vmin,vmax,nt,dt,jnt)
c  ***  tdata     real   input (nt x jnt) tau deviation data 
c  ***  vdata     real   output (nt x jnt) estimated velocity  
c  ***  vel       real   input rms velocity sampled at (j-1)*dt
c  ***  vmin      real   minimum velocity to map
c  ***  vmax      real   maximum velocity to map
c  ***  nt        int    number of time samples
c  ***  dt        real   time sample interval in seconds
c  ***  jnt       int    number of tau deviation data traces
      real vdata(nt,jnt), tdata(nt,jnt),vel(nt)
c   ***   map tau-deviation data (dtau) to corrected rms velocity (vc)
c   ***   Corrected rms velocity is estimated from tau deviation using:
c   ***   VC = VM * (1. - DTAU/TAU)
c   ***   DTAU = TAU * (1. - (VC/VM)**2 )
      dv=(vmax-vmin)/(jnt-1)
c   ***   dtau=0 at tau-deviaton trace  mjnt = jnt/2+1
      mjnt=jnt/2+1
      do j=1,jnt
      do n=1,nt
      vc=vmin+(j-1)*dv
      vm=vel(n)
      tau=(n-1)*dt
      dtau=tau*(1.-(vc/vm)**2)
c   ***  set rj to trace # for this dtau
      rj=dtau/dt+mjnt
      jj=int(rj)
      pctdv=abs((vm-vc)/vm)*100.
      f1=float(n)/float(nt)
      pct=f1*40+(1.-f1)*20
      if((1.le.jj .and. jj.lt.jnt).and.(pctdv.le.pct))then
      r1=rj-jj
      r2=1.-r1
      vdata(n,j)=r1*tdata(n,jj+1)+r2*tdata(n,jj)
      else
      vdata(n,j)=0.
      end if
c   ***   set dtau=0 location to zero
      if(jj .eq. mjnt)vdata(n,j)=0.
      if( jj .eq. mjnt-1 .and. dt .lt. .007)vdata(n,j)=0
      end do
      end do
      return
      end
