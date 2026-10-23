c-----------------------------------------------------------------------
c     Rotinas copiadas sem alteracao de test3dbcpmig.f (SEISPAK, programa
c     bcpmig) e spl1.f, para o modulo SPK_3DBCPMIG:
c     xxmvphs2d, xpspivels1, xpgsarr3d1, xpgsarr3d2, dfdx, boxf, zmutes,
c     embed, extractv, xxnmo, spl3, spl1.
c     (sgimrg2r, bndry, fftl2d, dwnn, fftc, vshift1, sgivshift1, fac235:
c      ../../spk_dvsmig3d/fortran)
c-----------------------------------------------------------------------
      subroutine xxmvphs2d(y,a,b,rs,nx,ny,dx,dy,dt,w,shft,vv,sgn,nvels)
      complex y(*),a(*),b(*) 
      real shft(*),vv(*),rs(*)
      real s1(20)
c   ***   find min and max shifts 
      nxy=nx*ny
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=0.
      do 5 i=1,nxy
      v1=amax1(v1,vv(i))
      v2=amin1(v2,vv(i))
      shft1=amin1(shft1,shft(i))
      shft2=amax1(shft2,shft(i))
5     continue
      vmax=1./sqrt(v2)
      vmin=1./sqrt(v1)
      shftmin=shft1
      shftmax=shft2
c   ***   fft, downward continue, and inverse fft for vmin
c   ***   fft, downward continue, and inverse fft for shift=vmin
      call fftl2d(y,b,nx,ny,+1)
c   ***   copy field to scratch array 
      do 7 i=1,nxy
7     a(i)=y(i) 
c   ***   if nvels = 1 the v1 = ( .5*(1./vmin+1./vmax) )**2
c   ***                shft1 = .5*(shftmin+shftmax)
      if(nvels .eq. 1)then
      v1=(.5*(1./vmin + 1./vmax) )**2
      shft1=.5*(shftmin+shftmax)
      end if
      call dwnn(y,nx,ny,dx,dy,w,shft1,v1,sgn)
      call fftl2d(y,b,nx,ny,-1)
c   ***   fft, downward continue, and inverse fft for shift=vmax
      if(vmax-vmin.lt..02*vmax .or. nvels .eq. 1)return
c   ***   velocity non-constant in this interval
c   ***   apply multi-velocity pspi
c   ***   make vmax slightly  larger than max of 1./sqrt(vv)
c   ***   make vmin slightly smaller than min of 1./sqrt(vv)
         dxz=vmin/shftmin
        vmax=1.0001*vmax
        vmin=.9999*vmin
      ichk=0
      dz=nint( shft1*dx*sqrt(v1))
      call xpspivels1(shft,rs,dz,dt,shftmin,shftmax,nxy,nvels1,ichk)
      shftinc=(shftmax-shftmin)/(nvels-1)
       if(nvels1.le.nvels)then
       do j=1,nvels1
       s1(j)=rs(j)
       end do
       else
       nvels1=nvels
       do j=1,nvels1
       s1(j)=shftmin+(j-1)*shftinc
       end do
       end if
      nsteps=nvels1-1
      call movlev(y,b,2*nxy)
      call clear(y,2*nxy)
c   ***   start multiple velocity loop
        vmin=.999*dxz*shftmin
        vmax=1.001*dxz*shftmax
      do istep=1,nsteps
      xs1=(s1(istep)-s1(1))/(shftmax-shftmin)
      xs2=(s1(istep+1)-s1(1))/(shftmax-shftmin)
      shft1=s1(istep)
      shft2=s1(istep+1)
      vv1=vmin+xs1*(vmax-vmin)
      vv2=vmin+xs2*(vmax-vmin)
      v1 = (1./vv1)**2
      v2 = (1./vv2)**2
c   ***   compute contribution of field downward continued using vv1
      do k=1,nxy
      rs(k) = ( vv2-1./sqrt(vv(k)) )/(vv2-vv1)
      rs(k) = amax1(0.,rs(k))
      if(rs(k) .gt. 1.)rs(k)=0.
      end do
c   ***   add contrubution of field downward continued using vv1
      do k=1,nxy
      y(k) = y(k) + rs(k)*b(k)
      end do
c   ***
      call movlev(a,b,2*nxy)
      call dwnn(b,nx,ny,dx,dy,w,shft2,v2,sgn)
c   ***   inverse fft data and scale
      call fftl2d(b,rs,nx,ny,-1)
c   ***   compute contribution of field downward continued using vv2
      do k=1,nxy
      rs(k)=(1./sqrt(vv(k))-vv1)/(vv2-vv1)
      rs(k)=amax1(0.,rs(k))
      if(rs(k) .ge. 1.)rs(k)=0.
      end do
c   ***   add contribution of field downward continued using vv2
      do k=1,nxy
      y(k) = y(k) + rs(k)*b(k)
      end do
c   ***
      end do
c   ***
      return
      end 
      subroutine xpspivels1(s1,s2,dz,dt,smin,smax,nx,nvels,icheck)
      real s1(*),s2(*)
      real vlist(2001)
      integer is(2001)
c   ***   clear  index array
c   ***   if smin nearly equal smax then set
c   ***   nvels=1    and   s2(1)=.5*(smin+smax)
         if(smax-smin .lt. .01*smax )then
         nvels=1
         s2(1)=.5*(smin+smax)
         return
         end if
c   ***
      ns=2001
      do j=1,ns
      is(j)=0
      vlist(j)=20000.-10.*(j-1)
      end do
       dzt=dz/dt
c   ***   set index array for histogram of velocity profile
      do j=1,nx
      vv=dzt/s1(j)
       ii=1+nint( (20000.-vv)*.1 )
       is(ii)=is(ii)+1
      end do
c   ***   compute histogram of velocity profile
      jj=0
       imax=0
      do j=2,ns
      imax=max0(is(j),imax)
      end do
      npct=imax/10
      vmn=10.*nint(.1*dzt/smax)
      vmx=10.*nint(.1*dzt/smin)
      jvmn=1+nint( (20000.-vmn)*.1)
      jvmx=1+nint( (20000.-vmx)*.1)
      do j=2,ns
      if(is(j).ge.npct .or. j.eq.jvmn .or. j.eq.jvmx)then
      jj=jj+1
      s2(jj)=dzt/vlist(j)
      end if
      end do
      nvels=jj
      s2(1)=smin
      s2(nvels)=smax
102   format(2i4,8f8.0)
103   format(8i10)
      return
      end
      subroutine xpgsarr3d1(t1,t2,s,v,xk,
     1                  nx,itzr,dt,ddz,dx,sxpos,ny,dy,sypos)
      real t1(nx,ny),t2(nx,ny),s(nx,ny),v(nx,ny)
      real s1(1000),s2(1000),s3(1000)
      save xdz,vdt,vmax,vmin
c
c   ***   compute arrival times at level k by downward continuation
c   ***   using Eichonal equation.
c
c   From an original by Dwight Sukup of Texaco E&P Research  
c
c   t1    - input arrival time in samples at level k-1
c   t2    - output arrival time in samples at level k
c   s     - scratch array
c   v     - velocity at level k
c   xk    - depth of shot in intervals of dz; i.e., depth samples
c   nx    - number of x (xline) samples in slice
c   itzr  - number of interpolation steps per downward-continuation step
c   dt    - delta t in seconds
c   ddz   - delta z in meters
c   dx    - delta x in meters
c   sxpos - x (xline) trace number of source position
c   ny    - number of y (subline) samples in slice
c   dy    - delta y in meters
c   sypos - y (subline) trace number of source position
c
      k=xk
      dz=itzr*ddz
      if(k.lt.1)return
      if(k.eq.1)then 
        xdz=(xk-k)*dz
        isx=sxpos 
        isy=sypos
        vdt=v(isx,isy)*dt
        vmax=0.  
        vmin=1000000. 
      end if
      do 10 i=1,nx  
        do 10 j=1,ny  
          vmin=amin1(vmin,v(i,j)) 
          vmax=amax1(vmax,v(i,j)) 
10    continue
      if(vmax-vmin .lt. .01*vmax)then
c
c   ***   compute arrival times for z=(k-1)*dz and z=k*dz
c
        c1=( (k-1)*dz+xdz)**2
        c2=(     k*dz+xdz)**2
        c3=1./vdt
        do 63 j=1,ny  
          do 63 i=1,nx
            r=((i-sxpos)*dx)**2 + ((j-sypos)*dy)**2
            t1(i,j)=c3*sqrt(r+c1)
            t2(i,j)=c3*sqrt(r+c2)
63      continue
      else
c   
c   ***   downward continue arrival times in samples by Eichonal. 
c  
        do 71 j=1,ny  
          do 71 i=1,nx 
            t1(i,j)=t2(i,j)  
71        continue  
       nitzr=2*itzr
       ritzr=1./nitzr
       do n=1,nitzr
c         
c   accumulate (dt/dx)**2 + (dt/dy)**2  into s(i,j)        
c
c   ***   compute dt/dx  
          do j=1,ny
            do i=1,nx
            s1(i)=t2(i,j)
            end do
            call dfdx(s1,s2,s3,dx,nx)
            do i=1,nx
            s(i,j)=s2(i)*s2(i)
            end do
          end do
c   ***   compute dt/dy 
          if(ny .gt. 1)then
          do i=1,nx
            do j=1,ny
            s1(j)=t2(i,j)
            end do
            call dfdx(s1,s2,s3,dy,ny)
            do j=1,ny
            s(i,j)=s(i,j)+s2(j)*s2(j)
            end do
          end do
          end if
c
c   ***   t2(i,j)=t2(i,j)+ritzr*dz*sqrt( (1/v)**2-(dt/dx)**2-(dt/dy)**2 )
c
          do j=1,ny
            do i=1,nx
            vdtr=1./(v(i,j)*dt)
            r=vdtr*vdtr-s(i,j)
            t2(i,j)=t2(i,j) + ritzr*dz*sign(1.,r)*sqrt(abs(r))
            end do
          end do
      end do
c   ***   smooth t2
      call boxf(t2,s,nx,ny,32.)
      call movlev(s,t2,nx*ny)
c   ***
      end if
      return
      end 
      subroutine xpgsarr3d2(t1,t2,s,v,xk,
     1                  nx,itzr,dt,ddz,dx,sxpos,ny,dy,sypos)
      real t1(nx,ny),t2(nx,ny),s(nx,ny),v(nx,ny)
      real s1(1000),s2(1000),s3(1000)
      save xdz,vdt,vmax,vmin
c
c   ***   compute arrival times at level k by downward continuation
c   ***   using Eichonal equation.
c
c   From an original by Dwight Sukup of Texaco E&P Research  
c
c   t1    - input arrival time in samples at level k-1
c   t2    - output arrival time in samples at level k
c   s     - scratch array
c   v     - velocity at level k
c   xk    - depth of shot in intervals of dz; i.e., depth samples
c   nx    - number of x (xline) samples in slice
c   itzr  - number of interpolation steps per downward-continuation step
c   dt    - delta t in seconds
c   ddz   - delta z in meters
c   dx    - delta x in meters
c   sxpos - x (xline) trace number of source position
c   ny    - number of y (subline) samples in slice
c   dy    - delta y in meters
c   sypos - y (subline) trace number of source position
c
      k=xk
      dz=itzr*ddz
      if(k.lt.1)return
      if(k.eq.1)then 
        xdz=(xk-k)*dz
        isx=sxpos 
        isy=sypos

        vdt=v(isx,isy)*dt
        vmax=0.  
        vmin=1000000. 
        c1=( (k-1)*dz+xdz)**2
        c2=(     k*dz+xdz)**2
        c3=1./vdt
         do  j=1,ny  
          do  i=1,nx
            r=((i-sxpos)*dx)**2 + ((j-sypos)*dy)**2
            t1(i,j)=c3*sqrt(r+c1)
            t2(i,j)=c3*sqrt(r+c2)
        enddo
        enddo
      end if
        isx=sxpos 
        isy=sypos
      kx1=max0(isx-k,1)
      kx2=min0(isx+k,nx)
      ky1=max0(1,isy-k)
      ky2=min0(ny,isy+k)
      do j=ky1,ky2
      do i=kx1,kx2
          vmin=amin1(vmin,v(i,j)) 
          vmax=amax1(vmax,v(i,j)) 
      enddo
      enddo
c     ivij=nint(v(isx,isy))
c     ivmin=nint(vmin)
c     ivmax=nint(vmax)
c     write(*,*)'xpgs ',k,kx1,kx2,ky1,ky2,ivij,ivmin,ivmax
      if(vmax-vmin .lt. .01*vmax)then
c
c   ***   compute arrival times for z=(k-1)*dz and z=k*dz
c
        c1=( (k-1)*dz+xdz)**2
        c2=(     k*dz+xdz)**2
        c3=1./vdt
        do 63 j=1,ny  
          do 63 i=1,nx
            r=((i-sxpos)*dx)**2 + ((j-sypos)*dy)**2
            t1(i,j)=c3*sqrt(r+c1)
            t2(i,j)=c3*sqrt(r+c2)
63      continue
      else
c   
c   ***   downward continue arrival times in samples by Eichonal. 
c  
        do 71 j=1,ny  
          do 71 i=1,nx 
            t1(i,j)=t2(i,j)  
71        continue  
       nitzr=2*itzr
       ritzr=1./nitzr
       do n=1,nitzr
c         
c   accumulate (dt/dx)**2 + (dt/dy)**2  into s(i,j)        
c
c   ***   compute dt/dx  
          do j=1,ny
            do i=1,nx
            s1(i)=t2(i,j)
            end do
            call dfdx(s1,s2,s3,dx,nx)
            do i=1,nx
            s(i,j)=s2(i)*s2(i)
            end do
          end do
c   ***   compute dt/dy 
          if(ny .gt. 1)then
          do i=1,nx
            do j=1,ny
            s1(j)=t2(i,j)
            end do
            call dfdx(s1,s2,s3,dy,ny)
            do j=1,ny
            s(i,j)=s(i,j)+s2(j)*s2(j)
            end do
          end do
          end if
c
c   ***   t2(i,j)=t2(i,j)+ritzr*dz*sqrt( (1/v)**2-(dt/dx)**2-(dt/dy)**2 )
c
          do j=1,ny
            do i=1,nx
            vdtr=1./(v(i,j)*dt)
            r=vdtr*vdtr-s(i,j)
            t2(i,j)=t2(i,j) + ritzr*dz*sign(1.,r)*sqrt(abs(r))
            end do
          end do
      end do
c   ***   smooth t2
c      call boxf(t2,s,nx,ny,32.)
c      call movlev(s,t2,nx*ny)
c   ***
      end if
      return
      end 
       subroutine dfdx(s1,s2,s3,dx,nx)
       real s1(*),s2(*),s3(*)
       dxr=1./dx
       s2(1)=(s1(2)-s1(1))*dxr
       s2(nx)=(s1(nx)-s1(nx-1))*dxr
       do j=2,nx-1
       dt1=s1(j)-s1(j-1)
       dt2=s1(j+1)-s1(j)
       dt3=.5*(dt1+dt2)
       if(dt1.ge.0.)dt=dt1
       if(dt2.le.0.)dt=dt2
       if(dt1*dt2 .lt. 0.)dt=dt3
       s2(j)=dt*dxr
       end do
c   ***   smooth
c      do j=1,nx
c      s3(j)=s2(j)
c      end do
c      s2(1)=s3(1)
c      s2(nx)=s3(nx)
c      do j=2,nx-1
c      s2(j)=.25*(s3(j+1)+s3(j-1))+.5*s3(j)
c      end do
       return
       end
      subroutine boxf(u,s,nx,ny,a)
      real u(nx,ny),s(nx,ny)
      d=a+8
      c1=1./d
      c2=a/d
      do j=2,ny-1
      do i=2,nx-1
      s(i,j)=c1*(u(i+1,j-1)+u(i+1,j)+u(i+1,j+1)+u(i,j-1)+u(i,j+1)
     : +u(i-1,j-1)+u(i-1,j)+u(i-1,j+1))+c2*u(i,j)
      end do
      end do
      d=a+5
      c1=1/d
      c2=a/d
      do i=2,nx-1
      s(i,1)=c1*(u(i+1,1)+u(i-1,1)+u(i+1,2)+u(i,2)+u(i-1,2))
     :      +c2*u(i,1)
      s(i,ny)=c1*(u(i+1,ny)+u(i-1,ny)+u(i+1,ny-1)+u(i,ny-1)
     :      +u(i-1,ny-1)) + c2*u(i,ny)
      end do
      do j=2,ny-1
      s(1,j)=c1*(u(1,j+1)+u(1,j-1)+u(2,j-1)+u(2,j)+u(2,j+1))
     :      +c2*u(1,j)
      s(nx,j)=c1*(u(nx,j+1)+u(nx,j-1)+u(nx-1,j-1)+u(nx-1,j)
     :       +u(nx-1,j+1))+c2*u(nx,j)
      end do
      d=a+3
      c1=1./d
      c2=a/d
      s(1,1)=c1*(u(1,2)+u(2,2)+u(2,1))+c2*u(1,1)
      s(1,ny)=c1*(u(1,ny-1)+u(2,ny-1)+u(2,ny))+c2*u(1,ny)
      s(nx,1)=c1*(u(nx,2)+u(nx-1,1)+u(nx-1,2))+c2*u(nx,1)
      s(nx,ny)=c1*(u(nx,ny-1)+u(nx-1,ny-1)+u(nx-1,ny))+c2*u(nx,ny)
      do j=1,ny
      do i=1,nx
      u(i,j)=s(i,j)
      end do
      end do
      return
      end
       subroutine zmutes(data,zmut,nx,ny,dx,dy,ddz,itzr,z,al1,al2,k)
       real data(*),zmut(*)
c   ***   apply amplitude restoration and zmutes
       nxy=nx*ny
       dx10=10.*dx
       do j = 1,itzr
       f1=float(j-1)/float(itzr)
       f2=1.-f1
       rmute=f1*zmut(k+1)+f2*zmut(k)
       zz=z+(j-1)*ddz
        do jy=1,ny
        do jx=1,nx
        jxy=(jy-1)*nx+jx
        xx=(al1-jx)*dx
        yy=(al2-jy)*dy
        r=sqrt(xx**2+yy**2+zz**2)
        rxy=sqrt( ( (jx-al1)*dx)**2 + ( (jy-al2)*dy)**2)
        amp=r*amax1(0. ,amin1(1., (rmute+dx10-rxy)/dx10 ) )
        data((j-1)*nxy+jxy)= amp*data((j-1)*nxy+jxy)
c        if(rxy.gt.rmute)data(jj+jxy)=0.
        end do
        end do
        end do
        end 
      subroutine embed(yin,yout,nx,ny,nxout,nyout,ixout,iyout)
      real yin(nx,ny)
      real yout(nxout,nyout)
c   ***   embed the nx by ny data plate contained in yin into
c   ***   the nxout by nyout data plate contained in yout
c   ***   starting at x=ixout , y=iyout
      do iy=1,ny
      iy1=iy+iyout-1
      if(1 .le. iy1 .and. iy1 .le. nyout)then
c   ***   case for ixout .ge. 1
      if(ixout .ge. 1)then
      ix1=ixout
      ix2=min0(ixout+nx-1,nxout)
      nn=ix2-ix1+1
      if(nn .ge. 1)then
      do k=1,nn
      yout(ixout+k-1,iy1)=yout(ixout+k-1,iy1)+yin(k,iy)
      end do
      end if
      end if
c   ***   case for ixout .lt. 1
      if(ixout .lt. 1)then
      ix1=1
      ix2=min0(ixout+nx-1,nxout)
      nn=ix2-ix1+1
      if(nn .ge. 1)then
      do k=1,nn
      yout(k,iy1)=yout(k,iy1)+yin(1-ixout+k,iy)
      end do
      end if
      end if
c   ***
      end if
      end do
      return
      end
      subroutine extractv(vin,vout,nxin,nyin,nxout,nyout,ix0,iy0)
c   ***   extract submatrix of values
      real vin(nxin,nyin),vout(nxout,nyout)
      do iy=1,nyout
      iy1=max0(iy0+iy-1,1)
      iy1=min0(iy1,nyin)
      do ix=1,nxout
      ix1=max0(ix0+ix-1,1)
      ix1=min0(ix1,nxin)
      vout(ix,iy)=vin(ix1,iy1)
      end do
      end do
      return
      end
       subroutine xxnmo(u,s,nt,dt,x1,x2,v)
c   ***  apply (normal) moveout from offset x1 to offset x2
c   ***  if x2=0 then normal movout is applied
c   ***  if x1=0 then normal movein is applied
c   ***  if x1>0 and x2>0 then data is moved out at x1 and moved in at x2
       real u(*),s(*)
       real v(*)
       t0=0.
       i1=1
       i2=i1+nt
       i3=i2+nt
c   ***   compute nmo times
       eps=.01*dt
       t1=0.
       t2=0.
       do j=1,nt
       t11=sqrt( t0**2 + (x1/v(j))**2 )
       t22=sqrt( t0**2 + (x2/v(j))**2 )
       t1=amax1(t1+eps,t11)
       t2=amax1(t2+eps,t22)
       s(i1+j-1)=t1
       s(i2+j-1)=t2
       t0=t0+dt
       end do
c   ***   call spl3
      call spl3(s(i2),u,s(i1),s(i3),nt,nt)
       do j=1,nt
       u(j)=s(i3+j-1)
       end do
       return
       end
      subroutine spl3(x,y,xx,f,n,m)
      real a(5000),b(5000),c(5000),d(5000)
      real y(*),x(*),f(*),xx(*)
c       cubic spline interpolation. input points are (x(i),y(i)),i=1,n.  
c     output points are (xx(i),f(i)),i=1,m. output independent variables 
c     xx(i) must safisfy x(1)<xx(i). if xx(i)>x(n)then f(i)=0. 
      n1=n-1 
      n2=n-2 
      n3=n-3 
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
      f(i)=0.
13     continue  
      del1=xx(i)-x(j)  
      del2=x(j+1)-xx(i)  
      if(del1.ge.0..and.del2.ge.0.)go to 12  
      if(del1 .lt. 0.)goto 14
      if(j.eq.n1)go to 12  
      j=j+1  
      go to 13 
12     continue  
      f(i)=y(j)+del1*(c(j)+del1*(b(j)+del1*a(j)))  
14     continue
11     continue  
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
