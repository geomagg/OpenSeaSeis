c-----------------------------------------------------------------------
c     Rotinas do caminho 3D de dvsmig.f (SEISPAK) usadas por SPK_DVSMIG3D:
c     sgimrg2r, cosphs3d, dwncos3d, fftl3d, bndry (sem alteracao)
c-----------------------------------------------------------------------
      subroutine sgimrg2r(z,u,y,shft,t1,t2,nl,nh,nt1,nx,itzr) 
      complex a,b,c,s1,s2,s3,s4,s5,s6
      real u(*),y(*),z(*),shft(*),t1(*),t2(*)
      pi=3.1415927
      rz=1./itzr
      nws=nh-nl+1 
      delw=2.*pi/nt1
      nl1=nl+1
      wf=(nl-1)*delw
      do 200 i=1,nx*itzr
      z(i)=0.
 200  continue
c    ***    initialize frequencies
      do 100 i=1,nx 
      sr1=cos(delw*t1(i))
      sr2=cos(wf*t1(i)) 
      sr3=cos(delw*(t1(i)-shft(i))) 
      sr4=cos(wf*(t1(i)-shft(i))) 
      sr5=cos(delw*rz*(t2(i)-t1(i)+shft(i)))  
      sr6=cos(wf*rz*(t2(i)-t1(i)+shft(i))) 
      si1=sin(delw*t1(i))
      si2=sin(wf*t1(i)) 
      si3=sin(delw*(t1(i)-shft(i))) 
      si4=sin(wf*(t1(i)-shft(i))) 
      si5=sin(delw*rz*(t2(i)-t1(i)+shft(i)))  
      si6=sin(wf*rz*(t2(i)-t1(i)+shft(i))) 
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
      do 10 l=0,itzr-1
      f1=(itzr-l)*rz
      f2=1.-f1
      z(ll+i)=z(ll+i)+f1*(ar*u(jr)-ai*u(ji))+f2*(br*y(jr)-bi*y(ji))
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
      subroutine cosphs3d(y,a,b,rs,nx,nh,ny,dx,dh,dy,w,shft,vv,sgn,nv)
      complex y(*),a(*),b(*) 
      real shft(*),vv(*),rs(*)
c   ***   find min and max shifts 
      nxhy=nx*nh*ny
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=0.
      do 5 i=1,nxhy
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
      call fftl3d(y,a,nx,nh,ny,+1)
c   ***   copy field to scratch array 
      do 7 i=1,nxhy
7     a(i)=y(i) 
c   ***   if nv = 1 the v1 = ( .5*(1./vmin+1./vmax) )**2
c   ***                shft1 = .5*(shftmin+shftmax)
      if(vmax-vmin.lt..005*vmax .or. nv .eq. 1)then
      v1=(.5*(1./vmin + 1./vmax) )**2
      shft1=.5*(shftmin+shftmax)
      call dwncos3d(y,nx,nh,ny,dx,dh,dy,w,shft1,v1,sgn)
      call fftl3d(y,b,nx,nh,ny,-1)
      return
      endif
c   ***   fft, downward continue, and inverse fft for shift=vmax
c   ***   velocity non-constant in this interval
c   ***   apply multi-velocity pspi
c   ***   make vmax slightly  larger than max of 1./sqrt(vv)
c   ***   make vmin slightly smaller than min of 1./sqrt(vv)
        vmax=1.0001*vmax
        vmin=.9999*vmin
      nsteps=(vmax-vmin)/(.2*vmax) +1
      nsteps=nv-1
      vvinc= (vmax-vmin)/nsteps
      shftinc=(shftmax-shftmin)/nsteps
      call movlev(y,b,2*nxhy)
      call clear(y,2*nxhy)
c   ***   start multiple velocity loop
        do k=1,nxhy
        rs(k)=1./sqrt(vv(k))
        end do
      do istep=1,nsteps
      vv1=vmin+(istep-1)*vvinc
      vv2=vv1+vvinc
      v1 = (1./vv1)**2
      v2 = (1./vv2)**2
      shft1=shftmin+(istep-1)*shftinc
      shft2=shft1+shftinc
      rdv=1./(vv2-vv1)
c   ***   compute contribution of field downward continued using vv1
      do k=1,nxhy
      rss = (vv2-rs(k))*rdv
      if(0. .lt. rss .and. rss .le. 1.)then
      xr=real(y(k))+rss*real(b(k))
      xi=aimag(y(k))+rss*aimag(b(k))
      y(k)=cmplx(xr,xi)
      end if
      end do
c   ***
      call movlev(a,b,2*nxhy)
      call dwncos3d(b,nx,nh,ny,dx,dh,dy,w,shft2,v2,sgn)
c   ***   inverse fft data and scale
      call fftl3d(b,rs(nxhy+1),nx,nh,ny,-1)
c   ***   compute contribution of field downward continued using vv2
      do k=1,nxhy
      rss = (rs(k)-vv1)*rdv
      if(0. .lt. rss .and. rss .lt. 1.)then
      xr=real(y(k))+rss*real(b(k))
      xi=aimag(y(k))+rss*aimag(b(k))
      y(k)=cmplx(xr,xi)
      end if
      end do
c   ***
      end do
c   ***
      return
      end 
      subroutine dwncos3d(b,nx,nh,ny,dx,dh,dy,w,shft,vv,sgn)
*     b has length 2*nx*nh*ny
      complex b(nx,nh,ny),e
      real kx,kh,ky
      if (w .eq. 0.) then
         do iy=1,ny
         do ih=1,nh
         do ix=1,nx 
         b(ix,ih,iy)=0.
         enddo
         enddo
         enddo
         return
      end if
      dkx = 2.*3.141593/(nx*dx)
      dkh = 2.*3.141593/(nh*dh)
      dky = 2.*3.141593/(ny*dy)
      t = w*shft
      vvrw2 = vv/w**2 
      do iy=1,ny
      ky=min0(iy-1,ny-iy+1)*dky
      do ih=1,nh
      kh=min0(ih-1,nh-ih+1)*dkh
      do ix=1,nx
      kx=min0(ix-1,nx-ix+1)*dkx
      r=1.-.25*vvrw2*((kx+sgn*kh)**2 + ky**2)
            if (r .lt. .001) then 
               b(ix,ih,iy)=0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               e=cmplx(er,ei)
               b(ix,ih,iy)=e*b(ix,ih,iy)
            endif
      enddo
      enddo
      enddo
      return
      end 
      subroutine fftl3d(y,a,nz,nx,ny,isgn)
      complex y(1),a(1)
c   ***   complex 3d matrix y(ny,nx,nz) in
c   ***   forward 3d fft returned in y(ny,nx,nz) if isgn=1
c   ***   inverse 3d fft returned in y(ny,nx,nz) if isgn=-1
      if(isgn.eq.1)then
c   ***   fft in z direction
      call cffti(nz,a)
      do k=1,nx*ny
      ii=(k-1)*nz+1
      call cfftf(nz,y(ii),a)
      enddo
c   ***   fft in x direction
      if(nx.gt.1)then
      do k=1,ny
      ii=(k-1)*nx*nz+1
      call ctran(y(ii),a,nz,nx)
      call cffti(nx,y(ii))
      do j=1,nz
      iii=(j-1)*nx+1
      call cfftf(nx,a(iii),y(ii))
      enddo
      call ctran(a,y(ii),nx,nz)
      enddo
      endif
c   ***   fft in y direction
      if(ny.gt.1)then
      call cffti(ny,a(ny+1))
      do k=1,nx*nz
      do j=1,ny
      a(j)=y((j-1)*nx*nz+k)
      enddo
      call cfftf(ny,a,a(ny+1))
      do j=1,ny
      y((j-1)*nx*nz+k)=a(j)
      enddo
      enddo
      endif
      endif
      if(isgn.eq.-1)then
c   ***   fft in z direction
      call cffti(nz,a)
      do k=1,nx*ny
      ii=(k-1)*nz+1
      call cfftb(nz,y(ii),a)
      enddo
c   ***   fft in x direction
      if(nx.gt.1)then
      do k=1,ny
      ii=(k-1)*nx*nz+1
      call ctran(y(ii),a,nz,nx)
      call cffti(nx,y(ii))
      do j=1,nz
      iii=(j-1)*nx+1
      call cfftb(nx,a(iii),y(ii))
      enddo
      call ctran(a,y(ii),nx,nz)
      enddo
      endif
c   ***   fft in y direction
      if(ny.gt.1)then
      call cffti(ny,a(ny+1))
      do k=1,nx*nz
      do j=1,ny
      a(j)=y((j-1)*nx*nz+k)
      enddo
      call cfftb(ny,a,a(ny+1))
      do j=1,ny
      y((j-1)*nx*nz+k)=a(j)
      enddo
      enddo
      endif
      endif
      if(isgn.eq.-1)then
      scale=1./(float(nx)*float(ny)*float(nz))
      do j=1,nx*ny*nz
      y(j)=scale*y(j)
      end do
      end if
      return
      end
      subroutine bndry(y,nx,ny,lbx,lby)
      complex y(*)
c   ***   apply attenuating absorbing boundary
c   ***   y     nx by ny complex field for one frequency 
c   ***   nx    # xlines
c   ***   ny    # ilines
c   ***   lbx   # of traces to attenuate at each end of inlines
c   ***   lby   # of traces to attenuate at each end of xlines
c   ***   lby trace attenuation region
c   ***   attenuate x boundary regions
      do iy=1,ny
         ii=(iy-1)*nx
         do j=1,lbx
            xx=(99.-lbx+j)/100.
            y(ii+j)=xx*y(ii+j)
         end do
         ii=ii+nx-lbx
         do j=1,lbx
            xx=(100.-j)/100.
            y(ii+j)=xx*y(ii+j)
         end do
      end do
c   ***   end x boundary attenuation
c   ***   attenuate y boundary region
      if(ny .gt. 30)then
      do ix=1,nx
         ii=ix
         do j=1,lby
            xx=(99.-lby+j)/100.
            y(ii+(j-1)*nx)=xx*y(ii+(j-1)*nx)
         end do
         ii=(ny-lby)*nx+ix
         do j=1,lby
            xx=(100.-j)/100
            y(ii+(j-1)*nx)=xx*y(ii+(j-1)*nx)
         end do
      end do
      end if
      return
      end
