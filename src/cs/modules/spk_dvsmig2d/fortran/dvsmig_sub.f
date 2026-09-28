c-----------------------------------------------------------------------
c     Rotinas copiadas sem alteracao de dvsmig.f (SEISPAK), exceto onde
c     marcado com "SEASEIS:". Usadas pelo modulo SPK_DVSMIG2D.
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
            r1=.01
            r2=.04
            if (r .le.r1)then 
               b(ind1) = 0.
               b(ind1+1) = 0.
               b(ind2) = 0.
               b(ind2+1) = 0.
            else

            if (r1.lt.r.and.r.lt.r2)scl= (r-r1)/(r2-r1)
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
