c-----------------------------------------------------------------------
c     Rotinas copiadas de xdvssyn2drs.f (SEISPAK) sem alteracao, exceto
c     onde marcado com "SEASEIS:". Usadas pelo modulo SPK_XDVSSYN2D.
c     Compilar com -fno-automatic (variaveis locais estaticas, como g77).
c-----------------------------------------------------------------------
      subroutine absorb1(u1,u2,v,nz,nx,ny,dz,dx,dy,xl,fs,tap1,tap2)
      real u1(*),u2(*),v(*),b1(100),b2(100)
c     SEASEIS: itime num common para o driver reiniciar a cada execucao
      common/spkrsit/itime,itdzz
      save
c   ***   compute koslov absorbing boundary taper
      if(itime.eq.0)then 
      if(nint(xl).eq.0)then
      nn=20
      else
      nn=nint(xl)
      endif
      write(*,*)'absorbing boundary attenuation : ',nn,' traces'
      nxz=nx*nz
      nxy=nx*ny
      nn2=2*nn
      itime=itime+1
      x1=nn*tap1
      x2=nn*tap2
      dw=.5*3.1415927/nn
      write(*,*)'cos taper: ',1.-x1,'+',x1,'*cos(i*',dw,')'
c   
      do 4 i=1,nn 
c      b1(i)=1.-i*tap1
c      b2(i)=1.-i*tap2
      b1(i)=1.-x1+x1*cos(i*dw)
      b2(i)=1.-x2+x2*cos(i*dw)
      b1(nn2+1-i)=b1(i)
      b2(nn2+1-i)=b2(i)
4     continue
c   *** 
c   *** 
      write(*,*)'inside absorb : nx ny nz nxz nn ',nx,ny,nz,nxz,nn
      endif
c   ***   apply Koslov absorbing boundary on x boundarys
      do 25 k=1,ny
      do 25 l=1,nn
      ii1=(k-1)*nxz+(l-1)*nz
      ii3=(k-1)*nxz+(nx-l)*nz
      do 25 j=1,nz
      u1(ii1+j)=b1(nn+l)*u1(ii1+j)
      u1(ii3+j)=b1(nn+l)*u1(ii3+j)
      u2(ii1+j)=b1(nn+l)*u2(ii1+j)
      u2(ii3+j)=b1(nn+l)*u2(ii3+j)
25    continue
c   ***   apply Koslov absorbing boundary on y boundarys
       if(ny.eq.1)goto 1234
      do l=1,nn
      ii2=(l-1)*nxz
      ii4=(ny-l)*nxz
      do j=1,nxz
      u1(ii2+j)=b1(nn+l)*u1(ii2+j)
      u1(ii4+j)=b1(nn+l)*u1(ii4+j)
      u2(ii2+j)=b1(nn+l)*u2(ii2+j)
      u2(ii4+j)=b1(nn+l)*u2(ii4+j)
      end do
      end do
1234   continue
c   ***   apply Koslov absorbing boundary on z boundarys
      do k=1,ny
      do j=1,nx
      if(fs .eq. 0.)then
      ii1=(k-1)*nxz+(j-1)*nz
      ii3=(k-1)*nxz+j*nz-nn
      do  l=1,nn
      u1(ii1+l)=b1(nn+l)*u1(ii1+l)
      u1(ii3+l)=b2(l)*u1(ii3+l)
      u2(ii1+l)=b1(nn+l)*u2(ii1+l)
      u2(ii3+l)=b2(l)*u2(ii3+l)
      end do 
      else
      ii3=(k-1)*nxz+j*nz-nn
      do  l=1,nn     
      u1(ii3+l)=b2(l)*u1(ii3+l)
      u2(ii3+l)=b2(l)*u2(ii3+l)
      enddo
      endif
 111  continue
c   ***   kill data at bottom
      ii3=(k-1)*nxz+j*nz-nn
      do l=0,3
      u1(ii3+nn-l)=0.
      u2(ii3+nn-l)=0.
      enddo
      enddo
      enddo
      return
      end
      subroutine xrhsfft(u,uz,ux,s,v,dz,dx,nz,nx) 
      real u(*),uz(*),ux(*),v(*)
      real s(*)
      nxz=nx*nz
c   ***   compute Uxx
      call rtran(u,s,nz,nx)
      call xfftddz(s,ux,nx,dx,nz)
      call rtran(s,ux,nx,nz)
c   ***   compute Uzz
      do j=1,nxz
      uz(j)=u(j)
      enddo
      call xfftddz(uz,s,nz,dz,nx)
c   ***   multiply by vdt
      do  j=1,nxz
      uz(j)=v(j)*(ux(j)+uz(j))
      enddo
      return
      end
      subroutine xfftddz(u,scr,nz,dz,nx)
      complex scr(*)
      real u(*)
      scale=1./nz
      call cffti(nz,scr(nz+1))
      pi=3.1415927 
      dkz=2.*pi/(nz*dz) 
      do k=1,nx/2
      k1=(2*k-2)*nz
      k2=k1+nz
      do j=1,nz
      scr(j)=cmplx(u(k1+j),u(k2+j))
      enddo
      call cfftf(nz,scr,scr(nz+1))
      do j=1,nz
      akz=min0(j-1,nz-j+1)*dkz
      scr(j)= -akz*akz*scr(j)
      enddo
      call cfftb(nz,scr,scr(nz+1))
      do j=1,nz
      u(k1+j)=scale*real(scr(j))
      u(k2+j)=scale*aimag(scr(j))
      enddo
      enddo      
      end
      subroutine rhsfft(u,r,s,v,dz,dx,nz,nx) 
      real u(1),v(1),r(1)
      complex s(1)
      nxz=nx*nz
      do 500 j=1,nxz
      s(j)=u(j)
500   continue
      call xfftl2d(s,s(nxz+1),nz,nx,1)
c   ***   multiply s(kz,kx) by -(kx**2 + kz**2)
      pi=3.1415927
      dkx=2.*pi/(nx*dx) 
      dkz=2.*pi/(nz*dz) 
      do 70 i=1,nx
      do 80 j=1,nz
      akx=min0(i-1,nx-i+1)*dkx
      akz=min0(j-1,nz-j+1)*dkz
      s((i-1)*nz+j)=-(akx**2+akz**2)*s((i-1)*nz+j)
80    continue
70    continue
      call xfftl2d(s,s(nxz+1),nz,nx,-1)
      do 700 j=1,nxz
      r(j)=v(j)*s(j)
700   continue
      return
         end
      subroutine rhsfftr(u,r,s,rho,v,dz,dx,nz,nx) 
c   ***   compute (rho*v**2)*(dx+dz)*(1/rho)*(dx+dz) U
      real u(*),r(*),rho(*),v(*)
      complex s(*)
      nxz=nx*nz
      do j=1,nxz
      s(j)=u(j)
      enddo
      call fftl2d(s,s(nxz+1),nz,nx,1)
c   ***   multiply s(kz,kx) by  kx - ikz
      pi=3.1415927
      dkx=2.*pi/(nx*dx) 
      dkz=2.*pi/(nz*dz) 
      do i=1,nx
      do j=1,nz
      if(i.le.nx/2)then
      akx=(i-1)*dkx
      else
      akx=(i-1-nx)*dkx
      endif
      if(j.le.nz/2)then
      akz=(j-1)*dkz
      else
      akz=(j-1-nz)*dkz
      endif
c      s((i-1)*nz+j)=cmplx(akx,-akz)*s((i-1)*nz+j)
      s((i-1)*nz+j)= cmplx(0.,sqrt(akx*akx+akz*akz))*s((i-1)*nz+j)
      enddo
      enddo
      call fftl2d(s,s(nxz+1),nz,nx,-1)
c   ***   multiply by 1/rho
      do j=1,nxz
      s(j)=rho(j)*s(j)
      enddo
      call fftl2d(s,s(nxz+1),nz,nx,1)
      do i=1,nx
      do j=1,nz
      if(i.le.nx/2)then
      akx=(i-1)*dkx
      else
      akx=(i-1-nx)*dkx
      endif
      if(j.le.nz/2)then
      akz=(j-1)*dkz
      else
      akz=(j-1-nz)*dkz
      endif
c      s((i-1)*nz+j)=cmplx(-akx,-akz)*s((i-1)*nz+j)
      s((i-1)*nz+j)= cmplx(0.,sqrt(akx*akx+akz*akz))*s((i-1)*nz+j)
      enddo
      enddo

      call fftl2d(s,s(nxz+1),nz,nx,-1)
c   ***   multiply by rho*(v*dt)**2 stored in v
      do j=1,nxz
      r(j)=v(j)*real(s(j))
      enddo
      return
         end
      subroutine test2dr(u,v,rho,s,nz,nx,dz,dx)
c     SEASEIS: w e scr de 10000 para 40000 (permite nx, nz ate ~4800)
      real u(*),v(*),rho(*),s(*),w(40000),scr(40000)
      nxz=nx*nz
c  Dx(rho*v**2)Dx(u)
      do j=1,nz,2
      do k=1,nx
      k1=(k-1)*nz
      w(k)=u(k1+j)
      w(nx+k)=u(k1+j+1)
      enddo
      call xfftdz(w,scr,nx,dx,2)
      do k=1,nx
      k1=(k-1)*nz
      w(k)=   v(k1+j  )*w(k)
      w(nx+k)=v(k1+j+1)*w(nx+k)
      enddo
      call xfftdz(w,scr,nx,dx,2)
      do k=1,nx
      k1=(k-1)*nz
      s(k1+j)=w(k)
      s(k1+j+1)=w(nx+k)
      enddo
      enddo
c  Dz(rho*v**2)Dz(u)
      do k=1,nx,2
      k1=(k-1)*nz+1
      call xfftdz(u(k1),scr,nz,dz,2)
      do j=1,2*nz
      u(k1+j-1)=v(k1+j-1)*u(k1+j-1)
      enddo
      call xfftdz(u(k1),scr,nz,dz,2)
      enddo
c  multiply by 1./rho (stored in rho)
      do j=1,nxz
      u(j)=rho(j)*(u(j)+s(j))
      enddo
c   finished
      end
      subroutine xfftdz(u,scr,nz,dz,nx)
      complex scr(*)
      real u(*)
      data nn/0/
      save nn
      scale=1./nz
      if(nn.ne.nz)then
      call cffti(nz,scr(nz+1))
      nn=nz
      endif
      pi=3.1415927 
      dkz=2.*pi/(nz*dz) 
      do k=1,nx/2
      k1=(2*k-2)*nz
      k2=k1+nz
      do j=1,nz
      scr(j)=cmplx(u(k1+j),u(k2+j))
      enddo
      call cfftf(nz,scr,scr(nz+1))
      do j=1,nz
      if(j.le.nz/2)then
      akz=(j-1)*dkz
      else
      akz=(j-1-nz)*dkz
      endif
c      akz=min0(j-1,nz-j+1)*dkz
      scr(j)= cmplx(0.,akz)*scr(j)
      enddo
      call cfftb(nz,scr,scr(nz+1))
      do j=1,nz
      u(k1+j)=scale*real(scr(j))
      u(k2+j)=scale*aimag(scr(j))
      enddo
      enddo      
      end
      subroutine xtest3drs(u,s,v,rho,nz,nx,ny,dz,dx,dy)
      real u(nz,nx,ny),v(nz,nx,ny),rho(nz,nx,ny)
c     SEASEIS: w e scr de 10000 para 40000 (permite nx, nz ate ~4800)
      real s(nz,nx,ny),w(40000),scr(40000)
c   ***   Dx( (1/rho)*Dx(u) )
      do iy=1,ny
      do iz=1,nz,2
      do ix=1,nx
      w(ix   )=u(iz  ,ix,iy)
      w(nx+ix)=u(iz+1,ix,iy)
      enddo
c     call xfftdz(w,scr,nx,dx,2)
      call fftdzsh(w,w,scr,dx,nx,2,1)
      do ix=1,nx
      ix1=min0(ix+1,nx)
      rho1=.5*(rho(iz,ix,iy)+rho(iz,ix1,iy))
      rho2=.5*(rho(iz+1,ix,iy)+rho(iz+1,ix1,iy))
      w(ix   )=rho1*w(ix)
      w(nx+ix)=rho2*w(nx+ix)
      enddo
c     call xfftdz(w,scr,nx,dx,2)
      call fftdzsh(w,w,scr,dx,nx,2,-1)
      do ix=1,nx
      s(iz  ,ix,iy)=w(ix)
      s(iz+1,ix,iy)=w(nx+ix)
      enddo
      enddo
      enddo
c   ***   Dy( (1/rho)*Dy(u) )
      if(ny.gt.1)then
      do ix=1,nx
      do iz=1,nz,2
      do iy=1,ny
      w(iy   )=u(iz  ,ix,iy)
      w(ny+iy)=u(iz+1,ix,iy)
      enddo
c     call xfftdz(w,scr,ny,dy,2)
      call fftdzsh(w,w,scr,dy,ny,2,1)
      do iy=1,ny
      iy1=min0(iy+1,ny)
      rho1=.5*(rho(iz,ix,iy)+rho(iz,ix,iy1))
      rho2=.5*(rho(iz+1,ix,iy)+rho(iz+1,ix,iy1))
      w(iy   )=rho1*w(iy   )
      w(ny+iy)=rho2*w(ny+iy)
      enddo
c     call xfftdz(w,scr,ny,dy,2)
      call fftdzsh(w,w,scr,dy,ny,2,-1)
      do iy=1,ny
      s(iz  ,ix,iy)=w(iy   )+s(iz  ,ix,iy)
      s(iz+1,ix,iy)=w(ny+iy)+s(iz+1,ix,iy)
      enddo
      enddo
      enddo
      endif
c   ***   Dz( (1/rho)*Dz(u) )
      do iy=1,ny
      do ix=1,nx,2
c     call xfftdz(u(1,ix,iy),scr,nz,dz,2)
      call fftdzsh(u(1,ix,iy),u(1,ix,iy),scr,dz,nz,2,1)
      do iz=1,nz
      iz1=min0(iz+1,nz)
      rho1=.5*(rho(iz,ix,iy)+rho(iz1,ix,iy))
      rho2=.5*(rho(iz,ix+1,iy)+rho(iz1,ix+1,iy))
c     SEASEIS: usa 1/rho no meio da celula (rho1, rho2), como em x e y;
c              o original calculava rho1/rho2 mas usava rho no no
      u(iz,ix  ,iy)=rho1*u(iz,ix  ,iy)
      u(iz,ix+1,iy)=rho2*u(iz,ix+1,iy)
      enddo
c     call xfftdz(u(1,ix,iy),scr,nz,dz,2)
      call fftdzsh(u(1,ix,iy),u(1,ix,iy),scr,dz,nz,2,-1)
      enddo
      enddo
c   ***   multiply by rho*(v*dt)**2
      do iy=1,ny
      do ix=1,nx
      do iz=1,nz
      u(iz,ix,iy)=v(iz,ix,iy)*(u(iz,ix,iy)+s(iz,ix,iy))
      enddo
      enddo
      enddo
c   finished
      end
      SUBROUTINE FFTDZSH(U,R,X,DZ,NZ,NX,ISHT)
      DIMENSION U(*),R(*),X(*)
      data nz0/0/
      save
C                                       ISHT= 1 MEANS SHIFT +HALF UNIT
C                                           =-1 MEANS SHIFT -HALF UNIT
      KZ   =NZ*2
      NN   =NZ/2-1
      CONST=2.*3.1415927/(DZ*NZ)
      NNP  =NN+1

      X(1)=0.
      X(2)=0.
      J =3
      JM=KZ-1
C                                   CONSTX=-CONST*DZ/2.*ISHT
      CONSTX= CONST*DZ/2.*ISHT
      COSX=COS(CONSTX)
      SINX=SIN(CONSTX)
      DO 40 I=1,NNP
        FACTOR=I*CONST
        IF (I.EQ.1) THEN
          COSNX=COSX
          SINNX=SINX
        ELSE
          CO=COSNX*COSX-SINNX*SINX
          SI=SINNX*COSX+COSNX*SINX
          COSNX=CO
          SINNX=SI
        ENDIF
          X(J)  =-SINNX*FACTOR
          X(J+1)= COSNX*FACTOR
        IF (I.LE.NN) THEN
          X(JM)  = X(J)
          X(JM+1)=-X(J+1)
        ENDIF
        J =J +2
        JM=JM-2
 40   CONTINUE
 60   CONTINUE
C
                      kz2=kz*2
      if(nz.ne.nz0)then
      call cffti(nz,x(kz2+1))
      nz0=nz
      endif
      DO 1 I=1,NX,2
        J=(I-1)*NZ+1
        CALL COMB(U(J),U(J+NZ),X(KZ+1),NZ, 1)
        call cfftf(nz,x(kz+1),x(kz2+1))
        CALL RFDERCC(X(KZ+1),X,KZ,1)
         call cfftb(nz,x(kz+1),x(kz2+1))
          call fsmpyhv(kz,1./real(nz),1,x,kz+1,x,kz+1)
        CALL COMB(R(J),R(J+NZ),X(KZ+1),NZ,-1)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE RFDERCC(R,FACT,KZ,ISHT)
      DIMENSION R(1),FACT(1)
      IF (ISHT.NE.0) THEN
        DO I=1,KZ,2
          D     =FACT(I)*R(I)  -FACT(I+1)*R(I+1)
          R(I+1)=FACT(I)*R(I+1)+FACT(I+1)*R(I)
          R(I)  =D
        ENDDO
      ELSE
        DO I=1,KZ,2
          D     =              -FACT(I+1)*R(I+1)
          R(I+1)=              +FACT(I+1)*R(I)
          R(I)  =D
        ENDDO
      ENDIF
      RETURN
      END
      SUBROUTINE COMB(A,B,C,NZ,ISW)
      DIMENSION A(1),B(1),C(1)
      IF (ISW.EQ.1) THEN
           J=1
        DO I=1,NZ
          C(J  )=A(I)
          C(J+1)=B(I)
          J=J+2
        ENDDO
      ELSEIF (ISW.EQ.-1) THEN
           J=1
        DO I=1,NZ
          A(I)=C(J  )
          B(I)=C(J+1)
          J=J+2
        ENDDO
      ENDIF
      RETURN
      END
      SUBROUTINE FSMPYHV(N,A,IA,B,IB,C,IC)
      DIMENSION A(1),B(1),C(1)
CC    HALF PRECISION B,C
      DO 1 I=1,N
      J=I-1
      C(IC+J)=A(IA)*B(IB+J)
    1 CONTINUE
      RETURN
      END
      subroutine xrhscon(u,uz,ux,s,v,dz,dx,nz,nx,n) 
      real u(*),uz(*),ux(*),v(*)
      real s(*)
      nxz=nx*nz
c   ***   compute Uxx
      call rtran(u,s,nz,nx)
      call xfftddz(s,ux,nx,dx,nz)
      call rtran(s,ux,nx,nz)
c   ***   compute Uzz
c     call xfftddz(uz,s,nz,dz,nx)
      do j=1,nxz
      s(j)=u(j)
      enddo
      call ddzcon2(s,uz,dz,nz,nx,n)
c     call dzz(s,uz,nz,nx,dz)
c   ***   multiply by vdt
      do  j=1,nxz
      uz(j)=v(j)*(ux(j)+uz(j))
      enddo
      return
      end
      subroutine ddzcon2(u,r,h,nz,nx,m)
      real u(1),r(1),s(nz+2*m)
      real b(7)
c     SEASEIS: itime (itdzz) num common para o driver reiniciar
      common/spkrsit/itabs,itime
      save
      itime=itime+1
      if(itime.eq.1)then
      id2=2
      write(*,*)'type ',id2,' second derivative filter length ',2*m-1
      call d2filt(b,id2,m)
      write(*,11)(b(j),j=1,m)
11    format(7f11.7)
      rh=1./h**2
      do j=1,m
      b(j)=rh*b(j)   
      enddo
      endif
      nzm=nz+m-1
      do 20 k=1,nx
      kk=(k-1)*nz
      m1=m-1
      do 10 j=1,m1
      s(j)=u(kk+nz+1-j)
      s(nzm+j)= u(kk+j)
10    continue
      do 15 j=1,nz
      s(m1+j)=u(kk+j)
15    continue
      if(m.eq.2)then
      do j=1,nz
      r(kk+j)=b(1)*s(m1+j)
     1    +b(2)*(s(m1+1+j)+s(m1-1+j))
      enddo
      elseif(m.eq.3)then
      do j=1,nz
      r(kk+j)=b(1)*s(m1+j)
     1    +b(2)*(s(m1+1+j)+s(m1-1+j))
     2    +b(3)*(s(m1+2+j)+s(m1-2+j))
       enddo
       elseif(m.eq.4)then
      do j=1,nz
      r(kk+j)=b(1)*s(m1+j)
     1    +b(2)*(s(m1+1+j)+s(m1-1+j))
     2    +b(3)*(s(m1+2+j)+s(m1-2+j))
     3    +b(4)*(s(m1+3+j)+s(m1-3+j))
       enddo
       elseif(m.eq.5)then
       do j=1,nz
      r(kk+j)=b(1)*s(m1+j)
     1    +b(2)*(s(m1+1+j)+s(m1-1+j))
     2    +b(3)*(s(m1+2+j)+s(m1-2+j))
     3    +b(4)*(s(m1+3+j)+s(m1-3+j))
     4    +b(5)*(s(m1+4+j)+s(m1-4+j))
      enddo
      elseif(m.eq.6)then
      do j=1,nz
      r(kk+j)=b(1)*s(m1+j)
     1    +b(2)*(s(m1+1+j)+s(m1-1+j))
     2    +b(3)*(s(m1+2+j)+s(m1-2+j))
     3    +b(4)*(s(m1+3+j)+s(m1-3+j))
     4    +b(5)*(s(m1+4+j)+s(m1-4+j))
     5    +b(6)*(s(m1+5+j)+s(m1-5+j))
      enddo
      elseif(m.eq.7)then
      do j=1,nz
      r(kk+j)=b(1)*s(m1+j)
     1    +b(2)*(s(m1+1+j)+s(m1-1+j))
     2    +b(3)*(s(m1+2+j)+s(m1-2+j))
     3    +b(4)*(s(m1+3+j)+s(m1-3+j))
     4    +b(5)*(s(m1+4+j)+s(m1-4+j))
     5    +b(6)*(s(m1+5+j)+s(m1-5+j))
     6    +b(7)*(s(m1+6+j)+s(m1-6+j))
      enddo
      endif
20    continue
      return
      end
      subroutine d2filt(b,itype,iord)
      real*8 alm(27,3)
      real*4 b(7)
c                                                 Taylor second derivatives
      ALM(1,1)=-2.982778
      ALM(2,1)= 1.714286
      ALM(3,1)=-0.267857
      ALM(4,1)= 0.052910
      ALM(5,1)=-0.008929
      ALM(6,1)= 0.001039
      ALM(7,1)=-0.000060
      ALM(8,1)=-2.927222
      ALM(9,1)= 1.666667
      ALM(10,1)=-0.238095
      ALM(11,1)= 0.039683
      ALM(12,1)=-0.004960
      ALM(13,1)= 0.000316
      ALM(14,1)=-2.847222
      ALM(15,1)= 1.600000
      ALM(16,1)=-0.200000
      ALM(17,1)= 0.025397
      ALM(18,1)=-0.001786
      ALM(19,1)=-2.722222
      ALM(20,1)= 1.500000
      ALM(21,1)=-0.150000
      ALM(22,1)= 0.011111
      ALM(23,1)=-2.500000
      ALM(24,1)= 1.333333
      ALM(25,1)=-0.083333
      ALM(26,1)=-2.0
      ALM(27,1)= 1.0
c						Holberg second derivatives
c						from Brian Sumner
      alm(1,2) = -3.0701709436264
      alm(2,2) = 1.7905029259120
      alm(3,2) = -0.31819137690324
      alm(4,2) = 7.7239746463727e-02
      alm(5,2) = -1.6861161667794e-02
      alm(6,2) = 2.3953380084827e-03
      alm(7,2) = 0.0
      alm(8,2) = -3.0701709436264
      alm(9,2) = 1.7905029259120
      alm(10,2) = -0.31819137690324
      alm(11,2) = 7.7239746463727e-02
      alm(12,2) = -1.6861161667794e-02
      alm(13,2) = 2.3953380084827e-03
      alm(14,2) = -2.9814394744996
      alm(15,2) = 1.7117485544487
      alm(16,2) = -0.26307621118030
      alm(17,2) = 4.7886001918499e-02
      alm(18,2) = -5.8386079371258e-03
      alm(19,2) = -2.8322153632841
      alm(20,2) = 1.5851861265830
      alm(21,2) = -0.18790397557714
      alm(22,2) = 1.8825530636174e-02
      alm(23,2) = -2.5597771680781
      alm(24,2) = 1.3741847787188
      alm(25,2) = -9.4296194679689e-02
      alm(26,2) = -2.0
      alm(27,2) =  1.0
c						Zhou et. al. 2nd deriv.
      alm(1,3) = -3.1378
      alm(2,3) =  1.8507
      alm(3,3) = -0.3643
      alm(4,3) = 0.1062
      alm(5,3) = -0.0313
      alm(6,3) = 0.0076
      alm(7,3) = 0.0
      alm(8,3) = -3.1378
      alm(9,3) =  1.8507
      alm(10,3) =  -0.3643
      alm(11,3) = 0.1062
      alm(12,3) = -0.0313
      alm(13,3) = 0.0076
      alm(14,3) = -3.0828
      alm(15,3) = 1.8068
      alm(16,3) = -0.3295
      alm(17,3) = 0.0830
      alm(18,3) = -0.0189
      alm(19,3) = -2.8826
      alm(20,3) = 1.6244
      alm(21,3) = -0.2109
      alm(22,3) = 0.0278
      alm(23,3) = -2.6796
      alm(24,3) = 1.4800
      alm(25,3) = -0.1402
      alm(26,3) = -2.1232
      alm(27,3) =  1.0616
       isum=0
       do i=2,iord
       isum=isum+i
       enddo
       index=28-isum
       do j=1,iord
       b(j)=alm(index+j-1,itype)
       enddo
       end
      subroutine fftfld2(u,s,vel,dz,dx,nz,nx) 
      real u(1),r(1)
      complex s(1)
      nxz=nx*nz
      do 500 j=1,nxz
      s(j)=u(j)
500   continue
      call fftl2d(s,u,nz,nx,1)
c   ***   apply frequency domain filter 
      pi=3.1415927
      dkx=2.*pi/(nx*dx) 
      dkz=2.*pi/(nz*dz) 
      do 50 i=1,nx
      do 60 j=1,nz
      ir=(i-1)*nz+j 
      akx=min0(i-1,nx-i+1)*dkx
      akz=min0(j-1,nz-j+1)*dkz
      aa=sqrt((akx*dx)**2+(akz*dz)**2)
      tl=3*pi/4.
      tr=pi
      rr=pi/(tr-tl)
      factor=.5*(1.+cos(rr*(aa-tl)))
      if(aa.le.tl)factor=1. 
      if(aa.ge.tr)factor=0. 
      s(ir)=factor*s(ir)
60    continue
50    continue
      call fftl2d(s,u,nz,nx,-1)
      do 600 j=1,nxz
      u(j)=s(j)
600   continue
      return
         end
      subroutine xyzfilt(a,b,ix0,iy0,iz0,nx,ny,nz)
      real a(*),b(*)
	iii=(iy0-1)*nx*nz+(ix0-1)*nz+1
	call xfilt1(a(iii),b,nz,1)
	write(*,*)'data filtered in z direction'
	ii=(iy0-1)*nx*nz
	do i=1,nz
	iii=ii+i
	call xfilt1(a(iii),b,nx,nz)
	end do
	write(*,*)'data filtered in x direction'
        if(ny .eq. 1)return
	do j=1,nx*nz
	iii=j
	call xfilt1(a(iii),b,ny,nx*nz)
	end do
	write(*,*)'data filtered in y direction'
	ii=(iy0-2)*nx*nz+(ix0-2)*nz+iz0-2
	do k=1,3
	ii1=ii+(k-1)*nx*nz
	do i=1,3
	iii=ii1+(i-1)*nz
	write(*,100)(int(.001*a(j)),j=iii+1,iii+3)
	end do
	write(*,*)
      end do
100   format(20i4)
	return
	end
      subroutine xfilt1(u,s,nx,inc)
      real u(*),s(*)
      real f(3)
      data f/.5,.25,0./
      s(1)=f(1)*u(1)+f(2)*(u(1+inc)-u(1))
      s(nx)=f(1)*u(1+(nx-1)*inc)+2.*f(2)*u(1+(nx-2)*inc)
      do i=2,nx-1
      s(i)=f(1)*u(1+(i-1)*inc)+f(2)*(u(1+(i-2)*inc)+u(1+i*inc))
      end do
      do i=1,nx
      u(1+(i-1)*inc)=s(i)
      end do
      return
      end
      SUBROUTINE SPGAUSSr(NX0,NZ0,AMP,U,DZ,DX,NZ,NX,ISWS,WF)
cc    call spgaussu(nx00,nz0,4.*scal,u(i1),dz,dx,nz,nx,isws)
C               
C         WF(4) IS NO OF PT AWAY FROM CENTER WHERE X=4 IN
C         EXP(-X**2/2.)=0.000335462
C
      DIMENSION WF(4)
      DIMENSION U(1)
      NXH=NX/2
      NZH=NZ/2
      NXZ=NX*NZ
      HFL=WF(4)
      print *,' spacial gauss with x=4 is ',hfl
      DI =4./HFL
c      CALL FSTOHV(NXZ*2,0.,1,U,NXZ*2+1)
      DO 1 IX=1,NX
         IXI   = (IX-NX0)
         IF     (IXI.LT.-NXH) THEN
                 IXI=IXI+NX
         ELSEIF (IXI.GE. NXH) THEN
                 IXI=IXI-NX
         ENDIF
         XI=IXI
         CONX =-(DI*XI)**2/2.
         TEMPX=-DI**2/DX*XI
         XI=.5+XI
         CONXU =-(DI*XI)**2/2.
         TEMPXU=-DI**2/DX*XI
      DO 1 IZ=1,NZ
         IZI   = (IZ-NZ0)
         IF     (IZI.LT.-NZH) THEN
                 IZI=IZI+NZ
         ELSEIF (IZI.GE. NZH) THEN
                 IZI=IZI-NZ
         ENDIF
         ZI=IZI
         CONZ =-(DI*ZI)**2/2.
         TEMPZ=-DI**2/DZ*ZI
         ZI=.5+ZI
         CONZU =-(DI*ZI)**2/2.

         TEMPZU=-DI**2/DZ*ZI
         IF (ISWS.EQ.1) THEN
            U(    (IX-1)*NZ+IZ)=TEMPX*EXP(CONX+CONZ)
            U(NXZ+(IX-1)*NZ+IZ)=TEMPZ*EXP(CONX+CONZ)
         ELSEIF (ISWS.EQ.2) THEN
            U(    (IX-1)*NZ+IZ)=TEMPZ*EXP(CONX+CONZ)
            U(NXZ+(IX-1)*NZ+IZ)=-TEMPX*EXP(CONX+CONZ)
         ELSEIF (ISWS.EQ.3) THEN
            U(    (IX-1)*NZ+IZ)=(TEMPX+TEMPZ)*EXP(CONX+CONZ)
            U(NXZ+(IX-1)*NZ+IZ)=(TEMPZ-TEMPX)*EXP(CONX+CONZ)
         ELSEIF (ISWS.EQ.-1) THEN
            U(    (IX-1)*NZ+IZ)=2.*TEMPX*EXP(2.*(CONX+CONZ))
            U(NXZ+(IX-1)*NZ+IZ)=2.*TEMPZ*EXP(2.*(CONX+CONZ))
         ELSEIF (ISWS.EQ.-2) THEN
            U(    (IX-1)*NZ+IZ)=2.*TEMPZ*EXP(2.*(CONX+CONZ))
            U(NXZ+(IX-1)*NZ+IZ)=-2.*TEMPX*EXP(2.*(CONX+CONZ))
         ELSEIF (ISWS.EQ.-3) THEN
            U(    (IX-1)*NZ+IZ)=2.*(TEMPX+TEMPZ)*EXP(2.*(CONX+CONZ))
            U(NXZ+(IX-1)*NZ+IZ)=2.*(TEMPZ-TEMPX)*EXP(2.*(CONX+CONZ))
         ENDIF
 1    CONTINUE
      if(nz0.eq.1)then
      do k=1,nx
      do j=0,10
      u((k-1)*nz+nz-j)=-u((k-1)*nz+j+2)
      u(nxz+(k-1)*nz+nz-j)=-u(nxz+(k-1)*nz+j+2)
      enddo
      enddo
      endif

      CALL FSMPYHV(NXZ*2,AMP,1,U,1,U,1)
 
      RETURN
      END
      SUBROUTINE SPGAUSSu(NX0,NZ0,AMP,U,DZ,DX,NZ,NX,ISWS,WF)
cc    call spgaussu(nx00,nz0,4.*scal,u(i1),dz,dx,nz,nx,isws)
C               
C         WF(4) IS NO OF PT AWAY FROM CENTER WHERE X=4 IN
C         EXP(-X**2/2.)=0.000335462
C
      DIMENSION WF(4)
      DIMENSION U(1)
      NXH=NX/2
      NZH=NZ/2
      NXZ=NX*NZ
      HFL=WF(4)
      print *,' spacial gauss with x=4 is ',hfl
      DI =4./HFL
c      CALL FSTOHV(NXZ*2,0.,1,U,NXZ*2+1)
      DO 1 IX=1,NX
         IXI   = (IX-NX0)
         IF     (IXI.LT.-NXH) THEN
                 IXI=IXI+NX
         ELSEIF (IXI.GE. NXH) THEN
                 IXI=IXI-NX
         ENDIF
         XI=IXI
         CONX =-(DI*XI)**2/2.
         TEMPX=-DI**2/DX*XI
         XI=.5+XI
         CONXU =-(DI*XI)**2/2.
         TEMPXU=-DI**2/DX*XI
      DO 1 IZ=1,NZ
         IZI   = (IZ-NZ0)
         IF     (IZI.LT.-NZH) THEN
                 IZI=IZI+NZ
         ELSEIF (IZI.GE. NZH) THEN
                 IZI=IZI-NZ
         ENDIF
         ZI=IZI
         CONZ =-(DI*ZI)**2/2.
         TEMPZ=-DI**2/DZ*ZI
         ZI=.5+ZI
         CONZU =-(DI*ZI)**2/2.
         TEMPZU=-DI**2/DZ*ZI
         IF (ISWS.EQ.1) THEN
            U(    (IX-1)*NZ+IZ)=TEMPXU*EXP(CONXU+CONZU)
            U(NXZ+(IX-1)*NZ+IZ)=TEMPZ*EXP(CONX+CONZ)
         ELSEIF (ISWS.EQ.2) THEN
            U(    (IX-1)*NZ+IZ)=TEMPZU*EXP(CONXU+CONZU)
            U(NXZ+(IX-1)*NZ+IZ)=-TEMPX*EXP(CONX+CONZ)
         ELSEIF (ISWS.EQ.3) THEN
            U(    (IX-1)*NZ+IZ)=(TEMPXU+TEMPZU)*EXP(CONXU+CONZU)
            U(NXZ+(IX-1)*NZ+IZ)=(TEMPZ-TEMPX)*EXP(CONX+CONZ)
         ELSEIF (ISWS.EQ.-1) THEN
            U(    (IX-1)*NZ+IZ)=2.*TEMPXU*EXP(2.*(CONXU+CONZU))
            U(NXZ+(IX-1)*NZ+IZ)=2.*TEMPZ*EXP(2.*(CONX+CONZ))
         ELSEIF (ISWS.EQ.-2) THEN
            U(    (IX-1)*NZ+IZ)=2.*TEMPZU*EXP(2.*(CONXU+CONZU))
            U(NXZ+(IX-1)*NZ+IZ)=-2.*TEMPX*EXP(2.*(CONX+CONZ))
         ELSEIF (ISWS.EQ.-3) THEN
            U(    (IX-1)*NZ+IZ)=2.*(TEMPXU+TEMPZU)*EXP(2.*(CONXU+CONZU))
            U(NXZ+(IX-1)*NZ+IZ)=2.*(TEMPZ-TEMPX)*EXP(2.*(CONX+CONZ))
         ENDIF
 1    CONTINUE
      if(nz0.eq.1)then
      do k=1,nx
      do j=0,10
      u((k-1)*nz+nz-j)=-u((k-1)*nz+j+2)
      u(nxz+(k-1)*nz+nz-j)=-u(nxz+(k-1)*nz+j+2)
      enddo
      enddo
      endif

      CALL FSMPYHV(NXZ*2,AMP,1,U,1,U,1)
 
      RETURN
      END
