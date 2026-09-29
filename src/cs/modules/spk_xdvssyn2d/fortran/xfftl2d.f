c-----------------------------------------------------------------------
c     xfftl2d.f (SEISPAK ice-pak), sem alteracao: FFT 2D com FFTPACK.
c     Usada por rhsfft (APX 1) no modulo SPK_XDVSSYN2D.
c-----------------------------------------------------------------------
      subroutine xfftl2d(y,a,nx,ny,isgn)
      complex y(*),a(*)
c      call fftc(y,a,2,2*nx,nx,ny,-isgn)
c      call fftc(y,a,2*nx,2,ny,nx,-isgn)
      if(isgn.eq.1)then
      call cffti(nx,a)
      do k=1,ny
      ii=(k-1)*nx+1
      call cfftf(nx,y(ii),a)
      enddo
      call ctran(y,a,nx,ny)
      call cffti(ny,y)
      do k=1,nx
      ii=(k-1)*ny+1
      call cfftf(ny,a(ii),y)
      enddo
      call ctran(a,y,ny,nx)
      endif
      if(isgn.eq.-1)then
      call cffti(nx,a)
      do k=1,ny
      ii=(k-1)*nx+1
      call cfftb(nx,y(ii),a)
      enddo
      call ctran(y,a,nx,ny)
      call cffti(ny,y)
      do k=1,nx
      ii=(k-1)*ny+1
      call cfftb(ny,a(ii),y)
      enddo
      call ctran(a,y,ny,nx)
      endif
      scale=1./(real(nx)*real(ny))
      if(isgn.eq.-1)then
         do 1 j=1,nx*ny
             y(j)=scale*y(j)
    1    continue
      end if
      return
      end
