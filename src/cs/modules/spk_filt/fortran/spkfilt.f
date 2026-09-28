c-----------------------------------------------------------------------
c     Filtro trapezoidal do programa xxfilt.f (SEISPAK), separado em
c     rotinas para uso traco a traco no SeaSeis.
c
c     O desenho do filtro (ifrq1..ifrq4, rampas lineares, fator 1/nt1)
c     e o tamanho da FFT (fac235) sao os mesmos do xxfilt.f.
c     A aplicacao usa o caminho de "1 traco" do xxfilt.f: o filtro so
c     e aplicado nas frequencias positivas e o resultado e 2*real().
c     O dado complexo e guardado como pares (real,imag) num vetor real,
c     para chamar CFFTF/CFFTB (fftpack.f) sem mismatch de tipo.
c-----------------------------------------------------------------------
c
c     SPKFNT: tamanho da FFT (menor par >= nt com fatores 2, 3 e 5)
c
      subroutine spkfnt(nt,nt1)
      integer ifax(5)
      call fac235(nt,ifax)
      nt1=ifax(1)
      return
      end
c
c     SPKFDES: desenho do filtro trapezoidal e inicializacao da FFT
c       nt        amostras do traco
c       sr        intervalo de amostragem [ms]
c       frq1..4   frequencias do trapezio [Hz]
c       nt1       tamanho da FFT (de SPKFNT)
c       filt(nt1) filtro (saida, ja multiplicado por 1/nt1)
c       wsave     vetor de trabalho da FFT (4*nt1+15)
c       ifrq(4)   indices das frequencias (saida)
c       ierr      0 = ok, 1 = frequencias invalidas para esta FFT
c
      subroutine spkfdes(nt,sr,frq1,frq2,frq3,frq4,nt1,filt,wsave,
     :                   ifrq,ierr)
      real filt(*),wsave(*)
      integer ifrq(4)
      nt2=nt1/2+1
      dt=.001*sr
      ifrq1=nint(2.*dt*nt2*frq1)+1
      ifrq2=nint(2.*dt*nt2*frq2)+1
      ifrq3=nint(2.*dt*nt2*frq3)+1
      ifrq4=nint(2.*dt*nt2*frq4)+1
      ifrq(1)=ifrq1
      ifrq(2)=ifrq2
      ifrq(3)=ifrq3
      ifrq(4)=ifrq4
      ierr=0
      if(ifrq2.le.ifrq1 .or. ifrq3.lt.ifrq2 .or. ifrq4.le.ifrq3
     :   .or. ifrq4.gt.nt2)then
        ierr=1
        return
      endif
      do k=1,ifrq1
      filt(k)=0
      enddo
      den=ifrq2-ifrq1
      do k=ifrq1,ifrq2
      xnum=k-ifrq1
      filt(k)=xnum/den
      enddo
      do k=ifrq2,ifrq3
      filt(k)=1.
      enddo
      den=ifrq4-ifrq3
      do k=ifrq3,ifrq4
      xnum=ifrq4-k
      filt(k)=xnum/den
      enddo
      do k=ifrq4,nt1
      filt(k)=0.
      enddo
      factor=1./nt1
      do n=1,nt1
      filt(n)=factor*filt(n)
      enddo
      call cffti(nt1,wsave)
      return
      end
c
c     SPKFAPP: aplica o filtro a um traco
c       s(nt)     traco (entrada/saida)
c       a(2*nt1)  vetor de trabalho (complexo como pares real,imag)
c
      subroutine spkfapp(s,nt,nt1,filt,a,wsave)
      real s(*),filt(*),a(*),wsave(*)
      do n=1,nt1
      a(2*n-1)=0.
      a(2*n)=0.
      enddo
      do n=1,nt
      a(2*n-1)=s(n)
      enddo
      call cfftf(nt1,a,wsave)
      do n=1,nt1
      a(2*n-1)=filt(n)*a(2*n-1)
      a(2*n)  =filt(n)*a(2*n)
      enddo
      call cfftb(nt1,a,wsave)
c     double scale when processing 1 trace (como no xxfilt.f)
      do n=1,nt
      s(n)=2.*a(2*n-1)
      enddo
      return
      end
