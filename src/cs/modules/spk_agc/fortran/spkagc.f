c-----------------------------------------------------------------------
c     SPKAGC - automatic gain control (running mean absolute amplitude)
c     Rotina AGC do programa xagc.f (SEISPAK), sem alteracoes no
c     algoritmo; apenas renomeada de AGC para SPKAGC para evitar
c     conflito de simbolos com outras bibliotecas.
c
c     s       (in/out) amostras do traco (nt)
c     buf1    (scratch) vetor de trabalho (nt)
c     nt      numero de amostras
c     sr      intervalo de amostragem [ms]
c     windowl comprimento da janela [ms]
c     valmn   limiar: amostras com |s| <= valmn no inicio/fim sao
c             ignoradas; tambem usado como epsilon no denominador
c     ampmax  amplitude media de saida
c-----------------------------------------------------------------------
      subroutine spkagc(s,buf1,nt,sr,windowl,valmn,ampmax)
      real s(*),buf1(*)
      LL=.5*WINDOWL/SR
      eps=valmn
      smax=0.
      s(1)=0.
      s(2)=0.
      s(nt-1)=0.
      s(nt)=0.
      do jj=1,nt
      if(abs(s(jj)).gt.valmn)goto 123
      end do
123   continue
      i1=jj
      do jj=nt,1,-1
      if(abs(s(jj)).gt.valmn)goto 124
      end do
124   continue
      i2=jj
       do jj=1,nt
       buf1(jj)=0.
       end do
      if(i2 .gt. i1)then
      do jj=i2,i1,-1
      j1=max0(i1,jj-ll)
      j2=min0(i2,jj+ll)
      if(jj.eq.i2)then
      do j=j1+1,j2
      smax=smax+abs(s(j))
      end do
      end if
      if(j2.lt.i2)smax=smax-abs(s(j2+1))
      if(j1.gt.i1 )smax=smax+abs(s(j1))
      avg=smax/(j2-j1+1)
      buf1(jj)=ampmax*s(jj)/(avg+eps)
      end do
      end if
      do  jj=1,nt
      s(jj)=buf1(jj)
      end do
       return
       end
