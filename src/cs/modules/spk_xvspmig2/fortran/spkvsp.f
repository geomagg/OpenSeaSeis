c-----------------------------------------------------------------------
c     Migracao pre-empilhamento em profundidade por registro do programa
c     XVSPMIG2 (SEISPAK, "migc3"/SHOTMIG), para o modulo SPK_XVSPMIG2.
c
c     Cada registro: campo da fonte (posicao ALPH1, profundidade GZ,
c     com imagem-espelho acima da superficie: amplitude negativa) e campo
c     registrado nos tracos (superficie), ambos continuados com PSPI de
c     multiplas velocidades (multpspi); imagem por correlacao
c     (xsgimrg2r), correcao de amplitude, mute por angulo e soma na secao
c     empilhada.
c
c     SPKVSP  : prepara tamanhos, velocidades (cacomt2) e memoria
c     SPKVSPW : laco 1001 (registros) e 2000 (profundidade) do programa
c               principal. Diferencas marcadas com SEASEIS:
c               - dado, velocidade e saidas em memoria (antes: ONLINE1,
c                 OFFLINE1/2, DVSVELS / cvevent)
c               - sem o modo de continuacao para baixo (APX 8, mig1/mig2)
c                 e sem a analise de coerencia (OFFLINE3/4/5, vtps, VAVG)
c               - geometria de cada registro (ALPH1, GZ, SPJN) vinda do
c                 chamador (headers ou parametros)
c     SPKVXT2 : lista VXT -> funcoes densas (VEVENTS)
c-----------------------------------------------------------------------
      SUBROUTINE SPKVSP(D,NT0,NNTR,NRECS0,SR,AFIT,ALPHR,GZR,SPJNR,
     :                  NNZ0,ITZR0,DDZ0,DX0,APX0,FRQ10,FRQ20,OTY,
     :                  AMUTE,LTM0,LEV10,STK,NSTK,OUTR,IOUTR,INFO,IERR)
c     D(NT0,NNTR,NRECS0)  registros de entrada
c     AFIT                funcoes de velocidade: afit(1..3)=nxv,nyv,np;
c                         por posicao 2*np+3 palavras: iy, ix, (v,z)...
c     ALPHR,GZR,SPJNR     por registro: traco (real, 1=primeiro) da fonte,
c                         profundidade da fonte, deslocamento (em tracos)
c                         do traco 1 em relacao a origem da secao
c     STK(NNZ0+LTM0,NSTK) secao empilhada
c     OUTR(NNZ0+LTM0,NNTR,NRECS0) registros migrados (se IOUTR=1)
c     INFO(8)             nx(FFT), nt1, nlow, nhigh, nz2, nxv, MB, kshot max
c     IERR                0 ok, 1 banda vazia, 2 memoria, 3 fonte abaixo
c                         do modelo de velocidade
      REAL D(NT0,NNTR,*),AFIT(*),ALPHR(*),GZR(*),SPJNR(*),STK(*),OUTR(*)
      INTEGER INFO(8),IFAX(5)
      REAL, ALLOCATABLE :: WK(:),VVV(:),DAT(:),TRC(:)
      character*4 vint,rfc,lunin,lunout
      common/velneed/nrecs,irec,lev,vinc,ipass,npass,nxv,lag,locv,ispj
      common/migneed/nx,nt,nz,itzr,vint,dt,nnz,ddz,beta,dx,apx,maxscr,
     *rfc,alph1,gam1,frq,rfstir,itype,cyc,raa,padl,padr,zzy,ltm,ny,dy,
     *frq1,frq2,alph2,rni,alph3,
     *vmin,vmax,seps,lev1,lev2,ifor,ilor,window,lunin,lunout
      IERR=0
c   ***   parametros no common, como o programa principal
      NT=NT0
      NRECS=NRECS0
      NNZ=NNZ0
      ITZR=ITZR0
      DDZ=DDZ0
      DX=DX0
      DY=DX0
      NY=1
      APX=APX0
      FRQ1=FRQ10
      FRQ2=FRQ20
      BETA=OTY
      LTM=LTM0
      LEV1=LEV10
      LEV=LEV1
      VINT='no'
      RFC='YES'
      ITYPE=1
      CYC=0.
      RAA=0.
      GAM1=0.
      ALPH2=0.
      ALPH3=0.
      IFOR=0
      ILOR=-1
c     SEASEIS: RFSTIR=1 -> posicao i da grade de velocidade = DPN i+OTY
      RFSTIR=1.
      NPASS=1
      IPASS=1
      DT=.001*SR
c   ***   nx = tracos por registro (padding FFT para APX >= 6)
c     SEASEIS: o original so fazia o padding para APX 6..8, mas
c              multpspi (APX >= 11) usa fftc no comprimento nx
      NX=NNTR
      IF(APX.GE.6.)THEN
        CALL FAC235(NX,IFAX)
        NX=IFAX(1)
      ENDIF
      CALL FAC235(NT,IFAX)
      NT1=IFAX(1)
      NT2=NT1/2+1
      IFRQ1=2.*DT*NT2*FRQ1+1
      IFRQ2=2.*DT*NT2*FRQ2+1
      NL=MAX0(IFRQ1,1)
      NH=MIN0(NT2,IFRQ2)
      INFO(1)=NX
      INFO(2)=NT1
      INFO(3)=NL
      INFO(4)=NH
      IF(NH.LT.NL)THEN
        IERR=1
        RETURN
      ENDIF
      NWS=NH-NL+1
      NXWS=NX*NWS
c   ***   ZZY (ispj): deslocamento maximo por registro, para cobrir a
c         secao inteira com a grade de velocidade (nxv)
      SPMAX=0.
      DO I=1,NRECS
        SPMAX=AMAX1(SPMAX,SPJNR(I))
      ENDDO
      IF(NRECS.GT.1)THEN
        ZZY=AINT(SPMAX/(NRECS-1))+1.
      ELSE
        ZZY=1.
      ENDIF
      NZ=(NNZ-1)/ITZR+1
      LLLTM=LTM/ITZR
      NZ1=NZ+LLLTM+2
      NZ2=NZ+LLLTM
      NXVMAX=NX+(NRECS-1)*NINT(ZZY)
      INFO(5)=NZ2
      INFO(6)=NXVMAX
c   ***   fonte abaixo do modelo de velocidade?
      KMAX=0
      DO I=1,NRECS
        KMAX=MAX0(KMAX,INT(GZR(I)/(ITZR*DDZ)))
      ENDDO
      INFO(8)=KMAX
      IF(KMAX+1.GT.NZ1)THEN
        IERR=3
        RETURN
      ENDIF
c   ***   memoria
      IY=4*MAX0(NX,NWS)+1
      IZ=IY+NXWS
      IY1=IZ+NXWS
      IZ1=IY1+NXWS
      IA=IZ1+NXWS
      IB=IA+NXWS
      IC=IB+NXWS
      ID=IC+NXWS
      IS=2*(ID+4*NXWS)
      NWK=MAX0(2*(IY+4*NT1*NX),IS+8*NX)+100000
      NDAT=(NZ2+2)*ITZR*NX+NX*(NNZ+LTM)+1000
      NVV=NXVMAX*NZ1+1000
      INFO(7)=(4.*(REAL(NWK)+REAL(NDAT)+REAL(NVV)))/1.E6+1
      ALLOCATE(WK(NWK),VVV(NVV),DAT(NDAT),TRC(NT+10),STAT=IST)
      IF(IST.NE.0)THEN
        IERR=2
        RETURN
      ENDIF
c   ***   velocidades por passo (cacomt2, pares velocidade-profundidade)
      CALL CACOMT2(AFIT,VVV,WK)
      CALL SPKVSPW(WK,WK,VVV,DAT,TRC,D,NNTR,NT1,NL,NH,ALPHR,GZR,SPJNR,
     :             AMUTE,STK,NSTK,OUTR,IOUTR)
      DEALLOCATE(WK,VVV,DAT,TRC)
      RETURN
      END
c
      SUBROUTINE SPKVSPW(Y,S,IFITVV,DATA,VEL,D,NNTR,NT1,NL,NH,ALPHR,
     :                   GZR,SPJNR,AMUTE,STKDAT,NSTK,OUTR,IOUTR)
c     laco principal de xvspmig2.f; Y (complexo) e S (real) sao o mesmo
c     vetor, como a EQUIVALENCE (s(1),y(1)) do original
      COMPLEX Y(*)
      REAL S(*),IFITVV(*),DATA(*),VEL(*),D(NT,NNTR,*),ALPHR(*),GZR(*)
      REAL SPJNR(*),STKDAT(*),OUTR(*)
      character*4 vint,rfc,lunin,lunout
      common/velneed/nrecs,irec,lev,vinc,ipass,npass,nxv,lag,locv,ispj
      common/migneed/nx,nt,nz,itzr,vint,dt,nnz,ddz,beta,dx,apx,maxscr,
     *rfc,alph1,gam1,frq,rfstir,itype,cyc,raa,padl,padr,zzy,ltm,ny,dy,
     *frq1,frq2,alph2,rni,alph3,
     *vmin,vmax,seps,lev1,lev2,ifor,ilor,window,lunin,lunout
      sr=1000.*dt
      nz=(nnz-1)/itzr+1
      dz=itzr*ddz
      delw=2.*3.1415927/nt1
      dtt=dt
      con1=dz/dtt
      con2=(dtt/dx)**2
      nwdst=nx*itzr
      nws=nh-nl+1
      nxws=nx*nws
c   ***   compute y (complex) indicies
c   ***   s (real) indicies start at end of complex
      ishft=1
      iv=ishft+nx
      it=iv+nx
      it1=it+max0(nx,nws)
      it2=it1+nx
      is1=it2+nx
      is2=is1+nx
      iy=4*max0(nx,nws)+1
      iz=iy+nxws
      iy1=iz+nxws
      iz1=iy1+nxws
      ia=iz1+nxws
      ib=ia+nxws
      ic=ib+nxws
      id=ic+nxws
      is= 2*(id+4*nxws)
      nnzt=nnz+ltm
c   ***   initialize stacked section
      do 220 j=1,nnzt*nstk
220   stkdat(j)=0.
      do 1001 i=1,nrecs
      irec=i
      call clear(y(iy),4*nt1*nx)
c     SEASEIS: geometria do registro vinda do chamador
      alph1=alphr(i)
      gz2=gzr(i)
      spjn=spjnr(i)
      do 1 k=1,nntr
c   ***   fft in t direction
c   ***   use negative amplitude for mirror migration
      factor=-2./nt1
      i2=iz+k-1
      do 2 j=1,nt
      y(i2+(j-1)*nx)= factor*d(j,k,i)
2     continue
1     continue
      call fftc(y(iz),y(iz+nt1*nx),2*nx,2,nt1,nx,-1)
      i3=iz+(nl-1)*nx
      call movlev(y(i3),y(iy),2*nx*nws)
c     SEASEIS: zera o campo da fonte (antes ficava com restos da FFT no
c              primeiro passo)
      call clear(y(iy1),2*nxws)
      locv=1+spjn+.5
      locv=min0(locv,nxv-nx+1)
      kshot=-gz2/dz
      kstart=2
      llltm=ltm/itzr
      nz2=nz+llltm
      iuu=1
      ialph1=alph1
      do 2000 k=1,nz2-kshot
c   ***   get velocities for this level
      if(k+kshot .le. 0)then
      k0=1-kshot-k
      else
      k0=k+kshot
      end if
      kk=(k0-1)*nxv+locv
      do j=1,nx
      s(ishft+j-1)=con1/ifitvv(kk+j-1)
      s(iv+j-1)=con2*ifitvv(kk+j-1)**2
      end do
      sgn=1.
c   ***   compute arrival time for shot migration.
c   ***   for vsp migration shot is fired at depth level kshot
      xkk=gz2/dz +kshot+k
      call arriv6(s(it1),s(it2),s(iv),s(ishft),s(is),s(is1),s(is2),
     * s(it),jltm,lltm,nx,xkk)
c   ***   compute initial field for downward propagation
       if(k.eq.kstart)then
       z=(k-1)*dz
c   ***   compute amplitude factors
       j1=max0(1,ialph1-16)
       j2=max0(1,ialph1-8)
       j4=min0(nx,ialph1+16)
       j3=min0(nx,ialph1+8)
       den1=j2-j1
       den2=j4-j3
       do j=1,nx
       x=(alph1-j)*dx
       r=sqrt(2*dz**2+x**2)
       if( 1.le.j .and. j.lt.j1)amp=0.
       if(j1.le.j .and. j.lt.j2)amp=(j-j1)/den1
       if(j2.le.j .and. j.lt.j3)amp=1.
       if(j3.le.j .and. j.lt.j4)amp=(j4-j)/den2
       if(j4.le.j .and. j.le.nx)amp=0.
       s(is+j-1)=1000.*amp/r
       end do
       do iw=1,nws
       w=(nl+iw-2)*delw
       iiy1=(iw-1)*nx+iy1
       do j=1,nx
       y(iiy1+j-1)=s(is+j-1)*cexp(cmplx(0.,-w*s(it1+j-1)))
       end do
       end do
       endif
c   ***   write complex field to working array
      do  j=0,nxws-1
      y(iz+j)=y(iy+j)
      y(iz1+j)=y(iy1+j)
      end do
       nvels=nint(apx-10.)
c   ***   downward continue downgoing field
       if(k.ge.kstart)then
       call multpspi(y(iy1),y(ia),y(ib),y(ic),s(is),nl,nh,nt1,nx,
     * s(ishft),s(iv),-sgn,nvels)
        endif
c   ***   downward continue upcomming field
        if(k+kshot .gt. 0)then
      call multpspi(y(iy),y(ia),y(ib),y(ic),s(is),nl,nh,nt1,nx,
     * s(ishft),s(iv),sgn,nvels)
        endif
c   ***   apply koslov absorbing boundary for pspi
      lb=10
      do ii=1,nws
      do j=1,lb
      xx=(99-lb+j)/100.
      y(iy+(ii-1)*nx+j-1)=y(iy+(ii-1)*nx+j-1)*xx
      y(iy1+(ii-1)*nx+j-1)=y(iy1+(ii-1)*nx+j-1)*xx
      end do
      do j=1,lb
      xx=(100-j)/100.
      y(iy+ii*nx-lb-1+j)=y(iy+ii*nx-lb-1+j)*xx
      y(iy1+ii*nx-lb-1+j)=y(iy1+ii*nx-lb-1+j)*xx
      end do
      end do
c   ***   image: multiply upcomming field by conjugate of downgoing field
        if(k+kshot.gt.0)then
        call xsgimrg2r(data(iuu),y(iz),y(iz1),y(iy),y(iy1),s(ishft),
     :  nl,nh,nt1,nx,itzr)
c   ***   apply amplitude restoration
        do jj=1,itzr
        i1=iuu+(jj-1)*nx-1
        z=(k-1)*dz+(jj-1)*ddz
        do j=1,nx
        x=(alph1-j)*dx
        r=sqrt(z**2+x**2)
        data(i1+j)=-r*data(i1+j)
        end do
        end do
      iuu=iuu+nwdst
        endif
2000   continue
c   ***   write out nnz+ltm depth samples; apply mute (angle AMUTE)
      call rtran(data,s,nx,nnzt)
      call movlev(s,data,nntr*nnzt)
c   ***   accumulate data to stacked section
      iws=aint(spjn+.5)*nnzt
      rng1=amax1( abs( (alph1-1)*dx),abs( (alph1-nntr)*dx ) )
      iisr=1
      do 210 j=1,nntr
      x=abs((alph1-j)*dx)
      ii=x*tan(3.14159*amute/180.)/ddz
      ii=(x/rng1)*ii
      ii=min0(nnzt,max0(ii,iisr))
      j1=(j-1)*nnzt+1
      j2=j1+nnzt-1
      call clear(data(j1),ii)
      do 210 jj=j1,j2
      stkdat(iws+jj)=stkdat(iws+jj)+data(jj)
210   continue
c     SEASEIS: registro migrado devolvido em OUTR (antes: OFFLINE1)
      if(ioutr.eq.1)then
        call movlev(data,outr(1+(i-1)*nntr*nnzt),nntr*nnzt)
      endif
1001   continue
      return
      end
c
      SUBROUTINE SPKVXT2(VEL,DX,NXV,NYV,NP,AFIT,EE,X,BFIT,XX)
c     Lista VXT (horizontes) -> funcoes de velocidade densas, com duas
c     chamadas de VEVENTS (como em SPK_DVSMIG2D).
      REAL VEL(*),AFIT(*),EE(*),X(*),BFIT(*),XX(*)
      CHARACTER*4 VINT
      NXX=0
      NY=1
      NT=0
      NWDS=50000
      CALL VEVENTS(EE,VEL,XX,X,NXX,NY,NT,VINT,NWDS,AFIT,BFIT,DX)
      CALL VEVENTS(EE,VEL,XX,X,NXX,NY,NT,VINT,NWDS,AFIT,BFIT,DX)
      NXV=NINT(AFIT(1))
      NYV=NINT(AFIT(2))
      NP=NINT(AFIT(3))
      RETURN
      END
