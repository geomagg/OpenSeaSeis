c-----------------------------------------------------------------------
c     Migracao pre-empilhamento em profundidade por registro do programa
c     XVSPMIG1 (SEISPAK, "sgivspmig1"/SHOTMIG), para o modulo SPK_XVSPMIG1.
c
c     Registro de no (OBN) com o campo ASCENDENTE: pela reciprocidade os
c     tracos sao fontes na superficie (z = 0) e o no e a fonte pontual na
c     profundidade GZ. O campo registrado e continuado para baixo desde a
c     superficie; o campo da fonte comeca no passo KSHOT+2 abaixo do no
c     (KSHOT = GZ/DZ), ambos com PSPI de multiplas velocidades (multpspi).
c     A imagem (xsgimrg2r) so existe abaixo do no, sem imagem-espelho.
c     Correcao de amplitude r, mute por angulo a partir da profundidade
c     do no com rampa de 20 amostras, soma na secao empilhada.
c
c     SPKVSP1 : prepara tamanhos, velocidades (cacomt2) e memoria
c               (igual ao SPKVSP do SPK_XVSPMIG2)
c     SPKVSP1W: laco 1001 (registros) e 2000 (profundidade) de xvspmig1.f
c               (APX > 10). Diferencas marcadas com SEASEIS:
c               - dado, velocidade e saidas em memoria (antes: ONLINE1,
c                 OFFLINE1/2, DVSVELS / cvevent)
c               - so APX > 10 (multpspi); sem APX 4/5/6/8, sem a
c                 continuacao para baixo (mig1/mig2) e sem a analise de
c                 coerencia (OFFLINE3/4/5, vtps, VAVG)
c               - geometria de cada registro (ALPH1, GZ, SPJN) vinda do
c                 chamador (headers ou parametros)
c               - nx com padding FFT (fac235), como no SPK_XVSPMIG2
c     Rotinas (multpspi, arriv6, cacomt2, xsgimrg2r, fftc ...) sao as do
c     SPK_XVSPMIG2 (../spk_xvspmig2/fortran), iguais as do xvspmig1.f.
c-----------------------------------------------------------------------
      SUBROUTINE SPKVSP1(D,NT0,NNTR,NRECS0,SR,AFIT,ALPHR,GZR,SPJNR,
     :                  NNZ0,ITZR0,DDZ0,DX0,APX0,FRQ10,FRQ20,OTY,
     :                  AMUTE,LTM0,LEV10,STK,NSTK,OUTR,IOUTR,INFO,IERR,
     :                  IFOC,RFOC,TDAT,OTAU,OVC,OVAVG,PICKS)
c     D(NT0,NNTR,NRECS0)  registros de entrada (campo ascendente)
c     AFIT                funcoes de velocidade: afit(1..3)=nxv,nyv,np;
c                         por posicao 2*np+3 palavras: iy, ix, (v,z)...
c     ALPHR,GZR,SPJNR     por registro: traco (real, 1=primeiro) da fonte,
c                         profundidade da fonte, deslocamento (em tracos)
c                         do traco 1 em relacao a origem da secao
c     STK(NNZ0+LTM0,NSTK) secao empilhada
c     OUTR(NNZ0+LTM0,NNTR,NRECS0) registros migrados (se IOUTR=1)
c     INFO(8)             nx(FFT), nt1, nlow, nhigh, nz2, nxv, MB, kshot max
c     IERR                0 ok, 1 banda vazia, 2 memoria, 3 fonte abaixo
c                         da imagem (NNZ)
      REAL D(NT0,NNTR,*),AFIT(*),ALPHR(*),GZR(*),SPJNR(*),STK(*),OUTR(*)
      INTEGER INFO(8),IFAX(5)
c     SEASEIS: analise de coerencia (FOR LOR RNI WINDOW VMIN VMAX EPS)
c     IFOC: ifor ilor irni jnt (ilor < ifor: sem analise)
c     RFOC: vmin vmax eps ; TDAT(jnt,nz2,ngsn) ; saidas OTAU, OVC
c     (nt x jnt1 x ngsn), OVAVG (nt x ngsn), PICKS (201 x ngsn)
      INTEGER IFOC(*)
      REAL RFOC(*),TDAT(*),OTAU(*),OVC(*),OVAVG(*),PICKS(*)
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
      IFOR=IFOC(1)
      ILOR=IFOC(2)
      RNI=IFOC(3)
      WINDOW=(IFOC(4)-1)*SR
      VMIN=RFOC(1)
      VMAX=RFOC(2)
      SEPS=RFOC(3)
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
c     SEASEIS: a imagem comeca no passo KSHOT+2, abaixo do no
      IF(KMAX+2.GT.NZ2)THEN
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
      CALL SPKVSP1W(WK,WK,VVV,DAT,TRC,D,NNTR,NT1,NL,NH,ALPHR,GZR,SPJNR,
     :             AMUTE,STK,NSTK,OUTR,IOUTR,TDAT)
      DEALLOCATE(WK,VVV,DAT,TRC)
      IF(ILOR.GE.IFOR)THEN
        NGSN=(ILOR-IFOR)/NINT(RNI)+1
c       SEASEIS: VMIN/VMAX do common sao usados como rascunho na
c       migracao (como no original, que guardava vmn=vmin, vmx=vmax)
        CALL SPKV1FOC(TDAT,AFIT,NT,SR,NZ2,ITZR*DDZ,IFOR,NINT(RNI),
     :       IFOC(4),NGSN,RFOC(1),RFOC(2),RFOC(3),BETA,NNZ+LTM,OTAU,OVC,
     :       OVAVG,
     :       PICKS,1)
      ENDIF
      RETURN
      END
c
c
      SUBROUTINE SPKVSP1W(Y,S,IFITVV,DATA,VEL,D,NNTR,NT1,NL,NH,ALPHR,
     :                   GZR,SPJNR,AMUTE,STKDAT,NSTK,OUTR,IOUTR,TDAT)
c     laco principal de xvspmig1.f; Y (complexo) e S (real) sao o mesmo
c     vetor, como a EQUIVALENCE (s(1),y(1)) do original
      COMPLEX Y(*)
      REAL S(*),IFITVV(*),DATA(*),VEL(*),D(NT,NNTR,*),ALPHR(*),GZR(*)
      REAL SPJNR(*),STKDAT(*),OUTR(*),TDAT(*)
      COMPLEX, ALLOCATABLE :: WF(:),WB(:)
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
      nvels=nint(apx-10.)
c     SEASEIS: analise de coerencia; WF = y(ia), WB = y(ib)/y(ic)
      irni=nint(rni)
      twindow=window
      jnt=twindow/sr+1
      if(ilor.ge.ifor)then
      allocate(wf(22*nt1+10),wb(nws+12))
      ngsn=(ilor-ifor)/irni+1
      llltm=ltm/itzr
      do j=1,jnt*(nz+llltm)*ngsn
      tdat(j)=0.
      enddo
      endif
c   ***   initialize stacked section
      do 220 j=1,nnzt*nstk
220   stkdat(j)=0.
      do 1001 i=1,nrecs
      irec=i
      call clear(y(iy),4*nt1*nx)
c     SEASEIS: geometria do registro vinda do chamador
      alph1=alphr(i)
      gz2=abs(gzr(i))
      spjn=spjnr(i)
      do 1 k=1,nntr
c   ***   fft in t direction (no mirror: positive amplitude)
      factor=2./nt1
c.  fix factor
      factor=factor/nt1
c.  end factor fix
      i2=iz+k-1
      do 2 j=1,nt
      y(i2+(j-1)*nx)= factor*d(j,k,i)
2     continue
1     continue
      call fftc(y(iz),y(iz+nt1*nx),2*nx,2,nt1,nx,-1)
      i3=iz+(nl-1)*nx
      call movlev(y(i3),y(iy),2*nx*nws)
c     SEASEIS: zera o campo da fonte (antes ficava com restos da FFT)
      call clear(y(iy1),2*nxws)
      locv=1+spjn+.5
      locv=min0(locv,nxv-nx+1)
      kshot=abs(gz2/dz)
      llltm=ltm/itzr
      nz2=nz+llltm
      iuu=1
      kstart=kshot+2
      ialph1=alph1
      do 2000 k=1,nz2
      k0=k
c   ***   get velocities for this level
      kk=(k-1)*nxv+locv
      do 300 j=0,nx-1
      s(ishft+j)=con1/ifitvv(kk+j)
300   s(iv+j)=con2*ifitvv(kk+j)**2
      sgn=1.
c   ***   compute arrival time for shot migration.
c   ***   for vsp migration shot is fired at depth level kshot
      xkk=k-gz2/dz
      kshot=xkk
      if(kshot .ge. 1)then
      call arriv6(s(it1),s(it2),s(iv),s(ishft),s(is),s(is1),s(is2),
     * s(it),jltm,lltm,nx,xkk)
c   ***   compute initial field for downward propagation
       if(kstart .eq. k0)then
c   ***   compute amplitude factors
       z=xkk*dz
       j1=max0(1,ialph1-8)
       j2=max0(1,ialph1-4)
       j4=min0(nx,ialph1+8)
       j3=min0(nx,ialph1+4)
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
c   ***   process frequencies
       do iw=1,nws
       w=(nl+iw-2)*delw
       iiy1=(iw-1)*nx+iy1
       do j=1,nx
       y(iiy1+j-1)=s(is+j-1)*cexp(cmplx(0.,-w*s(it1+j-1)))
       end do
       end do
       end if
      end if
c   ***   write complex field to working array
      do  j=0,nxws-1
      y(iz+j)=y(iy+j)
      y(iz1+j)=y(iy1+j)
      end do
c   ***   above the node: continue only the recorded (upcoming) field
      if(k.lt.kstart)then
      call multpspi(y(iy),y(ia),y(ib),y(ic),s(is),nl,nh,nt1,nx,
     * s(ishft),s(iv),sgn,nvels)
       do j=1,nxws
       y(iz+j-1)=y(iy+j-1)
       end do
      else
c   ***   downward continue upcomming field
      call multpspi(y(iy),y(ia),y(ib),y(ic),s(is),nl,nh,nt1,nx,
     * s(ishft),s(iv),sgn,nvels)
c   ***   downward continue downgoing (source) field
      call multpspi(y(iy1),y(ia),y(ib),y(ic),s(is),nl,nh,nt1,nx,
     * s(ishft),s(iv),-sgn,nvels)
      end if
c   ***   apply koslov absorbing boundary for pspi
      lb=30
c     SEASEIS: lb limitado a nx/4 (registros curtos)
      lb=min0(lb,nx/4)
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
c   ***   image
      if(k .ge. kstart)then
c   ***   multiply upcomming field by conjugate of downgoing field
        call xsgimrg2r(data(iuu),y(iz),y(iz1),y(iy),y(iy1),s(ishft),
     :  nl,nh,nt1,nx,itzr)
c   ***   apply amplitude restoration
        do jj=1,itzr
        i1=iuu+(jj-1)*nx-1
        z=(k-1)*dz+(jj-1)*ddz
        do j=1,nx
        x=(alph1-j)*dx
        r=sqrt(z**2+x**2)
        data(i1+j)=r*data(i1+j)
        end do
        end do
c   ***   SEASEIS: analise de coerencia de xvspmig1.f: janela de
c   ***   twindow ms do campo continuado (y(iz)*conjg(y(iz1))) em cada
c   ***   posicao FOR..LOR (passo RNI) e cada passo de profundidade
      if(ilor.ge.ifor)then
      do j=1,nx
      igsn=j+spjn+.5
      if(igsn.ge.ifor.and.igsn.le.ilor.and.mod(igsn-ifor,irni).eq.0)then
      call clear(wf,2*nt1)
      ii=nl
      do jj=0,nws-1
      wf(ii+jj)=y(iz+j-1+jj*nx)*conjg(y(iz1+j-1+jj*nx))
      end do
c   *** shift data back in time by .5*twindow
      shf=-.5*twindow/sr
c   ***   shift complex data by this amount
      call vshift1(wf(ii),wb,wb,wb(nws+3),shf,delw,nl,nh,1,1.)
c   ***   inverse transform
      call fftc(wf,wf(nt1+1),2,2*nt1,nt1,1,+1)
c   ***   sum time data to appropriate gsn locations
      jgsn=(igsn-ifor)/irni +1
      itdat=((jgsn-1)*nz2+k0-1)*jnt+1
      do jj=1,jnt
      tdat(itdat+jj-1)=tdat(itdat+jj-1)+real(wf(jj))
      end do
      end if
      end do
      end if
      else
        call clear(data(iuu),nwdst)
      end if
      iuu=iuu+nwdst
2000   continue
c   ***   write out nnz+ltm depth samples to depth migrated file
      call rtran(data,s,nx,nnzt)
      call movlev(s,data,nntr*nnzt)
c     SEASEIS: registro migrado (sem mute) devolvido em OUTR (OFFLINE1)
      if(ioutr.eq.1)then
        call movlev(data,outr(1+(i-1)*nntr*nnzt),nntr*nnzt)
      endif
c   ***   accumulate data to stacked section
c   ***   mute below node depth + angle, 20-sample ramp
      iws=aint(spjn+.5)*nnzt
      rng1=amax1( abs( (alph1-1)*dx),abs( (alph1-nntr)*dx ) )
      iisr=1
      do 210 j=1,nntr
      x=abs((alph1-j)*dx)
      ii=x*tan(3.14159*amute/180.)/ddz
      ii=(gz2/ddz)+(x/rng1)*ii
      ii=min0(nnzt,max0(ii,iisr))
      j1=(j-1)*nnzt+1
      j2=j1+nnzt-1
      do 210 jj=j1,j2
      jjj=jj-j1-ii
      xx=.05*min0(20,max0(0,jjj))
      stkdat(iws+jj)=stkdat(iws+jj)+xx*data(jj)
210   continue
1001   continue
      return
      end
