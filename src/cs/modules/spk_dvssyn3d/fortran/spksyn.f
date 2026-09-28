c-----------------------------------------------------------------------
c     Sintetico de afastamento nulo 2D/3D do programa DVSSYN (SEISPAK)
c     (refletor explosivo), para o modulo SPK_DVSSYN3D.
c
c     SPKSYN    : prepara parametros (common /migneed/), memoria, chama
c                 SPKSYNMOD (modelo) e SPKSYNW (continuacao para cima)
c     SPKSYNW   : programa principal sync3d de dvssyn.f.
c                 Diferencas (SEASEIS:): modelo, velocidade e espectro em
c                 memoria (antes: arquivos indata, velscr, velscr1, scr,
c                 outdata e transposicoes em disco); todas as
c                 frequencias num unico passo (nstepw=1); saida em float
c                 (antes: 16 bits com pi2).
c     SPKSYNMOD : xmodeli2 de dvssyn.f com o modelo em memoria
c-----------------------------------------------------------------------
      SUBROUTINE SPKSYN(VEL,NVELW,AFIT,NAFIT,IVT,NXIN,NYIN,NT0,SR0,
     :                  NNZ0,ITZR0,DX0,DY0,DDZ0,APX0,ALPH10,OTY0,
     :                  FRQ10,FRQ20,IRFC,OUT,INFO,IERR)
c     VEL(NVELW)   lista VXT (IVT=1); tambem usado como rascunho
c     AFIT(NAFIT)  funcoes de velocidade (IVT=3: preenchido pelo chamador)
c     OUT(NT0,NXIN*NYIN) tracos sinteticos (y externo, x interno)
c     INFO(8)      nx, ny (FFT), nt1, nlow, nhigh, nz, memoria MB
c     IERR         0 ok, 1 banda vazia, 2 sem memoria, 3 NT*ALPH1 > 9000
      REAL VEL(*),AFIT(*),OUT(NT0,*)
      INTEGER INFO(8),IFAX(5)
      REAL, ALLOCATABLE :: WK(:),EE(:),XS(:),BFIT(:),REFL(:),VELZ(:)
      common/velneed/nrecs,irec,lev,vinc,ipass,npass,nxv,lag,locv,ispj
      character*4 vint,rfc
      common/migneed/nx,nt,nz,itzr,vint,dt,nnz,ddz,beta,dx,apx,maxscr,
     *rfc,alph1,gam1,frq,rfstir,itype,cyc,raa,padl,padr,zzy,ltm,ny,dy,
     *frq1,frq2,alph2,rni,alph3,
     *vmin,vmax,seps,lev1,lev2,ifor,ilor,window,lunin,lunout
      IERR=0
      NX=NXIN
      NY=NYIN
      IF(NY.EQ.0)NY=1
      NT=NT0
      DT=.001*SR0
      NNZ=NNZ0
      ITZR=ITZR0
      DX=DX0
      DY=DY0
      IF(DY.EQ.0.)DY=DX
      DDZ=DDZ0
      APX=APX0
      ALPH1=ALPH10
      BETA=OTY0
      FRQ1=FRQ10
      FRQ2=FRQ20
      RNI=1.
      RFC='YES'
      IF(IRFC.EQ.0)RFC='NO'
      IF(NINT(ALPH1)*NT.GT.9000)THEN
        IERR=3
        RETURN
      ENDIF
      NZ=NNZ/ITZR
c   ***   pad traces out to powers of 2,3,5
      NXSAVE=NX
      NYSAVE=NY
      IF(APX .GE. 6.)THEN
        CALL FAC235(NX,IFAX)
        NX=IFAX(1)
        IF(NY .GT. 1)THEN
          CALL FAC235(NY,IFAX)
          NY=IFAX(1)
        ENDIF
      ENDIF
      NXY=NX*NY
      CALL FAC235(NT,IFAX)
      NT1=IFAX(1)
      NT2=NT1/2+1
      IFRQ1=2.*DT*NT2*FRQ1+1
      IFRQ2=2.*DT*NT2*FRQ2+1
      NLOW=MAX0(IFRQ1,1)
      NHIGH=MIN0(NT2,IFRQ2)
      INFO(1)=NX
      INFO(2)=NY
      INFO(3)=NT1
      INFO(4)=NLOW
      INFO(5)=NHIGH
      INFO(6)=NZ
      IF(NHIGH.LT.NLOW)THEN
        IERR=1
        RETURN
      ENDIF
      NWS=NHIGH-NLOW+1
      NXYWS=NXY*NWS
      NWDST=NXY*ITZR
c                                    modelo (xmodeli2)
      IPRM=NINT(ALPH1)
      NSC=MAX0(500000,20*IPRM*(NT+NNZ)+(NZ+2)*ITZR*IPRM+1000)
      NW=2*(6*NXYWS+2*NXY)+10*NXY+2*NWDST+4*MAX0(NXY,NWS)+4*NT1+10000
      NBYT=4*(NW+2*NSC+NXY*(NNZ+NZ)+NAFIT)
      INFO(7)=NBYT/1000000+1
      ALLOCATE(EE(NSC),XS(NSC),BFIT(NAFIT),REFL(NXY*NNZ),
     :         VELZ(NXY*(NZ+1)),STAT=IST)
      IF(IST.NE.0)THEN
        IERR=2
        RETURN
      ENDIF
      EE=0.
      XS=0.
      BFIT=0.
      REFL=0.
      VELZ=0.
      CALL SPKSYNMOD(VEL,EE,XS,BFIT,AFIT,REFL,VELZ,IVT)
      DEALLOCATE(EE,XS,BFIT)
      ALLOCATE(WK(NW),STAT=IST)
      IF(IST.NE.0)THEN
        IERR=2
        DEALLOCATE(REFL,VELZ)
        RETURN
      ENDIF
      WK=0.
      CALL SPKSYNW(WK,WK,REFL,VELZ,OUT,NX,NY,NXSAVE,NYSAVE,NT,NT1,
     :             NZ,NNZ,ITZR,DT,DX,DY,DDZ,APX,RFC,NLOW,NHIGH)
      DEALLOCATE(WK,REFL,VELZ)
      RETURN
      END
c
      SUBROUTINE SPKSYNW(Y,S,REFL,VELZ,OUT,NX,NY,NXSAVE,NYSAVE,NT,NT1,
     :                   NZ,NNZ,ITZR,DT,DX,DY,DDZ,APX,RFC,NLOW,NHIGH)
c     Y (complexo) e S (real) sao o MESMO vetor, como a
c     EQUIVALENCE (y(1),s(1)) de dvssyn.f
      COMPLEX Y(*),A(20000),B(20000)
      REAL S(*),REFL(NNZ,*),VELZ(NZ,*),OUT(NT,*)
      CHARACTER*4 RFC
      NXY=NX*NY
      NWDST=NXY*ITZR
      DZ=ITZR*DDZ
      DELW=2.*3.1415927/NT1
      NT2=NT1/2+1
      NL=NLOW
      NH=NHIGH
      NWS=NH-NL+1
      NXYWS=NWS*NXY
      NWWS=NWS
c   ***   compute y (complex) indicies (como em dvssyn.f)
      IF(NY.GT.1 .AND. NX.GT.1)THEN
        IY=1
        IA=IY+NXYWS
        IB=IA+NXY
        IC=IB+NXY
        ID=IC+NXY
        IU=2*(ID+NXY)+1
        ISX=IU+NWDST
        ISY=ISX+NXY
        IVX=ISY+NXY
        IVY=IVX+NXY
        IS=IVY+NXY
      ENDIF
      IF(NX.EQ.1 .OR. NY.EQ.1)THEN
        IY=1
        IA=IY+NXYWS
        IB=IA+NXYWS
        IC=IB+NXYWS
        ID=IC+NXYWS
        ISX=2*(ID+NXYWS-1)+1
        IVX=ISX+NXY
        IVY=IVX+NXY
        ISY=IVY+NXY
        IT=ISY+NXY
        IS=IT+MAX0(NXY,NWS)
        IU=IS+NXY
      ENDIF
c   ***   initialize complex field to zero.
      DO 5 J=1,NXYWS
5     Y(J)=0.
      CON1=DZ/(.5*DT)
      CON2=(.5*DT/DY)**2
      CON3=(DY/DX)**2
      SGN=-1.
      XSHFT=DZ/(2500.*DT)
      XVV=((2500.*DT)/DX)**2
      IF(RFC.EQ.'NO')XSHFT=ITZR
      IF(APX.EQ.10.)
     *CALL XXCONV(Y(IY),Y(IA),Y(IB),Y(IC),S(IS),NL,NH,NT1,NX,XSHFT,
     * XVV,DX,DT,DZ,SGN)
c   ***  upward continue and add model, from the bottom
      DO 2000 K=1,NZ
c     SEASEIS: velocidade do nivel nz-k+1 vem da memoria (antes: gi2)
      DO 299 J=0,NXY-1
299   S(IVY+J)=VELZ(NZ-K+1,J+1)
      DO 300 J=0,NXY-1
      S(ISY+J)=CON1/S(IVY+J)
300   S(IVY+J)=CON2*S(IVY+J)**2
      CALL RTRAN(S(ISY),S(ISX),NX,NY)
      CALL RTRAN(S(IVY),S(IVX),NX,NY)
      DO 21 J=1,NXY
21    S(IVX+J-1)=CON3*S(IVX+J-1)
c     SEASEIS: refletividade do nivel vem da memoria (antes: gi2 'scr')
      DO 601 JJ=1,ITZR
      IREC=(NZ-K)*ITZR+JJ
      IUU=(JJ-1)*NXY+IU
      DO 602 J=1,NXY
602   S(IUU+J-1)=REFL(IREC,J)
601   CONTINUE
c   ***   write x-y-z data into x-y-w space
      CALL DEPTH2W(Y,S(IU),Y(IA),Y(IB),S(ISY),NL,NH,NT1,ITZR,NXY)
      IF(NY.GT.1 .AND. NX.GT.1)THEN
        DO 100 I=1,NWS
        W=(NL+I-2)*DELW
        I1=(I-1)*NXY+1
        IF(APX .EQ. 4. .OR. APX .EQ. 5.)THEN
          IF(W.GT.0.)THEN
            CALL CTRAN(Y(I1),Y(IA),NX,NY)
            CALL DWN3D(Y(IA),Y(IA),Y(I1),Y(IB),Y(IC),Y(ID),W,S(IVX),
     *                 S(ISX),NX,NY,SGN)
            CALL CTRAN(Y(IA),Y(I1),NY,NX)
            CALL DWN3D(Y(I1),Y(I1),Y(IA),Y(IB),Y(IC),Y(ID),W,S(IVY),
     *                 S(ISY),NY,NX,SGN)
          ENDIF
        ENDIF
        IF(APX.EQ.6. .OR. APX.EQ.7.)THEN
          CALL PHS2D(Y(I1),Y(IA),Y(IB),Y(IC),NX,NY,DX,DY,W,
     :               S(ISY),S(IVY),SGN)
        ENDIF
        IF(APX .GE. 11.)THEN
          NVELS=NINT(APX-10.)
          CALL MVPHS2D(Y(I1),Y(IA),Y(IB),Y(IC),NX,NY,DX,DY,W,
     :                 S(ISY),S(IVY),SGN,NVELS)
        ENDIF
        IF(APX .GE. 13.)THEN
          DO IIY=1,NY
            II=I1+(IIY-1)*NX-1
            DO J=1,20
              XX=(79.+J)/100.
              Y(II+J)=XX*Y(II+J)
            END DO
            II=II+NX-20
            DO J=1,20
              XX=(100.-J)/100.
              Y(II+J)=XX*Y(II+J)
            END DO
          END DO
          DO IIX=1,NX
            II=I1+IIX-1
            DO J=1,20
              XX=(79.+J)/100.
              Y(II+(J-1)*NX)=XX*Y(II+(J-1)*NX)
            END DO
            II=I1+(NY-20)*NX+IIX-1
            DO J=1,20
              XX=(100.-J)/100
              Y(II+(J-1)*NX)=XX*Y(II+(J-1)*NX)
            END DO
          END DO
        ENDIF
100     CONTINUE
        CALL VSHIFT1(Y,Y(IA),Y(IA),Y(IB),S(ISY),DELW,NL,NH,NXY,SGN)
      ENDIF
      IF(NY.EQ.1 .OR. NX.EQ.1)THEN
        IF(APX .EQ. 1.)THEN
          CALL VSHIFT1(Y,Y(IA),Y(IA),Y(IB),S(ISX),DELW,NL,NH,NXY,SGN)
        ENDIF
        IF(APX .EQ. 4. .OR. APX .EQ. 5.)THEN
          CALL CTRAN(Y(IY),Y(IA),NXY,NWS)
          CALL DWN45B(Y(IA),Y(IA),Y(IY),Y(IB),Y(IC),Y(ID),S(IT),S(IS),
     *                S(IVX),S(ISX),NXY,NL,NH,DELW,SGN)
          CALL CTRAN(Y(IA),Y(IY),NWS,NXY)
        ENDIF
        IF(APX .EQ. 6. .OR. APX .EQ. 7.)
     *  CALL PSPI(Y(IY),Y(IA),Y(IB),Y(IC),S(IS),NL,NH,NT1,NXY,S(ISX),
     *  S(IVX),SGN)
        IF(APX.EQ.10.)
     *  CALL XXCONV(Y(IY),Y(IA),Y(IB),Y(IC),S(IS),NL,NH,NT1,NXY,S(ISX),
     *  S(IVX),DX,DT,DZ,SGN)
        IF(APX .GE. 11.)THEN
          NVELS=NINT(APX-10.)
          CALL XMULTPSPI(Y(IY),Y(IA),Y(IB),Y(IC),S(IS),NL,NH,NT1,NX,
     *                   S(ISX),S(IVX),SGN,NVELS)
        ENDIF
        IF(APX .GE. 13.)THEN
          DO I=1,NWS
            DO J=1,20
              XX=(79+J)/100.
              Y(IY+(I-1)*NX+J-1)=Y(IY+(I-1)*NX+J-1)*XX
            END DO
            DO J=1,20
              XX=(100-J)/100.
              Y(IY+I*NX-21+J)=Y(IY+I*NX-21+J)*XX
            END DO
          END DO
        ENDIF
      ENDIF
2000  CONTINUE
c   ***   inverse transform and output (como no ramo ONLINE1 de dvssyn.f)
c     SEASEIS: espectro lido direto de Y (antes: outdata/xctrnsp2/gi2)
      DO 10 I=1,NYSAVE
      DO 10 K=1,NXSAVE
      ITR=(I-1)*NX+K
      DO 20 J=1,NT1
20    A(J)=0.
      DO 30 J=NLOW,NHIGH
30    A(J)=Y((J-NLOW)*NXY+ITR)
      DO 50 J=2,NT2-1
50    A(NT1+2-J)=CONJG(A(J))
      CALL FFTC(A,B,2,2*NT1,NT1,1,+1)
      XFAC=1./NT1
      IOUT=(I-1)*NXSAVE+K
      DO 60 J=1,NT
60    OUT(J,IOUT)=XFAC*REAL(A(J))
10    CONTINUE
      RETURN
      END
      subroutine spksynmod(vel,e,x,bfit,afit,refl,velz,ivt)
c     xmodeli2 de dvssyn.f. SEASEIS: modelo em memoria (REFL, VELZ) em
c     vez dos arquivos 'indata' e 'velscr'; IVT=1 lista VXT (xvevents),
c     IVT=3 AFIT ja preenchido pelo chamador (arquivo ONL de v(z)).
      real refl(*),velz(*)
      character*4 vxtc
      real vel(*),xx(10 000),filt(10 000)  
      integer ifax(4)
      real e(*),x(*) 
      real afit(*),bfit(*)
      complex a(10000),b(10000),c(10000)  
      common/velneed/nrecs,irec,lev,vinc,ipass,npass,nxv,lag,locv,ispj
      common/migneed/nx,nt,nz,itzr,vint,dt,nnz,ddz,beta,dx,apx,maxscr, 
     *rfc,alph1,gam1,frq,rfstir,itype,cyc,raa,padl,padr,zzy,ltm,ny,dy,  
     *frq1,frq2,alph2,rni,alph3,
     *vmin,vmax,seps,lev1,lev2,ifor,ilor,window,lunin,lunout
      character*4 vint,rfc
      dz=itzr*ddz  
      nxx=nx 
      sr=1000*dt 
      nwds=5000 
      ia=1 
      ib=2 
      if(ivt.ge.1)then
      if(ivt.eq.1)then
      call xvevents(e,vel,xx,x,nxx,ny,nt,vxtc,nwds,afit,bfit,dx)  
      call xvevents(e,vel,xx,x,nxx,ny,nt,vxtc,nwds,afit,bfit,dx)  
      end if
c     SEASEIS: ramos DVSVELS (VEL 0) e VELOCITY.VXT (VEL -1) removidos
      end if
      iprm=alph1 
      nnnz=iprm*nnz  
      nnt=iprm*nt
      dddz=ddz/iprm  
      ddt=dt/iprm
      rddt=1./ddt
      fnyq=1./(2.*ddt)  
      if(iprm.gt.1.or.frq2.le.fnyq)then  
      ant=nnt-1 
       call fac235(nnt,ifax)
       nt1=ifax(1)
      nt2=nt1/2+1  
      xx(1)=0. 
      xx(2)=frq1*ddt*nt1
      xx(3)=2.*xx(2) 
      xx(5)=frq2*ddt*nt1
      xx(4)=.667*xx(5) 
      if(xx(3).gt.xx(4))then 
      xx(3)=.5*(xx(3)+xx(4)) 
      xx(4)=xx(3)  
      end if 
      write(*,301)frq1,2.*frq1,.667*frq2,frq2
301   format(4f6.0,'   filter applied')
      vel(1)=0.  
      vel(2)=0.  
      vel(3)=1.  
      vel(4)=1.  
      vel(5)=0.  
      if(xx(2).eq.0.)then  
      vel(1)=1.  
      vel(2)=1.  
      end if 
      call xint2(filt,xx,vel,0.,1.,nt2,5)  
      end if 
c   ***   write velociies to file "velscr" for use in syn3d.
c   ***   writes to "velscr" are via pi2 (2 byte data storage)
c     SEASEIS: arquivos 'velscr' e 'indata' substituidos por VELZ e REFL
      call movlev(afit,xx,5)
      nxv=xx(1)
      nyv=xx(2)
      np=xx(3) 
      idpn1=xx(5)
      write(*,*)'nxv nyv np idpn1 ',nxv,nyv,np,idpn1
      nwds=2*np+3 
      irec=1
      do 10 i=1,ny
      do 10 k=1,nx
c   ***   create dense (nt*iprm) reflection coefficient series trace in 
c   ***   time and velocity function (nz) in depth.
c   ***   time trace will be filtered in time, converted to depth, and
c   ***   decimated to nnz samples.
      iy=min0(i,nyv)
      ix=k-idpn1+beta+1
      ix=max0(1, min0(nxv,ix) )
      iwa=4+((iy-1)*nxv+ix-1)*nwds
      call movlev(afit(iwa),xx,nwds)
      call clear(x,nnt)  
      ii=3 
111   continue
      if(xx(ii).eq.0. )then
      ii=ii+2 
      goto 111
      end if
      v1=xx(ii)
      v2=v1 
      z=0.
      z1=0. 
      z2=xx(ii+1)
      t=0.
      eps=.1
      itprm=itzr*iprm
      do 101 n=1,nz+1
      rvbar=0.
      do 104 l=1,itprm
102   continue
      if(ii.ge.2*np+2 .or. z .le. z2)goto 103
      if(xx(ii).eq.0. .or. z .gt. xx(ii+1))then
      ii=ii+2
      goto 102
      end if
      v1=v2 
      v2=xx(ii)
      z1=z2 
      z2=xx(ii+1)
      it=t*rddt+1
      x(it)=10000.*(v2-v1)/(v2+v1+eps)
103   continue
      vvel=v2
      rvel=1./vvel
      z=z+dddz 
      t=t+2.*dddz*rvel
      rvbar=rvbar+rvel
      vel((n-1)*itprm+l)=vvel
c   *** 
104   continue
      e(n)=itprm/rvbar
101   continue
      do 164 jv=1,nz
164   velz((irec-1)*nz+jv)=e(jv)
c   ***   filter reflection coefficient series 
      if(iprm.gt.1.or.frq2.le.fnyq)then  
c   ***   filter trace if iprm > 1  or frq2 < fnyq . 
c   ***   complex data 
      do 400 n=1,nnt  
400   a(n)=cmplx(x(n),0.)  
      do 410 n=nnt+1,nt1  
410   a(n)=0.  
c   ***   take fft 
      factor=-1. 
c     call dvsft(a,b,c,nt1,factor) 
      log2=alog(real(nt1))/alog(2.)
c     call fftl(a,log2,-1)
       call fftc(a,b,2,2*nt1,nt1,1,+1)
c   ***   apply filter 
      b(1)=filt(1)*a(1) 
      do 420 n=2,nt2 
      b(n)=filt(n)*a(n) 
420   b(nt1-n+2)=conjg(b(n)) 
c   ***   take inverse fft 
      factor=1./nt1  
c     call dvsft(b,a,c,nt1,factor) 
c     call fftl(b,log2,+1)
      call fftc(b,a,2,2*nt1,nt1,1,-1)
      do 425 n=1,nt1
      b(n)=factor*b(n)
425   continue
      end if
c   ***   convert time data trace to depth
      z=0.
      t=0.
      x(1)=real(b(1))
      it=1
      do 500 n=2,nnnz
      if(it.lt.nnt)then
      t=t+2.*dddz/vel(n-1)
      rt=t/ddt+1.
      it=rt
      f1=iprm*(rt-it)
      f2=iprm-f1
      x(n)=f2*real(b(it))+f1*real(b(it+1))
      else
      x(n)=0.
      end if
500   continue
c   ***   decimate data to original sample rate and put away
      do 502 j=1,nnz
      e(j)=x((j-1)*iprm+1)
502   continue
      do 214 jv=1,nnz
214   refl((irec-1)*nnz+jv)=e(jv)
      irec=irec+1
10     continue  
      return
      end
c
      SUBROUTINE SPKSYNVXT(VEL,DX,NXV,NP,AFIT,EE,X,BFIT,XX)
c     SEASEIS: uma linha y da lista VXT -> funcoes de velocidade densas
c     em x, com XVEVENTS (duas chamadas, como em xmodeli2).
c     Saida em AFIT: afit(1..3) = nxv, 1, np; depois, para cada DPN,
c     2*np+3 palavras: iy, ix, pares (v, z).
      REAL VEL(*),AFIT(*),EE(*),X(*),BFIT(*),XX(*)
      CHARACTER*4 VXTC
      NXX=0
      NY=1
      NT=0
      NWDS=50000
      CALL XVEVENTS(EE,VEL,XX,X,NXX,NY,NT,VXTC,NWDS,AFIT,BFIT,DX)
      CALL XVEVENTS(EE,VEL,XX,X,NXX,NY,NT,VXTC,NWDS,AFIT,BFIT,DX)
      NXV=NINT(AFIT(1))
      NP=NINT(AFIT(3))
      RETURN
      END
