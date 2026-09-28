c-----------------------------------------------------------------------
c     Migracao 2D do programa DVSMIG (SEISPAK), para o modulo SPK_DVSMIG2D.
c
c     SPKMIG2D : prepara tamanhos e memoria e chama SPKMIGW
c     SPKMIGW  : programa principal de dvsmig.f, caminho ny=1 (2D):
c                entrada/FFT, continuacao para baixo (APX 1, 4/5, 6/7,
c                >=11) e imageamento (mrg3 / sgimrg2r).
c                Diferencas em relacao ao original (marcadas SEASEIS:):
c                - dado, velocidade e imagem em memoria (antes: arquivos
c                  'indata', 'velscr' e transposicao trni2da)
c                - todas as frequencias num unico passo (nstepw=1)
c                - imagem em float (antes: gravada em 16 bits com pi2)
c     SPKVCOL  : velocidade por passo de profundidade para uma posicao,
c                trecho de xvel3d (dvsmig.f) sem alteracao
c-----------------------------------------------------------------------
      SUBROUTINE SPKMIG2D(D,NXSAVE,NT,SR,VXZ,NZ,ITZR,NNZ,DDZ,DX,APX,
     :                    IRFC,FRQ1,FRQ2,OUT,INFO,IERR)
c     D(NT,NXSAVE)    secao de entrada
c     VXZ(NZ,NXSAVE)  velocidade intervalar por passo (de SPKVCOL)
c     OUT(NNZ,NXSAVE) secao migrada
c     IRFC            0 = tempo (RFC NO), 1 = profundidade (RFC YES)
c     INFO(6)         saida: nx (FFT), nt1, nlow, nhigh, nz, memoria(MB)
c     IERR            0 ok, 1 banda de frequencias vazia, 2 sem memoria
      REAL D(NT,*),VXZ(NZ,*),OUT(NNZ,*)
      INTEGER INFO(6),IFAX(5)
      REAL, ALLOCATABLE :: WK(:),VEL(:)
      IERR=0
      NX=NXSAVE
      IF(APX .GT. 5. .AND. NXSAVE .GT. 4)THEN
        CALL FAC235(NX,IFAX)
        NX=IFAX(1)
      ENDIF
      DT=.001*SR
      NTFFT=1.25*NT
      CALL FAC235(NTFFT,IFAX)
      NT1=IFAX(1)
      NT2=NT1/2+1
      IFRQ1=2.*DT*NT2*FRQ1+1
      IFRQ2=2.*DT*NT2*FRQ2+1
      NLOW=MAX0(IFRQ1,1)
      NHIGH=MIN0(NT2,IFRQ2)
      INFO(1)=NX
      INFO(2)=NT1
      INFO(3)=NLOW
      INFO(4)=NHIGH
      INFO(5)=NZ
      IF(NHIGH.LT.NLOW)THEN
        IERR=1
        RETURN
      ENDIF
      NWS=NHIGH-NLOW+1
      NXWS=NX*NWS
c                                  mesmo layout de dvsmig.f (ny=1)
      ID=1+5*NXWS
      NWK=2*(ID+NXWS+1)+8*MAX0(NX,NWS)+2*NXWS+2*(NT1+NX)+4096
      NWK=MAX0(NWK,2*(2*NXWS+1+NT1)+4*NT1+4096)
      NVEL=MAX0(ITZR*NX,8*NT1+1024)
      INFO(6)=(4*(NWK+NVEL))/1000000+1
      ALLOCATE(WK(NWK),VEL(NVEL),STAT=IST)
      IF(IST.NE.0)THEN
        IERR=2
        RETURN
      ENDIF
      CALL SPKMIGW(WK,WK,VEL,D,VXZ,OUT,NXSAVE,NX,NT,NT1,SR,NZ,ITZR,
     :             NNZ,DDZ,DX,APX,IRFC,NLOW,NHIGH)
      DEALLOCATE(WK,VEL)
      RETURN
      END
c
      SUBROUTINE SPKMIGW(Y,S,VEL,D,VXZ,OUT,NXSAVE,NX,NT,NT1,SR,NZ,ITZR,
     :                   NNZ,DDZ,DX,APX,IRFC,NLOW,NHIGH)
c     Y (complexo) e S (real) sao o MESMO vetor de trabalho, como a
c     EQUIVALENCE (y(1),s(1)) de dvsmig.f
      COMPLEX Y(*)
      REAL S(*),VEL(*),D(NT,*),VXZ(NZ,*),OUT(NNZ,*)
      NY=1
      NXY=NX*NY
      DT=.001*SR
      DY=DX
      DZ=ITZR*DDZ
      FACTOR=2./NT1
      DELW=2.*3.1415927/NT1
      CON1=DZ/(.5*DT)
      CON2=(.5*DT/DY)**2
      CON4=ITZR
      NL=NLOW
      NH=NHIGH
      NWS=NH-NL+1
      NXWS=NX*NWS
      NXYWS=NXY*NWS
c   ***   compute y (complex) indicies  (ny=1, como em dvsmig.f)
      IY=1
      IZ=IY+NXWS
      IA=IZ+NXWS
      IB=IA+NXWS
      IC=IB+NXWS
      ID=IC+NXWS
      ISX=2*(ID+NXWS+1)+1
      IVX=ISX+NX
      IVY=IVX+NX
      ISY=IVY+NX
      IT1=ISY+NX
      IT2=IT1+NX
      IT=IT2+NX
      IS=IT+MAX0(NX,NWS)
c   ***   read data, fft and keep only frequencies nl..nh
      CALL CLEAR(Y,2*NXYWS)
      DO 1 K=1,NXSAVE
      I1=K
      DO 70 J=1,NT1
70    VEL(J)=0.
c     SEASEIS: traco vem da memoria (antes: gettrc)
      DO 69 J=1,NT
69    VEL(J)=D(J,K)
      IF(APX.LT.6. .AND. (K.LE.2 .OR. K.GE.NX-1))
     :   CALL CLEAR(VEL,NT)
      DO 71 J=1,NT1
71    Y(IA+J)=FACTOR*VEL(J)
      CALL FFTC(Y(IA+1),VEL,2,2*NT1,NT1,1,-1)
      DO 73 J=1,NWS
73    Y(I1+(J-1)*NXY)=Y(IA+NL-1+J)
1     CONTINUE
c   ***   begin downward continuation loop
      SGN=+1.
      DO 85 J=1,NXY
      S(IT1+J-1)=0.
85    S(IT2+J-1)=0.
      DO 2000 K=1,NZ
c     SEASEIS: velocidade do passo k vem da memoria (antes: gi2 'velscr')
c              tracos de padding usam a velocidade do ultimo traco
      DO 86 J=1,NXY
      JJ=MIN0(J,NXSAVE)
86    S(IS+J-1)=VXZ(K,JJ)
      DO 87 J=1,NXY
87    S(IVY+J-1)=CON2*S(IS+J-1)**2
      IF(IRFC.EQ.1)THEN
        DO 88 J=1,NXY
88      S(ISY+J-1)=CON1/S(IS+J-1)
      ELSE
        DO 89 J=1,NXY
89      S(ISY+J-1)=CON4
      ENDIF
      CALL RTRAN(S(ISY),S(ISX),NX,NY)
      CALL RTRAN(S(IVY),S(IVX),NX,NY)
c   ***   save complex field from step k-1 in y(iz)
      DO 80 J=1,NXYWS
80    Y(IZ+J-1)=Y(IY+J-1)
c   ***   process 2d migration
      IF(APX .EQ. 6. .OR. APX .EQ. 7.)THEN
        CALL PSPI(Y(IY),Y(IA),Y(IB),Y(IC),S(IS),NL,NH,NT1,NX,
     :            S(ISX),S(IVX),SGN)
      ENDIF
c   ***   if apx > 10 then use multi-velocity pspi, nvels = apx-10
      IF(APX .GE. 11.)THEN
        NVELS=NINT(APX-10.)
        CALL XMULTPSPI(Y(IY),Y(IA),Y(IB),Y(IC),S(IS),NL,NH,NT1,NX,
     :                 S(ISX),S(IVX),SGN,NVELS)
c   ***   apply koslov absorbing boundary for pspi
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
      IF(APX.EQ.4. .OR. APX.EQ.5.)THEN
        CALL CTRAN(Y(IY),Y(IA),NX,NWS)
        CALL DWN45B(Y(IA),Y(IA),Y(IY),Y(IB),Y(IC),Y(ID),S(IT),S(IS),
     :              S(IVX),S(ISX),NX,NL,NH,DELW,SGN)
        CALL CTRAN(Y(IA),Y(IY),NWS,NX)
      ENDIF
      IF(APX .EQ. 1)THEN
        CALL VSHIFT1(Y(IY),Y(IA),Y(IA),Y(IB),S(ISY),DELW,NL,NH,NX,SGN)
      ENDIF
c   ***   sgimrg2r images for depth migration, mrg3 for time migration
      IF(IRFC.EQ.1)THEN
        IF(APX.EQ.6..OR.APX.EQ.7.)THEN
          DO 90 J=1,NXY
90        S(ISY+J-1)=2.*S(ISY+J-1)
        ENDIF
        CALL SGIMRG2R(VEL,Y(IZ),Y(IY),S(ISY),S(IT1),S(IT2),
     :                NL,NH,NT1,NXY,ITZR)
      ELSE
        CALL MRG3(VEL,Y(IZ),Y(IY),S(ISY),NL,NH,NT1,NXY,ITZR)
      ENDIF
c     SEASEIS: imagem vai direto para OUT (antes: pi2 'indata' + trni2da)
      DO 901 JJ=1,ITZR
        M=(K-1)*ITZR+JJ
        IF(M.LE.NNZ)THEN
          DO 902 J=1,NXSAVE
902       OUT(M,J)=VEL((JJ-1)*NXY+J)
        ENDIF
901   CONTINUE
2000  CONTINUE
      RETURN
      END
c
      SUBROUTINE SPKVCOL(XX,NP,NZ,ITZR,DDZ,DT,IRFC,E)
c     XX(2*NP+3): xx(1),xx(2) = posicao; pares (velocidade, profundidade)
c                 a partir de xx(3), como em xvel3d / DVSVELS
c     E(NZ)       saida: velocidade intervalar media em cada passo
      REAL XX(*),E(*)
      ii=3
      eps=.1
111   continue
c     SEASEIS: evita laco infinito se todas as velocidades forem < eps
      if(ii .gt. 2*np+1)then
        do 112 n=1,nz
112     e(n)=0.
        return
      endif
      if(xx(ii) .lt. eps)then
      ii=ii+2
      goto 111
      end if
      v1=xx(ii)
      v2=v1
      z=0.
      z1=0.
      z2=xx(ii+1)
      eps=.1
      do 101 n=1,nz
      rvbar=0.
      do 104 l=1,itzr
102   continue
      if(ii.gt.2*np+2.or.z.le.z2)goto 103
      if(xx(ii).lt.eps .or. (z .gt. xx(ii+1).and.ii .lt. 2*np))then
      ii=ii+2
      goto 102
      end if
      v1=v2
      v2=xx(ii)
      z1=z2
      z2=xx(ii+1)
103   continue
      vvel=v2
      if(irfc.eq.1)then
      z=z+ddz
      else
      z=z+.5*vvel*dt
      end if
      rvbar=rvbar+1./vvel
104   continue
      e(n)=itzr/rvbar
101   continue
      return
      end
c
      SUBROUTINE SPKVXT(VEL,DX,NXV,NYV,NP,AFIT,EE,X,BFIT,XX)
c     Lista VXT (horizontes) -> funcoes de velocidade densas, como no
c     ramo VEL(1) > 0 de xvel3d: duas chamadas de VEVENTS.
c     Saida em AFIT: afit(1..3) = nxv, nyv, np; depois, para cada
c     posicao (iy,ix), 2*np+3 palavras: iy, ix, pares (v, z).
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
