c-----------------------------------------------------------------------
c     Migracao 3D do programa DVSMIG (SEISPAK), para o modulo SPK_DVSMIG3D.
c
c     SPKMIG3D : prepara tamanhos e memoria e chama SPKMIG3W
c     SPKMIG3W : programa principal de dvsmig.f, caminho ny>1 (3D):
c                entrada/FFT, continuacao para baixo (APX 4/5 dwn3d em
c                x e em y, APX 6/7 cosphs3d+refrac, APX >=11 mvphs2d),
c                bordas absorventes (bndry), deslocamento vertical
c                (vshift1) e imageamento (sgimrg2r / mrg3).
c                Diferencas em relacao ao original (marcadas SEASEIS:):
c                - dado, velocidade e imagem em memoria (antes: arquivos
c                  'indata', 'velscr' e transposicao trni2da)
c                - todas as frequencias num unico passo (nstepw=1)
c                - imagem em float (antes: gravada em 16 bits com pi2)
c     SPKVCOL  : velocidade por passo de profundidade para uma posicao,
c                trecho de xvel3d (dvsmig.f) sem alteracao
c     SPKM3VXT : uma linha y da lista VXT -> funcoes densas (XVEVENTS)
c-----------------------------------------------------------------------
      SUBROUTINE SPKMIG3D(D,NXSAVE,NYSAVE,NT,SR,VXZ,NZ,ITZR,NNZ,DDZ,
     :                    DX,DY,APX,IRFC,FRQ1,FRQ2,OUT,INFO,IERR)
c     D(NT,NXSAVE*NYSAVE)    dado de entrada (x rapido, y lento)
c     VXZ(NZ,NXSAVE*NYSAVE)  velocidade intervalar por passo (SPKVCOL)
c     OUT(NNZ,NXSAVE*NYSAVE) volume migrado
c     IRFC            0 = tempo (RFC NO), 1 = profundidade (RFC YES)
c     INFO(7)         saida: nx, ny (FFT), nt1, nlow, nhigh, nz, MB
c     IERR            0 ok, 1 banda de frequencias vazia, 2 sem memoria
      REAL D(NT,*),VXZ(NZ,*),OUT(NNZ,*)
      INTEGER INFO(7),IFAX(5)
      REAL, ALLOCATABLE :: WK(:),VEL(:)
      IERR=0
      NX=NXSAVE
      NY=NYSAVE
c   ***   set up size for pspi migration (nx,ny = 2**p1 * 3**p2 * 5**p3)
      IF(APX .GT. 5. .AND. NXSAVE .GT. 4)THEN
        CALL FAC235(NX,IFAX)
        NX=IFAX(1)
      ENDIF
      IF(APX .GT. 5. .AND. NYSAVE .GT. 4)THEN
        CALL FAC235(NY,IFAX)
        NY=IFAX(1)
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
      NXY=NX*NY
      NXYWS=NXY*NWS
c                                  mesmo layout de dvsmig.f (ny>1)
      IA=1+2*NXYWS
      ID=IA+3*NXY
      IS=2*(ID+NXY)+1+6*NXY
      NWK=IS+8*NXY+4*NT1+8192
      NWK=MAX0(NWK,2*(IA+NT1+1)+4*NT1+8192)
      NVEL=MAX0(ITZR*NXY,8*NT1+1024)
      INFO(7)=(4.*(REAL(NWK)+REAL(NVEL)))/1.E6+1
      ALLOCATE(WK(NWK),VEL(NVEL),STAT=IST)
      IF(IST.NE.0)THEN
        IERR=2
        RETURN
      ENDIF
      CALL SPKMIG3W(WK,WK,VEL,D,VXZ,OUT,NXSAVE,NYSAVE,NX,NY,NT,NT1,SR,
     :              NZ,ITZR,NNZ,DDZ,DX,DY,APX,IRFC,NLOW,NHIGH)
      DEALLOCATE(WK,VEL)
      RETURN
      END
c
      SUBROUTINE SPKMIG3W(Y,S,VEL,D,VXZ,OUT,NXSAVE,NYSAVE,NX,NY,NT,NT1,
     :                    SR,NZ,ITZR,NNZ,DDZ,DX,DY,APX,IRFC,NLOW,NHIGH)
c     Y (complexo) e S (real) sao o MESMO vetor de trabalho, como a
c     EQUIVALENCE (y(1),s(1)) de dvsmig.f
      COMPLEX Y(*)
      REAL S(*),VEL(*),D(NT,*),VXZ(NZ,*),OUT(NNZ,*)
      NXY=NX*NY
      DT=.001*SR
      IF(DY.EQ.0.)DY=DX
      DZ=ITZR*DDZ
      FACTOR=2./NT1
      DELW=2.*3.1415927/NT1
      CON1=DZ/(.5*DT)
      CON2=(.5*DT/DY)**2
      CON3=(DY/DX)**2
      CON4=ITZR
      NL=NLOW
      NH=NHIGH
      NWS=NH-NL+1
      NXYWS=NXY*NWS
c   ***   compute y (complex) indicies (ny>1, como em dvsmig.f)
      IY=1
      IZ=IY+NXYWS
      IA=IZ+NXYWS
      IB=IA+NXY
      IC=IB+NXY
      ID=IC+NXY
c   ***   s (real) indicies start at end of complex
      ISX=2*(ID+NXY)+1
      ISY=ISX+NXY
      IVX=ISY+NXY
      IVY=IVX+NXY
      IT1=IVY+NXY
      IT2=IT1+NXY
      IS=IT2+NXY
c   ***   read data, fft and keep only frequencies nl..nh
      CALL CLEAR(Y,2*NXYWS)
      DO 1 I=1,NYSAVE
      DO 1 K=1,NXSAVE
      I1=(I-1)*NX+K
      DO 70 J=1,NT1
70    VEL(J)=0.
c     SEASEIS: traco vem da memoria (antes: gettrc)
      DO 69 J=1,NT
69    VEL(J)=D(J,(I-1)*NXSAVE+K)
      IF(APX.LT.6. .AND. (K.LE.2 .OR. K.GE.NX-1))
     :   CALL CLEAR(VEL,NT)
      IF(APX.LT.6. .AND. NY.GE.5 .AND. (I.LE.2 .OR. I.GE.NY-1))
     :   CALL CLEAR(VEL,NT)
      DO 71 J=1,NT1
71    Y(IA+J)=FACTOR*VEL(J)
      CALL FFTC(Y(IA+1),VEL,2,2*NT1,NT1,1,-1)
      DO 73 J=1,NWS
73    Y(I1+(J-1)*NXY)=Y(IA+NL-1+J)
1     CONTINUE
c   ***   begin downward continuation loop
      SGN=+1.
c   ***   preset arrival times to zero for post stack migration.
      DO 85 J=1,NXY
      S(IT1+J-1)=0.
85    S(IT2+J-1)=0.
      DO 2000 K=1,NZ
c     SEASEIS: velocidade do passo k vem da memoria (antes: gi2 'velscr');
c              tracos de padding usam a velocidade do ultimo traco/linha
      CALL SPKM3V(S(IS),VXZ,K,NZ,NX,NY,NXSAVE,NYSAVE)
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
      DO 21 J=1,NXY
21    S(IVX+J-1)=CON3*S(IVX+J-1)
c   ***   save complex field from step k-1 in y(iz)
      DO 80 J=1,NXYWS
80    Y(IZ+J-1)=Y(IY+J-1)
      IF(APX.EQ.6. .OR. APX.EQ.7.)THEN
        CALL SPKM3V(S(IS),VXZ,K,NZ,NX,NY,NXSAVE,NYSAVE)
        DO J=1,NXY
          S(IVY+J-1)=(S(IS+J-1)*DT)**2
          S(ISY+J-1)=DZ/(S(IS+J-1)*DT)
          S(IVX+J-1)=(S(IS+J-1)*DT)**2
          S(ISX+J-1)=DZ/(S(IS+J-1)*DT)
        ENDDO
      ENDIF
      DO 100 I=1,NWS
      W=(NL+I-2)*DELW
      I1=(I-1)*NXY+1
      IF(APX.EQ.4. .OR. APX.EQ.5.)THEN
        IF(W.GT.0.)THEN
c   ***   downward continue for x-direction.
          CALL CTRAN(Y(I1),Y(IA),NX,NY)
          CALL DWN3D(Y(IA),Y(IA),Y(I1),Y(IB),Y(IC),Y(ID),W,S(IVX),
     :               S(ISX),NX,NY,SGN)
          CALL CTRAN(Y(IA),Y(I1),NY,NX)
c   ***   downward continue for y-direction.
          CALL DWN3D(Y(I1),Y(I1),Y(IA),Y(IB),Y(IC),Y(ID),W,S(IVY),
     :               S(ISY),NY,NX,SGN)
        ENDIF
      ENDIF
      IF(APX .EQ. 6. .OR. APX .EQ. 7.)THEN
        NVELS=1
c   ***   source term
        CALL COSPHS3D(Y(I1),Y(IA),Y(IB),S(IS),NX,1,NY,DX,DX,
     :                DY,W,S(ISY),S(IVY),-SGN,NVELS)
        CALL REFRAC(Y(I1),S(ISY),W,NXY,SGN)
c   ***   receiver term
        CALL COSPHS3D(Y(I1),Y(IA),Y(IB),S(IS),NX,1,NY,DX,DX,
     :                DY,W,S(ISX),S(IVX),+SGN,NVELS)
        CALL REFRAC(Y(I1),S(ISX),W,NXY,SGN)
      ENDIF
      IF(APX .GE. 11.)THEN
        NVELS=NINT(APX-10.)
        CALL MVPHS2D(Y(I1),Y(IB),Y(IC),S(IS),NX,NY,DX,DY,
     :               W,S(ISY),S(IVY),SGN,NVELS)
      ENDIF
c   ***   apply koslov attenuating absorbing boundary
      LBX=10
      LBY=10
      CALL BNDRY(Y(I1),NX,NY,LBX,LBY)
100   CONTINUE
c     SEASEIS: em APX 6/7 o deslocamento vertical ja foi aplicado por
c              REFRAC (fonte + receptor); o VSHIFT1 do original aplicava
c              o deslocamento mais uma vez e a imagem saia rasa (2/3).
      IF(.NOT.(APX.EQ.6. .OR. APX.EQ.7.))THEN
        CALL VSHIFT1(Y(IY),Y(IA),Y(IB),Y(IC),S(ISY),DELW,NL,NH,NXY,SGN)
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
          DO 902 I=1,NYSAVE
          DO 902 J=1,NXSAVE
902       OUT(M,(I-1)*NXSAVE+J)=VEL((JJ-1)*NXY+(I-1)*NX+J)
        ENDIF
901   CONTINUE
2000  CONTINUE
      RETURN
      END
c
      SUBROUTINE SPKM3V(V,VXZ,K,NZ,NX,NY,NXSAVE,NYSAVE)
c     SEASEIS: velocidade do passo K em V(NX,NY); padding repete a borda
      REAL V(NX,NY),VXZ(NZ,*)
      DO 10 I=1,NY
      II=MIN0(I,NYSAVE)
      DO 10 J=1,NX
      JJ=MIN0(J,NXSAVE)
10    V(J,I)=VXZ(K,(II-1)*NXSAVE+JJ)
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
      SUBROUTINE SPKM3VXT(VEL,DX,NXV,NP,AFIT,EE,X,BFIT,XX)
c     SEASEIS: uma linha y da lista VXT -> funcoes de velocidade densas
c     em x, com XVEVENTS (duas chamadas, como em xmodeli2 / SPK_DVSSYN3D).
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
